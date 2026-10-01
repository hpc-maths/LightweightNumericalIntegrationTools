#ifndef LNIT_VECTOR_ADAPTIVE_INTEGRATION_HPP
#define LNIT_VECTOR_ADAPTIVE_INTEGRATION_HPP

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace LNIT
{

/**
 * @brief Options of integrateVector().
 */
template<typename LongScalar>
struct VectorIntegrationOptions
{
	LongScalar  relativeTol  = LongScalar(1e-13); ///< target: error_j <= relativeTol * absIntegral_j for every component j
	std::size_t maxIntervals = 20000;             ///< refinement stops (not converged) when this many intervals are alive
};

/**
 * @brief Result of integrateVector().
 */
template<typename Scalar, typename LongScalar>
struct VectorIntegrationResult
{
	std::vector<LongScalar> integrals;    ///< sum over the intervals of the integral of each component (compensated sum)
	std::vector<LongScalar> absIntegrals; ///< sum over the intervals of |integral| of each component (normalisation of the tolerance)
	std::vector<LongScalar> errors;       ///< estimated error of each component: quadrature errors plus roundoff floors
	std::vector<std::pair<Scalar, Scalar>> intervals; ///< final mesh, sorted
	std::size_t nRuleCalls = 0;           ///< number of estimateIntegrals() calls
	bool converged = false;               ///< true if the tolerance is met for every component
};

/**
 * @brief Globally adaptive integration of a vector-valued function over a given initial mesh.
 *
 * Every interval is integrated with rule.estimateIntegrals(), so the integrand is evaluated once per node
 * for all the components. The interval with the largest normalised error max_j E_j / A_j is bisected until
 *   sum_intervals E_j <= relativeTol * A_j   for every component j,
 * where A_j is the sum over the intervals of |integral of component j|. Normalising by A_j rather than
 * by |integral_j| keeps components whose integral nearly cancels (odd moments) under control.
 *
 * Roundoff floor: roundoffFloor(a, b, j, K) returns the roundoff level of the integral K of component j
 * over [a, b]. An error estimate at or below this level cannot be reduced by refinement: it is not
 * refined any further, but the floor is added to the reported error.
 *
 * @tparam Rule Adaptive quadrature providing estimateIntegrals().
 * @param rule Rule used on every interval.
 * @param f Callable with signature void f(const Scalar& x, std::span<LongScalar> values).
 * @param nComponents Number of components of f.
 * @param initialMesh Initial intervals (a < b); they do not need to be sorted.
 * @param options Tolerance and maximum number of intervals.
 * @param roundoffFloor Callable with signature LongScalar(Scalar a, Scalar b, std::size_t j, LongScalar K).
 *
 * @throw std::invalid_argument if the initial mesh is empty or contains an interval with a >= b.
 */
template<CAdaptiveQuadrature Rule, class Function, class RoundoffFloor>
VectorIntegrationResult<ScalarFor<Rule>, LongScalarFor<Rule>> integrateVector(
	Rule& rule, const Function& f, const std::size_t nComponents,
	std::span<const std::pair<ScalarFor<Rule>, ScalarFor<Rule>>> initialMesh,
	const VectorIntegrationOptions<LongScalarFor<Rule>>& options,
	const RoundoffFloor& roundoffFloor)
{
	using Scalar     = ScalarFor<Rule>;
	using LongScalar = LongScalarFor<Rule>;
	using std::abs;
	using std::max;

	if (initialMesh.empty()) { throw std::invalid_argument("LNIT::integrateVector: empty initial mesh"); }
	for (const auto& [a, b] : initialMesh)
	{
		if (!(a < b)) { throw std::invalid_argument("LNIT::integrateVector: every initial interval must satisfy a < b"); }
	}

	struct Piece
	{
		Scalar a, b;
		std::vector<LongScalar> K, E, N; // integral, error (0 when under the floor), roundoff floor
		bool alive;
	};

	VectorIntegrationResult<Scalar, LongScalar> result;
	std::vector<Piece> pieces;
	std::vector<LongScalar> Etot(nComponents, LongScalar{}), Atot(nComponents, LongScalar{});

	const auto evaluate = [&](const Scalar a, const Scalar b) -> Piece
	{
		Piece piece{a, b, std::vector<LongScalar>(nComponents), std::vector<LongScalar>(nComponents), std::vector<LongScalar>(nComponents), true};
		rule.estimateIntegrals(f, a, b, std::span<LongScalar>(piece.K), std::span<LongScalar>(piece.E));
		++result.nRuleCalls;
		for (std::size_t j = 0; j < nComponents; ++j)
		{
			piece.N[j] = roundoffFloor(a, b, j, piece.K[j]);
			if (piece.E[j] <= piece.N[j]) { piece.E[j] = LongScalar{}; }
		}
		return piece;
	};
	const auto account = [&](const Piece& piece, const LongScalar sign)
	{
		for (std::size_t j = 0; j < nComponents; ++j)
		{
			Etot[j] += sign*piece.E[j];
			Atot[j] += sign*abs(piece.K[j]);
		}
	};
	const auto reference = [&](const std::size_t j) { return max(Atot[j], std::numeric_limits<LongScalar>::min()); };
	const auto priority = [&](const Piece& piece)
	{
		LongScalar key{};
		for (std::size_t j = 0; j < nComponents; ++j) { key = max(key, piece.E[j]/reference(j)); }
		return key;
	};

	for (const auto& [a, b] : initialMesh)
	{
		pieces.push_back(evaluate(a, b));
		account(pieces.back(), LongScalar(1));
	}
	// Max-heap on the normalised error. Keys are computed with the A_j known at push time: they go
	// slightly stale as A_j evolves, which only changes the refinement order, never the stopping test.
	std::priority_queue<std::pair<LongScalar, std::size_t>> heap;
	for (std::size_t i = 0; i < pieces.size(); ++i) { heap.emplace(priority(pieces[i]), i); }

	std::size_t nAlive = pieces.size();
	while (true)
	{
		bool done = true;
		for (std::size_t j = 0; j < nComponents; ++j)
		{
			if (Etot[j] > options.relativeTol*reference(j)) { done = false; break; }
		}
		if (done) { result.converged = true; break; }
		if (nAlive >= options.maxIntervals || heap.empty()) { break; }
		const auto [key, index] = heap.top();
		heap.pop();
		if (key == LongScalar{}) { break; } // every remaining error is under its roundoff floor
		if (!pieces[index].alive) { continue; }

		pieces[index].alive = false;
		--nAlive;
		account(pieces[index], LongScalar(-1));
		const Scalar a = pieces[index].a, b = pieces[index].b, m = Scalar(0.5)*(a + b);
		for (const auto& [x0, x1] : {std::pair<Scalar, Scalar>{a, m}, std::pair<Scalar, Scalar>{m, b}})
		{
			pieces.push_back(evaluate(x0, x1));
			account(pieces.back(), LongScalar(1));
			heap.emplace(priority(pieces.back()), pieces.size() - 1);
			++nAlive;
		}
	}

	// Final sums from scratch: Kahan summation of the integrals, the running totals are not reused.
	result.integrals.assign(nComponents, LongScalar{});
	result.absIntegrals.assign(nComponents, LongScalar{});
	result.errors.assign(nComponents, LongScalar{});
	std::vector<LongScalar> compensation(nComponents, LongScalar{});
	for (const Piece& piece : pieces)
	{
		if (!piece.alive) { continue; }
		result.intervals.emplace_back(piece.a, piece.b);
		for (std::size_t j = 0; j < nComponents; ++j)
		{
			const LongScalar y = piece.K[j] - compensation[j];
			const LongScalar t = result.integrals[j] + y;
			compensation[j] = (t - result.integrals[j]) - y;
			result.integrals[j] = t;
			result.absIntegrals[j] += abs(piece.K[j]);
			result.errors[j] += piece.E[j] + piece.N[j];
		}
	}
	std::sort(result.intervals.begin(), result.intervals.end());
	return result;
}

/**
 * @brief integrateVector() without roundoff floor.
 */
template<CAdaptiveQuadrature Rule, class Function>
VectorIntegrationResult<ScalarFor<Rule>, LongScalarFor<Rule>> integrateVector(
	Rule& rule, const Function& f, const std::size_t nComponents,
	std::span<const std::pair<ScalarFor<Rule>, ScalarFor<Rule>>> initialMesh,
	const VectorIntegrationOptions<LongScalarFor<Rule>>& options = {})
{
	using Scalar     = ScalarFor<Rule>;
	using LongScalar = LongScalarFor<Rule>;
	return integrateVector(rule, f, nComponents, initialMesh, options,
	                       [](Scalar, Scalar, std::size_t, LongScalar) { return LongScalar{}; });
}

} // namespace LNIT

#endif // LNIT_VECTOR_ADAPTIVE_INTEGRATION_HPP
