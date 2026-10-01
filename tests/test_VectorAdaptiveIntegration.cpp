// Vector-valued estimation (estimateIntegrals) and the globally adaptive vector driver (integrateVector).

#include "Check.hpp"
#include "Tolerances.hpp"

#include <LNIT/AdaptiveQuadratures.hpp>

#include <cmath>
#include <numbers>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace LNIT::tests;

namespace
{

using LD = long double;

/// x^j exp(-x^2/2), j = 0..M-1
void gaussianMoments(const LD& x, std::span<LD> values)
{
	LD xj = 1;
	const LD g = std::exp(-x*x/2);
	for (auto& v : values) { v = xj*g; xj *= x; }
}

/// exact integral of x^j exp(-x^2/2) over R: sqrt(2 pi) (j-1)!! for even j, 0 for odd j
LD gaussianMomentExact(const std::size_t j)
{
	if (j % 2 == 1) { return 0; }
	LD m = std::sqrt(2*std::numbers::pi_v<LD>);
	for (std::size_t k = 1; k < j; k += 2) { m *= LD(k); }
	return m;
}

// estimateIntegrals must give exactly what one estimateIntegral call per component gives
template<class Quadrature>
void checkBitwiseAgainstScalar(const char* name)
{
	std::printf("-- estimateIntegrals vs estimateIntegral, %s\n", name);
	constexpr std::size_t M = 9;
	Quadrature quad;
	for (const auto& [a, b] : {std::pair<LD, LD>{-1, 2}, {0.25L, 0.5L}, {-40, -3}, {1e3L, 1e3L + 1e-3L}})
	{
		std::vector<LD> I(M), E(M);
		quad.estimateIntegrals(gaussianMoments, a, b, std::span<LD>(I), std::span<LD>(E));
		for (std::size_t j = 0; j < M; ++j)
		{
			const auto fj = [j](const LD& x) { std::vector<LD> v(M); gaussianMoments(x, v); return v[j]; };
			const auto [Ij, Ej] = quad.estimateIntegral(fj, a, b);
			CHECK(I[j] == Ij);
			CHECK(E[j] == Ej);
		}
	}
}

template<class Quadrature>
void checkEstimateIntegralsContract(const char* name)
{
	std::printf("-- estimateIntegrals contract, %s\n", name);
	Quadrature quad;
	// the integrand is evaluated once per node, whatever the number of components
	std::size_t calls = 0;
	const auto counting = [&calls](const LD& x, std::span<LD> v) { ++calls; gaussianMoments(x, v); };
	std::vector<LD> I1(1), E1(1), I9(9), E9(9);
	quad.estimateIntegrals(counting, LD(0), LD(1), std::span<LD>(I1), std::span<LD>(E1));
	const std::size_t callsForOne = calls;
	calls = 0;
	quad.estimateIntegrals(counting, LD(0), LD(1), std::span<LD>(I9), std::span<LD>(E9));
	CHECK(calls == callsForOne);
	CHECK(I9[0] == I1[0]);

	// zero components: nothing is evaluated
	calls = 0;
	quad.estimateIntegrals(counting, LD(0), LD(1), std::span<LD>(), std::span<LD>());
	CHECK(calls == 0);

	// mismatched output sizes are rejected
	std::vector<LD> I(3), E(2);
	bool thrown = false;
	try { quad.estimateIntegrals(gaussianMoments, LD(0), LD(1), std::span<LD>(I), std::span<LD>(E)); }
	catch (const std::invalid_argument&) { thrown = true; }
	CHECK(thrown);
}

template<class Quadrature>
void checkDriverGaussianMoments(const char* name, const bool conservativeEstimate)
{
	std::printf("-- integrateVector, gaussian moments, %s\n", name);
	constexpr std::size_t M = 11;
	Quadrature quad;
	// a single initial interval: the driver has to find the peak by itself
	const std::vector<std::pair<LD, LD>> mesh{{-40, 40}};
	LNIT::VectorIntegrationOptions<LD> opt;
	opt.relativeTol = 1e-14L;
	const auto res = LNIT::integrateVector(quad, gaussianMoments, M, std::span(mesh), opt);
	CHECK(res.converged);
	CHECK(res.intervals.size() > 1);
	CHECK(res.nRuleCalls == 2*res.intervals.size() - 1); // one initial interval, each bisection adds two
	for (std::size_t j = 0; j < M; ++j)
	{
		const LD exact = gaussianMomentExact(j);
		// the truncation of R to [-40, 40] is far below long double precision
		const LD trueError = std::abs(res.integrals[j] - exact);
		CHECK(trueError <= tol::TOL_USER_FACTOR*opt.relativeTol*res.absIntegrals[j]);
		// The rules that estimate the error as the difference of two rules (Hybrid, GLCC) are
		// conservative. GaussLegendre and ClenshawCurtis extrapolate it, err1 (err1/err2)^2, and
		// can underestimate it (by 4x for x^10 exp(-x^2/2) with ClenshawCurtis): for them, only
		// the accuracy above is required.
		if (conservativeEstimate) { CHECK(trueError <= res.errors[j] + 64*std::numeric_limits<LD>::epsilon()*res.absIntegrals[j]); }
		CHECK(res.errors[j] <= opt.relativeTol*res.absIntegrals[j]);
	}
	// the final mesh is a sorted partition of [-40, 40]
	CHECK(res.intervals.front().first == -40);
	CHECK(res.intervals.back().second == 40);
	for (std::size_t i = 0; i + 1 < res.intervals.size(); ++i) { CHECK(res.intervals[i].second == res.intervals[i + 1].first); }
}

void checkDriverOptions()
{
	std::printf("-- integrateVector options and errors\n");
	using Hyb = LNIT::ClenshawCurtisHybridAdaptiveQuadrature<LD, LD>;
	Hyb quad;
	const std::vector<std::pair<LD, LD>> mesh{{-40, -10}, {-10, 10}, {10, 40}};
	LNIT::VectorIntegrationOptions<LD> opt;
	opt.relativeTol = 1e-15L;

	// maxIntervals caps the refinement and reports non-convergence
	opt.maxIntervals = 4;
	const auto capped = LNIT::integrateVector(quad, gaussianMoments, 5, std::span(mesh), opt);
	CHECK(!capped.converged);
	CHECK(capped.intervals.size() <= 5);

	// a roundoff floor above every error estimate: nothing is refined, the floors are reported as error
	opt.maxIntervals = 20000;
	const auto floorAll = [](LD, LD, std::size_t, LD) { return LD(100); }; // coarse errors reach ~8 on [-10, 10]
	const auto floored = LNIT::integrateVector(quad, gaussianMoments, 5, std::span(mesh), opt, floorAll);
	CHECK(floored.converged);
	CHECK(floored.intervals.size() == mesh.size());
	CHECK(floored.nRuleCalls == mesh.size());
	for (std::size_t j = 0; j < 5; ++j) { CHECK(floored.errors[j] == LD(100*mesh.size())); }

	// invalid meshes are rejected
	bool thrown = false;
	try { LNIT::integrateVector(quad, gaussianMoments, 5, std::span<const std::pair<LD, LD>>(), opt); }
	catch (const std::invalid_argument&) { thrown = true; }
	CHECK(thrown);
	thrown = false;
	const std::vector<std::pair<LD, LD>> reversed{{1, 0}};
	try { LNIT::integrateVector(quad, gaussianMoments, 5, std::span(reversed), opt); }
	catch (const std::invalid_argument&) { thrown = true; }
	CHECK(thrown);
}

// A narrow peak away from the origin. The initial mesh must have a breakpoint at the peak: on
// [-10, 10] alone, no node of the first interval sees a peak of width 1e-3, the estimated error is 0
// and the driver would stop at once with a zero integral. Placing the initial mesh is the caller's job.
// The error normalisation per component keeps the odd moment, whose integral is small compared with
// the integral of its absolute value, under control.
void checkDriverShiftedPeak()
{
	std::printf("-- integrateVector, narrow shifted peak\n");
	using Hyb = LNIT::ClenshawCurtisHybridAdaptiveQuadrature<LD, LD>;
	Hyb quad;
	constexpr LD mu = 3, s = 1e-3L;
	const auto f = [](const LD& x, std::span<LD> v)
	{
		const LD g = std::exp(-(x - mu)*(x - mu)/(2*s*s));
		v[0] = g; v[1] = (x - mu)*g; v[2] = (x - mu)*(x - mu)*g;
	};
	const std::vector<std::pair<LD, LD>> mesh{{-10, mu}, {mu, 10}};
	LNIT::VectorIntegrationOptions<LD> opt;
	opt.relativeTol = 1e-13L;
	const auto res = LNIT::integrateVector(quad, f, 3, std::span(mesh), opt);
	const LD m0 = s*std::sqrt(2*std::numbers::pi_v<LD>);
	CHECK(res.converged);
	CHECK_CLOSE(double(res.integrals[0]), double(m0), tol::REF, 0.);
	CHECK(std::abs(res.integrals[1]) <= tol::TOL_USER_FACTOR*opt.relativeTol*res.absIntegrals[1]);
	CHECK_CLOSE(double(res.integrals[2]), double(s*s*m0), tol::REF, 0.);
}

} // namespace

int main()
{
	using GL   = LNIT::GaussLegendreAdaptiveQuadrature<LD, LD>;
	using CC   = LNIT::ClenshawCurtisAdaptiveQuadrature<LD, LD>;
	using Hyb  = LNIT::ClenshawCurtisHybridAdaptiveQuadrature<LD, LD>;
	using GLCC = LNIT::GLCCAdaptiveQuadrature<LD, LD>;

	checkBitwiseAgainstScalar<GL>("GaussLegendre");
	checkBitwiseAgainstScalar<CC>("ClenshawCurtis");
	checkBitwiseAgainstScalar<Hyb>("ClenshawCurtisHybrid");
	checkBitwiseAgainstScalar<GLCC>("GLCC");

	checkEstimateIntegralsContract<GL>("GaussLegendre");
	checkEstimateIntegralsContract<CC>("ClenshawCurtis");
	checkEstimateIntegralsContract<Hyb>("ClenshawCurtisHybrid");
	checkEstimateIntegralsContract<GLCC>("GLCC");

	checkDriverGaussianMoments<GL>("GaussLegendre", false);
	checkDriverGaussianMoments<CC>("ClenshawCurtis", false);
	checkDriverGaussianMoments<Hyb>("ClenshawCurtisHybrid", true);
	checkDriverGaussianMoments<GLCC>("GLCC", true);

	checkDriverOptions();
	checkDriverShiftedPeak();

	return report("test_VectorAdaptiveIntegration");
}
