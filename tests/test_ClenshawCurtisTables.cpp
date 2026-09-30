// Clenshaw-Curtis tables (SEN-04, TAB-01..04, TAB-06 of the test plan).
//
// A closed Clenshaw-Curtis rule with n+1 nodes cos(k pi / n), n even, integrates
// every polynomial of degree <= n+1 exactly on [-1, 1] (degree n by construction,
// n+1 because odd monomials vanish by symmetry). The strongest oracle is therefore
// polynomial exactness at the declared degree: it needs no external reference and
// catches both wrong weights and indexing mistakes.

#include "Check.hpp"
#include "Tolerances.hpp"

#include <LNIT/AdaptiveQuadratures.hpp>

#include <array>
#include <cmath>
#include <span>
#include <string>

using namespace LNIT::tests;

namespace
{

/// exact value of the integral of x^k over [-1, 1]
constexpr double monomialIntegral(const unsigned k) { return (k % 2 == 1) ? 0. : 2. / double(k + 1); }

/// exact value of the integral of the Chebyshev polynomial T_n over [-1, 1], n even
constexpr double chebyshevIntegral(const unsigned n) { return 2. / (1. - double(n)*double(n)); }

double chebyshev(const unsigned n, const double x) { return std::cos(double(n)*std::acos(x)); }

template<class F>
double apply(std::span<const double> xi, std::span<const double> wi, const F& f)
{
	long double s = 0;
	for (std::size_t i = 0; i != xi.size(); ++i) { s += static_cast<long double>(wi[i])*static_cast<long double>(f(xi[i])); }
	return static_cast<double>(s);
}

/// Checks one rule against its declared exactness degree.
/// @param degree  highest polynomial degree the rule must integrate exactly
/// @param mainRule  false for error-estimation companion rules whose zero weights are intentional
void checkRule(const std::string& name, std::span<const double> xi, std::span<const double> wi, const unsigned degree, const bool mainRule = true)
{
	std::printf("-- %s (%zu nodes, declared exactness degree %u)\n", name.c_str(), xi.size(), degree);
	const std::size_t n = xi.size();

	// TAB-03: sum of weights
	CHECK_CLOSE(apply(xi, wi, [](double) { return 1.; }), 2., tol::EXACT, 0.);

	// TAB-04: symmetry, sorted nodes in [-1, 1]
	for (std::size_t i = 0; i != n; ++i)
	{
		CHECK(wi[i] == wi[n - 1 - i]);
		CHECK(xi[i] == -xi[n - 1 - i]);
		CHECK(xi[i] >= -1. && xi[i] <= 1.);
		if (i + 1 != n) { CHECK(xi[i] > xi[i + 1]); }
	}

	// TAB-06: positivity of the main rules
	if (mainRule) { for (const double w : wi) { CHECK(w > 0.); } }

	// TAB-01: exactness for every monomial up to the declared degree
	for (unsigned k = 0; k <= degree; ++k)
	{
		const double rtol = (k <= 20) ? tol::EXACT : tol::EXACT_HI;
		const double got  = apply(xi, wi, [k](double x) { return std::pow(x, double(k)); });
		if (!closeTo(got, monomialIntegral(k), rtol, rtol))
		{
			std::fprintf(stderr, "%s: x^%u integrates to %.17g, exact %.17g\n", name.c_str(), k, got, monomialIntegral(k));
			++g_failures;
		}
	}

	// SEN-04: the last Chebyshev mode is where a missing 1/2 factor shows up
	if (degree % 2 == 1 && mainRule)
	{
		const unsigned nCheb = degree - 1; // even
		CHECK_CLOSE(apply(xi, wi, [nCheb](double x) { return chebyshev(nCheb, x); }), chebyshevIntegral(nCheb), tol::EXACT_HI, tol::EXACT_HI);
	}

	// TAB-02: the declared degree is not underestimated -> the next even degree must fail clearly.
	// Chebyshev polynomials are used because the leading Chebyshev coefficient of x^k is 2^(1-k):
	// for the 33-node rules the monomial error would be of the order of the round-off.
	{
		const unsigned k   = degree + 1 + (degree % 2); // next even degree above `degree`
		const double   got = apply(xi, wi, [k](double x) { return chebyshev(k, x); });
		CHECK(std::abs(got - chebyshevIntegral(k)) > 100.*tol::EXACT_HI);
	}
}

template<class Array>
std::span<const double> view(const Array& a) { return std::span<const double>(a.data(), a.size()); }

template<class Array, std::size_t Step>
std::array<double, (std::tuple_size_v<Array> - 1) / Step + 1> subsample(const Array& a)
{
	std::array<double, (std::tuple_size_v<Array> - 1) / Step + 1> out{};
	for (std::size_t i = 0; i != out.size(); ++i) { out[i] = a[i*Step]; }
	return out;
}

} // namespace

int main()
{
	using CC   = LNIT::ClenshawCurtisAdaptiveQuadrature<double, long double>;
	using Hyb  = LNIT::ClenshawCurtisHybridAdaptiveQuadrature<double, long double>;
	using GLCC = LNIT::GLCCAdaptiveQuadrature<double, long double>;

	// nested rules of ClenshawCurtisAdaptiveQuadrature: 33, 17 and 9 nodes
	const auto xi17 = subsample<decltype(CC::s_xi), 2>(CC::s_xi);
	const auto xi09 = subsample<decltype(CC::s_xi), 4>(CC::s_xi);
	checkRule("CC33", view(CC::s_xi), view(CC::s_wi33), 33);
	checkRule("CC17", view(xi17),     view(CC::s_wi17), 17);
	checkRule("CC09", view(xi09),     view(CC::s_wi09), 9);

	// hybrid: 13-node CC rule and its alternate companion rule
	checkRule("Hybrid13",    view(Hyb::s_xi), view(Hyb::s_wi),          13);
	checkRule("HybridAlt",   view(Hyb::s_xi), view(Hyb::s_alternateWi), 9, /*mainRule=*/false);

	// GLCC keeps its own copy of the 33-node CC table
	checkRule("GLCC-CC33", view(GLCC::s_xi_cc), view(GLCC::s_wi_cc), 33);
	for (std::size_t i = 0; i != 33; ++i) { CHECK(GLCC::s_wi_cc[i] == CC::s_wi33[i]); CHECK(GLCC::s_xi_cc[i] == CC::s_xi[i]); }

	return report("test_ClenshawCurtisTables");
}
