// Precision of the rule tables in every Scalar type (TAB-07 of the test plan).
//
// The tables are long double literals with 36 significant digits. A rule instantiated
// with Scalar = T must then be exact to a few epsilons of T, including long double, where
// double-precision literals would stop at ~1e-17. Checked through the public API:
// polynomial exactness on [-1, 1] for the adaptive rules, closed-form integrals for
// Gauss-Hermite and Gauss-Laguerre. The exactness of the error-estimation companions is
// checked with 80 digits by tools/generate_quadrature_tables.py.

#include "Check.hpp"

#include <LNIT/AdaptiveQuadratures.hpp>
#include <LNIT/GaussHermiteQuadrature.hpp>
#include <LNIT/GaussLaguerreQuadrature.hpp>

#include <cmath>
#include <limits>
#include <numbers>

using namespace LNIT::tests;

namespace
{

template<class T>
const char* typeName() { return sizeof(T) == sizeof(float) ? "float" : sizeof(T) == sizeof(double) ? "double" : "long double"; }

/// Integral of x^k over [-1, 1]
template<class T>
T monomialIntegral(const unsigned k) { return (k % 2 == 1) ? T(0) : T(2)/T(k + 1); }

/// @param degree exactness degree of the rule
template<class Quadrature>
void checkAdaptiveRule(const char* name, const unsigned degree)
{
	using T = LNIT::ScalarFor<Quadrature>;
	std::printf("-- %s<%s>\n", name, typeName<T>());
	const T eps = std::numeric_limits<T>::epsilon();
	Quadrature quad;
	for (unsigned k = 0; k <= degree; ++k)
	{
		const auto f = [k](const T x) { return std::pow(x, T(k)); };
		const T I = T(quad.estimateIntegral(f, T(-1), T(1)).first);
		const T exact = monomialIntegral<T>(k);
		// a few epsilons of T, relative to the integral of |x^k| (odd monomials cancel)
		const T scale = T(2)/T(k + 1);
		if (!(std::abs(I - exact) <= T(64)*eps*scale))
		{
			std::fprintf(stderr, "%s<%s>: x^%u integrates to %.21Lg, exact %.21Lg\n", name, typeName<T>(), k,
			             static_cast<long double>(I), static_cast<long double>(exact));
			++g_failures;
		}
	}
}

template<class T>
void checkGaussHermite()
{
	std::printf("-- GaussHermiteQuadrature<%s>\n", typeName<T>());
	const T eps = std::numeric_limits<T>::epsilon();
	const LNIT::GaussHermiteQuadrature<T, T> quad;
	// int x^{2m} exp(-x^2) = Gamma(m + 1/2)
	for (unsigned m = 0; m <= 10; ++m)
	{
		const T I = quad.integrate([m](const T x) { return std::pow(x, T(2*m))*std::exp(-x*x); });
		const T exact = std::tgamma(T(m) + T(0.5));
		CHECK(std::abs(I - exact) <= T(256)*eps*exact);
	}
	// shifted Gaussian: int exp(-(x - 1/2)^2) = sqrt(pi)
	const T I = quad.integrate([](const T x) { return std::exp(-(x - T(0.5))*(x - T(0.5))); });
	CHECK(std::abs(I - std::sqrt(std::numbers::pi_v<T>)) <= T(256)*eps*std::sqrt(std::numbers::pi_v<T>));
}

template<class T>
void checkGaussLaguerre()
{
	std::printf("-- GaussLaguerreQuadrature<%s>\n", typeName<T>());
	const T eps = std::numeric_limits<T>::epsilon();
	const LNIT::GaussLaguerreQuadrature<T, T> quad;
	// int_0^inf x^k exp(-x) = k!
	T factorial = 1;
	for (unsigned k = 0; k <= 20; ++k)
	{
		if (k > 0) { factorial *= T(k); }
		// x^k exp(-x) as exp(k log x - x): x^k alone overflows float at the last nodes (x ~ 115)
		const T I = quad.integrateRightInfinite([k](const T x) { return std::exp(T(k)*std::log(x) - x); }, T(0));
		CHECK(std::abs(I - factorial) <= T(256)*eps*factorial);
	}
	// int_{-inf}^0 exp(x) = 1
	CHECK(std::abs(quad.integrateLeftInfinite([](const T x) { return std::exp(x); }, T(0)) - T(1)) <= T(256)*eps);
}

template<class T>
void checkAll()
{
	checkAdaptiveRule<LNIT::GaussLegendreAdaptiveQuadrature<T, T>>("GaussLegendreAdaptiveQuadrature", 29);        // 15 Gauss nodes
	checkAdaptiveRule<LNIT::ClenshawCurtisAdaptiveQuadrature<T, T>>("ClenshawCurtisAdaptiveQuadrature", 33);      // 33 nodes
	checkAdaptiveRule<LNIT::ClenshawCurtisHybridAdaptiveQuadrature<T, T>>("ClenshawCurtisHybridAdaptiveQuadrature", 13); // 13 nodes
	checkAdaptiveRule<LNIT::GLCCAdaptiveQuadrature<T, T>>("GLCCAdaptiveQuadrature", 29);                          // Gauss-Legendre 15
	checkGaussHermite<T>();
	checkGaussLaguerre<T>();
}

} // namespace

int main()
{
	checkAll<float>();
	checkAll<double>();
	checkAll<long double>();
	return report("test_TablePrecision");
}
