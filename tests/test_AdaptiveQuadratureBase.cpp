// API contract of AdaptiveQuadratureBase (SEN-05, SEN-07, ROB-01 of the test plan).

#include "Check.hpp"
#include "Tolerances.hpp"

#include <LNIT/AdaptiveQuadratures.hpp>

#include <cmath>
#include <limits>

using namespace LNIT::tests;

namespace
{

template<class Quadrature>
void checkFreshObject(const char* name)
{
	std::printf("-- fresh %s\n", name);
	// SEN-05: state of a quadrature object before any integration
	const Quadrature quad{};
	CHECK(quad.hasConverged() == false);
	CHECK(quad.getNits() == 0);
}

template<class Quadrature>
void checkBounds(const char* name)
{
	std::printf("-- bounds %s\n", name);
	constexpr double nan = std::numeric_limits<double>::quiet_NaN();
	const auto f = [](const double x) { return std::exp(x); };
	const double exact = std::exp(1.) - 1.; // integral of exp over [0, 1]

	Quadrature quad;
	quad.setTol(1.e-12);

	const double forward = double(quad.integrate(f, 0., 1.));
	CHECK_CLOSE(forward, exact, tol::REF, 0.);
	CHECK(quad.hasConverged());

	// SEN-07 / ROB-01: reversed bounds follow the orientation convention of the integral
	const double backward = double(quad.integrate(f, 1., 0.));
	CHECK_CLOSE(backward, -exact, tol::REF, 0.);
	CHECK(quad.hasConverged());

	// equal bounds: zero, converged
	CHECK(double(quad.integrate(f, 0.5, 0.5)) == 0.);
	CHECK(quad.hasConverged());

	// infinite bounds delegate to the infinite-domain methods (exp integrates to 1 over ]-inf, 0])
	constexpr double inf = std::numeric_limits<double>::infinity();
	CHECK_CLOSE(double(quad.integrate(f, -inf, 0.)),  1., 1.e-12*tol::TOL_USER_FACTOR, 0.);
	CHECK(quad.hasConverged());
	CHECK_CLOSE(double(quad.integrate(f, 0., -inf)), -1., 1.e-12*tol::TOL_USER_FACTOR, 0.);
	CHECK(quad.hasConverged());
	CHECK(double(quad.integrate(f, inf, inf)) == 0.);

	// NaN bounds propagate, no plausible value is returned silently
	CHECK(std::isnan(double(quad.integrate(f, nan, 1.))));
	CHECK(!quad.hasConverged());
	CHECK(std::isnan(double(quad.integrate(f, 0., nan))));
	CHECK(!quad.hasConverged());
}

} // namespace

int main()
{
	using GL   = LNIT::GaussLegendreAdaptiveQuadrature<double, long double>;
	using CC   = LNIT::ClenshawCurtisAdaptiveQuadrature<double, long double>;
	using Hyb  = LNIT::ClenshawCurtisHybridAdaptiveQuadrature<double, long double>;
	using GLCC = LNIT::GLCCAdaptiveQuadrature<double, long double>;

	checkFreshObject<GL>("GaussLegendre");
	checkFreshObject<CC>("ClenshawCurtis");
	checkFreshObject<Hyb>("ClenshawCurtisHybrid");
	checkFreshObject<GLCC>("GLCC");

	checkBounds<GL>("GaussLegendre");
	checkBounds<CC>("ClenshawCurtis");
	checkBounds<Hyb>("ClenshawCurtisHybrid");
	checkBounds<GLCC>("GLCC");

	return report("test_AdaptiveQuadratureBase");
}
