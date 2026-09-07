#ifndef LNIT_ADAPTIVE_QUADRATURE_BASE_IMPL_HPP
#define LNIT_ADAPTIVE_QUADRATURE_BASE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase.hpp>
#include <LNIT/GaussLaguerreQuadrature.hpp>
#include <LNIT/misc/Numeric.hpp>

#include <cassert>
#include <algorithm>
#include <ranges>
#include <cmath>

#include <fmt/core.h>

namespace LNIT
{

template<class Derived>
AdaptiveQuadratureBase<Derived>::AdaptiveQuadratureBase(const Size& maxIt, const Scalar& relativeTol, const Scalar& absoluteTol) 
	: m_maxIt(maxIt)
	, m_relativeTol(relativeTol)
	, m_absoluteTol(absoluteTol) 
{ 
	m_intervals.reserve(maxIt); 
	m_subIntergrals.reserve(maxIt); 
	m_subIntergralsErr.reserve(maxIt);
}

template<class Derived>
constexpr void AdaptiveQuadratureBase<Derived>::resetState()
{
	m_hasConverged = false;
	m_it = 0;
	m_intervals.clear();
	m_subIntergrals.clear();
	m_subIntergralsErr.clear();
}

template<class Derived> template<class Function> 
auto AdaptiveQuadratureBase<Derived>::integrate(const Function& f, const Scalar& xmin, const Scalar& xmax) -> LongScalar
{	
	using std::ceil;
	using std::abs;
	using std::isfinite;
	using std::isnan;
	using std::isinf;
	
	using const_Iterator = typename std::vector<LongScalar>::const_iterator;
	
	// bounds: see the documentation of this method
	if (isnan(xmin) or isnan(xmax))  { resetState(); return NumTraits<LongScalar>::NaN; }
	if (xmin > xmax)                 { return -integrate(f, xmax, xmin); }
	if (xmin == xmax)                { resetState(); m_hasConverged = true; return LongScalar{}; }
	if (isinf(xmin) and isinf(xmax)) { return integrate(f); }
	if (isinf(xmin))                 { return integrateLeftInfinite(f, xmax); }
	if (isinf(xmax))                 { return integrateRightInfinite(f, xmin); }
	
	resetState();

	LongScalar res;
	LongScalar estimatedErr;
	
	if (m_out) { fmt::print(m_out, "#NumericalIntegrator addapting quadrature over [{}, {}]\n", xmin, xmax); }
	if (m_out) { fmt::print(m_out, "#Iteration integral estimated_error relative_tol absolute_tol\n"); }
	
	const Size N = Size(ceil(getMaxDeltaX(xmin, xmax)));
	
	m_intervals.reserve(N);
	m_subIntergrals.reserve(N);
	m_subIntergralsErr.reserve(N);
	
	for (Size i=0; i!=N; ++i)
	{
		const Scalar x_i   = xmin + Scalar(i)*(xmax - xmin) / Scalar(N);
		const Scalar x_ip1 = xmin + Scalar(i+1)*(xmax - xmin) / Scalar(N);
		
		std::tie(res, estimatedErr) = estimateIntegral(f, x_i, x_ip1);
		
		m_intervals.emplace_back(x_i, x_ip1);
		m_subIntergrals.push_back(res);
		m_subIntergralsErr.push_back(estimatedErr);
	}
	
	for (m_it=0; m_it!=m_maxIt; ++m_it)
	{
		const LongScalar I   = getEstimatedIntegral();
		const LongScalar err = getEstimatedError();

		if (m_out) { fmt::print(m_out, "{} {:10.4e} {:10.4e} {:10.4e} {:10.4e}\n", m_it, Scalar(I), Scalar(err), Scalar(abs(I))*m_relativeTol, Scalar(m_absoluteTol)); }

		if (not isfinite(I)) { return I; }
		if (err < abs(I)*LongScalar(m_relativeTol) or err < LongScalar(m_absoluteTol)) { m_hasConverged = true; return I; }

		// we find the interval over which the integral is the least accurate
		const const_Iterator maxErrIt = std::ranges::max_element(m_subIntergralsErr);
		const Size maxErrIdx = Size(std::ranges::distance(m_subIntergralsErr.begin(), maxErrIt));
		// we split it in two
		const auto [a, b] = m_intervals[maxErrIdx]; // copy: the entry is overwritten just below
		
		Scalar midPoint = std::midpoint(a, b); // non-const because I want to move it when I do not need it.
		// first interval
		m_intervals[maxErrIdx] = Interval(a, midPoint);
		std::tie(m_subIntergrals[maxErrIdx], m_subIntergralsErr[maxErrIdx]) = estimateIntegral(f, a, midPoint);
		// second interval
		std::tie(res, estimatedErr) = estimateIntegral(f, midPoint, b);
		m_intervals.emplace_back(std::move(midPoint), b);
		m_subIntergrals.push_back(res);
		m_subIntergralsErr.push_back(estimatedErr);
	}
	return getEstimatedIntegral();
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::integrateLeftInfinite(const Function& f, const Scalar& xmax) -> LongScalar
{
	using std::isfinite;
	using std::abs;
	
	
	GaussLaguerreQuadrature<Scalar,LongScalar> gLaguerreQuad;
	
	Scalar xmin = -1;
	LongScalar leftIntegral = gLaguerreQuad.integrateLeftInfinite(f, xmin);
	while (isfinite(leftIntegral) and abs(leftIntegral) >= NumTraits<LongScalar>::epsilon)
	{
		xmin *= 2;
		leftIntegral = gLaguerreQuad.integrateLeftInfinite(f, xmin);
	}
	if (not isfinite(leftIntegral)) { resetState(); return NumTraits<LongScalar>::NaN; }
	return integrate(f, xmin, xmax);
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::integrateRightInfinite(const Function& f, const Scalar& xmin) -> LongScalar
{
	using std::isfinite;
	using std::abs;
	
	GaussLaguerreQuadrature<Scalar,LongScalar> gLaguerreQuad;
	
	Scalar xmax = 1;
	LongScalar rightIntegral = gLaguerreQuad.integrateRightInfinite(f, xmax);
	while (isfinite(rightIntegral) and abs(rightIntegral) >= NumTraits<LongScalar>::epsilon)
	{
		xmax *= 2;
		rightIntegral = gLaguerreQuad.integrateRightInfinite(f, xmax);
	}
	if (not isfinite(rightIntegral)) { resetState(); return NumTraits<LongScalar>::NaN; }
	return integrate(f, xmin, xmax);
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::integrate(const Function& f) -> LongScalar
{
	using std::isfinite;
	using std::abs;
	
	GaussLaguerreQuadrature<Scalar,LongScalar> gLaguerreQuad;
	
	Scalar xmin = -1;
	LongScalar leftIntegral = gLaguerreQuad.integrateLeftInfinite(f, xmin);
	while (isfinite(leftIntegral) and abs(leftIntegral) >= NumTraits<LongScalar>::epsilon)
	{
		xmin *= 2;
		leftIntegral = gLaguerreQuad.integrateLeftInfinite(f, xmin);
	}
	
	if (not isfinite(leftIntegral)) { resetState(); return NumTraits<LongScalar>::NaN; }
	return integrateRightInfinite(f, xmin);
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::remapAndIntegrate(const Function& f) -> LongScalar
{		
	using std::isnan;
	
	constexpr Scalar eps = {};
	
	const auto fref = [&f](const Scalar t) -> LongScalar
	{			
		const LongScalar fx = f(t / (1 - t*t));
		const Scalar dxdt = (1 + t*t) / ((1 - t*t)*(1 - t*t));
			
		return isnan(fx*dxdt)
			? LongScalar{}
			: fx*dxdt;	
	};
	
	return integrate(fref, -1 + eps, 1 - eps);
}

template<class Derived>  template<class Function> 
auto AdaptiveQuadratureBase<Derived>::integrateWithoutAdaptation(const Function& f) const -> std::invoke_result_t<Function, Scalar>
{
	const auto localIntegrals = m_intervals | std::views::transform([&self = derived(), &f](const Interval& interval) -> std::invoke_result_t<Function, Scalar>
	{
		return self.integrateImpl(f, interval.first, interval.second);
	});
	
	return std::reduce(std::ranges::begin(localIntegrals), std::ranges::end(localIntegrals));
}

} // namespace LNIT

#endif // LNIT_ADAPTIVE_QUADRATURE_BASE_IMPL_HPP
