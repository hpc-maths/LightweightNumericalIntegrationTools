#ifndef LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_IMPL_HPP
#define LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/GaussLegendreAdaptiveQuadrature.hpp>

namespace LNIT
{

//// explicit template instanciations ////

extern template class GaussLegendreAdaptiveQuadrature<float, float>;
extern template class GaussLegendreAdaptiveQuadrature<float, double>;
extern template class GaussLegendreAdaptiveQuadrature<float, long double>;

extern template class GaussLegendreAdaptiveQuadrature<double, double>;
extern template class GaussLegendreAdaptiveQuadrature<double, long double>;

extern template class GaussLegendreAdaptiveQuadrature<long double, long double>;

//// method implementations ////

template<typename T, typename TT> template<class Function>
constexpr auto GaussLegendreAdaptiveQuadrature<T,TT>::estimateIntegralImpl(const Function& f, const Scalar& xmin, const Scalar& xmax) -> std::pair<LongScalar, LongScalar>
{
	using std::abs;
	
	const auto fx = s_xi | std::views::transform([&f, &xmin, &xmax](const Scalar& xi) -> LongScalar
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return static_cast<LongScalar>(f(x)); 
	});
	std::ranges::copy(fx, std::begin(m_fx15));
	
	const LongScalar I15 = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(std::ranges::begin(s_wi15), std::ranges::end(s_wi15), std::ranges::begin(m_fx15), LongScalar{});
	const LongScalar I14 = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(std::ranges::begin(s_wi14), std::ranges::end(s_wi14), std::ranges::begin(m_fx15), LongScalar{});
	const LongScalar I06 = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(std::ranges::begin(s_wi06), std::ranges::end(s_wi06), std::ranges::begin(m_fx15), LongScalar{});

	const LongScalar err1 = abs(I15 - I14);
	const LongScalar err2 = abs(I15 - I06);

	return std::make_pair(I15, err2 == LongScalar{} ? LongScalar{} : err1*(err1 / err2)*(err1 / err2));
	//~ return std::make_pair(I15, err2 == LongScalar{} ? LongScalar{} : err2*(err1 / err2)*(err1 / err2));
}

template<typename T, typename TT> template<class Function>
constexpr auto GaussLegendreAdaptiveQuadrature<T,TT>::integrateImpl(const Function& f, const Scalar& xmin, const Scalar& xmax) const -> std::invoke_result_t<Function, Scalar>
{
	const auto fx = s_xi | std::views::transform([&f, &xmin, &xmax](const Scalar& xi) -> std::invoke_result_t<Function, Scalar>
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return f(x); 
	});
	
	return LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(std::ranges::begin(s_wi15), std::ranges::end(s_wi15), fx.begin(), std::invoke_result_t<Function, Scalar>{});
}

} // namespace LNIT

#endif // LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_IMPL_HPP

