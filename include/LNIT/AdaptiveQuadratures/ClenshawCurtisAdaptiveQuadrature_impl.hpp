#ifndef LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_IMPL_HPP
#define LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/ClenshawCurtisAdaptiveQuadrature.hpp>

#include <LNIT/misc/StridedView.hpp>

namespace LNIT
{

//// explicit template instanciations ////
extern template class ClenshawCurtisAdaptiveQuadrature<float, float>;
extern template class ClenshawCurtisAdaptiveQuadrature<float, double>;
extern template class ClenshawCurtisAdaptiveQuadrature<float, long double>;

extern template class ClenshawCurtisAdaptiveQuadrature<double, double>;
extern template class ClenshawCurtisAdaptiveQuadrature<double, long double>;

extern template class ClenshawCurtisAdaptiveQuadrature<long double, long double>;

//// method implementations ////

template<typename T, typename TT> template<class Function>
constexpr auto ClenshawCurtisAdaptiveQuadrature<T,TT>::estimateIntegralImpl(const Function& f, const Scalar& xmin, const Scalar& xmax) -> std::pair<LongScalar, LongScalar>
{	
	using std::abs;
	
	const auto fx = s_xi | std::views::transform([&f, &xmin, &xmax](const Scalar& xi) -> LongScalar
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return static_cast<LongScalar>(f(x)); 
	});
	std::ranges::copy(fx, std::begin(m_fx33));
	
#if ((defined(_MSVC_LANG) && _MSVC_LANG >= 202302L) || __cplusplus >= 202302L)
	const auto fx17 = m_fx33 | std::views::stride(2);
	const auto fx09 = m_fx33 | std::views::stride(4);
#else
	const auto fx17 = misc::stride<2>(m_fx33);
	const auto fx09 = misc::stride<4>(m_fx33);
#endif // if using c++23

	const LongScalar I33 = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(s_wi33.begin(), s_wi33.end(), m_fx33.begin(), LongScalar{});
	const LongScalar I17 = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(s_wi17.begin(), s_wi17.end(),   fx17.begin(), LongScalar{});
	const LongScalar I09 = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(s_wi09.begin(), s_wi09.end(),   fx09.begin(), LongScalar{});
	
	const LongScalar err1 = abs(I33 - I17);
	const LongScalar err2 = abs(I33 - I09);
	
	return std::make_pair(I33, err2 == LongScalar{} ? LongScalar{} : err1*(err1 / err2)*(err1 / err2));
	//~ return std::make_pair(I33, err2 == LongScalar{} ? LongScalar{} : err2*(err1 / err2)*(err1 / err2));
}

template<typename T, typename TT> template<class Function>
constexpr auto ClenshawCurtisAdaptiveQuadrature<T,TT>::integrateImpl(const Function& f, const Scalar& xmin, const Scalar& xmax) const -> std::invoke_result_t<Function, Scalar>
{
	const auto fx = s_xi | std::views::transform([&f, &xmin, &xmax](const Scalar& xi) -> std::invoke_result_t<Function, Scalar>
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return f(x); 
	});
	
	return LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(s_wi33.begin(), s_wi33.end(), fx.begin(), std::invoke_result_t<Function, Scalar>{});
}

} // namespace LNIT

#endif // LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_IMPL_HPP
