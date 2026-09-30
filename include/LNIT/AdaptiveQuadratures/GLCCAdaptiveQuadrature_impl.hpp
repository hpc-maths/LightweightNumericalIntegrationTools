#ifndef LNIT_GL_CC_ADAPTIVE_QUADRATURE_IMPL_HPP
#define LNIT_GL_CC_ADAPTIVE_QUADRATURE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/GLCCAdaptiveQuadrature.hpp>

namespace LNIT
{

//// explicit template instanciations ////

extern template class GLCCAdaptiveQuadrature<float, float>;
extern template class GLCCAdaptiveQuadrature<float, double>;
extern template class GLCCAdaptiveQuadrature<float, long double>;

extern template class GLCCAdaptiveQuadrature<double, double>;
extern template class GLCCAdaptiveQuadrature<double, long double>;

extern template class GLCCAdaptiveQuadrature<long double, long double>;

//// method implementations ////

template<typename T, typename TT> template<class Function>
constexpr auto GLCCAdaptiveQuadrature<T,TT>::estimateIntegralImpl(const Function& f, const Scalar& xmin, const Scalar& xmax) const -> std::pair<LongScalar, LongScalar>
{
	using std::abs;
	
	const auto fx_gl = s_xi_gl | std::views::transform([&f, &xmin, &xmax](const Scalar& xi) -> LongScalar
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return static_cast<LongScalar>(f(x)); 
	});
	
	const auto fx_cc = s_xi_cc | std::views::transform([&f, &xmin, &xmax](const Scalar& xi) -> LongScalar
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return static_cast<LongScalar>(f(x)); 
	});
	
	const LongScalar I_gl = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(s_wi_gl.begin(), s_wi_gl.end(), fx_gl.begin(), LongScalar{});
	const LongScalar I_cc = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(s_wi_cc.begin(), s_wi_cc.end(), fx_cc.begin(), LongScalar{});
	
	return std::make_pair(I_gl, abs(I_cc - I_gl));
}

template<typename T, typename TT> template<class Function>
constexpr auto GLCCAdaptiveQuadrature<T,TT>::integrateImpl(const Function& f, const Scalar& xmin, const Scalar& xmax) const -> std::invoke_result_t<Function, Scalar>
{
	const auto fx = s_xi_gl | std::views::transform([&f, &xmin, &xmax](const Scalar& xi) -> std::invoke_result_t<Function, Scalar>
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return f(x); 
	});
	
	return LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(s_wi_gl.begin(), s_wi_gl.end(), fx.begin(), std::invoke_result_t<Function, Scalar>{});
}

} // namespace LNIT

#endif // LNIT_GL_CC_ADAPTIVE_QUADRATURE_IMPL_HPP

