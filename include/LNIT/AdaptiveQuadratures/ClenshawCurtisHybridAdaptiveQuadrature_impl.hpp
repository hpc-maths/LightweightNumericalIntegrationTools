#ifndef LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_IMPL_HPP
#define LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/ClenshawCurtisHybridAdaptiveQuadrature.hpp>

namespace LNIT
{

//// explicit template instanciations ////

extern template class ClenshawCurtisHybridAdaptiveQuadrature<float, float>;
extern template class ClenshawCurtisHybridAdaptiveQuadrature<float, double>;
extern template class ClenshawCurtisHybridAdaptiveQuadrature<float, long double>;

extern template class ClenshawCurtisHybridAdaptiveQuadrature<double, double>;
extern template class ClenshawCurtisHybridAdaptiveQuadrature<double, long double>;

extern template class ClenshawCurtisHybridAdaptiveQuadrature<long double, long double>;

//// method implementations ////

template<typename T, typename TT> template<class Function>
constexpr auto ClenshawCurtisHybridAdaptiveQuadrature<T,TT>::estimateIntegralImpl(const Function& f, const Scalar& xmin, const Scalar& xmax) -> std::pair<LongScalar, LongScalar>
{	
	using std::abs;
	
	const auto fx = s_xi | std::views::transform([&f, &xmin, &xmax](const Scalar& xi) -> LongScalar
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return static_cast<LongScalar>(f(x)); 
	});
	std::ranges::copy(fx, std::begin(m_fx));
	
	const LongScalar I1 = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(s_wi.begin(),		    s_wi.end(),		     m_fx.begin(), LongScalar{});
	const LongScalar I2 = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(s_alternateWi.begin(), s_alternateWi.end(), m_fx.begin(), LongScalar{});

	return std::make_pair(I1, abs(I1 - I2));
}

template<typename T, typename TT> template<class Function>
constexpr auto ClenshawCurtisHybridAdaptiveQuadrature<T,TT>::integrateImpl(const Function& f, const Scalar& xmin, const Scalar& xmax) const -> std::invoke_result_t<Function, Scalar>
{
	const auto fx = s_xi | std::views::transform([&f, &xmin, &xmax](const Scalar& xi) -> std::invoke_result_t<Function, Scalar>
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return f(x); 
	});
	
	return LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(s_wi.begin(), s_wi.end(), fx.begin(), std::invoke_result_t<Function, Scalar>{});
}

} // namespace LNIT

#endif // LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_IMPL_HPP
