#ifndef LNIT_GAUSS_LAGUERRE_QUADRATURE_IMPL_HPP
#define LNIT_GAUSS_LAGUERRE_QUADRATURE_IMPL_HPP

#include <numeric>
#include <fmt/core.h>
#include <ranges>

namespace LNIT
{

//// explicit template instanciations ////

extern template class GaussLaguerreQuadrature<double, double>;
extern template class GaussLaguerreQuadrature<double, long double>;

extern template class GaussLaguerreQuadrature<long double, long double>;

//// method implementations ////

template<typename Scalar, typename LongScalar> template<class Function> 
constexpr LongScalar GaussLaguerreQuadrature<Scalar, LongScalar>::integrateLeftInfinite(const Function& f, const Scalar& a) const
{
	const auto fx = s_xi | std::views::transform([&f, &a](const Scalar& x) -> LongScalar
	{
		return f(a - x); 
	});

	return std::inner_product(s_wi.begin(), s_wi.end(), fx.begin(), LongScalar{});	
}

template<typename Scalar, typename LongScalar> template<class Function> 
constexpr LongScalar GaussLaguerreQuadrature<Scalar, LongScalar>::integrateRightInfinite(const Function& f, const Scalar& a) const
{
	const auto fx = s_xi | std::views::transform([&f, &a](const Scalar& x) -> LongScalar
	{
		return f(x + a); 
	});

	return std::inner_product(s_wi.begin(), s_wi.end(), fx.begin(), LongScalar{});	
}

} // namespace LNIT

#endif // LNIT_GAUSS_LAGUERRE_QUADRATURE_IMPL_HPP
