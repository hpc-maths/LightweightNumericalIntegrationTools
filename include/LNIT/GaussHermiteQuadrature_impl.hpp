#ifndef LNIT_GAUSS_HERMITE_QUADRATURE_IMPL_HPP
#define LNIT_GAUSS_HERMITE_QUADRATURE_IMPL_HPP

#include <numeric>
#include <fmt/core.h>
#include <ranges>

namespace LNIT
{

//// explicit template instanciations ////

extern template class GaussHermiteQuadrature<float, float>;
extern template class GaussHermiteQuadrature<float, double>;
extern template class GaussHermiteQuadrature<float, long double>;

extern template class GaussHermiteQuadrature<double, double>;
extern template class GaussHermiteQuadrature<double, long double>;

extern template class GaussHermiteQuadrature<long double, long double>;

//// method implementations ////

template<typename Scalar, typename LongScalar> template<class Function> 
constexpr LongScalar GaussHermiteQuadrature<Scalar, LongScalar>::integrate(const Function& f) const
{
	const auto fx = s_xi | std::views::transform([&f](const Scalar& x) -> LongScalar
	{
		return static_cast<LongScalar>(f(x)); 
	});

	return std::inner_product(s_wi.begin(), s_wi.end(), fx.begin(), LongScalar{});	
}

} // namespace LNIT

#endif // LNIT_GAUSS_HERMITE_QUADRATURE_IMPL_HPP
