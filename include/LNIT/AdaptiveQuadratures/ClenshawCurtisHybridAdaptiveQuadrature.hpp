#ifndef LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_HPP
#define LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_HPP

#include <cmath>

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase.hpp>
#include <LNIT/misc/Numeric.hpp>

namespace LNIT
{

template<typename T, typename TT> class ClenshawCurtisHybridAdaptiveQuadrature;

template<typename T, typename TT> 
struct AdaptiveQuadratureTraits< ClenshawCurtisHybridAdaptiveQuadrature<T, TT> >
{
	using Size       = unsigned int;
	using Scalar     = T;
	using LongScalar = TT;
};

/**
 * @brief Adaptive quadrature using Clenshaw–Curtis rules and another quadrature for error estimation.
 *
 * This class implements adaptive quadrature using Clenshaw–Curtis (CC) rules with 13 nodes to evalutate the integral
 * and a 11 nodes quadrature to estimate the error.
 * More precisely, the 11 nodes quadrature is constructed by removing the leftmost, rightmost and mid-point CC node (-1, 1 and 0)
 * 
 * The error estimate is computed as \f$|I_{\text{CC}} - I_{11}|\f$ and is of order 11.
 * 
 * @tparam T Floating point type for integration (e.g., double).
 * @tparam TT Higher precision type for accumulation (e.g., long double).
 */
template<typename T, typename TT=T>
class ClenshawCurtisHybridAdaptiveQuadrature : public AdaptiveQuadratureBase< ClenshawCurtisHybridAdaptiveQuadrature<T,TT> >
{
	using Base = AdaptiveQuadratureBase< ClenshawCurtisHybridAdaptiveQuadrature<T,TT> >;
public:
	using Size       = Base::Size;       ///<  @brief Type for iteration counters.
	using Scalar     = Base::Scalar;     ///<  @brief Floating point type for integration (e.g., double).
	using LongScalar = Base::LongScalar; ///<  @brief Higher precision type for accumulation (e.g., long double).

	/**
	 * @brief Estimate integral and error on [xmin, xmax].
	 * @tparam Function Callable with signature Scalar f(Scalar).
	 * @param f Function to integrate.
	 * @param xmin Lower bound.
	 * @param xmax Upper bound.
	 * @return Pair (integral, estimated error).
	 */
	template<class Function> constexpr std::pair<LongScalar, LongScalar> estimateIntegralImpl(const Function& f, const Scalar& xmin, const Scalar& xmax);
	
	template<class Function> constexpr std::invoke_result_t<Function, Scalar> integrateImpl(const Function& f, const Scalar& xmin, const Scalar& xmax) const;
	
	constexpr Scalar getMaxDeltaXImpl(const Scalar& xmin, const Scalar& xmax) const { return (xmax - xmin)*misc::maxDiff(std::span{s_xi}); } 
#ifdef LNIT_TESTING
public: // quadrature tables exposed to the unit tests only
#else
private:
#endif // LNIT_TESTING
	std::array<LongScalar, 13> m_fx;

	static constexpr std::array<Scalar, 13> s_wi = {
		Scalar(0.0069930069930069930070), Scalar(0.066057424952074394517 ), Scalar(0.13154253154253154253  ), Scalar(0.18476338476338476338  ), Scalar(0.22697302697302697303  ), 
		Scalar(0.25267569378104433860  ), Scalar(0.26198986198986198986  ), Scalar(0.25267569378104433860  ), Scalar(0.22697302697302697303  ), Scalar(0.18476338476338476338  ), 
		Scalar(0.13154253154253154253  ), Scalar(0.066057424952074394517 ), Scalar(0.0069930069930069930070)};
	
	static constexpr std::array<Scalar, 13> s_alternateWi = {
		Scalar{}, Scalar(0.096656546636198874454), Scalar(0.055555555555555555556), Scalar(0.32275132275132275132), Scalar(0.026984126984126984127), Scalar(0.49805244807279583454 ), 
		Scalar{}, Scalar(0.49805244807279583454 ), Scalar(0.026984126984126984127), Scalar(0.32275132275132275132), Scalar(0.055555555555555555556), Scalar(0.096656546636198874454), 
		Scalar{}};
		
	static constexpr std::array<Scalar, 13> s_xi = {
		 Scalar( 1.0000000000000000000 ), Scalar( 0.96592582628906828675), Scalar( 0.86602540378443864676), Scalar( 0.70710678118654752440), Scalar( 0.50000000000000000000), 
		 Scalar( 0.25881904510252076235), Scalar(                       ), Scalar(-0.25881904510252076235), Scalar(-0.50000000000000000000), Scalar(-0.70710678118654752440), 
		 Scalar(-0.86602540378443864676), Scalar(-0.96592582628906828675), Scalar(-1.0000000000000000000 )};
};

} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/ClenshawCurtisHybridAdaptiveQuadrature_impl.hpp>

#endif // LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_HPP
