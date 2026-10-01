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
		Scalar(0.00699300699300699300699300699300699301L),  Scalar(0.0660574249520743945174839281743238182L),   Scalar(0.131542531542531542531542531542531543L),
		Scalar(0.184763384763384763384763384763384763L),    Scalar(0.226973026973026973026973026973026973L),    Scalar(0.252675693781044338601249190558794915L),
		Scalar(0.261989861989861989861989861989861990L),    Scalar(0.252675693781044338601249190558794915L),    Scalar(0.226973026973026973026973026973026973L),
		Scalar(0.184763384763384763384763384763384763L),    Scalar(0.131542531542531542531542531542531543L),    Scalar(0.0660574249520743945174839281743238182L),
		Scalar(0.00699300699300699300699300699300699301L)};
	
	static constexpr std::array<Scalar, 13> s_alternateWi = {
		Scalar{},                                         Scalar(0.0966565466361988744536980165133407152L), Scalar(0.0555555555555555555555555555555555556L),
		Scalar(0.322751322751322751322751322751322751L),  Scalar(0.0269841269841269841269841269841269841L), Scalar(0.498052448072795834541010978195653994L),
		Scalar{},                                         Scalar(0.498052448072795834541010978195653994L),  Scalar(0.0269841269841269841269841269841269841L),
		Scalar(0.322751322751322751322751322751322751L),  Scalar(0.0555555555555555555555555555555555556L), Scalar(0.0966565466361988744536980165133407152L),
		Scalar{}};
		
	static constexpr std::array<Scalar, 13> s_xi = {
		Scalar(1.00000000000000000000000000000000000L),   Scalar(0.965925826289068286749743199728897368L),  Scalar(0.866025403784438646763723170752936183L),
		Scalar(0.707106781186547524400844362104849039L),  Scalar(0.500000000000000000000000000000000000L),  Scalar(0.258819045102520762348898837624048328L),
		Scalar{},                                         Scalar(-0.258819045102520762348898837624048328L), Scalar(-0.500000000000000000000000000000000000L),
		Scalar(-0.707106781186547524400844362104849039L), Scalar(-0.866025403784438646763723170752936183L), Scalar(-0.965925826289068286749743199728897368L),
		Scalar(-1.00000000000000000000000000000000000L)};
};

} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/ClenshawCurtisHybridAdaptiveQuadrature_impl.hpp>

#endif // LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_HPP
