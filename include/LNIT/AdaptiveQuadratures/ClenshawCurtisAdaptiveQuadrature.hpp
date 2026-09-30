#ifndef LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_HPP
#define LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_HPP

#include <cmath>

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase.hpp>
#include <LNIT/misc/Numeric.hpp>

namespace LNIT
{

template<typename T, typename TT> class ClenshawCurtisAdaptiveQuadrature;

template<typename T, typename TT> 
struct AdaptiveQuadratureTraits< ClenshawCurtisAdaptiveQuadrature<T, TT> >
{
	using Size       = unsigned int;
	using Scalar     = T;
	using LongScalar = TT;
};

/**
 * @brief Adaptive Clenshaw–Curtis quadrature for numerical integration.
 *
 * This class implements adaptive quadrature using Clenshaw–Curtis rules
 * with nested evaluation points (9, 17, and 33 point formulas).
 *
 * @tparam T Floating point type for integration (e.g., double).
 * @tparam TT Higher precision type for accumulation (e.g., long double).
 */
template<typename T, typename TT=T>
class ClenshawCurtisAdaptiveQuadrature : public AdaptiveQuadratureBase< ClenshawCurtisAdaptiveQuadrature<T,TT> >
{
	using Base = AdaptiveQuadratureBase< ClenshawCurtisAdaptiveQuadrature<T,TT> >;
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
	std::array<LongScalar, 33> m_fx33;

	static constexpr std::array<Scalar, 9> s_wi09 = {
		Scalar(0.015873015873015873016), Scalar(0.14621864921601815501 ), Scalar(0.27936507936507936508 ), Scalar(0.36171785872048978150 ), Scalar(0.39365079365079365079 ), 
		Scalar(0.36171785872048978150 ), Scalar(0.27936507936507936508 ), Scalar(0.14621864921601815501 ), Scalar(0.015873015873015873016)};
	
	static constexpr std::array<Scalar, 17> s_wi17 = {
		Scalar(0.0039215686274509803922), Scalar(0.037368702837205610321 ), Scalar(0.075482331543151834413 ), Scalar(0.10890555258189093044  ), Scalar(0.13895646836823307412  ), 
		Scalar(0.16317266428170330256  ), Scalar(0.18147378423649335700  ), Scalar(0.19251386461292564687  ), Scalar(0.19641012582189052777  ), Scalar(0.19251386461292564687  ), 
		Scalar(0.18147378423649335700  ), Scalar(0.16317266428170330256  ), Scalar(0.13895646836823307412  ), Scalar(0.10890555258189093044  ), Scalar(0.075482331543151834413 ), 
		Scalar(0.037368702837205610321 ), Scalar(0.0039215686274509803922)};
	
	static constexpr std::array<Scalar, 33> s_wi33 = {
		Scalar(0.00097751710654936461388), Scalar(0.0093931979629550147012 ), Scalar(0.019234245132681149183  ), Scalar(0.028457916677233690094  ), Scalar(0.037594341914047206016  ), 
		Scalar(0.046262762837751749492  ), Scalar(0.054555016303980310438  ), Scalar(0.062272109545294004553  ), Scalar(0.069427575630435450900  ), Scalar(0.075883800441388470480  ), 
		Scalar(0.081634817654938510229  ), Scalar(0.086577538441827435439  ), Scalar(0.090706112867720998737  ), Scalar(0.093943244438768735729  ), Scalar(0.096292325945488179193  ), 
		Scalar(0.097698188208055581820  ), Scalar(0.098178577781768296767  ), Scalar(0.097698188208055581820  ), Scalar(0.096292325945488179193  ), Scalar(0.093943244438768735729  ), 
		Scalar(0.090706112867720998737  ), Scalar(0.086577538441827435439  ), Scalar(0.081634817654938510229  ), Scalar(0.075883800441388470480  ), Scalar(0.069427575630435450900  ), 
		Scalar(0.062272109545294004553  ), Scalar(0.054555016303980310438  ), Scalar(0.046262762837751749492  ), Scalar(0.037594341914047206016  ), Scalar(0.028457916677233690094  ), 
		Scalar(0.019234245132681149183  ), Scalar(0.0093931979629550147012 ), Scalar(0.00097751710654936461388)};
		
	static constexpr std::array<Scalar, 33> s_xi = {
		 Scalar( 1.0000000000000000000  ), Scalar( 0.99518472667219688624), Scalar( 0.98078528040323044913 ), Scalar( 0.95694033573220886494), Scalar( 0.92387953251128675613), 
		 Scalar( 0.88192126434835502971 ), Scalar( 0.83146961230254523708), Scalar( 0.77301045336273696081 ), Scalar( 0.70710678118654752440), Scalar( 0.63439328416364549822), 
		 Scalar( 0.55557023301960222474 ), Scalar( 0.47139673682599764856), Scalar( 0.38268343236508977173 ), Scalar( 0.29028467725446236764), Scalar( 0.19509032201612826785), 
		 Scalar( 0.098017140329560601994), Scalar(                       ), Scalar(-0.098017140329560601994), Scalar(-0.19509032201612826785), Scalar(-0.29028467725446236764), 
		 Scalar(-0.38268343236508977173 ), Scalar(-0.47139673682599764856), Scalar(-0.55557023301960222474 ), Scalar(-0.63439328416364549822), Scalar(-0.70710678118654752440),
		 Scalar(-0.77301045336273696081 ), Scalar(-0.83146961230254523708), Scalar(-0.88192126434835502971 ), Scalar(-0.92387953251128675613), Scalar(-0.95694033573220886494), 
		 Scalar(-0.98078528040323044913 ), Scalar(-0.99518472667219688624), Scalar(-1.0000000000000000000  )};
};
	
} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/ClenshawCurtisAdaptiveQuadrature_impl.hpp>

#endif // LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_HPP
