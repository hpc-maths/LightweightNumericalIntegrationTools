#ifndef LNIT_GL_CC_ADAPTIVE_QUADRATURE_HPP
#define LNIT_GL_CC_ADAPTIVE_QUADRATURE_HPP

#include <array>

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase.hpp>
#include <LNIT/misc/Numeric.hpp>

namespace LNIT
{

template<typename T, typename TT> class GLCCAdaptiveQuadrature;

template<typename T, typename TT> 
struct AdaptiveQuadratureTraits< GLCCAdaptiveQuadrature<T, TT> >
{
	using Size       = unsigned int;
	using Scalar     = T;
	using LongScalar = TT;
};

/**
 * @class GLCCAdaptiveQuadrature
 * @brief Adaptive quadrature using combined Gauss–Legendre (GL) and Clenshaw–Curtis (CC) rules.
 *
 * This class computes two high-order estimates of the same integral — a
 * 15-point Gauss–Legendre and a 33-point Clenshaw–Curtis rule — which both are
 * of order 31. Because the two rules share the same order, their
 * difference provides a very accurate estimate of the integration error;
 * that error estimate is used to drive recursive adaptive subdivision.
 *
 * @tparam T  Scalar floating-point type (e.g., double).
 * @tparam TT Extended precision accumulator type (e.g., long double).
 */
template<typename T, typename TT=T> 
class GLCCAdaptiveQuadrature : public AdaptiveQuadratureBase< GLCCAdaptiveQuadrature<T, TT> >
{
	using Base = AdaptiveQuadratureBase< GLCCAdaptiveQuadrature<T,TT> >;
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
	template<class Function> constexpr std::pair<LongScalar, LongScalar> estimateIntegralImpl(const Function& f, const Scalar& xmin, const Scalar& xmax) const;
	
	template<class Function> constexpr std::invoke_result_t<Function, Scalar> integrateImpl(const Function& f, const Scalar& xmin, const Scalar& xmax) const;
	
	constexpr Scalar getMaxDeltaXImpl(const Scalar xmin, const Scalar xmax) const { return (xmax - xmin)*misc::maxDiff(std::span{s_xi_cc}); } 
#ifdef LNIT_TESTING
public: // quadrature tables exposed to the unit tests only
#else
private:
#endif // LNIT_TESTING
	static constexpr std::array<Scalar, 15> s_xi_gl{
		Scalar(-0.98799251802048542849), -Scalar(0.93727339240070590431), -Scalar(0.84820658341042721620), Scalar(-0.72441773136017004742), Scalar(-0.57097217260853884754), 
		Scalar(-0.39415134707756336990), -Scalar(0.20119409399743452230),  Scalar(                      ), Scalar( 0.20119409399743452230), Scalar( 0.39415134707756336990), 
		Scalar( 0.57097217260853884754),  Scalar(0.72441773136017004742),  Scalar(0.84820658341042721620), Scalar( 0.93727339240070590431), Scalar( 0.98799251802048542849)};
	
	static constexpr std::array<Scalar, 15> s_wi_gl{
		Scalar(0.030753241996117268355), Scalar(0.070366047488108124709), Scalar(0.10715922046717193501), Scalar(0.13957067792615431445 ), Scalar(0.16626920581699393355 ), 
		Scalar(0.18616100001556221103 ), Scalar(0.19843148532711157646 ), Scalar(0.20257824192556127288), Scalar(0.19843148532711157646 ), Scalar(0.18616100001556221103 ), 
		Scalar(0.16626920581699393355 ), Scalar(0.13957067792615431445 ), Scalar(0.10715922046717193501), Scalar(0.070366047488108124709), Scalar(0.030753241996117268355)};
	
	static constexpr std::array<Scalar, 33> s_xi_cc = {
		 Scalar( 1.0000000000000000000  ), Scalar( 0.99518472667219688624), Scalar( 0.98078528040323044913 ), Scalar( 0.95694033573220886494), Scalar( 0.92387953251128675613), 
		 Scalar( 0.88192126434835502971 ), Scalar( 0.83146961230254523708), Scalar( 0.77301045336273696081 ), Scalar( 0.70710678118654752440), Scalar( 0.63439328416364549822), 
		 Scalar( 0.55557023301960222474 ), Scalar( 0.47139673682599764856), Scalar( 0.38268343236508977173 ), Scalar( 0.29028467725446236764), Scalar( 0.19509032201612826785), 
		 Scalar( 0.098017140329560601994), Scalar(                       ), Scalar(-0.098017140329560601994), Scalar(-0.19509032201612826785), Scalar(-0.29028467725446236764), 
		 Scalar(-0.38268343236508977173 ), Scalar(-0.47139673682599764856), Scalar(-0.55557023301960222474 ), Scalar(-0.63439328416364549822), Scalar(-0.70710678118654752440),
		 Scalar(-0.77301045336273696081 ), Scalar(-0.83146961230254523708), Scalar(-0.88192126434835502971 ), Scalar(-0.92387953251128675613), Scalar(-0.95694033573220886494), 
		 Scalar(-0.98078528040323044913 ), Scalar(-0.99518472667219688624), Scalar(-1.0000000000000000000  )};
		
	static constexpr std::array<Scalar, 33> s_wi_cc = {
		Scalar(0.00097751710654936461388), Scalar(0.0093931979629550147012 ), Scalar(0.019234245132681149183  ), Scalar(0.028457916677233690094  ), Scalar(0.037594341914047206016  ), 
		Scalar(0.046262762837751749492  ), Scalar(0.054555016303980310438  ), Scalar(0.062272109545294004553  ), Scalar(0.069427575630435450900  ), Scalar(0.075883800441388470480  ), 
		Scalar(0.081634817654938510229  ), Scalar(0.086577538441827435439  ), Scalar(0.090706112867720998737  ), Scalar(0.093943244438768735729  ), Scalar(0.096292325945488179193  ), 
		Scalar(0.097698188208055581820  ), Scalar(0.098178577781768296767  ), Scalar(0.097698188208055581820  ), Scalar(0.096292325945488179193  ), Scalar(0.093943244438768735729  ), 
		Scalar(0.090706112867720998737  ), Scalar(0.086577538441827435439  ), Scalar(0.081634817654938510229  ), Scalar(0.075883800441388470480  ), Scalar(0.069427575630435450900  ), 
		Scalar(0.062272109545294004553  ), Scalar(0.054555016303980310438  ), Scalar(0.046262762837751749492  ), Scalar(0.037594341914047206016  ), Scalar(0.028457916677233690094  ), 
		Scalar(0.019234245132681149183  ), Scalar(0.0093931979629550147012 ), Scalar(0.00097751710654936461388)};

};

} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/GLCCAdaptiveQuadrature_impl.hpp>

#endif // LNIT_GL_CC_ADAPTIVE_QUADRATURE_HPP
