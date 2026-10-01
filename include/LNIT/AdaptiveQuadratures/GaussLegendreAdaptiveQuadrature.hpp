#ifndef LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_HPP
#define LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_HPP

#include <array>

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase.hpp>
#include <LNIT/misc/Numeric.hpp>

namespace LNIT
{

template<typename T, typename TT> class GaussLegendreAdaptiveQuadrature;

template<typename T, typename TT> 
struct AdaptiveQuadratureTraits< GaussLegendreAdaptiveQuadrature<T, TT> >
{
	using Size       = unsigned int;
	using Scalar     = T;
	using LongScalar = TT;
};

template<typename T, typename TT=T> 
class GaussLegendreAdaptiveQuadrature : public AdaptiveQuadratureBase< GaussLegendreAdaptiveQuadrature<T, TT> >
{
	using Base = AdaptiveQuadratureBase< GaussLegendreAdaptiveQuadrature<T,TT> >;
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
	
	inline constexpr Scalar getMaxDeltaXImpl(const Scalar xmin, const Scalar xmax) const { return (xmax - xmin)*misc::maxDiff(std::span{s_xi}); } 
private:
	std::array<LongScalar, 15> m_fx15;
	
	static constexpr std::array<Scalar, 15> s_xi = {
		Scalar(-0.987992518020485428489565718586612581L), Scalar(-0.937273392400705904307758947710209471L), Scalar(-0.848206583410427216200648320774216851L),
		Scalar(-0.724417731360170047416186054613938010L), Scalar(-0.570972172608538847537226737253910641L), Scalar(-0.394151347077563369897207370981045468L),
		Scalar(-0.201194093997434522300628303394596208L), Scalar{},                                         Scalar(0.201194093997434522300628303394596208L),
		Scalar(0.394151347077563369897207370981045468L),  Scalar(0.570972172608538847537226737253910641L),  Scalar(0.724417731360170047416186054613938010L),
		Scalar(0.848206583410427216200648320774216851L),  Scalar(0.937273392400705904307758947710209471L),  Scalar(0.987992518020485428489565718586612581L)};
	
	static constexpr std::array<Scalar, 15> s_wi15 = {
		Scalar(0.0307532419961172683546283935772044177L),  Scalar(0.0703660474881081247092674164506673385L),  Scalar(0.107159220467171935011869546685869303L),
		Scalar(0.139570677926154314447804794511028323L),   Scalar(0.166269205816993933553200860481208811L),   Scalar(0.186161000015562211026800561866422825L),
		Scalar(0.198431485327111576456118326443839325L),   Scalar(0.202578241925561272880620199967519315L),   Scalar(0.198431485327111576456118326443839325L),
		Scalar(0.186161000015562211026800561866422825L),   Scalar(0.166269205816993933553200860481208811L),   Scalar(0.139570677926154314447804794511028323L),
		Scalar(0.107159220467171935011869546685869303L),   Scalar(0.0703660474881081247092674164506673385L),  Scalar(0.0307532419961172683546283935772044177L)};
		
	static constexpr std::array<Scalar, 15> s_wi14 = {
		Scalar(0.0429480564346795140110031852849063425L),  Scalar(0.0287463102008375282036713620303994957L),  Scalar(0.185198436210474185217363929710859293L),
		Scalar(0.0236554834186314199653732569241027923L),  Scalar(0.316940072793589169562543091505843745L),   Scalar(0.00768583788397500322224645308828128539L),
		Scalar(0.394825803057813179817798721455607046L),   Scalar{},                                          Scalar(0.394825803057813179817798721455607046L),
		Scalar(0.00768583788397500322224645308828128539L), Scalar(0.316940072793589169562543091505843745L),   Scalar(0.0236554834186314199653732569241027923L),
		Scalar(0.185198436210474185217363929710859293L),   Scalar(0.0287463102008375282036713620303994957L),  Scalar(0.0429480564346795140110031852849063425L)};
		
	static constexpr std::array<Scalar, 15> s_wi06 = {
		Scalar{},                                         Scalar(0.214315218972431542643713155033582112L),  Scalar{},
		Scalar(0.0622618038596276360658817591999255133L), Scalar{},                                         Scalar(0.723422977167940821290405085766492375L),
		Scalar{},                                         Scalar{},                                         Scalar{},
		Scalar(0.723422977167940821290405085766492375L),  Scalar{},                                         Scalar(0.0622618038596276360658817591999255133L),
		Scalar{},                                         Scalar(0.214315218972431542643713155033582112L),  Scalar{}};	
};

} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/GaussLegendreAdaptiveQuadrature_impl.hpp>

#endif // LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_HPP
