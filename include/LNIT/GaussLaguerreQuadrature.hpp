#ifndef LNIT_GAUSS_LAGUERRE_QUADRATURE_HPP
#define LNIT_GAUSS_LAGUERRE_QUADRATURE_HPP

#include <array>

namespace LNIT
{

/**
 * @class GaussLaguerreQuadrature
 * @brief Template class implementing Gauss-Laguerre quadrature.
 *
 * @tparam Scalar Numeric type for computation (e.g., float, double).
 * @tparam LongScalar Optional extended precision type for accumulation (defaults to Scalar).
 *
 * This class provides methods to numerically approximate integrals of functions defined
 * on semi-infinite intervals using Gauss-Laguerre quadrature rules.
 */
template<typename Scalar, typename LongScalar=Scalar>
class GaussLaguerreQuadrature
{
public:
	using Size = unsigned int; ///< Type alias for sizes/indices.

	constexpr GaussLaguerreQuadrature() {} ///<  @brief Default constructor.

	/**
	 * @brief Approximate integral over the left semi-infinite interval.
	 *
	 * Computes an approximation of the integral:
	 * \f[
	 * \int_{-\infty}^{a} f(x) dx
	 * \f]
	 * using Gauss-Laguerre quadrature.
	 *
	 * @tparam Function Type of the callable object (e.g., lambda, functor).
	 * @param f Integrand function \f$ f(x) \f$.
	 * @param a Upper bound of integration (default = 0).
	 * @return Approximated integral value.
	 */
	template<class Function> constexpr LongScalar integrateLeftInfinite  (const Function& f, const Scalar& a = Scalar{}) const;
	/**
	 * @brief Approximate integral over the right semi-infinite interval.
	 *
	 * Computes an approximation of the integral:
	 * \f[
	 * \int_{a}^{+\infty} f(x) dx
	 * \f]
	 * using Gauss-Laguerre quadrature.
	 *
	 * @tparam Function Type of the callable object (e.g., lambda, functor).
	 * @param f Integrand function \f$ f(x) \f$.
	 * @param a Lower bound of integration (default = 0).
	 * @return Approximated integral value.
	 */
	template<class Function> constexpr LongScalar integrateRightInfinite (const Function& f, const Scalar& a = Scalar{}) const;
private:
	static constexpr std::array<Scalar, 33> s_wi = {
		Scalar(0.110777305873207582811369398847383061L), Scalar(0.258105281281894765775426191201727153L), Scalar(0.406221768684373689182722487278337766L),
		Scalar(0.555262309599223033614149156625912369L), Scalar(0.705557387659582894269149721469723421L), Scalar(0.857465745197485609229412709557854856L),
		Scalar(1.01136895238079744580316305006241479L),  Scalar(1.16767566070191120284646309660091419L),  Scalar(1.32682791710863578123173794592947218L),
		Scalar(1.48930888165730039036528338741349956L),  Scalar(1.65565213521125686550559232766859788L),  Scalar(1.82645300470519336825076148878133150L),
		Scalar(2.00238251829042438705175117055488475L),  Scalar(2.18420483134563109633037959094957254L),  Scalar(2.37279928409113989563411721117838668L),
		Scalar(2.56918871553133835523706366327238758L),  Scalar(2.77457634678431634819664960331107515L),  Scalar(2.99039458946476195387029818967807993L),
		Scalar(3.21837074971698879981007017130934926L),  Scalar(3.46061716286773553012855282842708437L),  Scalar(3.71975748100558828635181744528749409L),
		Scalar(3.99910789299344150282858719087152052L),  Scalar(4.30294438794235735170168286496576500L),  Scalar(4.63690963770119300191804212610809531L),
		Scalar(5.00865603705641702982530424948045510L),  Scalar(5.42890847474740860685712622963535001L),  Scalar(5.91331937782164966954548794814728759L),
		Scalar(6.48593541397783535743845885932661548L),  Scalar(7.18627247285383848215059528921236267L),  Scalar(8.08557165264589566710884825274056954L),
		Scalar(9.33105925146152805395606955428187483L),  Scalar(11.3036226565875150444131601076159571L),  Scalar(15.5662914941639307905211633161087023L)};
	static constexpr std::array<Scalar, 33> s_xi = {
		Scalar(0.0431611356173268921917334738205544067L), Scalar(0.227517802803371123850290226912660576L),  Scalar(0.559616655851539887586282303915739251L),
		Scalar(1.04026850775100205382209621926814908L),   Scalar(1.67055919607571519092562973256704091L),   Scalar(2.45192079589763054651073898192462989L),
		Scalar(3.38615533758800483230187851831538199L),   Scalar(4.47545949839977145702059137904646435L),   Scalar(5.72245472027210352266790817932965387L),
		Scalar(7.13022434440010801631414039533509932L),   Scalar(8.70235923062140624893696399458808474L),   Scalar(10.4430136502059824268455293838802160L),
		Scalar(12.3569737593502859624441255236481692L),   Scalar(14.4497416815855402377145121177738701L),   Scalar(16.7276392186383223532615229941783570L),
		Scalar(19.1979365872124466372283088222023836L),   Scalar(21.8690135249281898713512287042778311L),   Scalar(24.7505629061577956433730931987332932L),
		Scalar(27.8538511114133567797747375537319113L),   Scalar(31.1920555455751298677734295988879550L),   Scalar(34.7807091535383377002292521852636020L),
		Scalar(38.6382967177740302250360622751110273L),   Scalar(42.7870720782534794879639219926530226L),   Scalar(47.2542066029932658172690829766684347L),
		Scalar(52.0734519015142202671640200482332471L),   Scalar(57.2876345410929400754514841077799217L),   Scalar(62.9525659469066302071906336861063195L),
		Scalar(69.1435133801098924457366348146804507L),   Scalar(75.9666870142470623437939790249571279L),   Scalar(83.5816372232708807614192336050103633L),
		Scalar(92.2511394441351012341481184391346172L),   Scalar(102.477844336823322575825984749546501L),   Scalar(115.554756448995807306876850792517920L)};
};

} // namespace LNIT

#include <LNIT/GaussLaguerreQuadrature_impl.hpp>

#endif // LNIT_GAUSS_LAGUERRE_QUADRATURE_HPP
