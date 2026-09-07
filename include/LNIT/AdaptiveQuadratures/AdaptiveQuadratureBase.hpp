#ifndef LNIT_ADAPTIVE_QUADRATURE_BASE_HPP
#define LNIT_ADAPTIVE_QUADRATURE_BASE_HPP

#include <LNIT/NumTraits.hpp>

#include <cstdio>
#include <array>
#include <vector>
#include <span>
#include <numeric>
#include <type_traits>

namespace LNIT
{

template<class Derived> struct AdaptiveQuadratureTraits;

/**
 * @class AdaptiveQuadratureBase
 * @brief Base class for all addaptive quadratures classes (implements the Curiously Recurring Template Pattern (CRTP))
 * @tparam Derived The derived quadrature implementation.
 * 
 * This class provides a generic framework for adaptive quadrature methods, a class of methods 
 * that tries to addapt the quadrature point to the integrated function. 
 * More precisely, an addaptive quadrature method tries to find an subdivision of the integration domain, 
 * \f$[a,b] = \cup_{i=1} [x_{i}, x_{i+1}]\f$, such that the estimated integration error on each subdomains
 * is bellow a given threshold. 
 * 
 * If the integration domain is infinite, the methods tries to find an interval \f$[x_\min, x_\max]\f$ such that
 * \f$\int_{-\infty}^{x_\min} f(x)\,\mathrm{d}x \approx 0\f$ and \f$\int_{x_\max}^{+\infty} f(x)\,\mathrm{d}x \approx 0\f$
 * using the Gauss-Laguerre quadrature to estimate the indefinite integrals. And then uses the addaptive quadrature method
 * to compute \f$\int_{x_\min}^{x_\max} f(x)\,\mathrm{d}x \approx \int_{-\infty}^{+\infty} f(x)\,\mathrm{d}x\f$.
 * 
 * The Derived class must implement the function `estimateIntegralImpl` that returns both the value of the 
 * integral iver a subdomain and the estimated error.
 */
template<class Derived>
class AdaptiveQuadratureBase
{
	using DerivedTraits = AdaptiveQuadratureTraits<Derived>;
public:
	using Size       = typename DerivedTraits::Size;       ///<  @brief Type for iteration counters.
	using Scalar     = typename DerivedTraits::Scalar;     ///<  @brief Floating point type for integration (e.g., double).
	using LongScalar = typename DerivedTraits::LongScalar; ///<  @brief Higher precision type for accumulation (e.g., long double).
	using Interval = std::pair<Scalar, Scalar>;            ///<  @brief Type representing an integration interval.
	
	      Derived& derived()       { return static_cast<      Derived&>(*this); }
	const Derived& derived() const { return static_cast<const Derived&>(*this); }
	
	/**
	 * @brief Construct quadrature object.
	 * @param maxIt Maximum number of iterations.
	 * @param tol Relative tolerance for convergence.
	 */
	AdaptiveQuadratureBase(const Size& maxIt = 20000, const Scalar& relativeTol = NumTraits<Scalar>::epsilon, const Scalar& absoluteTol = NumTraits<Scalar>::epsilon);

	/**
	 * @brief Estimate integral and error on [xmin, xmax].
	 * @tparam Function Callable with signature Scalar f(Scalar).
	 * @param f Function to integrate.
	 * @param xmin Lower bound.
	 * @param xmax Upper bound.
	 * 
	 * Delegates to Derived::estimateIntegralImpl.
	 * 
	 * @return Pair (integral, estimated error).
	 */
	template<class Function> constexpr std::pair<LongScalar, LongScalar> estimateIntegral(const Function& f, const Scalar& xmin, const Scalar& xmax) { return derived().estimateIntegralImpl(f, xmin, xmax); }

	/**
	 * @brief Perform adaptive quadrature on [xmin, xmax].
	 * @tparam Function Callable with signature Scalar f(Scalar).
	 * @param f Function to integrate.
	 * @param xmin Lower bound of interval.
	 * @param xmax Upper bound of interval.
	 * @return Approximation of the integral.
	 *
	 * Bounds follow the usual conventions of the definite integral:
	 * - xmin > xmax returns the opposite of the integral over [xmax, xmin];
	 * - xmin == xmax returns 0 (converged);
	 * - a NaN bound returns NaN (not converged);
	 * - an infinite bound delegates to integrate(f), integrateLeftInfinite() or integrateRightInfinite().
	 */
	template<class Function> LongScalar integrate(const Function& f, const Scalar& xmin, const Scalar& xmax);
	
	/**
	 * @brief Perform adaptive quadrature on (-inf, xmax].
	 * 
	 * First tries to find xmin such that:
	 * \f[
	 * \int_{-\infty}^{xmin} f(x) dx \approx 0.
	 * \f]
	 * And then addapt the quadrature over [xmin, xmax].
	 * 
	 * @tparam Function Callable with signature Scalar f(Scalar).
	 * @param f Function to integrate.
	 * @param xmin Lower bound of interval.
	 * @param xmax Upper bound of interval.
	 * @return Approximation of the integral.
	 */
	template<class Function> LongScalar integrateLeftInfinite(const Function& f, const Scalar& xmax);
	
	/**
	 * @brief Perform adaptive quadrature on [xmin, inf).
	 * 
	 * First tries to find xmax such that:
	 * \f[
	 * \int_{xmax}^{\infty} f(x) dx \approx 0.
	 * \f]
	 * And then addapt the quadrature over [xmin, xmax].
	 * 
	 * @tparam Function Callable with signature Scalar f(Scalar).
	 * @param f Function to integrate.
	 * @param xmin Lower bound of interval.
	 * @param xmax Upper bound of interval.
	 * @return Approximation of the integral.
	 */
	template<class Function> LongScalar integrateRightInfinite(const Function& f, const Scalar& xmin);
	
	/**
     * @brief Perform adaptive quadrature on (-inf, inf).
     * 
     * First tries to find xmin such that:
	 * \f[
	 * \int_{-\infty}^{xmin} f(x) dx \approx 0.
	 * \f]
	 * Then tries to find xmax such that:
	 * \f[
	 * \int_{xmax}^{\infty} f(x) dx \approx 0.
	 * \f]
	 * Finally addapt the quadrature over [xmin, xmax].
     */
	template<class Function> LongScalar integrate(const Function& f);
	
	/**
	 * @brief Perform adaptive quadrature on (-inf, inf) using a coordinate-remapping technique.
	 * 
	 * Uses the following coordinate transformation:
	 * \f[
	 *  x : t \to \frac{t}{1 - t^2}.
	 * \f]
	 * 
	 * So that the integral can be computed as:
	 * \f[
	 * \int_{-1}^{1} f(x(t))\frac{1 + t^2}{(1 - t^2)^2} dt.
	 * \f]
	 */
	template<class Function> LongScalar remapAndIntegrate(const Function& f);

	template<class Function> std::invoke_result_t<Function, Scalar> integrateWithoutAdaptation(const Function& f) const;
	
	constexpr Size   getMaxIt()       const { return m_maxIt; }       ///<  @brief Maximum iterations allowed.
	constexpr Size   getNits()        const { return m_it; }          ///<  @brief Number of iterations performed.
	constexpr Scalar getRelativeTol() const { return m_relativeTol; } ///<  @brief Relative tolerance of the quadrature.
	constexpr Scalar getAbsoluteTol() const { return m_absoluteTol; } ///<  @brief Absolute tolerance of the quadrature.
	
	constexpr bool hasConverged() const { return m_hasConverged; } ///<  @brief Whether convergence was achieved.
	
	constexpr void setMaxIt(const Size& maxIt)       { m_maxIt = maxIt; }                          ///<  @brief Set maximum number of iterations.
	constexpr void setTol(const Scalar& tol)         { m_absoluteTol = tol; m_relativeTol = tol; } ///<  @brief Set the tolerance of the quadrature.
	constexpr void setRelativeTol(const Scalar& tol) { m_relativeTol = tol; }                      ///<  @brief Set the relative tolerance of the quadrature.
	constexpr void setAbsoluteTol(const Scalar& tol) { m_absoluteTol = tol; }                      ///<  @brief Set the absolute tolerance of the quadrature.
	
	constexpr void setMaxIt(Size&& maxIt)       { m_maxIt = std::move(maxIt); }     ///<  @brief Set maximum number of iterations.
	constexpr void setRelativeTol(Scalar&& tol) { m_relativeTol = std::move(tol); } ///<  @brief Set the relative tolerance of the quadrature.
	constexpr void setAbsoluteTol(Scalar&& tol) { m_absoluteTol = std::move(tol); } ///<  @brief Set the absolute tolerance of the quadrature.
	
	constexpr void setOutput(std::FILE* out) { m_out = out; } ///<  @brief Redirect output to a file for logging progress.

	constexpr Scalar getMaxDeltaX(const Scalar& xmin, const Scalar& xmax) const { return derived().getMaxDeltaXImpl(xmin, xmax); } 
	
	constexpr LongScalar getEstimatedIntegral() const { return std::reduce(m_subIntergrals.begin(), m_subIntergrals.end()); }
	constexpr LongScalar getEstimatedError()    const { return std::reduce(m_subIntergralsErr.begin(), m_subIntergralsErr.end()); }
	
	constexpr std::span<const Interval> getSubIntervals() const { return m_intervals; }
private:
	/// Forgets the previous integration: no interval, no iteration, not converged.
	constexpr void resetState();

	std::vector<Interval>   m_intervals;
	std::vector<LongScalar> m_subIntergrals;
	std::vector<LongScalar> m_subIntergralsErr;

	Size   m_maxIt;
	Size   m_it           = 0;
	Scalar m_relativeTol;
	Scalar m_absoluteTol;
	bool   m_hasConverged = false;
	
	std::FILE* m_out = nullptr;
};

template<class T> inline constexpr bool isAdaptiveQuadrature = std::derived_from<T, AdaptiveQuadratureBase<T>>;

template<class T> concept CAdaptiveQuadrature = isAdaptiveQuadrature<T>;

template<CAdaptiveQuadrature Quadrature> using SizeFor       = typename Quadrature::Size;
template<CAdaptiveQuadrature Quadrature> using ScalarFor     = typename Quadrature::Scalar;
template<CAdaptiveQuadrature Quadrature> using LongScalarFor = typename Quadrature::LongScalar;

} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase_impl.hpp>

#endif // LNIT_ADAPTIVE_QUADRATURE_BASE_HPP
