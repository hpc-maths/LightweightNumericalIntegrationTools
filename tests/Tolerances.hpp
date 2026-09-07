#ifndef LNIT_TESTS_TOLERANCES_HPP
#define LNIT_TESTS_TOLERANCES_HPP

// Named tolerance tiers used by every test (see the LNIT test plan, section 1).
// Never write a magic tolerance literal in a test case: pick a tier here and,
// when a looser tier is needed, justify it with a comment at the call site.

namespace LNIT::tests::tol
{

/// Polynomial exactness of a rule, degree <= 20 on [-1, 1].
inline constexpr double EXACT    = 2.e-14;
/// Polynomial exactness, degrees 21 to 33 (conditioning of the weighted sums).
inline constexpr double EXACT_HI = 1.e-12;
/// Comparison with closed forms or high-precision references.
inline constexpr double REF      = 1.e-13;
/// Adaptive contract "the true error respects the requested tolerance": requested tolerance times this factor.
inline constexpr double TOL_USER_FACTOR = 10.;

} // namespace LNIT::tests::tol

#endif // LNIT_TESTS_TOLERANCES_HPP
