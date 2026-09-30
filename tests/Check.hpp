#ifndef LNIT_TESTS_CHECK_HPP
#define LNIT_TESTS_CHECK_HPP

// Minimal framework-free checks shared by the unit tests.
// CHECK() is not compiled out in Release builds, unlike assert().

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace LNIT::tests
{

inline int g_failures = 0;

/// Mixed absolute/relative comparison: |a - b| <= atol + rtol * max(|a|, |b|).
template<class T>
constexpr bool closeTo(const T a, const T b, const T rtol, const T atol = T{})
{
	using std::abs;
	using std::max;
	return abs(a - b) <= atol + rtol*max(abs(a), abs(b));
}

inline int report(const char* testName)
{
	if (g_failures != 0)
	{
		std::fprintf(stderr, "%s: %d check(s) failed\n", testName, g_failures);
		return EXIT_FAILURE;
	}
	std::printf("%s: all checks passed\n", testName);
	return EXIT_SUCCESS;
}

} // namespace LNIT::tests

#define CHECK(cond)                                                                            \
	do                                                                                         \
	{                                                                                          \
		if (!(cond))                                                                           \
		{                                                                                      \
			std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond);      \
			++LNIT::tests::g_failures;                                                         \
		}                                                                                      \
	} while (false)

#define CHECK_CLOSE(a, b, rtol, atol)                                                                             \
	do                                                                                                            \
	{                                                                                                             \
		const auto lnit_a = (a);                                                                                  \
		const auto lnit_b = (b);                                                                                  \
		if (!LNIT::tests::closeTo(lnit_a, lnit_b, rtol, atol))                                                    \
		{                                                                                                         \
			std::fprintf(stderr, "%s:%d: CHECK_CLOSE failed: %s = %.17g, %s = %.17g, |diff| = %.3g\n",             \
			             __FILE__, __LINE__, #a, double(lnit_a), #b, double(lnit_b), double(std::abs(lnit_a - lnit_b))); \
			++LNIT::tests::g_failures;                                                                            \
		}                                                                                                         \
	} while (false)

#endif // LNIT_TESTS_CHECK_HPP
