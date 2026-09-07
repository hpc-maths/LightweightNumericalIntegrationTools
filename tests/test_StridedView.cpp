// Test of LNIT::misc::StridedView: iterator/sentinel model, C++20 range concepts
// and iteration results over several kinds of underlying ranges.
//

#include "Check.hpp"

#include <LNIT/misc/StridedView.hpp>

#include <array>
#include <cstdio>
#include <cstdlib>
#include <forward_list>
#include <list>
#include <sstream>
#include <vector>

using namespace LNIT::tests;

namespace
{

using LNIT::misc::stride;

template<class V>
void checkConcepts()
{
	using Iterator      = decltype(std::declval<V&>().begin());
	using ConstIterator = decltype(std::declval<const V&>().begin());
	using Sentinel      = typename V::Sentinel;

	static_assert(std::input_iterator<Iterator>);
	static_assert(std::input_iterator<ConstIterator>);
	static_assert(std::sentinel_for<Sentinel, Iterator>);
	static_assert(std::sentinel_for<Sentinel, ConstIterator>);
	static_assert(std::ranges::input_range<V>);
	static_assert(std::ranges::input_range<const V>);
	static_assert(std::ranges::view<V>);
}

template<class Range>
std::vector<double> collect(Range&& range)
{
	std::vector<double> out;
	for (auto&& x : range) { out.push_back(static_cast<double>(x)); }
	return out;
}

using Values = std::vector<double>;

void testRandomAccess()
{
	Values v{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

	auto sv = stride<2>(v);
	checkConcepts<decltype(sv)>();
	CHECK((collect(sv) == Values{0, 2, 4, 6, 8}));

	const auto& csv = sv;
	CHECK((collect(csv) == Values{0, 2, 4, 6, 8}));

	// length not a multiple of the stride: the last step must stop at end()
	CHECK((collect(stride<3>(v))  == Values{0, 3, 6, 9}));
	CHECK((collect(stride<4>(v))  == Values{0, 4, 8}));
	CHECK((collect(stride<7>(v))  == Values{0, 7}));
	CHECK((collect(stride<10>(v)) == Values{0}));
	CHECK((collect(stride<11>(v)) == Values{0}));
	CHECK((collect(stride<1>(v))  == v));

	Values empty;
	CHECK(collect(stride<2>(empty)).empty());

	// std::array from a const source
	const std::array<int, 5> a{1, 2, 3, 4, 5};
	CHECK((collect(stride<2>(a)) == Values{1, 3, 5}));

	// the non-const iterator must yield mutable references
	for (auto& x : stride<2>(v)) { x = -x; }
	CHECK((v == Values{0, 1, -2, 3, -4, 5, -6, 7, -8, 9}));
}

void testBidirectional()
{
	std::list<int> l{0, 1, 2, 3, 4, 5, 6};

	auto sl = stride<3>(l);
	checkConcepts<decltype(sl)>();
	CHECK((collect(sl) == Values{0, 3, 6}));
	CHECK((collect(stride<2>(l)) == Values{0, 2, 4, 6}));
	CHECK((collect(stride<4>(l)) == Values{0, 4}));
}

void testForwardOnly()
{
	std::forward_list<int> f{0, 1, 2, 3, 4};

	auto sf = stride<2>(f);
	checkConcepts<decltype(sf)>();
	CHECK((collect(sf) == Values{0, 2, 4}));
	CHECK((collect(stride<3>(f)) == Values{0, 3}));
}

void testSinglePassInput()
{
	// a pure input range: the view must not consume it before iteration
	// (e.g. by calling std::ranges::distance in its constructor)
	std::istringstream in("0 1 2 3 4 5 6");
	auto si = stride<3>(std::views::istream<int>(in));
	static_assert(std::ranges::input_range<decltype(si)>);
	static_assert(std::ranges::view<decltype(si)>);
	CHECK((collect(si) == Values{0, 3, 6}));

	std::istringstream in2("0 1 2 3 4");
	CHECK((collect(stride<2>(std::views::istream<int>(in2))) == Values{0, 2, 4}));
}

void testIteratorAndSentinel()
{
	Values v{0, 1, 2, 3};
	auto sv = stride<2>(v);

	auto it  = sv.begin();
	auto old = it++;
	CHECK(*old == 0);
	CHECK(*it == 2);

	// comparison with the sentinel, both operand orders, == and !=
	CHECK(!(it == std::default_sentinel));
	CHECK(it != std::default_sentinel);
	CHECK(!(std::default_sentinel == it));
	CHECK(std::default_sentinel != it);

	++it;
	CHECK(it == std::default_sentinel);
	CHECK(std::default_sentinel == it);

	static_assert(std::ranges::borrowed_range<decltype(stride<2>(v))>);
}

} // namespace

int main()
{
	testRandomAccess();
	testBidirectional();
	testForwardOnly();
	testSinglePassInput();
	testIteratorAndSentinel();

	return report("test_StridedView");
}
