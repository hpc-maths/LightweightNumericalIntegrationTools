#ifndef LNIT_MISC_STRIDED_VIEW_HPP
#define LNIT_MISC_STRIDED_VIEW_HPP

#include <ranges>

namespace LNIT::misc
{

template<std::ptrdiff_t stride, std::ranges::view View>
class StridedView : public std::ranges::view_interface<StridedView<stride, View>>
{
	static_assert(std::ranges::input_range<View>);
public:
	using Sentinel = std::default_sentinel_t;

	template<bool isConst>
	class Iterator
	{
		template<typename T>
		using MaybeConst = std::conditional_t<isConst, const T, T>;

		using InnerIterator = std::ranges::iterator_t<MaybeConst<View>>;
		using InnerSentinel = std::ranges::sentinel_t<MaybeConst<View>>;

		using       Reference = std::iter_reference_t<InnerIterator>;
		using const_Reference = std::common_reference_t<
			const std::iter_value_t<InnerIterator>&&, 
			std::iter_reference_t<InnerIterator>>;
		
	public:	
		using difference_type   = std::iter_difference_t<InnerIterator>;
		using value_type        = std::iter_value_t<InnerIterator>;
		using reference         = std::conditional_t<isConst, const_Reference, Reference>;
		using iterator_category = std::input_iterator_tag;	

		Iterator(const InnerIterator first, const InnerSentinel bound) : m_current(first), m_bound(bound) {}

		constexpr reference operator*() const { return *m_current; }

		constexpr Iterator& operator++() { next(); return *this; }
		
		constexpr Iterator operator++(int) { Iterator ret(*this); ++(*this); return ret; }

		friend constexpr bool operator==(const Iterator& it, const Sentinel /* sentinel */) { return it.isDone(); }
	private:
		constexpr bool isDone() const;

		constexpr void next();

		InnerIterator m_current;
		InnerSentinel m_bound;
	};

	constexpr StridedView(View&& view) : m_view(std::forward<View>(view)), m_bound(stride*(std::ranges::distance(m_view) / stride)) {}

	constexpr Iterator<false> begin() { return Iterator<false>(std::ranges::begin(m_view), std::ranges::end(m_view)); }

	constexpr Sentinel end() { return {}; }

	constexpr Iterator<true> begin() const { return Iterator<true>(std::ranges::begin(m_view), std::ranges::end(m_view)); }

	constexpr Sentinel end() const { return {}; }
private:
	View m_view;
	std::ptrdiff_t m_bound;
};

template<std::ptrdiff_t Tstride, std::ranges::input_range Range>
constexpr StridedView<Tstride, std::views::all_t<Range>> stride(Range&& range) { return StridedView<Tstride, std::views::all_t<Range>>(std::views::all(range)); }
	
} // namespace LNIT::misc

template<std::ptrdiff_t Tstride, std::ranges::view View>
constexpr bool std::ranges::enable_borrowed_range< LNIT::misc::StridedView<Tstride, View> > = std::ranges::enable_borrowed_range< View >;

#include <LNIT/misc/StridedView_impl.hpp>

#endif // LNIT_MISC_STRIDED_VIEW_HPP
