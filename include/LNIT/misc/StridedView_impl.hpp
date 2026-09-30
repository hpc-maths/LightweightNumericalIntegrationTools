#ifndef LNIT_MISC_STRIDED_VIEW_IMPL_HPP
#define LNIT_MISC_STRIDED_VIEW_IMPL_HPP

#include <algorithm>

namespace LNIT::misc
{

template<std::ptrdiff_t stride, std::ranges::view View> template<bool isConst>
constexpr bool StridedView<stride,View>::Iterator<isConst>::isDone() const
{
	if constexpr (std::random_access_iterator<InnerIterator>)
	{
		return m_current >= m_bound;
	}
	else
	{
		return m_current == m_bound;
	}
} 

template<std::ptrdiff_t stride, std::ranges::view View> template<bool isConst>
constexpr void StridedView<stride,View>::Iterator<isConst>::next()
{
	if constexpr (std::random_access_iterator<InnerIterator>)
	{
		m_current += std::min<difference_type>(stride, m_bound - m_current);
	}
	else
	{
		std::ranges::advance(m_current, stride, m_bound);
	}
}

} // namespace LNIT::misc

#endif // LNIT_MISC_STRIDED_VIEW_IMPL_HPP
