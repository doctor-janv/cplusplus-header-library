#ifndef DRJANV_UTILITIES_RANGE_H
#define DRJANV_UTILITIES_RANGE_H

#include <type_traits>
#include <vector>

namespace drjanv
{
/**Returns a range of number according to the logic of the parameters.
 *
 * \param start First number in the sequence.
 * \param end   Termination criteria. If the delta is positive then
 *              the sequence will terminate if i>=end, otherwise if the
 *              delta is negative the sequence will terminate if i<=end
 * \param delta Cannot be 0. Default 1. Can be negative.
 * 
 * \return A `std::vector` of template type `T` containing the range according 
 *         to the logic.
 * 
 * Example:
 * \code {.cpp}
 * auto iorder = drjanv::Range<int>(0, 10); // 0,1,...,9
 * auto iorder = drjanv::Range<int>(9,-1,-1); //9,8,...,0
 * \endcode
 * 
 * */
template<typename T, typename D = int>
std::vector<T> Range(T start, T end, D delta=1)
{
  static_assert(std::is_signed<D>::value,
                "chi_math::Range delta parameter must be signed");
  const bool forward = (delta > 0);

  std::vector<T> sequence = {};

  if (    forward and start >= end) return sequence;
  if (not forward and start <= end) return sequence;

  T i = start;
  bool terminate = false;
  while (not terminate)
  {
    sequence.push_back(i);
    i += delta;

    if (    forward and i >= end) terminate = true;
    if (not forward and i <= end) terminate = true;

    if (not forward and i > start) terminate = true; //Wrap-around check
  }

  return sequence;
}

} // namespace drjanv


#endif //#DRJANV_UTILITIES_RANGE_H