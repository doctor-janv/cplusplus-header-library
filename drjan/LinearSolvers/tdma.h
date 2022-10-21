#ifndef DRJANV_TDMA_H
#define DRJANV_TDMA_H

#include <vector>
#include <stddef.h>

namespace drjanv
{
/**Tri-diagonal matrix algorithm.
 * 
 \param A STL-vector of doubles. The lower diagonal. a_0 lies on the first row 
        but is not used.
 \param B STL-vector of doubles. The main diagonal. Size N.
 \param C STL-vector of doubles. The upper diagonal. c_0 lies on the first row but
        c_(N-1) lies on the last row and is not used.
 \param D STL-vector of doubles. The right-hand-side. Size N.
 return std::vector<double> The solution.
 */
inline
std::vector<double> TDMA(std::vector<double>& A,
                         std::vector<double>& B,
                         std::vector<double>& C,
                         std::vector<double>& D)
{
    const size_t N = D.size();
    std::vector<double> X(N,0.0);

    for (int i=1; i<N; ++i)
    {
        double W = A[i] / B[i-1];
        B[i] = B[i] - W*C[i-1];
        D[i] = D[i] - W*D[i-1];
    }

    X[N-1] = D[N-1] / B[N-1];
    for (int i=(N-2); i>=0; --i)
    {
        X[i] = (D[i] - C[i]*X[i+1])/B[i];
    }

    return X;
}

}//namespace drjanv

#endif //DRJANV_TDMA_H