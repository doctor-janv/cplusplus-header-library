#ifndef DRJANV_QUADRATURES_LEGENDRE_H
#define DRJANV_QUADRATURES_LEGENDRE_H

#include <vector>
#include <cmath>
#include <algorithm>

namespace drjanv
{
  /**Provides the function evaluation of the Legendre polynomial
   * \f$P_N\f$ at value x.

  \param N int Order of the Legendre polynomial.
  \param x double The evaluation point.*/
  inline
  double Legendre(int N, double x)
  {
    double Pnm1 = 1;
    double Pn   = x;
    double Pnp1 = 0;

    if (N==0) {return 1;}

    if (N==1) {return x;}

    for (int n=2;n<=N; n++)
    {
      int ns = n-1;
      Pnp1 = ((2.0*ns+1)/(ns+1.0))*x*Pn - (ns/(ns+1.0))*Pnm1;
      Pnm1 = Pn;
      Pn = Pnp1;
    }

    return Pnp1;
  }

  /**Provides the function evaluation of the derivative of the Legendre 
   * polynomial \f$\frac{dP_N}{dx}\f$ at value x.

  \param N int Order of the Legendre polynomial.
  \param x double The evaluation point.*/
  inline
  double dLegendredx(int N, double x)
  {
    if (N==0) {return 0;}

    if (N==1) {return 1;}

    double retval = (N*x/(x*x-1))*Legendre(N,x);
          retval-= (N/(x*x-1))*Legendre(N-1,x);

    return retval;
  }

  /**Provides the function evaluation of the second derivative of the Legendre 
   * polynomial \f$\frac{d^2 P_N}{dx^2}\f$ at value x.

  \param N int Order of the Legendre polynomial.
  \param x double The evaluation point.*/
  inline
  double d2Legendredx2(int N, double x)
  {
    double epsilon = 1.0e-8;
    if (N==0) {return 0.0;}

    if (N==1) {return 0.0;}

    double xpos = std::min(x+epsilon, 1.0-1.0e-10);
    double xneg = std::max(x-epsilon,-1.0+1.0e-10);
    double dx = xpos - xneg;

    double dPdx_pos = dLegendredx(N,xpos);
    double dPdx_neg = dLegendredx(N,xneg);

    return (dPdx_pos - dPdx_neg)/dx;
  }

  /** Finds the roots of the Legendre polynomial.
   *
   * The algorithm is that depicted in:
   *
   * [1] Barrera-Figueroa, et al., "Multiple root finder algorithm for Legendre
   *     and Chebyshev polynomials via Newton's method", Annales Mathematicae et
   *     Informaticae, 33 (2006) pp. 3-13.
   *
   * \param N Is the order of the polynomial.
   * \param roots Is a reference to the roots.
   * \param max_iters Maximum newton iterations to perform for each root.
   *        Default: 1000.
   * \param tol Tolerance at which the newton iteration will be terminated.
   *        Default: 1.0e-12.
   * 
   * \return A `std::vector<double>` containg a sorted list of roots in the 
   *         interval [-1,1].
   *
   * \author Jan*/
  inline
  std::vector<double> LegendreRoots(unsigned int N, 
                                    unsigned int max_iters=1000, 
                                    double tol=1.0e-12)
{
  //======================================== Populate init guess
  //This initial guess proved to be quite important
  //at higher N since the roots start to get
  //squeezed to -1 and 1.
  int num_search_intvls = 1000;
  if (N>64)
    num_search_intvls *= 10;
  if (N>256)
    num_search_intvls *= 10;
  if (N>768)
    num_search_intvls *= 10;

  if (N>2056)
  {
    num_search_intvls *= 10;
    std::cout 
      << "drjan::LegendreRoots: "
      << "The order of the polynomial for which to find the roots is "
      << "greater than 2056. Accuracy of the root finder will be diminished "
      << "along with a reduction in stability.";
  }

  // For this code we simply check to see where the
  // polynomial changes sign.
  double delta = 2.0/num_search_intvls;
  std::vector<double> xk(N, 0.0);
  int counter = -1;
  for(size_t i=0; i<num_search_intvls; i++)
  {
    double x_i = -1.0 + i*delta;
    double x_ip1 = x_i + delta;

    if (Legendre(N,x_i)*Legendre(N,x_ip1) < 0.0)
      xk[++counter] = (x_ip1 + x_i) / 2.0;
  }

  //======================================== Apply algorithm
  // Refer to equation 4.3 in [1]. Sum 1 (S1) is used in the
  // computation of B at x_k. Sum 2 (S2) is used in equation 4.3.
  // Equation 4.3 is broken up into pieces as follows:
  //  - a = block bracket containing the second derivative
  //  - b = denominator
  //  - c = everything but xold
  for (int k=0; k<N; k++)
  {
    for (size_t iteration=0; iteration<max_iters; iteration++)
    {
      double xold = xk[k];
      double f   = Legendre(N,xold);      //Function evaluation
      double fp  = dLegendredx(N,xold);   //First derivative
      double fpp = d2Legendredx2(N,xold); //Second derivative

      //===================== Compute sum 1
      double S1 = 0.0;
      for (int i=0; i<=(k-1); i++)
        S1 += 1.0/(xk[k] - xk[i]);

      //===================== Compute B at x_k
      double B_xk = fp - f*S1;

      //===================== Compute sum 2
      double S2 = 0.0;
      for (int i=0; i<=(k-1); i++)
        S2 += 1.0 / (xk[k] - xk[i]) / (xk[k] - xk[i]);

      //===================== Compute final formula
      double a    = fpp + f*S2;
      double b    = B_xk*B_xk + fp*fp - f*a;
      double c    = 2.0*f*B_xk/b;

      xk[k] = xold - c;

      if (std::fabs(xk[k] - xold) < tol)
        break;
    }//for iteration
  }//for k

  std::stable_sort(xk.begin(), xk.end());

  return xk;
}

/**Populates the abscissae and weights for a Gauss-Legendre
* quadrature given the number of desired quadrature points.
* 
* \param N Is the number of quadrature points.
* \param roots Is a reference to the roots.
* \param max_iters Maximum newton iterations to perform for each root.
*        Default: 1000.
* \param tol Tolerance at which the newton iteration will be terminated.
*        Default: 1.0e-12.
*
* \return A pair with each part of type `std::vector<double>` and equal in size. 
*         The first part is a vector of quadrature points and the second part 
*         is a vector of weights. 
*
* \author Jan*/
inline 
std::pair<std::vector<double>, std::vector<double>> 
GaussLegendreQuadrature(unsigned int N, 
                        bool verbose=false,
                        unsigned int max_iters=1000, 
                        double tol=1.0e-12)
{
  if (verbose)
    std::cout << "Initializing Gauss-Legendre Quadrature "
                      "with " << N << " q-points\n";

  //========================= Compute the roots
  // which are the qpoints
  auto qpoints = drjan::LegendreRoots(N, max_iters, tol);

  //========================= Compute the weights
  std::vector<double> weights(N,1.0);
  for (size_t k=0; k < qpoints.size(); k++)
  {
    weights[k] =
      2.0 * (1.0 - qpoints[k] * qpoints[k]) /
      ((N + 1) * (N + 1) *
        Legendre(N+1, qpoints[k]) * Legendre(N + 1, qpoints[k]) );

    if (verbose)
      std::cout
        << "root[" << k << "]=" << qpoints[k]
        << ", weight=" << weights[k] 
        << "\n";
  }//for abscissae

  return {qpoints, weights};
}

}//namespace drjanv

#endif //DRJANV_QUADRATURES_LEGENDRE_H