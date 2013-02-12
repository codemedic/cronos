/*
 * -------------------------------------------------------------------
 *
 * Copyright 2004 Anthony Brockwell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * -------------------------------------------------------------------
 *
 *
 * This module finds the roots
 * of polynomials.  It was buggy, so I replaced
 * the core of it with calls to gsl routines on Jun, 2004.
 * Anthony Brockwell.
 *
 * polynomial format:
 *  a_0 + a_1 z + a_2 z^2 + ...
 *  a_0 = first coefficient, a_1 = second coefficient, ...
 */

#include <math.h>
#include <gsl/gsl_poly.h>
#include "poly.h"

namespace mslib {

// The following procedure creates a new
// complex array containing all the roots of
// the specified polynomial.

complex *rootsof(int order, double *coeffs)
{
  int i;
  complex *roots = new complex[order];
  double z[order*2];
  gsl_poly_complex_workspace *gwp =
    gsl_poly_complex_workspace_alloc(order+1);
  gsl_poly_complex_solve(coeffs, order+1, gwp, z);
  gsl_poly_complex_workspace_free(gwp);
  for (i=0 ; i<order ; ++i)
    roots[i] = complex(z[2*i],z[2*i+1]);
  return roots;
}

// The next procedure does the opposite.  Given
// the roots of a polynomial, it constructs the coefficients
// of the polynomial in which the $a_p$ coefficient is
// equal to 1.

complex *coeffsof(int order, complex *roots, complex normalize)
// This procedure works recursively
{
  complex *coeffs = new complex[order+1];
  if (order==1)
    {
      coeffs[0] = -roots[0];
      coeffs[1] = 1.0;
    }
  else
    {
      complex *sub = coeffsof(order-1, roots), c1, c2, c3;
      c1 = -1.0;
      c2 = sub[0];
      c3 = roots[order-1];
      coeffs[0] = c1*c2*c3;
      for (int i=1 ; i<order ; ++i)
        coeffs[i] = c1*sub[i]*c3 + sub[i-1];
      coeffs[order] = 1.0;
      delete[] sub;
    }
  for (int i=0 ; i<=order ; ++i)
    coeffs[i] *= normalize;
  return coeffs;
}

} // end of mslib namespace
