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
 */

// This module handles operations with polynomials.
// Polynomials are specified by two parameters:
//  (a) The order is p.
//  (b) $a_0, a_1, \ldots, a_p$ are coefficients of the
//      constant term, the $x$ term, the $x^2$ term, up to the $x^p$ term.

#ifndef poly_h
#define poly_h

#include <complex>
#include "matrix.h"

#define complex std::complex<double>

namespace mslib {

complex *rootsof(int order, double *coeffs);
complex *coeffsof(int order, complex *root, complex normalize=1.0);

}

#endif
