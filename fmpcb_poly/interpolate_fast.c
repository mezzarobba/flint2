/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmpcb_poly.h"

void
_fmpcb_poly_interpolation_weights(fmpcb_struct * w,
    fmpcb_struct ** tree, long len, long prec)
{
    fmpcb_struct * tmp;
    long i, n, height;

    if (len == 0)
        return;

    if (len == 1)
    {
        fmpcb_one(w);
        return;
    }

    tmp = _fmpcb_vec_init(len + 1);
    height = FLINT_CLOG2(len);
    n = 1L << (height - 1);

    _fmpcb_poly_mul_monic(tmp, tree[height-1], n + 1,
                        tree[height-1] + (n + 1), (len - n + 1), prec);

    _fmpcb_poly_derivative(tmp, tmp, len + 1, prec);
    _fmpcb_poly_evaluate_vec_fast_precomp(w, tmp, len, tree, len, prec);

    for (i = 0; i < len; i++)
        fmpcb_inv(w + i, w + i, prec);

    _fmpcb_vec_clear(tmp, len + 1);
}

void
_fmpcb_poly_interpolate_fast_precomp(fmpcb_struct * poly,
    const fmpcb_struct * ys, fmpcb_struct ** tree, const fmpcb_struct * weights,
    long len, long prec)
{
    fmpcb_struct *t, *u, *pa, *pb;
    long i, pow, left;

    if (len == 0)
        return;

    t = _fmpcb_vec_init(len);
    u = _fmpcb_vec_init(len);

    for (i = 0; i < len; i++)
        fmpcb_mul(poly + i, weights + i, ys + i, prec);

    for (i = 0; i < FLINT_CLOG2(len); i++)
    {
        pow = (1L << i);
        pa = tree[i];
        pb = poly;
        left = len;

        while (left >= 2 * pow)
        {
            _fmpcb_poly_mul(t, pa, pow + 1, pb + pow, pow, prec);
            _fmpcb_poly_mul(u, pa + pow + 1, pow + 1, pb, pow, prec);
            _fmpcb_vec_add(pb, t, u, 2 * pow, prec);

            left -= 2 * pow;
            pa += 2 * pow + 2;
            pb += 2 * pow;
        }

        if (left > pow)
        {
            _fmpcb_poly_mul(t, pa, pow + 1, pb + pow, left - pow, prec);
            _fmpcb_poly_mul(u, pb, pow, pa + pow + 1, left - pow + 1, prec);
            _fmpcb_vec_add(pb, t, u, left, prec);
        }
    }

    _fmpcb_vec_clear(t, len);
    _fmpcb_vec_clear(u, len);
}

void
_fmpcb_poly_interpolate_fast(fmpcb_struct * poly,
    const fmpcb_struct * xs, const fmpcb_struct * ys, long len, long prec)
{
    fmpcb_struct ** tree;
    fmpcb_struct *w;

    tree = _fmpcb_poly_tree_alloc(len);
    _fmpcb_poly_tree_build(tree, xs, len, prec);

    w = _fmpcb_vec_init(len);
    _fmpcb_poly_interpolation_weights(w, tree, len, prec);

    _fmpcb_poly_interpolate_fast_precomp(poly, ys, tree, w, len, prec);

    _fmpcb_vec_clear(w, len);
    _fmpcb_poly_tree_free(tree, len);
}

void
fmpcb_poly_interpolate_fast(fmpcb_poly_t poly,
        const fmpcb_struct * xs, const fmpcb_struct * ys, long n, long prec)
{
    if (n == 0)
    {
        fmpcb_poly_zero(poly);
    }
    else
    {
        fmpcb_poly_fit_length(poly, n);
        _fmpcb_poly_set_length(poly, n);
        _fmpcb_poly_interpolate_fast(poly->coeffs, xs, ys, n, prec);
        _fmpcb_poly_normalise(poly);
    }
}
