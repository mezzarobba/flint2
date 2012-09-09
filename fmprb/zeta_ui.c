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

#include <math.h>
#include "arith.h"
#include "fmprb.h"

void
fmprb_zeta_ui(fmprb_t x, ulong n, long prec)
{
    if (n == 0)
    {
        fmprb_set_si(x, -1);
        fmprb_mul_2exp_si(x, x, -1);
    }
    else if (n == 1)
    {
        printf("exception: zeta_ui(1)\n");
        abort();
    }
    /* asymptotic case */
    else if (n > 6 && n > 0.7 * prec)
    {
        fmprb_zeta_ui_euler_product(x, n, prec);
    }
    else if (n == 3)
    {
        fmprb_const_zeta3_bsplit(x, prec);
    }
    /* small even n */
    else if ((n % 2 == 0) && (n < 40 + 0.11*prec))
    {
        fmprb_zeta_ui_bernoulli(x, n, prec);
    }
    /* small odd n, extremely high precision */
    else if (n < prec * 0.0006)
    {
        fmprb_zeta_ui_bsplit(x, n, prec);
    }
    /* large n */
    else if (prec > 20 && n > 0.4 * pow(prec, 0.8))
    {
        fmprb_zeta_ui_euler_product(x, n, prec);
    }
    /* fallback */
    else
    {
        fmprb_zeta_ui_vec_borwein(x, n, 1, 0, prec);
    }
}
