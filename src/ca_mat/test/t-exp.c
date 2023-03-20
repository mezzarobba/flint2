/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

/*

iter = 49
n = 2
rand
1
success
s1
in small check
not QQ
not special
not repr
not algebraic number field
d1
d2
d3
d4
d5
d6
d7
s2
s3
first mul
second mul
inside ca_mat_check_equal
i = 0, j = 
in small check
not QQ
not special
after small check
i = 0, j = 
in small check
not QQ
not special
after small check
i = 1, j = 
in small check
after small check
i = 1, j = 
in small check
not QQ
not special
after small check
after success
in small check
not QQ
not special
not repr
not algebraic number field
d1
d2
d3
d4
d5
d7
in small check
not QQ
not special
not repr
not algebraic number field
d1
d2
d3
d4
d5
d7
true
inside ca_mat_check_equal
i = 0, j = 
in small check
not QQ
not special
not repr
not algebraic number field
d1
d2
d3
d4
d5
d6
d7
after small check
i = 0, j = 
in small check
not QQ
not special
not repr
not algebraic number field
d1
d2
d3
d4
d5
d6

*/

int main()
{
    slong iter;
    flint_rand_t state;

    ca_ctx_t ctx;
    ca_mat_t A, B, P, Q, C, D;
    ca_t d, t;
    slong n;
    int success;

    flint_printf("exp...");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        if (iter == 49)
        {
            flint_printf("\n\n\n\n\n\n==========================================================================================\n\n\n\n\n\n\n");
        }

        flint_printf("iter = %d\n", iter);
        fflush(stdout);

        ca_ctx_init(ctx);
        n = n_randint(state, 3);

        flint_printf("n = %wd\n", n);
        fflush(stdout);

        ca_mat_init(A, n, n, ctx);
        ca_mat_init(B, n, n, ctx);
        ca_mat_init(P, n, n, ctx);
        ca_mat_init(Q, n, n, ctx);
        ca_mat_init(C, n, n, ctx);
        ca_mat_init(D, n, n, ctx);
        ca_init(d, ctx);
        ca_init(t, ctx);

        flint_printf("rand\n");
        fflush(stdout);

        if (n <= 2 && n_randint(state, 5) == 0)
            ca_mat_randtest(A, state, 1, 5, ctx);
        else
            ca_mat_randtest_rational(A, state, 5, ctx);

        if (n_randint(state, 2))
        {
            flint_printf("1\n");
            fflush(stdout);

            success = ca_mat_exp(B, A, ctx);
        }
        else
        {
            flint_printf("0\n");
            fflush(stdout);

            ca_mat_set(B, A, ctx);
            success = ca_mat_exp(B, B, ctx);
        }

        if (success)
        {
            flint_printf("success\n");
            fflush(stdout);

            ca_mat_det(d, B, ctx);

            ca_mat_trace(t, A, ctx);
            ca_exp(t, t, ctx);

            flint_printf("s1\n");
            fflush(stdout);

            if (ca_check_equal(d, t, ctx) == T_FALSE)
            {
                flint_printf("FAIL (trace)\n");
                flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                flint_printf("B: "); ca_mat_print(B, ctx); flint_printf("\n");
                flint_printf("d: "); ca_print(d, ctx); flint_printf("\n");
                flint_printf("t: "); ca_print(t, ctx); flint_printf("\n");
                flint_abort();
            }

            do {
                ca_mat_randtest_rational(P, state, 1, ctx);
            } while (ca_mat_inv(Q, P, ctx) != T_TRUE);

            flint_printf("s2\n");
            fflush(stdout);

            ca_mat_mul(C, P, A, ctx);
            ca_mat_mul(C, C, Q, ctx);
            success = ca_mat_exp(C, C, ctx);

            flint_printf("s3\n");
            fflush(stdout);

            if (success)
            {
                flint_printf("first mul\n");
                fflush(stdout);

                ca_mat_mul(D, P, B, ctx);
                flint_printf("second mul\n");
                fflush(stdout);

                ca_mat_mul(D, D, Q, ctx);

                if (ca_mat_check_equal(C, D, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (similarity transform)\n");
                    flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("B: "); ca_mat_print(B, ctx); flint_printf("\n");
                    flint_printf("C: "); ca_mat_print(C, ctx); flint_printf("\n");
                    flint_printf("D: "); ca_mat_print(D, ctx); flint_printf("\n");
                    flint_printf("P: "); ca_mat_print(P, ctx); flint_printf("\n");
                    flint_printf("Q: "); ca_mat_print(Q, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        flint_printf("after success\n");
        fflush(stdout);

        if (ca_mat_log(B, A, ctx) == T_TRUE)
        {
            success = ca_mat_exp(C, B, ctx);

            flint_printf("true\n");
            fflush(stdout);

            if (success)
            {
                if (ca_mat_check_equal(A, C, ctx) == T_FALSE)
                {
                    flint_printf("going to fail\n");
                    fflush(stdout);

                    flint_printf("FAIL (logarithm)\n");
                    flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("B: "); ca_mat_print(B, ctx); flint_printf("\n");
                    flint_printf("C: "); ca_mat_print(C, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        flint_printf("clear\n");
        fflush(stdout);

        ca_mat_clear(A, ctx);
        ca_mat_clear(B, ctx);
        ca_mat_clear(P, ctx);
        ca_mat_clear(Q, ctx);
        ca_mat_clear(C, ctx);
        ca_mat_clear(D, ctx);
        flint_printf("mid clear\n");
        fflush(stdout);

        ca_clear(d, ctx);
        ca_clear(t, ctx);

        ca_ctx_clear(ctx);

        if (iter == 49)
            break;
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
