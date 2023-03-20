/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

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
        flint_printf("iter = %d", iter);
        fflush(stdout);

        ca_ctx_init(ctx);
        n = n_randint(state, 3);

        flint_printf("init", iter);
        fflush(stdout);

        ca_mat_init(A, n, n, ctx);
        ca_mat_init(B, n, n, ctx);
        ca_mat_init(P, n, n, ctx);
        ca_mat_init(Q, n, n, ctx);
        ca_mat_init(C, n, n, ctx);
        ca_mat_init(D, n, n, ctx);
        ca_init(d, ctx);
        ca_init(t, ctx);

        flint_printf("rand", iter);
        fflush(stdout);

        if (n <= 2 && n_randint(state, 5) == 0)
            ca_mat_randtest(A, state, 1, 5, ctx);
        else
            ca_mat_randtest_rational(A, state, 5, ctx);

        if (n_randint(state, 2))
        {
            flint_printf("1", iter);
            fflush(stdout);

            success = ca_mat_exp(B, A, ctx);
        }
        else
        {
            flint_printf("0", iter);
            fflush(stdout);

            ca_mat_set(B, A, ctx);
            success = ca_mat_exp(B, B, ctx);
        }

        if (success)
        {
            flint_printf("success", iter);
            fflush(stdout);

            ca_mat_det(d, B, ctx);

            ca_mat_trace(t, A, ctx);
            ca_exp(t, t, ctx);

            flint_printf("s1", iter);
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

            flint_printf("s2", iter);
            fflush(stdout);

            ca_mat_mul(C, P, A, ctx);
            ca_mat_mul(C, C, Q, ctx);
            success = ca_mat_exp(C, C, ctx);

            flint_printf("s3", iter);
            fflush(stdout);

            if (success)
            {
                ca_mat_mul(D, P, B, ctx);
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

        flint_printf("after success", iter);
        fflush(stdout);

        if (ca_mat_log(B, A, ctx) == T_TRUE)
        {
            success = ca_mat_exp(C, B, ctx);

            flint_printf("true", iter);
            fflush(stdout);

            if (success)
            {
                if (ca_mat_check_equal(A, C, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (logarithm)\n");
                    flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("B: "); ca_mat_print(B, ctx); flint_printf("\n");
                    flint_printf("C: "); ca_mat_print(C, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        flint_printf("clear", iter);
        fflush(stdout);

        ca_mat_clear(A, ctx);
        ca_mat_clear(B, ctx);
        ca_mat_clear(P, ctx);
        ca_mat_clear(Q, ctx);
        ca_mat_clear(C, ctx);
        ca_mat_clear(D, ctx);
        flint_printf("mid clear", iter);
        fflush(stdout);

        ca_clear(d, ctx);
        ca_clear(t, ctx);

        ca_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
