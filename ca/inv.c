/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_inv(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    truth_t is_zero;
    slong field_index;
    ca_field_type_t type;
    ulong xfield;

    xfield = x->field;

    if (xfield == CA_FIELD_ID_QQ)
    {
        if (fmpq_is_zero(CA_FMPQ(x)))
        {
            ca_uinf(res, ctx);
        }
        else
        {
            _ca_make_fmpq(res, ctx);
            fmpq_inv(CA_FMPQ(res), CA_FMPQ(x));
        }
        return;
    }

    if (CA_IS_SPECIAL(x))
    {
        if ((xfield & CA_SIGNED_INF) || (xfield & CA_UNSIGNED_INF))
            ca_zero(res, ctx);
        else
            ca_set(res, x, ctx);
        return;
    }

    is_zero = ca_check_is_zero(x, ctx);

    if (is_zero == T_TRUE)
    {
        ca_uinf(res, ctx);
        return;
    }
    else if (is_zero == T_UNKNOWN)
    {
        ca_unknown(res, ctx);
        return;
    }

    field_index = xfield;
    type = ctx->fields[field_index].type;

    _ca_make_field_element(res, field_index, ctx);

    if (type == CA_FIELD_TYPE_QQ)
    {
        fmpq_inv(CA_FMPQ(res), CA_FMPQ(x));
    }
    else if (type == CA_FIELD_TYPE_NF)
    {
        nf_elem_inv(CA_NF_ELEM(res), CA_NF_ELEM(x), CA_FIELD_NF(ctx->fields + field_index));
    }
    else if (type == CA_FIELD_TYPE_FUNC)
    {
        fmpz_mpoly_q_inv(CA_MPOLY_Q(res), CA_MPOLY_Q(x), ctx->mctx + 0);
    }
    else if (type == CA_FIELD_TYPE_MULTI)
    {
        fmpz_mpoly_q_inv(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_FIELD_MCTX(ctx->fields + field_index, ctx));
    }
}

