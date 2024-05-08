#pragma once
#ifndef PF_WIRGER_H
#define PF_WIRGER_H

#include "pf_rational.h"
#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
using _Bool = bool;
#endif

// judge if a number is a odd number
static inline _Bool isodd(int x) { return x % 2 != 0; }
// judge if a number is a even number
static _Bool iseven(int x) { return x % 2 == 0; }
// judge if two number are same odd or same even
static _Bool is_same_parity(int x, int y) { return iseven(x ^ y); }
// return (-1)^n
static int iphase(int x) { return iseven(x) - isodd(x); }
// check if m-quantum number if one of the components of a the j-quantum number
static _Bool check_jm(int dj, int dm) { return is_same_parity(dj, dm) && (abs(dm) <= dj); }
// judge if three angular momentum can couple
static _Bool check_couple(int dj1, int dj2, int dj3)
{
    return dj1 >= 0 && dj2 >= 0 && is_same_parity(dj1 + dj2, dj3) && (dj3 <= (dj1 + dj2)) && (dj3 >= abs(dj1 - dj2));
}

// \sum (-1)^i xs[i]
// result is `cf * sum`
void stagger_sum(mpz_t sum, pf_rational_t *cf, pf_rational_t *xs, int size);

void pf_binomial(pf_rational_t *r, int n, int k);

void pf_simplity(pf_rational_t *s, pf_rational_t *r, pf_rational_t *x);

typedef struct _sqrt_rational
{
    // sn/sd * \sqrt{rn/rd}
    mpz_t sn, rn, sd, rd;
} sqrt_rational_t;

void sqrt_rational_init(sqrt_rational_t *x);
void sqrt_rational_clear(sqrt_rational_t *x);
void print_sqrt_rational(sqrt_rational_t *x);
double sqrt_rational_to_double(sqrt_rational_t *x);

sqrt_rational_t pf_CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
sqrt_rational_t pf_3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
sqrt_rational_t pf_6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
sqrt_rational_t pf_Racah(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
sqrt_rational_t pf_9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9);

double dpf_CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
double dpf_3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
double dpf_6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
double dpf_Racah(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
double dpf_9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9);

#ifdef __cplusplus
}
#endif

#endif // PF_WIRGER_H