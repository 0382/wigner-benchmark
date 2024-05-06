#pragma once
#ifndef _PH_RATIONAL_H_
#define _PH_RATIONAL_H_

#include <gmp.h>
#include <stdint.h>
#include <stdio.h>

#define U16_PRIME_NUMBER 6542

#ifdef __cplusplus
extern "C" {
#endif

extern uint16_t primes_u16[U16_PRIME_NUMBER];

// find the index of the first prime number greater than n
int upper_bound_primes_u16(uint16_t n);

typedef struct _pf_rational
{
    int32_t size;
    int16_t data[];
} pf_rational_t;

extern int pf_integer_number;
extern pf_rational_t **pf_integers;
// initialize pf_integers with the first n integers
void pf_init_integers(int n);
// free pf_integers
void pf_free_integers();

// create a list of `pf_rational_t` representation of integers from `1` to `n`.
pf_rational_t **pf_range_u16(uint16_t n);

// directly init pf_rational_t with uninitialized `size` length data.
pf_rational_t *pf_init(int size);
// init pf_rational_t with `size` length data, and repeat `num` times.
pf_rational_t *pf_init_vec(int size, int num);

// free pf_rational_t, this will also free `pf_init_vec`, since they are allocated together.
void pf_free(pf_rational_t *r);

// pf_rational_t has no uniform size, so use `pf_next` to get the next pf_rational_t.
pf_rational_t *pf_next(pf_rational_t *r);

// make a minimum `pf_rational_t` representation of integer `n`.
pf_rational_t *pf_from_u16(uint16_t n);
// init a `pf_rational_t`, which the maximum prime factor is less `max_integer`.
pf_rational_t *pf_init_max(int max_integer);
// same as `pf_init_max`, but repeat `num` times.
pf_rational_t *pf_init_max_vec(int max_integer, int num);

// get numerator of `pf_rational_t`.
void pf_numerator(mpz_t n, pf_rational_t *r);
// get denominator of `pf_rational_t`.
void pf_denominator(mpz_t n, pf_rational_t *r);
// get numerator and denominator of `pf_rational_t`.
void pf_num_den(mpz_t n, mpz_t d, pf_rational_t *r);

// these functions will not check size
void unsafe_pf_copy(pf_rational_t *r, pf_rational_t *a);
void unsafe_pf_mul(pf_rational_t *r, pf_rational_t *a);
void unsafe_pf_div(pf_rational_t *r, pf_rational_t *a);
void unsafe_pf_pow(pf_rational_t *r, int16_t n);
void unsafe_pf_square(pf_rational_t *r);
void unsafe_pf_inv(pf_rational_t *r);
// pf_gcd means gcd of numerator ans lcm of denominator
void unsafe_pf_gcd(pf_rational_t *r, pf_rational_t *a);
// pf_sgcd means gcd of numerator and denominator
void unsafe_pf_sgcd(pf_rational_t *r, pf_rational_t *a);
// pf_lcm means lcm of numerator and gcd of denominator
void unsafe_pf_lcm(pf_rational_t *r, pf_rational_t *a);

void print_pf_rational(pf_rational_t *r);

void extract_to(pf_rational_t *r, mpz_t n);

#ifdef __cplusplus
}
#endif

#endif // _PH_RATIONAL_H_