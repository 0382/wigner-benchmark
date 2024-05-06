#include "pf_rational.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(b, msg)                        \
    if (!(b))                                \
    {                                        \
        fprintf(stderr, "Error: %s\n", msg); \
        exit(1);                             \
    }

static inline int min_int(int a, int b) { return a < b ? a : b; }
static inline int max_int(int a, int b) { return a > b ? a : b; }
static inline int min_i16(int16_t a, int16_t b) { return a < b ? a : b; }
static inline int max_i16(int16_t a, int16_t b) { return a > b ? a : b; }

uint16_t primes_u16[U16_PRIME_NUMBER] = {
#include "primes_u16.h"
};

int pf_integer_number = 0;
pf_rational_t **pf_integers = NULL;

int upper_bound_primes_u16(uint16_t n)
{
    assert(pf_integers != NULL);
    int l = 0, r = U16_PRIME_NUMBER - 1;
    while (l < r)
    {
        int m = (l + r) / 2;
        if (primes_u16[m] <= n)
        {
            l = m + 1;
        }
        else
        {
            r = m;
        }
    }
    return l;
}

void pf_init_integers(int n)
{
    pf_integer_number = n;
    pf_integers = pf_range_u16(pf_integer_number);
}

void pf_free_integers()
{
    free(pf_integers[0]);
    free(pf_integers);
}

pf_rational_t **pf_range_u16(uint16_t n)
{
    CHECK(n > 0, "pf_range_u16(0) is undefined");
    int mp_idx = upper_bound_primes_u16(n);
    int16_t *data = malloc(mp_idx * sizeof(int16_t));
    CHECK(data != NULL, "could not allocate memory for pf_range_u16()");
    pf_rational_t **r = malloc(n * sizeof(pf_rational_t *));
    CHECK(r != NULL, "could not allocate memory for pf_range_u16()");
    // size_t size = n * (sizeof(pf_rational_t) + mp_idx * sizeof(int16_t));
    // this is not the exact approx function, but it works
    size_t size = n * sizeof(pf_rational_t) + floor(0.02 * pow(n, 1.5) * log(n) * log(n) + 6) * sizeof(int16_t);
    pf_rational_t *p0 = malloc(size);
    CHECK(p0 != NULL, "could not allocate memory for pf_range_u16()");
    pf_rational_t *p = p0;
    for (int m = 1; m <= n; ++m)
    {
        memset(data, 0, mp_idx * sizeof(int16_t));
        int i = 0;
        int t = m;
        while (t > 1)
        {
            while (t % primes_u16[i] == 0)
            {
                t /= primes_u16[i];
                ++data[i];
            }
            ++i;
        }
        p->size = i;
        memcpy(p->data, data, i * sizeof(int16_t));
        r[m - 1] = p;
        p = pf_next(p);
    }
    free(data);
    return r;
}

pf_rational_t *pf_init(int size)
{
    pf_rational_t *r = malloc(sizeof(pf_rational_t) + size * sizeof(int16_t));
    CHECK(r != NULL, "could not allocate memory for pf_init()");
    r->size = size;
    return r;
}

pf_rational_t *pf_init_vec(int size, int num)
{
    pf_rational_t *r = malloc(num * (sizeof(pf_rational_t) + size * sizeof(int16_t)));
    CHECK(r != NULL, "could not allocate memory for pf_init_vec()");
    pf_rational_t *p = r;
    for (int i = 0; i < num; ++i)
    {
        p->size = size;
        p = pf_next(p);
    }
    return r;
}

void pf_free(pf_rational_t *r) { free(r); }

pf_rational_t *pf_next(pf_rational_t *r)
{
    size_t size = sizeof(pf_rational_t) + r->size * sizeof(int16_t);
    void *p = (void *)r + size;
    return (pf_rational_t *)p;
}

pf_rational_t *pf_from_u16(uint16_t n)
{
    CHECK(n > 0, "pf_from_u16(0) is undefined");
    int mp_idx = upper_bound_primes_u16(n);
    int16_t *data = malloc(mp_idx * sizeof(int16_t));
    CHECK(data != NULL, "could not allocate memory for pf_from_u16()");
    memset(data, 0, mp_idx * sizeof(int16_t));
    int i = 0;
    while (n > 1)
    {
        while (n % primes_u16[i] == 0)
        {
            n /= primes_u16[i];
            ++data[i];
        }
        ++i;
    }
    pf_rational_t *r = pf_init(i);
    memcpy(r->data, data, i * sizeof(int16_t));
    free(data);
    return r;
}

pf_rational_t *pf_init_max(int max_integer)
{
    int size = upper_bound_primes_u16(max_integer);
    return pf_init(size);
}

pf_rational_t *pf_init_max_vec(int max_integer, int num)
{
    int size = upper_bound_primes_u16(max_integer);
    return pf_init_vec(size, num);
}

void pf_numerator(mpz_t n, pf_rational_t *r)
{
    __gmpz_set_ui(n, 1);
    mpz_t t;
    __gmpz_init(t);
    for (int i = 0; i < r->size; ++i)
    {
        if (r->data[i] > 0)
        {
            __gmpz_ui_pow_ui(t, primes_u16[i], r->data[i]);
            __gmpz_mul(n, n, t);
        }
    }
    __gmpz_clear(t);
}

void pf_denominator(mpz_t n, pf_rational_t *r)
{
    __gmpz_set_ui(n, 1);
    mpz_t t;
    __gmpz_init(t);
    for (int i = 0; i < r->size; ++i)
    {
        if (r->data[i] < 0)
        {
            __gmpz_ui_pow_ui(t, primes_u16[i], -r->data[i]);
            __gmpz_mul(n, n, t);
        }
    }
    __gmpz_clear(t);
}

void pf_num_den(mpz_t n, mpz_t d, pf_rational_t *r)
{
    __gmpz_set_ui(n, 1);
    __gmpz_set_ui(d, 1);
    mpz_t t;
    __gmpz_init(t);
    for (int i = 0; i < r->size; ++i)
    {
        if (r->data[i] > 0)
        {
            __gmpz_ui_pow_ui(t, primes_u16[i], r->data[i]);
            __gmpz_mul(n, n, t);
        }
        else if (r->data[i] < 0)
        {
            __gmpz_ui_pow_ui(t, primes_u16[i], -r->data[i]);
            __gmpz_mul(d, d, t);
        }
    }
    __gmpz_clear(t);
}

void unsafe_pf_copy(pf_rational_t *r, pf_rational_t *a)
{
    // assume a->size <= r->size
    memcpy(r->data, a->data, a->size * sizeof(int16_t));
    memset(r->data + a->size, 0, (r->size - a->size) * sizeof(int16_t));
}

void unsafe_pf_mul(pf_rational_t *r, pf_rational_t *a)
{
    // assume a->size <= r->size
    for (int i = 0; i < a->size; ++i)
    {
        r->data[i] += a->data[i];
    }
}

void unsafe_pf_div(pf_rational_t *r, pf_rational_t *a)
{
    // assume a->size <= r->size
    for (int i = 0; i < a->size; ++i)
    {
        r->data[i] -= a->data[i];
    }
}

void unsafe_pf_pow(pf_rational_t *r, int16_t n)
{
    for (int i = 0; i < r->size; ++i)
    {
        r->data[i] *= n;
    }
}

void unsafe_pf_square(pf_rational_t *r)
{
    for (int i = 0; i < r->size; ++i)
    {
        r->data[i] *= 2;
    }
}

void unsafe_pf_inv(pf_rational_t *r)
{
    for (int i = 0; i < r->size; ++i)
    {
        r->data[i] = -r->data[i];
    }
}

void unsafe_pf_gcd(pf_rational_t *r, pf_rational_t *a)
{
    // assume a->size <= r->size
    for (int i = 0; i < a->size; ++i)
    {
        r->data[i] = min_int(r->data[i], a->data[i]);
    }
    for (int i = a->size; i < r->size; ++i)
    {
        r->data[i] = min_int(r->data[i], 0);
    }
}

void unsafe_pf_sgcd(pf_rational_t *r, pf_rational_t *a)
{
    // assume a->size <= r->size
    for (int i = 0; i < a->size; ++i)
    {
        int16_t tmin = min_i16(r->data[i], a->data[i]);
        int16_t tmax = max_i16(r->data[i], a->data[i]);
        if (tmin > 0)
        {
            r->data[i] = tmin;
        }
        else if (tmax < 0)
        {
            r->data[i] = tmax;
        }
        else
        {
            r->data[i] = 0;
        }
    }
    for (int i = a->size; i < r->size; ++i)
    {
        r->data[i] = 0;
    }
}

void unsafe_pf_lcm(pf_rational_t *r, pf_rational_t *a)
{
    // assume a->size <= r->size
    for (int i = 0; i < a->size; ++i)
    {
        r->data[i] = max_i16(r->data[i], a->data[i]);
    }
    for (int i = a->size; i < r->size; ++i)
    {
        r->data[i] = max_i16(r->data[i], 0);
    }
}

void pf_binomial(pf_rational_t *r, int n, int k)
{
    CHECK((unsigned)n <= (unsigned)pf_integer_number, "pf_binomial() requires 0 <= n <= pf_integer_number");
    CHECK((unsigned)k <= (unsigned)n, "pf_binomial() requires 0 <= k <= n");
    memset(r->data, 0, r->size * sizeof(int16_t));
    k = min_int(n - k, k);
    for (int i = 0; i < k; ++i)
    {
        unsafe_pf_mul(r, pf_integers[n - i - 1]);
        unsafe_pf_div(r, pf_integers[i]);
    }
}

void extract_to(pf_rational_t *pf, mpz_t n)
{
    if (__gmpz_cmp_ui(n, 1) == 0)
    {
        return;
    }
    mpz_t t;
    __gmpz_init(t);
    for (int i = 0; i < pf->size; ++i)
    {
        unsigned long r = __gmpz_fdiv_q_ui(t, n, primes_u16[i]);
        while (r == 0)
        {
            ++pf->data[i];
            __gmpz_set(n, t);
            r = __gmpz_fdiv_q_ui(t, n, primes_u16[i]);
        }
    }
    __gmpz_clear(t);
}

void print_pf_rational(pf_rational_t *r)
{
    mpz_t n, d;
    __gmpz_init(n);
    __gmpz_init(d);
    pf_num_den(n, d, r);
    __gmp_printf("%Zd/%Zd\n", n, d);
    __gmpz_clear(n);
    __gmpz_clear(d);
}
