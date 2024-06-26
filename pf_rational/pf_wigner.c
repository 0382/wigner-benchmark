#include "pf_wigner.h"
#include <math.h>
#include <string.h>

#define CHECK(b, msg)                        \
    if (!(b))                                \
    {                                        \
        fprintf(stderr, "Error: %s\n", msg); \
        exit(1);                             \
    }

static inline int max(int a, int b) { return a > b ? a : b; }
static inline int min(int a, int b) { return a < b ? a : b; }

void pf_binomial(pf_rational_t *r, int n, int k)
{
    CHECK((unsigned)n <= (unsigned)pf_integer_number, "pf_binomial() requires 0 <= n <= pf_integer_number");
    CHECK((unsigned)k <= (unsigned)n, "pf_binomial() requires 0 <= k <= n");
    pf_set_one(r);
    k = min(n - k, k);
    for (int i = 0; i < k; ++i)
    {
        unsafe_pf_mul(r, pf_integers[n - i - 1]);
        unsafe_pf_div(r, pf_integers[i]);
    }
}

void stagger_sum(mpz_t sum, pf_rational_t *cf, pf_rational_t *xs, int size)
{
    __gmpz_set_ui(sum, 0);
    mpz_t t;
    __gmpz_init(t);
    pf_rational_t *p = xs;
    for (int i = 0; i < size; ++i)
    {
        unsafe_pf_gcd(cf, p);
        p = pf_next(p);
    }
    p = xs;
    for (int i = 0; i < size; ++i)
    {
        unsafe_pf_div(p, cf);
        pf_numerator(t, p);
        if (i % 2 == 0)
        {
            __gmpz_add(sum, sum, t);
        }
        else
        {
            __gmpz_sub(sum, sum, t);
        }
        p = pf_next(p);
    }
    __gmpz_clear(t);
}

void pf_simplity(pf_rational_t *s, pf_rational_t *r, pf_rational_t *x)
{
    CHECK(s->size >= x->size, "pf_simplity() requires s->size >= x->size");
    CHECK(r->size >= x->size, "pf_simplity() requires r->size >= x->size");
    if (s->size > x->size)
    {
        memset(s->data + x->size, 0, (s->size - x->size) * sizeof(int16_t));
    }
    if (r->size > x->size)
    {
        memset(r->data + x->size, 0, (r->size - x->size) * sizeof(int16_t));
    }
    for (int i = 0; i < x->size; ++i)
    {
        s->data[i] = x->data[i] / 2;
        r->data[i] = x->data[i] % 2;
    }
}

void sqrt_rational_init(sqrt_rational_t *x)
{
    __gmpz_init_set_ui(x->sn, 0);
    __gmpz_init_set_ui(x->sd, 1);
    __gmpz_init_set_ui(x->rn, 1);
    __gmpz_init_set_ui(x->rd, 1);
}

void sqrt_rational_clear(sqrt_rational_t *x)
{
    __gmpz_clear(x->sn);
    __gmpz_clear(x->rn);
    __gmpz_clear(x->sd);
    __gmpz_clear(x->rd);
}

void print_sqrt_rational(sqrt_rational_t *x) { __gmp_printf("%Zd/%Zd√(%Zd/%Zd)\n", x->sn, x->sd, x->rn, x->rd); }

double sqrt_rational_to_double(sqrt_rational_t *x)
{
    return __gmpz_get_d(x->sn) / __gmpz_get_d(x->sd) * sqrt(__gmpz_get_d(x->rn) / __gmpz_get_d(x->rd));
}

sqrt_rational_t pf_CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    sqrt_rational_t ans;
    sqrt_rational_init(&ans);
    if (!(check_jm(dj1, dm1) && check_jm(dj2, dm2) && check_jm(dj3, dm3)))
        return ans;
    if (!check_couple(dj1, dj2, dj3))
        return ans;
    if (dm1 + dm2 != dm3)
        return ans;
    const int J = (dj1 + dj2 + dj3) / 2;
    const int jm1 = J - dj1;
    const int jm2 = J - dj2;
    const int jm3 = J - dj3;
    const int j1mm1 = (dj1 - dm1) / 2;
    const int j2mm2 = (dj2 - dm2) / 2;
    const int j3mm3 = (dj3 - dm3) / 2;
    const int j2pm2 = (dj2 + dm2) / 2;
    const int Jsize = upper_bound_primes_u16(J + 1);
    pf_rational_t *t = pf_alloc(Jsize);
    const int low = max(0, max(j1mm1 - jm2, j2pm2 - jm1));
    const int high = min(jm3, min(j1mm1, j2pm2));
    const int Bs_size = high - low + 1;
    pf_rational_t *Bs = pf_alloc_vec(Jsize, Bs_size);
    pf_rational_t *pBs = Bs;
    for (int z = low; z <= high; ++z)
    {
        pf_binomial(pBs, jm3, z);
        pf_binomial(t, jm2, j1mm1 - z);
        unsafe_pf_mul(pBs, t);
        pf_binomial(t, jm1, j2pm2 - z);
        unsafe_pf_mul(pBs, t);
        pBs = pf_next(pBs);
    }
    pf_rational_t *cf = pf_alloc(Jsize);
    mpz_t B;
    __gmpz_init(B);
    stagger_sum(B, cf, Bs, Bs_size);
    if (mpz_sgn(B) == 0)
    {
        __gmpz_clear(B);
        pf_free(t);
        pf_free(Bs);
        pf_free(cf);
        return ans;
    }

    pf_rational_t *A = pf_alloc(Jsize);
    pf_binomial(A, dj1, jm2);
    pf_binomial(t, dj2, jm3);
    unsafe_pf_mul(A, t);
    pf_binomial(t, J + 1, jm3);
    unsafe_pf_div(A, t);
    pf_binomial(t, dj1, j1mm1);
    unsafe_pf_div(A, t);
    pf_binomial(t, dj2, j2mm2);
    unsafe_pf_div(A, t);
    pf_binomial(t, dj3, j3mm3);
    unsafe_pf_div(A, t);

    extract_to(cf, B);
    pf_rational_t *s = pf_alloc(Jsize);
    pf_rational_t *r = pf_alloc(Jsize);
    pf_simplity(s, r, A);
    unsafe_pf_mul(s, cf);

    pf_rational_t *g = t;
    unsafe_pf_copy(g, r);
    unsafe_pf_inv(g);
    unsafe_pf_sgcd(g, s);
    unsafe_pf_div(s, g);
    unsafe_pf_mul(r, g);
    unsafe_pf_mul(r, g);

    pf_num_den(ans.sn, ans.sd, s);
    pf_num_den(ans.rn, ans.rd, r);
    __gmpz_mul(ans.sn, ans.sn, B);
    if (isodd(low))
        __gmpz_neg(ans.sn, ans.sn);
    __gmpz_clear(B);
    pf_free(A);
    pf_free(t);
    pf_free(Bs);
    pf_free(cf);
    pf_free(s);
    pf_free(r);
    return ans;
}

sqrt_rational_t pf_3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    sqrt_rational_t ans = pf_CG(dj1, dj2, dj3, -dm1, -dm2, dm3);
    if (isodd((dj3 + dm3) / 2 + dj1))
    {
        __gmpz_neg(ans.sn, ans.sn);
    }
    __gmpz_mul_ui(ans.rd, ans.rd, dj3 + 1);
    return ans;
}

sqrt_rational_t pf_6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6)
{
    sqrt_rational_t ans;
    sqrt_rational_init(&ans);
    if (!(check_couple(dj1, dj2, dj3) && check_couple(dj1, dj5, dj6) && check_couple(dj4, dj2, dj6) &&
          check_couple(dj4, dj5, dj3)))
        return ans;
    const int j123 = (dj1 + dj2 + dj3) / 2;
    const int j156 = (dj1 + dj5 + dj6) / 2;
    const int j426 = (dj4 + dj2 + dj6) / 2;
    const int j453 = (dj4 + dj5 + dj3) / 2;
    const int jpm123 = (dj1 + dj2 - dj3) / 2;
    const int jpm132 = (dj1 + dj3 - dj2) / 2;
    const int jpm231 = (dj2 + dj3 - dj1) / 2;
    const int jpm156 = (dj1 + dj5 - dj6) / 2;
    const int jpm426 = (dj4 + dj2 - dj6) / 2;
    const int jpm453 = (dj4 + dj5 - dj3) / 2;

    const int low = max(max(j123, j156), max(j426, j453));
    const int high = min(jpm123 + j453, min(jpm132 + j426, jpm231 + j156));
    const int max_J = max(max(high + 1, jpm123), max(jpm132, jpm231));
    const int max_size = upper_bound_primes_u16(max_J);
    const int Bs_size = high - low + 1;
    pf_rational_t *Bs = pf_alloc_vec(max_size, Bs_size);
    pf_rational_t *pB = Bs;
    pf_rational_t *t = pf_alloc(max_size);
    for (int x = low; x <= high; ++x)
    {
        pf_binomial(pB, x + 1, j123 + 1);
        pf_binomial(t, jpm123, x - j453);
        unsafe_pf_mul(pB, t);
        pf_binomial(t, jpm132, x - j426);
        unsafe_pf_mul(pB, t);
        pf_binomial(t, jpm231, x - j156);
        unsafe_pf_mul(pB, t);
        pB = pf_next(pB);
    }
    pf_rational_t *cf = pf_alloc(max_size);
    mpz_t B;
    __gmpz_init(B);
    stagger_sum(B, cf, Bs, Bs_size);
    if (mpz_sgn(B) == 0)
    {
        __gmpz_clear(B);
        pf_free(Bs);
        pf_free(t);
        pf_free(cf);
        return ans;
    }
    extract_to(cf, B);

    pf_rational_t *A = pf_alloc(max_size);
    pf_binomial(A, j123 + 1, dj1 + 1);
    pf_binomial(t, dj1, jpm123);
    unsafe_pf_mul(A, t);
    pf_binomial(t, j156 + 1, dj1 + 1);
    unsafe_pf_div(A, t);
    pf_binomial(t, dj1, jpm156);
    unsafe_pf_div(A, t);
    pf_binomial(t, j453 + 1, dj4 + 1);
    unsafe_pf_div(A, t);
    pf_binomial(t, dj4, jpm453);
    unsafe_pf_div(A, t);
    pf_binomial(t, j426 + 1, dj4 + 1);
    unsafe_pf_div(A, t);
    pf_binomial(t, dj4, jpm426);
    unsafe_pf_div(A, t);

    pf_rational_t *s = pf_alloc(max_size);
    pf_rational_t *r = pf_alloc(max_size);
    pf_rational_t *pfdj4 = pf_integers[dj4];
    pf_simplity(s, r, A);
    unsafe_pf_mul(s, cf);
    unsafe_pf_div(s, pfdj4);

    pf_rational_t *g = t;
    unsafe_pf_copy(g, r);
    unsafe_pf_inv(g);
    unsafe_pf_sgcd(g, s);
    unsafe_pf_div(s, g);
    unsafe_pf_mul(r, g);
    unsafe_pf_mul(r, g);

    pf_num_den(ans.sn, ans.sd, s);
    pf_num_den(ans.rn, ans.rd, r);
    __gmpz_mul(ans.sn, ans.sn, B);
    if (isodd(low))
        __gmpz_neg(ans.sn, ans.sn);
    __gmpz_clear(B);
    pf_free(A);
    pf_free(Bs);
    pf_free(t);
    pf_free(cf);
    pf_free(s);
    pf_free(r);
    return ans;
}

double dpf_CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    sqrt_rational_t x = pf_CG(dj1, dj2, dj3, dm1, dm2, dm3);
    double ans = sqrt_rational_to_double(&x);
    sqrt_rational_clear(&x);
    return ans;
}

double dpf_3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    sqrt_rational_t x = pf_3j(dj1, dj2, dj3, dm1, dm2, dm3);
    double ans = sqrt_rational_to_double(&x);
    sqrt_rational_clear(&x);
    return ans;
}

double dpf_6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6)
{
    sqrt_rational_t x = pf_6j(dj1, dj2, dj3, dj4, dj5, dj6);
    double ans = sqrt_rational_to_double(&x);
    sqrt_rational_clear(&x);
    return ans;
}
