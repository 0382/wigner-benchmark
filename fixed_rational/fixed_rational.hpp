#pragma once
#ifndef _FIXED_RATIONAL_H_
#define _FIXED_RATIONAL_H_

#include <array>
#include <cmath>
#include <immintrin.h>
#include <iostream>
#include <vector>

namespace fr136
{

constexpr std::size_t PN = 32;
using fr136_t = std::array<int8_t, PN>;

// clang-format off
constexpr std::array<int, PN> primes32 = {
    2, 3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73, 79, 83, 89,
    97, 101, 103, 107, 109, 113, 127, 131
};
// clang-format on

using i128 = __int128_t;
using u128 = __uint128_t;

inline std::ostream &operator<<(std::ostream &os, i128 x)
{
    int64_t y = x / 10000000000000000;
    int64_t z = x % 10000000000000000;
    if (y != 0)
    {
        os << y;
    }
    os << z;
}

double as_double(i128 x)
{
    int64_t y = x / 10000000000000000;
    int64_t z = x % 10000000000000000;
    return z + y * 1e16;
}

constexpr auto _make_table()
{
    std::array<std::array<u128, 128>, PN> table{};
    for (int i = 0; i < PN; ++i)
    {
        table[i][0] = primes32[i];
        for (int j = 1; j < 128; ++j)
        {
            table[i][j] = table[i][j - 1] * primes32[i];
        }
    }
    return table;
}

inline constexpr std::array<std::array<u128, 128>, PN> _pow_table = _make_table();

inline fr136_t operator*(const fr136_t &a, const fr136_t &b)
{
    fr136_t c{};
    __m256i a_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(a.data()));
    __m256i b_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(b.data()));
    __m256i c_ = _mm256_add_epi8(a_, b_);
    _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(c.data()), c_);
    return c;
}

inline fr136_t operator/(const fr136_t &a, const fr136_t &b)
{
    fr136_t c{};
    __m256i a_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(a.data()));
    __m256i b_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(b.data()));
    __m256i c_ = _mm256_sub_epi8(a_, b_);
    _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(c.data()), c_);
    return c;
}

inline fr136_t gcd(const fr136_t &a, const fr136_t &b)
{
    fr136_t c{};
    __m256i a_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(a.data()));
    __m256i b_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(b.data()));
    __m256i c_ = _mm256_min_epi8(a_, b_);
    _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(c.data()), c_);
    return c;
}

inline fr136_t lcm(const fr136_t &a, const fr136_t &b)
{
    fr136_t c{};
    __m256i a_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(a.data()));
    __m256i b_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(b.data()));
    __m256i c_ = _mm256_max_epi8(a_, b_);
    _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(c.data()), c_);
    return c;
}

inline fr136_t sgcd(const fr136_t &a, const fr136_t &b)
{
    fr136_t c{};
    __m256i a_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(a.data()));
    __m256i b_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(b.data()));
    __m256i n_ = _mm256_min_epi8(a_, b_);
    __m256i nmask = _mm256_cmpgt_epi8(n_, _mm256_setzero_si256());
    __m256i d_ = _mm256_max_epi8(a_, b_);
    __m256i dmask = _mm256_cmpgt_epi8(_mm256_setzero_si256(), d_);
    n_ = _mm256_and_si256(n_, nmask);
    d_ = _mm256_and_si256(d_, dmask);
    __m256i c_ = _mm256_add_epi8(n_, d_);
    _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(c.data()), c_);
    return c;
}

inline fr136_t square(const fr136_t &x)
{
    fr136_t r{};
    __m256i x_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(x.data()));
    __m256i r_ = _mm256_add_epi8(x_, x_);
    _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(r.data()), r_);
    return r;
}

inline fr136_t inv(const fr136_t &x)
{
    fr136_t r{};
    __m256i x_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(x.data()));
    __m256i r_ = _mm256_sub_epi8(_mm256_setzero_si256(), x_);
    _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(r.data()), r_);
    return r;
}

int last_nonzero(const fr136_t &r)
{
    __m256i r_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(r.data()));
    __m256i mask = _mm256_cmpeq_epi8(r_, _mm256_setzero_si256());
    int mask_ = _mm256_movemask_epi8(mask);
    return 32 - __builtin_clz(~mask_);
}

inline std::ostream &operator<<(std::ostream &os, const fr136_t &r)
{
    bool is_first = true;
    for (int i = 0; i < last_nonzero(r); ++i)
    {
        if (r[i] != 0)
        {
            if (is_first)
            {
                is_first = false;
            }
            else
            {
                os << " * ";
            }
            os << primes32[i] << '^' << int(r[i]);
        }
    }
    if (is_first)
        os << '1';
    return os;
}

// x -> f * x
fr136_t extract(i128 &num)
{
    fr136_t f{};
    if (num >= -1 && num <= 1)
        return f;
    for (int i = 0; i < PN; ++i)
    {
        i128 r = num % primes32[i];
        while (r == 0)
        {
            ++f[i];
            num = num / primes32[i];
            r = num % primes32[i];
        }
    }
    return f;
}

void extract_to(fr136_t &f, i128 &num)
{
    if (num >= -1 && num <= 1)
        return;
    for (int i = 0; i < PN; ++i)
    {
        i128 r = num % primes32[i];
        while (r == 0)
        {
            ++f[i];
            num = num / primes32[i];
            r = num % primes32[i];
        }
    }
    return;
}

i128 numerator(const fr136_t &r)
{
    i128 n = 1;
    for (int i = 0; i < last_nonzero(r); ++i)
    {
        if (r[i] > 0)
        {
            n = n * _pow_table[i][r[i] - 1];
        }
    }
    return n;
}

i128 denominator(const fr136_t &r)
{
    i128 d = 1;
    for (int i = 0; i < last_nonzero(r); ++i)
    {
        if (r[i] < 0)
        {
            d = d * _pow_table[i][-int(r[i]) - 1];
        }
    }
}

std::pair<i128, i128> rational(const fr136_t &r)
{
    i128 n = 1, d = 1;
    for (int i = 0; i < last_nonzero(r); ++i)
    {
        if (r[i] > 0)
        {
            n = n * _pow_table[i][r[i] - 1];
        }
        else if (r[i] < 0)
        {
            d = d * _pow_table[i][-int(r[i]) - 1];
        }
    }
    return std::make_pair(n, d);
}

// x = s^2 * r
inline void simplify(fr136_t &s, fr136_t &r, const fr136_t &x)
{
    fr136_t tmp{};
    __m256i x_ = _mm256_loadu_si256(reinterpret_cast<const __m256i_u *>(x.data()));
    __m256i neg_mask_ = _mm256_cmpgt_epi8(_mm256_setzero_si256(), x_);
    __m256i px_ = _mm256_abs_epi8(x_);
    __m256i pr_ = _mm256_and_si256(px_, _mm256_set1_epi8(1));
    __m256i ps_ = _mm256_sub_epi8(px_, pr_);
    __m256i pr_mask_ = _mm256_sub_epi8(_mm256_setzero_si256(), pr_);
    pr_mask_ = _mm256_and_si256(pr_mask_, neg_mask_);
    __m256i r_ = _mm256_or_si256(pr_, pr_mask_);
    ps_ = _mm256_srli_epi16(ps_, 1);
    __m256i ns_ = _mm256_and_si256(ps_, neg_mask_);
    __m256i s_ = _mm256_sub_epi8(ps_, ns_);
    s_ = _mm256_sub_epi8(s_, ns_);
    _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(r.data()), r_);
    _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(s.data()), s_);
}

// fr136_t representation from 1 to 136
// clang-format off
constexpr fr136_t rational136[136] = {
    {}, // 1 = 1
    {1}, // 2 = 2
    {0, 1}, // 3 = 3
    {2}, // 4 = 2^2
    {0, 0, 1}, // 5 = 5
    {1, 1}, // 6 = 2 * 3
    {0, 0, 0, 1}, // 7 = 7
    {3}, // 8 = 2^3
    {0, 2}, // 9 = 3^2
    {1, 0, 1}, // 10 = 2 * 5
    {0, 0, 0, 0, 1}, // 11 = 11
    {2, 1}, // 12 = 2^2 * 3
    {0, 0, 0, 0, 0, 1}, // 13 = 13
    {1, 0, 0, 1}, // 14 = 2 * 7
    {0, 1, 1}, // 15 = 3 * 5
    {4}, // 16 = 2^4
    {0, 0, 0, 0, 0, 0, 1}, // 17 = 17
    {1, 2}, // 18 = 2 * 3^2
    {0, 0, 0, 0, 0, 0, 0, 1}, // 19 = 19
    {2, 0, 1}, // 20 = 2^2 * 5
    {0, 1, 0, 1}, // 21 = 3 * 7
    {1, 0, 0, 0, 1}, // 22 = 2 * 11
    {0, 0, 0, 0, 0, 0, 0, 0, 1}, // 23 = 23
    {3, 1}, // 24 = 2^3 * 3
    {0, 0, 2}, // 25 = 5^2
    {1, 0, 0, 0, 0, 1}, // 26 = 2 * 13
    {0, 3}, // 27 = 3^3
    {2, 0, 0, 1}, // 28 = 2^2 * 7
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 29 = 29
    {1, 1, 1}, // 30 = 2 * 3 * 5
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 31 = 31
    {5}, // 32 = 2^5
    {0, 1, 0, 0, 1}, // 33 = 3 * 11
    {1, 0, 0, 0, 0, 0, 1}, // 34 = 2 * 17
    {0, 0, 1, 1}, // 35 = 5 * 7
    {2, 2}, // 36 = 2^2 * 3^2
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 37 = 37
    {1, 0, 0, 0, 0, 0, 0, 1}, // 38 = 2 * 19
    {0, 1, 0, 0, 0, 1}, // 39 = 3 * 13
    {3, 0, 1}, // 40 = 2^3 * 5
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 41 = 41
    {1, 1, 0, 1}, // 42 = 2 * 3 * 7
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 43 = 43
    {2, 0, 0, 0, 1}, // 44 = 2^2 * 11
    {0, 2, 1}, // 45 = 3^2 * 5
    {1, 0, 0, 0, 0, 0, 0, 0, 1}, // 46 = 2 * 23
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 47 = 47
    {4, 1}, // 48 = 2^4 * 3
    {0, 0, 0, 2}, // 49 = 7^2
    {1, 0, 2}, // 50 = 2 * 5^2
    {0, 1, 0, 0, 0, 0, 1}, // 51 = 3 * 17
    {2, 0, 0, 0, 0, 1}, // 52 = 2^2 * 13
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 53 = 53
    {1, 3}, // 54 = 2 * 3^3
    {0, 0, 1, 0, 1}, // 55 = 5 * 11
    {3, 0, 0, 1}, // 56 = 2^3 * 7
    {0, 1, 0, 0, 0, 0, 0, 1}, // 57 = 3 * 19
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 58 = 2 * 29
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 59 = 59
    {2, 1, 1}, // 60 = 2^2 * 3 * 5
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 61 = 61
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 62 = 2 * 31
    {0, 2, 0, 1}, // 63 = 3^2 * 7
    {6}, // 64 = 2^6
    {0, 0, 1, 0, 0, 1}, // 65 = 5 * 13
    {1, 1, 0, 0, 1}, // 66 = 2 * 3 * 11
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 67 = 67
    {2, 0, 0, 0, 0, 0, 1}, // 68 = 2^2 * 17
    {0, 1, 0, 0, 0, 0, 0, 0, 1}, // 69 = 3 * 23
    {1, 0, 1, 1}, // 70 = 2 * 5 * 7
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 71 = 71
    {3, 2}, // 72 = 2^3 * 3^2
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 73 = 73
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 74 = 2 * 37
    {0, 1, 2}, // 75 = 3 * 5^2
    {2, 0, 0, 0, 0, 0, 0, 1}, // 76 = 2^2 * 19
    {0, 0, 0, 1, 1}, // 77 = 7 * 11
    {1, 1, 0, 0, 0, 1}, // 78 = 2 * 3 * 13
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 79 = 79
    {4, 0, 1}, // 80 = 2^4 * 5
    {0, 4}, // 81 = 3^4
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 82 = 2 * 41
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 83 = 83
    {2, 1, 0, 1}, // 84 = 2^2 * 3 * 7
    {0, 0, 1, 0, 0, 0, 1}, // 85 = 5 * 17
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 86 = 2 * 43
    {0, 1, 0, 0, 0, 0, 0, 0, 0, 1}, // 87 = 3 * 29
    {3, 0, 0, 0, 1}, // 88 = 2^3 * 11
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 89 = 89
    {1, 2, 1}, // 90 = 2 * 3^2 * 5
    {0, 0, 0, 1, 0, 1}, // 91 = 7 * 13
    {2, 0, 0, 0, 0, 0, 0, 0, 1}, // 92 = 2^2 * 23
    {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 93 = 3 * 31
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 94 = 2 * 47
    {0, 0, 1, 0, 0, 0, 0, 1}, // 95 = 5 * 19
    {5, 1}, // 96 = 2^5 * 3
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 97 = 97
    {1, 0, 0, 2}, // 98 = 2 * 7^2
    {0, 2, 0, 0, 1}, // 99 = 3^2 * 11
    {2, 0, 2}, // 100 = 2^2 * 5^2
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 101 = 101
    {1, 1, 0, 0, 0, 0, 1}, // 102 = 2 * 3 * 17
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 103 = 103
    {3, 0, 0, 0, 0, 1}, // 104 = 2^3 * 13
    {0, 1, 1, 1}, // 105 = 3 * 5 * 7
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 106 = 2 * 53
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 107 = 107
    {2, 3}, // 108 = 2^2 * 3^3
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 109 = 109
    {1, 0, 1, 0, 1}, // 110 = 2 * 5 * 11
    {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 111 = 3 * 37
    {4, 0, 0, 1}, // 112 = 2^4 * 7
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 113 = 113
    {1, 1, 0, 0, 0, 0, 0, 1}, // 114 = 2 * 3 * 19
    {0, 0, 1, 0, 0, 0, 0, 0, 1}, // 115 = 5 * 23
    {2, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 116 = 2^2 * 29
    {0, 2, 0, 0, 0, 1}, // 117 = 3^2 * 13
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 118 = 2 * 59
    {0, 0, 0, 1, 0, 0, 1}, // 119 = 7 * 17
    {3, 1, 1}, // 120 = 2^3 * 3 * 5
    {0, 0, 0, 0, 2}, // 121 = 11^2
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 122 = 2 * 61
    {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 123 = 3 * 41
    {2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 124 = 2^2 * 31
    {0, 0, 3}, // 125 = 5^3
    {1, 2, 0, 1}, // 126 = 2 * 3^2 * 7
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 127 = 127
    {7}, // 128 = 2^7
    {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 129 = 3 * 43
    {1, 0, 1, 0, 0, 1}, // 130 = 2 * 5 * 13
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 131 = 131
    {2, 1, 0, 0, 1}, // 132 = 2^2 * 3 * 11
    {0, 0, 0, 1, 0, 0, 0, 1}, // 133 = 7 * 19
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 134 = 2 * 67
    {0, 3, 1}, // 135 = 3^3 * 5
    {3, 0, 0, 0, 0, 0, 1}, // 136 = 2^3 * 17
};
// clang-format on

fr136_t calc_fixed_binomial(int n, int k)
{
    if (unsigned(n) > 136 || unsigned(k) > unsigned(n))
    {
        throw std::runtime_error("n or k out of range");
    }
    fr136_t r{};
    k = std::min(k, n - k);
    for (int i = 1; i <= k; ++i)
    {
        r = r * rational136[n - i] / rational136[i - 1];
    }
    return r;
}

class Wigner
{
  public:
    Wigner()
    {
        m_data.resize(_binomial_data_size(136), fr136_t{});
        for (int n = 0; n <= 136; ++n)
        {
            for (int k = 0; k <= n / 2; ++k)
            {
                m_data[_binomial_index(n, k)] = calc_fixed_binomial(n, k);
            }
        }
    }
    // judge if a number is a odd number
    static bool isodd(int x) { return x % 2 != 0; }
    // judge if a number is a even number
    static bool iseven(int x) { return x % 2 == 0; }
    // judge if two number are same odd or same even
    static bool is_same_parity(int x, int y) { return iseven(x ^ y); }
    // return (-1)^n
    static int iphase(int x) { return iseven(x) - isodd(x); }
    // check if m-quantum number if one of the components of a the j-quantum number
    static bool check_jm(int dj, int dm) { return is_same_parity(dj, dm) && (std::abs(dm) <= dj); }
    // judge if three angular momentum can couple
    static bool check_couple(int dj1, int dj2, int dj3)
    {
        return dj1 >= 0 && dj2 >= 0 && is_same_parity(dj1 + dj2, dj3) && (dj3 <= (dj1 + dj2)) &&
               (dj3 >= std::abs(dj1 - dj2));
    }
    // only works for positive n
    static double quick_pow(double x, int n)
    {
        double ans = 1;
        while (n)
        {
            if (n & 1)
                ans = ans * x;
            n = n >> 1;
            x = x * x;
        }
        return ans;
    }

    fr136_t unsafe_binomial(int n, int k) const
    {
        k = std::min(k, n - k);
        return m_data[_binomial_index(n, k)];
    }

    double CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3) const
    {
        if (!(check_jm(dj1, dm1) && check_jm(dj2, dm2) && check_jm(dj3, dm3)))
            return 0;
        if (!check_couple(dj1, dj2, dj3))
            return 0;
        if (dm1 + dm2 != dm3)
            return 0;
        const int J = (dj1 + dj2 + dj3) / 2;
        const int jm1 = J - dj1;
        const int jm2 = J - dj2;
        const int jm3 = J - dj3;
        const int j1mm1 = (dj1 - dm1) / 2;
        const int j2mm2 = (dj2 - dm2) / 2;
        const int j3mm3 = (dj3 - dm3) / 2;
        const int j2pm2 = (dj2 + dm2) / 2;
        const auto A = (unsafe_binomial(dj1, jm2) * unsafe_binomial(dj2, jm3)) /
                       (unsafe_binomial(J + 1, jm3) * unsafe_binomial(dj1, j1mm1) * unsafe_binomial(dj2, j2mm2) *
                        unsafe_binomial(dj3, j3mm3));
        const int low = std::max(0, std::max(j1mm1 - jm2, j2pm2 - jm1));
        const int high = std::min(jm3, std::min(j1mm1, j2pm2));
        std::vector<fr136_t> Bs;
        Bs.reserve(high - low + 1);
        for (auto z = low; z <= high; ++z)
        {
            Bs.push_back(unsafe_binomial(jm3, z) * unsafe_binomial(jm2, j1mm1 - z) * unsafe_binomial(jm1, j2pm2 - z));
        }
        i128 B = 0;
        fr136_t cfB = stagger_sum(B, Bs);
        extract_to(cfB, B);
        fr136_t s, r;
        simplify(s, r, A);
        s = s * cfB;

        fr136_t g = sgcd(s, inv(r));
        s = s / g;
        r = r * square(g);

        auto [sn, sd] = rational(s);
        auto [rn, rd] = rational(r);
        return iphase(low) * as_double(sn) * as_double(B) / as_double(sd) * std::sqrt(as_double(rn) / as_double(rd));
    }

    double f3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
    {
        return iphase((dj3 + dm3) / 2 + dj1) / std::sqrt(dj3 + 1) * CG(dj1, dj2, dj3, -dm1, -dm2, dm3);
    }

  private:
    static fr136_t stagger_sum(i128 &sum, const std::vector<fr136_t> &xs)
    {
        fr136_t cf;
        cf.fill(std::numeric_limits<int8_t>::max());
        for (auto x : xs)
        {
            cf = gcd(cf, x);
        }
        int sign = 1;
        for (auto x : xs)
        {
            auto t = x / cf;
            i128 num = numerator(t);
            if (sign > 0)
                sum = sum + num;
            else
                sum = sum - num;
            sign = -sign;
        }
        return cf;
    }

  private:
    static std::size_t _binomial_data_size(int n)
    {
        std::size_t x = n / 2 + 1;
        return x * (x + isodd(n));
    }
    static std::size_t _binomial_index(int n, int k)
    {
        std::size_t x = n / 2 + 1;
        return x * (x - iseven(n)) + k;
    }
    std::vector<fr136_t> m_data;
};

inline Wigner _wigner{};

double CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3) { return _wigner.CG(dj1, dj2, dj3, dm1, dm2, dm3); }

double wigner_3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    return _wigner.f3j(dj1, dj2, dj3, dm1, dm2, dm3);
}

} // namespace fr136

#endif // _FIXED_RATIONAL_H_