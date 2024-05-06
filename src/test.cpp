#include "WignerSymbol.hpp"
#include "pf_wigner.h"
#include <gsl/gsl_sf_coupling.h>

using namespace util;

struct bench_stat_t
{
    std::size_t count;
    double max_diff;
    double sum_diff;
    double sum_diff2;
    double mean_diff;
    double std_dev;
};

void bench_3j(bench_stat_t *report, int N)
{
    std::fill(report, report + N, bench_stat_t());
    wigner_init(2 * N, "Jmax", 3);
    for (int dj1 = 1; dj1 <= N; ++dj1)
    {
        std::cout << "dj1: " << dj1 << std::endl;
        bench_stat_t &r = report[dj1 - 1];
        for (int dj2 = 1; dj2 <= dj1; ++dj2)
        {
            for (int dj3 = (dj1 + dj2) & 1; dj3 <= dj2; dj3 += 2)
            {
                for (int dm1 = -dj1; dm1 <= 0; dm1 += 2)
                {
                    for (int dm2 = -dj2; dm2 <= dj2; dm2 += 2)
                    {
                        const int dm3 = -dm1 - dm2;
                        sqrt_rational_t x = pf_3j(dj1, dj2, dj3, dm1, dm2, dm3);
                        double pf = sqrt_rational_to_double(&x);
                        sqrt_rational_clear(&x);
                        double w = wigner_3j(dj1, dj2, dj3, dm1, dm2, dm3);
                        double diff = std::abs(pf - w);
                        r.count += 1;
                        r.max_diff = std::max(r.max_diff, diff);
                        r.sum_diff += diff;
                        r.sum_diff2 += diff * diff;
                    }
                }
            }
        }
    }
}

void bench_3j_gsl(bench_stat_t *report, int N)
{
    std::fill(report, report + N, bench_stat_t());
    wigner_init(2 * N, "Jmax", 3);
    for (int dj1 = 1; dj1 <= N; ++dj1)
    {
        std::cout << "dj1: " << dj1 << std::endl;
        bench_stat_t &r = report[dj1 - 1];
        for (int dj2 = 1; dj2 <= dj1; ++dj2)
        {
            for (int dj3 = (dj1 + dj2) & 1; dj3 <= dj2; dj3 += 2)
            {
                for (int dm1 = -dj1; dm1 <= 0; dm1 += 2)
                {
                    for (int dm2 = -dj2; dm2 <= dj2; dm2 += 2)
                    {
                        const int dm3 = -dm1 - dm2;
                        sqrt_rational_t x = pf_3j(dj1, dj2, dj3, dm1, dm2, dm3);
                        double pf = sqrt_rational_to_double(&x);
                        sqrt_rational_clear(&x);
                        double w = gsl_sf_coupling_3j(dj1, dj2, dj3, dm1, dm2, dm3);
                        double diff = std::abs(pf - w);
                        r.count += 1;
                        r.max_diff = std::max(r.max_diff, diff);
                        r.sum_diff += diff;
                        r.sum_diff2 += diff * diff;
                    }
                }
            }
        }
    }
}

int main()
{
    pf_init_integers(1000);
    constexpr int N = 100;
    bench_stat_t reports[N];
    bench_3j_gsl(reports, N);
    std::printf("%3s %18s %18s %18s\n", "i", "mean", "std_dev", "max");
    for (int i = 0; i < N; ++i)
    {
        reports[i].mean_diff = reports[i].sum_diff / reports[i].count;
        reports[i].std_dev = sqrt(reports[i].sum_diff2 / reports[i].count);
        std::printf("%3d %18.6g %18.6g %18.6g\n", i, reports[i].mean_diff, reports[i].std_dev, reports[i].max_diff);
    }
    pf_free_integers();
    return 0;
}