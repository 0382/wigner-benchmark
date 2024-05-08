#include "WignerSymbol.hpp"
#include "wigxjpf.h"
#include <chrono>
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

using pfun6_t = double (*)(int, int, int, int, int, int);
using pfun9_t = double (*)(int, int, int, int, int, int, int, int, int);

void error_3j(int N, pfun6_t func, const char *name)
{
    std::cout << name << ":\n";
    auto t1 = std::chrono::high_resolution_clock::now();
    bench_stat_t *stats = new bench_stat_t[N];
    std::fill(stats, stats + N, bench_stat_t());
    wigner_init(N, "Jmax", 3);
    for (int dj1 = 1; dj1 <= N; ++dj1)
    {
        std::cout << dj1 << '\r' << std::flush;
        bench_stat_t &r = stats[dj1 - 1];
        double c1 = 0, c2 = 0;
        for (int dj2 = 0; dj2 <= dj1; ++dj2)
        {
            for (int dj3 = std::abs(dj1 - dj2); dj3 <= dj2; dj3 += 2)
            {
                for (int dm1 = -dj1; dm1 <= 0; dm1 += 2)
                {
                    for (int dm2 = -dj2; dm2 <= dj2; dm2 += 2)
                    {
                        const int dm3 = -dm1 - dm2;
                        double exact = wig3jj(dj1, dj2, dj3, dm1, dm2, dm3);
                        double w3j = func(dj1, dj2, dj3, dm1, dm2, dm3);
                        double diff = std::abs(exact - w3j);
                        r.count += 1;
                        r.max_diff = std::max(r.max_diff, diff);
                        double s1 = diff - c1;
                        double t1 = r.sum_diff + s1;
                        c1 = (t1 - r.sum_diff) - s1;
                        r.sum_diff = t1;
                        double s2 = diff * diff - c2;
                        double t2 = r.sum_diff2 + s2;
                        c2 = (t2 - r.sum_diff2) - s2;
                        r.sum_diff2 = t2;
                    }
                }
            }
        }
    }
    std::printf("\n%3s %10s %18s %18s %18s\n", "i", "count", "mean", "std_dev", "max");
    for (int i = 0; i < N; ++i)
    {
        bench_stat_t &r = stats[i];
        r.mean_diff = r.sum_diff / r.count;
        r.std_dev = sqrt(r.sum_diff2 / r.count);
        std::printf("%3d %10lu %18.6g %18.6g %18.6g\n", i + 1, r.count, r.mean_diff, r.std_dev, r.max_diff);
    }
    delete[] stats;
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Time: " << duration << " ms\n";
}

void error_6j(int N, pfun6_t func, const char *name)
{
    std::cout << name << ":\n";
    auto t1 = std::chrono::high_resolution_clock::now();
    bench_stat_t *stats = new bench_stat_t[N];
    std::fill(stats, stats + N, bench_stat_t());
    wigner_init(N, "Jmax", 6);
    for (int dj1 = 1; dj1 <= N; ++dj1)
    {
        std::cout << dj1 << '\r' << std::flush;
        bench_stat_t &r = stats[dj1 - 1];
        double c1 = 0, c2 = 0;
        for (int dj2 = 0; dj2 <= dj1; ++dj2)
        {
            for (int dj3 = std::abs(dj1 - dj2); dj3 <= dj2; dj3 += 2)
            {
                for (int dj4 = 0; dj4 <= dj1; ++dj4)
                {
                    for (int dj5 = std::abs(dj3 - dj4); dj5 <= std::min(dj1, dj3 + dj4); dj5 += 2)
                    {
                        const int dj6_min = std::max(std::abs(dj1 - dj5), std::abs(dj2 - dj4));
                        const int dj6_max = std::min(dj1, std::abs(dj2 + dj4));
                        for (int dj6 = dj6_min; dj6 <= dj6_max; dj6 += 2)
                        {
                            double exact = wig6jj(dj1, dj2, dj3, dj4, dj5, dj6);
                            double w6j = func(dj1, dj2, dj3, dj4, dj5, dj6);
                            double diff = std::abs(exact - w6j);
                            r.count += 1;
                            r.max_diff = std::max(r.max_diff, diff);
                            double s1 = diff - c1;
                            double t1 = r.sum_diff + s1;
                            c1 = (t1 - r.sum_diff) - s1;
                            r.sum_diff = t1;
                            double s2 = diff * diff - c2;
                            double t2 = r.sum_diff2 + s2;
                            c2 = (t2 - r.sum_diff2) - s2;
                            r.sum_diff2 = t2;
                        }
                    }
                }
            }
        }
    }
    std::printf("\n%3s %10s %18s %18s %18s\n", "i", "count", "mean", "std_dev", "max");
    for (int i = 0; i < N; ++i)
    {
        bench_stat_t &r = stats[i];
        r.mean_diff = r.sum_diff / r.count;
        r.std_dev = sqrt(r.sum_diff2 / r.count);
        std::printf("%3d %10lu %18.6g %18.6g %18.6g\n", i + 1, r.count, r.mean_diff, r.std_dev, r.max_diff);
    }
    delete[] stats;
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Time: " << duration << " ms\n";
}

void error_9j(int N, pfun9_t func, const char *name)
{
    std::cout << name << ":\n";
    auto t1 = std::chrono::high_resolution_clock::now();
    bench_stat_t *stats = new bench_stat_t[N];
    std::fill(stats, stats + N, bench_stat_t());
    wigner_init(N, "Jmax", 6);
    for (int dj1 = 1; dj1 <= N; ++dj1)
    {
        std::cout << dj1 << '\r' << std::flush;
        bench_stat_t &r = stats[dj1 - 1];
        double c1 = 0, c2 = 0;
        for (int dj2 = 0; dj2 <= dj1; ++dj2)
        {
            for (int dj3 = 0; dj3 <= dj2; ++dj3)
            {
                for (int dj4 = 0; dj4 <= dj3; ++dj4)
                {
                    for (int dj12 = std::abs(dj1 - dj2); dj12 <= dj1; dj12 += 2)
                    {
                        for (int dj34 = std::abs(dj3 - dj4); dj34 <= std::min(dj1, dj3 + dj4); dj34 += 2)
                        {
                            for (int dj13 = std::abs(dj1 - dj3); dj13 <= dj1; dj13 += 2)
                            {
                                for (int dj24 = std::abs(dj2 - dj4); dj24 <= std::min(dj1, dj2 + dj4); dj24 += 2)
                                {
                                    int Jmin = std::max(std::abs(dj12 - dj34), std::abs(dj13 - dj24));
                                    int Jmax = std::min(dj1, std::min(dj12 + dj34, dj13 + dj24));
                                    for (int J = Jmin; J <= Jmax; J += 2)
                                    {
                                        double exact = wig9jj(dj1, dj2, dj12, dj3, dj4, dj34, dj13, dj24, J);
                                        double w9j = func(dj1, dj2, dj12, dj3, dj4, dj34, dj13, dj24, J);
                                        double diff = std::abs(exact - w9j);
                                        r.count += 1;
                                        r.max_diff = std::max(r.max_diff, diff);
                                        double s1 = diff - c1;
                                        double t1 = r.sum_diff + s1;
                                        c1 = (t1 - r.sum_diff) - s1;
                                        r.sum_diff = t1;
                                        double s2 = diff * diff - c2;
                                        double t2 = r.sum_diff2 + s2;
                                        c2 = (t2 - r.sum_diff2) - s2;
                                        r.sum_diff2 = t2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::printf("\n%3s %10s %18s %18s %18s\n", "i", "count", "mean", "std_dev", "max");
    for (int i = 0; i < N; ++i)
    {
        bench_stat_t &r = stats[i];
        r.mean_diff = r.sum_diff / r.count;
        r.std_dev = sqrt(r.sum_diff2 / r.count);
        std::printf("%3d %10lu %18.6g %18.6g %18.6g\n", i + 1, r.count, r.mean_diff, r.std_dev, r.max_diff);
    }
    delete[] stats;
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Time: " << duration << " ms\n";
}

int main()
{
    wig_table_init(2 * 200, 9);
    wig_temp_init(2 * 200);
    // error_3j(150, wigner_3j, "wigner_3j");
    // error_3j(150, gsl_sf_coupling_3j, "gsl_3j");
    // error_6j(100, wigner_6j, "wigner_6j");
    // error_6j(85, gsl_sf_coupling_6j, "gsl_6j"); // gsl gamma overflows at 86
    error_9j(50, wigner_9j, "wigner_9j");
    error_9j(50, gsl_sf_coupling_9j, "gsl_9j");
    wig_temp_free();
    wig_table_free();
    return 0;
}