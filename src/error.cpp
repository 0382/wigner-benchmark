#include "WignerSymbol.hpp"
#include "wigxjpf.h"
#include <chrono>
#include <gsl/gsl_sf_coupling.h>

using namespace util;

struct Kahan_t
{
    double sum;
    double c;
    void add(double x)
    {
        double y = x - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    double value() const { return sum; }
};

struct bench_stat_t
{
    std::size_t count;
    double max_diff;
    double max_rel_diff;
    Kahan_t sum_diff2;
    Kahan_t sum_rel_diff2;
    double std_diff;
    double std_rel_diff;
};

void print_stats(std::vector<bench_stat_t> &stats)
{
    std::printf("\n%3s %10s %18s %18s %18s %18s\n", "i", "count", "max_diff", "max_rel_diff", "std_diff",
                "std_rel_diff");
    for (int i = 0; i < stats.size(); ++i)
    {
        bench_stat_t &r = stats[i];
        r.std_diff = sqrt(r.sum_diff2.value() / r.count);
        r.std_rel_diff = sqrt(r.sum_rel_diff2.value() / r.count);
        std::printf("%3d %10lu %18.6g %18.6g %18.6g %18.6g\n", i + 1, r.count, r.max_diff, r.max_rel_diff, r.std_diff,
                    r.std_rel_diff);
    }
}

using pfun6_t = double (*)(int, int, int, int, int, int);
using pfun9_t = double (*)(int, int, int, int, int, int, int, int, int);

void error_3j(int N, pfun6_t func, const char *name)
{
    std::cout << name << ":\n";
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<bench_stat_t> stats(N);
    std::fill(stats.begin(), stats.end(), bench_stat_t());
    wigner_init(N, "Jmax", 3);
    for (int dj1 = 1; dj1 <= N; ++dj1)
    {
        std::cout << '\r' << dj1 << std::flush;
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
                        if (exact == 0)
                            continue;
                        double w3j = func(dj1, dj2, dj3, dm1, dm2, dm3);
                        double diff = std::abs(exact - w3j);
                        double rel_diff = diff / std::abs(exact);
                        r.count += 1;
                        r.max_diff = std::max(r.max_diff, diff);
                        r.max_rel_diff = std::max(r.max_rel_diff, rel_diff);
                        r.sum_diff2.add(diff * diff);
                        r.sum_rel_diff2.add(rel_diff * rel_diff);
                    }
                }
            }
        }
    }
    print_stats(stats);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Time: " << duration << " ms\n";
}

void error_6j(int N, pfun6_t func, const char *name)
{
    std::cout << name << ":\n";
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<bench_stat_t> stats(N);
    std::fill(stats.begin(), stats.end(), bench_stat_t());
    wigner_init(N, "Jmax", 6);
    for (int dj1 = 1; dj1 <= N; ++dj1)
    {
        std::cout << '\r' << dj1 << std::flush;
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
                            if (exact == 0)
                                continue;
                            double w6j = func(dj1, dj2, dj3, dj4, dj5, dj6);
                            double diff = std::abs(exact - w6j);
                            double rel_diff = diff / std::abs(exact);
                            r.count += 1;
                            r.max_diff = std::max(r.max_diff, diff);
                            r.max_rel_diff = std::max(r.max_rel_diff, rel_diff);
                            r.sum_diff2.add(diff * diff);
                            r.sum_rel_diff2.add(rel_diff * rel_diff);
                        }
                    }
                }
            }
        }
    }
    print_stats(stats);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Time: " << duration << " ms\n";
}

void error_9j(int N, pfun9_t func, const char *name)
{
    std::cout << name << ":\n";
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<bench_stat_t> stats(N);
    std::fill(stats.begin(), stats.end(), bench_stat_t());
    wigner_init(N, "Jmax", 6);
    for (int dj1 = 1; dj1 <= N; ++dj1)
    {
        std::cout << '\r' << dj1 << std::flush;
        bench_stat_t &r = stats[dj1 - 1];
        double c1 = 0, c2 = 0;
        for (int dj2 = dj1 / 2; dj2 <= dj1; ++dj2)
        {
            for (int dj3 = dj1 / 2; dj3 <= dj2; ++dj3)
            {
                for (int dj4 = 0; dj4 <= dj1; ++dj4)
                {
                    const int dj12_odd = wigner.isodd(dj1 + dj2);
                    const int dj12 = ((dj2 / 2) * 2) | dj12_odd;
                    const int dj34_odd = wigner.isodd(dj3 + dj4);
                    const int dj34 = ((std::min(dj1, dj3 + dj4) / 2) * 2) | dj34_odd;
                    const int dj13_odd = wigner.isodd(dj1 + dj3);
                    const int dj13 = ((dj3 / 2) * 2) | dj13_odd;
                    const int dj24_odd = wigner.isodd(dj2 + dj4);
                    const int dj24 = ((std::min(dj1, dj2 + dj4) / 2) * 2) | dj24_odd;
                    const int Jmin = std::max(std::abs(dj12 - dj34), std::abs(dj13 - dj24));
                    const int Jmax = std::min(dj1, std::min(dj12 + dj34, dj13 + dj24));
                    for (int J = Jmin; J <= Jmax; J += 2)
                    {
                        double exact = wig9jj(dj1, dj2, dj12, dj3, dj4, dj34, dj13, dj24, J);
                        if (exact == 0)
                            continue;
                        double w9j = func(dj1, dj2, dj12, dj3, dj4, dj34, dj13, dj24, J);
                        double diff = std::abs(exact - w9j);
                        double rel_diff = diff / std::abs(exact);
                        r.count += 1;
                        r.max_diff = std::max(r.max_diff, diff);
                        r.max_rel_diff = std::max(r.max_rel_diff, rel_diff);
                        r.sum_diff2.add(diff * diff);
                        r.sum_rel_diff2.add(rel_diff * rel_diff);
                    }
                }
            }
        }
    }
    print_stats(stats);
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
    error_9j(100, wigner_9j, "wigner_9j");
    error_9j(84, gsl_sf_coupling_9j, "gsl_9j"); // gsl gamma overflows at 85
    wig_temp_free();
    wig_table_free();
    return 0;
}