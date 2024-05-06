#include "WignerSymbol.hpp"
#include "pf_wigner.h"

int main()
{
    pf_init_integers(1000);
    sqrt_rational_t x = pf_CG(24, 34, 56, 12, -32, -20);
    print_sqrt_rational(&x);
    sqrt_rational_clear(&x);
    x = pf_CG(2 * 100, 2 * 300, 2 * 285, 2 * 2, 2 * -2, 0);
    double rx = sqrt_rational_to_double(&x);
    print_sqrt_rational(&x);
    sqrt_rational_clear(&x);
    pf_free_integers();
    printf("%f\n", rx);
    return 0;
}