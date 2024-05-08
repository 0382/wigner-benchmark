# wigner-benchmark

Benchmark my [WignerSymbol](https://github.com/0382/WignerSymbol) library and `GNU-GSL` library's `gsl_sf_coupling_?j` functions.

The error is evaluated using [wigxjpf](https://fy.chalmers.se/subatom/wigxjpf/) library, which can evaluate  Wigner 3j, 6j and 9j symbols accurately using prime factorisation and multi-word integer arithmetic.

## Error

### 3j

![3j-error](data/bench_3j.svg)

### 6j

![6j-error](data/bench_6j.svg)

The `gsl_sf_coupling_6j` will throw error because `gamma` function overflow at `Jmax == 86`.

## Performance

See `src/bench.cpp` for details.

### 3j

Calculate all possible Wigner 3j symbol up to Jmax.

|        Jmax        |  10   |  15   |  20   |  25   |  30   |
| :----------------: | :---: | :---: | :---: | :---: | :---: |
|   wigner_3j(my)    | <1ms  |  6ms  | 20ms  | 64ms  | 144ms |
| gsl_sf_coupling_3j | 14ms  | 103ms | 446ms | 1.52s | 4.26s |
|  wig6jj(wigxjpf)   |  8ms  | 77ms  | 422ms | 1.55s | 4.72s |

### 6j

Calculate all possible Wigner 6j symbol up to Jmax. 

|        Jmax        |  10   |  15   |  20   |  25   |  30   |
| :----------------: | :---: | :---: | :---: | :---: | :---: |
|   wigner_6j(my)    |  7ms  | 75ms  | 327ms | 1.48s | 3.97s |
| gsl_sf_coupling_6j | 41ms  | 455ms | 2.43s | 9.84s | 32.4s |
|  wig6jj(wigxjpf)   | 82ms  | 1.01s | 6.90s | 32.1s | 118s  |