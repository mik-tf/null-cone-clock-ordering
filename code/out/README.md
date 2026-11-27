# Outputs of the Paper 9 numeric-table rerun

Everything is produced by `../run_table.sh PYTHON` (a Python with numpy, scipy, mpmath), run from `code/`. It calls the unchanged `../shear_clock_quantum.py` once per file, as

    python3 -I shear_clock_quantum.py N K n EPS 0.01,0.1,1

with threads limited to 3 (`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=3`) at `nice -n 19`. The last line of each file holds `runtime_seconds`, the return code and the thread count.

| file | command (N K n EPS) |
|---|---|
| N1_2_K{20,22,26}_n10_eps{0.1,0.01}.out | `1,2 K 10 EPS` |
| N3_K{8,10,12}_n10_eps{0.1,0.01}.out | `3 K 10 EPS` |
| N4_K{5,6}_n8_eps{0.1,0.01}.out | `4 K 8 EPS` |
| N4_K7_n8_eps0.1.out | `4 7 8 0.1` (killed by the 5400 s cap, rc=124; kappa = 0.01 and 0.1 only) |
| skipped.txt | items not started because the wall budget (7200 s) was used |

Software: numpy 2.5.5, scipy 1.18.1, mpmath 1.3.0. The script waits for a 1-minute load under 40 before starting.
