#!/bin/bash
# Reproduces the Section 5 numeric table of Paper 9 (referee r1 item M3). Usage: run_table.sh PYTHON
# (a python with numpy, scipy, mpmath). Run from this directory. One output file per (N, K, eps) in out/.
# Waits (max 10 h) for the 1-minute load to drop below 40, then runs at nice 19 on 3 threads.
# Items are ordered cheap first; no new item starts after the wall budget (WALL_BUDGET s, default 7200 = 6 core-h at 3 cores).
PY=${1:?python}
cd "$(dirname "$0")" || exit 2
mkdir -p out
T0=$(date +%s)
while :; do
  L=$(cut -d' ' -f1 /proc/loadavg)
  [ "${L%.*}" -lt 40 ] && break
  if [ $(( $(date +%s) - T0 )) -ge 36000 ]; then echo "ABORT: load stayed >= 40 for 10 h (last $L)"; exit 3; fi
  sleep 60
done
echo "start $(date -u +%FT%TZ) load $L"
export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3
S1=$(date +%s)
BUDGET=${WALL_BUDGET:-7200}
# spec: N K n
for spec in "1,2 20 10" "1,2 22 10" "1,2 26 10" \
            "3 8 10" "3 10 10" "4 5 8" "4 6 8" \
            "3 12 10" "4 7 8" "3 14 10" "4 8 8"; do
  set -- $spec
  for eps in 0.1 0.01; do
    f=out/N$(echo $1|tr , _)_K$2_n$3_eps$eps.out
    [ -s "$f.done" ] && continue
    if [ $(( $(date +%s) - S1 )) -ge $BUDGET ]; then echo "SKIP N=$1 K=$2 eps=$eps (wall budget used)" | tee -a out/skipped.txt; continue; fi
    a=$(date +%s)
    timeout 5400 nice -n 19 "$PY" -I shear_clock_quantum.py "$1" "$2" "$3" "$eps" 0.01,0.1,1 > "$f" 2>&1
    rc=$?
    echo "runtime_seconds=$(( $(date +%s) - a )) rc=$rc threads=3" >> "$f"
    touch "$f.done"; echo "done $f rc=$rc"
  done
done
echo "end $(date -u +%FT%TZ)"
