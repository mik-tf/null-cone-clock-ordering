"""Cutoff scan of the M=2 (N=1) weak-value margin for the coherent state a=(-1,0), eps=0.1, kcap=0.01.
Reuses q1_weak_value_clock.py (read only). Prints, for each K, the margin at every level with p>=1e-3 and the minimum.
Run: nice -n 19 python3 -I q1_k_scan.py   (from code/)"""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "2"
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import shear_clock_quantum as sc
import q1_weak_value_clock as q
eps, kc, Kb = 0.1, 0.01, 160
spec = [("disp", 0, (-1+0j)/np.sqrt(2))]
full = q.gaussian_state(1, Kb, spec)
Fc = sc.Fclock(eps); mod = sc.model(2, eps, kc*eps**2)
print(f"eps={eps} kcap={kc} state coh a=(-1,0) Kb={Kb}")
for K in (24, 30, 34, 40, 48, 52, 60, 70):
    psi, loss = q.truncate(full, 1, Kb, K)
    tau, ell, Fl, X = sc.clock_ops(mod, K, Fc)
    lam, V = np.linalg.eigh(ell); Fe = Fc(lam)
    c = V.conj().T @ psi; p = np.abs(c)**2; tc = V.conj().T @ (tau @ psi)
    marg = (tc/c).real - Fe
    ok = p >= 1e-3
    j = int(np.argmin(np.where(ok, marg, np.inf)))
    avg = (psi.conj() @ tau @ psi).real - float(np.sum(p*Fe))
    print(f"K={K} trunc_loss={loss:.2e} min margin(p>=1e-3)={marg[j]:.6e} at ell={lam[j]:.5f} p={p[j]:.3e} | <tau>-<F(ell)>={avg:.6e} | levels p>=1e-3: " +
          " ".join(f"({lam[i]:.4f},{marg[i]:+.2e})" for i in np.where(ok)[0]), flush=True)
