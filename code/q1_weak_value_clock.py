"""Q1 first test (born-phase/scoping.md): weak-value clock bound  Re tau_w(ell) >= F_eps(ell)  for Gaussian states.
Reuses paper 9's shear_clock_quantum.py (read only, imported; no bytecode written there).
tau_w(l_j) = <l_j|tau|psi>/<l_j|psi>, |l_j> eigenvectors of the (truncated) ell_hat; p_j = |<l_j|psi>|^2.
margin_j = Re tau_w(l_j) - F_eps(l_j), kept where p_j >= pmin.  Also the averaged margin <tau> - <F(ell_hat)>.
States: pure Gaussian = displaced, squeezed (and for N=2 beam-split / two-mode squeezed) vacua of the Williamson modes,
built in a large Fock space (Kb) then truncated to K and K+4 (same state, truncation loss reported).
Fixed seed 20261126. Run: nice -n 19 python3 -I q1_weak_value_clock.py N   (N = 1 -> M=2, N = 2 -> M=4)
"""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "2"
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply
import shear_clock_quantum as sc

SEED = 20261126
PMIN, PKILL = 1e-6, 1e-3

def ann(Kb):
    return sp.diags(np.sqrt(np.arange(1, Kb)), 1, format="csr")

def mode_op(a, N, m):
    I = sp.identity(a.shape[0], format="csr"); out = sp.identity(1, format="csr")
    for i in range(N): out = sp.kron(out, a if i == m else I, format="csr")
    return out

def gaussian_state(N, Kb, spec):
    """spec: list of ops applied in order to the vacuum. ('sq',m,z) ('disp',m,alpha) ('bs',theta,phi) ('tms',g,phi)."""
    a = ann(Kb); A = [mode_op(a, N, m) for m in range(N)]
    dim = Kb**N; psi = np.zeros(dim, complex); psi[0] = 1
    for op in spec:
        if op[0] == "sq":   G = 0.5*(np.conj(op[2])*A[op[1]]@A[op[1]] - op[2]*A[op[1]].T@A[op[1]].T)
        elif op[0] == "disp": G = op[2]*A[op[1]].T - np.conj(op[2])*A[op[1]]
        elif op[0] == "bs":  G = op[1]*(np.exp(1j*op[2])*A[0].T@A[1] - np.exp(-1j*op[2])*A[1].T@A[0])
        elif op[0] == "tms": G = op[1]*(np.exp(1j*op[2])*A[0].T@A[1].T - np.exp(-1j*op[2])*A[0]@A[1])
        psi = expm_multiply(G.tocsc(), psi)
    return psi

def truncate(psi, N, Kb, K):
    t = psi.reshape((Kb,)*N)[(slice(0, K),)*N].reshape(-1)
    n2 = np.vdot(t, t).real
    return t/np.sqrt(n2), 1 - n2

def states(N, rng):
    S = []   # (label, spec, is_coherent)
    if N == 1:
        for aq in np.arange(-3, 3.01, 1.0):
            for ap in np.arange(-3, 3.01, 1.0):
                S.append((f"coh a=({aq:g},{ap:g})", [("disp", 0, (aq+1j*ap)/np.sqrt(2))], True))
                for r in (0.3, 0.6, 1.0):
                    for ph in (0, np.pi/4, np.pi/2, 3*np.pi/4):
                        S.append((f"sq r={r:g} ph={ph:.3f} a=({aq:g},{ap:g})",
                                  [("sq", 0, r*np.exp(1j*ph)), ("disp", 0, (aq+1j*ap)/np.sqrt(2))], False))
        for i in range(2000):
            S.append((f"rand{i}", [("sq", 0, rng.uniform(0, 1.0)*np.exp(2j*np.pi*rng.random())),
                                   ("disp", 0, 2.5*np.sqrt(rng.random())*np.exp(2j*np.pi*rng.random()))], False))
    else:
        for aq in np.arange(-2, 2.01, 1.0):
            for ap in np.arange(-2, 2.01, 1.0):
                for m in (0, 1):
                    S.append((f"coh mode{m} a=({aq:g},{ap:g})", [("disp", m, (aq+1j*ap)/np.sqrt(2))], True))
        for a1 in (-1.5, 0, 1.5):
            for a2 in (-1.5, 0, 1.5):
                S.append((f"coh2 a1={a1:g} a2={a2:g}", [("disp", 0, a1/np.sqrt(2)), ("disp", 1, 1j*a2/np.sqrt(2))], True))
        for r1 in (0, 0.4, 0.8):
            for r2 in (0, 0.4, 0.8):
                for ph in (0, np.pi/2):
                    for g in (0, 0.3):
                        S.append((f"sqprod r=({r1:g},{r2:g}) ph={ph:.3f} bs={g:g}",
                                  [("sq", 0, r1), ("sq", 1, r2*np.exp(1j*ph)), ("bs", g, 0.7), ("disp", 0, 0.7)], False))
        for i in range(1500):
            spec = [("sq", 0, rng.uniform(0, .7)*np.exp(2j*np.pi*rng.random())),
                    ("sq", 1, rng.uniform(0, .7)*np.exp(2j*np.pi*rng.random())),
                    ("bs", rng.uniform(0, .6), 2*np.pi*rng.random()),
                    ("tms", rng.uniform(0, .3), 2*np.pi*rng.random()),
                    ("disp", 0, 1.8*np.sqrt(rng.random())*np.exp(2j*np.pi*rng.random())),
                    ("disp", 1, 1.8*np.sqrt(rng.random())*np.exp(2j*np.pi*rng.random()))]
            S.append((f"rand{i}", spec, False))
    return S

def analyse(psi, V, lam, tauV, Fl_eig, tau_exp):
    c = V.conj().T @ psi            # <l_j|psi>
    p = np.abs(c)**2
    tc = V.conj().T @ (tauV @ psi)  # <l_j|tau|psi>
    ok = p >= PMIN
    marg = np.full(len(lam), np.inf)
    marg[ok] = (tc[ok]/c[ok]).real - Fl_eig[ok]
    j = int(np.argmin(marg)); ok3 = p >= PKILL
    m3 = np.where(ok3, marg, np.inf); j3 = int(np.argmin(m3))
    avg = float(np.sum(p*np.where(np.isfinite(marg), marg, 0)) + 0)  # restricted sum, for diagnostics only
    full = (psi.conj() @ tauV @ psi).real - float(np.sum(p*Fl_eig))
    return marg[j], lam[j], p[j], m3[j3], lam[j3], p[j3], full

def main(N):
    rng = np.random.default_rng(SEED)
    M = 2*N
    Kbase, Kb = (48, 160) if N == 1 else (20, 44)
    Ks = (Kbase, Kbase+4)
    S = states(N, rng)
    cache = {}
    for (spec_label, spec, coh) in S:
        pass
    print(f"Q1 first test, M={M} (N={N}), seed={SEED}, K={Ks}, Kb={Kb}, {len(S)} states, pmin={PMIN:g}, pkill={PKILL:g}")
    psis = {}; losses = {}
    for lab, spec, coh in S:
        full = gaussian_state(N, Kb, spec)
        for K in Ks:
            t, loss = truncate(full, N, Kb, K); psis[(lab, K)] = t; losses[(lab, K)] = loss
    for eps in (0.1, 0.01):
        Fc = sc.Fclock(eps)
        for kc in (0.01, 0.1, 1.0):
            mod = sc.model(M, eps, kc*eps**2)
            res = {}
            for K in Ks:
                tau, ell, Fl, X = sc.clock_ops(mod, K, Fc)
                lam, V = np.linalg.eigh(ell)
                Fl_eig = Fc(lam)
                best = None; bestc = None; worstavg = None; n_used = 0
                for lab, spec, coh in S:
                    if losses[(lab, K)] > 1e-4: continue
                    n_used += 1
                    r = analyse(psis[(lab, K)], V, lam, tau, Fl_eig, None)
                    if best is None or r[0] < best[0][0]: best = (r, lab)
                    if coh and (bestc is None or r[3] < bestc[0][3]): bestc = (r, lab)
                    if worstavg is None or r[6] < worstavg[0]: worstavg = (r[6], lab)
                res[K] = (best, bestc, worstavg, n_used)
            for K in Ks:
                best, bestc, wa, n = res[K]
                r, lab = best
                print(f"eps={eps:g} kcap={kc:g} K={K} used={n}/{len(S)}: min margin(p>=1e-6)={r[0]:.6e} at ell={r[1]:.5g} p={r[2]:.3g} state[{lab}]"
                      f" | min coherent margin(p>=1e-3)={bestc[0][3]:.6e} at ell={bestc[0][4]:.5g} state[{bestc[1]}]"
                      f" | min <tau>-<F(ell)>={wa[0]:.6e} state[{wa[1]}]", flush=True)
            a, b = res[Ks[0]][0][0][0], res[Ks[1]][0][0][0]
            print(f"   stability of min margin K->K+4: {a:.4e} vs {b:.4e}, rel diff {abs(a-b)/max(abs(a),1e-300):.2e}", flush=True)

if __name__ == "__main__":
    main(int(sys.argv[1]))
