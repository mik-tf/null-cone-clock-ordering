"""Shear-clock quantum, operator inequality (2026-11-26): numbers for shear-clock-quantum-proof.md.

Theorem B (tangent family): for every function c(t), e^x >= e^c (1 + x - c) gives the operator bound
  tau_hat >= A0(c) + sum_k w_k(c) Delta_k^2 >= A0(c) + ell_hat^2 / W(c),   W = sum_k 1/w_k,
the second step by the weighted operator Cauchy-Schwarz for the noncommuting |Delta_k|. With c = X* of the discrete
classical optimum of length x the quadratic G_x(l) = A0 + l^2/W touches the discrete classical bound F^(M) at x.
Sections:
  1. discrete bound F^(M) vs continuum F_eps, envelope check G_x(x) = F^(M)(x);
  2. low window [0, l*] where the vacuum tangent G_0 >= F_eps (operator inequality on that spectral subspace);
  3. two-point obstruction (no tangent certifies a superposition of a small-l and a large-l eigenvector);
  4. operator checks at N = 1: lam_min(tau_hat - G_x(ell_hat)) (should be >= 0) and lam_min(tau_hat - F^(M)(ell_hat));
  5. normal-ordered clock: vacuum violates even the expectation inequality for small kappa (counterexample).
Run: python3 -I shear_clock_quantum_proof.py   (light: N = 1, 2 only)
"""
import sys, os, numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import shear_clock_quantum as Q

def quad_nodes(M, h, nq):
    xg, wg = np.polynomial.legendre.leggauss(nq)
    return [(k*h + (xg + 1)/2*h, wg*h/2) for k in range(M)]

def tau_disc(D, eps, M, nq=40):
    """Classical clock of increments D (same quadrature as the operator version)."""
    Lam = np.log(1/eps); h = Lam/M; tau = 0.0; Xa = 0.0
    for k, (t, w) in enumerate(quad_nodes(M, h, nq)):
        X = Xa + D[k]**2*(t - k*h)/(4*h*h)
        tau += eps*np.sum(w*np.exp(t + X)); Xa += D[k]**2/(4*h)
    return tau

def tangent(D, eps, M, nq=40):
    """A0, w_k, W of the linearisation at c = X[D]."""
    Lam = np.log(1/eps); h = Lam/M; Xa = 0.0; A0 = 0.0; w = np.zeros(M); nodes = quad_nodes(M, h, nq)
    for k, (t, wq) in enumerate(nodes):
        c = Xa + D[k]**2*(t - k*h)/(4*h*h); g = eps*wq*np.exp(t + c)
        A0 += np.sum(g*(1 - c))
        w[k] += np.sum(g*(t - k*h))/(4*h*h)       # coefficient of Delta_k^2 inside its own cell
        for j in range(k):                          # coefficient h/(4h^2) of Delta_j^2 for later cells
            w[j] += np.sum(g)*h/(4*h*h)
        Xa += D[k]**2/(4*h)
    return A0, w, np.sum(1/w)

def FM(x, eps, M, nq=40):
    """Discrete classical bound: min tau over D >= 0, sum D = x (convex program; log tau minimised)."""
    if x == 0: return tau_disc(np.zeros(M), eps, M, nq), np.zeros(M)
    best = None
    for start in (np.ones(M)/M, np.linspace(0.5, 1.5, M)/np.sum(np.linspace(0.5, 1.5, M))):
        f = lambda u: np.log(tau_disc(x*np.abs(u)/np.sum(np.abs(u)), eps, M, nq))
        r = minimize(f, start, method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-14, maxiter=40000))
        r = minimize(f, np.abs(r.x), method='BFGS', options=dict(gtol=1e-12))
        if best is None or r.fun < best.fun: best = r
    D = x*np.abs(best.x)/np.sum(np.abs(best.x))
    return np.exp(best.fun), D

def G(l, tg):
    A0, w, W = tg
    return A0 + np.asarray(l)**2/W

def main():
    np.set_printoptions(precision=4)
    out = []
    P = lambda *a: (print(*a), out.append(' '.join(str(x) for x in a)))
    Fcs = {eps: Q.Fclock(eps) for eps in (0.1, 0.01)}

    P("## 1. F^(M) vs F_eps and envelope check")
    for eps in (0.1, 0.01):
        Fc = Fcs[eps]
        for M in (2, 4, 8):
            for x in (0.1, 0.5, 1.0, 2.0, 4.0):
                fm, D = FM(x, eps, M); tg = tangent(D, eps, M)
                P(f"eps={eps} M={M} x={x}: F^M={fm:.6g} F_eps={float(Fc(x)):.6g} rel.gap={(fm-float(Fc(x)))/(fm-(1-eps)):.3f} "
                  f"G_x(x)-F^M={G(x, tg)-fm:.2e} Dk*w_k/mean={np.round(D*tg[1]/np.mean(D*tg[1]), 4)}")

    P("## 2. low window [0, l*]: G_0(l) = (1-eps) + l^2/W_0 >= F_eps(l)")
    for eps in (0.1, 0.01):
        Fc = Fcs[eps]
        for M in (2, 4, 8, 16, 32):
            tg0 = tangent(np.zeros(M), eps, M)
            ls = np.linspace(1e-4, 20, 200001); d = G(ls, tg0) - Fc(ls)
            neg = np.where(d < 0)[0]; lstar = ls[neg[0]] if len(neg) else np.inf
            P(f"eps={eps} M={M}: W_0={tg0[2]:.5g} l*={lstar:.4f}")

    P("## 3. two-point obstruction (M=2, eps=0.1): best tangent margin p(G-F)(s)+(1-p)(G-F)(L)")
    eps, M = 0.1, 2; Fc = Fcs[eps]
    xs = np.concatenate([[0.0], np.geomspace(0.01, 12, 120)])
    tgs = []
    for x in xs:
        fm, D = FM(x, eps, M); tgs.append(tangent(D, eps, M))
    for s, L in ((0.2, 3.0), (0.2, 6.0), (1.0, 6.0), (0.2, 10.0)):
        for p in (0.01, 0.1, 0.5, 0.9, 0.99):
            best = max(p*(G(s, tg) - float(Fc(s))) + (1-p)*(G(L, tg) - float(Fc(L))) for tg in tgs)
            P(f"s={s} L={L} p={p}: best margin {best:.4g}  (F(L)={float(Fc(L)):.4g})")

    P("## 4. operator checks at N = 1 (K = 26, full space)")
    for eps in (0.1,):
        Fc = Fcs[eps]
        for kcap in (0.1, 1.0):
            mod = Q.model(2, eps, kcap*eps**2); tau, ell, Fl, Xacc = Q.clock_ops(mod, 26, Fc)
            idx = Q.occ_index(1, 26, 20)
            lev = np.linalg.eigvalsh(ell[np.ix_(idx, idx)])
            for x in (0.0, 0.5, 2.0):
                fm, D = FM(x, eps, 2, nq=6); tg = tangent(D, eps, 2, nq=6)
                Op = tau - Q.fun(ell, lambda l: G(l, tg))
                P(f"eps={eps} kcap={kcap}: lam_min(tau - G_x(ell)) x={x}: {np.linalg.eigvalsh(Op[np.ix_(idx, idx)])[0]:.3e}")
            grid = np.linspace(0, max(lev.max(), 1)*1.2, 61)
            fms = np.array([FM(g, eps, 2, nq=6)[0] for g in grid])
            FMf = lambda l: np.exp(np.interp(l, grid, np.log(fms)))
            Op = tau - Q.fun(ell, FMf)
            P(f"eps={eps} kcap={kcap}: lam_min(tau - F^(M)(ell)) = {np.linalg.eigvalsh(Op[np.ix_(idx, idx)])[0]:.3e}; "
              f"lam_min(tau - F_eps(ell)) = {np.linalg.eigvalsh((tau-Fl)[np.ix_(idx, idx)])[0]:.3e}")

    P("## 5. normal-ordered clock: vacuum (N = 1, 2)")
    for eps in (0.1, 0.01):
        Fc = Fcs[eps]
        for N, K in ((1, 40), (2, 14)):
            for kcap in (1.0, 0.1, 0.01, 0.001):
                M = 2*N; mod = Q.model(M, eps, kcap*eps**2); h = mod['h']
                Dl, D2 = Q.operators(mod, K)
                z2 = [D2k[0, 0].real for D2k in D2]                    # vacuum <Delta_k^2>
                xg, wg = np.polynomial.legendre.leggauss(6)
                dim = Dl[0].shape[0]; I = np.eye(dim)
                tauN = np.zeros((dim, dim), complex); Xa = np.zeros((dim, dim), complex)
                for k in range(M):
                    for x, w in zip(xg, wg):
                        s = (x + 1)/2*h
                        Xt = Xa + (D2[k] - z2[k]*I)*s/(4*h*h)
                        tauN += eps*(w*h/2)*np.exp(k*h + s)*Q.fun(Xt, np.exp)
                    Xa = Xa + (D2[k] - z2[k]*I)/(4*h)
                ell = sum(Q.fun(D, np.abs) for D in Dl)
                idx = Q.occ_index(N, K, K-2 if N == 1 else 10)
                Fl = Q.fun(ell, Fc)
                e_tau, e_ell, e_F = tauN[0, 0].real, ell[0, 0].real, Fl[0, 0].real
                lam = np.linalg.eigvalsh((tauN - Fl)[np.ix_(idx, idx)])[0]
                P(f"eps={eps} N={N} kcap={kcap}: <tauN>-(1-eps)={e_tau-(1-eps):.3e}  F(<ell>)-(1-eps)={float(Fc(e_ell))-(1-eps):.3e}  "
                  f"<F(ell)>-(1-eps)={e_F-(1-eps):.3e}  lam_min(tauN-F(ell))={lam:.3e}")
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'shear_clock_quantum_proof.out'), 'w') as f:
        f.write('\n'.join(out) + '\n')

if __name__ == '__main__':
    main()
