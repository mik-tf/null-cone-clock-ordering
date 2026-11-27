"""Shear-clock quantum (2026-11-26): truncated model of tau_hat - F_eps(ell_hat) in the paper's representation
(arXiv 2611.34924, Section 4.2 Gaussian/Fock representation, (4.20)-(4.22) Weyl ordering, (4.26) clock ordering,
(A.34)-(A.36) bracket of the conformal metric along one null generator).

Model (one generator, one angular harmonic, polarization z = z^1; on the line z^2 = 0 the hyperbolic metric is dz^2,
so Xi = (1/4) int s (d_s z)^2 ds and ell = int |d_t z| dt hold exactly, and the classical extremal lives on this line):
  t = ln(v/eps) in [0, Lambda], M cells of width h, nodes v_k = eps e^{k h}, increments Delta_k = z(v_k) - z(v_{k-1}),
  speed a = |Delta_k|/h on cell k (piecewise linear z in t), bracket [z(v), z(w)] = i eta sgn(w - v)/(v w)  (A.36),
  eta = 4 pi G hbar / (rho_S^2 * angular area of the harmonic), dimensionless; kappa_cap := eta/eps^2.
  Xi_hat(t) = Weyl ordered quadratic = (1/4h^2)[sum_{j<k} Delta_j^2 h + Delta_k^2 (t - t_{k-1})]   ((4.22))
  tau_hat   = eps int_0^Lambda e^t exp(Xi_hat(t)) dt, exp by spectral calculus                     ((4.26)-(4.27))
  ell_hat   = sum_k |Delta_hat_k| = Op^W(ell) (Weyl quantization of a function of one linear form = functional calculus)
  F_eps(ell_hat) by spectral calculus.  N = M/2 canonical pairs (the increments' bracket matrix is nondegenerate).
Numerics: Fock basis of the Williamson modes of sum Delta_k^2, cutoff K per mode (big space), then Rayleigh-Ritz on
total occupation <= n (small space): every Ritz value is an upper bound on the true lowest eigenvalue up to the
big-space truncation error, which is checked by varying K.
Run: python3 -I shear_clock_quantum.py [N_list] [K] [n]   e.g.  python3 -I shear_clock_quantum.py 1,2 24 10
"""
import sys, itertools, numpy as np, mpmath as mp
from scipy.linalg import schur, eigh

mp.mp.dps = 40

def dmax(tau, eps):
    """Theorem 1 of article D (closed form), as in shear_clock_deeper.py."""
    tau, eps = mp.mpf(tau), mp.mpf(eps)
    lam0, Lam = tau/eps, mp.log(1/eps)
    def beta_of(al):
        r = lam0*(1-al)/al
        return r/(1-r)
    def f(al):
        b = beta_of(al)
        return (b*mp.log(1/(1-al)) - al*mp.log(1+b))/(al+b) - Lam
    lo, hi = lam0/(1+lam0)*(1+mp.mpf(10)**-30), 1-mp.mpf(10)**-30
    flo, fhi = f(lo), f(hi)
    assert flo*fhi < 0, (tau, eps, flo, fhi)
    for _ in range(200):            # bracketing bisection (real arithmetic only)
        al = (lo+hi)/2; fa = f(al)
        if fa*flo > 0: lo, flo = al, fa
        else: hi = al
        if hi-lo < mp.mpf(10)**-34*hi: break
    al = (lo+hi)/2
    b = beta_of(al)
    return 2*mp.sqrt(al*b)/(al+b)*mp.log((1+b)/(1-al))

class Fclock:
    """F_eps(D) = inverse of D_max, tabulated; log(tau-(1-eps)) interpolated in log D; beyond the table
    log tau extrapolated linearly in D^2 (Phi ~ 4 Lambda ln tau)."""
    def __init__(self, eps, npts=400):
        self.eps = eps
        xs = []; Ds = []
        for x in np.linspace(-9, 40, npts):       # stop where the 1e-30 bracket can no longer hold the root
            try: Ds.append(float(dmax((1-eps) + mp.mpf(10)**x, eps))); xs.append(x)
            except AssertionError: break
        xs = np.array(xs); taus = [(1-eps) + 10.0**x for x in xs]
        self.lD, self.lx = np.log(Ds), xs*np.log(10.0)
        from scipy.interpolate import CubicSpline
        self.spl = CubicSpline(self.lD, self.lx)
        D0 = Ds[0]; self.D0 = D0
        self.c_lo = 10.0**xs[0]/(D0**2/(4*np.log(1/(eps*D0**2))))   # small-excess law constant
        self.D2a, self.D2b = Ds[-2]**2, Ds[-1]**2
        self.lt2a, self.lt2b = np.log(taus[-2]), np.log(taus[-1])
    def __call__(self, D):
        D = np.asarray(D, float); out = np.empty_like(D)
        lo = D <= 1e-300
        out[lo] = 1 - self.eps
        Dm = np.where(lo, 1.0, D); lD = np.log(Dm)
        sm = (~lo) & (lD < self.lD[0])
        out[sm] = (1 - self.eps) + self.c_lo*Dm[sm]**2/(4*np.log(1/(self.eps*Dm[sm]**2)))
        mid = (~lo) & (~sm) & (lD <= self.lD[-1])
        out[mid] = (1 - self.eps) + np.exp(self.spl(lD[mid]))
        hi = (~lo) & (lD > self.lD[-1])
        sl = (self.lt2b - self.lt2a)/(self.D2b - self.D2a)
        out[hi] = np.exp(self.lt2b + sl*(Dm[hi]**2 - self.D2b))
        return out

def clock_classical(Delta, eps, M):
    """Classical tau and ell for increments Delta (piecewise constant speed |Delta_k|/h)."""
    Lam = np.log(1/eps); h = Lam/M
    a2 = (np.asarray(Delta)/h)**2
    tau, X = 0.0, 0.0
    for k in range(M):
        # int_{t_{k-1}}^{t_k} e^{t + X + a2 (t - t_{k-1})/4} dt
        c = 1 + a2[k]/4; t0 = k*h
        tau += np.exp(t0 + X)*(np.expm1(c*h)/c)
        X += a2[k]*h/4
    return eps*tau, np.sum(np.abs(Delta))

def model(M, eps, eta):
    """Increment bracket B ([Delta_j, Delta_k] = i B_jk) and L with Delta = L r, r = (q1,p1,q2,p2,...), [q,p] = i,
    r the Williamson modes of sum Delta_k^2."""
    Lam = np.log(1/eps); h = Lam/M
    v = eps*np.exp(h*np.arange(M+1))
    c = np.sign(np.subtract.outer(np.arange(M+1), np.arange(M+1)).T)/np.outer(v, v)   # sgn(k-j)/(v_j v_k)
    Dm = np.zeros((M, M+1))
    for k in range(M):
        Dm[k, k+1], Dm[k, k] = 1, -1
    B = eta*Dm @ c @ Dm.T
    T, O = schur(B, output='real')
    blocks, k = [], 0
    while k < M:
        if k+1 < M and abs(T[k+1, k]) > 1e-14*np.abs(T).max():
            s = T[k, k+1]
            cols = (O[:, k], O[:, k+1]) if s > 0 else (O[:, k+1], O[:, k])
            blocks.append((abs(s), cols)); k += 2
        else:
            blocks.append((0.0, (O[:, k],))); k += 1
    L = []
    for s, cols in blocks:
        if s == 0.0:
            raise RuntimeError("degenerate bracket: central element")
        L += [np.sqrt(s)*cols[0], np.sqrt(s)*cols[1]]
    L = np.array(L).T                     # M x 2N
    assert np.allclose(L @ np.kron(np.eye(M//2), [[0, 1], [-1, 0]]) @ L.T, B, atol=1e-12*np.abs(B).max())
    return dict(M=M, eps=eps, eta=eta, h=h, v=v, B=B, L=L, s=[b[0] for b in blocks])

def single_mode(K):
    n = np.arange(1, K+2)
    a = np.diag(np.sqrt(n), 1)                         # (K+2)x(K+2)
    q = (a + a.T)/np.sqrt(2); p = (a - a.T)/(1j*np.sqrt(2))
    tr = lambda X: X[:K, :K]
    return dict(q=tr(q), p=tr(p), qq=tr(q@q), pp=tr(p@p), qp=tr((q@p + p@q)/2))

def operators(mod, K):
    """Full-space matrices of Delta_k (linear) and Op^W(Delta_k^2) (quadratic, products truncated correctly)."""
    M = mod['M']; N = M//2; L = mod['L']; sm = single_mode(K); I = np.eye(K)
    def emb(ops):   # ops: dict mode->matrix
        out = np.array([[1.0+0j]])
        for m in range(N):
            out = np.kron(out, ops.get(m, I))
        return out
    lin = [emb({m: sm['q']}) for m in range(N)], [emb({m: sm['p']}) for m in range(N)]
    R = [None]*(2*N)
    for m in range(N):
        R[2*m], R[2*m+1] = lin[0][m], lin[1][m]
    def quad(i, j):   # Weyl ordered r_i r_j
        mi, mj = i//2, j//2
        if mi != mj:
            return R[i] @ R[j]
        key = {(0, 0): 'qq', (1, 1): 'pp', (0, 1): 'qp', (1, 0): 'qp'}[(i % 2, j % 2)]
        return emb({mi: sm[key]})
    Q = {}
    for i in range(2*N):
        for j in range(i, 2*N):
            Q[(i, j)] = quad(i, j)
    Dl, D2 = [], []
    for k in range(M):
        l = L[k]
        Dl.append(sum(l[i]*R[i] for i in range(2*N)))
        D2.append(sum((1 if i == j else 2)*l[i]*l[j]*Q[(i, j)] for i in range(2*N) for j in range(i, 2*N)))
    return Dl, D2

def fun(A, f):
    w, U = eigh(A)
    return (U*f(w)) @ U.conj().T

def clock_ops(mod, K, Fc, nq=6):
    M, h, eps = mod['M'], mod['h'], mod['eps']
    Dl, D2 = operators(mod, K)
    xg, wg = np.polynomial.legendre.leggauss(nq)
    dim = Dl[0].shape[0]
    tau = np.zeros((dim, dim), complex); Xacc = np.zeros((dim, dim), complex)
    for k in range(M):
        t0 = k*h
        for x, w in zip(xg, wg):
            s = (x + 1)/2*h
            Xt = Xacc + D2[k]*s/(4*h*h)
            tau += eps*(w*h/2)*np.exp(t0 + s)*fun(Xt, np.exp)
        Xacc = Xacc + D2[k]/(4*h)
    ell = sum(fun(D, np.abs) for D in Dl)
    Fl = fun(ell, Fc)
    return tau, ell, Fl, Xacc

def occ_index(N, K, nmax):
    idx = []
    for i, occ in enumerate(itertools.product(range(K), repeat=N)):
        if sum(occ) <= nmax:
            idx.append(i)
    return np.array(idx)

def run(N, eps, kcap, K, nmax, Fc, verbose=True):
    eta = kcap*eps**2
    mod = model(2*N, eps, eta)
    tau, ell, Fl, X = clock_ops(mod, K, Fc)
    Op = tau - Fl
    idx = occ_index(N, K, nmax)
    sub = Op[np.ix_(idx, idx)]
    w, U = eigh(sub)
    vac = 0
    res = dict(N=N, eps=eps, kcap=kcap, K=K, n=nmax, lam=w[0].real,
               vac_tau=tau[vac, vac].real, vac_F=Fl[vac, vac].real, vac_ell=ell[vac, vac].real,
               F_of_vac_ell=float(Fc(ell[vac, vac].real)), vac_X=X[vac, vac].real,
               gs_tau=(U[:, 0].conj() @ tau[np.ix_(idx, idx)] @ U[:, 0]).real,
               gs_F=(U[:, 0].conj() @ Fl[np.ix_(idx, idx)] @ U[:, 0]).real, s=mod['s'])
    return res, (mod, tau, ell, Fl, idx, U[:, 0])

def coherent_check(eps, Delta_cl, kcaps, K=None):
    """Coherent states (displaced Williamson vacuum, mean increments = Delta_cl) as kappa_cap -> 0, N = 1 (M = 2).
    Returns rows (kcap, <tau>, <F(ell_hat)>, <ell>) next to the classical tau_cl, F(ell_cl), ell_cl."""
    from scipy.linalg import expm
    Fc = Fclock(eps); M = len(Delta_cl)
    tcl, lcl = clock_classical(Delta_cl, eps, M)
    rows = []
    for kc in kcaps:
        mod = model(M, eps, kc*eps**2)
        r0 = np.linalg.lstsq(mod['L'], np.asarray(Delta_cl, float), rcond=None)[0]   # (q0, p0) of the single mode
        Kc = K or int(2.5/kc) + 70            # Fock cutoff grows like 1/kappa (|alpha|^2 ~ 1/kappa)
        tau, ell, Fl, X = clock_ops(mod, Kc, Fc)
        sm = single_mode(Kc)
        Dsp = expm(-1j*(r0[1]*sm['q'] - r0[0]*sm['p']))     # displacement; truncation error negligible for |alpha|^2 << K
        psi = Dsp[:, 0]; psi = psi/np.linalg.norm(psi)
        ev = lambda A: (psi.conj() @ A @ psi).real
        rows.append((kc, ev(tau), ev(Fl), ev(ell)))
    return (tcl, float(Fc(np.array([lcl]))[0]), lcl), rows

if __name__ == "__main__":
    Ns = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else "1,2").split(",")]
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    eps_list = [float(x) for x in (sys.argv[4] if len(sys.argv) > 4 else "0.1,0.01").split(",")]
    kc_list = [float(x) for x in (sys.argv[5] if len(sys.argv) > 5 else "0.01,0.1").split(",")]
    for eps in eps_list:
        Fc = Fclock(eps)
        for kc in kc_list:
            for N in Ns:
                r, _ = run(N, eps, kc, K, nmax, Fc)
                print(f"eps={eps:g} kcap={kc:g} N={N} K={K} n={nmax}: lam_min={r['lam']:.6g} "
                      f"| vac: tau={r['vac_tau']:.5g} F(ell)={r['vac_F']:.5g} ell={r['vac_ell']:.4g} "
                      f"F(<ell>)={r['F_of_vac_ell']:.5g} Xi(1)={r['vac_X']:.4g} | gs: tau={r['gs_tau']:.5g} "
                      f"F={r['gs_F']:.5g} | s={np.round(r['s'],5).tolist()}", flush=True)
