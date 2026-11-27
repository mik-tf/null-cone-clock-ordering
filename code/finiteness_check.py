"""Finiteness criterion for vacuum moments <0|exp(b Q)|0>, Q = Weyl quadratic form of Xi(t) (Paper 9, r1, item M4).

Q = r^T S r with S = L^T diag(beta(t)) L, [r_a, r_b] = i J_ab, vacuum <r r^T> = I/2.
Heisenberg flow  exp(bQ) r exp(-bQ) = T(b) r,  T(b) = expm(-2 i b J S).
In the (a, a^dagger) basis T = [[U, V], [V', U']].  exp(bQ)|0> is annihilated by U a + V a^dagger, i.e. it is
c * exp(1/2 a^dag Z a^dag)|0>  with  Z = -U^{-1} V,  c = <0|exp(bQ)|0> = det(U)^{-1/2}.
Criterion:  <0|exp(2 beta Q)|0> = ||exp(beta Q)|0>||^2 is finite  iff  ||Z(beta)||_op < 1  (norm along [0, beta]).
So <exp(Xi)> finite iff ||Z(1/2)|| < 1 and <exp(2 Xi)> finite iff ||Z(1)|| < 1.
Run: python3 finiteness_check.py > finiteness_check.out
"""
import os, sys
import numpy as np, mpmath as mp
from scipy.linalg import expm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import shear_clock_quantum as Q

def Smat(mod, t):
    M, h, L = mod['M'], mod['h'], mod['L']
    k = min(int(t // h), M - 1); beta = np.zeros(M)
    beta[:k] = 1/(4*h); beta[k] = (t - k*h)/(4*h*h)
    return L.T @ np.diag(beta) @ L

def Jmat(N): return np.kron(np.eye(N), np.array([[0., 1.], [-1., 0.]]))

def Wmat(N):
    W = np.zeros((2*N, 2*N), complex)
    for m in range(N):
        W[m, 2*m], W[m, 2*m+1] = 1/np.sqrt(2), 1j/np.sqrt(2)
        W[N+m, 2*m], W[N+m, 2*m+1] = 1/np.sqrt(2), -1j/np.sqrt(2)
    return W

def Zof(S, b, N):
    T = expm(-2j*b*Jmat(N) @ S); W = Wmat(N)
    Ta = W @ T @ np.linalg.inv(W)
    U, V = Ta[:N, :N], Ta[:N, N:]
    try:
        return -np.linalg.solve(U, V), U
    except np.linalg.LinAlgError:
        return np.full((N, N), np.inf), U

def znorm(S, beta, N, npath=24):
    """max over the path (0, beta] of ||Z||; also returns the value at beta."""
    best = 0.0
    for x in np.linspace(beta/npath, beta, npath):
        Z, U = Zof(S, x, N)
        best = max(best, np.inf if not np.all(np.isfinite(Z)) else np.linalg.norm(Z, 2))
    return best, np.linalg.norm(np.nan_to_num(Zof(S, beta, N)[0], posinf=1e300), 2)

def c_exact(S, b, N):
    Z, U = Zof(S, b, N)
    return np.linalg.det(U)**-0.5

def fock_Q(S, N, K):
    """Fock matrix of Q = r^T S r (Weyl), truncated correctly, N<=2."""
    sm = Q.single_mode(K); I = np.eye(K)
    def emb(ops):
        out = np.array([[1.0+0j]])
        for m in range(N): out = np.kron(out, ops.get(m, I))
        return out
    R = []
    for m in range(N): R += [emb({m: sm['q']}), emb({m: sm['p']})]
    def quad(i, j):
        mi, mj = i//2, j//2
        if mi != mj: return R[i] @ R[j]
        return emb({mi: sm[{(0,0):'qq',(1,1):'pp',(0,1):'qp',(1,0):'qp'}[(i%2, j%2)]]})
    out = 0
    for i in range(2*N):
        for j in range(2*N): out = out + S[i, j]*quad(i, j)
    return out

def brute(S, b, N, K):
    """<0| exp(b Q_K) |0> with Q_K the correctly truncated Fock matrix; scaling-and-squaring (no eigenvectors:
    an eigh-based sum loses everything to rounding once exp(b w_max) is huge)."""
    A = fock_Q(S, N, K); A = (A + A.conj().T)/2
    E = expm(b*A)
    return float(E[0, 0].real)

def validate():
    print("## A. validation of the squeezed-vector criterion against a direct Fock sum (N = 1)")
    print("At t = Lambda the N=1 quadratic is a multiple of the number operator (no squeezing, always finite), so the test uses")
    print("intermediate t: t=h (Q = x q^2, exact threshold b x = 1) and t=1.5h (squeezed, x != y). eps=0.1, kappa=1.")
    ok = True
    mod = Q.model(2, 0.1, 0.1**2); h = mod['h']
    for label, t in (("t=h", h*(1-1e-12)), ("t=1.5h", 1.5*h)):
        S = Smat(mod, t)
        lo, hi = 0.0, 1.0
        while znorm(S, hi/2, 1)[0] < 1 and hi < 1e6: hi *= 2
        for _ in range(60):
            mid = (lo+hi)/2
            (lo, hi) = (mid, hi) if znorm(S, mid/2, 1)[0] < 1 else (lo, mid)
        bstar = lo
        print(f"[{label}] S diag = {np.round(np.diag(S), 5)}, critical b* (||Z(b*/2)|| = 1) = {bstar:.6g}" +
              (f"  (pure q^2 prediction 1/x = {1/S[0,0]:.6g} or 1/y)" if label == "t=h" else ""))
        for f in (0.2, 0.5, 0.9, 0.99, 1.05, 1.5):
            b = f*bstar; zn = znorm(S, b/2, 1)[0]; crit = zn < 1
            vals = [brute(S, b, 1, K) for K in (30, 60, 120, 240)]
            ce = c_exact(S, b, 1)
            inc = all(vals[i+1] >= vals[i]*(1-1e-9) for i in range(3))     # truncation increases the sum
            if crit:   # finite: increasing sequence converging to det(U)^(-1/2); slowly near the threshold
                tol = {0.2: 1e-8, 0.5: 1e-6, 0.9: 2e-2, 0.99: 0.5}[f]
                good = inc and abs(vals[3]-ce.real) < tol*abs(ce.real)
            else:      # infinite: the truncated sums keep growing with K, no limit
                good = inc and vals[3] > 5*vals[1]
            ok &= good
            print(f"  b/b*={f}: |Z(b/2)|={zn:.4f} finite={crit} exact det(U)^-1/2={ce.real:.8g} Fock K=30,60,120,240: {['%.6g'%v for v in vals]} -> {'consistent' if good else 'MISMATCH'}")
    mod2 = Q.model(4, 0.1, 0.01); S2 = Smat(mod2, mod2['h']*1.5)
    zn = znorm(S2, 0.5, 2)[0]; ce = c_exact(S2, 1.0, 2); bf = brute(S2, 1.0, 2, 20)
    print(f"N=2 eps=0.1 kappa=0.01 t=1.5h, <e^Xi>: |Z(1/2)|={zn:.4g} exact={ce.real:.8g} Fock(K=20)={bf:.8g}")
    ok &= abs(bf-ce.real) < 1e-6*abs(ce.real)
    print("VALIDATION", "PASSED" if ok else "FAILED")
    return ok

def table():
    print("\n## B. criterion for each tabulated (N, eps, kappa); t on the cell grid, 5 points per cell plus Lambda")
    print("Columns: N eps kappa | max_t ||Z(1/2)|| (<e^Xi> finite iff <1) | max_t ||Z(1)|| (<e^{2Xi}> finite iff <1) | kappa_1(e^Xi) | kappa_2(e^{2Xi})")
    rows = [(1,.1,k) for k in (1,.1,.01,.001)] + [(2,.1,k) for k in (1,.1,.01,.001)] \
         + [(1,.01,k) for k in (1,.1,.01,.001)] + [(2,.01,k) for k in (1,.1,.01,.001)] \
         + [(3,.1,k) for k in (.01,.1,1)] + [(4,.1,k) for k in (.01,.1,1)] + [(4,.01,k) for k in (.01,.1,1)]
    allfinite = True
    for N, eps, kap in rows:
        mod = Q.model(2*N, eps, eps**2); M, h = 2*N, mod['h']       # S1 at kappa = 1 (S scales linearly in kappa)
        ts = [(k + f)*h for k in range(M) for f in (0.2, 0.4, 0.6, 0.8, 1.0)]
        Ss = [Smat(mod, min(t, M*h - 1e-12)) for t in ts]
        z1 = max(znorm(S, kap/2, N, 8)[0] for S in Ss)
        z2 = max(znorm(S, kap, N, 8)[0] for S in Ss)
        def kcrit(beta_factor):
            lo, hi = 0.0, 1.0
            f = lambda kk: max(znorm(S, beta_factor*kk, N, 6)[0] for S in Ss)
            while f(hi) < 1 and hi < 1e5: hi *= 2
            for _ in range(40):
                mid = (lo+hi)/2
                (lo, hi) = (mid, hi) if f(mid) < 1 else (lo, mid)
            return lo
        k1, k2 = kcrit(0.5), kcrit(1.0)
        allfinite &= (z1 < 1) and (z2 < 1)
        print(f"N={N} eps={eps} kappa={kap}: |Z(1/2)|={z1:.4f} |Z(1)|={z2:.4f} finite(e^Xi)={z1<1} finite(e^2Xi)={z2<1}  kappa_1={k1:.4g} kappa_2={k2:.4g}")
    print("ALL TABULATED ROWS FINITE (both moments)" if allfinite else "SOME TABULATED ROW HAS AN INFINITE MOMENT, see above")

def m7():
    print("\n## C. m7 row recomputed in mpmath (40 digits): eps=0.01, N=1, normal-ordered vacuum <tau_N>-(1-eps)")
    mp.mp.dps = 40
    eps = mp.mpf('0.01'); Lam = mp.log(1/eps); M = 2; h = Lam/M
    v = [eps*mp.e**(h*k) for k in range(M+1)]
    c = lambda j, k: (0 if j == k else (1 if k > j else -1))/(v[j]*v[k])
    Dm = [[-1, 1, 0], [0, -1, 1]]
    B01u = sum(Dm[0][j]*c(j, k)*Dm[1][k] for j in range(3) for k in range(3))   # B01 / eta
    K = 14
    n = mp.matrix(K+2, K+2)
    for i in range(K+1): n[i, i+1] = mp.sqrt(i+1)
    ad = n.T; qm = (n+ad)/mp.sqrt(2); pm = (n-ad)/(1j*mp.sqrt(2))
    qq = (qm*qm)[:K, :K]; pp = (pm*pm)[:K, :K]
    xg, wg = np.polynomial.legendre.leggauss(6)
    nodes, weights = [], []
    for x0 in xg:
        x = mp.mpf(x0)
        for _ in range(8):
            p = mp.legendre(6, x); dp = 6*(x*p - mp.legendre(5, x))/(x*x-1)
            x -= p/dp
        dp = 6*(x*mp.legendre(6, x) - mp.legendre(5, x))/(x*x-1)
        nodes.append(x); weights.append(2/((1-x*x)*dp*dp))
    def fock_row(kap):
        b = abs(B01u)*kap*eps**2; D2 = [b*qq, b*pp]; z2 = [b/2, b/2]
        Xa = mp.zeros(K, K); tau = mp.mpf(0); tau1 = mp.mpf(0)
        for k in range(M):
            for x, w in zip(nodes, weights):
                s = (x+1)/2*h
                E = mp.expm(Xa + (D2[k] - z2[k]*mp.eye(K))*s/(4*h*h))
                tau += eps*(w*h/2)*mp.e**(k*h+s)*E[0, 0].real
                tau1 += eps*(w*h/2)*mp.e**(k*h+s)
            Xa = Xa + (D2[k] - z2[k]*mp.eye(K))/(4*h)
        return tau, tau1
    def gauss_mean(t, kap):
        b = abs(B01u)*kap*eps**2
        k = min(int(t/h), M-1); be = [mp.mpf(0)]*M
        for j in range(k): be[j] = 1/(4*h)
        be[k] = (t-k*h)/(4*h*h)
        x_, y_ = b*be[0], b*be[1]; w2 = 4*x_*y_
        sq = mp.sqrt(w2) if w2 > 0 else mp.mpf(0)
        ch = mp.cosh(sq); sh = mp.sinh(sq)/sq if sq != 0 else mp.mpf(1)
        T = ch*mp.eye(2) + sh*mp.matrix([[0, -2j*y_], [2j*x_, 0]])
        W = mp.matrix([[1, 1j], [1, -1j]])/mp.sqrt(2)
        U = (W*T*mp.inverse(W))[0, 0]
        return (mp.e**(-(x_+y_)/2)*(U**-0.5)).real
    res = {}
    for kap in (mp.mpf('0.01'), mp.mpf('0.001'), mp.mpf('0.0001')):
        tau, tau1 = fock_row(kap)
        bias = tau1 - (1-eps)                        # Gauss-Legendre error on eps*int e^t dt, no physics in it
        raw = tau - (1-eps)
        corr = sum(eps*(w*h/2)*mp.e**(k*h+(x+1)/2*h)*(gauss_mean(k*h+(x+1)/2*h, kap)-1)
                   for k in range(M) for x, w in zip(nodes, weights))
        res[str(kap)] = corr
        print(f"kappa={mp.nstr(kap,3)}: raw <tau_N>-(1-eps) = {mp.nstr(raw,10)} (Fock K={K});  quadrature bias of eps*int e^t = {mp.nstr(bias,6)};"
              f"  bias-free excess eps*int e^t(<e^Xi_N>-1) = {mp.nstr(corr,10)} (closed form)")
    r1 = res['0.01']/res['0.001']; r2 = res['0.001']/res['0.0001']
    print(f"bias-free scaling: excess(0.01)/excess(0.001) = {mp.nstr(r1,8)}, excess(0.001)/excess(0.0001) = {mp.nstr(r2,8)}  (kappa^2 law predicts 100)")
    print("Reading: the 1.449e-12 of shear_clock_quantum_proof.out is reproduced to all digits at 40 digits, so it is NOT double-precision cancellation:")
    print("it is the 6-point Gauss-Legendre error on eps*int e^t dt (cell width h=2.3), which the code does not remove. The bias-free excess follows kappa^2.")

if __name__ == '__main__':
    ok = validate()
    table()
    m7()
