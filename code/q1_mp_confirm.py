"""Q1 first test: higher-precision (mpmath, 40 digits) confirmation of the M=2 coherent-state violation.
Case: eps=0.1, kcap=0.01, coherent state a=(-1,0) of the Williamson mode, Fock cutoff K=34.
Everything is rebuilt in mpmath (operators, ell_hat = sum |Delta_k| by eigh, tau_hat by Gauss-Legendre nq nodes with
mp exp(Xi), F_eps via the exact closed-form dmax root, not the spline) and compared with the float64 pipeline."""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ[v] = "2"
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, mpmath as mp
import shear_clock_quantum as sc
mp.mp.dps = 40   # the paper-9 module needs 40 for its closed form
eps, kc, K, NQ = 0.1, 0.01, 34, int(sys.argv[1]) if len(sys.argv) > 1 else 6
a_disp = (-1 + 0j)/np.sqrt(2)           # alpha for a=(-1,0)
M = 2; mod = sc.model(M, eps, kc*eps**2); h = mod['h']; L = mod['L']
Lm = mp.matrix(L.tolist())
def mat(f, n): return mp.matrix(n, n) if f is None else f
def build(K):
    n = K+2
    a = mp.zeros(n, n)
    for i in range(1, n): a[i-1, i] = mp.sqrt(i)
    ad = a.T
    q = (a+ad)/mp.sqrt(2); p = (a-ad)/(1j*mp.sqrt(2))
    cut = lambda X: X[:K, :K]
    return cut(q), cut(p), cut(q*q), cut(p*p), cut((q*p+p*q)/2)
q, p, qq, pp, qp = build(K)
Dl = [Lm[k, 0]*q + Lm[k, 1]*p for k in range(M)]
D2 = [Lm[k, 0]**2*qq + Lm[k, 1]**2*pp + 2*Lm[k, 0]*Lm[k, 1]*qp for k in range(M)]
def fun(A, f):
    E, U = mp.eigh(A)
    return U*mp.diag([f(e) for e in E])*U.H, E, U
def exp_h(A):
    E, U = mp.eigh(A); return U*mp.diag([mp.exp(e) for e in E])*U.H
# tau
xg, wg = np.polynomial.legendre.leggauss(NQ)
xg = [mp.mpf(float(x)) for x in xg]; wg = [mp.mpf(float(w)) for w in wg]
# refine nodes/weights in mp
xg = [mp.findroot(lambda x: mp.legendre(NQ, x), x) for x in xg]
wg = [2/((1-x**2)*mp.diff(lambda t: mp.legendre(NQ, t), x)**2) for x in xg]
hm = mp.mpf(h); epsm = mp.mpf(eps)
Id = mp.eye(K)
tau = mp.zeros(K, K); Xacc = mp.zeros(K, K)
for k in range(M):
    for x, w in zip(xg, wg):
        s = (x+1)/2*hm
        Xt = Xacc + D2[k]*s/(4*hm*hm)
        tau += epsm*(w*hm/2)*mp.exp(k*hm+s)*exp_h(Xt)
    Xacc = Xacc + D2[k]/(4*hm)
tau = (tau+tau.H)/2
ell = mp.zeros(K, K)
for D in Dl:
    A, _, _ = fun(D, abs); ell += A
ell = (ell+ell.H)/2
lam, V = mp.eigh(ell)
# coherent state: D(alpha)|0> amplitudes c_n = e^{-|a|^2/2} a^n/sqrt(n!)  (exact, truncated)
alpha = mp.mpc(float(a_disp.real), float(a_disp.imag))
c = mp.matrix([mp.exp(-abs(alpha)**2/2)*alpha**n/mp.sqrt(mp.factorial(n)) for n in range(K)])
nrm = mp.sqrt((c.H*c)[0]); c = c/nrm
amp = V.H*c; tc = V.H*(tau*c)
rows = []
for j in range(K):
    pj = abs(amp[j])**2
    if pj >= mp.mpf('1e-3'):
        Fj = sc.dmax  # exact F via root of dmax: solve dmax(tau)=ell for tau
        ell_j = lam[j]
        F_exact = mp.findroot(lambda t: sc.dmax(t, eps) - ell_j, (1-epsm) + mp.mpf('1e-4'), tol=mp.mpf(10)**-25) \
            if False else None
        wv = (tc[j]/amp[j]).real
        rows.append((j, ell_j, pj, wv))
# F by bisection on tau for dmax(tau)=ell (dmax increasing in tau)
def F_exact(l):
    lo, hi = (1-epsm) + mp.mpf(10)**-12, (1-epsm) + mp.mpf(10)**3
    for _ in range(160):
        mid = (lo+hi)/2
        if sc.dmax(mid, eps) < l: lo = mid
        else: hi = mid
    return (lo+hi)/2
Fc = sc.Fclock(eps)
print(f"mp dps=40, eps={eps}, kcap={kc}, K={K}, NQ={NQ}, state coherent a=(-1,0), dim {K}")
mn = None
for j, l, pj, wv in rows:
    Fx = F_exact(l); Fs = float(Fc(np.array([float(l)]))[0])
    m = wv - Fx
    print(f"  level {j}: ell={mp.nstr(l,12)} p={mp.nstr(pj,6)} Re tau_w={mp.nstr(wv,15)} F_exact={mp.nstr(Fx,15)} F_spline={Fs:.15f} margin={mp.nstr(m,8)}")
    if mn is None or m < mn[0]: mn = (m, j, l, pj)
print(f"MIN margin (mp, p>=1e-3) = {mp.nstr(mn[0],10)} at level {mn[1]}, ell={mp.nstr(mn[2],10)}, p={mp.nstr(mn[3],6)}")
# float64 reference, same K, same state
Fcf = sc.Fclock(eps)
tauf, ellf, Flf, Xf = sc.clock_ops(mod, K, Fcf, nq=NQ)
lf, Vf = np.linalg.eigh(ellf)
n = np.arange(K); import math
cf = np.array([np.exp(-abs(a_disp)**2/2)*a_disp**k/math.sqrt(math.factorial(k)) for k in n]); cf /= np.linalg.norm(cf)
af = Vf.conj().T@cf; tf = Vf.conj().T@(tauf@cf)
mf = [((tf[j]/af[j]).real - float(Fcf(np.array([lf[j]]))[0]), lf[j], abs(af[j])**2) for j in range(K) if abs(af[j])**2 >= 1e-3]
print("float64 same K: min margin", min(mf)[0], "at ell", min(mf)[1], "p", min(mf)[2])
