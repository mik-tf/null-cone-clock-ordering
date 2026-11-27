---
title: "A one-mode quantization of the clock bound on a null cone, and its sensitivity to operator ordering"
author: "Michaël Perreault"
date: "2026-11-26"
format: "article"
toc: false
header_footer_policy: "all"
footer: "© 2026 Michaël Perreault · CC BY 4.0"
pageof: true
date_footer: "YYYY-MM-DD"
header_title: "A one-mode quantization of the clock bound"
---

*Corresponding author: Michaël Perreault, independent researcher, ORCID 0009-0003-0017-7658.*

**Abstract.** On a vacuum null cone the affine length of a generator is a clock $\tau=\varepsilon\int e^{t}e^{\Xi}dt$, bounded below by a convex function of the conformal length, $\tau\ge F_\varepsilon(\ell)$ [Per26]. We ask whether this survives quantization in Chandrasekaran's null-cone scheme [Cha26], where the clock is the exponential of a Weyl-ordered shear. For a free system with $N$ canonical pairs we prove, for every $N$, an operator family bound $\hat\tau\ge A(c)+\hat\ell^{2}/W(c)$, the expectation bound $\langle\hat\tau\rangle\ge F_\varepsilon(\langle\hat\ell\rangle)$ in every state, and the operator inequality on a spectral window. Ordering matters at order $\hbar$: the expectation bound holds in the Weyl ordering, fails at the vacuum for small $\kappa=\eta/\varepsilon^{2}$ when the zero-point constant is removed from the shear $\Xi$ with $\hat\ell$ unchanged (Proposition 6), and survives, at the vacuum and for small $\kappa$, removal of any fraction $\lambda<1$ of it (Proposition 6a). The full operator inequality is open; we give its obstruction and small-$N$ numerics.

**Keywords.** null cone; quantization of gravity; Raychaudhuri equation; focusing; Weyl ordering; operator inequality; zero-point energy; operator ordering.

**MSC 2020.** 83C45, 81S05, 47A63, 83C30.

# Introduction {#p9.intro}

The classical focusing statement is the Raychaudhuri equation; its quantum versions are of two kinds. One keeps the equation and corrects it: Das replaces geodesics by Bohmian trajectories [Das14], and Bak, Parikh, Sarkar and Setti compute the effect of a quantized metric on a bundle of light rays and find a correction of the same sign as the classical terms for the vacuum gravitational field [BPSS24]. The other replaces area by generalised entropy: the Quantum Focussing Conjecture of Bousso, Fisher, Leichenauer and Wall [BFLW16]. Neither concerns an inequality between the affine length of a null generator and the shape change of its conformal metric. To our knowledge, in the sources we could consult, there is no earlier statement of such an inequality for a quantized clock.

Such an inequality exists classically. In [Per26] we proved that on a vacuum null cone the clock $\tau=\varepsilon\int_0^\Lambda e^{t+X(t)}dt$, with $\Lambda=\ln(1/\varepsilon)$ and $X'=a^{2}/4$ for a speed $a\ge0$, satisfies $\tau\ge F_\varepsilon(\ell)$, $\ell=\int a\,dt$, with $F_\varepsilon$ convex, increasing, explicit and sharp. The quantization of Chandrasekaran [Cha26] turns $\Xi$ into a Weyl-ordered operator (4.22) and the clock into an operator (4.26), and the bracket along a generator is (A.36). The question of this note is the one a reader of both papers would ask: does $\hat\tau\ge F_\varepsilon(\hat\ell)$ survive?

**What is new.** (i) A quantum model of the sector and the expectation bound for all states and all $N$, strengthened to the root-mean-square length and the discrete bound. (ii) An operator inequality for all $N$ with a free function $c$, by a tangent construction, and the operator inequality on a spectral window. (iii) The role of ordering: with the whole vacuum zero point removed the inequality fails already for expectations, while removing any fraction less than one does not, at the vacuum and for small $\kappa$ (Proposition 6a). (iv) The open problem, with its obstruction, and eigenvalue computations at small $N$. **Not new:** the Weyl ordering, the clock operator and the bracket, which are in [Cha26]; the classical bound, which is in [Per26]; the operator Cauchy-Schwarz inequality and Jensen's inequality for the spectral measure, which are standard.

**Prior art.** We searched for the quantum Raychaudhuri equation, quantum focusing, clock or length bounds in null-cone quantization, and operator Jensen inequalities for the exponential (Section 8). To our knowledge, in the sources we could consult, none of Theorem 1, Corollary 3 or Proposition 6 appears earlier. Our search was limited to what we could read freely; the sources we did not read in full text are named in Section 8 and in the references.

# The model {#p9.model}

The definitions are those of [Cha26] specialised to one sector; the reduction is ours and is spelled out in the subsection "Faithfulness to the source" below. Take one null generator, one angular harmonic and one polarisation, $z=z^{1}$, with $z^{2}=0$. On this line the hyperbolic metric is $dz^{2}$, so $\Xi=\tfrac14\int s\,(\partial_s z)^{2}ds$ and $\ell=\int|\partial_t z|\,dt$ are exact, and the sector is a free system. The boundary term of (A.36) is set to zero. This is a choice, which the source does not determine (see below), and the Limitations section lists it.

Put $t=\ln(v/\varepsilon)\in[0,\Lambda]$ and divide $[0,\Lambda]$ into $M=2N$ cells of width $h=\Lambda/M$, with nodes $v_k=\varepsilon e^{kh}$, node operators $z_k=z(v_k)$, increments $\Delta_k=z_k-z_{k-1}$, and $z$ piecewise linear in $t$. The bracket of (A.36) restricted to the cone reads $[z(v),z(w)]=i\eta\,\mathrm{sgn}(w-v)/(vw)$, with $\eta=4\pi G\hbar/(\rho_S^{2}\cdot\text{angular area})$. The matrix $[\Delta_j,\Delta_k]$ is a nondegenerate $c$-number matrix for even $M$, so there are $N$ canonical pairs and a Fock space on the Williamson modes of $\sum_k\Delta_k^{2}$. The dimensionless parameter at the cap is $\kappa=\eta/\varepsilon^{2}$ (the metric fluctuations scale like $1/\rho^{2}$ and are largest at $v=\varepsilon$).

On cell $k$ the Weyl-ordered $\Xi$ is
$$
\hat\Xi(t)=\sum_j\beta_j(t)\,\Delta_j^{2},\qquad \beta_j=\tfrac1{4h}\ (j<k),\quad \beta_k=\tfrac{t-t_{k-1}}{4h^{2}},\quad \beta_j=0\ (j>k),
$$
and, because $\Delta_j^{2}$ is the square of one linear form, its Weyl ordering is the plain square. The clock and length are
$$
\hat\tau=\varepsilon\int_0^\Lambda e^{t}\exp\bigl(\hat\Xi(t)\bigr)\,dt,\qquad \hat\ell=\sum_k|\Delta_k|,
$$
with the exponential and $|\cdot|$ by spectral calculus. This is the ordering $e^{\mathrm{Op}^{W}(\Xi)}$, which is the one [Cha26] prescribes in (C.10); the ordering $\mathrm{Op}^{W}(e^{\Xi})$ is not used (see Limitations). Here $\hat\ell\ge0$, $\hat\Xi(t)\ge0$ and $\hat\tau\ge1-\varepsilon$. Quadratic-form statements are on finite-particle vectors, with values in $[0,+\infty]$ (see "Domain and finiteness" below).

The **discrete classical bound** is $F^{(M)}(x)=\min\{\tau(\Delta):\Delta_k\ge0,\ \sum\Delta_k=x\}$, where $\tau(\Delta)$ is the classical clock of the piecewise-linear profile. Piecewise-constant speeds are admissible in [Per26], so $F^{(M)}\ge F_\varepsilon$. The continuum extremal has $\Xi(1)=\infty$ and is not piecewise constant, so $F^{(M)}(x)>F_\varepsilon(x)$ for $x>0$. *Proof of strictness.* In the integration by parts of [Per26, Theorem 1] it suffices that the profile $G$ be absolutely continuous, which a piecewise-constant speed gives; and $G(0)=\varepsilon e^{-\Xi(1)}>0$ because $\Xi(1)$ is finite for a piecewise-constant speed, whereas the sharp bound requires $G(0)=0$ with $\Xi(1)=\infty$. Computed relative gap $(F^{(M)}-F_\varepsilon)/(F^{(M)}-(1-\varepsilon))$: $0.25$ to $0.71$ for $x\in[0.1,4]$, $M\in\{2,4,8\}$, $\varepsilon=0.1$.

## Domain and finiteness {#p9.dom}

**Domain.** The clock is $\hat\tau=\varepsilon\int_0^\Lambda e^{t}e^{\hat\Xi(t)}dt$ with each $\hat\Xi(t)\ge0$ a positive self-adjoint operator (closure of its restriction to finite-particle vectors, Nelson's analytic vector theorem). We consider $\hat\tau$ only as a *quadratic form* on the dense space $\mathcal D$ of finite-particle vectors of the Williamson Fock space, with values in $[0,+\infty]$:
$$
\langle\psi,\hat\tau\psi\rangle=\varepsilon\int_0^\Lambda e^{t}\,\|e^{\hat\Xi(t)/2}\psi\|^{2}dt\in[0,+\infty].
$$
All inequalities between forms in this note (Lemmas 1 and 2, Theorem 1, Proposition 2, Corollaries 3 to 5) are inequalities of $[0,+\infty]$-valued forms on $\mathcal D$; where the left side is $+\infty$ they hold trivially, and Lemmas 1 and 2 need no finiteness. The form domain of $\hat\tau$ (used in Corollary 3) is $\{\psi:\langle\psi,\hat\tau\psi\rangle<\infty\}$. Whether $\mathcal D$ lies in it depends on $\kappa$ as follows.

**Lemma F (finiteness criterion, vacuum).** *Proved*, using Landau's theorem on Laplace transforms of positive measures [Wid41] for the converse direction (analytic continuation; numerically validated). Write $e^{b\hat Q}\,r\,e^{-b\hat Q}=T(b)\,r$ with $T(b)=\exp(-2ibJS)$ for $\hat Q=r^{T}Sr$, and in the $(a,a^\dagger)$ basis $T=\begin{pmatrix}U&V\\ \bar V&\bar U\end{pmatrix}$. Where $U(b)$ is invertible, $e^{b\hat Q}|0\rangle=\det U(b)^{-1/2}\exp(\tfrac12a^{\dagger}Z(b)a^{\dagger})|0\rangle$ with $Z=-U^{-1}V$ symmetric. Then
$$
\langle0|e^{2\beta\hat Q}|0\rangle=\|e^{\beta\hat Q}|0\rangle\|^{2}<\infty\iff\|Z(b)\|_{\rm op}<1\ \text{ for all }b\in(0,\beta],
$$
and then it equals $|\det U(\beta)|^{-1}\det(1-Z(\beta)\bar Z(\beta))^{-1/2}$. In particular, with $S=\kappa S_1(t)$ and $Z_t(\beta)$ the matrix for $\hat Q(t)$ at flow parameter $\beta$: $\langle e^{\hat\Xi(t)}\rangle_0<\infty$ iff $\|Z_t(\beta)\|<1$ for all $\beta\in(0,\kappa/2]$, and $\langle e^{2\hat\Xi(t)}\rangle_0<\infty$ iff $\|Z_t(\beta)\|<1$ for all $\beta\in(0,\kappa]$. The same criterion governs every finite-particle vector (a polynomial in $a^\dagger$ times the squeezed vector has finite norm iff $\|Z\|<1$; stated, not separately checked numerically).

*Proof.* The Heisenberg flow follows from $[\hat Q,r_i]=-2i(JSr)_i$. The vector $\psi(b)=\det U^{-1/2}e^{\frac12a^\dagger Za^\dagger}|0\rangle$ is annihilated by $Ua+Va^\dagger$, hence is the unique (up to scale, fixed by $\langle0|\psi\rangle=\det U^{-1/2}$) Gaussian solution of $\partial_b\psi=\hat Q\psi$ with $\psi(0)=|0\rangle$; it is a norm-analytic function of $b$ while $\|Z\|<1$. On the other side $b\mapsto e^{b\hat Q}|0\rangle$ is analytic in the strip $\mathrm{Re}\,b<b^{*}$ where $b^{*}=\sup\{b:\int e^{2b\lambda}d\mu_0(\lambda)<\infty\}$ (spectral theorem; $\hat Q\ge0$). Both agree for small $b$ (power series, analytic vectors), hence on the common connected domain. If $\|Z(b_1)\|\to1$ as $b\uparrow b_1$ with $U$ finite, then $\|\psi(b)\|^{2}=|\det U|^{-1}\det(1-Z\bar Z)^{-1/2}\to\infty$, while it is monotone in $b$ by the spectral theorem, so $b^{*}\le b_1$. This proves the direction "finite $\Rightarrow\|Z(b)\|<1$ on $(0,\beta]$". For the converse suppose $\|Z\|<1$ on $(0,\beta]$ but $b^{*}<\beta$. Then $f(b)=\int e^{2b\lambda}d\mu_0(\lambda)$ is the Laplace transform of a positive measure with abscissa of convergence $b^{*}$ (finite), so by Landau's theorem (Pringsheim's theorem for Dirichlet series) $f$ has a singularity at the real point $b^{*}$ [Wid41, Chapter II, Theorem 5b]. But $f(b)=\det U(b)^{-1}\det(1-Z(b)\bar Z(b))^{-1/2}$ on $(0,b^{*})$, and the right side is analytic in a neighbourhood of every real $b$ with $\|Z(b)\|<1$ ($T(b)$ is entire in $b$), so $f$ continues analytically through $b^{*}$, a contradiction. Hence the norm is finite exactly up to the first crossing of $\|Z\|=1$. $\square$ Proposition 6 does not depend on the converse direction: step (iii) of its proof uses the criterion only through the analytic-vector proof for small $\kappa$ and the computed $\kappa_2$.

*Remark (one mode).* For $\hat Q=xq^{2}$, $\|Z(\beta)\|=\beta x/(1-\beta x)$, which is below $1$ iff $2\beta x<1$, matching the direct Gaussian integral $\langle e^{2\beta xq^2}\rangle=(1-2\beta x)^{-1/2}$ (validated in item 1 below).

**Check (computed), `code/finiteness_check.py`, output `code/finiteness_check.out`.**

1. *Validation.* At $N=1$ the criterion was compared with a direct Fock sum $\langle0|e^{bQ_K}|0\rangle$ ($K=30,60,120,240$, scaling-and-squaring $\exp$; an eigendecomposition sum is useless once $e^{bw_{\max}}$ is huge) at $t=h$ (pure $q^2$, exact threshold $bx=1$, found $b^*=18.5804=1/x$) and $t=1.5h$ (squeezed, $b^*=23.1596$), at $b/b^*=0.2,0.5,0.9,0.99,1.05,1.5$: below threshold the Fock sums increase to $\det U^{-1/2}$ (slowly near $0.99$: $6.99,8.65,9.68,9.98$ against exact $10$), above threshold they grow without limit with $K$. An $N=2$ check at small $b$ agrees to $10^{-7}$. **Validation passed.** (At $t=\Lambda$ with $N=1$ all $\beta_j$ are equal and $\hat Q\propto\hat n+\tfrac12$, no squeezing, finite for all $\kappa$; the test therefore uses intermediate $t$.)
2. *Tabulated points.* For every $(N,\varepsilon,\kappa)$ of the Section 5 vacuum table and the Section 6 eigenvalue table, with $t$ on a grid of 5 points per cell, $\max_t\|Z\|<1$ for both moments: $\langle e^{\hat\Xi}\rangle_0$ and $\langle e^{2\hat\Xi}\rangle_0$ are finite, so the Fock-truncated numerics converge (they are sums of finite quantities) and $\hat\tau$ is finite on $\mathcal D$ at all tabulated points. Largest values at $\kappa=1$: $\|Z(1/2)\|,\|Z(1)\|$ = $0.0277,0.0569$ ($N=1,\varepsilon=0.1$); $0.0988,0.2192$ ($N=2$); $0.1960,0.4875$ ($N=3$); $0.3228,0.9532$ ($N=4$, the closest to a threshold); $0.0989,0.2195$ ($N=4$, $\varepsilon=0.01$).
3. *Thresholds* (computed, all $t$ on the grid; $\kappa_1$ for $\langle e^{\hat\Xi}\rangle$, $\kappa_2$ for $\langle e^{2\hat\Xi}\rangle$, the one used in the proof of Proposition 6): $\varepsilon=0.1$: $N=1$: $\kappa_1=18.58$, $\kappa_2=9.29$; $N=2$: $5.562$, $2.781$; $N=3$: $3.051$, $1.526$; $N=4$: $2.049$, $1.025$. $\varepsilon=0.01$: $N=1$: $101.2$, $50.61$; $N=2$: $18.32$, $9.158$; $N=4$: $5.556$, $2.778$. So $\hat\tau$ is finite on the vacuum iff $\kappa<\kappa_1(N,\varepsilon)$ (and, by the criterion, on every finite-particle vector), $\kappa_1$ decreasing with $N$ in these data; its behaviour for large $N$ is *not* established. For $\kappa\ge\kappa_1$ one has $\langle0,\hat\tau0\rangle=+\infty$ and the Fock-truncated numerics do not converge.

**Ordering of the exponential and a threshold.** The Weyl symbol of $e^{\Xi}$ diverges at a threshold of the symplectic eigenvalues, so $\mathrm{Op}^{W}(e^{\Xi})$ does not exist beyond a threshold in $N$ (computed from the model, not proved in general). The ordering used here, $e^{\mathrm{Op}^{W}(\Xi)}$, has a threshold of the same kind: its vacuum expectation is finite iff $\kappa<\kappa_1(N,\varepsilon)$ (Lemma F), and $\kappa_1$ is $2.05$ at $N=4$, $\varepsilon=0.1$; all tabulated points lie below it. Statements of this note are for $\kappa<\kappa_1(N,\varepsilon)$ where a finite value is asserted, and are trivial (form value $+\infty$) beyond.

## Faithfulness to the source {#p9.faith}

The source is [Cha26] ([arXiv:2611.34924](https://arxiv.org/abs/2611.34924)v1, 25 November 2026), read in full text. The five statements the model rests on are quoted verbatim (mathematical glyphs retyped from the PDF; page numbers are those of v1).

**(A.36), p. 58, with the continuation on p. 59 (the bracket).** "The inverse along a null generator involves parallel transport P(v, w) from z(w) to z(v). Integration by parts yields the boundary term -(1/(16 pi G_N)) [ rho^2 g_{alpha beta} delta z^beta delta' z^alpha ]_boundary, (A.35), which enters the boundary conditions together with delta Theta_boundary from the cuts. With these contributions included the kernel is (Delta_C^{-1})^{alpha beta}(v, w) = - (8 pi G_N sgn(v - w) / (2 rho(v) rho(w))) P^alpha_gamma(v, w) g^{gamma beta}(z(w)) + k_boundary^{alpha beta}(v, w). (A.36)" and, on p. 59, "The conditions at the cuts and corners in (A.14) fix k_boundary. On a null cone, rho^2 = rho_S(y)^2 v^2, so the coefficient of parallel transport is - sgn(v - w)/(2 K(y) v w)." Here $K(y)=\rho_S(y)^2/(8\pi G_N)$ (eq. (3.63b)).

**(4.22), p. 30 (the shear).** "Defining z^1 = Q(z^1), z^2 = Q(z^2) we find that Xi_y(v) = (1/4) integral_eps^v [ s (d_s z^1)^dagger d_s z^1 + (d_s z^2)^dagger e^{2 z^1} d_s z^2 ] ds, (4.22) which is nonnegative since e^{2 z^1} is positive."

**(4.26), p. 31 (the clock).** "Using the same angular average as in the classical expression, the clock operator takes the form a = integral_S rho_S(y) eps_0, u*bar = a^{-1/2} [ integral_S eps_0(y) integral_eps^1 dv rho_S(y) e^{Xi_y(v)} rho_S(y) ] a^{-1/2}. (4.26)"

**(C.9) and (C.10), p. 65 (the ordering of the exponential).** "For functions f, g of smeared canonical variables, Weyl ordering expresses the product of their operators as [61], Q(f) Q(g) = Q(f star g), f star g = f exp[ (i hbar / 2) (left-arrow delta_I) J^{IJ} (right-arrow delta_J) ] g. (C.9)" and "To evaluate e^{Xi_y(v)}, let x denote the phase space function obtained by multiplying the factors in (4.22) with the indicated ordering using (C.9), so that Xi_y(v) = Q(x). Order by order, e^{Xi_y(v)} = Q(1 + x + (1/2) x star x + (1/3!) x star x star x + ...). (C.10)"

**What these fix.** The shear is Weyl quantized, $\Xi_y(v)=Q(x)$, and $e^{\Xi_y(v)}$ in (4.26) is the exponential of that operator, expanded in (C.10): this is $e^{\mathrm{Op}^{W}(\Xi)}$, the ordering used here, and not $\mathrm{Op}^{W}(e^{\Xi})$. No normal ordering of $\Xi$ is prescribed anywhere in [Cha26], so the zero-point term is part of the source's operator; the normal-ordered variant of the headline section is our own comparison.

**Reductions made here.** Each is a restriction of the source, not something the source states.

1. *One generator.* One fixed $y$ on the corner; several generators are not studied.
2. *Frozen area.* The area operator $a$ and the density $\rho_S$ of (4.24) and (4.26) are held fixed. With $\rho_S$ constant, $\rho_S e^{\Xi}\rho_S=\rho_S^{2}e^{\Xi}$ and the normalisation $a^{-1/2}(\cdots)a^{-1/2}$ cancels, leaving the plain integral $\tau$. This agrees with the source's clock, which "preserves the corner area on incoming wavefunctions". The cross terms between $z$ and $\delta(\rho^{2})$ are dropped and not studied.
3. *One polarisation, $z^{2}=0$.* The second term of (4.22) and the coupling $e^{2z^{1}}$ drop, so $\Xi=\tfrac14\int s(\partial_sz)^{2}ds$.
4. *The $y$-independent mode.* $z(y,v)=a(v)$ with $a=(1/\mathrm{Area})\int\epsilon_0 z$. Smearing the kernel of (A.36) over the corner then gives $[a(v),a(w)]=\text{kernel}/\mathrm{Area}$; this is the origin of the angular area in $\eta$. A harmonic normalised by $\int\epsilon_0Y^{2}=1$ would give the factor $1$ instead.
5. *Trivial parallel transport.* Along $z^{2}=0$ with $z^{1}$ varying, the curve is a geodesic of $\mathrm{diag}(1,e^{2z^{1}})$ and $\partial_{z^{1}}$ is parallel along it ($\Gamma^{1}_{11}=\Gamma^{2}_{11}=0$), so $P(v,w)$ is the identity.
6. *The boundary term.* We set $k_{\rm boundary}=0$. The source says only that the conditions at the cuts and corners in (A.14) fix it and gives no formula in the text. What it does to the nondegeneracy of the bracket matrix for even $M$, or to the number of canonical pairs, is therefore not determined by the source, and we do not claim it is harmless.

*Bracket magnitude.* $\hbar$ times the coefficient $8\pi G_N/(2\rho_S^{2}vw)$ is $4\pi G\hbar/(\rho_S^{2}vw)$, which is the $\eta/(vw)$ of this note (before the division by the area); the sign, $\mathrm{sgn}(w-v)$ against $-\mathrm{sgn}(v-w)$, is a convention.

# The operator family bound {#p9.thmB}

**Lemma 1 (tangent of the exponential).** *Proved.* For real $c$ and self-adjoint $A$, $e^{A}\ge e^{c}(1+A-c)$. Indeed $e^{x}\ge e^{c}(1+x-c)$ for real $x$ by convexity, applied in the spectral calculus of $A$. The right side is affine in $A$, so it is monotone in $A$, while $e^{A}$ is not.

**Lemma 2 (weighted operator Cauchy-Schwarz).** *Proved.* For self-adjoint $X_1,\dots,X_M$, not necessarily commuting, and $w_k>0$: $(\sum_kX_k)^{2}\le(\sum_k1/w_k)(\sum_kw_kX_k^{2})$. The difference equals $\sum_{j<k}(\sqrt{w_j/w_k}\,X_j-\sqrt{w_k/w_j}\,X_k)^{2}\ge0$: expand $(aX-bY)^{2}=a^{2}X^{2}+b^{2}Y^{2}-ab(XY+YX)$ with $ab=1$ and compare with $(\sum X_k)^{2}=\sum X_k^{2}+\sum_{j<k}(X_jX_k+X_kX_j)$. Only squares of self-adjoint operators occur, never a product inequality, so noncommutativity costs nothing.

**Theorem 1 (operator family bound).** *Proved.* For every measurable $c:[0,\Lambda]\to\mathbb R$ put
$$
A(c)=\varepsilon\!\int e^{t+c(t)}\bigl(1-c(t)\bigr)dt,\quad w_k(c)=\varepsilon\!\int e^{t+c(t)}\beta_k(t)\,dt>0,\quad W(c)=\sum_k\frac1{w_k(c)}.
$$
Then for every $N$, $\varepsilon$ and $\eta$,
$$
\hat\tau\ \ge\ A(c)+\sum_kw_k(c)\,\Delta_k^{2}\ \ge\ A(c)+\frac{\hat\ell^{2}}{W(c)}=:G_c(\hat\ell).
$$

*Proof.* Lemma 1 with $A=\hat\Xi(t)$, $c=c(t)$, multiplied by $\varepsilon e^{t}>0$ and integrated (an integral of operator inequalities against a positive weight; likewise for any positive quadrature) gives the first inequality, since $\hat\Xi(t)$ is linear in the $\Delta_j^{2}$. Lemma 2 with $X_k=|\Delta_k|$, $X_k^{2}=\Delta_k^{2}$, gives the second. $\square$

**Proposition 2 (envelope).** *Proved.* For every $x\ge0$, $\sup_cG_c(x)=F^{(M)}(x)$, attained at $c=X^{*}_x$, the $\Xi$ of a minimiser of $F^{(M)}(x)$.

*Proof.* (i) $G_c(x)\le F^{(M)}(x)$: for classical increments summing to $x$ the same two steps hold in scalar form. (ii) Put $y_k=\Delta_k^{2}$. The classical clock is convex in $y$ (an exponential of a linear function, integrated) and increasing, and the constraint set $\{y\ge0:\sum\sqrt{y_k}\ge x\}$ is convex (a superlevel set of a concave function), so the minimum $F^{(M)}(x)$ is attained at some $y^{*}$ and the linearisation $L$ of the clock at $y^{*}$ satisfies $L\ge\tau(y^{*})$ on the constraint set (first-order optimality). $L$ is exactly $A(c)+\sum w_k(c)y_k$ at $c=X^{*}_x$; minimising it over $\sum\Delta_k=x$ by Cauchy-Schwarz gives $A+x^{2}/W=G_c(x)$. Hence $G_c(x)\ge F^{(M)}(x)$. $\square$

*Computed check:* $|G_x(x)-F^{(M)}(x)|\le2\times10^{-15}$ and $\Delta^{*}_kw_k$ constant in $k$ (the optimality condition) for all tested $(\varepsilon,M,x)$ (`code/shear_clock_quantum_proof.out`, Section 1). Write $G_x=G_{X^{*}_x}$.

*Computed (operator check, $N=1$, $\varepsilon=0.1$, cutoff 26).* $\lambda_{\min}(\hat\tau-G_x(\hat\ell))=3.9\times10^{-4},\,3.7\times10^{-4},\,2.6\times10^{-2}$ for $x=0,0.5,2$ at $\kappa=0.1$, and $7.9\times10^{-3},\,7.4\times10^{-3},\,2.5\times10^{-2}$ at $\kappa=1$. All positive, as Theorem 1 requires.

# Consequences {#p9.cons}

**Corollary 3 (expectation bound, strengthened).** *Proved.* For every $N$ and every normalised $\psi$ in the form domain of $\hat\tau$ (the vectors with $\langle\psi,\hat\tau\psi\rangle<\infty$, Section 2; for form value $+\infty$ the bound is trivial):
$$
\langle\hat\tau\rangle_\psi\ \ge\ F^{(M)}\bigl(\langle\hat\ell^{2}\rangle^{1/2}\bigr)\ \ge\ F^{(M)}\bigl(\langle\hat\ell\rangle\bigr)\ \ge\ F_\varepsilon\bigl(\langle\hat\ell\rangle\bigr).
$$
*Proof.* Take expectations in Theorem 1: $\langle\hat\tau\rangle\ge A+\langle\hat\ell^{2}\rangle/W=G_c(\langle\hat\ell^{2}\rangle^{1/2})$. Take $\sup_c$ and use Proposition 2. Then $\langle\hat\ell^{2}\rangle\ge\langle\hat\ell\rangle^{2}$ and $F^{(M)}$ is increasing, and $F^{(M)}\ge F_\varepsilon$. $\square$

This contains the weaker form, proved earlier by Jensen on the spectral measure of $\hat\Xi(t)$, then $\langle\Delta_k^{2}\rangle\ge\langle|\Delta_k|\rangle^{2}$, then the classical theorem. It is strictly stronger whenever $\hat\ell$ has variance in $\psi$ or $\langle\hat\ell\rangle_\psi>0$ (the latter through $F^{(M)}>F_\varepsilon$).

**Corollary 4 (spectral reduction and window).** *Proved; window values computed.* Let $\mu_\psi$ be the spectral measure of $\hat\ell$ in $\psi$. Then $\langle\hat\tau-F_\varepsilon(\hat\ell)\rangle_\psi\ge\sup_c\int(G_c-F_\varepsilon)\,d\mu_\psi$, so $\hat\tau$ drops out. Let $J_c=\{\ell\ge0:G_c(\ell)\ge F_\varepsilon(\ell)\}$. Then $\langle\hat\tau\rangle\ge\langle F_\varepsilon(\hat\ell)\rangle$ for every $\psi$ whose $\mu_\psi$ is supported in $J_c$; that is, the operator inequality holds in the quadratic-form sense on the spectral subspace $E_{\hat\ell}(J_c)\mathcal H$. Every $x>0$ lies in $J_{X^{*}_x}$ because $F^{(M)}(x)>F_\varepsilon(x)$. The vacuum tangent $c=0$ gives $G_0(\ell)=(1-\varepsilon)+\ell^{2}/W_0$ with $W_0$ finite for each $M$, while $F_\varepsilon(\ell)-(1-\varepsilon)\sim\ell^{2}/(4\ln(1/(\varepsilon\ell^{2})))=o(\ell^{2})$ [Per26, Proposition 2]. Hence $J_0\supset[0,\ell_*]$ for some $\ell_*>0$. *Computed:* $\ell_*=3.29,\,2.55,\,1.81$ for $M=2,8,32$ at $\varepsilon=0.1$, and $4.23,\,3.52,\,2.67$ at $\varepsilon=0.01$. *Proved:* $W_0=4\ln M+O(1)$ and $\ell_*\to0$ as $M\to\infty$. Indeed, with $j=M-k$ and $\varepsilon e^{kh}=e^{-jh}$, $w_k(0)=\varepsilon\int e^{t}\beta_k\,dt=(1-e^{-jh})/(4h)+r_k$, where $r_k$ is the contribution of cell $k$ itself and $0<r_k=O(1)$ uniformly. Since $x/(1-e^{-x})=1+O(x)$, $1/w_k=4/(j+O(1))\,(1+O(jh))$, so $W_0=\sum_j4/j+O(1)=4\ln M+O(1)$. Since $W_0\to\infty$ while $F_\varepsilon(\ell)-(1-\varepsilon)>0$ is fixed for each $\ell>0$, $G_0(\ell)<F_\varepsilon(\ell)$ at any fixed $\ell>0$ once $M$ is large, so $\ell_*\to0$. The computed increments of $W_0$ per doubling of $M$ ($2.31,\,2.54,\,2.65,\,2.71$) approach $4\ln2=2.77$. For small $\ell_*$ the small-$\ell$ law of [Per26] suggests $\ell_*^{2}$ of order $e^{-W_0/4}/\varepsilon$, that is, of order $1/(\varepsilon M)$ (a heuristic, not proved).

In the numerics of Section 6 the low-lying eigenvectors have $\langle\hat\ell\rangle$ of order one (the computed vacuum value is $\langle\hat\ell\rangle_0=1.04$ at $N=2$, $\kappa=1$, $\varepsilon=0.1$), with a spread that leaves part of the spectral weight outside $[0,\ell_*]$, so Corollary 4 alone does not certify those states.

**Corollary 5 (commuting states).** *Proved.* For every density matrix $\rho$ commuting with $\hat\ell$, $\mathrm{Tr}\,\rho\hat\tau\ge\mathrm{Tr}\,\rho F_\varepsilon(\hat\ell)$. *Proof.* No pure point spectrum of $\hat\ell$ is assumed. $\rho$ is trace class, so each eigenspace $V$ of $\rho$ for an eigenvalue $p>0$ is finite dimensional, and $\rho=\sum p\,P_V$. Because $\rho$ commutes with $\hat\ell$, it commutes with every spectral projection $E_{\hat\ell}(B)$, so $V$ is invariant under all of them. A finite-dimensional subspace invariant under the spectral resolution of a self-adjoint operator reduces it, and the restriction of $\hat\ell$ to $V$ is a self-adjoint operator on a finite-dimensional space; it has an orthonormal basis of eigenvectors, which lie in the domain of $\hat\ell$. Thus $\rho$ is a mixture of joint eigenvectors $\phi$, with $\hat\ell\phi=\ell\phi$, and for each Corollary 3 gives $\langle\hat\tau\rangle_\phi\ge F_\varepsilon(\ell)=\langle F_\varepsilon(\hat\ell)\rangle_\phi$. $\square$ This covers eigenstates and Gibbs states of $\hat\ell$.

# Ordering sensitivity at order $\hbar$ {#p9.zp}

**Headline: ordering sensitivity at order $\hbar$.** Under the source's own ordering (Weyl, as in (4.22), with the exponential taken as $e^{\mathrm{Op}^{W}(\Xi)}$ in the clock (4.26)) the expectation bound $\langle\hat\tau\rangle\ge F_\varepsilon(\langle\hat\ell\rangle)$ holds in the proved senses, in every state (Corollary 3). With full normal ordering, that is, with the whole zero-point constant $\langle0|\hat\Xi(t)|0\rangle$ removed from $\hat\Xi$ while $\hat\ell$ keeps its vacuum fluctuation, the bound fails at the vacuum for $\kappa<\kappa_0(N,\varepsilon)$ (Proposition 6, proved). Removing any fraction $\lambda<1$ of the constant is still consistent with the bound at small $\kappa$, with a vacuum margin that grows like $\ln(1/\kappa)$ relative to the length term (Proposition 6a, proved). So the bound is sensitive to the ordering at order $\hbar$: it holds with the Weyl zero-point lift and fails when that lift is removed completely. We read this as a statement about an inconsistent subtraction (the length operator keeps its zero-point fluctuation, the clock loses its lift, and normal ordering also destroys $\hat\Xi\ge0$) and not as a discovery about what carries positivity. The closest analogue in the literature is the zero-point mechanism of [BPSS24], where vacuum metric noise gives a correction to the null Raychaudhuri equation; it is perturbative, cutoff dependent and restricted to the vacuum, and it does not test an ordering (Section 8). Proposition 6 is the geometric analogue of a known fact: renormalized (normal-ordered) energy densities lose pointwise positivity, which the quantum energy inequalities answer with smeared lower bounds [FR97, Few12, FK20, FKK22]. We do not derive Proposition 6 from that literature, and we read none of it in full text (Section 8). What we state is a sharp geometric clock inequality in a quantized-gravity model, with the threshold $\lambda<1$ for partial removal.

Replace $\hat\Xi(t)$ by its normal-ordered form $\hat\Xi_N(t)=\hat\Xi(t)-\langle0|\hat\Xi(t)|0\rangle$ (Williamson vacuum $|0\rangle$), and $\hat\tau_N=\varepsilon\int e^{t}\exp(\hat\Xi_N(t))dt$. Theorem 1 still holds for $\hat\tau_N$ with a shift: $\hat\tau_N\ge G_c(\hat\ell)-\varepsilon\int e^{t+c(t)}\langle0|\hat\Xi(t)|0\rangle dt$, of order $\eta$ (*proved*, same two lemmas). This shift is what the vacuum can exploit.

**Scaling.** The bracket matrix is proportional to $\eta=\kappa\varepsilon^2$, hence $L=\sqrt\kappa\,L_1$ with $L_1=L_1(N,\varepsilon)$ independent of $\kappa$. Write $\hat\Xi(t)=\kappa\,\hat Q(t)$, $\hat Q(t)=r^{T}S_1(t)r$, $S_1(t)=L_1^{T}\mathrm{diag}(\beta(t))L_1\ge0$, with $\|S_1(t)\|\le C_S$ uniformly in $t$. The vacuum $|0\rangle$ is the same for every $\kappa$. Put $m_1(t)=\langle0|\hat Q(t)|0\rangle=\tfrac12\mathrm{tr}\,S_1(t)$, which is $>0$ for $t>0$, and $X_t=\hat\Xi_N(t)=\kappa(\hat Q(t)-m_1(t))$. Since $\hat Q\ge0$ (a sum of squares of self-adjoint operators), $\hat Q(t)$ is a positive self-adjoint operator (closure on finite-particle vectors, Nelson).

**Proposition 6 (normal-ordered vacuum).** *Proved* (the lower asymptotics of $F_\varepsilon$ is [Per26, Proposition 2]). For every $N$ and $\varepsilon\in(0,1)$ there are constants $c_N>0$, $C_1<\infty$ and $\kappa_0(N,\varepsilon)>0$ such that for $0<\kappa<\kappa_0$:

1. $\langle\hat\ell\rangle_0=c_N\sqrt\kappa$ exactly, with $c_N=\sqrt{2/\pi}\sum_k\sqrt{s_k}$, $s_k=\tfrac12|L_{1,k}|^2>0$;
2. $0\le\langle\hat\tau_N\rangle_0-(1-\varepsilon)\le C_1\kappa^{2}$;
3. $\langle\hat\tau_N\rangle_0<F_\varepsilon(\langle\hat\ell\rangle_0)$, so the expectation bound fails in the vacuum, and so does the operator inequality $\hat\tau_N\ge F_\varepsilon(\hat\ell)$ in the quadratic-form sense on finite-particle vectors.

*Proof.* (i) *Length.* For a real linear form $X=u\cdot r$ one has $\langle0|e^{isX}|0\rangle=e^{-s^{2}|u|^{2}/4}$ (Weyl operator in the Gaussian vacuum), so the spectral measure of $X$ in $|0\rangle$ is $N(0,|u|^{2}/2)$ and $\langle|X|\rangle_0=\sqrt{2/\pi}\,|u|/\sqrt2$. Apply this to each $\Delta_k=\sqrt\kappa L_{1,k}\cdot r$: $\langle|\Delta_k|\rangle_0=\sqrt{2/\pi}\sqrt{\kappa s_k}$, and $\langle\hat\ell\rangle_0=\sum_k\langle|\Delta_k|\rangle_0=c_N\sqrt\kappa$. If $s_k=0$ then $\Delta_k=0$ as an operator, contradicting the nondegenerate bracket $[\Delta_j,\Delta_k]=iB_{jk}$ (row $k$ of $B$ would vanish); hence every $s_k>0$ and $c_N>0$. This is exact, not leading order.

(ii) *Moments of $X_t$.* $\langle X_t\rangle_0=0$ by construction. $\hat Q|0\rangle$ has components only in the $0$- and $2$-particle sectors of the Williamson Fock space, so $(\hat Q-m_1)|0\rangle=\tfrac12 a^{\dagger T}Y(t)a^{\dagger}|0\rangle$ with $Y=Y(S_1)$ linear in $S_1$. Hence $\langle X_t^{2}\rangle_0=\kappa^{2}V(t)$, $V(t)=\|(\hat Q-m_1)|0\rangle\|^{2}=2\,\mathrm{tr}(S_1CS_1C^{T})$ with $C=\tfrac12(I+iJ)$ (Wick), finite, continuous in $t$, $V\le C_VC_S^{2}$. Likewise $X_t^{2}|0\rangle$ lies in the sectors $\le4$ with norm$^2$ homogeneous of degree 4 in $S_1$, so $\langle X_t^{4}\rangle_0=\kappa^{4}V_4(t)\le\kappa^{4}C_4C_S^{4}$.

(iii) *Exponential moment.* By the criterion of Lemma F (Section 2, proved), $\langle0|e^{2\kappa\hat Q(t)}|0\rangle<\infty$ iff $\|Z_t(\beta)\|<1$ for all $\beta\in(0,\kappa]$. Because $T(\beta)=I+O(\beta)$ uniformly in $t\in[0,\Lambda]$ (compact, $S_1$ bounded), $\|Z_t(\beta)\|\le C_Z\beta$ uniformly, so there is $\kappa_2(N,\varepsilon)>0$, explicit by the criterion and computed in Section 2, with $\sup_t\|Z_t(\beta)\|<1$ for $\beta<\kappa_2$. For $\kappa\le\kappa_2/2$ the quantity $E(\kappa)=\sup_t\langle e^{2X_t}\rangle_0=\sup_t e^{-2\kappa m_1(t)}\,|\det U_t(\kappa)|^{-1}\det\bigl(1-Z_t(\kappa)\bar Z_t(\kappa)\bigr)^{-1/2}$ (Lemma F, with $U_t,Z_t$ the matrices of $\hat Q(t)$ at flow parameter $\kappa$) is finite (continuous on a compact set) and bounded by a constant $E_*$ independent of $\kappa$. (An independent proof of finiteness for small $\kappa$, with a radius $1/(2C)$: finite-particle vectors are analytic vectors of the quadratic $\hat Q$, since $\hat Q$ raises the particle number by at most 2 and $\|\hat Q\phi\|\le C(n+2)\|\phi\|$ on the sectors $\le n$.)

(iv) *Taylor bound.* For real $x$, $e^{x}-1-x\le\tfrac12x^{2}(1+e^{x})$ (for $x\le0$ it is $\le x^{2}/2$; for $x\ge0$ it is $\le\tfrac12x^{2}e^{x}$). In the spectral calculus of $X_t$ with Cauchy-Schwarz,
$$
\langle e^{X_t}\rangle_0-1\le\tfrac12\langle X_t^{2}\rangle_0+\tfrac12\langle X_t^{4}\rangle_0^{1/2}\langle e^{2X_t}\rangle_0^{1/2}\le C_1\kappa^{2},\qquad C_1=\tfrac12\bigl(C_VC_S^{2}+\sqrt{C_4C_S^{4}E_*}\bigr).
$$
Jensen on the spectral measure gives $\langle e^{X_t}\rangle_0\ge e^{\langle X_t\rangle_0}=1$. Integrate against $\varepsilon e^{t}dt$ (nonnegative integrand, Tonelli; $\varepsilon\int_0^\Lambda e^{t}dt=1-\varepsilon$): $0\le\langle\hat\tau_N\rangle_0-(1-\varepsilon)\le C_1\kappa^{2}$. This is (2).

(v) *Comparison with [Per26, Proposition 2].* The proposition gives $F_\varepsilon(x)-(1-\varepsilon)\sim x^{2}/(4\ln(1/(\varepsilon x^{2})))$ as $x\to0$, a two-sided statement: there are $c,C,x_0>0$ with $c\,x^{2}/\ln(1/x)\le F_\varepsilon(x)-(1-\varepsilon)\le C\,x^{2}/\ln(1/x)$ for $0<x<x_0$. Put $x=\langle\hat\ell\rangle_0=c_N\sqrt\kappa$ and $L=\ln(1/\kappa)$. Then $\ln(1/x)=\tfrac12L-\ln c_N\le L$ as soon as $L\ge2|\ln c_N|$, so
$$
F_\varepsilon(\langle\hat\ell\rangle_0)-(1-\varepsilon)\ \ge\ \frac{c\,c_N^{2}\,\kappa}{\ln(1/\kappa)}\ >\ C_1\kappa^{2}\ \ge\ \langle\hat\tau_N\rangle_0-(1-\varepsilon)
$$
as soon as $\kappa\,\ln(1/\kappa)<c\,c_N^{2}/C_1$, which holds for all small $\kappa$ because $\kappa/\ln(1/\kappa)\gg\kappa^{2}$. Define $\kappa_0$ as the minimum of: (a) $\kappa_2/2$ (step (iii), so that $C_1$ is finite); (b) $x_0^{2}/c_N^{2}$ (so that $x<x_0$); (c) $\exp(-2|\ln c_N|)\le1$ (so that $0<\ln(1/x)\le\ln(1/\kappa)$); (d) $\kappa_{\log}$, the largest $\kappa'\le1/e$ with $\kappa'\ln(1/\kappa')\le c\,c_N^{2}/C_1$ (the function $\kappa\ln(1/\kappa)$ increases on $(0,1/e]$, so the inequality holds for all $\kappa<\kappa_{\log}$). For $0<\kappa<\kappa_0$ all four conditions hold. This is (3), first part.

(vi) *Operator version.* Suppose $\hat\tau_N\ge F_\varepsilon(\hat\ell)$ as forms on finite-particle vectors. Evaluate in $|0\rangle$: $\langle\hat\tau_N\rangle_0<\infty$ by (iv), so $\langle F_\varepsilon(\hat\ell)\rangle_0\le\langle\hat\tau_N\rangle_0<\infty$, and Jensen for the convex $F_\varepsilon$ on the spectral measure of $\hat\ell$ in $|0\rangle$ gives $\langle F_\varepsilon(\hat\ell)\rangle_0\ge F_\varepsilon(\langle\hat\ell\rangle_0)$. Hence $\langle\hat\tau_N\rangle_0\ge F_\varepsilon(\langle\hat\ell\rangle_0)$, contradicting (v). $\square$

**Remarks (asymmetry, positivity).** (1) The subtraction is asymmetric: $\hat\ell$ keeps its vacuum fluctuation ($\langle\hat\ell\rangle_0\propto\sqrt\kappa$) while the clock loses its vacuum lift, so the failure is the expected consequence of an inconsistent subtraction; see Proposition 6a for what is and is not shown. (2) $\hat\Xi_N(t)\ge0$ is lost: wherever the vacuum variance $V(t)>0$ the operator has mean $0$ and $X_t|0\rangle\ne0$, so its spectrum contains negative values; the classical $\Xi$ is never negative. (At $t=\Lambda$, $N=1$, $\hat\Xi$ is a multiple of $\hat n+\tfrac12$ and $V=0$, so the loss occurs at intermediate $t$.) (3) The fourth column of the table below is $\langle\hat\tau_N\rangle-(1-\varepsilon)$; the note after the table treats the one row of the output file that was contaminated by quadrature error.

*Computed* (vacuum; cutoff 40 for $N=1$, 14 for $N=2$):

| $\varepsilon$ | $N$ | $\kappa$ | $\langle\hat\tau_N\rangle-(1-\varepsilon)$ | $F_\varepsilon(\langle\hat\ell\rangle)-(1-\varepsilon)$ | $\lambda_{\min}(\hat\tau_N-F_\varepsilon(\hat\ell))$ |
|---|---|---|---|---|---|
| 0.1 | 1 | 1 | $1.96\times10^{-4}$ | $7.88\times10^{-3}$ | $-1.17\times10^{-2}$ |
| 0.1 | 1 | 0.01 | $1.88\times10^{-8}$ | $5.11\times10^{-5}$ | $-7.79\times10^{-5}$ |
| 0.1 | 2 | 1 | $7.02\times10^{-4}$ | $3.23\times10^{-2}$ | $-4.51\times10^{-2}$ |
| 0.1 | 2 | 0.01 | $6.14\times10^{-8}$ | $1.92\times10^{-4}$ | $-2.27\times10^{-4}$ |
| 0.01 | 1 | 0.1 | $4.88\times10^{-8}$ | $1.74\times10^{-4}$ | $-2.17\times10^{-4}$ |
| 0.01 | 2 | 0.001 | $2.01\times10^{-11}$ | $5.61\times10^{-6}$ | $-6.24\times10^{-6}$ |

The quadratic scaling of the fourth column and the linear-over-logarithm scaling of the fifth are visible (the fourth after removing the quadrature bias described next); all 16 rows (both $\varepsilon$, $N=1,2$, $\kappa=1,0.1,0.01,0.001$) are in `code/shear_clock_quantum_proof.out`, Section 5. The reading: in the Weyl ordering the zero-point energy of the Williamson modes lifts $\langle\hat\Xi\rangle$ by an amount linear in $\eta$, which is what covers the linear-in-$\kappa$ cost of a nonzero $\langle\hat\ell\rangle$. Remove it and the covering is gone. The ordering matters because of zero-point energy, not only because of the classical profile.

*Quadrature bias (one row of the output file).* The row $\varepsilon=0.01$, $N=1$, $\kappa=0.001$ of `code/shear_clock_quantum_proof.out` lists $1.449\times10^{-12}$ for the fourth column. This value is a 6-point Gauss-Legendre bias of $-3.43\times10^{-12}$ in $\varepsilon\int e^{t}dt$ (cell width $h=2.3$), not a double-precision cancellation: a 40-digit computation and the closed Gaussian formula reproduce it. The bias-free value is $4.874\times10^{-12}$, and the bias-free excess at $\kappa=10^{-2},10^{-3},10^{-4}$ is $4.8748\times10^{-10}$, $4.8744\times10^{-12}$, $4.874\times10^{-14}$, ratios $100.01$ and $100.01$ to the digits shown, so the $\kappa^{2}$ law holds. The other entries carry a quadrature bias below $4\times10^{-12}$ (relative $7\times10^{-3}$ at the $\varepsilon=0.01$, $N=1$, $\kappa=0.01$ entry, negligible at $N=2$ and at $\varepsilon=0.1$).

**Proposition 6a (partial zero-point removal).** *Proved, given [Per26, Proposition 2].* For $0\le\lambda<1$ let $\hat\Xi_\lambda(t)=\hat\Xi(t)-\lambda\langle0|\hat\Xi(t)|0\rangle$ and $\hat\tau_\lambda=\varepsilon\int e^{t}e^{\hat\Xi_\lambda(t)}dt$ ($\lambda=0$ is the Weyl ordering, $\lambda=1$ is $\hat\tau_N$). Then for every $N$, $\varepsilon$ and $\lambda<1$ there is $\kappa_\lambda(N,\varepsilon)>0$ such that for $0<\kappa<\kappa_\lambda$
$$
\langle\hat\tau_\lambda\rangle_0-(1-\varepsilon)\ \ge\ (1-\lambda)\,a\,\kappa,\qquad a=\varepsilon\!\int_0^\Lambda e^{t}m_1(t)\,dt>0,
$$
$$
F_\varepsilon(\langle\hat\ell\rangle_0)-(1-\varepsilon)\ \le\ \frac{3C\,c_N^{2}\,\kappa}{\ln(1/\kappa)},
$$
hence $\langle\hat\tau_\lambda\rangle_0\ge F_\varepsilon(\langle\hat\ell\rangle_0)$ with margin at least $\kappa\,[(1-\lambda)a-3Cc_N^{2}/\ln(1/\kappa)]$, and the ratio of the vacuum excess of the clock to that of $F_\varepsilon(\langle\hat\ell\rangle)$ grows at least like $(1-\lambda)\,a\,\ln(1/\kappa)/(3Cc_N^{2})$. Here $c,C,x_0$ are the constants of the two-sided law of [Per26, Proposition 2] (step (v) of the proof of Proposition 6) and $c_N$ is as in Proposition 6. One may take $\kappa_\lambda=\min\{\exp(-6Cc_N^{2}/((1-\lambda)a)),\ \exp(-6|\ln c_N|),\ x_0^{2}/c_N^{2},\ 1\}$.

*Proof.* In the vacuum spectral measure of $\hat\Xi(t)$, Jensen gives $\langle e^{\hat\Xi}\rangle_0\ge e^{m(t)}$, $m=\kappa m_1$; this holds in $[0,\infty]$ with no finiteness assumption. Then $\langle e^{\hat\Xi_\lambda}\rangle_0=e^{-\lambda m}\langle e^{\hat\Xi}\rangle_0\ge e^{(1-\lambda)m}\ge1+(1-\lambda)m$. Integrate against $\varepsilon e^{t}dt$. The upper bound is the upper side of the two-sided law of [Per26, Proposition 2], with $x=c_N\sqrt\kappa<x_0$. Put $L=\ln(1/\kappa)$; then $\ln(1/x)=\tfrac12L-\ln c_N$, and for $\kappa<\kappa_\lambda$ we have $L>6|\ln c_N|$, so $\ln(1/x)\ge\tfrac12L-\tfrac16L=\tfrac13L$ and $x^{2}/\ln(1/x)\le3c_N^{2}\kappa/L$ (the factor $\tfrac12$ in $\ln(1/x)$ is the reason for $3C$ and not $C$). Also for $\kappa<\kappa_\lambda$, $L>6Cc_N^{2}/((1-\lambda)a)$, so $3Cc_N^{2}/L<\tfrac12(1-\lambda)a$ and the margin is at least $\tfrac12(1-\lambda)a\kappa>0$. $\square$

So every fraction $\lambda<1$ of the Weyl zero-point constant may be removed and the expectation bound still holds in the vacuum for $\kappa$ small enough (the threshold $\kappa_\lambda$ shrinks like $\exp(-\mathrm{const}/(1-\lambda))$); only the complete removal $\lambda=1$ fails (Proposition 6). The Weyl ordering ($\lambda=0$) has a large margin, linear in $\kappa$, against a competitor of order $\kappa/\ln(1/\kappa)$. *What is not shown:* this concerns the vacuum expectation only; the operator inequality for $\lambda<1$ is not claimed (it remains Conjecture 2 for $\lambda=0$). Both claims are for small $\kappa$ at fixed $N$ and $\varepsilon$; we do not claim uniformity in $N$.

# Numerics for the full operator inequality {#p9.num}

Fock basis of the Williamson modes with per-mode cutoff $K$; operators built with correctly truncated quadratics; $\hat\tau$ by 6-point Gauss-Legendre per cell; exponentials, $|\cdot|$ and $F_\varepsilon$ by spectral calculus in the large space; then Rayleigh-Ritz on states of total occupation $\le n$. $F_\varepsilon$ is tabulated from its closed form [Per26] (relative error of $\tau-(1-\varepsilon)$ below $3\times10^{-6}$ against the closed form). Rayleigh-Ritz gives upper bounds on the lowest eigenvalue of the truncated problem, so a positive Ritz value does not by itself certify positivity of $\hat\tau-F_\varepsilon(\hat\ell)$; the numerics are evidence, not proof. Lowest Ritz value $\lambda_{\min}$ of $\hat\tau-F_\varepsilon(\hat\ell)$, $\varepsilon=0.1$ unless stated:

| $\varepsilon$ | $N$ | $K$ | $n$ | $\kappa=0.01$ | $\kappa=0.1$ | $\kappa=1$ |
|---|---|---|---|---|---|---|
| 0.1 | 1 | 20 | 10 | 0.000249305 | 0.00236075 | 0.0217855 |
| 0.1 | 1 | 22 | 10 | 0.000249442 | 0.00236226 | 0.0217997 |
| 0.1 | 1 | 26 | 10 | 0.00024965 | 0.00236456 | 0.0218211 |
| 0.1 | 2 | 20 | 10 | 0.00141947 | 0.0137592 | 0.141278 |
| 0.1 | 2 | 22 | 10 | 0.0014198 | 0.0137632 | 0.141297 |
| 0.1 | 2 | 26 | 10 | 0.0014203 | 0.0137692 | 0.141316 |
| 0.1 | 3 | 8 | 10 | 0.00332367 | 0.0320765 | 0.375172 |
| 0.1 | 3 | 10 | 10 | 0.00333055 | 0.0321351 | 0.375391 |
| 0.1 | 3 | 12 | 10 | 0.00333443 | 0.0322412 | 0.377224 |
| 0.1 | 4 (not converged) | 5 | 8 | 0.00617947 | 0.0616709 | 0.907962 |
| 0.1 | 4 (not converged) | 6 | 8 | 0.00602785 | 0.0595153 | 0.871667 |
| 0.01 | 1 | 20 | 10 | 5.85133e-05 | 0.00055678 | 0.00517063 |
| 0.01 | 1 | 22 | 10 | 5.85501e-05 | 0.000557177 | 0.00517471 |
| 0.01 | 1 | 26 | 10 | 5.86061e-05 | 0.000557781 | 0.00518094 |
| 0.01 | 2 | 20 | 10 | 0.000459228 | 0.00444177 | 0.0431199 |
| 0.01 | 2 | 22 | 10 | 0.000459413 | 0.00444402 | 0.0431476 |
| 0.01 | 2 | 26 | 10 | 0.000459695 | 0.00444746 | 0.04319 |
| 0.01 | 3 | 8 | 10 | 0.00106184 | 0.0103738 | 0.105916 |
| 0.01 | 3 | 10 | 10 | 0.00106192 | 0.0103692 | 0.105776 |
| 0.01 | 3 | 12 | 10 | 0.00106288 | 0.0103783 | 0.105858 |
| 0.01 | 4 (not converged) | 5 | 8 | 0.00195313 | 0.019351 | 0.209955 |
| 0.01 | 4 (not converged) | 6 | 8 | 0.00189848 | 0.0186815 | 0.201091 |

*Table 1. Lowest Ritz value $\lambda_{\min}(\hat\tau-F_\varepsilon(\hat\ell))$ for $N=1,2,3,4$, per-mode Fock cutoff $K$, restricted to total occupation $\le n$, values as stored. A positive Ritz value is evidence, not a certificate. The output file of each row is in `code/out/` and the command that produced it is in `code/out/README.md` (run from `code/` with `run_table.sh`). Rows $N=1$ and $N=2$ come from one run per $(K,\varepsilon)$.*

*Notes on convergence.* For $N=1,2$ the values grow with $K$ and move by at most 0.2 percent between $K=20$ and $K=26$ (largest case, $\varepsilon=0.01$, $N=1$, $\kappa=1$: $0.00517063$ to $0.00518094$; at $\varepsilon=0.1$ the largest move is 0.16 percent, $2.1786\times10^{-2}$ to $2.1821\times10^{-2}$ at $\kappa=1$). For $N=3$, $K=8,10,12$ agree within 0.33 percent at $\kappa=0.01$ and 0.55 percent at $\kappa=1$ ($\varepsilon=0.1$), and within 0.13 percent at $\varepsilon=0.01$; at $\varepsilon=0.01$, $\kappa=1$ the sequence $0.105916,\ 0.105776,\ 0.105858$ is not monotone. **$N=4$ is not converged.** At $\varepsilon=0.1$, $\kappa=0.01$ the values for $K=5,6,7$ are $6.17947\times10^{-3}$, $6.02785\times10^{-3}$, $6.12948\times10^{-3}$ (a 2.5 percent drop, then a 1.7 percent rise), and at $\kappa=1$ the step $K=5\to6$ drops 4 percent ($0.907962$ to $0.871667$). The $K=7$ run reached $\kappa=0.01$ and $0.1$ only ($\kappa=0.1$: $0.0608549$); it hit its time cap, so it is not in the table. The $N=4$ entries are accurate to a few percent at best, and no limit can be read off. Larger $K$ ($N=3$, $K=14$; $N=4$, $K=8$) was not reached. All stored values are positive. The growth with $N$ at $\varepsilon=0.1$, $\kappa=0.01$ ($2.5\times10^{-4}$, $1.4\times10^{-3}$, $3.3\times10^{-3}$, $6.2\times10^{-3}$) is a factor 25 from $N=1$ to $N=4$, between $N^{2}$ and $N^{3}$.

*Two nearby numbers, labelled.* At $N=1$, $\varepsilon=0.1$, $\kappa=1$, $K=26$ the table gives $\lambda_{\min}=2.182\times10^{-2}$ (more digits: $2.18211\times10^{-2}$), a Ritz value restricted to total occupation $n=10$. The full-space check in the Conjecture 1 paragraph below (`code/shear_clock_quantum_proof.out`, Section 4) gives $2.181\times10^{-2}$ at the same $K$, with no occupation cutoff. These are two different computations, not a rounding difference; the restricted value is the larger, as a Ritz value on a subspace must be.


The vacuum is not the minimiser: at $N=1$, $\varepsilon=0.1$, $\kappa=0.1$ the vacuum margin is $\langle\hat\tau\rangle_0-\langle F_\varepsilon(\hat\ell)\rangle_0=0.90329-0.90067=2.6\times10^{-3}$, and the ground Ritz vector gives slightly less.

*Classical limit, computed ($N=1$, $\varepsilon=0.1$).* Coherent states with mean increments $(0.3,0.5)$ have classical $\tau_{cl}=0.938951$, $F_\varepsilon(\ell_{cl})=0.917181$, gap $2.177\times10^{-2}$. The quantum gap $\langle\hat\tau\rangle-\langle F_\varepsilon(\hat\ell)\rangle$ tends to the classical gap from above as $\kappa\to0$, with errors linear in $\kappa$ ($0.036\kappa$ for $\hat\tau$, $0.0104\kappa$ for $F_\varepsilon(\hat\ell)$).

*Sharper conjecture, tested at $N=1$ only.* **Conjecture 1.** $\hat\tau\ge F^{(M)}(\hat\ell)$ for all $N$ in the Weyl ordering. At $N=1$, $\varepsilon=0.1$: $\lambda_{\min}(\hat\tau-F^{(M)}(\hat\ell))=1.70\times10^{-4}$ ($\kappa=0.1$) and $1.84\times10^{-3}$ ($\kappa=1$), versus $2.36\times10^{-3}$ and $2.18\times10^{-2}$ for $F_\varepsilon$ (full space, no occupation cutoff, $K=26$: $2.364\times10^{-3}$ and $2.181\times10^{-2}$): positive, about ten times smaller. Conjecture 1 implies the $F_\varepsilon$ statement and isolates the zero-point effect from the discretisation gap.

# The open problem and its obstruction {#p9.open}

**Conjecture 2 (the main open problem).** $\hat\tau\ge F_\varepsilon(\hat\ell)$ as an operator inequality for all $N$, $\varepsilon$, $\eta$, in the Weyl ordering $e^{\mathrm{Op}^{W}(\Xi)}$. *Status:* conjectured; supported by Corollary 3, Corollary 4, Corollary 5, and the numerics of Section 6; no counterexample found.

**Remark (pointwise weak-value form; why off-diagonal control is needed).** *Hypothesis:* $\hat\ell$ has an orthonormal eigenbasis $\{|\ell\rangle\}$ in the space considered (Corollary 5 deliberately assumes no pure point spectrum). This holds at $N=1$, where $\hat\ell$ is a multiple of $|q|+|p|$ in one canonical pair and has discrete spectrum; for $N\ge2$ it is not shown for the operator on the full Fock space. In every truncation used below $\hat\ell$ is a finite Hermitian matrix, hence diagonalisable, and all levels quoted are levels of the truncated matrix; they move with the cutoff $K$. For a state $\psi$ and such an eigenvector with $\langle\ell|\psi\rangle\ne0$, the weak value of $\hat\tau$ is $\tau_w(\ell)=\langle\ell|\hat\tau|\psi\rangle/\langle\ell|\psi\rangle$ [AAV88]. Since $\langle\hat\tau\rangle_\psi=\sum_\ell p(\ell)\,\mathrm{Re}\,\tau_w(\ell)$ with $p(\ell)=|\langle\ell|\psi\rangle|^{2}$, the averaged inequality $\langle\hat\tau\rangle_\psi\ge\langle F_\varepsilon(\hat\ell)\rangle_\psi$ would follow from the pointwise bound $\mathrm{Re}\,\tau_w(\ell)\ge F_\varepsilon(\ell)$ for every $\ell$. The pointwise bound is false in every truncation tested (computed). *Computed* ($M=2$, i.e. $N=1$, $\varepsilon=0.1$, $\kappa=0.01$, coherent state with Williamson displacement $a=(-1,0)$, levels with $p\ge10^{-3}$, `code/q1_k_scan.out`): the margin $\mathrm{Re}\,\tau_w(\ell)-F_\varepsilon(\ell)$ is negative at every cutoff tested, $K=24,30,34,40,48,52,60,70$. At the level $\ell\approx0.207$ ($p\approx2.6\times10^{-3}$) it is $-1.455\times10^{-3}$ ($K=24$), $-1.472\times10^{-3}$ ($K=30$), $-1.4825\times10^{-3}$ ($K=34$, with $\ell=0.2071$, $p=2.55\times10^{-3}$, $\mathrm{Re}\,\tau_w=0.899360<F_\varepsilon(\ell)=0.900842$) and about $-1.50\times10^{-3}$ for $K\ge48$. From $K=40$ on a second level at $\ell\approx0.187$ has weight just above the cut ($p=1.00\times10^{-3}$ at $K=40$, $1.09\times10^{-3}$ at $K=70$), with margin $-6.72\times10^{-3}$ ($K=40$), $-6.61\times10^{-3}$ ($K=48$), $-6.56\times10^{-3}$ ($K=52$), $-6.50\times10^{-3}$ ($K=60$), $-6.44\times10^{-3}$ ($K=70$); the minimum over levels with $p\ge10^{-3}$ therefore jumps by a factor of about 4.5 between $K=34$ and $K=40$ only because this level crosses the weight cut, and its magnitude depends on $K$. A 40-digit recomputation at $K=34$ (operators, spectrum of $\hat\ell$ and $F_\varepsilon$ from its closed form rebuilt in mpmath; `code/q1_mp_confirm_nq6.out`) gives $-1.482541001\times10^{-3}$ at $\ell=0.2071$, equal to double precision to nine digits. For $M=4$ ($N=2$, $K=24$, `code/q1_weak_value_clock_M4.out`) the minimum coherent-state margin at $p\ge10^{-3}$ is negative at all six pairs $(\varepsilon,\kappa)$, from $-6.7\times10^{-4}$ at $(0.1,0.01)$ to $-1.96$ at $(0.1,1)$; for $M=2$ ($K=52$, `code/q1_weak_value_clock_M2.out`) it runs from $-1.34\times10^{-3}$ at $(0.01,0.01)$ to $-1.44$ at $(0.1,1)$. The averaged form held on every tested state: over the grid ($M=2$ and $M=4$; coherent, squeezed and 2000 and 1500 random Gaussian states) the smallest value of $\langle\hat\tau\rangle-\langle F_\varepsilon(\hat\ell)\rangle$ is $+6.2\times10^{-5}$, at $(\varepsilon,\kappa)=(0.01,0.01)$, $M=2$, coherent state at the origin ($6.1665\times10^{-5}$ at $K=48$ and $6.1679\times10^{-5}$ at $K=52$, `code/q1_weak_value_clock_M2.out`, lines 11 and 12). The grid is small and the states are Gaussian with modest displacement, so this is evidence, not proof. Weak values of $\hat\tau$ at single levels are anomalous: levels with positive and negative margin alternate, so the average is positive only through cancellation between levels, which is controlled by the off-diagonal elements $\langle\ell|\hat\tau|\ell'\rangle$. This motivates why Conjecture 2 needs control of the off-diagonal elements in the eigenbasis of $\hat\ell$, and why a level-by-level argument cannot work. Scripts and outputs: `code/q1_weak_value_clock.py`, `code/q1_mp_confirm.py`, `code/q1_k_scan.py` and their `.out` files.

**Why the standard tools fail** (statements about the tools; the growth facts are quoted from the operator-monotone and operator-convex literature, [Bha97] Chapter V, whose theorem numbers we have not checked: full text not consulted).

1. *Operator Jensen* (Hansen-Pedersen [FK20] C. J. Fewster and E.-A. Kontou, *A new derivation of singularity theorems with weakened energy hypotheses*, Class. Quantum Grav. 37 (2020) 065010, [doi:10.1088/1361-6382/ab685b](https://doi.org/10.1088/1361-6382/ab685b), [arXiv:1907.13604](https://arxiv.org/abs/1907.13604) (arXiv record checked 2026-11-26; full text not consulted).

[Few12] C. J. Fewster, *Lectures on quantum energy inequalities*, [arXiv:1208.5399](https://arxiv.org/abs/1208.5399) (arXiv record checked 2026-11-26; full text not consulted).

[FKK22] B. Freivogel, E.-A. Kontou and D. Krommydas, *The return of the singularities: applications of the smeared null energy condition*, SciPost Phys. 13 (2022) 001, [doi:10.21468/SciPostPhys.13.1.001](https://doi.org/10.21468/SciPostPhys.13.1.001), [arXiv:2012.11569](https://arxiv.org/abs/2012.11569) (arXiv record checked 2026-11-26; full text not consulted).

[FR97] L. H. Ford and T. A. Roman, *Restrictions on negative energy density in flat spacetime*, Phys. Rev. D 55 (1997) 2082, [doi:10.1103/PhysRevD.55.2082](https://doi.org/10.1103/PhysRevD.55.2082), [arXiv:gr-qc/9607003](https://arxiv.org/abs/gr-qc/9607003) (Crossref and arXiv records checked 2026-11-26; full text not consulted).

[HP03], Davis [Dav57], Choi [Cho74]): $f(\Phi(A))\le\Phi(f(A))$ for all unital positive $\Phi$ holds if and only if $f$ is operator convex. An operator convex function on $[0,\infty)$ grows at most like $x^{2}$. $F_\varepsilon$ grows like $\exp(\ell^{2}/(4\Lambda))$ up to lower order (computed large-clock law of [Per26]), so it is not operator convex. The Berezin anti-Wick map is such a $\Phi$, so a coherent-state sandwich gives only trace inequalities (Berezin-Lieb [Lie73]), not quadratic-form ones.
2. *Löwner.* Comparing $\hat\Xi$ with a function of $\hat\ell$ needs $\exp$ or $F_\varepsilon$ operator monotone; operator monotone functions on $[0,\infty)$ grow at most linearly. Lemma 1 is the strongest monotone replacement of $\exp$, and Theorem 1 is what it yields.
3. *Golden-Thompson* [Gol65, Tho65] bounds $\mathrm{Tr}\,e^{A+B}$ from above, the wrong direction. Peierls-Bogoliubov is Jensen for Gibbs states and reproduces Corollary 3 for those states only.

**The precise obstruction (computed).** By Corollary 4 the tangent method certifies a state if and only if some $c$ gives $\int(G_c-F_\varepsilon)d\mu_\psi\ge0$. For a superposition of two $\hat\ell$ eigenvectors with eigenvalues $s<L$ and weights $p,1-p$, the best margin over 121 tangents ($M=2$, $\varepsilon=0.1$) is, with negative entries bold:

| $s$ | $L$ | $p=0.01$ | $0.1$ | $0.5$ | $0.9$ | $0.99$ |
|---|---|---|---|---|---|---|
| 0.2 | 3 | 0.318 | 0.262 | 0.086 | 0.0089 | 0.0021 |
| 0.2 | 6 | 6.32 | 2.92 | **-1.99** | **-0.66** | **-0.068** |
| 1.0 | 6 | 6.34 | 3.01 | **-1.85** | **-0.62** | **-0.036** |
| 0.2 | 10 | 3668 | **-449** | **-1928** | **-397** | **-39.7** |

($F_\varepsilon(6)=10.0$, $F_\varepsilon(10)=3979$.) Here $s$ and $L$ are prescribed values, not computed eigenvalues of $\hat\ell$ at a given $\kappa$, so the table tests the tangent method on two-point measures, not the operator. The tangent at $L$ is a parabola with a very negative intercept $A$ because $X^{*}_L$ is large, and the tangent at $s$ lies far below $F_\varepsilon(L)$. So the tangent family, the strongest bound obtained by linearising the exponential, cannot certify such states. On such a span the inequality is equivalent to
$$
|\langle\phi_s,\hat\tau\,\phi_L\rangle|^{2}\le\bigl(\langle\hat\tau\rangle_s-F_\varepsilon(s)\bigr)\bigl(\langle\hat\tau\rangle_L-F_\varepsilon(L)\bigr)
$$
(*proved* as a restatement: positivity of $\hat\tau-F_\varepsilon(\hat\ell)$ on $\mathrm{span}\{\phi_s,\phi_L\}$, that is for all weights and phases of the superposition, is equivalent to this $2\times2$ determinant condition, since $F_\varepsilon(\hat\ell)$ is diagonal on eigenvectors of distinct eigenvalues), and both diagonal factors are positive by Corollary 3. **The obstruction is therefore a bound on the off-diagonal matrix elements of $\hat\tau$ in the eigenbasis of $\hat\ell$, that is, on the Fock structure itself.** A proof would plausibly use the Gaussian (symplectic) covariance of the model; we have no such proof.

*Large $N$ (conjectured scaling).* The vacuum has $\langle\hat\ell\rangle$ of order $M\kappa^{1/2}$ and a spread of order $(M\kappa)^{1/2}$, while the window of any single tangent shrinks as $M$ grows. So Theorem 1 certifies no fixed state family uniformly as $N\to\infty$. The ground energy of $\hat\tau$ grows like $\exp(c\,\kappa M^{2})$ and $F_\varepsilon(\langle\hat\ell\rangle)$ like $\exp(M^{2}\eta/(4\Lambda\varepsilon^{2}))$: the exponents have the same form, so the sign at large $N$ is a genuine competition. In the computed data at small $N$ the margin grows with $N$ (an observation from computed data, see Section 6, not a trend we claim beyond them). This is a scaling argument, not a proof.

# Related work and prior art {#p9.prior}

Quantum corrections to the Raychaudhuri equation are in [Das14] (Bohmian trajectories on a fixed background; no null cone, no clock, no ordering question) and [BPSS24] (a quantized metric on a pencil of null geodesics in linearised gravity); the entropic version is [BFLW16]. We read [Das14] and [BPSS24] in full text. [BPSS24] is the closest analogue of the zero-point mechanism of the headline section: in the Poincaré-invariant vacuum the metric noise gives a correction to the null Raychaudhuri equation of the same sign as the classical terms (their eq. (3.10)), but perturbatively, dependent on a cutoff and restricted to the vacuum. They adopt Weyl ordering, by an anticommutator in the equations of motion (their Section 2.3, p. 8), but for the metric fluctuation operators and not for a clock, and they do not test the ordering. None of these papers states an inequality between affine length and conformal excursion, and none bounds a clock operator. Quantizations of null data are in [Wie24, Wie25], [CFL24] and, in July 2026, [CFL24b] L. Ciambelli, L. Freidel and R. G. Leigh, *Null Raychaudhuri: canonical structure and the dressing time*, JHEP 01 (2024) 166, [doi:10.1007/JHEP01(2024)166](https://doi.org/10.1007/JHEP01(2024)166), [arXiv:2309.03932](https://arxiv.org/abs/2309.03932) (Crossref and arXiv records checked 2026-11-26; full text not consulted).

[CK26], read by abstract only. The canonical structure of the null Raychaudhuri equation, with a dressing time, is in [CFL24b] (arXiv record checked, full text not consulted). The quantum energy inequalities, lower bounds on smeared normal-ordered energy densities, are in [FR97] (Ford and Roman), the lectures [Few12], and their use in singularity theorems and focusing in [FK20] and [FKK22] (arXiv records checked, full text not consulted). Proposition 6 is the geometric analogue, for a clock, of the known loss of pointwise positivity under normal ordering that these inequalities address; we claim the sharp clock inequality and its threshold in the one-mode model, not the loss of positivity itself. The scheme used here is [Cha26], read in full text. Operator Jensen inequalities are in [HP03] (abstract and statement of the characterisation read; the claim that $\exp$ is not operator convex rests on [Bha97], full text not consulted). To our knowledge, in the sources we could consult, there is no operator or expectation-value inequality bounding a quantized null-cone clock below by a function of the length operator, and no study of how its ordering matters. The sources we did not read in full text are [BFLW16], [CFL24], [CFL24b], [CK26], [Wie24], [Wie25], [FR97], [Few12], [FK20], [FKK22], [HP03] (beyond its characterisation statement) and [Bha97]; our novelty statements rest on titles and abstracts for these. The search covered arXiv and the open web through 2026-11-26; absence from our search is not proof of absence.

# Limitations {#p9.lim}

One generator, one harmonic and one polarisation (several generators add commuting copies only if cross terms vanish, untested); area and density operators frozen; boundary term $k_{\rm boundary}$ of (A.36) set to zero, and the source text does not determine its effect on the nondegeneracy of the bracket matrix or on the number of canonical pairs; the discretisation $M=2N$ is coarse, so convergence as $M\to\infty$ at fixed $\eta$ is untested; the numerics are Rayleigh-Ritz with truncation, hence upper bounds on the lowest eigenvalue of the truncated problem, and not a certificate of positivity; the ordering $\mathrm{Op}^{W}(e^{\Xi})$ is not used, and for the ordering used the finiteness of $\hat\tau$ needs $\kappa<\kappa_1(N,\varepsilon)$ (Section 2; the large-$N$ behaviour of $\kappa_1$ is observed for $N\le4$ only, and the extension of Lemma F from the vacuum to every finite-particle vector is argued, not separately checked numerically); constants $c,C,x_0$ of the small-$x$ law of [Per26, Proposition 2], used in Propositions 6 and 6a, are not made explicit, and Propositions 6 and 6a are for small $\kappa$ at fixed $N$ and $\varepsilon$, without uniformity in $N$. Cha26 is an unrefereed v1 and equation numbers are those of v1. No claim here is about the full quantum gravity theory.

# Code and data {#p9.code .unnumbered}

`code/shear_clock_quantum.py` (model, Fock space, Rayleigh-Ritz, coherent states; output `code/shear_clock_quantum.out` for $N=1,2$) `code/finiteness_check.py` (finiteness criterion of Lemma F, thresholds, validation against a direct Fock sum; output `code/finiteness_check.out`), and `code/shear_clock_quantum_proof.py` (Section 1 envelope check, windows, two-point obstruction, operator checks at $N=1$, normal-ordered vacuum; output `code/shear_clock_quantum_proof.out`). Python 3 with numpy, scipy and mpmath, one core. Run `python3 shear_clock_quantum_proof.py` from `code/`. The eigenvalue table of Section 6 is reproduced by `code/run_table.sh` (outputs in `code/out/`, commands in `code/out/README.md`). The weak-value test of the remark in Section 7 is `code/q1_weak_value_clock.py` (outputs `code/q1_weak_value_clock_M2.out`, `code/q1_weak_value_clock_M4.out`) with the 40-digit confirmation `code/q1_mp_confirm.py` (outputs `code/q1_mp_confirm_nq6.out`, `code/q1_mp_confirm_nq12.out`); run each with `python3 -I`, the first with the argument 1 or 2 for $N$. The cutoff scan of the remark in Section 7 is `code/q1_k_scan.py` (output `code/q1_k_scan.out`, run with `nice -n 19 python3 -I q1_k_scan.py` from `code/`, 3 seconds).

# Acknowledgements {#p9.ack .unnumbered}

We are grateful to Jason Glowney, whose geometric treatment of light, the light cone and hyperbolic curvature in *Arithmoí*, Book 1 [Glo23], led us to the hyperbolic-plane view of null data on which the classical paper [Per26], and hence this note, is built.

# Use of AI {#p9.ai .unnumbered}

Claude (Anthropic) was used in November 2026, under the author's direction, for exploratory computation, derivations, code, checks, literature search and drafting.

# References {#p9.refs .unnumbered}

Each entry ends with its read status: "read in full text", or "full text not consulted" when we did not read the text (then the entry says what was read, if anything). Records were checked against Crossref and the arXiv listing on 2026-11-26.

[AAV88] Y. Aharonov, D. Z. Albert and L. Vaidman, *How the result of a measurement of a component of the spin of a spin-1/2 particle can turn out to be 100*, Phys. Rev. Lett. 60 (1988) 1351, [doi:10.1103/PhysRevLett.60.1351](https://doi.org/10.1103/PhysRevLett.60.1351) (full text not consulted).

[Bha97] R. Bhatia, *Matrix Analysis*, Graduate Texts in Mathematics 169, Springer, New York, 1997, [doi:10.1007/978-1-4612-0653-8](https://doi.org/10.1007/978-1-4612-0653-8) (full text not consulted).

[BFLW16] R. Bousso, Z. Fisher, S. Leichenauer and A. C. Wall, *Quantum focusing conjecture* (arXiv title: *A quantum focussing conjecture*), Phys. Rev. D 93 (2016) 064044, [doi:10.1103/PhysRevD.93.064044](https://doi.org/10.1103/PhysRevD.93.064044), [arXiv:1506.02669](https://arxiv.org/abs/1506.02669) (full text not consulted).

[BPSS24] S.-E. Bak, M. Parikh, S. Sarkar and F. Setti, *Quantum-gravitational null Raychaudhuri equation*, JHEP 07 (2024) 214, [doi:10.1007/JHEP07(2024)214](https://doi.org/10.1007/JHEP07(2024)214), [arXiv:2312.17214](https://arxiv.org/abs/2312.17214) (v2, 25 Jul 2024; read in full text).

[CFL24] L. Ciambelli, L. Freidel and R. G. Leigh, *Quantum null geometry and gravity*, JHEP 12 (2024) 028, [doi:10.1007/JHEP12(2024)028](https://doi.org/10.1007/JHEP12(2024)028), [arXiv:2407.11132](https://arxiv.org/abs/2407.11132) (full text not consulted).

[Cha26] V. Chandrasekaran, *An algebraic quantization of gravity from the null cone*, [arXiv:2611.34924](https://arxiv.org/abs/2611.34924) (v1, 25 Nov 2026, hep-th; unrefereed, no journal reference; read in full text on 2026-11-26, equations (4.20), (4.22), (4.26), (A.36), (C.9), (C.10); equation numbers are those of v1).

[Cho74] M.-D. Choi, *A Schwarz inequality for positive linear maps on C\*-algebras*, Illinois J. Math. 18 (1974) 565-574, [doi:10.1215/ijm/1256051007](https://doi.org/10.1215/ijm/1256051007) (full text not consulted).

[CK26] L. Ciambelli and M. S. Klinger, *Quantization of gravity on null hypersurfaces*, [arXiv:2607.07785](https://arxiv.org/abs/2607.07785) (v1, 8 Jul 2026; no journal record found; abstract only; full text not consulted).

[Das14] S. Das, *Quantum Raychaudhuri equation*, Phys. Rev. D 89 (2014) 084068, [doi:10.1103/PhysRevD.89.084068](https://doi.org/10.1103/PhysRevD.89.084068), [arXiv:1311.6539](https://arxiv.org/abs/1311.6539) (read in full text).

[Dav57] C. Davis, *A Schwarz inequality for convex operator functions*, Proc. Amer. Math. Soc. 8 (1957) 42-44, [doi:10.1090/S0002-9939-1957-0084120-4](https://doi.org/10.1090/S0002-9939-1957-0084120-4) (full text not consulted).

[Glo23] J. Glowney, *Arithmoí: A Study of Mathematical Cohesiveness in Physics and Metaphysics*, Book 1 (Volume I, Part 1), First Edition, Athene Noctua, 2023; Sections 1.3.1 (Curvature and Light), 1.3.6 (Hyperbolic Curvature), 1.4.2 and 1.4.4 (read in full text).

[Gol65] S. Golden, *Lower bounds for the Helmholtz function*, Phys. Rev. 137 (1965) B1127-B1128, [doi:10.1103/PhysRev.137.B1127](https://doi.org/10.1103/PhysRev.137.B1127) (full text not consulted).

[HP03] F. Hansen and G. K. Pedersen, *Jensen's operator inequality*, Bull. London Math. Soc. 35 (2003) 553-564, [doi:10.1112/S0024609303002200](https://doi.org/10.1112/S0024609303002200), [arXiv:math/0204049](https://arxiv.org/abs/math/0204049) (abstract and characterisation statement read; full text not consulted).

[Lie73] E. H. Lieb, *The classical limit of quantum spin systems*, Commun. Math. Phys. 31 (1973) 327-340, [doi:10.1007/BF01646493](https://doi.org/10.1007/BF01646493) (full text not consulted).

[Per26] M. Perreault, *Affine length bounds the conformal excursion on a null cone: a sharp focusing inequality*, 2026 (preprint), [doi:10.5281/zenodo.24572220](https://doi.org/10.5281/zenodo.24572220) (Zenodo record; read in full text, our own paper).

[Tho65] C. J. Thompson, *Inequality with applications in statistical mechanics*, J. Math. Phys. 6 (1965) 1812-1813, [doi:10.1063/1.1704727](https://doi.org/10.1063/1.1704727) (full text not consulted).

[Wie24] W. Wieland, *Quantum geometry of the null cone*, Phys. Rev. D 112 (2025) 044042, [doi:10.1103/k5vj-61my](https://doi.org/10.1103/k5vj-61my), [arXiv:2401.17491](https://arxiv.org/abs/2401.17491) (abstract only; full text not consulted).

[Wie25] W. Wieland, *Quantum geometry of the light cone: Fock representation and spectrum of radiated power*, Class. Quantum Grav. 42 (2025) 195006, [doi:10.1088/1361-6382/ae0235](https://doi.org/10.1088/1361-6382/ae0235), [arXiv:2504.10802](https://arxiv.org/abs/2504.10802) (abstract only; full text not consulted).

[Wid41] D. V. Widder, *The Laplace Transform*, Princeton Mathematical Series 6, Princeton University Press, Princeton, 1941, Chapter II, Theorem 5b (full text not consulted; the statement used, that the abscissa of convergence of the Laplace-Stieltjes transform of a nondecreasing function is a singular point, is cited from memory of the standard form of Landau's theorem).
