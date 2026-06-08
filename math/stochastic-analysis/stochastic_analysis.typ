// ============================================================
//  Stochastic Analysis Summary — Typst Template
//  Based on Karatzas & Shreve, "Brownian Motion and
//  Stochastic Calculus", 2nd ed.
// ============================================================

// ── Page setup ──────────────────────────────────────────────
#set page(
  paper: "a4",
  margin: (top: 2.5cm, bottom: 2.5cm, left: 3cm, right: 3cm),
  numbering: "1",
  number-align: center,
)

#set text(font: "New Computer Modern", size: 11pt, lang: "en")
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "1.1")

// ── Theorem environments ─────────────────────────────────────
#let thmbox(title, color, body) = block(
  width: 100%,
  inset: (x: 12pt, y: 10pt),
  radius: 3pt,
  stroke: (left: 3pt + color, rest: 0.5pt + color.lighten(50%)),
  fill: color.lighten(92%),
)[
  #text(weight: "bold", fill: color)[#title]
  #h(0.5em)
  #body
]

#let theorem-counter  = counter("theorem")
#let def-counter      = counter("definition")
#let prop-counter     = counter("proposition")
#let lemma-counter    = counter("lemma")
#let remark-counter   = counter("remark")

#let theorem(title: none, body) = {
  theorem-counter.step()
  let num = theorem-counter
  let head = if title != none [Theorem #num (#title)] else [Theorem #num]
  thmbox(head, rgb("#1a4f8a"), body)
}

#let definition(title: none, body) = {
  def-counter.step()
  let num = def-counter
  let head = if title != none [Definition #num (#title)] else [Definition #num]
  thmbox(head, rgb("#276749"), body)
}

#let proposition(title: none, body) = {
  prop-counter.step()
  let num = prop-counter
  let head = if title != none [Proposition #num (#title)] else [Proposition #num]
  thmbox(head, rgb("#6b3fa0"), body)
}

#let lemma(title: none, body) = {
  lemma-counter.step()
  let num = lemma-counter
  let head = if title != none [Lemma #num (#title)] else [Lemma #num]
  thmbox(head, rgb("#b45309"), body)
}

#let remark(body) = {
  remark-counter.step()
  let num = remark-counter
  thmbox([Remark #num], rgb("#555555"), body)
}

#let proof(body) = block(width: 100%, inset: (left: 12pt))[
  #text(style: "italic")[Proof.] #body #h(1fr) $square$
]

// ── Code block ───────────────────────────────────────────────
// Usage: #code("python", `...`)
#let code(lang, body) = block(
  width: 100%,
  inset: 10pt,
  radius: 3pt,
  fill: rgb("#f5f5f5"),
  stroke: 0.5pt + rgb("#cccccc"),
)[
  #text(size: 8.5pt, fill: rgb("#444444"))[_#lang _]
  #v(2pt)
  #raw(body, lang: lang)
]

// ── Title page ───────────────────────────────────────────────
#align(center)[
  #v(3cm)
  #text(size: 22pt, weight: "bold")[
    Brownian Motion and Stochastic Calculus
  ]
  #v(0.4cm)
  #text(size: 14pt, fill: rgb("#555555"))[
    A Summary
  ]
  #v(0.6cm)
  #line(length: 60%, stroke: 0.8pt + rgb("#aaaaaa"))
  #v(0.6cm)
  #text(size: 12pt)[Your Name]
  #v(0.2cm)
  #text(size: 11pt, fill: rgb("#666666"))[
    Based on Karatzas & Shreve, _Brownian Motion and Stochastic Calculus_, 2nd ed.
  ]
  #v(0.3cm)
  #text(size: 11pt, fill: rgb("#666666"))[#datetime.today().display("[month repr:long] [year]")]
  #v(3cm)
]

#pagebreak()

// ── Table of contents ────────────────────────────────────────
#outline(indent: 1.5em, depth: 2)
#pagebreak()

// ============================================================
//  1. PROBABILITY SPACE AND FILTRATIONS
// ============================================================
= Probability Space and Filtrations

#definition(title: "Probability Space")[
  A _probability space_ is a triple $(Omega, cal(F), PP)$ where $Omega$ is a sample
  space, $cal(F)$ is a $sigma$-algebra on $Omega$, and $PP : cal(F) -> [0,1]$ is a
  probability measure with $PP(Omega) = 1$.
]

#definition(title: "Filtration")[
  A _filtration_ on $(Omega, cal(F), PP)$ is a family
  ${cal(F)_t}_{t >= 0}$ of sub-$sigma$-algebras of $cal(F)$ such that
  $cal(F)_s subset.eq cal(F)_t$ for all $s <= t$.
  The filtration satisfies the _usual conditions_ if it is right-continuous
  ($cal(F)_t = cal(F)_{t+}$) and $cal(F)_0$ contains all $PP$-null sets.
]

// ============================================================
//  2. BROWNIAN MOTION
// ============================================================
= Brownian Motion

#definition(title: "Standard Brownian Motion")[
  A stochastic process $W = {W_t}_{t >= 0}$ on $(Omega, cal(F), PP)$ is a
  _standard Brownian motion_ (Wiener process) if:
  + $W_0 = 0$ almost surely;
  + $W$ has independent increments: for $0 <= t_0 < t_1 < dots < t_n$,
    the r.v.'s $W_{t_1} - W_{t_0}, dots, W_{t_n} - W_{t_(n-1)}$ are independent;
  + $W_t - W_s tilde cal(N)(0, t - s)$ for all $0 <= s < t$;
  + $t |-> W_t$ is continuous $PP$-a.s.
]

#theorem(title: "Existence (Wiener, 1923)")[
  There exists a probability space $(Omega, cal(F), PP)$ carrying a standard
  Brownian motion.
]

#proof[
  See Karatzas & Shreve §2.2. The construction proceeds via the Kolmogorov
  extension theorem on the space $C[0,∞)$ equipped with the cylinder
  $sigma$-algebra, followed by Kolmogorov's continuity criterion.
]

#theorem(title: "Non-differentiability")[
  Almost surely, $t |-> W_t$ is nowhere differentiable.
]

#remark[
  Although Brownian paths are continuous, they have infinite variation on
  every interval, which motivates the need for a new integration theory.
]

// Simulation code example
#v(0.5em)
#raw(lang: "python",
```python
import numpy as np
import matplotlib.pyplot as plt

def simulate_bm(T: float = 1.0, n: int = 1000, seed: int = 42) -> tuple:
    """Simulate a standard Brownian motion on [0, T]."""
    rng  = np.random.default_rng(seed)
    dt   = T / n
    dW   = rng.normal(0, np.sqrt(dt), n)
    W    = np.concatenate([[0], np.cumsum(dW)])
    t    = np.linspace(0, T, n + 1)
    return t, W

t, W = simulate_bm(T=1.0, n=1000)
plt.plot(t, W)
plt.xlabel("$t$"); plt.ylabel("$W_t$")
plt.title("Sample path of standard Brownian motion")
plt.tight_layout(); plt.savefig("bm_path.pdf")
```.text
)

// ============================================================
//  3. MARTINGALES
// ============================================================
= Martingales

#definition(title: "Martingale")[
  An $cal(F)_t$-adapted process $M = {M_t}_{t >= 0}$ with $EE[|M_t|] < ∞$
  for all $t$ is a _martingale_ if
  $ EE[M_t | cal(F)_s] = M_s, quad forall 0 <= s <= t. $
  It is a _supermartingale_ (resp. _submartingale_) if "$=$" is replaced by
  "$<=$" (resp. "$>=$").
]

#theorem(title: "Doob's Optional Stopping")[
  Let $M$ be a continuous martingale and $tau$ a bounded stopping time.
  Then $EE[M_tau] = EE[M_0]$.
]

#proposition[
  $W_t$, $W_t^2 - t$, and $exp(sigma W_t - sigma^2 t \/ 2)$ are all
  martingales with respect to the natural filtration of $W$.
]

// ============================================================
//  4. STOCHASTIC INTEGRATION
// ============================================================
= Stochastic Integration

== The Itô Integral

#definition(title: "Itô Integral for Simple Processes")[
  For a simple adapted process
  $H_t = sum_(i=0)^(n-1) H_(t_i) bb(1)_{(t_i, t_(i+1)]}(t)$,
  the _Itô integral_ is defined as
  $ integral_0^T H_t dif W_t := sum_(i=0)^(n-1) H_(t_i)(W_(t_(i+1)) - W_(t_i)). $
]

#theorem(title: "Itô Isometry")[
  For $H in cal(L)^2([0,T] times Omega)$,
  $ EE lr([lr(|integral_0^T H_t dif W_t|)^2]) = EE[integral_0^T H_t^2 dif t]. $
]

#proof[
  For simple processes the identity follows by expanding the square and using
  the independence of non-overlapping increments. The general case follows by
  a density argument; see K&S Proposition 3.2.10.
]

== Itô's Formula

#theorem(title: "Itô's Formula (1-dimensional)")[
  Let $f in C^2(RR)$ and $X_t = X_0 + integral_0^t mu_s dif s + integral_0^t sigma_s dif W_s$.
  Then
  $ f(X_t) = f(X_0)
      + integral_0^t f'(X_s) mu_s dif s
      + integral_0^t f'(X_s) sigma_s dif W_s
      + 1/2 integral_0^t f''(X_s) sigma_s^2 dif s. $
]

#remark[
  The last term $(1\/2) f''(X_s) sigma_s^2 dif s$ is the _Itô correction_,
  absent in ordinary calculus. It arises because $(dif W_t)^2 = dif t$.
]

// ============================================================
//  5. STOCHASTIC DIFFERENTIAL EQUATIONS
// ============================================================
= Stochastic Differential Equations

#definition(title: "SDE")[
  A _stochastic differential equation_ (SDE) is written informally as
  $ dif X_t = mu(t, X_t) dif t + sigma(t, X_t) dif W_t, quad X_0 = x_0, $
  and rigorously means the integral equation
  $ X_t = x_0 + integral_0^t mu(s, X_s) dif s + integral_0^t sigma(s, X_s) dif W_s. $
]

#theorem(title: "Existence and Uniqueness (K&S Theorem 5.2.1)")[
  Suppose $mu, sigma : [0,T] times RR -> RR$ satisfy, for some $K > 0$:
  - *(Lipschitz)* $|mu(t,x) - mu(t,y)| + |sigma(t,x) - sigma(t,y)| <= K|x - y|$;
  - *(Linear growth)* $|mu(t,x)| + |sigma(t,x)| <= K(1 + |x|)$.
  Then there exists a unique (up to indistinguishability) continuous adapted
  process $X$ satisfying the SDE, with $EE[sup_{t in [0,T]} X_t^2] < ∞$.
]

#proof[
  Apply the Picard–Lindelöf iteration:
  set $X^0_t equiv x_0$ and
  $ X^(n+1)_t = x_0 + integral_0^t mu(s, X^n_s) dif s
                      + integral_0^t sigma(s, X^n_s) dif W_s. $
  Using Itô isometry and Gronwall's inequality one shows
  $EE[|X^(n+1)_t - X^n_t|^2] -> 0$,
  so the sequence converges in $L^2$ to a solution. Uniqueness follows
  similarly; see K&S §5.2 for details. $square$
]

== Geometric Brownian Motion

#definition(title: "Geometric Brownian Motion")[
  The SDE
  $ dif S_t = mu S_t dif t + sigma S_t dif W_t, quad S_0 = s_0 > 0 $
  has the explicit solution
  $ S_t = s_0 exp lr(((mu - sigma^2/2) t + sigma W_t)). $
]

#raw(lang: "python",
```python
def simulate_gbm(
    s0: float = 1.0, mu: float = 0.05, sigma: float = 0.2,
    T: float = 1.0, n: int = 252, seed: int = 0
) -> tuple:
    """Euler–Maruyama scheme for Geometric Brownian Motion."""
    rng  = np.random.default_rng(seed)
    dt   = T / n
    dW   = rng.normal(0, np.sqrt(dt), n)
    S    = np.empty(n + 1); S[0] = s0
    for i in range(n):
        S[i+1] = S[i] * (1 + mu * dt + sigma * dW[i])  # Euler–Maruyama
    t = np.linspace(0, T, n + 1)
    return t, S
```.text
)

// ============================================================
//  REFERENCES
// ============================================================
#pagebreak()
= References

#set par(hanging-indent: 1.5em)

[KS91] Karatzas, I. and Shreve, S. E. (1991). _Brownian Motion and Stochastic
Calculus_, 2nd ed. Springer, New York.

[RY99] Revuz, D. and Yor, M. (1999). _Continuous Martingales and Brownian
Motion_, 3rd ed. Springer, Berlin.

[Oks03] Øksendal, B. (2003). _Stochastic Differential Equations: An Introduction
with Applications_, 6th ed. Springer, Berlin.
