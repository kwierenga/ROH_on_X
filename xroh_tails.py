"""
xroh_tails.py — exact closed-form X-autozygosity tail probabilities + validation
================================================================================
The whole-chromosome outcomes (F_X = 1, fully autozygous; F_X = 0, fully outbred)
are NOT intrinsically Monte-Carlo quantities — they have closed forms read off the
recombination process. This module computes them exactly and validates the
gene-drop simulation in `xroh_sim.py` against them (the simulation must land within
Monte-Carlo error of the analytic value — the scientific unit test).

Model: Haldane (no interference). Crossovers on one female-X meiosis are a Poisson
process of mean Λ = L/100 Morgans (L ≈ 176.3 cM, Λ ≈ 1.763). A key fact used
throughout: the extreme events are "all-or-nothing along the whole chromosome", so
they depend ONLY on the total genetic length Λ, not on the map's spatial shape
(the map shape matters only for intermediate F_X and for classifier accuracies).

Derivations
-----------
FD — one governing meiosis. O_pat = grandfather X (intact). O_mat = a single
  meiosis of {G, V} (V = great-grandmother-derived). A mosaic with k crossovers
  alternates G/V, so it is "entirely G" (F_X=1) or "entirely V" (F_X=0) iff k=0,
  times the ½ start-strand choice:
      P(F_X=1) = P(F_X=0) = ½ e^{-Λ}.

MS — two independent meioses of the SAME mother on {m1,m2}. The XOR (disagreement)
  process flips at rate 2Λ; F_X=1 needs it 0 everywhere, F_X=0 needs it 1
  everywhere — each requires the XOR to start in that phase (½) and never flip
  (e^{-2Λ}):
      P(F_X=1) = P(F_X=0) = ½ e^{-2Λ}.

SS — three governing meioses: O_pat = meiosis of {gm1,gm2} (brother); O_mat =
  meiosis of {GF, U}, with U itself a meiosis of {gm1,gm2} (sister chain).
  • F_X=1 needs O_mat to be entirely U-side (k_s=0 starting on U: ½e^{-Λ}) AND
    U ≡ O_pat everywhere (both single strands, same one: e^{-Λ}·e^{-Λ}·½):
        P(F_X=1) = ¼ e^{-3Λ}.
  • F_X=0 (never autozygous) is a survival probability of a 3-telegraph CTMC
    (state = (O_pat, U, selector)); the two "autozygous" states (selector=U-side
    AND O_pat==U) are made absorbing and we compute the taboo probability
        P(F_X=0) = π0_allowed · exp(Q_allowed · Λ) · 1
    over the 6 non-autozygous states. Exact (map-independent), ≈0.203.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import expm
from xroh_sim import wilson_ci   # shared Wilson score interval helper


def _ss_fx0_generator():
    """6x6 restricted generator for the SS 'never autozygous' survival problem.

    State = (a, u, s): a=O_pat strand, u=U strand, s=selector (0=GF, 1=U-side),
    each a 2-state telegraph flipping at rate 1 per Morgan. Autozygous (bad) states
    are s=1 and a==u, i.e. (0,0,1) and (1,1,1); they are made absorbing by dropping
    their rows/cols, leaving the 6 allowed states."""
    states = [(a, u, s) for a in (0, 1) for u in (0, 1) for s in (0, 1)]
    bad = {(0, 0, 1), (1, 1, 1)}
    allowed = [st for st in states if st not in bad]
    idx = {st: i for i, st in enumerate(allowed)}
    Q = np.zeros((len(allowed), len(allowed)))
    for st in allowed:
        i = idx[st]
        Q[i, i] = -3.0                      # three components each leave at rate 1
        for c in range(3):                  # single-component flips
            nb = list(st); nb[c] ^= 1; nb = tuple(nb)
            if nb in idx:                   # flips into a bad state are absorbing (lost)
                Q[i, idx[nb]] += 1.0
    pi0 = np.full(len(allowed), 1.0 / 8.0)  # uniform over all 8; mass on allowed only
    return Q, pi0


def analytic_tails(L_cM=176.3):
    """Exact Haldane tail probabilities per union. Returns
    {union: {'fx1':.., 'fx0':..}}. Map-independent (depends only on Λ=L/100)."""
    lam = L_cM / 100.0
    Q, pi0 = _ss_fx0_generator()
    ss_fx0 = float(pi0 @ expm(Q * lam) @ np.ones(len(pi0)))
    return {
        "father-daughter": {"fx1": 0.5 * np.exp(-lam),    "fx0": 0.5 * np.exp(-lam)},
        "mother-son":      {"fx1": 0.5 * np.exp(-2 * lam), "fx0": 0.5 * np.exp(-2 * lam)},
        "brother-sister":  {"fx1": 0.25 * np.exp(-3 * lam), "fx0": ss_fx0},
    }


def validate_against_sim(n=300000, seed=11, tol=1e-6):
    """Compare exact analytic extremes to the gene-drop simulation, with Wilson 95%
    CIs. Measures the EXACT events F_X>=1-tol and F_X<=tol (not F_X>0.99, which also
    catches near-telomere slivers — a threshold artifact, see paper §4.2.1)."""
    import xroh_sim as X
    if X.XMAP is None:
        X.init()
    L = X.XMAP[4]
    ana = analytic_tails(L)
    X.reset(seed)
    peds = [("father-daughter", X.FD), ("mother-son", X.MS), ("brother-sister", X.SS)]
    rows = []
    sims = {nm: np.array([X.Fx(*fn()) for _ in range(n)]) for nm, fn in peds}
    for nm, _ in peds:
        fx = sims[nm]
        for kind in ("fx1", "fx0"):
            k = int(np.sum(fx >= 1 - tol) if kind == "fx1" else np.sum(fx <= tol))
            lo, hi = wilson_ci(k, n)
            a = ana[nm][kind]
            rows.append(dict(union=nm, event=("F_X=1" if kind == "fx1" else "F_X=0"),
                             analytic=a, sim=k / n, ci_lo=lo, ci_hi=hi,
                             inside=bool(lo <= a <= hi)))
    return dict(L_cM=float(L), n=n, rows=rows)


if __name__ == "__main__":
    import xroh_sim as X
    X.init()
    L = X.XMAP[4]
    print(f"Female X = {L:.1f} cM (Lambda = {L/100:.3f} Morgans).  Haldane, exact closed forms.\n")
    ana = analytic_tails(L)
    print("Exact tail probabilities (map-independent):")
    print(f"{'union':16s}{'P(F_X=1)':>12}{'P(F_X=0)':>12}")
    for nm, d in ana.items():
        print(f"{nm:16s}{d['fx1']:>12.5f}{d['fx0']:>12.5f}")
    print(f"\nFully-homozygous-X likelihood ratio  FD:MS:SS = "
          f"{ana['father-daughter']['fx1']/ana['brother-sister']['fx1']:.0f}"
          f":{ana['mother-son']['fx1']/ana['brother-sister']['fx1']:.0f}:1")

    print("\nValidation vs gene-drop simulation (exact events, Wilson 95% CI):")
    v = validate_against_sim(n=300000)
    print(f"  (n = {v['n']:,} per union)")
    print(f"{'union':16s}{'event':8s}{'analytic':>10}{'sim':>10}"
          f"{'95% CI':>20}{'  ok':>5}")
    for r in v["rows"]:
        ci = f"[{r['ci_lo']:.5f},{r['ci_hi']:.5f}]"
        print(f"{r['union']:16s}{r['event']:8s}{r['analytic']:>10.5f}{r['sim']:>10.5f}"
              f"{ci:>20}{'  yes' if r['inside'] else '  NO':>5}")
