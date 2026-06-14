"""
xroh_posterior.py — calibrated union-type posteriors + the maternal prior
=========================================================================
Forensic use needs a *calibrated* likelihood ratio / posterior, not just an
accuracy: when the method says "P(father-daughter) = 0.9", that should be right
~90% of the time. This module (i) turns the exact X likelihoods of
`xroh_likelihood.py` into posteriors, (ii) measures their calibration (expected
calibration error, ECE) including under model mismatch (data simulated WITH
crossover interference but scored with the Haldane likelihood — the realistic
worst case), and (iii) folds in the §6.1 maternal prior: because a mother-son
mother must already have a reproductive-age son, a YOUNG birth mother sets the
MS prior ~0, collapsing the problem to the well-separated FD-vs-SS pair.
"""
from __future__ import annotations
import numpy as np
import xroh_sim as S
import xroh_likelihood as L

NAMES = ["FD", "MS", "SS"]


def simulate_lls(mp, n=2000, nu_data=None, seed=7, bin_bp=0.5e6):
    """n cases/class: gene-drop (optionally with interference nu_data) -> the three
    Haldane log-likelihoods. Returns (truth_idx array, lls array n*3 x 3)."""
    S.XMAP = mp
    bins = L.make_bins(mp, bin_bp)
    old = S.INTERFERENCE_NU
    S.INTERFERENCE_NU = nu_data
    S.reset(seed)
    sims = [S.FD, S.MS, S.SS]
    truth, lls = [], []
    for ci, fn in enumerate(sims):
        for _ in range(n):
            A = L.autozygosity_track(*fn(), bins)
            lls.append([L.loglik_FD(A, bins), L.loglik_MS(A, bins), L.loglik_SS(A, bins)])
            truth.append(ci)
    S.INTERFERENCE_NU = old
    return np.array(truth), np.array(lls)


def posteriors(lls, priors=(1 / 3, 1 / 3, 1 / 3)):
    lp = np.log(np.array(priors) + 1e-300)
    x = lls + lp
    x = x - x.max(1, keepdims=True)
    p = np.exp(x)
    return p / p.sum(1, keepdims=True)


def ece(truth, post, n_bins=10):
    """Expected calibration error + reliability table on the top-class confidence."""
    pred = post.argmax(1); conf = post.max(1); correct = (pred == truth).astype(float)
    edges = np.linspace(0, 1, n_bins + 1)
    e, table = 0.0, []
    for i in range(n_bins):
        m = (conf >= edges[i]) & (conf <= edges[i + 1] if i == n_bins - 1 else conf < edges[i + 1])
        if m.sum() == 0:
            continue
        acc, cf, w = correct[m].mean(), conf[m].mean(), m.mean()
        e += w * abs(acc - cf)
        table.append((cf, acc, int(m.sum())))
    return e, table


def acc(truth, post):
    return float((post.argmax(1) == truth).mean())


if __name__ == "__main__":
    mp = L.bherer_map()
    print("Calibration of the exact-likelihood union posterior (3-way, equal priors)\n")

    # (a) matched: Haldane data, Haldane score
    t, ll = simulate_lls(mp, n=2500, nu_data=None, seed=11)
    post = posteriors(ll)
    e, tab = ece(t, post)
    print(f"(a) matched (Haldane/Haldane):   accuracy={acc(t,post):.3f}  ECE={e:.3f}")

    # (b) mismatch: interference data (nu=4.3), Haldane score
    t2, ll2 = simulate_lls(mp, n=2500, nu_data=4.3, seed=12)
    post2 = posteriors(ll2)
    e2, tab2 = ece(t2, post2)
    print(f"(b) mismatch (nu=4.3 data, Haldane score): accuracy={acc(t2,post2):.3f}  ECE={e2:.3f}")
    print("    reliability (confidence -> empirical accuracy):")
    for cf, ac, k in tab2:
        print(f"      conf~{cf:.2f}  acc={ac:.2f}  (n={k})")

    # (c) maternal prior: young mother => P(MS)=0  (collapses to FD vs SS)
    print("\nMaternal prior (§6.1): a young birth mother sets P(MS)=0.")
    for label, pri in [("flat prior        ", (1/3, 1/3, 1/3)),
                       ("young mother (MS=0)", (0.5, 0.0, 0.5))]:
        p = posteriors(ll, priors=pri)
        # restrict evaluation to the FD/SS cases that remain possible under this prior
        keep = (t == 0) | (t == 2) if pri[1] == 0 else np.ones(len(t), bool)
        print(f"  {label}: 3-way acc={acc(t,p):.3f}   FD-vs-SS subset acc={acc(t[keep],p[keep]):.3f}")
