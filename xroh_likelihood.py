"""
xroh_likelihood.py — exact map-based likelihood of the union type from a female
X-chromosome autozygosity (ROH) track.
==========================================================================
Given an observed autozygosity track A(x) on a female X (1 = the two X copies are
IBD/homozygous-by-descent at x, 0 otherwise), compute the likelihood of each
first-degree hypothesis directly from the empirical recombination map — no Monte
Carlo. This turns "which union?" into a closed-form likelihood ratio / posterior.

GENERATIVE MODEL (females recombine the non-PAR X; males pass it intact).
Discretize the X into bins; let p_i = 1/2(1 - e^{-2 mu_i}) be the probability
that ONE female meiosis switches parental homolog across bin i (mu_i = local
genetic length in Morgans). Each hypothesis maps the latent homolog telegraph(s)
onto the observed autozygosity track:

  FD (father-daughter): child's paternal X = grandfather's X, intact; the
    maternal X is a single-meiosis recombinant. A(x) is ONE symmetric telegraph
    that flips with prob p_i per bin. P(A=1) = 1/2  -> mean F_X = 1/2.

  MS (mother-son): both of the child's X's are independent single-meiosis
    recombinants of the SAME mother. A = (U1==U2) flips when EITHER flips, i.e.
    with prob q_i = 2 p_i (1-p_i) per bin. P(A=1)=1/2 -> mean F_X = 1/2.
    => L_FD/L_MS reduces to the junction COUNT (the map cancels); FD vs MS is
       the hard pair (~0.71 ceiling), unaffected by the spatial map.

  SS (sib-sib): father's X is one GM meiosis (V_B); mother carries the
    grandfather's X intact (g) plus one GM meiosis (V_Z), and transmits a
    recombinant of {g, V_Z} (selector W_Z). A=1 iff W_Z picks the GM side AND
    V_B==V_Z. P(A=1)=1/4 -> mean F_X = 1/4. A is NOT a simple telegraph, so the
    likelihood is an 8-state HMM over (V_B, V_Z, W_Z) -> the spatial map does NOT
    cancel. This is why FD vs SS (the forensically important pair) is separable.

Run:  python xroh_likelihood.py        # validate vs simulation, both maps
Deps: numpy, plus xroh_sim.py in the same folder.
"""
from __future__ import annotations
import os
import numpy as np
import xroh_sim as S


# --------------------------------------------------------------------------
# Maps -> (pos, cumulative_cM, P0, P1, L). Bherer female X (via xroh_sim) and the
# deCODE maternal chrX built from the UCSC recombMat track (GRCh38).
# --------------------------------------------------------------------------
def bherer_map():
    S.init()
    return S.XMAP


def decode_map(path="data/decode_chrX_maternal_hg38.txt"):
    import pandas as pd
    d = pd.read_csv(path, sep="\t")
    pos = d["pos"].to_numpy(float)
    cm = d["cM"].to_numpy(float) + np.arange(len(d)) * 1e-9   # strict monotone
    return pos, cm, pos[0], pos[-1], cm[-1]


# --------------------------------------------------------------------------
# Binning + per-bin single-meiosis switch probability from the map.
# --------------------------------------------------------------------------
def make_bins(mp, bin_bp=0.5e6):
    POS, CM, P0, P1, L = mp
    edges = np.arange(P0, P1 + bin_bp, bin_bp)
    edges[-1] = P1
    cm_edges = np.interp(edges, POS, CM)
    mu = np.diff(cm_edges) / 100.0                 # Morgans per bin
    p = 0.5 * (1.0 - np.exp(-2.0 * mu))            # 1-meiosis switch prob / bin
    mids = 0.5 * (edges[:-1] + edges[1:])
    return dict(edges=edges, mids=mids, p=p, n=len(p))


# --------------------------------------------------------------------------
# Observed autozygosity track on the bin grid from a (pat, mat) mosaic pair.
# --------------------------------------------------------------------------
def autozygosity_track(pat, mat, bins):
    mids = bins["mids"]

    def lab(h):
        st = np.array([s for s, _ in h]); lb = np.array([l for _, l in h])
        return lb[np.searchsorted(st, mids, "right") - 1]

    return (lab(pat) == lab(mat)).astype(np.int8)


# --------------------------------------------------------------------------
# Likelihoods.  All return log L of the observed track A under the hypothesis.
# FD / MS: closed-form telegraph.  SS: 8-state HMM forward (log-scaled).
# --------------------------------------------------------------------------
def _telegraph_loglik(A, flip_p):
    """log L of a binary track A whose state flips with prob flip_p[i] per bin
    and is observed directly (stationary symmetric start)."""
    flips = A[1:] != A[:-1]
    fp = np.clip(flip_p[1:], 1e-12, 1 - 1e-12)
    ll = np.log(0.5)                       # initial state prob
    ll += np.sum(np.where(flips, np.log(fp), np.log1p(-fp)))
    return ll


def loglik_FD(A, bins):
    return _telegraph_loglik(A, bins["p"])


def loglik_MS(A, bins):
    p = bins["p"]
    q = 2.0 * p * (1.0 - p)                # A flips if exactly one of 2 meioses flips
    return _telegraph_loglik(A, q)


# SS latent state = (V_B, V_Z, W_Z) in {0,1}^3, indexed 0..7 as bits b2 b1 b0.
_SS_STATES = np.arange(8)
_VB = (_SS_STATES >> 2) & 1
_VZ = (_SS_STATES >> 1) & 1
_WZ = _SS_STATES & 1
_SS_A = ((_WZ == 1) & (_VB == _VZ)).astype(np.int8)   # autozygous iff GM-side & V_B==V_Z
# Hamming distance between states (# of the 3 telegraphs that must flip)
_SS_HAM = (((_SS_STATES[:, None] ^ _SS_STATES[None, :]) >> np.array([0, 1, 2])[:, None, None]) & 1).sum(0)


def loglik_SS(A, bins):
    p = bins["p"]
    H = _SS_HAM                                   # (8,8) flips needed s'->s
    alpha = np.where(_SS_A == A[0], 1.0 / 8.0, 0.0)
    ll = 0.0
    for i in range(1, len(A)):
        pi = min(max(p[i], 1e-12), 1 - 1e-12)
        T = (pi ** H) * ((1 - pi) ** (3 - H))     # independent flips of 3 telegraphs
        nxt = alpha @ T
        nxt = np.where(_SS_A == A[i], nxt, 0.0)
        z = nxt.sum()
        if z <= 0:
            return -np.inf
        alpha = nxt / z
        ll += np.log(z)
    return ll


def posterior(A, bins, priors=(1 / 3, 1 / 3, 1 / 3)):
    lls = np.array([loglik_FD(A, bins), loglik_MS(A, bins), loglik_SS(A, bins)])
    logp = lls + np.log(np.array(priors))
    logp -= logp.max()
    post = np.exp(logp); post /= post.sum()
    return dict(zip(["FD", "MS", "SS"], post)), dict(zip(["FD", "MS", "SS"], lls))


# --------------------------------------------------------------------------
# Validation: simulate genomes, classify by the analytic likelihood, report
# confusion + the pairwise accuracies (esp. FD vs SS, the important pair).
# --------------------------------------------------------------------------
def validate(mp, name, n=2000, bin_bp=0.5e6, seed=2024):
    S.reset(seed)
    # the simulators use xroh_sim's XMAP; point it at this map for gene-drop
    S.XMAP = mp
    bins = make_bins(mp, bin_bp)
    sims = {"FD": S.FD, "MS": S.MS, "SS": S.SS}
    LL = {k: [] for k in sims}          # per true class: list of (llFD,llMS,llSS)
    for cls, fn in sims.items():
        for _ in range(n):
            A = autozygosity_track(*fn(), bins)
            LL[cls].append([loglik_FD(A, bins), loglik_MS(A, bins), loglik_SS(A, bins)])
        LL[cls] = np.array(LL[cls])

    names = ["FD", "MS", "SS"]
    print(f"\n=== Likelihood classifier on {name} map "
          f"({bins['n']} bins x {bin_bp/1e6:.2f} Mb, n={n}/class) ===")

    # 3-way (equal priors)
    correct = 0
    conf = np.zeros((3, 3))
    for ti, cls in enumerate(names):
        pred = LL[cls].argmax(1)
        for pj in range(3):
            conf[ti, pj] = np.mean(pred == pj)
        correct += np.mean(pred == ti)
    print(f"3-way accuracy = {correct/3:.3f}   (chance 0.333)")
    print("  rows=truth cols=pred ", names)
    for ti, cls in enumerate(names):
        print(f"   {cls}  {conf[ti].round(3)}")

    # pairwise (restrict to the two relevant log-liks)
    def pair(a, b):
        ia, ib = names.index(a), names.index(b)
        # truth a classified as a if ll_a > ll_b
        acc_a = np.mean(LL[a][:, ia] > LL[a][:, ib])
        acc_b = np.mean(LL[b][:, ib] > LL[b][:, ia])
        return 0.5 * (acc_a + acc_b)

    print(f"  pairwise  FD vs SS = {pair('FD','SS'):.3f}   "
          f"FD vs MS = {pair('FD','MS'):.3f}   MS vs SS = {pair('MS','SS'):.3f}")
    return LL


if __name__ == "__main__":
    bh = bherer_map()
    validate(bh, "Bherer female X (GRCh37)", n=2000)

    if os.path.exists("data/decode_chrX_maternal_hg38.txt"):
        dc = decode_map()
        print(f"\n[deCODE maternal chrX loaded: {dc[4]:.2f} cM, "
              f"{dc[2]/1e6:.1f}-{dc[3]/1e6:.1f} Mb]")
        validate(dc, "deCODE maternal X (GRCh38)", n=2000)
    else:
        print("\n(deCODE map file not found; run the fetch step first.)")
