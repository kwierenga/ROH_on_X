"""
xroh_sim.py — X-chromosomal autozygosity simulation toolkit
============================================================
Gene-dropping engine for runs-of-homozygosity (ROH) on the X chromosome in
incestuous / consanguineous matings, using empirical sex-specific recombination
maps. Reproduces all results in the two methods notes.

Core idea: in a FEMALE offspring the paternal X is transmitted intact (males do
not recombine the non-PAR X and pass no X to sons), so X autozygosity depends on
the SEX-PATH of the consanguinity loop, not just its depth. F_X differs by union
type even when autosomal F is identical.

Run:  python xroh_sim.py            # prints headline results + catalogue
Deps: numpy, pandas, scipy, matplotlib, scikit-learn   (see requirements.txt)

Data: female X recombination map from Bhérer et al. 2017 (auto-downloaded from
GitHub on first run). See data-provenance notes in PROJECT_HANDOFF.md.
"""
from __future__ import annotations
import os, tarfile, urllib.request
import numpy as np
import pandas as pd


def wilson_ci(k, n, z=1.96):
    """Wilson score 95% interval for a binomial proportion — correct in the
    small-Np tail where the normal approximation fails. Returns (lo, hi)."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))

# --------------------------------------------------------------------------
# Data: Bhérer et al. 2017 refined sex-specific map (GRCh37). chrX is FEMALE-only
# (no male/sexavg X exists, because the non-PAR X recombines only in females).
# --------------------------------------------------------------------------
BHERER_URL = ("https://raw.githubusercontent.com/cbherer/"
              "Bherer_etal_SexualDimorphismRecombination/master/"
              "Refined_genetic_map_b37.tar.gz")
DATA_DIR = os.environ.get("XROH_DATA", "data")
MAP_DIR = os.path.join(DATA_DIR, "Refined_genetic_map_b37")


def ensure_maps():
    """Download + extract the Bhérer map if not already present."""
    if os.path.isdir(MAP_DIR):
        return
    os.makedirs(DATA_DIR, exist_ok=True)
    tgz = os.path.join(DATA_DIR, "bherer_b37.tar.gz")
    if not os.path.exists(tgz):
        print("Downloading Bhérer map ...")
        urllib.request.urlretrieve(BHERER_URL, tgz)
    with tarfile.open(tgz) as t:
        t.extractall(DATA_DIR)


def load_map(fname):
    """Return (pos, cumulative_cM, P0, P1, L_cM) for one chromosome file."""
    d = pd.read_csv(os.path.join(MAP_DIR, fname), sep="\t")
    pos = d["pos"].to_numpy(float)
    cm = d["cM"].to_numpy(float) + np.arange(len(d)) * 1e-9  # strict monotone
    return pos, cm, pos[0], pos[-1], cm[-1]


# Module-level maps (set by init()).
XMAP = None          # female X map tuple
AUT = {}             # autosome -> {"M":map, "F":map}


def init(autosomes=False):
    ensure_maps()
    global XMAP, AUT
    XMAP = load_map("female_chrX.txt")
    if autosomes:
        AUT = {c: {"M": load_map(f"male_chr{c}.txt"),
                   "F": load_map(f"female_chr{c}.txt")} for c in range(1, 23)}
    return XMAP


# --------------------------------------------------------------------------
# Core gene-drop primitives. A haplotype mosaic = list of (start_bp, label),
# covering [P0, P1]; each unique integer label = one ancestral founder X.
# --------------------------------------------------------------------------
RNG = np.random.default_rng()
_LABEL = [0]

# Crossover interference. None or 1.0 => Haldane/Poisson (no interference, the
# validated default). nu > 1 => gamma-renewal (Housworth-Stahl) model: inter-
# crossover gaps ~ Gamma(shape=nu, mean=1 Morgan), so expected crossover count
# (= genetic length) is UNCHANGED and the mean F_X is preserved (½/½/¼ unit
# test still holds); only the *spacing variance* tightens. Set a literature nu
# (human ~ 2.6–7; verify against Campbell et al. 2015 before final use).
INTERFERENCE_NU = None


def reset(seed=None):
    global RNG
    RNG = np.random.default_rng(seed)
    _LABEL[0] = 0


def _xover_cM(L_cM):
    """Crossover positions (cM) on a gamete over a chromosome of length L_cM.

    Haldane (nu None/1): homogeneous Poisson. Otherwise a stationary gamma
    renewal process with shape nu and mean gap 1 Morgan (burn-in from -10 M so
    the boundary at 0 is forgotten). Both have expected count = L_cM/100.
    """
    nu = INTERFERENCE_NU
    LM = L_cM / 100.0
    if nu is None or nu == 1:
        k = RNG.poisson(LM)
        return np.sort(RNG.uniform(0, L_cM, k)) if k > 0 else np.empty(0)
    xs, pos = [], -10.0
    while True:
        pos += RNG.gamma(nu, 1.0 / nu)
        if pos > LM:
            break
        if pos > 0:
            xs.append(pos * 100.0)
    return np.asarray(xs)


def founder(mp):
    """A new founder haplotype (single mosaic) on map mp."""
    _LABEL[0] += 1
    return [(mp[2], _LABEL[0])]


def meiosis(hA, hB, mp):
    """Recombine two mosaics under map mp (interference per INTERFERENCE_NU)."""
    POS, CM, P0, P1, L = mp
    cm_x = _xover_cM(L)
    xs = sorted(np.interp(cm_x, CM, POS)) if len(cm_x) else []
    bnd = [P0] + list(xs) + [P1]
    par = [hA, hB]
    s = int(RNG.integers(2))
    seg = []
    for j in range(len(bnd) - 1):
        a, b = bnd[j], bnd[j + 1]
        h = par[(s + j) % 2]
        for i, (st, lab) in enumerate(h):
            e = h[i + 1][0] if i + 1 < len(h) else P1
            if e <= a or st >= b:
                continue
            if not seg or seg[-1][1] != lab:
                seg.append((max(st, a), lab))
    if seg[0][0] > P0:
        seg[0] = (P0, seg[0][1])
    return seg


def autozygosity(pat, mat, mp):
    """Fraction of the chromosome where the two copies share a founder label."""
    P0, P1 = mp[2], mp[3]
    bks = sorted(set([s for s, _ in pat] + [s for s, _ in mat] + [P1]))

    def lab(h, x):
        c = h[0][1]
        for s, l in h:
            if s <= x:
                c = l
            else:
                break
        return c

    t = 0.0
    for j in range(len(bks) - 1):
        a, b = bks[j], bks[j + 1]
        if lab(pat, a) == lab(mat, a):
            t += b - a
    return t / (P1 - P0)


# X-specific shorthands (XMAP must be init'd).
def Hx():               return founder(XMAP)
def son(motherpair):    return meiosis(motherpair[0], motherpair[1], XMAP)  # male's single X
def dau(fX, motherpair): return (list(fX), son(motherpair))                 # daughter's X pair
def female_child(fX, motherpair): return list(fX), son(motherpair)          # (O_pat, O_mat)
def Fx(pat, mat):       return autozygosity(pat, mat, XMAP)


# --------------------------------------------------------------------------
# Pedigree library. Each returns (O_pat, O_mat) for a FEMALE child.
# autoF = textbook autosomal inbreeding coefficient for that union.
# --------------------------------------------------------------------------
def FD():  ff = Hx(); W = (Hx(), Hx()); return female_child(ff, dau(ff, W))
def MS():  M = (Hx(), Hx()); return female_child(son(M), M)
def SS():  GF = Hx(); GM = (Hx(), Hx()); return female_child(son(GM), dau(GF, GM))

def matHS(): W = (Hx(), Hx()); return female_child(son(W), dau(Hx(), W))
def patHS(): H = Hx(); return female_child(son((Hx(), Hx())), dau(H, (Hx(), Hx())))

def unc_mom_bro():  GGF = Hx(); GGM = (Hx(), Hx()); return female_child(son(GGM), dau(Hx(), dau(GGF, GGM)))
def unc_dad_bro():  GGF = Hx(); GGM = (Hx(), Hx()); return female_child(son(GGM), dau(son(GGM), (Hx(), Hx())))
def aunt_via_PGM(): GGF = Hx(); GGM = (Hx(), Hx()); PGM = dau(GGF, GGM); return female_child(son(PGM), dau(GGF, GGM))
def aunt_via_PGF(): GGF = Hx(); GGM = (Hx(), Hx()); return female_child(son((Hx(), Hx())), dau(GGF, GGM))

def cous_moms_sis(): GGF = Hx(); GGM = (Hx(), Hx()); return female_child(son(dau(GGF, GGM)), dau(Hx(), dau(GGF, GGM)))
def cous_dads_bro(): GGF = Hx(); GGM = (Hx(), Hx()); return female_child(son((Hx(), Hx())), dau(son(GGM), (Hx(), Hx())))
def cous_PGM_MGF():  GGF = Hx(); GGM = (Hx(), Hx()); return female_child(son(dau(GGF, GGM)), dau(son(GGM), (Hx(), Hx())))
def cous_PGF_MGM():  GGF = Hx(); GGM = (Hx(), Hx()); _ = son(GGM); return female_child(son((Hx(), Hx())), dau(Hx(), dau(GGF, GGM)))

CATALOGUE = [
    ("First degree",      "father-daughter",                 FD,           0.25),
    ("First degree",      "mother-son",                      MS,           0.25),
    ("First degree",      "brother-sister",                  SS,           0.25),
    ("Half-sib (2nd)",    "maternal half-sibs",              matHS,        0.125),
    ("Half-sib (2nd)",    "paternal half-sibs",              patHS,        0.125),
    ("Avuncular (2nd)",   "uncle = mother's brother",        unc_mom_bro,  0.125),
    ("Avuncular (2nd)",   "uncle = father's brother",        unc_dad_bro,  0.125),
    ("Avuncular (2nd)",   "aunt = father's-mother's sister", aunt_via_PGM, 0.125),
    ("Avuncular (2nd)",   "aunt = father's-father's sister", aunt_via_PGF, 0.125),
    ("First cousin (3rd)", "mothers are sisters",            cous_moms_sis, 0.0625),
    ("First cousin (3rd)", "fathers are brothers",           cous_dads_bro, 0.0625),
    ("First cousin (3rd)", "PGM & MGF are sibs",             cous_PGM_MGF,  0.0625),
    ("First cousin (3rd)", "PGF & MGM are sibs",             cous_PGF_MGM,  0.0625),
]


def run_catalogue(n=6000, seed=2024):
    reset(seed)
    rows = []
    for cls, name, fn, fa in CATALOGUE:
        fx = float(np.mean([Fx(*fn()) for _ in range(n)]))
        rows.append((cls, name, fa, fx, (fx / fa if fa else 0.0)))
    return pd.DataFrame(rows, columns=["class", "configuration", "F_auto", "F_X", "ratio"])


# --------------------------------------------------------------------------
# SNP-genotype layer + PLINK-style sliding-window ROH caller (detection study).
# NOTE: 1000G VCFs were not reachable in the build sandbox, so the SNP panel is
# simulated at realistic density. In VS Code (open network) you can swap in real
# founder haplotypes from 1000G/HGDP.
# --------------------------------------------------------------------------
def labels_at(hap, snp):
    st = np.array([s for s, _ in hap]); lb = np.array([l for _, l in hap])
    return lb[np.searchsorted(st, snp, "right") - 1]


def call_roh(snp, het, W=50, maxhet=2, minlen=1.5e6, minSNP=50, span=None):
    """PLINK --homozyg-style caller. Returns (froh, n_segments)."""
    if span is None:
        span = XMAP[3] - XMAP[2]
    n = len(snp)
    c = np.concatenate([[0], np.cumsum(het.astype(int))])
    ends = np.minimum(np.arange(n) + W, n)
    win_pass = (c[ends] - c[np.arange(n)]) <= maxhet
    cp = np.concatenate([[0], np.cumsum(win_pass.astype(int))])
    lo = np.maximum(np.arange(n) - W + 1, 0)
    covered = (cp[np.arange(n) + 1] - cp[lo]) > 0
    segs, i = [], 0
    while i < n:
        if covered[i]:
            j = i
            while j + 1 < n and covered[j + 1]:
                j += 1
            a, b = snp[i], snp[j]
            if (b - a) >= minlen and (j - i + 1) >= minSNP:
                segs.append((a, b))
            i = j + 1
        else:
            i += 1
    return sum(b - a for a, b in segs) / span, len(segs)


def detection_study(ped_fn, k_founders, n=150, nsnp=12000, eps=0.005, seed=8):
    """Compare true vs detected F_X for one pedigree at a given SNP density."""
    reset(seed)
    P0, P1 = XMAP[2], XMAP[3]
    snp = np.sort(RNG.choice(np.arange(int(P0) + 1, int(P1) - 1), nsnp, replace=False)).astype(int)
    af = RNG.uniform(0.05, 0.5, nsnp)
    tru, det, nseg = [], [], []
    ar = np.arange(nsnp)
    for _ in range(n):
        pat, mat = ped_fn()
        A = (RNG.random((k_founders + 2, nsnp)) < af).astype(np.int8)  # +slack on labels
        l1 = labels_at(pat, snp) % A.shape[0]
        l2 = labels_at(mat, snp) % A.shape[0]
        a1, a2 = A[l1, ar], A[l2, ar]
        het = a1 != a2
        flip = (RNG.random(nsnp) < eps) & (~het)
        het = het | flip
        fr, ns = call_roh(snp, het, span=P1 - P0)
        tru.append(Fx(pat, mat)); det.append(fr); nseg.append(ns)
    return dict(true=np.mean(tru), detected=np.mean(det), segs=np.mean(nseg))


# --------------------------------------------------------------------------
# Bayes ceiling for father-daughter vs mother-son (the within-class hard pair).
# Junctions ~ Poisson(lambda) [FD] vs Poisson(2*lambda) [MS], lambda = L/100.
# --------------------------------------------------------------------------
def bayes_ceiling_fd_ms():
    from scipy.stats import poisson
    L = XMAP[4]; lam = L / 100.0
    ks = np.arange(0, 60)
    acc_count = 0.5 * np.sum(np.maximum(poisson.pmf(ks, lam), poisson.pmf(ks, 2 * lam)))
    return dict(lam=lam, acc_junction_count=float(acc_count))


def bayes_ceiling_fullinfo(n=12000, seed=21):
    """Full-information ceiling via a classifier on simulated noise-free features."""
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.model_selection import cross_val_score
    reset(seed)
    L = XMAP[4]

    def feats(pat, mat):
        # IBD segment lengths (cM) + summary stats
        bks = sorted(set([s for s, _ in pat] + [s for s, _ in mat] + [XMAP[3]]))

        def lab(h, x):
            c = h[0][1]
            for s, l in h:
                if s <= x: c = l
                else: break
            return c
        segs, cur = [], None
        for j in range(len(bks) - 1):
            a, b = bks[j], bks[j + 1]
            if lab(pat, a) == lab(mat, a):
                cur = [a, b] if cur is None else [cur[0], b]
            elif cur is not None:
                segs.append(tuple(cur)); cur = None
        if cur is not None:
            segs.append(tuple(cur))
        cl = [np.interp(b, XMAP[0], XMAP[1]) - np.interp(a, XMAP[0], XMAP[1]) for a, b in segs]
        return [len(segs), sum(cl) / L, (max(cl) if cl else 0), (np.std(cl) if len(cl) > 1 else 0)]

    X = [feats(*FD()) for _ in range(n)] + [feats(*MS()) for _ in range(n)]
    y = [0] * n + [1] * n
    acc = cross_val_score(GradientBoostingClassifier(max_depth=3, n_estimators=120),
                          np.array(X), np.array(y), cv=5).mean()
    return float(acc)


# --------------------------------------------------------------------------
# F_X tail statistics (the fully-homozygous-X tail differs sharply by type).
# --------------------------------------------------------------------------
def fx_tail_stats(n=40000, seed=7):
    reset(seed)
    out = {}
    for nm, fn in [("father-daughter", FD), ("mother-son", MS), ("brother-sister", SS)]:
        fx = np.array([Fx(*fn()) for _ in range(n)])
        out[nm] = dict(mean=float(fx.mean()), p_gt99=float(np.mean(fx > 0.99)),
                       p_gt95=float(np.mean(fx > 0.95)), p_lt05=float(np.mean(fx < 0.05)))
    return out


def fx_by_crossover_count(n=60000, seed=1):
    """Decompose F_X by the per-meiosis crossover count of each union.

    FD is governed by ONE female-X meiosis (maternal transmission, count k), MS by
    TWO meioses of the same mother (k1+k2), SS by THREE (the GF/grandmother mixing
    count ks is the structural switch). The whole-chromosome outcomes (F_X=1 fully
    homozygous, F_X=0 fully outbred) live almost entirely in the all-zero corner:
    FD needs e^{-λ}, MS needs e^{-2λ}, and SS can essentially never be F_X=1 (the
    grandfather's X is a one-sided intrusion that blocks homozygosity). Returns
    {union: {k: dict(p, mean, p_fx1, p_fx0)}} with the governing count as key."""
    POS, CM, P0, P1, Lm = XMAP

    def cmeiosis(hA, hB):                    # meiosis() that also returns crossover count
        cm_x = _xover_cM(Lm); k = len(cm_x)
        xs = sorted(np.interp(cm_x, CM, POS)) if k else []
        bnd = [P0] + list(xs) + [P1]; par = [hA, hB]; s = int(RNG.integers(2)); seg = []
        for j in range(len(bnd) - 1):
            a, b = bnd[j], bnd[j + 1]; h = par[(s + j) % 2]
            for i, (st, lab) in enumerate(h):
                e = h[i + 1][0] if i + 1 < len(h) else P1
                if e <= a or st >= b: continue
                if not seg or seg[-1][1] != lab: seg.append((max(st, a), lab))
        if seg[0][0] > P0: seg[0] = (P0, seg[0][1])
        return seg, k

    def FD_c():                              # governing count: maternal meiosis km
        G = Hx(); W = (Hx(), Hx()); Wson, _ = cmeiosis(*W)
        Omat, km = cmeiosis(G, Wson); return Fx(G, Omat), km

    def MS_c():                              # governing count: k1+k2 (same mother)
        M = (Hx(), Hx()); Op, k1 = cmeiosis(*M); Om, k2 = cmeiosis(*M); return Fx(Op, Om), k1 + k2

    def SS_c():                              # governing count: ks (GF/grandmother mixing)
        GF = Hx(); GM = (Hx(), Hx()); Op, _ = cmeiosis(*GM); U, _ = cmeiosis(*GM)
        Om, ks = cmeiosis(GF, U); return Fx(Op, Om), ks

    reset(seed)
    tol = 1e-6                               # exact extremes, not >0.99 (avoids telomere slivers)
    out = {}

    def cell(f, n_tot):
        n_c = len(f)
        k1 = int(np.sum(f >= 1 - tol)); k0 = int(np.sum(f <= tol))
        return dict(p=n_c / n_tot, n=n_c, mean=float(f.mean()),
                    p_fx1=k1 / n_c if n_c else 0.0, ci_fx1=wilson_ci(k1, n_c),
                    p_fx0=k0 / n_c if n_c else 0.0, ci_fx0=wilson_ci(k0, n_c))

    for nm, fn in [("father-daughter", FD_c), ("mother-son", MS_c), ("brother-sister", SS_c)]:
        fxv = np.empty(n); kv = np.empty(n, int)
        for i in range(n):
            fxv[i], kv[i] = fn()
        d = {}
        for k in range(0, 5):
            m = kv == k
            if m.any():
                d[k] = cell(fxv[m], n)
        m = kv >= 5
        if m.any():
            d["5+"] = cell(fxv[m], n)
        d["overall"] = cell(fxv, n)
        out[nm] = d
    return out


# --------------------------------------------------------------------------
# Genotyped-mother experiment: does observing the birth mother's X break the
# father-daughter vs mother-son degeneracy? Mechanism: the child's PATERNAL X
# equals an INTACT maternal homolog under father-daughter (father = mother's
# father, transmits X undivided), but is a RECOMBINANT MOSAIC of the mother's
# two homologs under mother-son (father = mother's son, got a recombined X from
# her). So the number of maternal-homolog switches in the child's paternal X is
# ~0 (FD) vs ~Poisson(L) (MS) — an observable that child-only data lack.
# --------------------------------------------------------------------------
def FD_duo():
    """Father-daughter; returns (child_pat, child_mat, mother_hapA, mother_hapB)."""
    ff = Hx(); W = (Hx(), Hx())
    D = dau(ff, W)                       # birth mother = father's daughter
    pat, mat = female_child(ff, D)
    return pat, mat, D[0], D[1]


def MS_duo():
    """Mother-son; returns (child_pat, child_mat, mother_hapA, mother_hapB)."""
    M = (Hx(), Hx())                     # birth mother
    pat, mat = female_child(son(M), M)
    return pat, mat, M[0], M[1]


def _ibd_segments(pat, mat):
    bks = sorted(set([s for s, _ in pat] + [s for s, _ in mat] + [XMAP[3]]))

    def lab(h, x):
        c = h[0][1]
        for s, l in h:
            if s <= x: c = l
            else: break
        return c
    segs, cur = [], None
    for j in range(len(bks) - 1):
        a, b = bks[j], bks[j + 1]
        if lab(pat, a) == lab(mat, a):
            cur = [a, b] if cur is None else [cur[0], b]
        elif cur is not None:
            segs.append(tuple(cur)); cur = None
    if cur is not None:
        segs.append(tuple(cur))
    return segs


def seg_features(pat, mat):
    """Child-only X features: [n_segments, F_X, max_seg_cM, std_seg_cM]."""
    L = XMAP[4]
    segs = _ibd_segments(pat, mat)
    cl = [np.interp(b, XMAP[0], XMAP[1]) - np.interp(a, XMAP[0], XMAP[1]) for a, b in segs]
    return [len(segs), sum(cl) / L, (max(cl) if cl else 0.0),
            (float(np.std(cl)) if len(cl) > 1 else 0.0)]


def _mother_features(pat, mhapA, mhapB):
    """Express the child's paternal X in terms of the two maternal homologs:
    [n_switches, homolog_purity, fraction_traced, is_intact(0/1)]."""
    SA = {l for _, l in mhapA}; SB = {l for _, l in mhapB}
    P0, P1 = XMAP[2], XMAP[3]
    sw, prev, lenA, lenB, unk = 0, None, 0.0, 0.0, 0.0
    for i, (st, l) in enumerate(pat):
        e = pat[i + 1][0] if i + 1 < len(pat) else P1
        idx = 0 if l in SA else (1 if l in SB else -1)
        if idx == 0: lenA += e - st
        elif idx == 1: lenB += e - st
        else: unk += e - st
        if idx != -1:
            if prev is not None and idx != prev: sw += 1
            prev = idx
    tot = lenA + lenB
    purity = (max(lenA, lenB) / tot) if tot > 0 else 1.0
    return [sw, purity, tot / (P1 - P0), 1.0 if sw == 0 else 0.0]


def genotyped_mother_experiment(n=4000, seed=23):
    """Full-info ceiling for FD vs MS, child-only vs child+genotyped-mother."""
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.model_selection import cross_val_score
    reset(seed)
    Xc, Xf, y = [], [], []
    for cls, fn in [(0, FD_duo), (1, MS_duo)]:
        for _ in range(n):
            pat, mat, mA, mB = fn()
            child = seg_features(pat, mat)
            Xc.append(child)
            Xf.append(child + _mother_features(pat, mA, mB))
            y.append(cls)
    y = np.array(y)

    def acc(X):
        return float(cross_val_score(
            GradientBoostingClassifier(max_depth=3, n_estimators=120),
            np.array(X), y, cv=5).mean())
    return dict(child_only=acc(Xc), with_mother=acc(Xf))


# --------------------------------------------------------------------------
# Multiple offspring. An incestuous union often yields several children. Full
# sisters share the father's X *identically* (a father transmits his single X
# intact to every daughter), but each gets an INDEPENDENT maternal recombination.
# So multiple daughters do not add independent paternal-X observations, but they
# collectively reveal the mother's two X homologs — which is what breaks FD vs MS
# (under FD the shared paternal X equals one intact maternal homolog; under MS it
# is a recombinant of both). Accuracy should climb from the single-child ceiling
# toward the genotyped-mother ceiling as the sibship grows.
# --------------------------------------------------------------------------
def FD_family(k):
    ff = Hx(); W = (Hx(), Hx()); M = (ff, son(W))      # mother = grandfather's daughter
    return [female_child(ff, M) for _ in range(k)]     # shared paternal ff, independent maternal

def MS_family(k):
    M = (Hx(), Hx()); fX = son(M)                      # father = the son; his X fixed
    return [female_child(fX, M) for _ in range(k)]

def SS_family(k):
    GF = Hx(); GM = (Hx(), Hx()); fX = son(GM); Z = (GF, son(GM))   # father = brother, mother = sister
    return [female_child(fX, Z) for _ in range(k)]


def _family_xfeatures(kids):
    """Permutation-invariant X summary over a sibship (each kid's [F_X, #seg, max_cM])."""
    fx, ns, mx = [], [], []
    for pat, mat in kids:
        sf = seg_features(pat, mat)                    # [n_seg, F_X, max_seg_cM, std_seg_cM]
        fx.append(sf[1]); ns.append(sf[0]); mx.append(sf[2])
    fx, ns, mx = np.array(fx), np.array(ns), np.array(mx)
    return [fx.mean(), fx.std(), fx.min(), fx.max(), ns.mean(), ns.max(), mx.mean(), mx.max()]


def multiple_offspring_typing(ks=(1, 2, 3, 5, 8), n=1200, seed=31):
    """Single-genome -> sibship: X-only FD/MS/SS accuracy as a function of the
    number of female children of the same union (shared paternal X, independent
    maternal X). Shows how many sisters it takes to break the FD-vs-MS hard pair."""
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.model_selection import cross_val_score
    fams = {"father-daughter": FD_family, "mother-son": MS_family, "brother-sister": SS_family}
    out = {}
    for k in ks:
        reset(seed)
        feat = {nm: np.array([_family_xfeatures(fn(k)) for _ in range(n)]) for nm, fn in fams.items()}

        def acc(a, b):
            Xm = np.vstack([feat[a], feat[b]]); y = np.array([0] * n + [1] * n)
            return float(cross_val_score(GradientBoostingClassifier(max_depth=3, n_estimators=120),
                                         Xm, y, cv=5).mean())
        out[k] = dict(FD_vs_MS=acc("father-daughter", "mother-son"),
                      FD_vs_SS=acc("father-daughter", "brother-sister"),
                      MS_vs_SS=acc("mother-son", "brother-sister"))
    return out


# --------------------------------------------------------------------------
# Autosomal gene-drop (diploid, sex-specific maps) + whole-genome union typing.
# All three first-degree unions have autosomal F=1/4, but the ROH SEGMENT
# architecture orders them FD < MS < SS: parent-child loops span g=3 meioses,
# sib-sib g=4 (shorter, more numerous ROH), and within g=3 the MS loop routes
# more FEMALE meioses (the female map is ~1.6x longer) so it fragments slightly
# more than FD. Requires init(autosomes=True). Reproduces the paper's autosomal
# (#ROH FD~26 < MS~30 < SS~35) and whole-genome (FD-vs-SS ~0.87) numbers.
# --------------------------------------------------------------------------
def _seg_ibd(pat, mat, P0, P1):
    """IBD segments (bp) where two diploid copies share a founder label."""
    bks = sorted(set([s for s, _ in pat] + [s for s, _ in mat] + [P1]))

    def lab(h, x):
        c = h[0][1]
        for s, l in h:
            if s <= x: c = l
            else: break
        return c
    segs, cur = [], None
    for j in range(len(bks) - 1):
        a, b = bks[j], bks[j + 1]
        if lab(pat, a) == lab(mat, a):
            cur = [a, b] if cur is None else [cur[0], b]
        elif cur is not None:
            segs.append(tuple(cur)); cur = None
    if cur is not None:
        segs.append(tuple(cur))
    return segs


def auto_FD(mpM, mpF):                       # common ancestor = grandfather G; g=3
    g1, g2 = founder(mpM), founder(mpM); w1, w2 = founder(mpF), founder(mpF)
    Mpat = meiosis(g1, g2, mpM); Mmat = meiosis(w1, w2, mpF)   # mother = G x W
    return meiosis(g1, g2, mpM), meiosis(Mpat, Mmat, mpF)      # child: from G, from mother


def auto_MS(mpM, mpF):                       # common ancestor = mother M; g=3
    f1, f2 = founder(mpF), founder(mpF); u = founder(mpM)      # mother M, unrelated grandfather
    Spat, Smat = u, meiosis(f1, f2, mpF)                       # son = u x M
    return meiosis(Spat, Smat, mpM), meiosis(f1, f2, mpF)      # child: from son, from M


def auto_SS(mpM, mpF):                       # common ancestors = both grandparents; g=4
    a1, a2 = founder(mpM), founder(mpM); b1, b2 = founder(mpF), founder(mpF)
    Bpat, Bmat = meiosis(a1, a2, mpM), meiosis(b1, b2, mpF)    # brother
    Zpat, Zmat = meiosis(a1, a2, mpM), meiosis(b1, b2, mpF)    # sister
    return meiosis(Bpat, Bmat, mpM), meiosis(Zpat, Zmat, mpF)  # child: from brother, from sister


_AUTO_PED = {"father-daughter": auto_FD, "mother-son": auto_MS, "brother-sister": auto_SS}
_X_PED = {"father-daughter": FD, "mother-son": MS, "brother-sister": SS}


def genome_features(name, minlen=1.5e6):
    """One female child's whole-genome feature vector
    [#autoROH, autoF, mean_ROH_Mb, max_ROH_Mb, F_X, X_nseg, X_max_seg_cM].
    The autosomes carry loop depth (indices 0-3); the X carries the sex-path via
    its autozygosity LEVEL (F_X, idx 4) and SEGMENT STRUCTURE (count + longest run,
    idx 5-6) — so the joint classifier uses the full X signal, not F_X alone.
    Requires init(autosomes=True)."""
    nroh, tot, span, lens = 0, 0.0, 0.0, []
    for c in range(1, 23):
        mpM, mpF = AUT[c]["M"], AUT[c]["F"]
        P0, P1 = mpF[2], mpF[3]; span += P1 - P0
        cp, cm = _AUTO_PED[name](mpM, mpF)
        for a, b in _seg_ibd(cp, cm, P0, P1):
            tot += b - a
            if b - a >= minlen:
                nroh += 1; lens.append(b - a)
    xp, xm = _X_PED[name]()
    xseg = seg_features(xp, xm)              # [n_seg, F_X, max_seg_cM, std_seg_cM]
    return [nroh, tot / span, (np.mean(lens) / 1e6 if lens else 0.0),
            (max(lens) / 1e6 if lens else 0.0), xseg[1], xseg[0], xseg[2]]


def autosomal_roh_by_union(n=250, seed=2024):
    """Mean autosomal #ROH(>=1.5 Mb), length, autosomal F and F_X per union, with
    the SD and a normal 95% CI on the mean #ROH."""
    reset(seed)
    out = {}
    for nm in _AUTO_PED:
        r = np.array([genome_features(nm) for _ in range(n)])
        nr = r[:, 0]; sd = float(nr.std(ddof=1)); se = sd / np.sqrt(n)
        out[nm] = dict(nroh=float(nr.mean()), nroh_sd=sd,
                       nroh_ci95=(float(nr.mean() - 1.96 * se), float(nr.mean() + 1.96 * se)),
                       autoF=float(r[:, 1].mean()), meanMb=float(r[:, 2].mean()),
                       fx=float(r[:, 4].mean()))
    return out


def _acc_se(Xm, y, cv=5):
    """5-fold CV accuracy with its cross-fold standard error (mean, SE)."""
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.model_selection import cross_val_score
    s = cross_val_score(GradientBoostingClassifier(max_depth=3, n_estimators=120), Xm, y, cv=cv)
    return float(s.mean()), float(s.std(ddof=1) / np.sqrt(len(s)))


def whole_genome_typing(n=600, seed=77):
    """FD-vs-SS single-genome accuracy from X only, autosomes only, and combined
    (5-fold CV gradient boosting, each with cross-fold SE). init(autosomes=True)."""
    reset(seed)
    feats = np.array([genome_features("father-daughter") for _ in range(n)] +
                     [genome_features("brother-sister") for _ in range(n)])
    y = np.array([0] * n + [1] * n)
    def run(idx):
        a, se = _acc_se(feats[:, idx], y)
        return dict(acc=a, se=se)
    return dict(x_only=run([4, 5, 6]), auto_only=run([0, 1, 2, 3]),
                whole_genome=run([0, 1, 2, 3, 4, 5, 6]))


def whole_genome_all_pairs(n=600, seed=77):
    """All three pairwise + the 3-way single-genome accuracies on the full
    whole-genome feature vector [#autoROH, autoF, meanMb, maxMb, F_X], each with a
    5-fold CV cross-fold SE. Replaces the d'-in-quadrature approximation for the
    MS-vs-SS and FD-vs-MS combined figures. Requires init(autosomes=True)."""
    reset(seed)
    feats = {nm: np.array([genome_features(nm) for _ in range(n)])
             for nm in ("father-daughter", "mother-son", "brother-sister")}
    out = {}
    for a, b, na, nb in [("FD", "MS", "father-daughter", "mother-son"),
                         ("FD", "SS", "father-daughter", "brother-sister"),
                         ("MS", "SS", "mother-son", "brother-sister")]:
        acc, se = _acc_se(np.vstack([feats[na], feats[nb]]), np.array([0] * n + [1] * n))
        out[f"{a}_vs_{b}"] = dict(acc=acc, se=se)
    acc, se = _acc_se(np.vstack([feats[k] for k in feats]), np.array([0] * n + [1] * n + [2] * n))
    out["three_way"] = dict(acc=acc, se=se)
    return out


# --------------------------------------------------------------------------
if __name__ == "__main__":
    init()
    print("Female X map length:", round(XMAP[4], 1), "cM\n")

    print("=== Catalogue: X autozygosity by mating type (female offspring) ===")
    df = run_catalogue(n=4000)
    print(df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    print("\n=== Detection (true vs detected F_X) ===")
    for nm, fn, k in [("sib-sib", SS, 3), ("father-daughter", FD, 3), ("mother-son", MS, 2)]:
        r = detection_study(fn, k, n=60, nsnp=12000)
        print(f"{nm:16s} true={r['true']:.3f} detected={r['detected']:.3f} segs={r['segs']:.2f}")

    print("\n=== Bayes ceiling: father-daughter vs mother-son ===")
    bc = bayes_ceiling_fd_ms()
    print(f"lambda={bc['lam']:.2f}  acc(junction count)={bc['acc_junction_count']:.3f}")
    print(f"acc(full info, sim+GBM)={bayes_ceiling_fullinfo(n=4000):.3f}  (chance=0.50)")

    print("\n=== F_X tail by type (the fully-homozygous-X signature) ===")
    for nm, s in fx_tail_stats(n=40000).items():
        print(f"{nm:16s} mean={s['mean']:.3f}  P(Fx>0.99)={s['p_gt99']:.3f}  "
              f"P(Fx>0.95)={s['p_gt95']:.3f}  P(Fx<0.05)={s['p_lt05']:.3f}")

    print("\n=== F_X decomposed by governing crossover count (exact extremes, Wilson 95% CI) ===")
    for nm, d in fx_by_crossover_count(n=60000).items():
        z = d.get(0, {}); ov = d["overall"]
        lo, hi = ov["ci_fx1"]
        print(f"{nm:16s} k=0: P={z.get('p',0):.3f} P(Fx=1|k0)={z.get('p_fx1',0):.3f}  |  "
              f"overall P(Fx=1)={ov['p_fx1']:.4f} [{lo:.4f},{hi:.4f}]  "
              f"P(Fx=0)={ov['p_fx0']:.4f}")

    print("\n=== Genotyped-mother experiment: FD vs MS (chance=0.50) ===")
    gm = genotyped_mother_experiment(n=4000)
    print(f"child only      acc={gm['child_only']:.3f}")
    print(f"+ genotyped mom acc={gm['with_mother']:.3f}")

    print("\n=== Whole genome: autosomal ROH by union (mean #ROH +/- 95% CI) ===")
    init(autosomes=True)
    for nm, s in autosomal_roh_by_union(n=200).items():
        lo, hi = s['nroh_ci95']
        print(f"{nm:16s} #ROH={s['nroh']:.1f} [{lo:.1f},{hi:.1f}]  meanMb={s['meanMb']:.1f}  "
              f"autoF={s['autoF']:.3f}  F_X={s['fx']:.3f}")
    wg = whole_genome_typing(n=500)
    print(f"FD-vs-SS  X-only={wg['x_only']['acc']:.3f}  autosomes={wg['auto_only']['acc']:.3f}  "
          f"combined={wg['whole_genome']['acc']:.3f}+/-{wg['whole_genome']['se']:.3f}")
    print("Joint classifier, all pairs + 3-way (acc +/- CV SE):")
    for k, v in whole_genome_all_pairs(n=500).items():
        print(f"  {k:10s} {v['acc']:.3f} +/- {v['se']:.3f}")

    print("\n=== Crossover-interference sensitivity (nu) ===")
    print(f"{'nu':>6} {'FD_meanFx':>10} {'FD_P(Fx>.99)':>13} {'FD/MS_fullinfo':>15}")
    for nu in (None, 2.6, 4.3):
        INTERFERENCE_NU = nu
        t = fx_tail_stats(n=20000)["father-daughter"]
        full = bayes_ceiling_fullinfo(n=4000)
        print(f"{str(nu):>6} {t['mean']:>10.3f} {t['p_gt99']:>13.3f} {full:>15.3f}")
    INTERFERENCE_NU = None
