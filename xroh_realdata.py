"""
xroh_realdata.py — validate the X-ROH engine on REAL 1000 Genomes haplotypes
============================================================================
The headline accuracies in `xroh_sim.py` are gene-drops on *true* autozygosity
tracks with a simulated SNP layer (no LD). This module checks that the result
survives real human haplotype structure: it drops REAL phased 1000G chrX
haplotypes through the FD/MS/SS pedigrees (so the union type is ground-truth),
genotypes the child at real SNP positions, and calls ROH with the validated
PLINK-style caller (`xroh_sim.call_roh`) under real LD and allele frequencies.

No plink/bcftools/tabix needed. The panel is built by STREAMING the position-
sorted phased VCF and stopping early (we only need a chrX region), so we never
download the whole 1.9 GB file.

Usage:
  # 1. build a real-haplotype panel over a chrX region (streamed, stops at END):
  curl -s <1000G chrX phased VCF URL> | python xroh_realdata.py build
  # 2. validate:
  python xroh_realdata.py validate
"""
from __future__ import annotations
import sys, gzip, numpy as np

PANEL_FILE = "data/chrX_realpanel.npz"
SAMPLES = "data/1kg_samples.panel"
START, END = 2_700_000, 25_000_000     # non-PAR chrX region (GRCh37)
THIN = 3_000                           # min bp between kept SNPs
MAF_MIN = 0.10
N_FEMALES = 200                        # -> 400 real haplotypes in the founder pool
MAX_SNP = 8000


def _female_samples():
    names = []
    with open(SAMPLES) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[3] == "female":
                names.append(p[0])
    return names[:N_FEMALES]


def build():
    """Stream a phased VCF from stdin; write a thinned real-haplotype panel."""
    want = set(_female_samples())
    gz = gzip.GzipFile(fileobj=sys.stdin.buffer)
    cols = None
    pos_list, hap_rows = [], []
    last = -10 ** 9
    kept = 0
    for raw in gz:
        line = raw.decode("ascii", "replace")
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            hdr = line.rstrip("\n").split("\t")
            cols = [i for i, s in enumerate(hdr) if s in want]
            continue
        t1 = line.index("\t"); t2 = line.index("\t", t1 + 1)
        pos = int(line[t1 + 1:t2])
        if pos < START:
            continue
        if pos > END or kept >= MAX_SNP:
            break
        if pos - last < THIN:
            continue
        p = line.split("\t")
        ref, alt = p[3], p[4]
        if len(ref) != 1 or len(alt) != 1:      # biallelic SNP only
            continue
        haps = []
        for c in cols:
            gt = p[c]
            # phased diploid female GT "a|b"; take both haplotypes
            a, b = gt[0], gt[2] if len(gt) >= 3 else gt[0]
            haps.append(1 if a == "1" else 0)
            haps.append(1 if b == "1" else 0)
        h = np.array(haps, dtype=np.int8)
        maf = h.mean()
        if maf < MAF_MIN or maf > 1 - MAF_MIN:
            continue
        pos_list.append(pos); hap_rows.append(h); last = pos; kept += 1
    H = np.array(hap_rows, dtype=np.int8).T          # n_hap x n_snp
    pos = np.array(pos_list, dtype=np.int64)
    np.savez_compressed(PANEL_FILE, pos=pos, H=H)
    print(f"panel: {H.shape[0]} real haplotypes x {H.shape[1]} SNPs over "
          f"{pos.min():,}-{pos.max():,} bp  -> {PANEL_FILE}")


def validate(n=150, seed=8):
    import xroh_sim as X
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.model_selection import cross_val_score
    X.init()
    d = np.load(PANEL_FILE)
    pos, H = d["pos"], d["H"]
    nsnp, nhap = len(pos), H.shape[0]
    span = float(pos.max() - pos.min())
    ar = np.arange(nsnp)
    X.reset(seed)

    def one(fn):
        pat, mat = fn()
        lp = X.labels_at(pat, pos); lm = X.labels_at(mat, pos)
        uniq = set(lp.tolist()) | set(lm.tolist())
        hmap = {l: int(X.RNG.integers(nhap)) for l in uniq}      # founder -> real haplotype
        a1 = H[[hmap[l] for l in lp], ar]
        a2 = H[[hmap[l] for l in lm], ar]
        het = a1 != a2
        froh, nseg = X.call_roh(pos.astype(int), het, span=span)
        true_reg = float(np.mean(lp == lm))                      # true autozygous SNP fraction (region)
        # longest detected ROH (cM-free, Mb) via a quick pass on het runs
        return true_reg, froh, nseg

    res, feats, ys = {}, [], []
    for cls, (nm, fn) in enumerate([("father-daughter", X.FD), ("mother-son", X.MS),
                                    ("brother-sister", X.SS)]):
        tr, fr, sg = [], [], []
        for _ in range(n):
            t, f, s = one(fn)
            tr.append(t); fr.append(f); sg.append(s)
            feats.append([f, s]); ys.append(cls)
        res[nm] = dict(true_region=float(np.mean(tr)), detected_froh=float(np.mean(fr)),
                       n_seg=float(np.mean(sg)))
    feats = np.array(feats); ys = np.array(ys)

    def acc(a, b):
        m = (ys == a) | (ys == b)
        return float(cross_val_score(GradientBoostingClassifier(max_depth=3, n_estimators=120),
                                     feats[m], (ys[m] == b).astype(int), cv=5).mean())

    print(f"Real 1000G chrX haplotypes: {nhap} haps x {nsnp} SNPs over "
          f"{pos.min():,}-{pos.max():,} bp ({span/1e6:.1f} Mb), n={n}/union\n")
    print(f"{'union':16s}{'true F(region)':>15}{'detected F_roh':>16}{'#ROH':>7}")
    for nm, s in res.items():
        print(f"{nm:16s}{s['true_region']:>15.3f}{s['detected_froh']:>16.3f}{s['n_seg']:>7.2f}")
    print("\nReal-data ROH classification (detected features only, 5-fold CV):")
    print(f"  FD vs SS = {acc(0,2):.3f}   MS vs SS = {acc(1,2):.3f}   FD vs MS = {acc(0,1):.3f}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "validate"
    (build if mode == "build" else validate)()
