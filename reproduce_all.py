"""
reproduce_all.py — regenerate every headline number in union_typing_paper.md
from code with fixed seeds, so the manuscript is mechanically linked to the
engine. Run: `python reproduce_all.py`. (For the real-data check see
xroh_realdata.py; for the full per-module output run each module's __main__.)
"""
import numpy as np
import xroh_sim as X
import xroh_tails as T


def main():
    X.init(autosomes=True)
    print("=" * 64)
    print("union_typing_paper.md — headline numbers reproduced")
    print("=" * 64)

    print("\n[§4.1] mean F_X (target 1/2, 1/2, 1/4):")
    X.reset(1)
    for nm, fn in [("FD", X.FD), ("MS", X.MS), ("SS", X.SS)]:
        print(f"   {nm} = {np.mean([X.Fx(*fn()) for _ in range(4000)]):.3f}")

    print("\n[§4.2.1] exact fully-homozygous-X tails P(F_X=1) and P(F_X=0):")
    a = T.analytic_tails(X.XMAP[4])
    for nm in a:
        print(f"   {nm:16s} fx1={a[nm]['fx1']:.5f}  fx0={a[nm]['fx0']:.5f}")
    fr = a["father-daughter"]["fx1"] / a["brother-sister"]["fx1"]
    mr = a["mother-son"]["fx1"] / a["brother-sister"]["fx1"]
    print(f"   fully-homozygous LR  FD:MS:SS = {fr:.0f}:{mr:.0f}:1")

    print("\n[§4.2.1] tails validated against simulation (Wilson 95% CI):")
    v = T.validate_against_sim(n=120000)
    print(f"   all {len(v['rows'])} analytic values inside CI: "
          f"{all(r['inside'] for r in v['rows'])}")

    print("\n[§4.3] autosomal #ROH by union (target FD<MS<SS):")
    for nm, s in X.autosomal_roh_by_union(n=200).items():
        lo, hi = s["nroh_ci95"]
        print(f"   {nm:16s} #ROH={s['nroh']:.1f} [{lo:.1f},{hi:.1f}]  autoF={s['autoF']:.3f}")

    print("\n[§4.4] whole-genome joint classifier (acc +/- CV SE):")
    for k, val in X.whole_genome_all_pairs(n=500).items():
        print(f"   {k:10s} {val['acc']:.3f} +/- {val['se']:.3f}")

    print("\n[§4.4] genotyped-mother lift for FD-vs-MS:")
    g = X.genotyped_mother_experiment(n=3000)
    print(f"   child only = {g['child_only']:.3f}   + genotyped mother = {g['with_mother']:.3f}")

    print("\nDone. (Sund-case posterior: see xroh_tails; calibration: xroh_posterior;")
    print(" real 1000G haplotypes: xroh_realdata.)")


if __name__ == "__main__":
    main()
