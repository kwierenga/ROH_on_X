"""
test_xroh.py — regression tests locking the manuscript invariants to code.

Runs under pytest OR standalone (`python test_xroh.py`), so it needs no extra
deps. Every assertion corresponds to a number claimed in union_typing_paper.md.
"""
import numpy as np
import xroh_sim as X
import xroh_tails as T

X.init()


def test_means_converge_to_half_half_quarter():
    """Unit test: simulated F_X -> 1/2, 1/2, 1/4 (the X/sex transmission is correct)."""
    X.reset(1)
    for fn, target in [(X.FD, 0.5), (X.MS, 0.5), (X.SS, 0.25)]:
        fx = np.mean([X.Fx(*fn()) for _ in range(4000)])
        assert abs(fx - target) < 0.02, (fn.__name__, fx, target)


def test_catalogue_ratios_0123():
    """X/autosome ratio takes values {0,1,2,3} by sex-path (paper §4.1, catalogue)."""
    df = X.run_catalogue(n=3000, seed=2)
    r = dict(zip(df["configuration"], df["ratio"]))
    assert abs(r["father-daughter"] - 2) < 0.15
    assert abs(r["brother-sister"] - 1) < 0.15
    assert r["paternal half-sibs"] < 0.10            # X=0 (loop reaches father via his father)
    assert r["fathers are brothers"] < 0.10
    assert abs(r["aunt = father's-mother's sister"] - 3) < 0.4   # maximal matrilineal


def test_analytic_tails_closed_form():
    """Exact Haldane tails (paper §4.2.1): 1/2 e^-L, 1/2 e^-2L, 1/4 e^-3L, SS survival."""
    a = T.analytic_tails(176.3)
    assert abs(a["father-daughter"]["fx1"] - 0.0858) < 0.002
    assert abs(a["mother-son"]["fx1"] - 0.0147) < 0.002
    assert abs(a["brother-sister"]["fx1"] - 0.00126) < 0.0005
    assert abs(a["brother-sister"]["fx0"] - 0.1785) < 0.010


def test_tails_match_simulation_within_wilson_ci():
    """Every closed-form tail lies inside the simulation's Wilson 95% CI (§4.2.1)."""
    v = T.validate_against_sim(n=80000, seed=5)
    for row in v["rows"]:
        assert row["inside"], row


def test_genotyped_mother_breaks_fd_ms():
    """Genotyped mother lifts FD-vs-MS above the single-child ceiling (paper §4.4)."""
    g = X.genotyped_mother_experiment(n=1500, seed=3)
    assert g["child_only"] < g["with_mother"]
    assert g["with_mother"] > 0.85


def test_wilson_ci_basic():
    lo, hi = X.wilson_ci(5, 1000)
    assert 0.0 < lo < 0.005 < hi < 0.02


def test_fx_decomposition_extremes_only_in_all_zero_corner():
    """Whole-X extremes live in the k=0 corner; SS fully-homozygous ~ 0 (§4.2.1)."""
    d = X.fx_by_crossover_count(n=40000, seed=4)
    assert d["brother-sister"]["overall"]["p_fx1"] < 0.004
    assert d["father-daughter"][0]["p"] > 0.15          # P(k=0) ~ 0.17


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    fails = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
        except AssertionError as e:
            fails += 1
            print(f"FAIL  {fn.__name__}: {e}")
    print(f"\n{len(fns) - fails}/{len(fns)} passed")
    raise SystemExit(1 if fails else 0)
