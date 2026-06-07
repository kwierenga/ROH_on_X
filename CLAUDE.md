# CLAUDE.md — agent orientation

Research project on **X-chromosome runs of homozygosity (ROH) and the type of consanguineous union**. Read `PROJECT_HANDOFF.md` first for the full picture.

## TL;DR of the science

In a **female** offspring, X autozygosity (F_X) depends on the *sex-path* of the consanguinity loop, because a father transmits his X intact to daughters and sons inherit no X from their father. So F_X distinguishes union types that are identical on the autosomes (all first-degree unions have autosomal F = ¼, but F_X = ½/½/¼ for father–daughter/mother–son/sib–sib).

## How to run

```bash
pip install -r requirements.txt
python xroh_sim.py        # prints catalogue + detection + Bayes ceiling
```

First run auto-downloads the Bhérer 2017 female X recombination map into `data/`.

## Where things are

- `xroh_sim.py` — the engine (gene-drop, pedigree catalogue, SNP+ROH detection, Bayes ceiling). All results derive from here.
- `PROJECT_HANDOFF.md` — full context, results, data provenance, TODO.
- `X_ROH_incest_type_methods_note.md` — Paper A (clinical/forensic: which first-degree union?).
- `X_autozygosity_mating_type_fingerprint.md` — Paper B (population: marriage-system fingerprint).
- `*.png` — figures.

## Key facts the agent must respect

- **Validation/unit test:** simulated mean F_X must converge to the analytic values (½, ½, ¼ for first degree; the {0,1,2,3}× catalogue ratios). If it doesn’t, the X/sex transmission logic is broken.
- **Network in the original sandbox was restricted** (no 1000G VCFs); the SNP layer is *simulated*. Here you have open network — swap in real 1000G/HGDP haplotypes where useful.
- **ROH calling:** use the built-in PLINK-style `call_roh`, or real tools (PLINK `--homozyg`, `bcftools roh`). Avoid `scikit-allel.roh_mhmm` (OOMs at chromosome scale).
- **Don’t fabricate** numbers or citations. Verify references (several are flagged “verify” in both notes). When citing PubMed articles, include DOIs.

## Top open tasks (see PROJECT_HANDOFF §6)

1. Add crossover interference; re-check the Bayes ceiling.
1. Add a genotyped parent to test whether it breaks father–daughter vs mother–son.
1. Swap in real 1000G/HGDP haplotypes for the detection study.
1. Paper B biobank scan with a sex-biased-demography control.
1. SLiM forward model under a fixed marriage rule.