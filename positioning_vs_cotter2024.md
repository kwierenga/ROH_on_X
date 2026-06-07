# Positioning vs. Cotter et al. 2024 (G3) — prior art for the X-ROH program

**Reference.** Cotter DJ, Severson AL, Kang JTL, Godrej HN, Carmi S, Rosenberg NA.
*Modeling the effects of consanguinity on autosomal and X-chromosomal runs of
homozygosity and identity-by-descent sharing.* **G3** 2024;14(2):jkad264.
doi:10.1093/g3journal/jkad264. PDF: `docs/cotter-rosenberg-2024-x-auto-roh-consanguinity.pdf`.

This is the closest prior art to **Paper B** (the population / marriage-system
paper) and must be cited and engaged. It does **not** touch **Paper A** (the
single-genome clinical/forensic union-typing paper). Below: exactly what they
own, where we overlap, and the defensible white space — updated with our own
new simulation results.

---

## 1. What Cotter et al. own (do not re-derive — cite)

- A **diploid coalescent model** giving the **expected population-mean** fraction
  of the X vs the autosomes in ROH and IBD, as a function of consanguinity rates,
  via TMRCA distributions (Severson 2019; Cotter 2021, 2022) plus a
  recombination/segment-length layer (Palamara 2012).
- The consanguinity axis is exactly our sex-path axis: the **four first-cousin
  types** — patrilateral-parallel (c_pp), patrilateral-cross (c_pc),
  matrilateral-parallel (c_mp), matrilateral-cross (c_mc) (their Fig. 1).
- Three headline predictions: (1) X-IBD rises with X-ROH; (2) **even with no
  consanguinity the X:autosome ratio ≈ 2** (smaller N_e,X, lower coalescence,
  ~2/3 recombination — their eq. 5); (3) **matrilateral consanguinity pushes the
  ratio above 2 (parallel > cross); patrilateral pushes it below 2**, and
  "patrilateral consanguinity produces no ROH on the X chromosome." The
  consanguinity effect is **larger for ROH than for IBD** (their Fig. 2).
- An **empirical check**: 13 Jewish populations, X genotypes from Behar 2013
  (1,647 individuals, **only 420 female**, 13,052 X SNPs), demographic cousin-type
  rates from Goldschmidt 1960 (Israeli births 1955–57). Result: each 1% rise in
  autosomal ROH ↔ **2.1% rise in X ROH** (close to the predicted ~2×).

So Cotter et al. own the **coalescent expectation of population-mean X/autosome
ROH-and-IBD as a function of first-cousin-marriage type**, qualitatively
validated in aggregate.

---

## 2. What they explicitly do NOT do (our white space)

Read against the paper's own Methods/Discussion:

1. **No single-individual inference.** Everything is a population *mean* total ROH
   proportion (regressions *across populations*). There is no posterior/likelihood
   for "what union produced **this** genome." → **Paper A is entirely open.**
2. **No first-degree or non-cousin loops.** Only the four first-cousin types, and
   only at infinite-N coalescent limits. → Our catalogue covers parent–child,
   full-sib, half-sib, avuncular, double-first-cousin, where the sex-path signal
   is **strongest** (X:auto ratio {0,1,2,3}).
3. **Expectations only — no realized distribution.** They give means; single-genome
   typing needs the **full sampling distribution and its tail** (see §3). They do
   not model the variance at all.
4. **No ROH-*calling* / detection study.** No treatment of SNP density, genotyping
   error, or caller behavior on the X. Their empirical X-ROH used **autosome-
   borrowed LOD cutoffs** (they flag this as uncertain), **no X length classes**,
   and a small, **female-poor** (n=420), dated panel. → A modern, **X-specific
   ROH-calling biobank scan** is unclaimed.
5. **No disentangling of the N_e,X / sex-biased-demography baseline from the
   consanguinity signal** within real data (their baseline ≈2 is assumed
   universal; Cox-Hammer 2008 shows N_e,X is strongly population-specific). →
   The identifiability problem (polygyny vs matrilateral marriage both raise X:A)
   is open.

---

## 3. Our distinct, already-demonstrated contributions

From `xroh_sim.py` (gene-drop on the Bhérer female X map; all simulation-confirmed):

- **Realized single-genome distribution & its tail** (Cotter give only means).
  The fully-homozygous-X tail is strongly type-specific and is itself a
  discriminator: P(F_X > 0.99) = **9.3% (father–daughter)** vs **2.0% (mother–son)**
  vs **0.2% (sib–sib)** under Haldane. (Mechanism: father–daughter needs only one
  non-recombinant maternal transmission; this is also the explanation for the
  Sund 2013 "entirely homozygous X" girl.)
- **Single-genome union-typing ceiling and how to break it** (Paper A — no Cotter
  analog). FD vs MS from one genome ≈ coin flip (full-info ceiling **0.69**); a
  **genotyped birth mother raises it to 0.91**, because the child's paternal X is
  an *intact* maternal homolog under FD but a *recombinant mosaic* of the mother's
  two homologs under MS (observable maternal-homolog switch count: ~0 vs ~Poisson(L)).
- **Recombination-model sensitivity.** Adding gamma-renewal crossover interference
  (ν) leaves the mean F_X at ½/½/¼ (unit test intact) but **shrinks the FD
  fully-homozygous tail** (P(F_X>0.99): 0.093 → 0.034 at ν=2.6 → 0.014 at ν=4.3)
  and **modestly raises** the FD/MS single-genome ceiling (0.69 → 0.73 → 0.76).
  So both the headline tail probability and the ceiling are interference-dependent
  — a methodological point Cotter's coalescent treatment cannot make.

---

## 4. One-paragraph "how we differ" (drop-in for both drafts)

> Cotter et al. (2024) derived the expected population-mean ratio of X-chromosomal
> to autosomal ROH and IBD under the four types of first-cousin marriage in a
> diploid coalescent framework, and confirmed the predicted ~2× X:autosome ROH
> relationship in aggregate genotype data from Jewish populations. Our work is
> complementary and non-overlapping in three respects. First, where Cotter et al.
> model population means, we characterize the **realized single-genome
> distribution** of X autozygosity — in particular its type-specific tail — which
> is the quantity relevant to clinical and forensic inference about an individual
> (Paper A). Second, we extend beyond first cousins to the **full catalogue of
> first- and second-degree loops**, where the sex-path signal is strongest and
> where forensic questions arise. Third, we add an explicit **ROH-detection layer**
> (SNP density, genotyping error, X-specific calling) and a recombination-
> interference model, enabling a real-data biobank scan with X-specific thresholds
> and an explicit treatment of the N_e,X/sex-biased-demography baseline that Cotter
> et al. assume as a fixed factor of two.

---

## 5. Action items this creates

- [ ] Paper B intro/discussion: cite Cotter 2024 as the population-level precedent;
      use the paragraph in §4. Drop the handoff §5 claim that "nobody has used ROH
      on the X to discriminate the type of union" — replace with the realized-vs-
      expected / single-genome / detection framing.
- [ ] Paper B method: commit to **X-specific LOD/length thresholds** (not
      autosome-borrowed), X length classes, and a population-calibrated X:A null —
      i.e., do the three things Cotter explicitly could not.
- [ ] Verify the KinInbcoefX/X-kinship lineage (Thornton 2012; Grossman & Eisen
      1989) and Campbell 2015 for the interference ν before final ν choice.
