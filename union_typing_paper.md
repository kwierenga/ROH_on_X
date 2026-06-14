# Autosomes give the loop depth, the X gives the sex-path: resolving the *type* of first-degree consanguinity from a single female genome

*Methods note — v1.0 (2026-06-13). All numbers are simulation results reproduced by `xroh_sim.py` and `xroh_likelihood.py` in this repository, on empirical sex-specific recombination maps (Bhérer 2017; deCODE/Halldorsson 2019). References verified against PubMed (PMIDs/DOIs given). This note extends and supersedes the earlier X-only draft (`X ROH incest type methods note.md`).*

---

## Abstract

Genome-wide runs of homozygosity (ROH) reliably reveal that a child’s parents are first-degree relatives, but the autosomal inbreeding coefficient is **F ≈ ¼ for all three first-degree unions** — father–daughter (FD), mother–son (MS), and brother–sister (SS) — so the autosomes cannot, on their face, say *which* union produced the child. We show that the genome nonetheless carries this information in two complementary places. (1) On the **autosomes**, the ROH *segment architecture* differs because the consanguinity loop spans a different number of meioses (g = 3 for parent–child, g = 4 for sib–sib) and a different mix of male and female meioses; this orders the unions by ROH count (FD < MS < SS) even at identical F. (2) On the **X chromosome** of a female offspring, autozygosity depends on the *sex-path* of the loop, because a father transmits his X intact to daughters and sons inherit no X from their father: expected X autozygosity is **½ (FD), ½ (MS), ¼ (SS)**. We formalise the X signal as an exact, closed-form likelihood read directly from the recombination map (no Monte Carlo), prove that the FD-vs-MS likelihood ratio collapses to the junction *count* (the spatial map cancels), and validate the whole framework against gene-dropping simulation on two independent maps. Combining autosomes and X in a single joint classifier, single-genome accuracies are **≈0.87 (FD vs SS), ≈0.77 (MS vs SS), and ≈0.74 (FD vs MS)**. The forensically and clinically decisive contrast — father–daughter versus sibling — is the *best* resolved; father–daughter versus mother–son is the irreducible hard core, breakable in practice only by genotyping a relative (which lifts it to ≈0.91). The discriminating signal lives partly on the one chromosome clinical pipelines routinely discard.

---

## 1. Background

SNP/chromosomal microarray and exome/genome sequencing are first-tier tests for developmental delay, intellectual disability and congenital anomalies, and they incidentally reveal parental relatedness through genome-wide ROH. Standard practice estimates the *degree* of relatedness from the autosomal ROH burden (~25% of the genome autozygous → first-degree) and, following ACMG guidance, deliberately does **not** assign the specific parental relationship, because array data are not a paternity test and total ROH burden maps imperfectly to a named relationship (Grote 2012; Fan 2013; Bennett 2021).

This conservatism is well founded for the autosomes alone, where all first-degree unions share F ≈ ¼. But two features of inheritance break the symmetry the autosomal *burden* hides:

- **Loop depth and meiosis sex (autosomes).** The autozygous tracts in an inbred individual are fragments of an ancestral chromosome broken by recombination accumulated over the meioses of the consanguinity loop. The number of those meioses, and whether each was male or female (female maps are ~1.6× longer), shapes the ROH *length spectrum* even when total F is fixed.
- **X transmission asymmetry (sex chromosome).** A father passes his single X **intact** to every daughter (no recombination); a son inherits **no** X from his father; a mother passes a recombined X. So in a female child, X autozygosity depends on the **sexes along the consanguinity loop**, not only its depth.

We quantify both, give the X signal an exact likelihood, and show the two sources are complementary.

That the sex-path leaves an X signature is established at the **population** scale: Cotter et al. (2024, *G3*) showed in a diploid coalescent model that, across the four first-cousin marriage types, **matrilateral** consanguinity raises the population-mean X-to-autosome ROH ratio while **patrilateral** consanguinity lowers it (“patrilateral consanguinity produces no ROH on the X chromosome”), and confirmed a ~2× X:autosome ROH relationship in aggregate human data. The present work moves the same sex-path principle from the **population mean of cousin marriages** to the **realised single genome of a first-degree union**, where the signal is strongest and where the clinical/forensic question of *type* actually arises.

---

## 2. Theory

### 2.1 The X: sex-path determines F_X

For a female child, the **paternal X is transmitted intact** and derives entirely from the paternal grandmother; the **maternal X** is a single-meiosis recombinant of the mother’s two X’s. Tracing the three first-degree pedigrees gives the expected X inbreeding coefficient (= the X-kinship coefficient of the parents; McPeek framework; Grossman & Eisen 1989):

| Union (female offspring) | Autosomal F | Expected F_X |
|---|---|---|
| Father × daughter (FD) | ¼ | **½** |
| Mother × son (MS) | ¼ | **½** |
| Brother × sister (SS) | ¼ | **¼** |

Mechanistically, the child’s X autozygosity track A(x) (1 = the two X copies are identical-by-descent at x) is:

- **FD:** one symmetric telegraph — the maternal recombinant alternates between the grandfather’s X (autozygous, because the paternal copy *is* that X, intact) and the grandmother’s X. Switches at the crossover intensity λ(x) of one female meiosis. P(A=1)=½.
- **MS:** both X copies are independent single-meiosis recombinants of the **same** mother; A=1 where they coincide. A flips when *either* meiosis crosses over → intensity 2λ(x). P(A=1)=½.
- **SS:** the grandfather’s X enters only on the maternal side and can never be autozygous, forcing guaranteed non-autozygous blocks; A=1 requires the maternal transmission to be on the grandmother-derived side *and* the two grandmother homologs to coincide. P(A=1)=¼.

### 2.2 An exact likelihood from the map (no simulation)

Crossovers along a female meiosis form an **inhomogeneous Poisson process** whose intensity λ(x) *is* the recombination map. The likelihood of an observed junction configuration {x₁…x_m} under a process of intensity c·λ(x) is the standard
L = exp(−c∫λ) · Π_i c·λ(x_i).
Hence:

- **FD vs MS likelihood ratio.** Both intensities are the *same* λ(x), scaled by 1 (FD) vs 2 (MS):
  **L_FD / L_MS = exp(Λ) / 2^m**, Λ = ∫λ ≈ 1.76 Morgans.
  The λ(xᵢ) **cancel**: for FD vs MS the fine-grained map adds nothing beyond the junction *count* m. This is an exact statement of why FD vs MS is irreducibly hard from a single genome — it is the discrimination of Poisson(1.76) from Poisson(3.53) on one integer.
- **FD/MS vs SS.** SS has a different total intensity *and* a different structure (it is an 8-state hidden process over the latent homolog telegraphs, not a simple telegraph). The map does **not** cancel; the likelihood uses the full spatial signal. This is the pair the X resolves well.

The X-only Bayes-optimal classifier is therefore an explicit likelihood, not a black box (implemented in `xroh_likelihood.py`).

### 2.3 The autosomes: loop depth and meiosis sex

An autozygous tract whose two copies coalesce through g meioses has mean genetic length ≈ 100/g cM. For first-degree unions:

- **FD** (common ancestor = the grandfather G): paternal allele 1 meiosis from G, maternal allele 2 → **g = 3**.
- **MS** (common ancestor = the mother M): maternal allele 1 meiosis, paternal 2 → **g = 3**.
- **SS** (common ancestors = both grandparents): each side 2 meioses → **g = 4**.

So sib–sib offspring have **shorter, more numerous** autosomal ROH than parent–child offspring, at the same total F=¼. Within the g=3 parent–child class, FD and MS differ in the **sexes** of the loop meioses: FD = 2 male + 1 female; MS = 2 female + 1 male. Because female autosomal maps are ~1.57× longer than male, the MS loop accumulates ~17% more crossovers → MS ROH are slightly shorter and more numerous than FD. This is a faint autosomal trace of the *same* sex-path asymmetry the X shows strongly.

### 2.4 The demographic baseline, and why it barely matters at first degree

A real X:autosome ROH comparison sits on top of a demographic floor. Cotter et al. (2024) show that **even with no recent consanguinity, the X carries ~2× the autosomal ROH/IBD level**, because the X has a smaller effective population size (N_e,X ≈ ¾ N_e,A, further shifted by sex-biased demography) and hence shorter coalescence times. This baseline is the central confound for *population* inference from cousin marriage, where the loop contribution is small and must be separated from N_e,X effects.

At **first degree it is essentially negligible**: the loop adds F ≈ ¼ on the autosomes and F_X up to ½ — one to two orders of magnitude above the background autozygosity that the N_e,X factor modulates — so the union-type signal is read against a floor it dwarfs. Two consequences follow. First, the demographic identifiability problem that complicates Paper-B-style population scans does not materially degrade single-genome union typing. Second, the per-genome F_X/F_auto ratios we predict (2 for parent–child, **1 for sib–sib**) straddle Cotter’s ~2 demographic baseline: a sib–sib offspring’s recent-loop ratio of ~1 lies *below* the no-consanguinity population expectation, a reminder that the sex-path — not just depth — sets the ratio.

---

## 3. Methods

**Gene-dropping.** Founder haplotypes were assigned unique labels and dropped through each pedigree with sex-specific recombination: female meioses used the empirical female map, male meioses the male map (X: males transmit intact; sons receive no paternal X). Crossovers were drawn as an (inhomogeneous) Poisson process in genetic distance (Haldane; an optional gamma-renewal interference model is available). Autozygosity was scored as the genome fraction where the two copies carried the same founder label. PAR1/PAR2 excluded. Maps: Bhérer 2017 refined female X (GRCh37, 176.3 cM) and sex-specific autosomes; cross-checked on the deCODE/Halldorsson 2019 **maternal** chrX (GRCh38), built from the UCSC `recombMat` track and integrated to 175.88 cM (agreeing with Bhérer to 0.4 cM).

**X likelihood (`xroh_likelihood.py`).** The X is binned (0.5 Mb); per bin the single-meiosis switch probability is p_i = ½(1−e^{−2μ_i}), μ_i the local genetic length in Morgans. FD and MS likelihoods are closed-form telegraph products (flip probability p_i and 2p_i(1−p_i) respectively); SS is an 8-state HMM forward pass over the latent homolog telegraphs (V_B, V_Z, W_Z), all exact given the binning. Classification uses the maximum-posterior union at equal priors.

**Autosomal ROH.** All 22 autosomes were gene-dropped through each pedigree as **diploid** chromosomes with sex-specific maps (male meioses through males, female through females); ROH segments ≥1.5 Mb were extracted per genome and counted/measured (`autosomal_roh_by_union()`). Whole-genome accuracy is the **direct joint gene-drop simulation**: one gradient-boosted classifier (5-fold CV, n=600/union) on the combined per-genome feature vector [# autosomal ROH, autosomal F, mean and max autosomal ROH length, F_X, X-ROH segment count, longest X-ROH] (`whole_genome_all_pairs()`), reported with its cross-fold standard error. This replaces the earlier d′-in-quadrature approximation (which it reproduces within ~0.01 for all three pairs).

**Tail probabilities.** The whole-chromosome events (F_X=1, F_X=0) are computed in **closed form** (`xroh_tails.py`): FD/MS from the no-crossover probability of one/two meioses (½ e^{−Λ}, ½ e^{−2Λ}), SS fully-homozygous as ¼ e^{−3Λ}, and SS fully-outbred as the taboo survival probability of a three-telegraph Markov chain (the two autozygous states absorbing, via matrix exponential). These are exact and map-shape-independent. The gene-drop simulation validates them: the analytic value lies inside the **Wilson 95% CI** of the simulated proportion (n=300,000) for all six quantities. Simulated proportions are reported with Wilson intervals throughout (rather than the normal approximation, which fails in the small-Np tail).

**Reproducibility.** `python xroh_sim.py` reproduces the catalogue, detection, crossover-count decomposition (with Wilson CIs), genotyped-relative experiment, the autosomal and whole-genome union typing, and the interference sweep; `python xroh_tails.py` reproduces the closed-form tails and the analytic-vs-simulation validation table; `python xroh_likelihood.py` reproduces the exact-likelihood classification on both maps; `xroh_realdata.py` builds a real 1000G chrX haplotype panel (streamed) and validates ROH calling on it (§4.6); `xroh_posterior.py` computes calibrated posteriors (§6.2). `python reproduce_all.py` regenerates every headline number with fixed seeds, and `python test_xroh.py` (also runs under pytest) asserts the invariants — means → ½/½/¼, the {0,1,2,3} catalogue ratios, the closed-form tails inside their Wilson 95% CIs, and the genotyped-mother lift.

---

## 4. Results

### 4.1 Means match theory (validation)
Simulated F_X = 0.50 / 0.50 / 0.25 (FD/MS/SS) and F_auto = 0.25 for all three (outbred control 0). The full first/second/third-degree catalogue reproduces the analytic X/autosome ratios {0,1,2,3}. This convergence is the unit test that the sex-specific X transmission is coded correctly.

### 4.2 The X separates the sex-path
The exact map-based likelihood classifier gives, on a single female X (0.5 Mb bins, n=2000/union):

| pair | Bhérer X | deCODE X |
|---|---|---|
| FD vs SS | 0.736 | 0.741 |
| MS vs SS | 0.715 | 0.713 |
| FD vs MS | 0.710 | 0.712 |
| 3-way | 0.573 | 0.570 |

Two independent maps agree to ≈0.005 — the conclusion is map-robust. FD vs MS reaches exactly the analytic count-only ceiling (0.71), confirming §2.2: the spatial map cannot help that pair. The full-information tail is itself diagnostic: at the operational threshold F_X > 0.99, **P = 9.3% (FD) vs 2.0% (MS) vs 0.2% (SS)**; the *exact* fully-homozygous probabilities are 8.6/1.5/0.13% (closed forms, §4.2.1). Either way a near-fully-homozygous X is a strong FD signature (and explains the Sund 2013 “entirely homozygous X” girl).

#### 4.2.1 Decomposition by crossover count: the all-zero corner

The whole-chromosome outcomes are governed entirely by the **number of crossovers in the female-X meioses of the loop**, of which each union has a different count: FD is set by **one** meiosis (the maternal transmission, crossover count k ~ Poisson(λ), λ = Λ = 1.763 Morgans), MS by **two** independent meioses of the same mother (k₁+k₂), and SS by **three** (the operative one being k_s, the grandfather/grandmother mixing in the sister→child transmission). Because both extremes are *all-or-nothing along the whole chromosome*, their probabilities are **exact closed forms that depend only on the total genetic length Λ, not on the map’s spatial shape** (the map shape matters only for intermediate F_X and for the classifier accuracies). They are therefore reported as analytic values, not Monte-Carlo estimates, and the gene-drop simulation is used only to *validate* them — every analytic value below lies inside the simulation’s Wilson 95% CI at n = 300,000 (`xroh_tails.py`):

| union | governing meioses | derivation | **P(F_X=1)** (exact) | **P(F_X=0)** (exact) |
|---|---|---|---|---|
| FD | 1 (k) | ½ e^{−Λ} | **0.0858** | 0.0858 |
| MS | 2 (k₁+k₂) | ½ e^{−2Λ} | 0.0147 | 0.0147 |
| SS | 3 (k_s, k_g, k_b) | ¼ e^{−3Λ}  /  8-state survival | **0.00126** | **0.1785** |

The FD and MS extremes are the no-crossover events of one and two meioses; the SS fully-homozygous form ¼ e^{−3Λ} requires all three meioses crossover-free *and* a strand coincidence, while the SS fully-outbred probability is the taboo (never-autozygous) survival of a three-telegraph Markov chain over Λ (the two autozygous states made absorbing). Three facts follow. (i) A 0-crossover maternal transmission under **FD** is a coin flip between an **entirely homozygous** X (grandfather strand) and an **entirely outbred-looking** X (great-grandmother strand): P(F_X=1)=P(F_X=0)=½ e^{−Λ}=8.6%. (ii) **MS** needs the *same* zero event in **two** meioses (e^{−2Λ}), making its homozygous tail ~5.8× rarer than FD’s — the quantitative core of the FD-vs-MS hard pair seen from the tail. (iii) **SS is asymmetric**: the grandfather’s X is a one-sided structural intrusion that can never be autozygous, so a fully homozygous X is essentially impossible (0.13%) while a fully *outbred* X is common (18%). The two extremes therefore **invert the diagnosis**: a fully homozygous X gives a likelihood ratio FD:MS:SS = **68:12:1** (toward father–daughter), whereas a fully outbred X gives SS:FD:MS = **12:6:1** (toward sibling). The mean F_X is invariant to crossover count (½/½/¼ within every count class) — the union-type information the mean discards is entirely in the variance and these tails.

A note on definitions (and why these differ from a naïve threshold). The values above are the **exact** events F_X = 1 and F_X = 0. An operational caller instead reports, e.g., F_X > 0.99; that thresholded probability is *higher* (FD ≈ 0.096 vs the exact 0.086) because a crossover landing within ~1% of a telomere leaves a >99% but not fully homozygous X. The exact value is the right quantity for a likelihood; the thresholded value is what a given ROH pipeline at a stated resolution would actually call, and should be computed under that same definition when comparing to real data.

The tabled values are Haldane, which is an **upper bound** on the fully-homozygous tail: real meiosis has crossover interference and an obligate chiasma, both of which suppress the zero-crossover class. Under the gamma-renewal model with the human shape parameter (ν ≈ 4; Housworth & Stahl 2003), the realised P(0 crossovers on the X) falls from 0.171 to ≈0.025, so the exact P(F_X=1) falls to FD 1.3%, MS 0.03%, SS ≈0 (`tails_from_p0()`). Crucially, interference *sharpens* the FD-vs-MS contrast for this observation: because MS requires the zero-crossover event in **two** meioses (½p₀²) and FD in **one** (½p₀), the FD:MS likelihood ratio for a fully-homozygous X is 1/p₀, rising from ≈6:1 (Haldane) to ≈40:1 (ν≈4). A fully-homozygous X is thus an even stronger father–daughter signal under realistic recombination. In the strict obligate-chiasma limit (p₀→0) the exact event vanishes — so a clinically reported "entirely homozygous X" is necessarily the *operational* event (a single near-telomeric crossover leaving the X homozygous at array resolution), not a literal zero-crossover transmission.

### 4.3 The autosomes give the loop depth
Autosomal ROH (≥1.5 Mb; n=250/union; sex-specific maps):

| union | g | F_auto | # ROH | mean length |
|---|---|---|---|---|
| FD | 3 | 0.252 | 26.3 | 26.8 Mb |
| MS | 3 | 0.253 | 30.0 | 23.6 Mb |
| SS | 4 | 0.253 | 34.8 | 20.3 Mb |

All three share F=¼, yet the ROH **count orders them FD < MS < SS**. The MS/FD count ratio is 1.14 (predicted 1.16 from the female:male autosomal map lengths 4087:2606 cM) — the sex-specific-recombination effect of §2.3, confirmed. ROH count alone gives FD vs SS = 0.82, MS vs SS = 0.69, FD vs MS = 0.65.

### 4.4 Whole genome: the decisive contrast is the best resolved
A **direct joint classifier** on the whole-genome feature vector [# autosomal ROH, autosomal F, mean/max autosomal ROH length, F_X, X-ROH segment count, longest X-ROH] — using the full X signal, not F_X alone — gives the single-genome accuracies (5-fold CV, n=600/union, ± cross-fold SE):

| pair | combined accuracy |
|---|---|
| FD vs SS | **0.87 ± 0.01** |
| MS vs SS | **0.77 ± 0.01** |
| FD vs MS | **0.74 ± 0.01** |
| 3-way | 0.67 ± 0.004 |

These supersede the earlier d′-in-quadrature estimates and agree with them within ~0.01 — the quadrature approximation is vindicated by the joint simulation across all three pairs, not just FD-vs-SS. The decomposition for FD-vs-SS is X-only 0.72, autosomes-only 0.85, combined 0.88, showing the autosomes carry the depth and the X adds the sex-path. Sib–sib is well separated from both parent–child types; **father–daughter versus mother–son is the irreducible hard core** (joint 0.74, barely above the X-only junction-count ceiling of 0.71, since the autosomes add only the weak sex-path #ROH trace). In casework a genotyped birth mother breaks even that pair, lifting FD vs MS to **0.91**, because the child’s paternal X is an *intact* maternal homolog under FD but a *recombinant mosaic* of the mother’s two homologs under MS (an observable switch count of ~0 vs ~Poisson). Crossover interference (gamma-renewal, ν=2.6–4.3) leaves the means at ½/½/¼ but tightens segment spacing, with three consequences: it modestly raises the FD/MS ceiling (0.69→0.76), shrinks the FD homozygous-X tail (P(F_X>0.99) 0.093→0.015), and — because sharper F_X distributions are easier to separate — *raises* the X-only FD-vs-SS accuracy (0.64→0.70). So the Haldane (no-interference) defaults used throughout are the **conservative** end: realistic interference improves resolution of every contrast except the fully-homozygous-X tail probability, which it lowers.

### 4.5 Consistency with reported cases

No published case names the *type* of a first-degree union (clinical reports deliberately stop at "first-degree"), so real cases cannot give a ground-truth classification accuracy. They can, however, test the model for **quantitative consistency** and show the posterior it would have produced.

The sharpest is **Sund et al. (2013)**, whose motivating example is a female patient with first-degree parents in whom **the entire X chromosome was homozygous**, raising total homozygosity from **25.2% to 29.2%**. Both numbers match the model. (i) *Magnitude:* adding a fully-autozygous X (the covered non-PAR X ≈ 152 Mb) to a first-degree autosomal F = 0.252 on a ≈2.85 Gb genome gives 25.2% → **29.2%**, reproducing the report exactly. (ii) *Interpretation:* the report could only confirm "first-degree", but the closed-form tails convert a fully-homozygous X into a union-type posterior. With equal priors, P(F_X=1 | FD/MS/SS) = 0.086/0.015/0.0013 yields **P(union | F_X=1) = 0.84 (FD), 0.15 (MS), 0.01 (SS)** — the observation nearly excludes a sibling union and points strongly to parent–child, consistent with Sund's note that the homozygous X "supported the classification." This is precisely the discrimination the autosome-only pipeline discarded. Under realistic interference the call gets *sharper*, not weaker: with ν≈4 the FD:MS likelihood ratio for a fully-homozygous X rises to ≈40:1 (§4.2.1), pushing the posterior to ≈0.98 FD; and because an obligate chiasma makes a literal zero-crossover X vanishingly rare, the Sund "entirely homozygous X" is best read as the resolution-limited operational event (a near-telomeric crossover), still overwhelmingly indicating parent–child over sibling. The same paper recommends *excluding* the X from homozygosity calculations — i.e. discarding the signal we exploit.

Two further reports corroborate the framing rather than the numbers. **Schaaf et al. (2011)** frame the forensic question around whether the mother was a minor or adult and a perpetrator who is "a father or brother" — supporting the maternal-identity/age prior of §6.1, and tellingly considering only father–daughter and sib–sib while **omitting mother–son**. **Chaves et al. (2024)**, detecting first-degree incest by autosomal homozygosity in a large cohort, **exclude the sex chromosomes** from the calculation — the standard practice that sets aside exactly the chromosome that carries the type signal.

### 4.6 Real-haplotype check (1000 Genomes)

To test that the result survives real human haplotype structure (linkage disequilibrium, the realistic allele-frequency spectrum), we dropped **real phased 1000 Genomes chrX haplotypes** through the FD/MS/SS pedigrees — so the union type is ground truth — genotyped the child at real SNP positions, and called ROH with the same PLINK-style caller (`xroh_realdata.py`; 400 phased haplotypes, 4,791 SNPs over a 22 Mb Xp region, n=150/union). The caller recovers true autozygosity essentially **unbiased**: detected F_roh vs true (region) = 0.579/0.575 (FD), 0.534/0.528 (MS), 0.272/0.249 (SS), and the **½/½/¼ ordering is preserved under real LD** (the small SS upward bias is the expected LD-driven short-ROH inflation). Classification over this single 22 Mb segment is underpowered — it carries only ~1/6 of the X's genetic length, so the regional F_X is high-variance (FD-vs-SS 0.63) — and reproducing the full-X accuracies of §4.4 requires a whole-chromosome panel (the pipeline streams it on demand). The point established here is that the **detection layer (true autozygosity → called ROH) holds on real genomes**, so the simulation's idealised tracks are a faithful proxy; what remains is to scale the real-haplotype panel to the whole X (and to autosomes) for an end-to-end real-data accuracy.

---

## 5. Relationship to prior work

The closest prior art is Cotter, Severson, Kang, Godrej, Carmi & Rosenberg (2024, *G3* 14(2):jkad264). In a **diploid coalescent model**, they relate ROH and IBD-sharing at a site to coalescence time (TMRCA), which they obtain analytically as a function of consanguinity rates, and compare the autosomal and X-chromosomal models across a population in which a fraction of matings are first-cousin unions of four types — patrilateral-parallel, patrilateral-cross, matrilateral-parallel, matrilateral-cross (their Fig. 1). They predict that (i) X-IBD rises with X-ROH as consanguinity increases; (ii) the X carries ~2× the autosomal ROH/IBD **even without consanguinity**, from its smaller N_e and shorter coalescence time; and (iii) **matrilateral** consanguinity raises the X:autosome ratio above this baseline while **patrilateral** consanguinity lowers it. In genome-wide SNP data from human populations with estimated consanguinity, they find each 1% rise in autosomal ROH is associated with a **2.1% rise in X-chromosomal ROH** (1.6% for IBD), close to the ~2× their model predicts.

We use these results in three ways. We adopt their **demographic baseline** (prediction ii) as the floor against which first-degree signals are read (§2.4) and show it is negligible at first degree. We treat their **matrilateral-vs-patrilateral population result** (prediction iii) as the cousin-level, population-mean instance of the same sex-path mechanism whose first/second-degree, single-genome extension is our catalogue (the discrete X/autosome ratios {0,1,2,3}; their “patrilateral → no X ROH” is the cousin analogue of our X=0 configurations, e.g. fathers-are-brothers). And we position our work as **complementary and non-overlapping**: it characterises the **realised single-genome distribution** (and its tail) rather than population means; treats **first-degree** unions, where the sex-path signal is strongest and where they do not go; provides an **exact map-based likelihood** with the FD/MS cancellation result, where their coalescent treatment models only expectations (no variance) and no ROH-*calling*; and adds the **autosomal loop-depth/meiosis-sex** axis and its joint use with the X. The X-kinship coefficients underlying §2.1 trace to Grossman & Eisen (1989) and the McPeek-group derivations (Thornton 2012). Clinically, the autosomal-degree literature (Fan 2013; Sund 2013; Grote 2012; Bennett 2021) and ACMG guidance explicitly *avoid* naming the relationship and, in Sund’s case, *recommend excluding the X* — the chromosome that carries the sex-path signal.

The autosomal “is this ROH real autozygosity?” calibration (length → recombination- and ancestry-aware Bayes factor) is a **separate** problem handled in a companion project; here the loop-depth result is used only to *interpret* called ROH, and that calibration should be cited, not re-derived.

### 5.1 What is new here, and what is inherited

We are explicit about the boundary, because two of the building blocks are not ours.

**Inherited (cite, do not claim):**
- *The expected values F_X = ½/½/¼.* These follow from established X-coancestry coefficients (Grossman & Eisen 1989; Thornton 2012) — they are an evaluation of known theory, not a new derivation.
- *The autosomal depth law (mean ROH ≈ 100/g cM; sib–sib g=4 vs parent–child g=3).* This is standard ROH/IBD-relationship inference (Thompson; Ceballos 2018; Ringbauer 2021) and underpins existing relationship estimators.
- *Population-mean X:autosome ROH under consanguinity type.* Owned by Cotter et al. 2024 for first cousins.

**New contributions (specific):**
1. **Single-genome union *typing* of first-degree unions from the X sex-path.** Prior X work types *degree* or, in forensics, detects incest while omitting mother–son; the population work gives means only. Casting FD/MS/SS as a single-genome classification with an explicit posterior is, to our search, unaddressed.
2. **An exact, closed-form union likelihood read directly from the recombination map** (inhomogeneous-Poisson/HMM; `xroh_likelihood.py`), replacing simulation or black-box classifiers — and reaching the analytic Bayes ceiling.
3. **The FD-vs-MS cancellation identity: L_FD/L_MS = e^Λ/2^m.** This *proves* the spatial map is uninformative for FD vs MS — converting the empirical observation “FD vs MS is hard” into a theorem and showing exactly where the single-genome information limit sits. We are not aware of this result in the literature.
4. **The autosomal sex-path trace.** A daughter’s *autosomal* ROH count weakly but measurably encodes whether the loop ran father-side (FD: 2 male + 1 female meiosis) or mother-side (MS: 2 female + 1 male), through the 1.57× female:male map-length ratio (predicted 1.16, observed 1.14 MS/FD count ratio). Distinguishing two unions of the *same* depth g via meiotic-sex composition is, to our knowledge, novel — and it is the autosomal echo of the X sex-path.
5. **The joint autosome–X decomposition** (“autosomes give depth, the X gives sex-path”) with whole-genome single-genome accuracies, identifying FD-vs-SS (the forensically decisive contrast) as the *best*-resolved pair (≈0.87) and FD-vs-MS as the irreducible core.
6. **The realised single-genome distribution and its tail as discriminators** — e.g. P(F_X>0.99) = 9.3/2.0/0.2% — versus the population means of prior work, with the mechanistic link to the Sund 2013 fully-homozygous-X case.
7. **The genotyped-relative power result** (FD vs MS 0.69→0.91) via the intact-paternal-X switch count, defining the minimal extra data that breaks the hard pair.

In short, the *framing* (single-genome typing) and three *technical* results — the cancellation identity (3), the autosomal sex-path trace (4), and the joint decomposition with the genotyped-relative fix (5,7) — are the genuinely new content; the coefficient values and the depth law are inherited and cited.

---

## 6. Discussion

The information needed to type a first-degree union is present in a single female genome, distributed across two axes that the conventional autosomal *burden* statistic averages away. The autosomes encode the **depth and meiotic sex** of the consanguinity loop in the ROH length spectrum; the X encodes the **sex-path** in its autozygosity level and junction structure. Read jointly, they make the forensically and clinically consequential contrast — **father–daughter (an adult perpetrator and a child) versus sibling union** — the best-resolved question (≈0.87 from one genome), while father–daughter versus mother–son remains near the information limit unless a relative is genotyped.

Two results are, to our knowledge, new: the **exact likelihood with the FD/MS spatial cancellation** (which converts “FD vs MS is hard” from an empirical observation into a theorem), and the **autosomal sex-path trace** (a daughter’s autosomal ROH count weakly encodes whether the loop ran father-side or mother-side, via sex-specific recombination).

### 6.1 The mother is almost always known — which is exactly where the hard pair gives way

The single-genome accuracies above are deliberately conservative: they use *only* the child’s genome. But there is a structural asymmetry in real casework that bears directly on the one pair the genome cannot resolve. **The mother’s identity is essentially always known** — she was pregnant and gave birth, so she is identified and, in most settings, available — whereas the **father is frequently unknown or absent**: he is the perpetrator under father–daughter and the son under mother–son. This matters because FD-vs-MS is the genetically irreducible pair (the X likelihood *cancels*, §2.2), yet it is precisely the pair that a **genotyped birth mother breaks** (FD vs MS → 0.91, §4.4), via the intact-vs-recombinant paternal-X signature. So the relative who would be most useful is the one routinely *present*, while the relative who is typically missing (the father) is the one the method does not need. The genotyped-mother result is therefore not a theoretical add-on but the **realistic operating mode** of the method: type the child, obtain the mother, resolve the hard pair.

The birth mother’s **age** points the same way as a weaker, secondary prior — under MS she necessarily already has a reproductive-age son and is thus older, whereas under FD/SS she is of the offspring generation — but we do not lean on it; it is asymmetric (a young mother weighs against MS; an older one is ambiguous) and easily confounded, and the substantive point is her *availability for genotyping*, not her age.

One caution: any such case-context fact must enter as an **explicit prior and not be double-counted** — if the mother’s circumstances are what raised suspicion, they cannot also serve as independent corroboration — and the genomic likelihood, conditionally independent of how the case came to attention, remains the evidentiary core.

### 6.2 Calibrated likelihood ratios for casework

A forensic report needs a *calibrated* likelihood ratio: when the method states P(FD)=0.9 it should be right about 90% of the time. The exact-likelihood posterior (§2.2) is well-calibrated **when the recombination model matches the data** — simulating and scoring under Haldane gives an expected calibration error (ECE) of **0.013** (`xroh_posterior.py`). It degrades under **model mismatch**: data generated with realistic crossover interference (ν=4.3) but scored with the Haldane likelihood is miscalibrated (ECE **0.16**, in the conservative/underconfident direction — a stated confidence of 0.54 corresponds to an empirical accuracy of 0.72). The practical implication is concrete and was not obvious a priori: the headline *accuracies* are robust to (indeed improved by) interference, but the *calibrated probabilities* are not — **court-grade LRs require the likelihood itself to use a realistic, interference-aware recombination map**, anchored to a measured ν (Housworth & Stahl 2003; Campbell et al. 2015) or, better, the empirical female-X crossover-count distribution. Folding in the §6.1 maternal prior is then immediate: a young birth mother sets the MS prior to ≈0, collapsing the three-way problem (accuracy 0.57) to the well-separated **FD-vs-SS** pair (0.74) — the genome answers precisely the question that remains once the non-genomic facts have removed mother–son.

---

## 7. Limitations

- Accuracies are for a **single** genome and are now from a direct joint gene-drop classifier for all three pairs (with 5-fold-CV standard errors), superseding the earlier d′-in-quadrature approximation, which they reproduce within ~0.01.
- The genotype layer is validated against real 1000G haplotypes for the *detection* step (§4.6: ROH calling recovers true autozygosity unbiased over a 22 Mb chrX segment); the remaining step is to scale the real-haplotype panel to the whole X and the autosomes for an end-to-end real-data classification accuracy (the `xroh_realdata.py` pipeline streams the whole chromosome on demand).
- Crossover interference is modelled as gamma-renewal; the human interference parameter (ν ≈ 4–11; Housworth & Stahl 2003) and its known increase in deregulation with maternal age (Campbell et al. 2015) should anchor the ν used for final tail/ceiling claims. An obligate chiasma on the X further suppresses the zero-crossover class below the Haldane value, lowering the exact fully-homozygous-X probabilities (these are reported as Haldane upper bounds).
- Elevating resolution from “first-degree” to a named union raises consent, reporting and psychosocial stakes; any application must engage the ELSI literature (Tarini 2013; Grote 2012/2014; Helm 2013; Bennett 2021).

---

## 8. Conclusion

All first-degree unions look alike on the autosomal inbreeding coefficient, but not on the genome. The autosomal ROH spectrum reveals the loop’s depth and meiotic sex; the X chromosome of a female offspring reveals its sex-path through an exactly computable likelihood. Together they type the union from one genome — separating sibling from parent–child with ~0.87 accuracy and bounding how far father–daughter versus mother–son can ever be pushed without a second sample. The decisive signal is partly carried by the one chromosome routine pipelines set aside.

---

## References (verify formatting before submission)

1. Bhérer C, Campbell CL, Auton A. Refined genetic maps reveal sexual dimorphism in human meiotic recombination at multiple scales. *Nat Commun* 2017;8:14994. doi:10.1038/ncomms14994
2. Halldorsson BV, et al. Characterizing mutagenic effects of recombination through a sequence-level genetic map. *Science* 2019;363:eaau1043. doi:10.1126/science.aau1043 (deCODE maternal chrX via UCSC `recombMat` track).
3. Cotter DJ, Severson AL, Kang JTL, Godrej HN, Carmi S, Rosenberg NA. Modeling the effects of consanguinity on autosomal and X-chromosomal runs of homozygosity and identity-by-descent sharing. *G3* 2024;14(2):jkad264. doi:10.1093/g3journal/jkad264
4. Grossman M, Eisen EJ. Inbreeding, coancestry, and covariance between relatives for X-chromosomal loci. *J Hered* 1989;80(2):137–142. doi:10.1093/oxfordjournals.jhered.a110812 (PMID 2926116).
5. Thornton T, Zhang Q, Cai X, Ober C, McPeek MS. XM: association testing on the X-chromosome in case-control samples with related individuals. *Genet Epidemiol* 2012;36(5):438–450. doi:10.1002/gepi.21638
6. Ceballos FC, Joshi PK, Clark DW, Ramsay M, Wilson JF. Runs of homozygosity: windows into population history and trait architecture. *Nat Rev Genet* 2018;19:220–234. doi:10.1038/nrg.2017.109
7. Fan YS, et al. Frequent detection of parental consanguinity in children with developmental disorders by a combined CGH and SNP microarray. *Mol Cytogenet* 2013;6:38. doi:10.1186/1755-8166-6-38
8. Sund KL, et al. Regions of homozygosity identified by SNP microarray analysis aid in the diagnosis of autosomal recessive disease and incidentally detect parental blood relationships. *Genet Med* 2013;15:70–78. doi:10.1038/gim.2012.94
9. Grote L, et al. Variability in laboratory reporting practices for regions of homozygosity. *Genet Med* 2012;14:971–976. doi:10.1038/gim.2012.83
10. Tarini BA, et al. The perils of SNP microarray testing: uncovering unexpected consanguinity. *Pediatr Neurol* 2013;49:50–53. doi:10.1016/j.pediatrneurol.2013.03.008
11. Helm BM, et al. Three clinical experiences with SNP array results consistent with parental incest. *J Genet Couns* 2013;23:489–495. doi:10.1007/s10897-013-9669-0
12. Bennett RL, et al. Genetic counseling and screening of consanguineous couples and their offspring: focused revision. *J Genet Couns* 2021;30:1354–1357. doi:10.1002/jgc4.1477
13. Ringbauer H, Novembre J, Steinrücken M. Parental relatedness through time revealed by runs of homozygosity in ancient DNA. *Nat Commun* 2021;12:5425. doi:10.1038/s41467-021-25289-w
14. Schaaf CP, Scott DA, Wiszniewska J, Beaudet AL. Identification of incestuous parental relationships by SNP-based DNA microarrays. *Lancet* 2011;377(9765):555–556. doi:10.1016/S0140-6736(11)60201-8 (PMID 21315943).
15. Palsson G, Hardarson MT, Jonsson H, et al. Complete human recombination maps. *Nature* 2025;639(8055):700–707. doi:10.1038/s41586-024-08450-5 (PMID 39843742).
16. Housworth EA, Stahl FW. Crossover interference in humans. *Am J Hum Genet* 2003;73(1):188–197. doi:10.1086/376610 (PMID 12772089).
17. Campbell CL, Furlotte NA, Eriksson N, Hinds D, Auton A. Escape from crossover interference increases with maternal age. *Nat Commun* 2015;6:6260. doi:10.1038/ncomms7260 (PMID 25695863).
18. Hammer MF, Mendez FL, Cox MP, Woerner AE, Wall JD. Sex-biased evolutionary forces shape genomic patterns of human diversity. *PLoS Genet* 2008;4(9):e1000202. doi:10.1371/journal.pgen.1000202 (PMID 18818765).
19. Wilson Sayres MA. Genetic diversity on the sex chromosomes. *Genome Biol Evol* 2018;10(4):1064–1078. doi:10.1093/gbe/evy039 (PMID 29635328).
