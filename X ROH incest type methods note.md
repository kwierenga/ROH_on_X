# Runs of homozygosity on the X chromosome partially resolve the *type* of first-degree incest in a female offspring

*Working methods-note draft — v0.1. Numbers are from the simulations described below; references marked “verify” need a final citation check.*

## Summary

Autosomal runs of homozygosity (ROH) detected on SNP microarrays reliably reveal that a child’s parents are first-degree relatives, but cannot say *which* first-degree union (father–daughter, mother–son, or brother–sister) produced the child: all three give an autosomal inbreeding coefficient F ≈ ¼. Because of the X chromosome’s sex-specific transmission, the X carries information the autosomes do not. Using gene-dropping simulation on empirical sex-specific genetic maps, we show that in a **female** offspring the expected X autozygosity is **½ for father–daughter and mother–son but only ¼ for brother–sister**, so the X-to-autosome homozygosity ratio cleanly separates sib–sib (≈1) from parent–child (≈2). Father–daughter versus mother–son, however, share the same expected X burden and differ only in ROH segment architecture; the Bayes-optimal accuracy for that pair from a single genome is ~0.68–0.71. Detection is robust: because first-degree ROH are chromosome-scale, the X signal survives sparse markers and genotyping error. We note that this is an underused signal — clinical laboratories frequently *exclude* the X from homozygosity calculations — and that suitable clinical data (SNP arrays on girls born of incest) already exist.

## Background

- SNP/chromosomal microarray is first-tier testing for developmental delay, intellectual disability, autism, and congenital anomalies, and incidentally reveals parental relatedness through genome-wide ROH.
- Incidental detection of first-degree incest on these arrays is documented, alongside the attendant consent, reporting, and psychosocial issues (Grote 2012, 2013; Tarini 2013; Helm 2013; Bennett 2021; Wang 2014; Fan 2013; Schaaf 2011).
- Standard practice estimates the *degree* of relatedness from autosomal ROH burden (~25% → first-degree) and, per ACMG guidance, deliberately does **not** assign the specific parental relationship, since array data is not a paternity test and ROH size maps imperfectly to relationship (Fan 2013).
- The X chromosome is transmitted asymmetrically: a father passes his single X intact to every daughter (no recombination); a mother passes a recombined X; males carry only a maternal X. This asymmetry makes X autozygosity depend on the *sexes in the consanguinity loop*, not just its depth.

## Theory

The X inbreeding coefficient of a female equals the X-kinship coefficient between her parents (McPeek framework). Evaluating the three first-degree pedigrees:

|Union (female offspring)|Autosomal F|Expected X F (F_X)|
|------------------------|-----------|------------------|
|Father × daughter       |¼          |**½**             |
|Mother × son            |¼          |**½**             |
|Brother × sister        |¼          |**¼**             |

Consequences:

- The autosomes are blind to type (all ¼).
- The **F_X / F_auto ratio** is ≈2 for parent–child and ≈1 for sib–sib → separates sib–sib from the parent–child pair.
- Father–daughter and mother–son share F_X ≈ ½ and differ only in segment structure: father–daughter has one recombination-free X copy (the intact paternal X), so its autozygous tracts are bounded by a *single* maternal meiosis; mother–son’s two copies are separated by *two* meioses.

## Methods

**Gene-dropping simulation.** Founder X (and autosomal) haplotypes were assigned unique labels and dropped through each pedigree with sex-specific recombination: female meioses used the empirical female map (Bhérer 2017 refined map, GRCh37; cross-checked against deCODE/Halldorsson 2019, GRCh38), males transmitted the X intact, and autosomes recombined in both sexes using sex-specific maps. Crossovers were drawn as a Poisson process in genetic distance (Haldane, no interference). Autozygosity was scored as the genome fraction where both copies carried the same founder label. PAR1/PAR2 were excluded (consistent with the maps). Female X genetic length ≈ 176 cM (λ ≈ 1.76 crossovers per female meiosis); the empirical female-X map and the deCODE chrX map agreed on X length to within ~0.4 cM.

**Realistic genotype layer.** SNP panels at array-like (~12k) and WGS-like (~50k) density across the X were simulated with a common-variant allele-frequency spectrum; founder haplotype alleles were dropped through the pedigree to produce offspring genotypes, with 0.5% genotyping error injected as false heterozygotes. ROH were called with a standard PLINK-style sliding-window algorithm (window 50 SNPs, ≤2 heterozygotes, ≥1.5 Mb, ≥50 SNPs).

**Bayes ceiling (father–daughter vs mother–son).** The number of IBD/non-IBD junctions on the X follows Poisson(λ) for father–daughter (one maternal meiosis) and Poisson(2λ) for mother–son (two female meioses). Optimal accuracy from the junction count is ½·Σ_k max[Poisson(k;λ), Poisson(k;2λ)]. A simulation classifier (gradient boosting, 5-fold CV) trained on full features (segment count, F_X, longest segment, length spread) gives the full-information ceiling.

## Results

- **Means match theory** (validation): F_X = 0.50 / 0.50 / 0.25 for father–daughter / mother–son / sib–sib, and F_auto ≈ 0.25 for all three (outbred control 0). *(Figure 2)*
- **Ratio separates sib–sib:** F_X/F_auto ≈ 2.0 (father–daughter), 2.1 (mother–son), 1.0 (sib–sib). *(Figure 3, panel 1)*
- **Detection is robust:** detected F_X tracked true F_X to within ~0.001 at both densities despite genotyping error — first-degree ROH are chromosome-scale, so the limit is statistical, not technical. *(Figure 3, panel 2)*
- **Father–daughter vs mother–son is the wall:** equal F_X; mean detected ROH segment counts ≈1.3 vs ≈2.1. Bayes accuracy from junction count = **0.71**; full-information classifier = **0.68** (chance = 0.50); segment lengths add essentially nothing beyond the count. *(Figure 3, panel 3)*

Net: the X turns “which first-degree union?” into one easy classification (sib–sib vs parent–child) and one near-intractable one (father–daughter vs mother–son) from a single female genome.

## Real-world anchor cases

Three published clinical cases anchor the simulation in girls actually born of incest:

- **Fan et al. 2013, Case 001** — a 6-year-old girl, ~598 Mb of homozygosity (F ≈ ¼), referred for intellectual disability; parents a **brother–sister** pair; later given a tentative diagnosis of α-ketoglutarate dehydrogenase deficiency uncovered through the homozygosity. A sib–sib girl — predicted lower X burden (F_X ≈ ¼); her X was not specifically flagged.
- **Fan et al. 2013, Case 004** — a 23-day-old female neonate with congenital anomalies, 345 Mb homozygosity (F ≈ 1/8); parents an **uncle–niece** pair (second-degree).
- **Sund et al. 2013** — among 25 female patients in whom the X was evaluated, one **girl with a suspected first-degree relationship had an entirely homozygous X**, raising total homozygosity by ~4% and supporting the first-degree call. A fully homozygous X (F_X ≈ 1) is the high tail that essentially only parent–child unions reach — the opposite end of the distribution from Fan’s sib–sib girl.

These two girls bracket the model: a sib–sib girl with low X burden and a parent–child girl with a fully homozygous X. Notably, Sund et al. recommended that laboratories consider **excluding** the X from homozygosity calculations — i.e., the one chromosome carrying the type-discriminating signal is the one routinely discarded.

## The gap

No published work appears to use X-ROH to discriminate the *type* of first-degree union. The X-kinship theory exists but has not been applied this way; clinical reports use autosomal ROH for the *degree* and (per ACMG) avoid assigning the relationship; forensic X-STR and machine-learning work treat the X as a complementary *detector* of incest and omit mother–son. The proposal here — read the X-ROH pattern in a female offspring to separate sib–sib from parent–child, and to bound how far father–daughter vs mother–son can ever be pushed — is, to our knowledge, unaddressed.

## Limitations

- The SNP genotypes are simulated at realistic density rather than drawn from a specific reference panel (1000G VCFs were not accessible in the compute environment); the detection conclusions concern density and het-noise, which this reproduces.
- Crossovers are modeled without interference (Haldane); adding a chiasma-interference model (e.g., the map’s native interference) would slightly tighten segment counts.
- PAR1/PAR2 are excluded.
- The accuracy bounds are for a **single** female genome. In casework, a genotyped parent or relative would likely break the father–daughter/mother–son degeneracy via the intact-paternal-X check — the obvious next experiment.
- Heightening the resolution from “first-degree” to a specific union raises the reporting and consent stakes; any application must engage the ELSI literature (Tarini 2013; Grote 2012/2013; Bennett 2021).

## Conclusion

The X chromosome adds genuine, quantifiable information about the type of first-degree incest in a female offspring: it cleanly separates brother–sister from parent–child via a ~2-fold difference in X autozygosity that is robust to realistic genotyping, while father–daughter versus mother–son remains near a coin flip from a single genome (~0.68–0.71 ceiling). The signal is currently discarded in clinical practice, and the data to test it — SNP arrays on girls born of incest — already exist.

## Figures

- **Figure 1.** Female recombination landscape of chrX (Bhérer refined map), with PAR1 excluded, the pericentromeric cold valley, and the distal rise toward PAR2. *(chrX_female_recomb.png)*
- **Figure 2.** Distribution of X autozygosity F_X by pedigree (means 0.50 / 0.50 / 0.25 / 0), and the father–daughter vs mother–son segment-count contrast. *(Fx_incest_sim.png)*
- **Figure 3.** Three-panel summary: (1) F_X/F_auto ratio splits sib–sib from parent–child; (2) detected vs true F_X (detection ≈ truth); (3) Poisson junction-count likelihoods and the father–daughter vs mother–son Bayes ceiling. *(incest_X_summary.png)*

## References (verify formatting before submission)

1. Bhérer C, Campbell CL, Auton A. Refined genetic maps reveal sexual dimorphism in human meiotic recombination at multiple scales. *Nat Commun* 2017;8:14994. doi:10.1038/ncomms14994
1. Halldorsson BV, et al. Characterizing mutagenic effects of recombination through a sequence-level genetic map. *Science* 2019;363:eaau1043. doi:10.1126/science.aau1043
1. Zhang Q, Bourgain C, McPeek MS. KinInbcoefX: calculation of X-chromosome kinship and inbreeding coefficients (software/definitions). [verify primary citation]
1. X-chromosomal inbreeding, coancestry, and covariance between relatives — explicit derivations for full sibs and other relatives. PubMed PMID 2926116. [verify author/journal]
1. Fan YS, et al. Frequent detection of parental consanguinity in children with developmental disorders by a combined CGH and SNP microarray. *Mol Cytogenet* 2013;6:38. doi:10.1186/1755-8166-6-38
1. Sund KL, et al. Regions of homozygosity identified by SNP microarray analysis aid in the diagnosis of autosomal recessive disease and incidentally detect parental blood relationships. *Genet Med* 2013;15:70–78. doi:10.1038/gim.2012.94
1. Tarini BA, et al. The perils of SNP microarray testing: uncovering unexpected consanguinity. *Pediatr Neurol* 2013;49:50–53. doi:10.1016/j.pediatrneurol.2013.03.008
1. Helm BM, et al. Three clinical experiences with SNP array results consistent with parental incest: a narrative with lessons learned. *J Genet Couns* 2013;23:489–495. doi:10.1007/s10897-013-9669-0
1. Wang JC, et al. Regions of homozygosity identified by oligonucleotide SNP arrays: evaluating the incidence and clinical utility. *Eur J Hum Genet* 2014;23:663–671. doi:10.1038/ejhg.2014.153
1. Grote L, et al. Variability in laboratory reporting practices for regions of homozygosity indicating parental relatedness as identified by SNP microarray testing. *Genet Med* 2012;14:971–976. doi:10.1038/gim.2012.83
1. Grote L, et al. Variable approaches to genetic counseling for microarray regions of homozygosity associated with parental relatedness. *Am J Med Genet A* 2014;164A:87–98. doi:10.1002/ajmg.a.36206
1. Bennett RL, et al. Genetic counseling and screening of consanguineous couples and their offspring: focused revision. *J Genet Couns* 2021;30:1354–1357. doi:10.1002/jgc4.1477
1. Schaaf CP, et al. Identification of incestuous parental relationships by SNP-based DNA microarrays. *Lancet* 2011. [verify authors/volume/pages]
1. Ringbauer H, Novembre J, Steinrücken M. Parental relatedness through time revealed by runs of homozygosity in ancient DNA. *Nat Commun* 2021;12:5425. doi:10.1038/s41467-021-25289-w
1. Miles A, et al. scikit-allel (software; roh caller).