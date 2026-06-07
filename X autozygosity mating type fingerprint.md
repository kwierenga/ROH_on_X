# X-chromosomal autozygosity as a fingerprint of consanguineous mating type

*Working methods-note draft — v0.1. A population-genetics / anthropological-genetics companion to the clinical incest-typing note. Numbers are simulation results; references marked “verify” need a citation check.*

## Summary

The degree of relatedness between a person’s parents is read from genome-wide autosomal runs of homozygosity (ROH), but autosomal ROH are blind to the *sex-structure* of the consanguinity loop: father–daughter, mother–son, and brother–sister unions all give F = ¼; maternal and paternal half-sib unions both give ⅛; the four wirings of first-cousin parents all give 1/16. Because the X chromosome is transmitted asymmetrically (a father passes his single X intact only to daughters, and a son receives no X from his father), X autozygosity in a **female** offspring depends on the sexes along the loop, not only its depth. We catalogue this with gene-dropping simulation on an empirical female recombination map and show that, holding autosomal F fixed within each relationship class, the ratio F_X/F_auto takes the discrete values **0, 1, 2, or 3** according to the sex-path. The X is therefore a genome-readable fingerprint of *which kind* of consanguineous union — and, at the population level, of which marriage system a society practises.

## The transmission rule

For a female child:

- Her **paternal X** is her father’s single X, transmitted intact; it derives entirely from her **paternal grandmother**.
- Her **maternal X** is a recombinant of her mother’s two X’s, deriving from the maternal grandfather (intact) and maternal grandmother.

Consequently a consanguinity loop contributes to her X autozygosity **only if it threads ancestral X-material into the paternal grandmother through females** (and likewise reaches the maternal grandparents). Any father→son step on the path carries no X and zeroes that branch. This single rule generates the entire catalogue below.

## Catalogue (female offspring; simulation, female X map ≈ 176 cM)

|Class       |Autosomal F|Configuration                  |F_X  |F_X/F_auto|
|------------|-----------|-------------------------------|-----|----------|
|First degree|¼          |father–daughter                |0.50 |2         |
|First degree|¼          |mother–son                     |0.50 |2         |
|First degree|¼          |brother–sister                 |0.25 |1         |
|Half-sib    |⅛          |maternal half-sibs             |0.25 |2         |
|Half-sib    |⅛          |paternal half-sibs             |0.00 |0         |
|Avuncular   |⅛          |uncle = mother’s brother       |0.125|1         |
|Avuncular   |⅛          |uncle = father’s brother       |0.25 |2         |
|Avuncular   |⅛          |aunt = father’s-mother’s sister|0.375|3         |
|Avuncular   |⅛          |aunt = father’s-father’s sister|0.00 |0         |
|First cousin|1/16       |mothers are sisters            |0.188|3         |
|First cousin|1/16       |fathers are brothers           |0.00 |0         |
|First cousin|1/16       |PGM & MGF are sibs             |0.125|2         |
|First cousin|1/16       |PGF & MGM are sibs             |0.00 |0         |

(PGM/PGF/MGF/MGM = paternal/maternal grand-mother/father.) Analytic values agree with simulation (e.g., mothers-are-sisters cousins F_X = 3/16; aunt-via-paternal-grandmother F_X = 3/8).

**Key regularities.**

- Within each autosomal-F class, the autosomes are constant but F_X spans 0 → 3× — the X resolves variation the autosomes cannot see.
- F_X = 0 whenever the loop reaches the father only through *his* father (paternal half-sib; fathers-are-brothers cousins; aunt via paternal grandfather; PGF–MGM cousins). These are genome configurations with autosomal homozygosity but **no** X homozygosity.
- F_X is maximised (3×) by **all-female / matrilineal** loops (mothers-are-sisters cousins; aunt = father’s-mother’s sister).
- Father–daughter and mother–son both sit at 2× — the within-class degeneracy explored in the companion clinical note.

## Why this is a population-genetic signal

Human societies have systematic, non-random marriage preferences, and the preferred type fixes the sex-path of the consanguinity loop:

- **Patrilateral parallel-cousin marriage** (father’s-brother’s-daughter), common in parts of the Middle East/North Africa, is the *fathers-are-brothers* configuration → autosomal ROH present but **X ROH ≈ 0**.
- **Matrilateral cross-cousin** and matrilineal preferences route the loop through females → **elevated** X ROH.

Therefore the population-mean **F_roh(X) / F_roh(autosome) ratio** is a genome-readable index of the descent/marriage system, distinguishable from a simple consanguinity-rate measure (which only sees autosomal F). This is the non-clinical, evolutionary-anthropology payoff: marriage systems leave a chromosome-specific autozygosity signature recoverable from genotype data alone, without genealogies.

## Methods

Gene-dropping simulation: founder X haplotypes were dropped through each pedigree with female-only X recombination (empirical female map; males transmit the X intact), and F_X scored as the fraction of the X where a female offspring’s two X’s carry the same founder haplotype (PAR excluded; Haldane crossovers). Autosomal F is the textbook value for each relationship. 6,000 replicates per configuration.

## Proposed empirical tests

1. **Cross-population biobank scan.** Compute F_roh(X)/F_roh(autosome) per individual in HGDP, 1000 Genomes, and biobanks; compare population means against ethnographically documented marriage systems. Prediction: populations practising patrilateral parallel-cousin marriage show the *lowest* X/autosome ROH ratio at matched autosomal consanguinity.
1. **Forward simulation under a fixed marriage rule.** Use a sex-explicit forward simulator (e.g., SLiM) to evolve a population under a prescribed cousin-marriage preference and read the equilibrium X/autosome ROH ratio, calibrating the genealogical prediction against drift and finite-size effects.
1. **Within-population validation.** Where pedigrees exist (founder populations, genealogical biobanks), test whether the X/autosome ROH ratio recovers the loop type case-by-case.

## Caveats

- **Demographic confounding.** The X/autosome ratio also responds to sex-biased effective population size (baseline Ne_X/Ne_A = ¾, shifted by sex-biased reproductive variance and migration). Disentangling the mating-type signal from sex-biased demography is essential — and is itself an interesting inference problem (the diversity-ratio literature provides priors).
- **Population-level statistic.** Per-individual F_X is high-variance (short, one-sex-recombining chromosome), so this is a distributional/population signal, not a single-genome call (except at first-degree).
- **ROH detection on the X** requires sex-aware calling; male hemizygosity and reduced Ne miscalibrate standard diploid callers.
- Crossovers modelled without interference; PAR excluded.

## Relationship to the clinical note

The companion note treats the *individual* question (which first-degree union produced this child, in casework). This note treats the *population* question (which marriage system a society practises), using the same transmission rule and simulation engine but a different unit of analysis and a non-clinical motivation.

## Figure

- **Figure 1.** X autozygosity F_X by mating type in a female offspring (bars), with autosomal F marked (ticks); colour = relationship class. Within each class the autosomal F is constant while F_X varies 0–3×. *(X_signature_catalogue.png)*

## References (verify formatting before submission)

1. Bhérer C, Campbell CL, Auton A. Refined genetic maps reveal sexual dimorphism in human meiotic recombination at multiple scales. *Nat Commun* 2017;8:14994. doi:10.1038/ncomms14994
1. Halldorsson BV, et al. Characterizing mutagenic effects of recombination through a sequence-level genetic map. *Science* 2019;363:eaau1043. doi:10.1126/science.aau1043
1. Zhang Q, Bourgain C, McPeek MS. X-chromosome kinship and inbreeding coefficients (KinInbcoefX). [verify primary citation]
1. X-chromosomal inbreeding and coancestry derivations for relatives. PubMed PMID 2926116. [verify author/journal]
1. Ceballos FC, Joshi PK, Clark DW, Ramsay M, Wilson JF. Runs of homozygosity: windows into population history and trait architecture. *Nat Rev Genet* 2018;19:220–234. doi:10.1038/nrg.2017.109
1. Bittles AH, Black ML. Consanguinity, human evolution, and complex diseases. *Proc Natl Acad Sci USA* 2010;107(Suppl 1):1779–1786. doi:10.1073/pnas.0906079106
1. Sex-biased demography and the X-to-autosome diversity ratio (e.g., Hammer et al. 2008 *Nat Genet*; Wilson Sayres 2018 *Genome Biol Evol*). [verify exact citations]
1. Ringbauer H, Novembre J, Steinrücken M. Parental relatedness through time revealed by runs of homozygosity in ancient DNA. *Nat Commun* 2021;12:5425. doi:10.1038/s41467-021-25289-w