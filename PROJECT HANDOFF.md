# PROJECT HANDOFF — X-chromosome autozygosity & the type of consanguineous union

Handoff for continuing this project in **VS Code with the Claude Code extension**. This captures the full exchange: the scientific idea, the code, the results, the literature, the two draft papers, data provenance, and the open threads.

-----

## 1. What this project is

A female child of an incestuous/consanguineous union carries information on the **X chromosome** that the autosomes do not, because of X transmission asymmetry: a father passes his single X **intact** to daughters (no recombination), a son gets **no** X from his father, and a mother passes a recombined X. As a result, **X autozygosity (F_X) in a female offspring depends on the sex-path of the consanguinity loop, not just its depth** — while the autosomal inbreeding coefficient F is the same for many distinct union types.

Two papers fall out of this:

- **Paper A (clinical / forensic):** can the X-ROH pattern tell *which* first-degree union (father–daughter / mother–son / brother–sister) produced a girl? Draft: `X_ROH_incest_type_methods_note.md`.
- **Paper B (population / anthropological):** the population-mean F_roh(X)/F_roh(autosome) ratio fingerprints the **marriage system** (e.g., patrilateral parallel-cousin marriage → near-zero X ROH). Draft: `X_autozygosity_mating_type_fingerprint.md`.

-----

## 2. Headline results (all simulation-confirmed)

**First-degree, female offspring:** F_X = **½ (father–daughter), ½ (mother–son), ¼ (brother–sister)**; F_auto = ¼ for all three. So the X/autosome ratio is ~2 for parent–child vs ~1 for sib–sib → **sib–sib is cleanly separable; father–daughter vs mother–son is the hard pair.**

**Detection robustness:** with a realistic SNP panel (array ~12k / WGS ~50k) + 0.5% genotyping error and a PLINK-style ROH caller, **detected F_X ≈ true F_X (bias ~0.001)** — first-degree ROH are chromosome-scale, so the limit is statistical, not technical.

**Bayes ceiling, father–daughter vs mother–son:** junctions ~ Poisson(1.76) vs Poisson(3.53); optimal accuracy **0.71 (junction count)**, **0.69 (full-info classifier)** — near a coin flip; segment lengths add ~nothing. From one genome this pair is close to intractable.

**Catalogue (X/autosome ratio by union type, female offspring):** within each autosomal-F class the ratio takes values **{0, 1, 2, 3}** by sex-path alone. Zero whenever the loop reaches the father only via *his* father (paternal half-sib; fathers-are-brothers cousins; aunt-via-paternal-grandfather; PGF–MGM cousins). Maximal (3×) for all-female matrilineal loops (mothers-are-sisters cousins; aunt = father’s-mother’s sister, F_X ≈ 0.375). Full table is in Paper B and `X_signature_catalogue.png`.

-----

## 3. Files in this folder

|File                                       |What it is                                                                                                                       |
|-------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------|
|`xroh_sim.py`                              |Consolidated, runnable engine: gene-drop, pedigree catalogue, SNP+ROH detection, Bayes ceiling. **Reproduces all numbers above.**|
|`requirements.txt`                         |Python deps.                                                                                                                     |
|`CLAUDE.md`                                |Orientation for the Claude Code agent.                                                                                           |
|`X_ROH_incest_type_methods_note.md`        |Paper A draft (clinical/forensic).                                                                                               |
|`X_autozygosity_mating_type_fingerprint.md`|Paper B draft (population/marriage-system).                                                                                      |
|`incest_X_summary.png`                     |Paper A 3-panel figure (ratio / detection / Bayes).                                                                              |
|`Fx_incest_sim.png`                        |F_X distributions by first-degree type.                                                                                          |
|`chrX_female_recomb.png`                   |Female X recombination landscape (PAR/centromere).                                                                               |
|`X_signature_catalogue.png`                |Paper B catalogue figure (13 configs).                                                                                           |

Run: `pip install -r requirements.txt && python xroh_sim.py`

-----

## 4. Data provenance & network notes (important)

- **Female X recombination map — Bhérer et al. 2017 (GRCh37).** Auto-downloaded by `xroh_sim.py` from
  `https://raw.githubusercontent.com/cbherer/Bherer_etal_SexualDimorphismRecombination/master/Refined_genetic_map_b37.tar.gz`.
  Contains `female_chrN.txt`, `male_chrN.txt`, `sexavg_chrN.txt`. **chrX has only `female_chrX.txt`** (no male/sexavg X — the non-PAR X recombines only in females). Columns: `chr, pos, rate(cM/Mb), cM(cumulative)`. Female X length ≈ **176.3 cM**.
- **deCODE 2019 map — Halldorsson et al. (GRCh38), cross-check.** Official files (`aau1043_datas3.gz`, sex-averaged; for chrX the average = the female rate) live on science.org / UCSC. We pulled a working copy from a GitHub mirror; it is **gzip-compressed despite the `.gz`-as-text claim in some repos** — decompress it. Columns: `Chr, Begin, End, cMperMb, cM`. deCODE chrX length 175.9 cM ≈ Bhérer 176.3 cM (validated).
- **1000 Genomes / HGDP VCFs were NOT reachable** from the original build sandbox (allowed domains were limited to GitHub/PyPI/etc.). The SNP genotype layer is therefore **simulated** at realistic density. **In VS Code you have open network access** — swap in real founder haplotypes from 1000G/HGDP for the detection study and for Paper B’s biobank scan.
- **ROH caller:** we use a transparent **PLINK-style sliding-window** caller (`call_roh`). We tried `scikit-allel`’s `roh_mhmm` (OOMs at chromosome scale — allocates per-bp arrays) and `roh_poissonhmm` (collapsed to one segment on this synthetic data); both are avoidable. For real data, PLINK `--homozyg` or `bcftools roh` are the standards.

-----

## 5. Literature status

**The key gap:** as far as PubMed + web searching reached, **nobody has used ROH on the X to discriminate the *type* of union.** The X-kinship theory exists (McPeek framework; X-coancestry derivations), clinical reports use autosomal ROH for the *degree* (and per ACMG avoid naming the relationship), and forensic X-STR/ML treats the X as a complementary *detector* (omitting mother–son). The closest is Sund et al. 2013, which only used X homozygosity to support first-degree-vs-not and recommended *excluding* the X.

**Real-world anchor cases (girls born of incest, SNP array):**

- Fan et al. 2013, Case 001 — 6-y girl, brother–sister parents, ~598 Mb ROH (F≈¼); tentative α-ketoglutarate dehydrogenase deficiency. doi:10.1186/1755-8166-6-38
- Fan et al. 2013, Case 004 — 23-day female neonate, uncle–niece parents (2nd degree), 345 Mb ROH. (same paper)
- Sund et al. 2013 — a first-degree girl with an **entirely homozygous X** (+~4% homozygosity). doi:10.1038/gim.2012.94

**Other clinical/ethics refs:** Helm 2013 (doi:10.1007/s10897-013-9669-0); Tarini 2013 (infant, sex not stated, half-sib; doi:10.1016/j.pediatrneurol.2013.03.008); Wang 2014 (doi:10.1038/ejhg.2014.153); Grote 2012 (doi:10.1038/gim.2012.83), 2013 (doi:10.1002/ajmg.a.36206); Bennett 2021 (doi:10.1002/jgc4.1477); Schaaf et al. Lancet 2011 (a boy; verify citation).

-----

## 6. Open threads / TODO (roughly prioritized)

1. **Add crossover interference** (replace Haldane/Poisson with the map’s Campbell-style interference) and re-check segment counts and the Bayes ceiling.
1. **Add a genotyped parent/relative** to the simulation — likely breaks the father–daughter/mother–son degeneracy via the intact-paternal-X check. The obvious power experiment for Paper A.
1. **Real-data swap:** pull 1000G/HGDP haplotypes; redo the detection study on real allele frequencies and LD.
1. **Paper B biobank scan:** compute F_roh(X)/F_roh(autosome) across populations with contrasting marriage systems; design the **sex-biased-demography control** (Ne_X/Ne_A also moves the ratio — must disentangle).
1. **SLiM forward model:** evolve a population under a fixed cousin-marriage rule; confirm the genealogical ratio prediction survives drift/finite-size.
1. **`ped-sim` cross-check:** reproduce the catalogue with the Williams-lab `ped-sim` (uses the same Bhérer map natively) as an independent engine.
1. **Reference cleanup (both papers):** verify the items below.

-----

## 7. References to verify before submission

- McPeek X-kinship / `KinInbcoefX` — find the correct primary citation (software vs paper).
- PubMed **PMID 2926116** (X-chromosomal inbreeding/coancestry derivations) — confirm authors/journal.
- Schaaf et al., *Lancet* 2011 — confirm authors/volume/pages.
- X-to-autosome diversity-ratio / sex-biased-demography refs (Hammer 2008 *Nat Genet*; Wilson Sayres 2018 *GBE*) — confirm exact citations.

-----

## 8. Working norms (carry into Claude Code)

- **Don’t fabricate** results or citations; verify present-day facts and references rather than asserting from memory.
- When citing PubMed-derived articles, attribute and include DOI links.
- Prefer **simulation-confirmed** numbers; the gene-drop means should converge to the analytic F_X (½/½/¼, the {0,1,2,3} catalogue ratios) — that convergence is the unit test that the X/sex transmission is coded correctly.