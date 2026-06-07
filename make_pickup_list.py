"""Generate docs/papers_to_obtain.pdf — a one-page pickup list of the
paywalled references still needed for the ROH_on_X project, to fetch via
institutional access. Clickable DOI links. Everything else is already in docs/.
"""
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

GENERATED = "2026-06-07"

PAPERS = [
    dict(
        tag="X-kinship theory anchor (HIGH)",
        authors="Grossman M, Eisen EJ.",
        title="Inbreeding, coancestry, and covariance between relatives for X-chromosomal loci.",
        venue="J Hered. 1989;80(2):137–142.",
        doi="10.1093/oxfordjournals.jhered.a110812",
        pmid="2926116",
        why="Primary derivation behind the 1/2, 1/2, 1/4 F_X result. Note: animal-breeding journal -- cite as such.",
    ),
    dict(
        tag="KinInbcoefX primary citation (HIGH)",
        authors="Thornton T, Zhang Q, Cai X, Ober C, McPeek MS.",
        title="XM: Association testing on the X-chromosome in case-control samples with related individuals.",
        venue="Genet Epidemiol. 2012;36(5):438–450.",
        doi="10.1002/gepi.21638",
        pmid="22552845",
        why="The McPeek-group paper defining the X-kinship/inbreeding coefficients that KinInbcoefX computes.",
    ),
    dict(
        tag="Clinical case set (completeness)",
        authors="Tarini BA, et al.",
        title="The perils of SNP microarray testing: uncovering unexpected consanguinity.",
        venue="Pediatr Neurol. 2013;49(1):50–53.",
        doi="10.1016/j.pediatrneurol.2013.03.008",
        pmid="23827427",
        why="Clinical/ethics case in the incest-by-array literature cited in Paper A.",
    ),
    dict(
        tag="Clinical case set (completeness)",
        authors="Grote L, et al.",
        title="Variable approaches to genetic counseling for microarray regions of homozygosity associated with parental relatedness.",
        venue="Am J Med Genet A. 2014;164A(1):87–98.",
        doi="10.1002/ajmg.a.36206",
        pmid="",
        why="The 2014 Grote (you already have the 2012 Grote). Completes the lab-reporting-practice pair.",
    ),
]

fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
fig.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.05)
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")

y = 0.965
ax.text(0.07, y, "ROH_on_X — Papers still to obtain", fontsize=19, fontweight="bold")
y -= 0.028
ax.text(0.07, y, f"Fetch via institutional access at work.  Generated {GENERATED}.",
        fontsize=10, color="#444")
y -= 0.016
ax.text(0.07, y, "Drop the PDFs into  ROH_on_X/docs/  and they'll be renamed/verified. "
                 "All other refs are already in docs/ (20 papers).",
        fontsize=9, color="#666")
y -= 0.022
ax.plot([0.07, 0.95], [y, y], color="#999", lw=0.8)
y -= 0.030

for i, p in enumerate(PAPERS, 1):
    ax.text(0.07, y, f"{i}.", fontsize=12, fontweight="bold")
    ax.text(0.11, y, p["tag"], fontsize=9.5, fontweight="bold", color="#1a5fb4")
    y -= 0.024
    ax.text(0.11, y, p["authors"], fontsize=10.5)
    y -= 0.022
    # title may wrap; matplotlib won't auto-wrap, so wrap manually
    import textwrap
    for line in textwrap.wrap(p["title"], width=92):
        ax.text(0.11, y, line, fontsize=10.5, style="italic")
        y -= 0.020
    ax.text(0.11, y, p["venue"], fontsize=10.5)
    y -= 0.022
    doi_url = f"https://doi.org/{p['doi']}"
    ax.text(0.11, y, f"DOI: {p['doi']}", fontsize=10, color="#0b6e0b", url=doi_url)
    if p["pmid"]:
        ax.text(0.62, y, f"PMID: {p['pmid']}", fontsize=10, color="#0b6e0b",
                url=f"https://pubmed.ncbi.nlm.nih.gov/{p['pmid']}/")
    y -= 0.020
    ax.text(0.11, y, f"Get: {doi_url}", fontsize=9.5, color="#1a5fb4", url=doi_url)
    y -= 0.022
    for line in textwrap.wrap("Why: " + p["why"], width=95):
        ax.text(0.11, y, line, fontsize=9, color="#444")
        y -= 0.018
    y -= 0.012
    ax.plot([0.07, 0.95], [y, y], color="#ddd", lw=0.6)
    y -= 0.026

_footer = ("Already obtained & in docs/ (20): bherer-2017, halldorsson-2019, ceballos-2018, "
           "ringbauer-2021, fan-2013, sund-2013, schaaf-2011, helm-2013, grote-2012, wang-2014, "
           "bennett-2021, delgado-2015, frontiers-2025-ml, bittles-black-2010, cox-hammer-2008, "
           "wilson-sayres-2018, chaves-2024, cotter-rosenberg-2024, hammer-2008, genes-2022-x-str, "
           "coancestry-inbreeding-lecture.")
_fy = 0.075
for _line in textwrap.wrap(_footer, width=115):
    ax.text(0.07, _fy, _line, fontsize=7.2, color="#888")
    _fy -= 0.013

with PdfPages("docs/papers_to_obtain.pdf") as pdf:
    pdf.savefig(fig)
plt.close(fig)
print("wrote docs/papers_to_obtain.pdf")
