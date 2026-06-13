# How a single genome can tell *which kind* of first-degree incest produced a girl — step by step

*A plain-language companion to the methods note. Robust but readable: every step is something the simulation in this repository actually shows. “First-degree” means parents as closely related as possible short of being the same person: father–daughter (FD), mother–son (MS), or brother–sister (SS).*

---

## Step 0 — The puzzle

When a child is conceived by close relatives, large stretches of their two genome copies are identical, because both copies descend from the same recent ancestor. These stretches are **runs of homozygosity (ROH)**. Clinical genetic tests see them and can say *“the parents are first-degree relatives.”*

But here is the catch. The standard summary number — the fraction of the genome that is autozygous, called **F** — is **¼ for all three** first-degree unions. Father–daughter, mother–son, brother–sister: identical on that number. So the usual reading of the test **cannot tell which kind of union it was.** That distinction matters enormously: a father–daughter case means an adult abusing a child; a sibling case is a different situation legally, clinically, and ethically.

The question of this work: **is the missing information really gone, or is it hiding somewhere the standard summary throws away?** It’s hiding — in two places.

---

## Step 1 — Why the autosomes seem blind

The 22 non-sex chromosomes (autosomes) are inherited symmetrically: you get one copy from each parent, and both parents shuffle (recombine) their two copies before passing one on. For all three first-degree unions, the total amount of autozygous genome works out to the same ¼. That’s why the headline F number is uninformative. Hold that thought — the autosomes are *not* actually blind, but you have to look past the total (Step 6).

---

## Step 2 — The X chromosome breaks the symmetry

The X is inherited in a lopsided way:

- A **father gives his single X to his daughters completely intact** — no shuffling. (He has only one X; he hands it over whole.)
- A **son gets no X at all from his father** (he gets the Y instead).
- A **mother shuffles her two X’s** and passes one mixed copy to each child.

Because of this, in a **girl**, how much of her two X’s match depends on **which sexes lie along the path connecting her parents back to their shared ancestor** — not just how close that ancestor is. We call that path the **sex-path**.

---

## Step 3 — The three unions give three different X signatures

Run the inheritance forward for a daughter:

- **Father–daughter (FD):** her father hands her his X *intact*. That same X also sits, intact, inside her mother (because the mother is the father’s daughter). So one of the girl’s X’s is an unshuffled copy that the other X is trying to match. Result: on average **half** the X is autozygous (**F_X = ½**), in one or a few long blocks.
- **Mother–son (MS):** both of the girl’s X’s are shuffled copies of the *same* grandmother-mother’s two X’s — but shuffled in two separate events. They still match **half** the time on average (**F_X = ½**), but broken into *more* pieces, because two shufflings make more seams.
- **Brother–sister (SS):** the matching can only happen through the shared grandmother’s X, and the grandfather’s X (which rides along on one side) can never match. Result: only **a quarter** of the X is autozygous (**F_X = ¼**).

So just from the *level* of X matching: **sib–sib (¼) stands apart from the two parent–child types (½).** That’s the first real signal the standard F number hid.

---

## Step 4 — Reading the pattern as “seams”

Picture the girl’s X as a strip that is either “matching” or “not matching” as you walk along it. The points where it flips are **seams**, and each seam is a recombination (shuffling) event. The three stories predict different numbers of seams:

- **FD:** seams come from **one** shuffling → few seams.
- **MS:** seams come from **two** shufflings → about twice as many seams.
- **SS:** fewer matching regions overall, with built-in non-matching blocks.

The recombination **map** tells us exactly how likely a seam is at each spot on the X. So instead of guessing, we can **calculate the probability of the observed seam-pattern under each story** and compare. That probability is the *likelihood*, and the ratio of two likelihoods is the weight of evidence for one story over another. No guessing, no simulation — a formula read off the map.

---

## Step 5 — Why father–daughter vs mother–son is genuinely stuck

When you write down that calculation for FD vs MS, something clean happens: the **map cancels out**. The only thing that ends up mattering is the **number of seams** — FD averages about 1.8, MS about 3.5. From a single genome you see *one* number of seams, and “about 1.8” versus “about 3.5” overlap a lot. So FD vs MS is **near a coin flip (~0.71 at best)** — and we can prove it’s not a tooling limitation; it’s the information limit. The fine detail of *where* the seams fall, which helps elsewhere, gives nothing here.

For **FD vs SS** (and MS vs SS) the map does **not** cancel — the stories differ in the total amount of matching and in their structure — so the full pattern helps. This is the good news, because…

---

## Step 6 — The autosomes are not blind after all (two hidden signals)

Go back to the autosomes. The *total* is ¼ for everyone, but the **sizes and number of the ROH pieces differ**, for two reasons:

**(a) Depth.** An autozygous piece gets chopped up by every shuffling between the child and the shared ancestor. Parent–child loops pass through **3** shufflings; sibling loops pass through **4** (because siblings share *two* grandparents, one extra step on each side). More shufflings → **shorter, more numerous** pieces. So:

| union | typical # of ROH pieces | typical size |
|---|---|---|
| father–daughter | ~26 | larger |
| mother–son | ~30 | medium |
| brother–sister | ~35 | smaller |

All add up to ¼ — but the **count** cleanly separates sib–sib from parent–child. And because this is averaged over 22 chromosomes, it’s a **reliable** signal, unlike the single, noisy X.

**(b) Meiosis sex (the subtle one).** Female shuffling produces ~1.6× more seams than male shuffling. A father–daughter loop runs mostly through *male* steps; a mother–son loop runs mostly through *female* steps. So mother–son fragments a bit more than father–daughter even though both have the same depth — giving the autosomes a faint echo of the *same* sex-path story the X tells loudly. It’s weak, but real (we predicted a 16% effect and measured 14%).

---

## Step 7 — Putting both chromosomes together: the scoreboard

Combine the autosomes (depth) and the X (sex-path), and ask how well a single genome can tell each pair apart:

| question | answer from one genome |
|---|---|
| **father–daughter vs sibling** | **~0.87** — well resolved |
| mother–son vs sibling | ~0.77 |
| father–daughter vs mother–son | ~0.75 — the stubborn pair |

The headline: the **most consequential question — father–daughter versus sibling — is the best-answered one.** The autosomes nail the depth; the X confirms the sex-path. The only pair that stays hard is father–daughter vs mother–son, because they are identical in depth *and* in X level, differing only in the noisy seam-count.

---

## Step 8 — How to crack the last hard pair

Father–daughter vs mother–son is stuck **from the child alone**. But add **one genotyped relative — the birth mother — and it jumps to ~0.91.** The reason is elegant: under father–daughter, the girl’s intact paternal X is a perfect, unshuffled copy of one of her mother’s X’s; under mother–son, it’s a *shuffled mix* of both. Counting how many times the girl’s paternal X “switches” between her mother’s two X’s gives the answer: about zero (father–daughter) versus several (mother–son). The child’s genome can’t show this; the child plus mother can.

---

## Step 9 — Why this matters, and the irony

- **Clinically and forensically**, separating father–daughter from sibling is the question that changes what happens next — and it’s the one a single genome answers best.
- **The irony:** the X chromosome carries the sex-path signal, yet clinical labs routinely **exclude the X** from homozygosity calculations. One published case even describes a girl from a first-degree union whose X was *entirely* homozygous — exactly the father–daughter fingerprint — noted only in passing. The discriminating chromosome is the one usually set aside.

---

## Step 10 — How solid is this?

- The core inheritance numbers (½, ½, ¼) **match theory exactly** in simulation — that agreement is the built-in correctness check.
- Everything was checked on **two independent recombination maps** (Bhérer and deCODE), which agree to within half a percent — so the conclusions don’t depend on map choice.
- The honest limits: these are single-genome accuracies; the genotype data here are realistically simulated, not yet real patient data; and turning “first-degree” into a *named* relationship raises real consent and reporting questions that any use must respect.

**In one sentence:** the autosomes tell you *how deep* the family loop is, the X tells you *which sexes* it ran through, and together — from a single girl’s genome — they reveal not just *that* the parents were first-degree relatives, but *which kind*, with the most important distinction being the clearest.
