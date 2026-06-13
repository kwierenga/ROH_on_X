"""
fetch_decode_xmap.py — build the deCODE maternal chrX genetic map (GRCh38).
============================================================================
Independent cross-check map for xroh_likelihood.py. The non-PAR X recombines
only in females, so the deCODE *maternal* recombination rate is the relevant X
map. UCSC hosts it as the `recombMat` bigWig track; we pull it over the UCSC
REST API in chunks, integrate the per-interval rate (cM/Mb) to cumulative cM,
and write data/decode_chrX_maternal_hg38.txt (columns: pos<TAB>cM).

Validated total length ~175.9 cM, matching Bhérer 2017 female X (176.3 cM).

Source: Halldorsson et al. 2019 Science (deCODE), via UCSC hg38 recombRate track.
Run:    python fetch_decode_xmap.py
"""
import json, os, time, urllib.request

CHROM, SIZE, STEP = "chrX", 156_040_895, 20_000_000
BASE = ("https://api.genome.ucsc.edu/getData/track"
        "?genome=hg38;track=recombMat;chrom=chrX")
OUT = os.path.join("data", "decode_chrX_maternal_hg38.txt")


def main():
    ivals = []
    for s in range(0, SIZE, STEP):
        e = min(s + STEP, SIZE)
        url = f"{BASE};start={s};end={e}"
        for attempt in range(3):
            try:
                d = json.load(urllib.request.urlopen(url, timeout=120))
                break
            except Exception as ex:
                print("retry", s, ex); time.sleep(2)
        chunk = d["recombMat"]
        ivals.extend(chunk)
        print(f"{s:>11}-{e:<11} intervals={len(chunk):>5} total={len(ivals)}")

    ivals = [iv for iv in ivals if iv["end"] > iv["start"]]
    ivals.sort(key=lambda x: x["start"])
    pos, cm, cum = [ivals[0]["start"]], [0.0], 0.0
    for iv in ivals:                                  # rate (cM/Mb) -> cumulative cM
        cum += iv["value"] * (iv["end"] - iv["start"]) / 1e6
        pos.append(iv["end"]); cm.append(cum)

    print(f"\nchrX maternal span {pos[0]/1e6:.1f}-{pos[-1]/1e6:.1f} Mb, "
          f"total {cum:.2f} cM (Bhérer female X = 176.3 cM)")
    os.makedirs("data", exist_ok=True)
    with open(OUT, "w") as f:
        f.write("pos\tcM\n")
        for p, c in zip(pos, cm):
            f.write(f"{int(p)}\t{c:.6f}\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
