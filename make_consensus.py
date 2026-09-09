#!/usr/bin/env python3
import sys
from collections import Counter
import numpy as np

BASES = "ACGT"
LUT = np.full(256, 4, dtype=np.uint8)
for i, b in enumerate(BASES):
    LUT[ord(b)] = i
    LUT[ord(b.lower())] = i


def iter_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0] if len(line) > 1 else "", []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def main():
    if len(sys.argv) != 4:
        sys.exit("usage: make_consensus.py <cluster.fasta> <name> <out.fasta>")
    path, label, out_path = sys.argv[1:4]

    lengths = Counter(len(s) for _, s in iter_fasta(path))
    if not lengths:
        sys.exit("no sequences found")
    total = sum(lengths.values())
    length, count = lengths.most_common(1)[0]
    print("length distribution (top 6):")
    for L, c in lengths.most_common(6):
        print("   {:>9,} bp   {:>6} seqs   {:5.1f}%".format(L, c, 100.0 * c / total))
    print("\nusing modal length {:,} bp: {} of {} sequences ({:.1f}%)".format(
        length, count, total, 100.0 * count / total))
    if length == 0:
        sys.exit("modal length is zero, check the FASTA")

    counts = np.zeros((length, 5), dtype=np.int32)
    rows = np.arange(length)
    n = 0
    for sid, seq in iter_fasta(path):
        if len(seq) != length:
            continue
        counts[rows, LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]] += 1
        n += 1
        if n % 250 == 0:
            print("  {} tallied".format(n))

    acgt = counts[:, :4]
    depth = acgt.sum(axis=1)
    top = acgt.max(axis=1)
    consensus = np.array(list(BASES), dtype="<U1")[acgt.argmax(axis=1)]
    consensus[depth == 0] = "N"

    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(depth > 0, top / depth, 1.0)

    poly = int(((acgt > 0).sum(axis=1) > 1).sum())
    print("\n{} sequences used".format(n))
    print("polymorphic positions : {:,}  ({:.4f}%)".format(poly, 100.0 * poly / length))
    print("majority base under 90% : {:,}".format(int((frac < 0.9).sum())))
    print("majority base under 60% : {:,}".format(int((frac < 0.6).sum())))

    seq = "".join(consensus)
    with open(out_path, "w") as fh:
        fh.write(">{}_consensus_n{}\n".format(label, n))
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    print("\nwrote {}".format(out_path))


main()
