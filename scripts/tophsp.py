import re
from collections import defaultdict

CRE = re.compile(r"^(BIG|SMALL|SINGLE)_.*?_size(\d+)_(.+)$")

m2r = {}
for p in ("full_big_cluster_members.tsv", "full_small_cluster_members.tsv"):
    try:
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 6:
                m2r[f[5]] = f[4]; m2r.setdefault(f[4], f[4])
    except FileNotFoundError:
        pass

rep2c = {}
for line in open("consensus_samples.fasta"):
    if line.startswith(">"):
        cid = line[1:].split()[0]
        m = CRE.match(cid)
        if m:
            rep2c[m.group(3)] = cid

summed = defaultdict(lambda: defaultdict(float))
best   = defaultdict(lambda: defaultdict(float))
for line in open("blast_work/forward.tsv"):
    q, s, pid, ln, b = line.rstrip("\n").split("\t")
    b = float(b)
    summed[q][s] += b
    best[q][s] = max(best[q][s], b)

def score(table, label):
    stats = defaultdict(lambda: [0, 0])
    for q, d in table.items():
        rep = m2r.get(q)
        truth = rep2c.get(rep) if rep else None
        if not truth:
            continue
        cls = CRE.match(truth).group(1)
        top = max(d.items(), key=lambda kv: kv[1])[0]
        stats[cls][0] += 1
        stats[cls][1] += (top == truth)
    print(f"\n{label}")
    for cls in ("BIG", "SMALL", "SINGLE"):
        if cls in stats:
            n, c = stats[cls]
            print(f"  {cls:<8} {c}/{n}  ({100*c/n:.1f}%)")

score(summed, "SUMMED bitscore across all HSPs (what I reported)")
score(best,   "BEST single HSP only")

print("\nshare of summed score coming from the top HSP:")
fr = []
for q in summed:
    for s in summed[q]:
        if summed[q][s] > 0:
            fr.append(best[q][s] / summed[q][s])
fr.sort()
print(f"  median {100*fr[len(fr)//2]:.1f}%   min {100*fr[0]:.1f}%   max {100*fr[-1]:.1f}%")
