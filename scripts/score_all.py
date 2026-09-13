from collections import defaultdict, Counter
import statistics as st, glob

cls_of = {}
for pat, c in (("*big_members.tsv","BIG"), ("*medium_members.tsv","MEDIUM"),
               ("*small_members.tsv","SMALL")):
    for f in glob.glob(pat):
        for line in open(f):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 6:
                cls_of[p[5]] = c
                cls_of.setdefault(p[4], c)
print("members labelled from TSVs:", len(cls_of))

sc = defaultdict(dict)
for line in open("all_4116.tsv"):
    q, s, pid, ln, b = line.rstrip("\n").split("\t")
    k = s.split("_")[0]
    sc[q][k] = sc[q].get(k, 0.0) + float(b)
print("queries with BLAST hits:", len(sc))

conf = Counter()
margins = defaultdict(list)
unlabelled = 0
for q, d in sc.items():
    truth = cls_of.get(q)
    if not truth:
        unlabelled += 1
        continue
    ranked = sorted(d.items(), key=lambda kv: -kv[1])
    conf[(truth, ranked[0][0])] += 1
    if len(ranked) > 1:
        margins[truth].append(ranked[0][1] - ranked[1][1])
if unlabelled:
    print("queries with no class label:", unlabelled)

order = ("BIG", "MEDIUM", "SMALL")
print("\nCONFUSION MATRIX  (rows = truth, cols = predicted)")
print("%-8s %8s %8s %8s %10s" % ("truth", *order, "n"))
for t in order:
    n = sum(conf[(t, p)] for p in order)
    if n:
        print("%-8s %8d %8d %8d %10d" % (t, *[conf[(t, p)] for p in order], n))

print("\nACCURACY AND MARGIN")
print("%-8s %10s %10s %12s" % ("class", "n", "correct", "med margin"))
for t in order:
    n = sum(conf[(t, p)] for p in order)
    if not n:
        continue
    m = st.median(margins[t]) if margins[t] else 0
    print("%-8s %10d %9.2f%% %12.0f" % (t, n, 100.0 * conf[(t, t)] / n, m))
