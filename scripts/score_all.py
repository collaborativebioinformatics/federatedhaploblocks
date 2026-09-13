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

stats = defaultdict(lambda: {"n":0, "uc":0, "uw":0, "tie_inc":0, "tie_exc":0})
conf = Counter(); margins = defaultdict(list); unlabelled = 0

for q, d in sc.items():
    truth = cls_of.get(q)
    if not truth:
        unlabelled += 1; continue
    s = stats[truth]; s["n"] += 1
    top = max(d.values())
    winners = [k for k, v in d.items() if v == top]
    if len(winners) == 1:
        pred = winners[0]
        conf[(truth, pred)] += 1
        s["uc" if pred == truth else "uw"] += 1
        rest = [v for v in d.values() if v != top]
        if rest: margins[truth].append(top - max(rest))
    else:
        s["tie_inc" if truth in winners else "tie_exc"] += 1

if unlabelled: print("queries with no class label:", unlabelled)
order = ("BIG", "MEDIUM", "SMALL")

print("\nUNAMBIGUOUS ASSIGNMENTS ONLY (ties excluded)")
print("%-8s %7s %9s %9s %8s %8s %11s" % (
    "class","n","unique ok","unique no","tie+own","tie-own","med margin"))
for t in order:
    s = stats.get(t)
    if not s or not s["n"]: continue
    m = st.median(margins[t]) if margins[t] else 0
    print("%-8s %7d %9d %9d %8d %8d %11.0f" % (
        t, s["n"], s["uc"], s["uw"], s["tie_inc"], s["tie_exc"], m))

print("\nACCURACY BOUNDS")
print("%-8s %7s %14s %14s %12s" % ("class","n","lower","upper","decided only"))
for t in order:
    s = stats.get(t)
    if not s or not s["n"]: continue
    n = s["n"]; ties = s["tie_inc"] + s["tie_exc"]; dec = s["uc"] + s["uw"]
    print("%-8s %7d %13.2f%% %13.2f%% %11.2f%%" % (
        t, n, 100.0*s["uc"]/n, 100.0*(s["uc"]+ties)/n,
        100.0*s["uc"]/dec if dec else 0))
print("\nlower = ties counted wrong, upper = ties counted right,")
print("decided only = accuracy over unambiguous assignments")

print("\nCONFUSION MATRIX, unambiguous only (rows = truth)")
print("%-8s %8s %8s %8s" % ("truth", *order))
for t in order:
    if sum(conf[(t, p)] for p in order):
        print("%-8s %8d %8d %8d" % (t, *[conf[(t, p)] for p in order]))
