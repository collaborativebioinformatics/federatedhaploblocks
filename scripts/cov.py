from collections import defaultdict

QLEN = 152010  # chr1:10023878-10175888

pairs = defaultdict(lambda: {"n": 0, "total": 0, "longest": 0})
for line in open("blast_work/forward.tsv"):
    q, s, pid, ln, b = line.rstrip("\n").split("\t")
    if "chr1" not in q or "chr1" not in s:
        continue
    k = (q, s.split("_")[0])
    p = pairs[k]
    p["n"] += 1
    p["total"] += int(ln)
    p["longest"] = max(p["longest"], int(ln))

by_class = defaultdict(list)
for (q, cls), p in pairs.items():
    by_class[cls].append(p)

print("per query-subject pair, chr1 only")
print("%-8s %8s %14s %14s" % ("class", "med HSPs", "med total cov", "med longest HSP"))
for cls in ("BIG", "SMALL", "SINGLE"):
    v = by_class.get(cls)
    if not v:
        continue
    med = lambda key: sorted(x[key] for x in v)[len(v) // 2]
    print("%-8s %8d %13.1f%% %13.1f%%" % (
        cls, med("n"),
        100.0 * med("total") / QLEN,
        100.0 * med("longest") / QLEN))
