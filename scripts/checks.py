from collections import defaultdict
sc = defaultdict(dict)
for line in open("all_4116.tsv"):
    q, s, pid, ln, b = line.rstrip("\n").split("\t")
    k = s.split("_")[0]
    sc[q][k] = sc[q].get(k, 0.0) + float(b)
ties = near = clear = 0
for q, d in sc.items():
    v = sorted(d.values(), reverse=True)
    if len(v) < 2: continue
    gap = v[0] - v[1]
    if gap == 0: ties += 1
    elif gap < 50: near += 1
    else: clear += 1
print("queries:", len(sc))
print("  exact ties (top two identical):", ties)
print("  gap under 50 bits            :", near)
print("  clearly separated            :", clear)
