from collections import defaultdict

tot = defaultdict(lambda: defaultdict(float))
for line in open("blast_work/forward.tsv"):
    q, s, pid, ln, b = line.rstrip("\n").split("\t")
    tot[q][s] += float(b)

print("WHERE QUERIES LAND")
counts = defaultdict(int)
for q, d in tot.items():
    counts[max(d.items(), key=lambda kv: kv[1])[0]] += 1
for s, c in sorted(counts.items(), key=lambda kv: -kv[1]):
    print("  %4d -> %s" % (c, s[:75]))

print("\nRELATIVE MARGIN (top vs second, as %% of top)")
rels = []
for q, d in tot.items():
    r = sorted(d.values(), reverse=True)
    if len(r) > 1 and r[0] > 0:
        rels.append(100.0 * (r[0] - r[1]) / r[0])
rels.sort()
if rels:
    print("  median %.3f%%   min %.3f%%   max %.3f%%" % (
        rels[len(rels)//2], rels[0], rels[-1]))

print("\nEXAMPLE: full ranking for 2 queries")
for q in list(tot)[:2]:
    print("  " + q[:60])
    for s, b in sorted(tot[q].items(), key=lambda kv: -kv[1]):
        print("      %10.0f  %s" % (b, s[:65]))
