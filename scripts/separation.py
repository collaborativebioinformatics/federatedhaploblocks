import re, subprocess
from collections import defaultdict

CRE = re.compile(r"^(BIG|SMALL|SINGLE)_(chr\d+_[\d-]+)_size(\d+)_")

meta = {}
for line in open("consensus_samples.fasta"):
    if line.startswith(">"):
        cid = line[1:].split()[0]
        m = CRE.match(cid)
        if m:
            meta[cid] = (m.group(2), m.group(1), int(m.group(3)))

subprocess.run(["makeblastdb", "-in", "consensus_samples.fasta",
                "-dbtype", "nucl", "-out", "/tmp/allcons"],
               capture_output=True)
r = subprocess.run(["blastn", "-query", "consensus_samples.fasta",
                    "-db", "/tmp/allcons",
                    "-outfmt", "6 qseqid sseqid pident length mismatch gaps",
                    "-dust", "no", "-qcov_hsp_perc", "80"],
                   capture_output=True, text=True)

seen = set()
rows = defaultdict(list)
for line in r.stdout.splitlines():
    q, s, pid, ln, mm, gaps = line.split("\t")
    if q == s or q not in meta or s not in meta:
        continue
    if meta[q][0] != meta[s][0]:
        continue
    key = tuple(sorted([q, s]))
    if key in seen:
        continue
    seen.add(key)
    rows[meta[q][0]].append((meta[q][1], meta[s][1], float(pid),
                             int(ln), int(mm), int(gaps)))

print("%-26s %-16s %9s %10s %6s" % ("haploblock", "pair", "identity", "aligned", "mism"))
print("-" * 74)
for region in sorted(rows):
    for a, b, pid, ln, mm, gaps in sorted(rows[region]):
        print("%-26s %-16s %8.3f%% %10s %6d" % (
            region, a + " vs " + b, pid, "{:,}".format(ln), mm))
