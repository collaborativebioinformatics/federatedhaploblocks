#!/usr/bin/env python3
import random, re, subprocess, sys
from collections import defaultdict
from pathlib import Path

CONSENSUS = sys.argv[1] if len(sys.argv) > 1 else "consensus_samples.fasta"
QUERIES   = "full_clusters.fasta"
TSVS      = ["full_big_cluster_members.tsv", "full_small_cluster_members.tsv"]
SUBSAMPLE = 40
EVALUE    = "1e-20"
DUST      = "no"
WORKDIR   = Path("blast_work")
SEED      = 0

def run(cmd):
    print("  $", " ".join(str(c) for c in cmd))
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[:2000]); sys.exit("FAILED: " + str(cmd[0]))
    return r

def read_fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"): ids.append(line[1:].split()[0])
    return ids

def write_subset(src, keep, dest):
    keep = set(keep); out = open(dest, "w"); writing = False; n = 0
    with open(src) as fh:
        for line in fh:
            if line.startswith(">"):
                writing = line[1:].split()[0] in keep; n += writing
            if writing: out.write(line)
    out.close(); return n

def load_truth(paths):
    m2r = {}; size = defaultdict(int)
    for p in paths:
        if not Path(p).exists():
            print("  (missing " + p + ")"); continue
        with open(p) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 6: continue
                rep, member = f[4], f[5]
                m2r[member] = rep; size[rep] += 1; m2r.setdefault(rep, rep)
    return m2r, size

CRE = re.compile(r"^(BIG|SMALL|SINGLE)_.*?_size(\d+)_(.+)$")
def parse_consensus(cid):
    m = CRE.match(cid)
    return (m.group(1), int(m.group(2)), m.group(3)) if m else None

def blast(q, db, out):
    run(["blastn","-query",q,"-db",db,"-out",out,
         "-outfmt","6 qseqid sseqid pident length bitscore",
         "-evalue",EVALUE,"-dust",DUST,"-max_target_seqs","20","-qcov_hsp_perc","80","-num_threads","4"])

def parse_hits(path):
    tot = defaultdict(lambda: defaultdict(float))
    with open(path) as fh:
        for line in fh:
            q,s,pid,ln,bits = line.rstrip("\n").split("\t")
            tot[q][s] += float(bits)
    return {q: sorted(d.items(), key=lambda kv: -kv[1]) for q,d in tot.items()}

def main():
    WORKDIR.mkdir(exist_ok=True); random.seed(SEED)
    meta = {}; rep2c = {}
    for cid in read_fasta_ids(CONSENSUS):
        p = parse_consensus(cid)
        if not p: print("  ! unparsed header: " + cid); continue
        meta[cid] = p; rep2c[p[2]] = cid
    print(str(len(meta)) + " consensus sequences:")
    for cid,(cls,size,_) in meta.items(): print("    " + cls.ljust(7) + " size=" + str(size))

    m2r, _ = load_truth(TSVS)
    print("\n" + str(len(m2r)) + " members mapped from TSV")

    allq = read_fasta_ids(QUERIES); byc = defaultdict(list)
    for q in allq:
        cid = rep2c.get(m2r.get(q))
        if cid: byc[meta[cid][0]].append(q)
    print(str(len(allq)) + " queries, " + str(len(allq)-sum(len(v) for v in byc.values())) + " with no consensus")
    for cls,qs in byc.items(): print("    " + cls.ljust(7) + " " + str(len(qs)) + " available")

    chosen = []
    for cls,qs in byc.items():
        chosen += random.sample(qs, min(SUBSAMPLE,len(qs))) if SUBSAMPLE else qs
    qf = WORKDIR/"queries_subset.fasta"
    print("\nwrote " + str(write_subset(QUERIES, chosen, qf)) + " queries")

    print("\nFORWARD"); cdb = WORKDIR/"cdb"
    run(["makeblastdb","-in",CONSENSUS,"-dbtype","nucl","-out",cdb])
    fwd = WORKDIR/"forward.tsv"; blast(qf, cdb, fwd); forward = parse_hits(fwd)

    print("\nREVERSE"); qdb = WORKDIR/"qdb"
    run(["makeblastdb","-in",qf,"-dbtype","nucl","-out",qdb])
    rev = WORKDIR/"reverse.tsv"; blast(CONSENSUS, qdb, rev); reverse = parse_hits(rev)
    rbest = {c:h[0][0] for c,h in reverse.items() if h}

    st = defaultdict(lambda: {"n":0,"correct":0,"mutual":0,"none":0,"margins":[]})
    for q in chosen:
        truth = rep2c[m2r[q]]; s = st[meta[truth][0]]; s["n"] += 1
        h = forward.get(q, [])
        if not h: s["none"] += 1; continue
        top, tb = h[0]; second = h[1][1] if len(h) > 1 else 0.0
        s["margins"].append(tb - second)
        if top == truth: s["correct"] += 1
        if rbest.get(top) == q: s["mutual"] += 1

    print("\n" + "="*74)
    print("class      n   correct   mutual   no hit   med margin")
    print("-"*74)
    for cls in ("BIG","SMALL","SINGLE"):
        s = st.get(cls)
        if not s or not s["n"]: continue
        m = sorted(s["margins"]); med = m[len(m)//2] if m else 0.0
        print(cls.ljust(8) + str(s["n"]).rjust(4) + ("%8.1f%%" % (100*s["correct"]/s["n"]))
              + str(s["mutual"]).rjust(9) + str(s["none"]).rjust(9) + ("%13.1f" % med))
    print("="*74)

main()
