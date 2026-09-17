import pandas as pd

df = pd.read_csv("block_stats.tsv", sep="\t")
N = int(df["n_haplotypes"].max())

print("blocks: {:,}     haplotypes: {:,}  ({:,} individuals)".format(
    len(df), N, N // 2))

tot = int(df["singleton_count"].sum())
print("\nSINGLETON EXPOSURE")
print("  total singleton clusters genome-wide : {:,}".format(tot))
print("  per haplotype, mean blocks as singleton: {:.1f}".format(tot / N))
print("  per individual (2 haplotypes)          : {:.1f}".format(2 * tot / N))

has = (df["singleton_count"] > 0).sum()
print("  blocks containing at least one singleton: {:,} of {:,}  ({:.1f}%)".format(
    has, len(df), 100 * has / len(df)))

print("\n  singleton_count per block:")
for q in (0.5, 0.75, 0.9, 0.99, 1.0):
    print("    {:>5.0%}  {:>8.0f}".format(q, df["singleton_count"].quantile(q)))

d = df.sort_values("singleton_count", ascending=False).reset_index(drop=True)
d["cum"] = d["singleton_count"].cumsum()
print("\nSUPPRESSING WHOLE BLOCKS, worst first")
print("  {:>8}  {:>14}  {:>10}".format("blocks", "singletons gone", "% removed"))
for k in (10, 100, 500, 1000, 5000, 10000):
    if k <= len(d):
        c = int(d.loc[k - 1, "cum"])
        print("  {:>8,}  {:>14,}  {:>9.1f}%".format(k, c, 100 * c / tot))

need = int((d["cum"] >= tot).idxmax()) + 1
print("  blocks needed to remove every singleton: {:,} ({:.1f}% of all blocks)".format(
    need, 100 * need / len(d)))

print("\nSTRUCTURE")
print("  dominance      median {:.3f}   min {:.4f}   max {:.4f}".format(
    df["dominance"].median(), df["dominance"].min(), df["dominance"].max()))
print("  n_clusters     median {:.0f}   max {:,}".format(
    df["n_clusters"].median(), int(df["n_clusters"].max())))
print("  entropy        median {:.2f}   max {:.2f}".format(
    df["shannon_entropy"].median(), df["shannon_entropy"].max()))

df["singleton_frac"] = df["singleton_count"] / df["n_haplotypes"]
w = df.nlargest(8, "singleton_frac")[
    ["block", "n_clusters", "max_cluster_size", "singleton_count",
     "dominance", "shannon_entropy"]]
print("\nMOST EXPOSED BLOCKS (highest singleton fraction)")
print(w.to_string(index=False))

corr = df[["singleton_count", "dominance", "shannon_entropy",
           "n_clusters", "block_length"]].corr()["singleton_count"]
print("\nCORRELATION with singleton_count")
for k, v in corr.items():
    if k != "singleton_count":
        print("  {:<16} {:+.3f}".format(k, v))
