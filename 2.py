import pandas as pd
from rpy2 import robjects
from rpy2.robjects import Formula

from rpy2.robjects import pandas2ri
pandas2ri.activate()

from rpy2.robjects.packages import importr

base = importr("base")
stats = importr("stats")
DESeq2 = importr("DESeq2")

# Load read counts table
counts = pd.read_csv("colon_cancer_tumor_vs_normal_paired_counts.tsv", sep="\t", index_col=0)

# Define meta
meta = pd.DataFrame({"Tissue": ["Tumor"]*5 + ["Normal"]*5, "Patient": list(range(1, 6))*2}, index=counts.columns)
meta["Tissue"] = stats.relevel(robjects.vectors.FactorVector(meta["Tissue"]), ref="Normal")

# Define meta
meta2 = pd.DataFrame({"Tissue": ["Tumor"]*5 + ["Normal"]*5}, index=counts.columns)
meta2["Tissue"] = stats.relevel(robjects.vectors.FactorVector(meta2["Tissue"]), ref="Normal")

# Calculate normalization factors
dds = DESeq2.DESeqDataSetFromMatrix(countData=counts, colData=meta, design=Formula("~ Patient + Tissue"))
dds = DESeq2.DESeq(dds)

# Calculate normalization factors
dds2 = DESeq2.DESeqDataSetFromMatrix(countData=counts, colData=meta2, design=Formula("~ Tissue"))
dds2 = DESeq2.DESeq(dds2)

res = DESeq2.results(dds, name="Tissue_Tumor_vs_Normal")
res = DESeq2.lfcShrink(dds, coef="Tissue_Tumor_vs_Normal", type="apeglm")
res = pd.DataFrame(base.as_data_frame(res))
res.index = counts.index
res = res.sort_values("padj")
res = res.loc[res["padj"] < 0.05]
#res = res.loc[res["log2FoldChange"].abs() >= 1]

res2 = DESeq2.results(dds2, name="Tissue_Tumor_vs_Normal")
res2 = DESeq2.lfcShrink(dds2, coef="Tissue_Tumor_vs_Normal", type="apeglm")
res2 = pd.DataFrame(base.as_data_frame(res2))
res2.index = counts.index
res2 = res2.sort_values("padj")
res2 = res2.loc[res2["padj"] < 0.05]
#res = res.loc[res["log2FoldChange"].abs() >= 1]

print(len(res))
print(res.iloc[:10].index)
res.to_csv("DESeq2_results_paired.tsv", sep="\t")

print(len(res2))
print(res2.iloc[:10].index)
res2.to_csv("DESeq2_results_unpaired.tsv", sep="\t")