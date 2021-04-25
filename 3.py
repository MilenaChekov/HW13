import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu

from rpy2 import robjects
from rpy2.robjects import Formula

from rpy2.robjects import pandas2ri
pandas2ri.activate()

from rpy2.robjects.packages import importr

base = importr("base")
stats = importr("stats")
DESeq2 = importr("DESeq2")

# Load read counts table
counts = pd.read_csv("colon_cancer_tumor_vs_normal_unpaired_counts.tsv", sep="\t", index_col=0)

# Define meta
meta = pd.DataFrame({"Tissue": ["Tumor"]*5 + ["Normal"]*5}, index=counts.columns)
meta["Tissue"] = stats.relevel(robjects.vectors.FactorVector(meta["Tissue"]), ref="Normal")

# Calculate normalization factors
dds = DESeq2.DESeqDataSetFromMatrix(countData=counts, colData=meta, design=Formula("~ Tissue"))
dds = DESeq2.DESeq(dds)

res = DESeq2.results(dds, name="Tissue_Tumor_vs_Normal")
res = DESeq2.lfcShrink(dds, coef="Tissue_Tumor_vs_Normal", type="apeglm")
res = pd.DataFrame(base.as_data_frame(res))
res.index = counts.index
res = res.sort_values("padj")
res = res.loc[res["padj"] < 0.05]
#res = res.loc[res["log2FoldChange"].abs() >= 1]

print(res.iloc[:10].index)

'''Index(['FABP6', 'ETV4', 'IER5L', 'KRT80', 'FUT1', 'C17orf96', 'CLDN1', 'ATG9B',
       'KIAA1257', 'SLC51B'],
      dtype='object')'''

df = pd.read_csv("colon_cancer_tumor_vs_normal_unpaired_FPKM.tsv", sep="\t", index_col=0)
df["ttest"] = [ttest_ind(df.loc[gene].iloc[0:5], df.loc[gene].iloc[5:10])[1] for gene in df.index]
df = df.sort_values("ttest")
print('ttest:', df.iloc[:10].index)

df = pd.read_csv("colon_cancer_tumor_vs_normal_unpaired_FPKM.tsv", sep="\t", index_col=0)
df["mannwhitneyu"] = [mannwhitneyu(df.loc[gene].iloc[0:5], df.loc[gene].iloc[5:10])[1] for gene in df.index]
df = df.sort_values("mannwhitneyu")
print('mannwhitneyu:', df.iloc[:10].index)

'''ttest: Index(['C17orf96', 'IER5L', 'FUT1', 'CDH3', 'FXYD5', 'ZNHIT2', 'CLCA4',
       'ACADSB', 'MT1F', 'PIGN'],
      dtype='object')
mannwhitneyu: Index(['SFTA2', 'RAET1L', 'CTD-2147F2.1', 'LINC00460', 'AC007128.1',
       'RP5-884M6.1', 'VAC14-AS1', 'RP11-399O19.9', 'CST1', 'LINC00858'],
      dtype='object')'''

deseq2=set(['FABP6', 'ETV4', 'IER5L', 'KRT80', 'FUT1', 'C17orf96', 'CLDN1', 'ATG9B', 'KIAA1257', 'SLC51B'])
ttest=set(['C17orf96', 'IER5L', 'FUT1', 'CDH3', 'FXYD5', 'ZNHIT2', 'CLCA4', 'ACADSB', 'MT1F', 'PIGN'])
mann=set(['SFTA2', 'RAET1L', 'CTD-2147F2.1', 'LINC00460', 'AC007128.1', 'RP5-884M6.1', 'VAC14-AS1', 'RP11-399O19.9', 'CST1', 'LINC00858'])

print('deseq2+ttest', deseq2.intersection(ttest))
print('deseq2+mann', deseq2.intersection(mann))
print('ttest+mann', ttest.intersection(mann))

'''deseq2+ttest {'C17orf96', 'IER5L', 'FUT1'}
deseq2+mann set()
ttest+mann set()'''
