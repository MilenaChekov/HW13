import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

counts = pd.read_csv("colon_cancer_tumor_vs_normal_unpaired_counts.tsv", sep="\t", index_col=0)
deseq2 = counts.loc[['FABP6', 'ETV4', 'IER5L', 'KRT80', 'FUT1', 'C17orf96', 'CLDN1', 'ATG9B', 'KIAA1257', 'SLC51B']].T
print(deseq2)
sns.violinplot(data=deseq2)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("violin_deseq2.png")
plt.close()

fpkm = pd.read_csv("colon_cancer_tumor_vs_normal_unpaired_FPKM.tsv", sep="\t", index_col=0)

ttest = fpkm.loc[['C17orf96', 'IER5L', 'FUT1', 'CDH3', 'FXYD5', 'ZNHIT2', 'CLCA4', 'ACADSB', 'MT1F', 'PIGN']].T
sns.violinplot(data=ttest)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("violin_ttest.png")
plt.close()

mann = fpkm.loc[['SFTA2', 'RAET1L', 'CTD-2147F2.1', 'LINC00460', 'AC007128.1', 'RP5-884M6.1', 'VAC14-AS1', 'RP11-399O19.9', 'CST1', 'LINC00858']].T
sns.violinplot(data=mann)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("violin_mann.png")
plt.close()
