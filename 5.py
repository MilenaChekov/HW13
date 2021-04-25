import random
import pandas as pd

all_genes = pd.read_csv('all_genes.txt', sep='\t', header=None, index_col=0)
print(all_genes)

for i in range(5):
    all_genes.sample(1000).to_csv('all_genes'+str(i)+'.txt', header=None, sep='\t')

