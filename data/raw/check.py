import pandas as pd
ppi = pd.read_csv("BIOGRID-ORGANISM-Homo_sapiens-4.4.248.tab3.txt", sep='\t')

protein = ppi[['Official Symbol Interactor A','Official Symbol Interactor B',
			   'SWISS-PROT Accessions Interactor A','TREMBL Accessions Interactor A',
			   'SWISS-PROT Accessions Interactor B','TREMBL Accessions Interactor B']]

protein.drop_duplicates().to_csv('protein.csv',index=False)