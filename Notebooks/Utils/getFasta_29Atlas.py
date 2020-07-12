#!/usr/bin/env python
# coding: utf-8

# In[41]:


import requests, sys
import pandas as pd

data = pd.read_csv('../../Data/29Atlas/QuantificationTables/Table_EV3/Table_EV3.tsv',sep='\t')

server = "https://rest.ensembl.org"

missing = []

#pd.DataFrame(columns=['gene','status'])

for gene in data['EnsemblGeneID']:
    
    last = open("lastgene.txt", "w")
    last.write(gene)
    last.close()
    
    ext = "/sequence/id/"+gene+"?type=genomic;expand_5prime=2000;expand_3prime=2000"
    r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
    open(gene+'.fasta', 'wb').write(r.content)
    
    print(gene)
    print(gene,r.status_code)
    if not r.ok:
        missing.append(gene,r.status_code)

if len(missing)>0:
    missingDF = pd.DataFrame(missing, columns=['gene','status'])
    missingDF.set_index('gene',inplace=True)
    missingDF.to_csv('missing_29Atlas.tsv',sep='\t')
        
        
        
        
#r.raise_for_status()
#sys.exit() 
#print(r.text)

