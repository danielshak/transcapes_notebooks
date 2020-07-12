-data.ipynb
This file takes in the raw rna_seq reads as well as the proteomics tables. This does some filtering to remove missing values etc.

-gene2protein.ipynb
In the gene2protein mappings, a list of protein names is generated. This list of protein names is used to fetch protein sequences from uniprot (such that the first name in the uniprot.txt file corresponds to the first sequence from the protein sequence file). However the order of the names generated may change if the mappings (gene2protein.ipynb) is run again due to the nature of set comparisons/sorting etc... If the file is run again the protein sequences should be refetched from the uniprot database to match the names in the uniprot.txt file since the fetched sequences do not contain identifiers that match the names in the uniprot.txt file.

This file is used to fetch majority protein identifier for each gene and generate a gene to protein mapping for both AT2 and AM cell lines. This file generates a list of protein names that can be used to fetch protein sequences. 