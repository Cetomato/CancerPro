CancerPro: Deciphering the Pan-Cancer Prognostic Landscape through Combinatorial Enrichment Analysis and Knowledge Network Insights

CancerPro knowledge network insight platform
(https://medcode.link/cancerpro). This platform integrates extensive biomedical data on drugs,
genes, pathways, and diseases, along with interaction information such as protein-protein
interactions, gene-phenotype relations, and disease associations.

We specifically designed three functional modules to analyze pan-cancer prognosis genes comprehensively to meet our research objectives. 
The first module, named X-enrich, focuses on functional enrichment analysis. 
We use the Over-Representation Analysis method to test whether known biological functions or processes are overrepresented in a given gene list. 
Input genes are first mapped to the selected annotation sets. 
Subsequently, a hypergeometric test is conducted to identify overrepresented terms within these sets, 
using all genes associated with the selected annotations as a reference. 
Since the annotation information encompasses not only functional annotations such as GO and pathways but also includes regulations of drugs, diseases, 
phenotypes, etc., we can analyze whether a gene list is enriched in these diverse entities. 

The second module focuses on a novel perspective on drug called Drug Clue. By considering drug target proteins, genes affected by drugs, 
genes upregulated or downregulated by drugs, Drug Clue can identify the biological pathways, diseases, and phenotypes associated with a particular drug. 
It can also identify drugs that may have similar effects. This information can be used to better understand how drugs work and to identify potential 
side effects. Using this Drug Clue module, we gain insight into drugs that significantly downregulate pan-cancer prognostic unfavorable genes. 

The third section introduces a new method for examining gene lists known as Gene List (GL) Insight. Through this module, genes are interconnected using 
highly reliable protein-protein interaction information or protein regulation information. Subsequently, we extract drugs associated with at least two or 
more genes, which can target these genes or affect/regulate their expression. Additionally, the module can retrieve biological pathways from KEGG or 
Reactome, GO annotations, and disease/phenotype information linked to multiple genes. Furthermore, the module offers support for selecting critical nodes 
based on the 4 centrality measurements. We conducted a detailed analysis of prognostically unfavorable genes specific to prostate cancer 
in this module. Using HINT, we connected them and retrieved Reactome pathway information. 

DOI: 10.5281/zenodo.14010461

Email: wangzg@pumc.edu.cn
