# trait-associated eQTL (treQTL) mapping analysis module
-------
### YoSon Park
-------
## Citations and notes

manuscript citation:
Park Y, Voight BF, Engelhardt BE and Brown CD (2016), Mendelian randomization to study genetic mechanisms of complex traits, (in preparation)  

code citation:
http://github.com/ypar/treqtl.git

-------
Examples analysis flows involving publicly available datasets as well as useful relevant scripts are provided in the notebook  
<http://nbviewer.jupyter.org/gist/ypar/06c185759643c436cce5>  

-------
## Installation

Scripts in the source code can be used without installation  

wget https://github.com/ypar/treqtl/archive/master.zip  

once downloaded, you may call the main prompt by indicating the path  

\<dir path to treqtl\>/treqtl  

If preferred, export path to the treqtl module to access from anywhere  

export PATH=$PATH:\<dir path to treqtl\>  

If necessary, the repo can also be cloned  

git clone http://github.com/ypar/treqtl.git  


This module may also be used as a part of analyses pipeline by importing to python 3 environment. However, this is still under active development and users may use these at their own risks.    

To install  

python3 setup.py install  
Or  
pip install treqtl  



## Input files


### Input file option 1: merged GWAS X eQTL summary statistics

This method is a lot faster than the option 2 and is recommended for large-ish datasets.  
The xtreqtl input file example  

GENE	RSNUM	CHR	POS_HG19	A1_e	A2_e	INC_AFRQ_e	BETA_e	SE_e	PVAL_e	SNPID	A1_g	A2_g	INC_ALLELE	INC_AFRQ_g	BETA_g	SE_g	PVAL_g  
ENSG00000256448.1	rs12421053	11	73276900	T	C	0.2781065	0.196160052625968	0.043988824	1.1687065744466196e-05	chr11:72954548	T	C	T	0.2428	0.02979999966919422	0.01600000075995922	0.06267999857664108  
ENSG00000214826.4	rs2003608	12	9519429	T	C	0.42159763	0.212636884718404	0.03152041	7.98540667976702e-11	chr129410696	T	C	T	0.5787	0.027699999511241913  


The input data is presumed to be matched per row and that header is consistent with the indicated file example above  


### Input file option 2: separate GWAS and eQTL summary statistics files

eqtl input file example   

GENE	RSNUM	CHR	POS_HG19	A1	A2	INC_ALLELE	INC_AFRQ	BETA	SE	PVAL  
GPR17   rs17262104      chr2    128747549       G       T       G       0.06456953      0.736983583560685       0.11432743      5.541546921310881e-10  
HAX1    rs12749691      chr1    154251259       T       A       T       0.27152318      0.280817746771117       0.05622703      1.08387876813775e-06  


gwas input file example  

RSNUM	SNPID	CHR	POS_HG19	A1	A2	INC_ALLELE	INC_AFRQ	BETA	SE	PVAL  
rs1172982	chr1:100230111	chr1	100230111	T	C	T	0.3219	0.0043	0.0055	0.4689  
rs1172981	chr1:100230197	chr1	100230197	T	C	T	0.06069	0.0057	0.0103	0.7688  



The pvalue column is not used in the analysis but retained for convenience of results interpretations (see relevant notes in the ipython notebook). In practice, this column may be used to contain notes and comments regarding the region or additional annotations such as alternate gene and transcript names, etc.  


## Usage

example usage:  
for data trimmed to contain only independent variants per transcript, --treqtl will sum up treqtl results per transcript  
./treqtl -e example/Liver_eqtl.txt -g example/GLGC_LDLC_chr1.txt --treqtl  

for data with non-independent set of variants, each variant may be used for analysis using --itreqtl  


## Other options

for users' convenience, select publicly availble datasets can be downloaded using --getgwas argument for treqtl analyses.  

currently only following options are available in treqtl  
- rheumatoid arthritis gwas dataset (stahl et al 2010)  
- anorexia nervosa gwas dataset (GCAN)  
- major depressive disorder gwas dataset (CONVERGE)  
- neuroticism gwas dataset (GPC phase 2)  
- extraversion gwas dataset (GPC phase 2)  
- heart rate gwas dataset (den hoed et al 2013)  
- low-density lipid cholesterol levels gwas dataset (GLGC)  
- high-density lipid choleterol levels gwas dataset (GLGC)  
- total cholesterol levels gwas dataset (GLGC)  
- tryglyceride levels gwas dataset (GLGC)  

While these files are formatted for treqtl analysis, input files may require additional formatting/filtering for each analysis specifically  



