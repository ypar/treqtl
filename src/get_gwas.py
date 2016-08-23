#!/usr/bin/env python

###
# ypar
# @ treqtl/get_gwas.py
# format gwas summary statistics for mr
###

# mreqtl header
# RSNUM SNPID CHR POS_HG19 A1 A2 INC_ALLELE INC_AFRQ BETA SE PVAL

import pandas as pd
import numpy as np
import subprocess
import sys, os
from sys import argv



# CONVERGE MDD header (hg19)
# CHR BP.GRCh37 RSID REF_ALLELE ALT_ALLELE FRQ_ALT.1KGP_ASN_n286 INFO.plink OR.logistic SE P.lmm

def format_converge():
  
  subprocess.call('wget http://www.med.unc.edu/pgc/files/resultfiles/converge.MDD.summary_stats.2Sep2015.tbl.gz', shell=True)
  infile = 'converge.MDD.summary_stats.2Sep2015.tbl.gz'
  
  df = pd.read_csv(infile, sep=' ', header=0, names=['CHR', 'POS_HG19', 'RSNUM','A2','A1', 'INC_AFRQ', 'INFO', 'OR', 'SE', 'PVAL'], na_values='NA', dtype=object, compression='gzip')
  
  df['CHR'] = 'chr' + df['CHR'].astype(str)
  df['SNPID'] = df['CHR'].astype(str) + ':' + df['POS_HG19'].astype(str)
  df.loc[df['INC_AFRQ'] < 0, 'INC_AFRQ'] = 'NA'
  df = df[np.isfinite(df['OR'].astype(float))]
  df['BETA'] = np.log(df['OR'].astype(float))
  df = df[['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'A1', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']]
  df.names = ['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'INC_ALLELE', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']
  
  df.to_csv('converge_for_treqtl.txt', sep='\t', index=False)


# GCAN header (hg18)
# chromosome	position	SNP	reference_allele	other_allele	eaf	OR	OR_se	OR_95L	OR_95U	z	p_sanger	_-log10_p-value	q_statistic	q_p-value	i2	n_studies	n_samples	effects

def format_gcan():
  
  subprocess.call('wget http://www.med.unc.edu/pgc/files/resultfiles/gcan_meta-out.gz', shell=True)
  infile = 'gcan_meta-out.gz'
  
  # function to match snpIDs is deprecated in treqtl.
  df = pd.read_csv(infile, sep='\t', header=0, names=['CHR', 'POS_HG18', 'RSNUM','A2','A1', 'INC_AFRQ', 'OR', 'OR_se', 'OR_95L', 'OR_95U', 'z', 'PVAL', '_-log10_p-value', 'q_statistic', 'q_p-value', 'i2', 'n_studies', 'n_samples', 'effects'], compression='gzip')
  df['CHR'] = 'chr' + df['CHR'].astype(str)
  df['SNPID'] = df['CHR'].astype(str) + ':' + df['POS_HG18'].astype(str)
  df.loc[df['INC_AFRQ'] < 0, 'INC_AFRQ'] = 'NA'
  df['BETA'] = np.log(df['OR'].astype(float))
  df['SE'] = np.log(df['OR'].astype(float))/df['z'].astype(float)
  df = df[['RSNUM', 'SNPID', 'CHR', 'POS_HG18', 'A1', 'A2', 'A1', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']]
  df.names = ['RSNUM', 'SNPID', 'CHR', 'POS_HG18', 'A1', 'A2', 'INC_ALLELE', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']
  
  df.to_csv('gcan_for_treqtl.txt', sep='\t', index=False)



# RA header (hg18)
# CHR SNP POS major_al minor_al wtccc_info narac1_info narac2_info eira_info canada_info brass_info wtccc_z narac1_z narac2_z eira_z canada_z brass_z cases_MM cases_Mm cases_mm controls_MM controls_Mm controls_mm meta_OR OR_95%CI_lo OR_95%CI_up meta_z meta_2tP CochranQ Q_Pval

def format_ra():
  
  subprocess.call('wget http://www.broadinstitute.org/ftp/pub/rheumatoid_arthritis/Stahl_etal_2010NG/RA_GWASmeta2_20090505-results.txt', shell=True)
  infile = 'RA_GWASmeta2_20090505-results.txt'
  
  df = pd.read_csv(infile, sep=' ', header=0, names=['CHR', 'RSNUM', 'POS_HG18', 'A2', 'A1', 'wtccc_info', 'narac1_info', 'narac2_info', 'eira_info', 'canada_info', 'brass_info', 'wtccc_z', 'narac1_z', 'narac2_z', 'eira_z', 'canada_z', 'brass_z', 'cases_MM', 'cases_Mm', 'cases_mm', 'controls_MM', 'controls_Mm', 'controls_mm', 'meta_OR', 'OR_95%CI_lo', 'OR_95%CI_up', 'meta_z', 'PVAL', 'CochranQ', 'Q_Pval'], na_values='NA', dtype=object)
  df['CHR'] = 'chr' + df['CHR'].astype(str)
  df['SNPID'] = df['CHR'].astype(str) + ':' + df['POS_HG18'].astype(str)
  df['INC_AFRQ'] = 'NA'
  df = df[df['meta_OR']!='1.00']
  df['BETA'] = np.log(df['meta_OR'].astype(float))
  df['SE'] = np.log(df['meta_OR'].astype(float))/df['meta_z'].astype(float)
  df = df[['RSNUM', 'SNPID', 'CHR', 'POS_HG18', 'A1', 'A2', 'A1', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']]
  df.columns = ['RSNUM', 'SNPID', 'CHR', 'POS_HG18', 'A1', 'A2', 'INC_ALLELE', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']
  
  df.to_csv('ra_for_treqtl.txt', sep='\t', index=False)



# GPC header
# rsid chr bp a1 a2 beta se pvalue info ncoh maf

def format_neuroticism2():
  
  subprocess.call('wget https://www.dropbox.com/s/ql1en6s0kramxne/GPC-2.NEUROTICISM.zip', shell=True)
  subprocess.call('unzip GPC-2.NEUROTICISM.zip', shell=True)
  infile = 'GPC-2.NEUROTICISM.full.txt'
  
  df = pd.read_csv(infile, sep='\t', header=None, names=['RSNUM', 'CHR', 'POS_HG19', 'A1', 'A2', 'BETA', 'SE', 'PVAL', 'INFO', 'NCOH', 'INC_AFRQ'])
  df['CHR'] = 'chr' + df['CHR'].astype(str)
  df['SNPID'] = df['CHR'].astype(str) + ':' + df['POS_HG19'].astype(str)
  df['A1'] = df['A1'].str.upper()
  df['A2'] = df['A2'].str.upper()
  df = df[['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'A1', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']]
  df.names = ['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'INC_ALLELE', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']
  
  df.to_csv('neuroticism2_for_treqtl.txt', sep='\t', index=False)


def format_extraversion2():
  
  subprocess.call('wget https://www.dropbox.com/s/bk2jn41vrfl3zna/GPC-2.EXTRAVERSION.zip', shell=True)
  subprocess.call('unzip GPC-2.EXTRAVERSION.zip', shell=True)
  infile = 'GPC-2.EXTRAVERSION.full.txt'
  
  df = pd.read_csv(infile, sep='\t', header=None, names=['RSNUM', 'CHR', 'POS_HG19', 'A1', 'A2', 'BETA', 'SE', 'PVAL', 'INFO', 'NCOH', 'INC_AFRQ'])
  df['CHR'] = 'chr' + df['CHR'].astype(str)
  df['SNPID'] = df['CHR'].astype(str) + ':' + df['POS_HG19'].astype(str)
  df['A1'] = df['A1'].str.upper()
  df['A2'] = df['A2'].str.upper()
  df = df[['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'A1', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']]
  df.names = ['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'INC_ALLELE', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']
  
  df.to_csv('extraversion2_for_treqtl.txt', sep='\t', index=False)


# heart rate header
# rsid	CHR	POS	EFFECT_ALLELE	OTHER_ALLELE	EAF	beta	SE	N_TOTAL	P_VALUE	P_VALUE_GCADJ

def format_heartrate():
  
  subprocess.call('wget http://walker05.u.hpc.mssm.edu/META_STAGE1_GWASHR_SUMSTATS.txt', shell=True)
  infile = 'META_STAGE1_GWASHR_SUMSTATS.txt'

  df = pd.read_csv(infile, sep='\t', header=0, names=['RSNUM', 'CHR', 'POS_HG18', 'A1', 'A2', 'INC_AFRQ', 'BETA', 'SE', 'N', 'PVAL', 'PVAL_GCADJ'], dtype=object)
  df['CHR'] = 'chr' + df['CHR'].astype(str)
  df['SNPID'] = df['CHR'].astype(str) + ':' + df['POS_HG18'].astype(str)
  df['A1'] = df['A1'].str.upper()
  df['A2'] = df['A2'].str.upper()
  df = df[['RSNUM', 'SNPID', 'CHR', 'POS_HG18', 'A1', 'A2', 'A1', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']]
  df.names = ['RSNUM', 'SNPID', 'CHR', 'POS_HG18', 'A1', 'A2', 'INC_ALLELE', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']
  
  outfile = 'heartrate_for_treqtl.txt'
  df.to_csv(outfile, sep='\t', index=False)
  

# glgc header
# SNP_hg18	SNP_hg19	rsid	A1	A2	beta	se	N	P-value	Freq.A1.1000G.EUR

# mreqtl header
# RSNUM SNPID CHR POS_HG19 A1 A2 INC_ALLELE INC_AFRQ BETA SE PVAL

def format_ldlc():

  subprocess.call('wget http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_LDL.txt.gz', shell=True)
  infile = 'jointGwasMc_LDL.txt.gz'
  
  df = pd.read_csv(infile, sep='\t', header=0, names=['SNP_HG18', 'SNP_HG19', 'RSNUM', 'A1', 'A2', 'BETA', 'SE', 'N', 'PVAL', 'INC_AFRQ'], compression='gzip')
  add = pd.DataFrame(df.SNP_HG19.str.split(':',1).tolist(), columns = ['CHR','POS_HG19'])

  df['CHR'] = add['CHR']
  df['POS_HG19'] = add['POS_HG19']
  df['SNPID'] = df['SNP_HG19']
  df['A1'] = df['A1'].str.upper()
  df['A2'] = df['A2'].str.upper()
  df = df[['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'A1', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']]
  df.names = ['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'INC_ALLELE', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']
  
  outfile = 'ldlc_for_treqtl.txt'
  df.to_csv(outfile, sep='\t', index=False)
  
  
def format_hdlc():

  subprocess.call('wget http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_HDL.txt.gz', shell=True)
  infile = 'jointGwasMc_HDL.txt.gz'
  
  df = pd.read_csv(infile, sep='\t', header=0, names=['SNP_HG18', 'SNP_HG19', 'RSNUM', 'A1', 'A2', 'BETA', 'SE', 'N', 'PVAL', 'INC_AFRQ'], compression='gzip')
  add = pd.DataFrame(df.SNP_HG19.str.split(':',1).tolist(), columns = ['CHR','POS_HG19'])

  df['CHR'] = add['CHR']
  df['POS_HG19'] = add['POS_HG19']
  df['SNPID'] = df['SNP_HG19']
  df['A1'] = df['A1'].str.upper()
  df['A2'] = df['A2'].str.upper()
  df = df[['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'A1', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']]
  df.names = ['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'INC_ALLELE', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']
  
  outfile = 'hdlc_for_treqtl.txt'
  df.to_csv(outfile, sep='\t', index=False)


def format_tg():

  subprocess.call('wget http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_TG.txt.gz', shell=True)
  infile = 'jointGwasMc_HDL.txt.gz'
  
  df = pd.read_csv(infile, sep='\t', header=0, names=['SNP_HG18', 'SNP_HG19', 'RSNUM', 'A1', 'A2', 'BETA', 'SE', 'N', 'PVAL', 'INC_AFRQ'], compression='gzip')
  add = pd.DataFrame(df.SNP_HG19.str.split(':',1).tolist(), columns = ['CHR','POS_HG19'])

  df['CHR'] = add['CHR']
  df['POS_HG19'] = add['POS_HG19']
  df['SNPID'] = df['SNP_HG19']
  df['A1'] = df['A1'].str.upper()
  df['A2'] = df['A2'].str.upper()
  df = df[['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'A1', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']]
  df.names = ['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'INC_ALLELE', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']
  
  outfile = 'tg_for_treqtl.txt'
  df.to_csv(outfile, sep='\t', index=False)
  

def format_tc():

  subprocess.call('wget http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_TC.txt.gz', shell=True)
  infile = 'jointGwasMc_HDL.txt.gz'
  
  df = pd.read_csv(infile, sep='\t', header=0, names=['SNP_HG18', 'SNP_HG19', 'RSNUM', 'A1', 'A2', 'BETA', 'SE', 'N', 'PVAL', 'INC_AFRQ'], compression='gzip')
  add = pd.DataFrame(df.SNP_HG19.str.split(':',1).tolist(), columns = ['CHR','POS_HG19'])

  df['CHR'] = add['CHR']
  df['POS_HG19'] = add['POS_HG19']
  df['SNPID'] = df['SNP_HG19']
  df['A1'] = df['A1'].str.upper()
  df['A2'] = df['A2'].str.upper()
  df = df[['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'A1', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']]
  df.names = ['RSNUM', 'SNPID', 'CHR', 'POS_HG19', 'A1', 'A2', 'INC_ALLELE', 'INC_AFRQ', 'BETA', 'SE', 'PVAL']
  
  outfile = 'tc_for_treqtl.txt'
  df.to_csv(outfile, sep='\t', index=False)


def parse_getgwas(dataset):
  
  dataset = str(dataset)
  
  if dataset == 'converge':
    format_converge()
  elif dataset == 'gcan':
    format_gcan()
  elif dataset == 'ra':
    format_ra()
  elif dataset == 'neuroticism2':
    format_neuroticism2()
  elif dataset == 'extraversion2':
    format_extraversion2()
  elif dataset == 'heartrate':
    format_heartrate()
  elif dataset == 'ldlc':
    format_ldlc()
  elif dataset == 'hdlc':
    format_hdlc()
  elif dataset == 'tg':
    format_tg()
  elif dataset == 'tc':
    format_tc()
  
  else:
    logmessage = str(dataset) + ' is not a recognized dataset'
    print(logmessage)
    sys.exit()
    



if __name__ == '__main__':
  
  if len(sys.argv) < 2:
    logmessage = 'analysis terminated. input files are not specified'
    print(logmessage)
    
  else:
    parse_getgwas(argv[1])

  
