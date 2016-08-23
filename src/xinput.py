#!/usr/bin/env python3

###
# YoSon
# @treqtl/xinput.py
# produce xtreqtl input files by matching rs numbers from trait and iv summary statistics
###


import pandas as pd
import numpy as np
import os, sys
from sys import argv
from os import walk

from treqtl_input import read_dir


def xinput(ewkdir, efile, gwkdir):
  
  ewkdirs = ewkdir + '/'
  gwkdirs = gwkdir + '/'
  
  edf = pd.read_csv(efile, delim_whitespace=True, header=0, na_values='NA')
  efilename = efile.replace(ewkdirs, '')
  efilename = efilename.replace('.txt', '_w_')
  
  gfilelist = read_dir(gwkdir)
  
  for gfile in gfilelist:
    
    gdf = pd.read_csv(gfile, delim_whitespace=True, header=0, na_values='NA', names= ['RSNUM', 'SNPID', 'CHR', 'POS', 'A1', 'A2', 'INC_ALLELE', 'INC_AFRQ', 'BETA', 'SE', 'PVAL'], dtype = {'RSNUM': str, 'SNPID': str, 'CHR': str, 'POS': str, 'A1': str, 'A2': str, 'INC_ALLELE': str, 'INC_AFRQ': str, 'BETA': np.float32, 'SE': np.float32, 'PVAL': np.float32})
    gdf['CHR'] = gdf['CHR'].replace('.0', '')
    gdf['POS'] = gdf['POS'].replace('.0', '')
    gfilename = gfile.replace(gwkdirs, '')
    
    outfile = str(efilename) + str(gfilename)
    
    merged = edf.merge(gdf, suffixes=('_e', '_g'), how='inner', left_on='RSNUM', right_on='RSNUM', left_index=False, right_index=False)
    merged.columns = ['GENE', 'i', 'tier', 'RSNUM', 'CHR', 'POS_HG19', 'A1_e', 'A2_e', 'INC_AFRQ_e', 'BETA_e', 'SE_e', 'LOGP', 'SNPID', 'CHR_g', 'POS_g', 'A1_g', 'A2_g', 'INC_ALLELE', 'INC_AFRQ_g', 'BETA_g', 'SE_g', 'PVAL']
    merged = merged[['GENE', 'RSNUM', 'CHR', 'POS_HG19', 'A1_e', 'A2_e', 'INC_AFRQ_e', 'BETA_e', 'SE_e', 'LOGP', 'SNPID', 'A1_g', 'A2_g', 'INC_ALLELE', 'INC_AFRQ_g', 'BETA_g', 'SE_g', 'PVAL']]
    merged['A1_g'] = merged['A1_g'].str.upper()
    merged['A2_g'] = merged['A2_g'].str.upper()
    merged['INC_ALLELE'] = merged['INC_ALLELE'].str.upper()
    
    mdf = merged[(merged.A1_e == merged.A1_g) & (merged.A2_e == merged.A2_g) & (merged.A1_e == merged.INC_ALLELE)]
    mdf = mdf.reset_index(drop=True)
    
    mdf0 = merged[(merged.A1_e == merged.A2_g) & (merged.A2_e == merged.A1_g) & (merged.A1_e == merged.INC_ALLELE)]
    mdf0 = mdf0.reset_index(drop=True)
    mdf0 = mdf0[['GENE', 'RSNUM', 'CHR', 'POS_HG19', 'A1_e', 'A2_e', 'INC_AFRQ_e', 'BETA_e', 'SE_e', 'LOGP', 'SNPID', 'A2_g', 'A1_g', 'INC_ALLELE', 'INC_AFRQ_g', 'BETA_g', 'SE_g', 'PVAL']]
    mdf0.columns = ['GENE', 'RSNUM', 'CHR', 'POS_HG19', 'A1_e', 'A2_e', 'INC_AFRQ_e', 'BETA_e', 'SE_e', 'LOGP', 'SNPID', 'A1_g', 'A2_g', 'INC_ALLELE', 'INC_AFRQ_g', 'BETA_g', 'SE_g', 'PVAL']
    
    # some summary statistics have minor/major coding, etc. rather than alt/ref
    # match alleles and compare to inc_allele (effective allele reported)
    # if effective allele in a1 or a2, flip beta and reorder columns accordingly
    mdf1 = merged[(merged.A1_e == merged.A1_g) & (merged.A2_e == merged.A2_g) & (merged.A2_e == merged.INC_ALLELE)]
    mdf1 = mdf1.dropna(subset=['BETA_g'])
    mdf1 = mdf1.reset_index(drop=True)
    mdf1['BETA_g_adj'] = mdf1['BETA_g'] * -1
    mdf1 = mdf1[['GENE', 'RSNUM', 'CHR', 'POS_HG19', 'A1_e', 'A2_e', 'INC_AFRQ_e', 'BETA_e', 'SE_e', 'LOGP', 'SNPID', 'A1_g', 'A2_g', 'A1_g', 'INC_AFRQ_g', 'BETA_g_adj', 'SE_g', 'PVAL']]
    mdf1.columns = ['GENE', 'RSNUM', 'CHR', 'POS_HG19', 'A1_e', 'A2_e', 'INC_AFRQ_e', 'BETA_e', 'SE_e', 'LOGP', 'SNPID', 'A1_g', 'A2_g', 'INC_ALLELE', 'INC_AFRQ_g', 'BETA_g', 'SE_g', 'PVAL']
    
    mdf2 = merged[(merged.A1_e == merged.A2_g) & (merged.A2_e == merged.A1_g) & (merged.A2_e == merged.INC_ALLELE)]
    mdf2 = mdf2.dropna(subset=['BETA_g'])
    mdf2 = mdf2.reset_index(drop=True)
    #mdf2 = mdf2[np.isfinite(mdf2['BETA_g'])]
    mdf2['BETA_g_adj'] = mdf2['BETA_g'] * -1
    mdf2 = mdf2[['GENE', 'RSNUM', 'CHR', 'POS_HG19', 'A1_e', 'A2_e', 'INC_AFRQ_e', 'BETA_e', 'SE_e', 'LOGP', 'SNPID', 'A2_g', 'A1_g', 'A2_g', 'INC_AFRQ_g', 'BETA_g_adj', 'SE_g', 'PVAL']]
    mdf2.columns = ['GENE', 'RSNUM', 'CHR', 'POS_HG19', 'A1_e', 'A2_e', 'INC_AFRQ_e', 'BETA_e', 'SE_e', 'LOGP', 'SNPID', 'A1_g', 'A2_g', 'INC_ALLELE', 'INC_AFRQ_g', 'BETA_g', 'SE_g', 'PVAL']
    
    temp0 = mdf.append(mdf0, ignore_index=True)
    temp1 = temp0.append(mdf1, ignore_index=True)
    temp2 = temp1.append(mdf2, ignore_index=True)
    outdf = temp2.reset_index(drop=True)
    del temp0, temp1, temp2, mdf, mdf0, mdf1, mdf2
    outdf.to_csv(outfile, sep='\t', index=False, na_rep='NA')
    
    return(outdf)


if __name__ == '__main__':
  
  # if main, run test with input files
  ewkdir = argv[1]
  efile = argv[2]
  gwkdir = argv[3]
  
  outdf = xinput(ewkdir, efile, gwkdir)
  
  
  