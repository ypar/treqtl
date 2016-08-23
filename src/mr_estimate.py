#!/usr/bin/env python3

###
# ypar
# @ treqtl/mr_estimate.py
# core calculations
###


import sys, os
import pandas as pd
import numpy as np
from sys import argv
from math import sqrt
from scipy import log10
from scipy.stats import chisqprob

import treqtl_input
import treqtl_analyze
from run_mp import *


def prep_ivstat(eqtl, gwas):
  
  """prep_ivstat
  
  parameters
  -----
  eqtl : pandas dataframe for iv summary statistics file
  gwas : pandas dataframe for trait summary statistics file
  
  action
  -----
  checks and matches rs numbers of variants in both input data
  returns coinciding subset
  """
  
  # check gwas
  gwastemp = pd.DataFrame()
      
  # filter vars to include ones with both eqtl and gwas results
  snps = pd.unique(eqtl.RSNUM.ravel())
      
  for snp in snps:
    if str(snp) != '.':
      temp = gwas[gwas['RSNUM'] == snp]
        
      if not temp.empty:
        gwastemp = gwastemp.append(temp, ignore_index=True)
          
  
  gwasvar = len(gwas)
  gwassel = len(gwastemp)
  
  logmessage = 'of ' + str(gwasvar) + ' variants in the gwas input,\n' + str(gwassel) + ' variants are mapped to the eqtl input and will remain for analysis'
  print(logmessage)
  
  # check eqtl
  eqtltemp = pd.DataFrame()
      
  # filter vars to include ones with both eqtl and gwas results
  snps = pd.unique(gwastemp.RSNUM.ravel())
      
  for snp in snps:
    if str(snp) != '.':
      temp = eqtl[eqtl['RSNUM'] == snp]
        
      if not temp.empty:
        eqtltemp = eqtltemp.append(temp, ignore_index=True)
          
  
  eqtlvar = len(eqtl)
  eqtlsel = len(eqtltemp)
  
  logmessage = 'of ' + str(eqtlvar) + ' variants in the eqtl input,\n' + str(eqtlsel) + ' variants are mapped to the gwas input and will remain for analysis'
  print(logmessage)
  
  return(eqtltemp, gwastemp) 


def printinit(eqtl, gwas, outfile):

  # note: functions related to verbose log deleted 01/07/2016

  # initialize outfile the lazy way
  colnames = 'GENE\tRSNUM\tCHR\tPOS\tAHAT\tSE\tCHISQ\tPVAL\tIV_VA\tIV_NSNPS\tIV_BETA\tIV_SE\tIV_P\tENDP_BETA\tENDP_SE\tENDP_P\n'
  out = open(outfile, 'w')
  out.write(colnames)
  out.close()
  
  # filter input files for coinciding vars
  eqtldf, gwasdf = prep_ivstat(eqtl, gwas)
  
  #genes = pd.unique(eqtldf.GENE.ravel())
  genes = eqtldf.GENE.ravel()

  # count genes in input
  ngenes = 0

  return(eqtldf, gwasdf, genes, ngenes)


def pervariant(eqtldf, gwasdf, snp, etemp, nsnps, esnps, gene, outfile, logfile):
  """pervariant

  parameters
  -----
  eqtldf : iv summary statistics in pandas dataframe
  gwasdf : trait summary statistics in pandas dataframe
  etemp : eqtldf with variants in gwasdf
  nsnps : number of variants per transcript
  esnps : variant in etemp
  gene : gene/transcript unit annotated to variant in analysis
  outfile : output file name
  logfile : log file name 

  action
  -----
  runs analysis for each variant given 

  """

  # initialize
  if nsnps == 0:
    va = 0
    numer = 0
    denom = 0
  
  # select correspoding rows from two input files
  if isinstance(etemp, pd.DataFrame):
    erow = etemp[etemp['RSNUM'] == snp]
    # pick out variables from eqtl
    gi = erow['GENE'].values[0]
    chrom = erow['CHR'].values[0]
    pos = erow['POS_HG19'].values[0]
    rsnum = erow['RSNUM'].values[0]

  else:
    erow = etemp
    gi = gene
    chrom = erow['CHR']
    pos = erow['POS_HG19']
    rsnum = erow['RSNUM']

  grow = gwasdf[gwasdf['RSNUM'] == snp]

  try:
    weight = abs(float(erow['BETA']))
  except ValueError:
    logmessage = 'check for value. possible duplicates causing type conversion error for lines in\n' + erow['BETA']
    print(logmessage)
  
  iv_se = float(erow['SE'])
  iv_alt = str(erow['A1'])
  iv_afrq = float(erow['INC_AFRQ'])
  
  if 'PVAL' in erow:
    iv_p = float( erow['PVAL'])
  elif 'LOGP' in erow:
    iv_p = float( erow['LOGP'])
  else:
    iv_p = erow.iloc[:,len(erow.columns)-1]
  
  # check using snpid has been deprecated 01/07/2016
  #snpid = erow['SNPID']
  
  # pick out variables from gwas
  beta = float(grow['BETA'])
  se = float(grow['SE'])
  endp_alt = str(grow['A1'])
  endp_p = float(grow['PVAL'])
  #endp_snpid = grow['SNPID']
  
  # check if erow.A1 and grow.A1 match
  if treqtl_analyze.check_pal(iv_alt, endp_alt):
    pass
  else:
    # note: add checking and matching procedure back into this
    logmessage = iv_alt + ' and ' + endp_alt + ' do not match. check strands'
    pass
  
  numer = numer + weight * beta * se**-2
  denom = denom + weight**2 * se**-2
  
  va = va + 2 * iv_afrq * (1 - iv_afrq) * weight**2
  
  # done with snp
  nsnps = nsnps + 1

  return(snp, erow, grow, gi, chrom, pos, rsnum, weight, iv_se, iv_alt, iv_afrq, iv_p, beta, se, endp_alt, endp_p, numer, denom, va, nsnps)


def pergene(eqtldf, gwasdf, gene, ngenes, nsnps, etemp, esnps, outfile, logfile):
  """pergene

  action
  ------
  analyze independent variants per transcript

  """
  
  # calculate stats for each snp in gene
  for snp in esnps:
    
    snp, erow, grow, gi, chrom, pos, rsnum, weight, iv_se, iv_alt, iv_afrq, iv_p, beta, se, endp_alt, endp_p, numer, denom, va, nsnps = pervariant(eqtldf, gwasdf, snp, etemp, nsnps, esnps, gene, outfile, logfile)
    # done with all snps in gene (if best snp per gene, these are equivalent)
    # note: add ld back for penalized gene-based sum

  ngenes = printsum(snp, erow, grow, gi, chrom, pos, rsnum, weight, iv_se, iv_alt, iv_afrq, iv_p, beta, se, endp_alt, endp_p, numer, denom, va, nsnps, ngenes, outfile, logfile)
  
  return(ngenes)
  

def printsum(snp, erow, grow, gi, chrom, pos, rsnum, weight, iv_se, iv_alt, iv_afrq, iv_p, beta, se, endp_alt, endp_p, numer, denom, va, nsnps, ngenes, outfile, logfile):
  """printsum

  action
  -----
  followup on pervariant - redefined from class object to test multiproc.
  equivalent for treqtl and itreqtl
  """

  # finish calculations
  if int(nsnps) != 0:
    ahat = numer / denom
    se_ahat = sqrt((1/denom))
    chisqstat = (ahat/se_ahat)**2
    pval = chisqprob(chisqstat, 1)
    ngenes = ngenes + 1
    
  else:
    ahat = 'NA'
    se_ahat = 'NA'
    chisqstat = 'NA'
    pval = 'NA'
    
  # write output
  outrow = str(gi) + '\t' + str(rsnum) + '\t' + str(chrom) + '\t' + str(pos) + '\t' + str(ahat) + '\t' + str(se_ahat) + '\t' + str(chisqstat) + '\t' + str(pval) + '\t' + str(va) + '\t' + str(nsnps) + '\t' + str(weight) + '\t' + str(iv_se) + '\t' + str(iv_p) + '\t' + str(beta) + '\t' + str(se) + '\t' + str(endp_p) + '\n'

  out = open(outfile, 'a')
  out.write(outrow)
  out.close()

  return(ngenes)


if __name__ == '__main__':
  
  print('check document')



