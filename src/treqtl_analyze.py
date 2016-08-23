#!/usr/bin/env python

###
# ypar
# @ treqtl/treqtl_analyze.py
# read eqtl and gwas summary stats files and analyze
###


import sys, os
import pandas as pd
import numpy as np
from sys import argv
from math import sqrt
from scipy import log10
from scipy.stats import chisqprob

import treqtl_input
import mr_estimate


def check_pal(ea1, ga1):
  # all input should be checked prior to analysis
  # function deleted 01/07/2016
  return True


def calculate_ivstat(eqtl, gwas, outfile, logfile):
  
  """calculate_ivstat
  
  parameters
  -----
  eqtl : pandas dataframe for iv summary statistics file
  gwas : pandas dataframe for gwas summary statistics file
  outfile : name of output file
  logfile : name of log file
  
  action
  -----
  runs treqtl analysis
  prints output
  """

  eqtldf, gwasdf, genes, ngenes = mr_estimate.printinit(eqtl, gwas, outfile)

  for gene in genes:
    # reset snp per gene count
    nsnps = 0
    # pick out relevant rows
    etemp = eqtldf[eqtldf['GENE'] == gene]
    # make a list of snps in gene
    #esnps = pd.unique(etemp.RSNUM.ravel())
    esnps = etemp.RSNUM.ravel()
    ngenes = mr_estimate.pergene(eqtldf, gwasdf, gene, ngenes, nsnps, etemp, esnps, outfile, logfile)

  # per input bonf pval
  bonfp = 1/(ngenes * 20)
  logmessage = '\n-------\na total of ' + str(ngenes) + ' variants are tested\n' + 'bonf nominal p = ' + str(bonfp)
  print(logmessage)


def calculate_ivstat_itreqtl(eqtl, gwas, outfile, logfile):
  
  """calculate_ivstat
  
  parameters
  -----
  eqtl : pandas dataframe for iv summary statistics file
  gwas : pandas dataframe for gwas summary statistics file
  outfile : name of output file
  logfile : name of log file
  
  action
  -----
  runs itreqtl analysis
  prints output
  """

  eqtldf, gwasdf, genes, ngenes = mr_estimate.printinit(eqtl, gwas, outfile)

  # for gene in genes:
  #   # reset snp per gene count
  #   nsnps = 0
  #   # pick out relevant rows
  #   etemp = eqtldf[eqtldf['GENE'] == gene]
  #   # make a list of snps in gene
  #   esnps = pd.unique(etemp.RSNUM.ravel())
  #   ngenes = mr_estimate.pergene(eqtldf, gwasdf, gene, ngenes, nsnps, etemp, esnps, outfile, logfile)

  nsnps = 0

  for index, row in eqtldf.iterrows():

    etemp = row
    esnps = etemp.RSNUM
    gene = etemp.GENE

    #print(index)
    #print(row)

    snp, erow, grow, gi, chrom, pos, rsnum, weight, iv_se, iv_alt, iv_afrq, iv_p, beta, se, endp_alt, endp_p, numer, denom, va, nsnps = mr_estimate.pervariant(eqtldf, gwasdf, etemp, nsnps, esnps, gene, outfile, logfile)
  
    ngenes = mr_estimate.printsum(snp, erow, grow, gi, chrom, pos, rsnum, weight, iv_se, iv_alt, iv_afrq, iv_p, beta, se, endp_alt, endp_p, numer, denom, va, nsnps, ngenes, outfile, logfile)

  # per input bonf pval
  bonfp = 1/(ngenes * 20)
  logmessage = '\n-------\na total of ' + str(ngenes) + ' variants are tested\n' + 'bonf nominal p = ' + str(bonfp)
  print(logmessage)


def run_treqtl(eqtldf, gwasdf, outfile, logfile, treqtl, itreqtl):
  
  """run_treqtl
  
  parameters
  -----
  eqtldf : pandas dataframe for iv summary statistics file
  gwasdf : pandas dataframe for gwas summary statistics file
  outfile : name of output file
  logfile : name of log file
  treqtl : boolean variable indicating whether to run treqtl analysis 
  itreqtl : boolean variable indicating whether to run itreqtl analysis
  
  action
  -----
  if treqtl = True, invokes calculate_ivstats function to run treqtl analysis
  """
  
  if isinstance(eqtldf, pd.DataFrame):
    logmessage = 'good to go! eqtl dataframe is loaded!'
    print(logmessage)
    
    if isinstance(gwasdf, pd.DataFrame):
      logmessage = 'good to go! gwas dataframe is loaded!'
      print(logmessage)
      
      if treqtl:
        # now analyze!
        calculate_ivstat(eqtldf, gwasdf, outfile, logfile)
      elif itreqtl:
        calculate_ivstat_itreqtl(eqtldf, gwasdf, outfile, logfile)
      else:
        #nothing to do
        pass
      
    else:
      logmessage = 'something is wrong with gwas data loading'
      print(logmessage)
      sys.exit()
      
  else:
    logmessage = 'something is wrong with eqtl data loading'
    print(logmessage)
    sys.exit()



if __name__ == '__main__':
  
  if len(sys.argv) < 2:
    logmessage = 'analysis terminated. input files are not specified'
    print(logmessage)
    
  else:
    # without main prompt, input files are taken as positional arguments...
    eqtlfile = argv[1]
    gwasfile = argv[2]
    
    eqtldf = treqtl_input.read_eqtl(eqtlfile)
    gwasdf = treqtl_input.read_gwas(gwasfile)
    
    calculate_ivstat(eqtldf, gwasdf, 'treqtl.out', 'treqtl.log')



