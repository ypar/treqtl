#!/usr/bin/env python3

###
# ypar
# @ treqtl/xtreqtl.py
# take combined input of iv + trait
###


import sys, os
import pandas as pd
import numpy as np
from os import walk
from sys import argv
from math import sqrt
from scipy import log10
from scipy.stats import chisqprob

import treqtl_input
import mr_estimate
import treqtl_analyze


# header xinput
# mdf2.columns = ['GENE', 'RSNUM', 'CHR', 'POS_HG19', 'A1_e', 'A2_e', 'INC_AFRQ_e', 'BETA_e', 'SE_e', 'LOGP', 'SNPID', 'A1_g', 'A2_g', 'INC_ALLELE', 'INC_AFRQ_g', 'BETA_g', 'SE_g', 'PVAL']


def read_xfile(infile):
  
  """read_xfile
  
  parameters
  -----
  infile : input summary statistics file in the following format (tab or space delimited)
  # column headings
  # GENE RSNUM CHR POS_HG19 A1_e A2_e INC_AFRQ_e BETA_e SE_e LOGP SNPID A1_g A2_g INC_ALLELE INC_AFRQ_g BETA_g SE_g PVAL
  
  action
  -----
  reads contents of the file into a pandas dataframe
  
  """
  df = pd.read_csv(infile, delim_whitespace=True, header=0, na_values='NA')
  logmessage = str(len(df.RSNUM)) + ' variants are found in file ' + infile
  print(logmessage)
  return(df)
  


def calculate_xstat(xfile, outfile, logfile):
  
  """calculate_ivstat
  
  parameters
  -----
  xfile : pandas dataframe for iv and trait summary statistics file
  outfile : name of output file
  logfile : name of log file
  
  action
  -----
  runs xtreqtl analysis
  prints output
  """

  xdf = read_xfile(xfile)
  
  if isinstance(xdf, pd.DataFrame):
    logmessage = 'good to go! eqtl x gwas dataframe is loaded!'
    print(logmessage)
      
  else:
    logmessage = 'something is wrong with eqtl x gwas data loading'
    print(logmessage)
    sys.exit()
  
  nsnps = 0
  ngenes = 0
  
  if not xdf.empty:
    outdf = pervariant(xdf, outfile, logfile)
    nsnps = len(outdf.RSNUM)
    #ngenes = mr_estimate.pergene(eqtldf, gwasdf, gene, ngenes, nsnps, etemp, esnps, outfile, logfile)
    # per input bonf pval
    bonfp = 1/(nsnps * 20)
    logmessage = '\n-------\na total of ' + str(nsnps) + ' variants are tested\n' + 'bonf nominal p = ' + str(bonfp) + '\n-------\n'
    print(logmessage)
  else:
    logmessage = 'given input file is empty\n-------\n'
    print(logmessage)
    sys.exit()


def pervariant(xdf, outfile, logfile):
  """pervariant

  parameters
  -----
  xdf : iv and trait summary statistics in pandas dataframe
  outfile : output file name
  logfile : log file name 

  action
  -----
  runs analysis for each variant given 

  """
  
  temp = xdf
  
  # initialize
  nsnps = 0
  ngenes = 0
  temp['va'] = 0
  temp['numer'] = 0
  temp['denom'] = 0
  
  # select correspoding rows from two input files
  if isinstance(xdf, pd.DataFrame):
    pass

  else:
    logmessage = 'check the input'
    print(logmessage)
    sys.exit()

  # header xinput
  # mdf2.columns = ['GENE', 'RSNUM', 'CHR', 'POS_HG19', 'A1_e', 'A2_e', 'INC_AFRQ_e', 'BETA_e', 'SE_e', 'LOGP', 'SNPID', 'A1_g', 'A2_g', 'INC_ALLELE', 'INC_AFRQ_g', 'BETA_g', 'SE_g', 'PVAL']

  try:
    temp['weight'] = abs(xdf['BETA_e'].convert_objects(convert_numeric=True))
  except ValueError:
    logmessage = 'check for value. possible duplicates causing type conversion error'
    print(logmessage)
  
  temp['iv_se'] = xdf['SE_e'].convert_objects(convert_numeric=True)
  temp['iv_alt'] = str(xdf['A1_e'])
  temp['iv_afrq'] = xdf['INC_AFRQ_e'].convert_objects(convert_numeric=True)
  
  if 'PVAL_e' in xdf.columns:
    temp['iv_p'] = xdf['PVAL_e'].convert_objects(convert_numeric=True)
  elif 'LOGP' in xdf.columns:
    temp['iv_p'] = xdf['LOGP'].convert_objects(convert_numeric=True)
  else:
    temp['iv_p'] = xdf.iloc[:,len(xdf.columns)-1]
  
  # pick out variables from gwas
  temp['beta'] = xdf['BETA_g'].convert_objects(convert_numeric=True)
  temp['se'] = xdf['SE_g'].convert_objects(convert_numeric=True)
  temp['endp_alt'] = str(xdf['A1_g'])
  
  if 'PVAL_g' in xdf.columns:
    temp['endp_p'] = xdf['PVAL_g'].convert_objects(convert_numeric=True)
  else:
    temp['endp_p'] = xdf['PVAL'].convert_objects(convert_numeric=True)
  
  temp['numer'] = temp['numer'].convert_objects(convert_numeric=True)
  temp['denom'] = temp['denom'].convert_objects(convert_numeric=True)
  
  temp['numer'] = temp['numer'] + temp['weight'] * temp['beta'] * temp['se']**-2
  temp['denom'] = temp['denom'] + temp['weight']**2 * temp['se']**-2
  
  #print(xdf.numer)
  #print(xdf.denom)
  
  temp['va'] = temp['va'] + 2 * temp['iv_afrq'] * (1 - temp['iv_afrq']) * temp['weight']**2
  
  # done with snp
  nsnps = len(xdf.RSNUM)

  # finish calculations
  temp['ahat'] = temp['numer'] / temp['denom']
  temp['overdenom'] = 1 / temp['denom']
  temp['se_ahat'] = temp['overdenom'].apply(np.sqrt)
  #temp['se_ahat'] = sqrt((1/temp['denom']))
  temp['chisqstat'] = (temp['ahat']/temp['se_ahat'])**2
  temp['chisqstat'] = temp['chisqstat'].convert_objects(convert_numeric=True)
  temp['pval'] = chisqprob(temp['chisqstat'], 1)
  ngenes = ngenes + 1
  
  outdf = temp[['GENE', 'RSNUM', 'CHR', 'POS_HG19', 'A1_e', 'A2_e', 'INC_AFRQ_e', 'ahat', 'se_ahat', 'chisqstat', 'pval', 'weight', 'iv_se', 'iv_p', 'BETA_g', 'SE_g', 'endp_p']]
  
  outdf.to_csv(outfile, sep='\t', na_rep='NA', index=False)
  # write output
  #'GENE\tRSNUM\tCHR\tPOS\tAHAT\tSE\tCHISQ\tPVAL\tIV_VA\tIV_NSNPS\tIV_BETA\tIV_SE\tIV_P\tENDP_BETA\tENDP_SE\tENDP_P\n'
  #outrow = str(gi) + '\t' + str(rsnum) + '\t' + str(chrom) + '\t' + str(pos) + '\t' + str(ahat) + '\t' + str(se_ahat) + '\t' + str(chisqstat) + '\t' + str(pval) + '\t' + str(va) + '\t' + str(nsnps) + '\t' + str(weight) + '\t' + str(iv_se) + '\t' + str(iv_p) + '\t' + str(beta) + '\t' + str(se) + '\t' + str(endp_p) + '\n'

  #out = open(outfile, 'a')
  #out.write(outrow)
  #out.close()

  return(outdf)


def iter_dir(inputdir):
  
  """iter_dir
  
  parameters
  -----
  dirname : input file dir
  treqtl : boolean variable indicating whether to run treqtl analysis 
  itreqtl : boolean variable indicating whether to run itreqtl analysis
  
  action
  -----
  initializes a parameter setup for each pairwise combinations of efiles and gfiles
  invokes run_treqtl function
  """
  
  os.chdir(inputdir)
  files = []
  
  for (dirpath, dirnames, filenames) in walk(inputdir):
    files.extend(filenames)
    files = [i if i.startswith('/') else dirpath + '/' + i for i in files]
    break
  
  if type(files) is list:
    i = 0
    for xfile in files:
      outfile = str(xfile).replace('xinput', 'xoutput')
      outfile = str(outfile).replace('.txt', '_xtreqtl.txt')
      logfile = str(outfile).replace('.txt', '.log')
      print('python3 xtreqtl.py ', xfile, outfile, logfile)
    return(xfile, outfile, logfile)
  else:
    pass
  pass


if __name__ == '__main__':
  
  if len(argv) == 2:
    inputdir = argv[1]
    iter_dir(inputdir)
  else:
    pass
  
  if len(argv) == 4:
    #example = '/project/chrbrolab/gtex/mreqtl/xinput/Adipose_Subcutaneous_w_CARDIoGRAM_CHD_results.txt'
    #outfile = '/project/chrbrolab/gtex/mreqtl/xtreqtl.out'
    #logfile = '/project/chrbrolab/gtex/mreqtl/xtreqtl.log'
    xfile = argv[1]
    outfile = argv[2]
    logfile = argv[3]
    
    calculate_xstat(xfile, outfile, logfile)
  else:
    pass


