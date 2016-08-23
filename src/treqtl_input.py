#!/usr/bin/env python

###
# ypar
# @ treqtl/treqtl_input.py
# read input files
###

import sys, os, types
import pandas as pd
import numpy as np
from os import walk
from sys import argv
from itertools import product, chain, combinations

import treqtl_analyze
import count_vars


def check_inputs(infile):
  
  """check_inputs
  
  parameters
  -----
  infile: input file name
  
  action
  -----
  check if the file can be opened. return boolean
  """
  
  try:
    f = open(infile, 'r')
  
  except IOError:
    logmessage = 'check ' + infile + '. the file could not be opened.'
    print(logmessage)
    sys.exit()
    
  else:
    f.close()
    return True
  

def read_eqtl(infile):
  
  """read_eqtl
  
  parameters
  -----
  infile : input iv summary statistics file in the following format (tab or space delimited)
  # column headings
  # GENE	RSNUM	CHR	POS_HG19	A1	A2  INC_ALLELE  INC_AFRQ	BETA	SE	LOGP
  
  action
  -----
  reads contents of the file into pandas dataframe
  
  """
  if check_inputs(infile):
    df = pd.read_csv(infile, delim_whitespace=True, header=0, na_values='NA')
    logmessage = str(len(df)) + ' eqtl variants are found in file ' + infile
    print(logmessage)
  
    return(df)
  
  else:
    sys.exit()


def read_gwas(infile):
  
  """read_gwas
  
  parameters
  -----
  infile : input trait summary statistics file in the following format (tab or space delimited)
  # column headings
  # RSNUM SNPID CHR POS_HG18|19 A1 A2 INC_ALLELE INC_AFRQ BETA SE PVAL
  
  action
  -----
  reads contents of the file into pandas dataframe
  
  """
  if check_inputs(infile):
    df = pd.read_csv(infile, delim_whitespace=True, header=0, na_values='NA', low_memory=False)
    logmessage = str(len(df)) + ' gwas variants are found in file ' + infile
    print(logmessage)
    
    # no longer supports hg18...
    # checks deleted 01/07/2016
    # will fix this to be able to distinguish hg18/19/20 in the future though
    #if df.columns[3] == 'POS_HG19':
    #  pass
    #elif df.columns[3] == 'POS_HG18':
    #  treqtl_mapping.find_snp(df)
    #else:
    #  sys.exit()
    
    return(df)
    
  else:
    sys.exit()
    

def read_dir(inputdir):
  
  """read_dir
  
  parameters
  -----
  inputdir : path to directory with input files (string)
  
  action
  -----
  given an input directory path,
  iterates and returns a list of file names in the directory
  """
  
  os.chdir(inputdir)
  files = []
  
  for (dirpath, dirnames, filenames) in walk(inputdir):
    files.extend(filenames)
    files = [i if i.startswith('/') else dirpath + '/' + i for i in files]
    break
  
  return files


def iter_dir(efiles, gfiles, treqtl, itreqtl):
  
  """iter_dir
  
  parameters
  -----
  efiles : list containing names of iv files (e.g. eqtls)
  gfiles : list containing names of trait files (e.g. gwas)
  treqtl : boolean variable indicating whether to run treqtl analysis 
  itreqtl : boolean variable indicating whether to run itreqtl analysis
  
  action
  -----
  initializes a parameter setup for each pairwise combinations of efiles and gfiles
  invokes run_treqtl function
  """
  
  efiles = list(efiles)
  gfiles = list(gfiles)
  if type(efiles) is list and type(gfiles) is list:
    i = 0
    for r in chain(product(efiles, gfiles)):
      efile = r[0]
      gfile = r[1]
      outfile = 'treqtl_results_' + str(i) + '.txt'
      logfile = 'treqtl_results_' + str(i) + '.log'
      edf = read_eqtl(efile)
      gdf = read_gwas(gfile)
      treqtl_analyze.run_treqtl(edf, gdf, outfile, logfile, treqtl, itreqtl)
      i = i + 1
    pass
  else:
    pass
  pass



if __name__ == '__main__':
  
  # if main, run test with input files
  eqtldf = read_eqtl(argv[1])
  gwasdf = read_gwas(argv[2])




