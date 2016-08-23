#!/usr/bin/env pytyhon3

###
# YoSon
# @treqtl/permutation.py
###

import pandas as pd
import numpy as np
from sys import argv

def shuffvars(_infile):
  """
  shuffle beta_eqtl from treqtl xinput file
  output a shuffled data
  """
  _outfile = _infile + '.shuffvars'
  _df = pd.read_csv(_infile, sep='\t', header=0)
  _betae = _df['BETA_e']
  _vars = len(_betae)
  _shuff = _betae.loc[np.random.permutation(len(_betae))]
  _shuffled = _shuff.reset_index(drop=True)
  _df['BETA_e'] = _shuffled
  return _df, _outfile

def _writedf(_df, _outfile):
  _df.to_csv(_outfile, sep='\t', index=False, na_rep='NA')
  

if __name__ == '__main__':
  _infile = argv[1]

  _df = shuffvars(_infile)
  _writedf(_df, _outfile)




