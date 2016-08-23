#!/usr/bin/env python

###
# ypar
# @ treqtl/count_vars.py
# check number of variants in files
###

import os


def count_lines(files):

  lines = []
  totalfi = 0
  totalfl = 0
  
  for fi in files:
    fl = sum(1 for line in open(fi))
    totalfl = totalfl + fl
    lines.append(fl)
    totalfi = totalfi + 1
    print(str(fi) + ' includes ' + str(fl) + ' variants')
    
  print('-----')
  print('total ' + str(totalfi) + ' files are in the directory')
  print('total ' + str(totalfl) + ' variants are included in input files')
  
  return(lines)


def count_vars(files):
  
  """count_vars
  
  parameters
  ------
  files : list of files
  
  action
  -----
  iterates through files in the list and count number of lines (variants)
  """
  
  if type(files) is list:
    files = list(files)
    count_lines(files)

  elif isinstance(files, str) and os.path.isfile(filelist):
    filelist = [line.strip() for line in open(files, 'r')]
    count_lines(filelist)

  else:
    print('input argument must be a list object or a filename')
    sys.exit()


if __name__ == '__main__':
  
  # if main, run test with input files
  count_vars(argv[1])

