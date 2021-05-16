#!/usr/bin/env python3

import argparse
import pandas as pd
import functools
import time
import os
import sys

# Define timer decorator
def timer(func):
  '''Time the execution time of a function'''
  @functools.wraps(func)
  def wrapper(*args, **kwargs):
    start = time.time()
    res = func(*args, **kwargs)
    end = time.time()
    diff = end - start
    print("{} took {} seconds to complete".format(func.__name__, diff))
    return res
  return wrapper

def make_outdir(outdir):
  if not os.path.exists(outdir):
    os.makedirs(outdir)

@timer
def get_representatives(ST_to_HierCC):
  df = pd.read_csv(ST_to_HierCC, sep = '\t')
#  df = df.rename(columns = {'HC400 (cgST Cplx)': 'HC400'})
  HC400_representative = pd.DataFrame()
  HC400 = df['HC400'].unique()
  for HC in HC400:
    ST = df[df['HC400'] == HC].sort_values('ST').reset_index().loc[0, 'ST']
    HC400_representative = HC400_representative.append({'ST': ST, 'HC400': HC}, ignore_index=True)
  HC400_representative = HC400_representative.astype(int)
  return HC400_representative

@timer
def process_ST(df_out, ST, HC400, df_profiles):
  ST_row = df_profiles.query('ST == @ST')
  if ST_row.shape[0] > 0:
    ST_row.insert(0, 'HC400', HC400)
    df_out = pd.concat([df_out, ST_row])
  return df_out

@timer
def process_HC400(profiles_outdir, HC400, df_profiles, ST_to_HierCC):
  df = pd.read_csv(ST_to_HierCC, sep = '\t')
  df = df.rename(columns = {'HC400 (cgST Cplx)': 'HC400'})
  HC400_STs = df.query('HC400 == @HC400')['ST']
  if len(HC400_STs) > 0:
    list_to_concat = []
    for ST in HC400_STs:
      HC400_ST_rows = df_profiles.query('ST == @ST')
      list_to_concat.append(HC400_ST_rows)
    HC_rows = pd.concat(list_to_concat)
    filestring = profiles_outdir + '/HC400_' + str(HC400) + '.tsv'
    HC_rows.to_csv(filestring, sep = '\t', index=False)

@timer
def load_profiles(path_to_profiles):
  df_profiles = pd.read_csv(path_to_profiles, sep = '\t')
  return df_profiles

@timer
def main(args):
  make_outdir(args.profiles_outdir)
  HC400_representative = get_representatives(args.ST_to_HierCC)
  df_profiles = load_profiles(args.profiles)
  df_out = pd.DataFrame()
  for index, row in HC400_representative.iterrows():
    df_out = process_ST(df_out, row['ST'], row['HC400'], df_profiles)
    process_HC400(args.profiles_outdir, row['HC400'], df_profiles, args.ST_to_HierCC)
  df_out.to_csv(args.output, sep = '\t', index=False)

if __name__ == '__main__':
  # Parse arguments
  parser = argparse.ArgumentParser(description='Link chewBBACA typing results to Enterobase cgMLST and corresponding hierCCs.')

  parser.add_argument('-o', '--output', dest='profiles_outdir', help="Directory where profiles are stored", type=str, required=True)
  parser.add_argument('-i', '--input', dest='input', help="Input typing from chewBBACA", type=str, required=True)
  parser.add_argument('-s', '--st-to-hiercc', dest='st_to_hiercc', help="Table relating ST to hierCC levels (default: ST_to_HierCC.tsv)", type=str, default="ST_to_HierCC.tsv")
  parser.add_argument('-r', '--representatives', dest='output', help="Output table with HC400 representatives", type=str, required=True)
  parser.add_argument('-p', '--profiles', dest='profiles', help="Table listing cgMLSTs with profiles (default: profiles.list)", type=str, default="profiles.list")

  args = parser.parse_args()

  main(args)
