#!/usr/bin/env python3

import argparse
import pandas as pd
import time
import functools
import os
import sys

# Suppress warning raised in process_chunk from https://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas
pd.options.mode.chained_assignment = None

# Define timer decorator
def timer(func):
  '''
  Time the execution time of a function.
  '''
  @functools.wraps(func)
  def wrapper(*args, **kwargs):
    start = time.time()
    res = func(*args, **kwargs)
    end = time.time()
    diff = end - start
    print("{} took {} seconds to complete".format(func.__name__, diff))
    return res
  return wrapper

@timer
def load_results(HC1100_representative_file, chewbbaca_output_file):
  '''
  Load results of cgMLST typing with correct header order.

  Parameters
  ----------
  profiles_file : str
    Path to profiles.list file from Enterobase.
  chewbbaca_output_file : str
    Path to chewBBACA results_alleles.tsv file.

  Returns
  -------
  chewbbaca_output : pandas.DataFrame
    DataFrame containing the loaded chewBBACA output, with columns ordered the same as in profiles_file.
  nr_loci : int
    Number of loci present in profiles_file.

  '''
  # Extract header, loci, number of loci
  HC1100_representative = pd.read_csv(HC1100_representative_file, sep = '\t', dtype=str)
  loci_header = HC1100_representative.columns[2:]
  nr_loci = len(HC1100_representative.columns[2:])
  # Init header for chewbbaca output
  chewbbaca_output_header = ['FILE']
  # Append .fasta to loci names to correspond to chewBBACA output
  for locus in loci_header:
    chewbbaca_output_header.append(locus + '.fasta')
  # Read chewBBACA output file with correct column names and order
  chewbbaca_output = pd.read_csv(chewbbaca_output_file, sep = '\t', usecols=chewbbaca_output_header, dtype=str)[chewbbaca_output_header]
  chewbbaca_output.columns = ['FILE'] + list(loci_header)
  HC1100_representative.drop(columns = 'ST', inplace = True)
  return HC1100_representative, chewbbaca_output, nr_loci

@timer
def process_chunk(chunk, level, chewbbaca_output_isolate, i):
  '''
  Compare a chunk against chewBBACA output and return cgMLSTs.

  Return cgMLSTs together with a count of matching alleles, which can be used to select the most similar cgMLST.

  Parameters
  ------
  chunk : pandas.DataFrame
    Generated chunk of pandas.DataFrame containing cgMLSTs and corresponding allelic profiles.
  chewbbaca_output_isolate : pandas.DataFrame
    DataFrame consisting of a single row, where columns represent allelic profiles.
  i : int
    Used to keep track of the for loop.

  Returns
  ------
  df_ST : pandas DataFrame
    DataFrame containing two columns: cgMLST, number of alleles matching that cgMLST
  nr_rows_chunk: int
    Number of rows of the chunk.

  '''
  # Read STs and shape from chunk
  df_ST = chunk[[level]]
  nr_rows_chunk, nr_cols_chunk = chunk.shape
  # Recreate df with same size as chunk, filled with isolate allelic profile
  chewbbaca_output_comparator = pd.DataFrame(index=range(nr_rows_chunk))
  for locus in chewbbaca_output_isolate.columns:
    chewbbaca_output_comparator[locus] = chewbbaca_output_isolate.at[i,locus]
  chewbbaca_output_comparator.index = chunk.index
  # Compare dfs and sum the number of matching cells
  comparison = chewbbaca_output_comparator.iloc[:,1:] == chunk.iloc[:,1:]
  df_ST['matches'] = comparison.sum(axis=1)
  return df_ST

@timer
def identify_cgMLST_isolate(HC1100_representative, profiles_outdir, chewbbaca_output, i):
  '''
  Find cgMLST that most closely matches the allelic profile from chewBBACA output.

  Process the profiles from profiles_file in chunks and call process_chunks() on every chunk. Subsequently select the highest matching cgMLST.

  Parameters
  ----------
  profiles_file : str
    Path to file containing profiles.list from Enterobase.
  chewbbaca_output : pandas.DataFrame
    chewBBACA output loaded by load_results function..
  i : int
    Used to keep track of the for loop.

  Returns
  -------
  ST : str
    cgMLST that shows most matches with allelic profile in chewbbaca output.
  matching_alleles : int
    Number of alleles from chewbbaca_output that match the printed cgMLST.
  isolate_name : str
    Isolate name.

  '''
  chewbbaca_output_isolate = chewbbaca_output.iloc[[i]]
  isolate_name = chewbbaca_output_isolate.iloc[0,0].rstrip('.fasta')
  df_HC1100 = process_chunk(HC1100_representative, 'HC1100', chewbbaca_output_isolate, i)
  df_HC1100 = df_HC1100.sort_values('matches', ascending=False)
  selected_HC1100 = df_HC1100.iloc[0,:]['HC1100']

  filepath = profiles_outdir + '/HC1100_' + str(selected_HC1100) + '.tsv'
  chunk = pd.read_csv(filepath, sep = '\t', dtype=str)
  df_ST = process_chunk(chunk, 'ST', chewbbaca_output_isolate, i)
  df_ST = df_ST.sort_values('matches', ascending=False)
  selected_ST = df_ST.iloc[0,:]['ST']
  matching_alleles = int(df_ST.iloc[0,:]['matches'])
  return selected_ST, matching_alleles, isolate_name

def select_hierCC(df_ST_to_hierCC, ST):
  '''
  Find corresponding hierarchical CC from cgMLST.

  Looks up the provided cgMLST in a table relating cgMLST to hierarchical clonal complexes.

  Parameters
  ----------
  df_ST_to_hierCC : pandas.DataFrame
    Lookup table containing cgMLST and corresponding hierarchical CCs.
  ST : str
    cgMLST to query for.

  Returns
  -------
  pandas.DataFrame
    DataFrame containing a single row containing the provided cgMLST and corresponding hierarchical CCs.

  Raises
  ------
  ExceptionError
    If the supplied cgMLST cannot be found in the ST to hierCC lookup table (df_ST_to_hierCC).

  Examples
  --------
  >>> import pandas as pd
  >>> df = pd.DataFrame([['10', '10', '10'], ['131', '131', '10'], ['38', '131', '10']], index=[0, 1, 2], columns=['ST', 'HC10', 'HC200'])
  >>> ST = '38'
  >>> select_hierCC(df, ST)
     ST HC10 HC200
  2  38  131    10
  '''
  # Query input dataframe
  df_selected_ST_hierCC = df_ST_to_hierCC.query('ST == @ST')
  # Check whether resulting dataframe contains zero rows
  if df_selected_ST_hierCC.shape[0] == 0:
    raise Exception('cgMLST {} was not found in the ST to hierCC lookup table. Please check whether these files were produced from the same database version'.format(ST))
  else:
    return df_selected_ST_hierCC

def collect_output_data(selected_df, isolate_name, matching_alleles, nr_loci):
  '''
  Collect all output data and return pandas.DataFrame.

  Collects input data (isolate name, number of matching alleles, total number of loci), computes one extra column (maximum number of mismatches) and returns edited DataFrame.

  Parameters
  ------
  selected_df : pandas.DataFrame
    DataFrame containing cgMLST and corresponding hierCCs.
  isolate_name : str
    Isolate name.
  matching_alleles : int
    Number of loci where the alleles match between cgMLST and chewBBACA results.
  nr_loci : int
    Number of loci in cgMLST typing scheme.

  Returns
  -------
  pandas.DataFrame
    Same as input DataFrame, but with extra data added.

  Raises
  ------
  ExceptionError
    If fewer or more than four columns are added to the supplied DataFrame.

  '''
  hierCC_threshold = nr_loci - matching_alleles
  levels_dict = {'HC20': 20, 'HC50': 50, 'HC100': 100, 'HC200': 200, 'HC400': 400, 'HC1100': 1100, 'unreliable': 9999}
  for confidence_level, threshold in levels_dict.items():
    if hierCC_threshold <= threshold:
      break
  nr_columns_original = len(selected_df.columns)
  selected_df.insert(0, 'isolate_name', [isolate_name])
  selected_df.insert(1, 'matching_alleles', [matching_alleles])
  selected_df.insert(2, 'max_mismatches', [hierCC_threshold])
  selected_df.insert(3, 'confidence_level', [confidence_level])
  selected_df.insert(4, 'nr_loci', [nr_loci])
  nr_columns_final = len(selected_df.columns)
  nr_columns_difference = nr_columns_final - nr_columns_original
  if nr_columns_difference == 5:
    return selected_df
  else:
    raise Exception('Expected 5 columns would be inserted, but {} columns were inserted'.format(nr_columns_differences))

@timer
def main(args):
  '''
  Matches chewBBACA allele calls to Enterobase cgMLST scheme

  The file results_alleles.tsv can be provided to this script and will be matched to Enterobase profiles (linking cgMLST to alleles) and subsequently to hierarchical clonal complexes (hierCC).

  Parameters
  ------
  -i, --input: results_alleles.tsv file from chewBBACA
  -p, --profiles: table with in the first column a cgMLST and in the other columns the alleles the correspond to that cgMLST
  -s, --st-to-hiercc: table describing to which hierarchical clonal complexes cgMLSTs are assigned

  Returns
  ------
  -o, --output: table containing cgMLST, matching alleles and various HCs for an isolate

  Notes
  ------
  chewBBACA excludes a small percent of the Enterobase cgMLST alleles. This results that chewBBACA profiles will never fully match the Enterobase cgMLST scheme, but can come pretty close.
  '''
  # Load chewbbaca output and table linking ST to hierCC
  HC1100_representative, chewbbaca_output, nr_loci = load_results(args.HC1100_representative, args.input)
  df_ST_to_hierCC = pd.read_csv(args.st_to_hiercc, sep = '\t', dtype=str)
  # Init df which will be written to output
  output_df = pd.DataFrame()
  # Loop over chewbbaca output
  for i in range(chewbbaca_output.shape[0]):
    # identify most closely matching cgMLST
    selected_ST, matching_alleles, isolate_name = identify_cgMLST_isolate(HC1100_representative, args.profiles_outdir, chewbbaca_output, i)
    # Select corresponding hierCCs
    selected_df = select_hierCC(df_ST_to_hierCC, selected_ST)
    # Collect output data
    selected_df_added = collect_output_data(selected_df, isolate_name, matching_alleles, nr_loci)
    # Concat results to output df
    output_df = pd.concat([output_df, selected_df_added])

  # Save output to csv file
  output_df.to_csv(args.output, index=False)

if __name__ == "__main__":
  # Parse arguments
  parser = argparse.ArgumentParser(description='Link chewBBACA typing results to Enterobase cgMLST and corresponding hierCCs.')

  parser.add_argument('-p', '--profiles', dest='profiles_outdir', help="Directory where profiles are stored", type=str, required=True)
  parser.add_argument('-i', '--input', dest='input', help="Input typing from chewBBACA", type=str, required=True)
  parser.add_argument('-o', '--output', dest='output', help="Output file", type=str, required=True)
  parser.add_argument('-s', '--st-to-hiercc', dest='st_to_hiercc', help="Table relating ST to hierCC levels", type=str, required=True)
  parser.add_argument('-r', '--representatives', dest='HC1100_representative', help="Table with HC1100 representatives", type=str, required=True)

  args = parser.parse_args()

  main(args)
