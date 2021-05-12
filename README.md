# Using chewBBACA typing output to find corresponding Enterobase cgMLST and hierCC

https://github.com/boasvdp/match_chewBBACA_to_Enterobase/actions/workflows/match_chewBBACA_to_Enterobase.yml/badge.svg

## Introduction

[ChewBBACA](https://github.com/B-UMMI/chewBBACA) is a "comprehensive pipeline including a set of functions for the creation and validation of whole genome and core genome MultiLocus Sequence Typing (wg/cgMLST) schemas". The Enterobase cgMLST + hierCC scheme can be downloaded through the PubMLST API and used in chewBBACA. However, to match the chewBBACA output to the Enterobase typing scheme, some extra scripts are needed which are provided (experimentally) here.

## Methods

The main script loads chewBBACA typing output (which may contain allelic profiles of multiple isolates) and matches these against the Enterobase cgMLST allelic profiles. Because there are >160.000 cgMLSTs with allelic profiles, the profiles file is read in chunks (fixed at 10,000 lines now). The scripts takes approximately 4 per isolate.

## Results

The script writes its output to a csv file, containing the columns:

- Isolate name
- Number of matching alleles, number of loci in cgMLST typing scheme and the difference between those, which represents the maximum number of mismatches
- cgMLST
- Various HierCC levels (HC0, HC2, HC5, HC10, HC20, HC50, HC100, HC200, HC400, HC1100, HC1500, HC2000, HC2350)

## What I learned or plan to learn

- [x] Writing numpy style docstrings
- [x] Writing a simple decorator to find out which functions were slowest
- [x] Using GitHub Actions
- [x] Some additional experience with writing efficient functions
- [ ] How to handle PubMLST API efficiently
- [ ] Unit testing

## To do

- [ ] Add scripts to download Enterobase cgMLST alleles through the PubMLST API
- [ ] Check whether there is an appropriate way to update profiles and ST to hierCC table from Enterobase (keeping in mind fair usage)
- [ ] Implement pytest or other testing framework
- [ ] Improve the correspondence between Enterobase cgMLST typing scheme and chewBBACA. Currently, a lot of alleles are removed in the PrepExternalSchema step of chewBBACA.
- [ ] Improve the speed of comparisons. This is currently done by comparing chunks of the profiles file against the isolate's allelic profile, but this is by far the slowest step in the script.
