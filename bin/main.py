import argparse
import pandas as pd
import numpy as np
from scipy import stats

def load_sumstats(manifest, pheno_category, phenocode, sex_set):
  # load manifest tsv
  manifest_df = pd.read_csv(manifest, sep="\t")

  # filter for given trait_type, pheno_code, and pheno_sex
  manifest_df = manifest_df[(manifest_df["trait_type"] == pheno_category) &
                            (manifest_df["pheno_code"] == phenocode) &
                            (manifest_df["pheno_sex"] == sex_set)]

  # get aws_link and aws_link_tabix
  aws_link = manifest_df["aws_link"].iloc[0]
  aws_link_tabix = manifest_df["aws_link_tabix"].iloc[0]

  # load data for those AWS paths as file_sumstats and file_sumstats_index
  file_sumstats = pd.read_csv(aws_link)
  file_sumstats_index = pd.read_csv(aws_link_tabix)

  return file_sumstats, file_sumstats_index

if __name__ == "__main__":
  # take arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("--manifest", default="/seq/vgb/dd/gwas/annotate/gwa/panUKBB/phenotype_manifest.tsv")
  parser.add_argument("--pheno_category", default="biomarkers")
  parser.add_argument("--phenocode", default="30600")
  parser.add_argument("--sex_set", default="both_sexes")
  args = parser.parse_args()

  file_sumstats, file_sumstats_index = load_sumstats(args.manifest,
                                                     args.pheno_category,
                                                     args.phenocode,
                                                     args.sex_set)
