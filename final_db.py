# Set up
import pandas as pd
%matplotlib inline
import matplotlib.pyplot as plt
import re
import random
import numpy as np
import os
from functools import reduce

# open the initial database
wd = "/Volumes/matthews/Guy/Raw_data/monomorph/data/final_db/"
df = pd.read_csv(wd + 'gene_data.csv', sep=",", na_values=['-'], low_memory=False)
df.columns = [c.lower() for c in df.columns]
df = df.replace('\*','', regex=True)


df ['difn-bfd3'] = pd.to_numeric(df['difn/tet- (cds only)']) - pd.to_numeric(df['bfd3/tet- (cds only)']) 
df ['difn-bfd6'] = pd.to_numeric(df['difn/tet- (cds only)']) - pd.to_numeric(df['bfd6/tet- (cds only)']) 
df ['difn-pf'] = pd.to_numeric(df['difn/tet- (cds only)']) - pd.to_numeric(df['pf/tet- (cds only)']) 

#Add list of sif genes
sifs = {'gene id':("Tb927.2.1810","Tb927.2.2720","Tb927.2.4020","Tb927.3.4560","Tb927.4.670","Tb927.4.3620","Tb927.4.3630","Tb927.4.3640","Tb927.4.3650","Tb927.5.3580","Tb927.6.2300","Tb927.6.2360","Tb927.7.2100","Tb927.7.7160","Tb927.8.2860","Tb927.9.4080","Tb927.9.7550","Tb927.9.13530","Tb927.10.5930","Tb927.10.5940","Tb927.10.5950","Tb927.10.12100","Tb927.10.15020","Tb927.10.16120","Tb927.11.290","Tb927.11.300","Tb927.11.750","Tb927.11.760","Tb927.11.1640","Tb927.11.2250","Tb927.11.3650","Tb927.11.4610","Tb927.11.6600","Tb927.11.11470","Tb927.11.11480","Tb927.10.12090","Tb927.1.1930","Tb927.11.9270","Tb927.6.4220","Tb927.8.1530","Tb927.10.1740","Tb927.10.1750","Tb927.10.2030","Tb927.10.2530","Tb927.10.12110","Tb927.11.1480","Tb927.11.6610","Tb927.7.2660","Tb927.4.5390","Tb927.8.7020","Tb927.11.2500","Tb927.8.8330","Tb927.11.3570","Tb927.6.400","Tb927.11.6590","Tb927.3.2090","Tb927.3.3410","Tb927.11.12850","Tb927.3.4750","Tb927.10.12260","Tb927.1.2100","Tb927.8.2780")}

sifs = pd.DataFrame(sifs)
sifs['QS_Pathway'] = "Y"
df = pd.merge(df, sifs, how="left", on=["gene id"])

#Add list of sif genes
target = {'gene id':("Tb927.8.1530", "Tb927.5.2580", "Tb927.2.4020", "Tb927.11.6600", "Tb927.4.3650", "Tb927.11.3400")}
target = pd.DataFrame(target)
target['target'] = "Y"
df = pd.merge(df, target, how="left", on=["gene id"])

#Add dnds results
code = pd.read_table(directory + "codeml.txt", sep = "_", names = ["source_id", "txt", "dn", "dn/ds"])
dnds = code['dn'].str.split(' ',expand=True)[[5,6,7,8,9,10]]
dnds["source_id"] = code["source_id"]
dnds = dnds.rename({5:'background', 6:'ev_a', 7:'ev_b', 8:'ivmt1', 9:'ovi', 10:'botat'}, axis='columns')
dnds[["background", "ev_a", "ev_b", "ivmt1", "ovi", "botat"]] = dnds[["background", "ev_a", "ev_b", "ivmt1", "ovi", "botat"]].apply(pd.to_numeric, errors='coerce')
dnds[['gene id', 'type']] = dnds['source_id'].str.split(':', n=1, expand=True)
df = pd.merge(df, dnds, how="left", on=["gene id"])

#Add snpeff results
def process_file(path, suffix):
    df = pd.read_table(path, sep=" ", na_values=['-'], low_memory=False, header=None).iloc[1:][0].str.split('\t', expand=True)
    df.columns = df.iloc[0]
    df = df[1:].add_suffix('_' + suffix)
    df[['gene id', 'type']] = df['TranscriptId_' + suffix].str.split(':', n=1, expand=True)
    df.drop('type', axis=1, inplace=True)
    return df

wd = "/Volumes/matthews/Guy/Raw_data/monomorph/data/snpeff/lineage_specific"
files = ['T.b.evansitypeA.vcf.dir/snpEff_genes.txt', 'T.b.evansitypeB.vcf.dir/snpEff_genes.txt', 'T.b.evansitypeC.vcf.dir/snpEff_genes.txt', 'T.b.equitypeOVI.vcf.dir/snpEff_genes.txt', 'T.b.equitypeBOTAT.vcf.dir/snpEff_genes.txt']
suffixes = ['eva', 'evb', 'evc', 'ovi', 'botat']

dfs = {suffix: process_file(wd + '/' + file, suffix) for file, suffix in zip(files, suffixes)}

# Filter with pleomorphic data
pleo = pd.read_table('/Volumes/matthews/Guy/Raw_data/monomorph/data/snpeff/pleomorph/snpEff_genes.txt', sep=" ", na_values=['-'], low_memory=False, header = None).iloc[1:][0].str.split('\t', expand=True)
pleo.columns = pleo.iloc[0]
pleo = pleo[1:].add_suffix('_pleo')
pleo[['gene id', 'type']] = pleo['TranscriptId_pleo'].str.split(':', n=1, expand=True)
pleo.drop('type', axis=1, inplace=True)

dfList = [df, eva, evb, evc, ovi, botat, pleo]
df = reduce(lambda left,right: left.merge(right, how="left", on='gene id'), dfList)

df[['background','ev_a', 'ev_b', 'ivmt1', 'ovi', 'botat']] = df[['background','ev_a', 'ev_b', 'ivmt1', 'ovi', 'botat']].apply(pd.to_numeric)
# Remove pseudo gene
filtered = df[(df['is pseudo'] != "Yes")]
# Remove VSG
filtered = filtered[~filtered['product description'].str.contains('variant surface')]
# Keep only megabase chromosome
filtered = filtered[filtered['chromosome'].isin(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10','11'])] 

# Filter by predicted mutation impact in pleo
filtered[['ev_a', 'variants_impact_HIGH_pleo', 'variants_impact_MODERATE_eva', 'variants_impact_MODERATE_evb', 'variants_impact_MODERATE_evc', 'variants_impact_MODERATE_ovi', 'variants_impact_MODERATE_botat']] = filtered[['ev_a', 'variants_impact_HIGH_pleo', 'variants_impact_MODERATE_eva', 'variants_impact_MODERATE_evb', 'variants_impact_MODERATE_evc', 'variants_impact_MODERATE_ovi', 'variants_impact_MODERATE_botat']].apply(pd.to_numeric)
dnds_filtered = filtered[(filtered['background'] < 1) & (filtered['variants_impact_HIGH_pleo'] == 0 )]
dnds_filtered[['variants_impact_HIGH_eva', 'variants_impact_HIGH_evb', 'variants_impact_HIGH_evc', 'variants_impact_HIGH_ovi', 'variants_impact_HIGH_botat']] = dnds_filtered[['variants_impact_HIGH_eva', 'variants_impact_HIGH_evb', 'variants_impact_HIGH_evc', 'variants_impact_HIGH_ovi', 'variants_impact_HIGH_botat']].apply(pd.to_numeric)

# Extract positive dnds in monomorphic
dnds_filtered = dnds_filtered[(dnds_filtered['ev_a'] > 1) | (dnds_filtered['ev_b'] > 1) | (dnds_filtered['ivmt1'] > 1) | (dnds_filtered['ovi'] > 1) | (dnds_filtered['botat'] > 1)]
dnds_filtered = dnds_filtered[(dnds_filtered['variants_impact_MODERATE_eva'] > 0) | (dnds_filtered['variants_impact_MODERATE_evb'] > 0) | (dnds_filtered['variants_impact_MODERATE_evc'] > 0) | (dnds_filtered['variants_impact_MODERATE_ovi'] > 0) | (dnds_filtered['variants_impact_MODERATE_botat'] > 0) | (dnds_filtered['variants_impact_HIGH_eva'] > 0) | (dnds_filtered['variants_impact_HIGH_evb'] > 0) | (dnds_filtered['variants_impact_HIGH_evc'] > 0) | (dnds_filtered['variants_impact_HIGH_ovi'] > 0) | (dnds_filtered['variants_impact_HIGH_botat'] > 0)]

# Filter by RITseq data
dnds_filtered = dnds_filtered[(dnds_filtered['difn-bfd3'] <= 0) & (dnds_filtered['difn-bfd6'] <= 0) & (dnds_filtered['difn-pf'] <= 0)]
dnds_filtered[['difn/tet- (cds only)']] = dnds_filtered[['difn/tet- (cds only)']].apply(pd.to_numeric)
dnds_filtered = dnds_filtered[(dnds_filtered['difn/tet- (cds only)'] < -1.5)]

# Create a second category based on QS genes
qs_filtered = filtered

# Extract genes with no high impact mutation in pleomorphic cell line but moderate or high in monomorphic
qs_filtered[['ev_a', 'variants_impact_HIGH_pleo', 'variants_impact_MODERATE_eva', 'variants_impact_MODERATE_evb', 'variants_impact_MODERATE_evc', 'variants_impact_MODERATE_ovi', 'variants_impact_MODERATE_botat']] = qs_filtered[['ev_a', 'variants_impact_HIGH_pleo', 'variants_impact_MODERATE_eva', 'variants_impact_MODERATE_evb', 'variants_impact_MODERATE_evc', 'variants_impact_MODERATE_ovi', 'variants_impact_MODERATE_botat']].apply(pd.to_numeric)
qs_filtered[['variants_impact_HIGH_eva', 'variants_impact_HIGH_evb', 'variants_impact_HIGH_evc', 'variants_impact_HIGH_ovi', 'variants_impact_HIGH_botat']] = qs_filtered[['variants_impact_HIGH_eva', 'variants_impact_HIGH_evb', 'variants_impact_HIGH_evc', 'variants_impact_HIGH_ovi', 'variants_impact_HIGH_botat']].apply(pd.to_numeric)
qs_filtered = qs_filtered[(qs_filtered['variants_impact_HIGH_pleo'] == 0 )]
qs_filtered = filtered[(filtered['variants_impact_MODERATE_eva'] > 0) | (filtered['variants_impact_MODERATE_evb'] > 0) | (filtered['variants_impact_MODERATE_evc'] > 0) | (filtered['variants_impact_MODERATE_ovi'] > 0) | (filtered['variants_impact_MODERATE_botat'] > 0) | (filtered['variants_impact_HIGH_eva'] > 0) | (dnds_filtered['variants_impact_HIGH_evb'] > 0) | (dnds_filtered['variants_impact_HIGH_evc'] > 0) | (dnds_filtered['variants_impact_HIGH_ovi'] > 0) | (filtered['variants_impact_HIGH_botat'] > 0)]
# Filter by predefined QS pathway list
qs_filtered = qs_filtered[(qs_filtered['QS_Pathway'] == "Y")]

# Export
dnds_filtered.to_csv('dnds_targets.csv', index = False)
qs_filtered.to_csv('qs_targets.csv', index = False)
