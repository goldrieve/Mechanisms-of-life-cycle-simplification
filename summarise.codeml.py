import pandas as pd
import math
import numpy as np
import re

#Load files and preprocess
directory = "/Users/s1886853/datastore/Raw_data/monomorph/data/codeml/"
code = pd.read_table(directory + "codeml/codeml.txt", sep = "_", names = ["source_id", "txt", "dn", "dn/ds"])
dnds = code['dn'].str.split(' ',expand=True)[[5,6,7,8,9,10]]
dnds["source_id"] = code["source_id"]
dnds = dnds.rename({5:'background', 6:'ev_a', 7:'ev_b', 8:'ivmt1', 9:'ovi', 10:'botat'}, axis='columns')
dnds[["background", "ev_a", "ev_b", "ivmt1", "ovi", "botat"]] = dnds[["background", "ev_a", "ev_b", "ivmt1", "ovi", "botat"]].apply(pd.to_numeric, errors='coerce')

cds = pd.read_table(directory + "gene_lists/genome.csv", sep = ",")
dnds = dnds.merge(cds)

MG_PV = pd.read_table(directory + "gene_lists/Naguleswaran/MG_vs_PV.csv", sep = ",")
PV_SG = pd.read_table(directory + "gene_lists/Naguleswaran/PV_SG.csv", sep = ",")

D3 = pd.read_table(directory + "gene_lists/ritseq/D3.csv", sep = ",")[["Gene ID"]]
D6 = pd.read_table(directory +  "gene_lists/ritseq/D6.csv", sep = ",")[["Gene ID"]]
PS = pd.read_table(directory + "gene_lists/ritseq/PS.csv", sep = ",")[["Gene ID"]]
DIF = pd.read_table(directory + "gene_lists/ritseq/DIF.csv", sep = ",")[["Gene ID"]]
D3_D6_PS_DIF = pd.read_table(directory + "gene_lists/ritseq/D3_D6_PS_DIF.csv", sep = ",")[["Gene ID"]]

#SIF list from lncRNA paper...update
SIF = ["Tb927.2.1810","Tb927.2.2720 ","Tb927.2.4020","Tb927.3.4560","Tb927.4.670 ","Tb927.4.3620","Tb927.4.3630","Tb927.4.3640","Tb927.4.3650","Tb927.5.3580","Tb927.5.5380","Tb927.6.2300","Tb927.6.2360","Tb927.7.2100","Tb927.7.7160","Tb927.8.2860","Tb927.9.4080","Tb927.9.7550","Tb927.9.13530","Tb927.10.1740","Tb927.10.1750","Tb927.10.2030","Tb927.10.2530","Tb927.10.5930 ","Tb927.10.5940","Tb927.10.5950","Tb927.10.12100","Tb927.10.12110","Tb927.10.15020","Tb927.10.16120","Tb927.11.290","Tb927.11.300","Tb927.11.750","Tb927.11.760","Tb927.11.1480","Tb927.11.1640","Tb927.11.2250","Tb927.11.3650","Tb927.11.4610","Tb927.11.6600","Tb927.11.6610","Tb927.11.11470","Tb927.11.11480"]

#Calculate Log2FC
dnds['T.b.evansi type A'] = np.log2(dnds['ev_a']) - np.log2(dnds['background']) 
dnds['T.b.evansi type B'] = np.log2(dnds['ev_b']) - np.log2(dnds['background']) 
dnds['T.b.evansi type IVM-t1'] = np.log2(dnds['ivmt1']) - np.log2(dnds['background']) 
dnds['T.b.equiperdum type OVI'] = np.log2(dnds['ovi']) - np.log2(dnds['background']) 
dnds['T.b.equiperdum type BoTat'] = np.log2(dnds['botat']) - np.log2(dnds['background']) 

#Create dnds lists of interest
D3 = D3.merge(dnds, how = "inner")
D6 = D6.merge(dnds, how = "inner")
PS = PS.merge(dnds, how = "inner")
DIF = DIF.merge(dnds, how = "inner")
D3_D6_PS_DIF = D3_D6_PS_DIF.merge(dnds, how = "inner")
SIF = pd.DataFrame(SIF, columns=['Gene ID'])
MG_PV_sig = MG_PV.query('padj<0.01 & Log2FoldChange>1.5 | padj<0.01 & Log2FoldChange<1.5')[["Gene ID"]]
PV_SG_sig = PV_SG.query('padj<0.01 & Log2FoldChange>1.5 | padj<0.01 & Log2FoldChange<1.5')[["Gene ID"]]
MG_PV_dnds = MG_PV_sig.merge(dnds, on = "Gene ID", how = "inner")
PV_SG_dnds = PV_SG_sig.merge(dnds, on = "Gene ID", how = "inner")
SIF_dnds = SIF.merge(dnds, on = "Gene ID", how = "inner")
D3_D6_PS_DIF['Category'] = 'Essential'
MG_PV_dnds['Category'] = 'Tsetse'
PV_SG_dnds['Category'] = 'Tsetse'
SIF_dnds['Category'] = 'SIF'

#Combine and save
ridges_dnds = pd.concat([MG_PV_dnds, PV_SG_dnds, D3_D6_PS_DIF, SIF_dnds])
ridges_dnds.to_csv(directory + 'dnds.txt', index=False, header = True, sep = ',')

#Clade dnds

dnds.query('ev_a>1 & background <1').replace(':', ' ', regex=True).replace(',', '', regex=True).replace('\.\.+', ' ', regex=True)['Genomic Location (Gene)'].str.split("(", expand = True)[0].to_csv(directory + '../circos/evansi_a_dnds.txt', index=False, header = False, sep = ',')
dnds.query('ev_b>1 & background <1').replace(':', ' ', regex=True).replace(',', '', regex=True).replace('\.\.+', ' ', regex=True)['Genomic Location (Gene)'].str.split("(", expand = True)[0].to_csv(directory + '../circos/evansi_b_dnds.txt', index=False, header = False, sep = ',')
dnds.query('ivmt1>1 & background <1').replace(':', ' ', regex=True).replace(',', '', regex=True).replace('\.\.+', ' ', regex=True)['Genomic Location (Gene)'].str.split("(", expand = True)[0].to_csv(directory + '../circos/evansi_c_dnds.txt', index=False, header = False, sep = ',')
dnds.query('ovi>1 & background <1').replace(':', ' ', regex=True).replace(',', '', regex=True).replace('\.\.+', ' ', regex=True)['Genomic Location (Gene)'].str.split("(", expand = True)[0].to_csv(directory + '../circos/equi_ovi_dnds.txt', index=False, header = False, sep = ',')
dnds.query('botat>1 & background <1').replace(':', ' ', regex=True).replace(',', '', regex=True).replace('\.\.+', ' ', regex=True)['Genomic Location (Gene)'].str.split("(", expand = True)[0].to_csv(directory + '../circos/equi_botat_dnds.txt', index=False, header = False, sep = ',')
dnds.query('ev_a>1 & ev_b>1 & ivmt1>1 & ovi>1 & botat>1 & background <1').to_csv(directory + '../circos/common_dnds.txt', index=False, header = True, sep = ',')

dnds.query('ev_a>1 & background <1')['Gene ID'].to_csv(directory + '../codeml/clade_dnds/evansi_a_dnds.txt', index=False, header = False, sep = ',')
dnds.query('ev_b>1 & background <1')['Gene ID'].to_csv(directory + '../codeml/clade_dnds/evansi_b_dnds.txt', index=False, header = False, sep = ',')
dnds.query('ivmt1>1 & background <1')['Gene ID'].to_csv(directory + '../codeml/clade_dnds/evansi_c_dnds.txt', index=False, header = False, sep = ',')
dnds.query('ovi>1 & background <1')['Gene ID'].to_csv(directory + '../codeml/clade_dnds/equi_ovi_dnds.txt', index=False, header = False, sep = ',')
dnds.query('botat>1 & background <1')['Gene ID'].to_csv(directory + '../codeml/clade_dnds/equi_botat_dnds.txt', index=False, header = False, sep = ',')
dnds.query('ev_a>1 & ev_b>1 & ivmt1>1 & ovi>1 & botat>1 & background <1')['Gene ID'].to_csv(directory + '../codeml/clade_dnds/common_dnds.txt', index=False, header = False, sep = ',')
