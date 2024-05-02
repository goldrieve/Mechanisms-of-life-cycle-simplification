import pandas as pd

#Read in data
columns = ['CHROM', 'START', 'END', 'ID', 'MU10', 'Bosendja', 'STIB818', 'Zagora-I-17', 'MHOM-ZM-80-TRPZ-23', 'Alfort', 'AnTat-3-1', 'NKOUA', 'ATCC30019', 'HAMBURG', 'GMOM-ZM-83-TRPZ-317', 'UPA_PO', 'LOGRA', 'MSUS-CI-78-TSW-157', 'MHOM-ZM-83-TRPZ-349', 'MSUS-CI-82-TSW62', 'MBO-NG-74-R10', 'Pakwa', 'MCAM-ET-2013-MU-02', 'MCAM-ET-2013-MU-14', 'MCAM-ET-2013-MU-05', 'E28', 'FEO-AnTat-16-1', 'Kenya', 'FEO', 'American-Strain', 'ROUPO-VAVOUA--80-MURAZ-14', 'MU09', 'Nabe', 'Canadian-Strain', 'Rumphi', 'Philippines', 'SVP', 'MSUS-CI-78-TSW382', 'LiTat-1-5-P9', 'BIM-AnTat-8-1-P8', 'GPAP-CI-82-KP10-29', '348BT', 'MCAP-CI-91-BALEA-2', 'Jua', 'IVMt1', 'MBA', 'MCAM-ET-2013-MU-01', 'STIB-816', '940', 'MBOT-GM-77-GB2', 'STIB-851', 'ATCC-30019', 'MSUS-CI-83-TSW-11', 'STIB247', 'MHOM-SD-82-MUSIKIA-cloneA', 'BoTat-1-1', 'Te-Ap-N-D1', 'STIB805', 'NDMI', 'MCAM-ET-2013-MU-09', '15BT-relapse', 'MCAM-ET-2013-MU-04', 'STIB386', 'ATCC-30023', 'RoTat_1.2', 'Merzouga-56', 'LIGO', 'MHOM-SD-82-BIYAMINA', 'Kazakstan', '57AT', 'EATRO3', 'OVI', 'AnTat-12-1S', 'MCAM-ET-2013-MU-17', 'STIB900', '108AT', '280104', '108BT', '927', 'Colombia', 'EATRO2340', '340AT', 'AnTat-3-3', 'OUSOU', 'Dodola_943', 'ABBA', 'AGAL-CI-78-TCH312']
df = pd.read_table('data/delly/merged/summary.tsv', header=None,names=columns)

# Define clades
evansi_a = ['280104', 'Alfort', 'American-Strain', 'AnTat-3-1', 'AnTat-3-3', 'ATCC-30019', 'ATCC-30023', 'ATCC30019', 'Canadian-Strain', 'Colombia', 'HAMBURG', 'Kazakstan', 'Kenya', 'MCAM-ET-2013-MU-01', 'MCAM-ET-2013-MU-02', 'MCAM-ET-2013-MU-04', 'MCAM-ET-2013-MU-05', 'MCAM-ET-2013-MU-09', 'MCAM-ET-2013-MU-17', 'Merzouga-56', 'MU09', 'Philippines', 'RoTat_1.2', 'STIB-816', 'STIB805', 'STIB818', 'SVP', 'Zagora-I-17']
evansi_b = ['MCAM-ET-2013-MU-14', 'MU10']
evansi_c = ['IVMt1']
equi_ovi = ['940', 'Dodola_943', 'OVI', 'Te-Ap-N-D1']
equi_botat = ['BoTat-1-1', 'E28']
pleo = ['927', '108AT', '108BT', '15BT-relapse', '340AT', '348BT', '57AT', 'ABBA', 'AGAL-CI-78-TCH312', 'AnTat-12-1S', 'BIM-AnTat-8-1-P8', 'Bosendja', 'EATRO2340', 'EATRO3', 'FEO', 'FEO-AnTat-16-1', 'GMOM-ZM-83-TRPZ-317', 'GPAP-CI-82-KP10-29', 'Jua', 'LIGO', 'LiTat-1-5-P9', 'LOGRA', 'MBA', 'MBO-NG-74-R10', 'MBOT-GM-77-GB2', 'MCAP-CI-91-BALEA-2', 'MHOM-SD-82-BIYAMINA', 'MHOM-SD-82-MUSIKIA-cloneA', 'MHOM-ZM-80-TRPZ-23', 'MHOM-ZM-83-TRPZ-349', 'MSUS-CI-78-TSW-157', 'MSUS-CI-78-TSW382', 'MSUS-CI-82-TSW62', 'MSUS-CI-83-TSW-11', 'Nabe', 'NDMI', 'NKOUA', 'OUSOU', 'Pakwa', 'ROUPO-VAVOUA--80-MURAZ-14', 'Rumphi', 'STIB-851', 'STIB247', 'STIB386', 'STIB900', 'UPA_PO']
mask = ((df[pleo] > 1) & (df[pleo] < 3)).all(1)
pleo_df = df[mask]

# Create csv for clade specific CNVs
def create_clade_csv(df, clade, output_file):
    mask = ((df[clade] <= 1) | (df[clade] >= 3)).all(1)
    clade_df = df[mask]
    clade_df['mean'] = clade_df[clade].mean(axis=1)
    clade_df[['CHROM', 'START', 'END', 'mean']].to_csv(output_file, sep=' ', index=False, header=False)

create_clade_csv(pleo_df, evansi_a, "data/circos/evansi_a_cnv.txt")
create_clade_csv(pleo_df, evansi_b, "data/circos/evansi_b_cnv.txt")
create_clade_csv(pleo_df, evansi_c, "data/circos/evansi_c_cnv.txt")
create_clade_csv(pleo_df, equi_ovi, "data/circos/equi_ovi_cnv.txt")
create_clade_csv(pleo_df, equi_botat, "data/circos/equi_botat_cnv.txt")