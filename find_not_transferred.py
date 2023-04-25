import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Finds a list of samples in HiFi_Libraries-Summary.tsv that have not been transferred.')
parser.add_argument('--hifi','-hifi', type=str,required=True,
                    help='Path to HiFi_Libraries-Summary table')

parser.add_argument('--clinical_table','-clinical', type=str,
					default='/net/eichler/vol28/projects/long_read_archive/nobackups/clinical/cell_table.tsv',
                    help='Path to Clinical Cell table')

parser.add_argument('--nhp_table', '-nhp', type=str,
					default='/net/eichler/vol28/projects/long_read_archive/nobackups/nhp/cell_table.tsv',
                    help='Path to NHP Cell table')

parser.add_argument('--pop_table', '-pop', type=str,
					default='/net/eichler/vol28/projects/long_read_archive/nobackups/pop/cell_table.tsv',
                    help='Path to pop Cell table')

args = parser.parse_args()

#hifi table
df = pd.read_csv(args.hifi, sep='\t')
hifi_df = df.dropna(subset=['Sample_Cell','Run_ID','CCS Job ID'])
hifi_df['Barcode'].fillna(1, inplace=True)

# transfer done cell tables
clinical_table = pd.read_csv(args.clinical_table,  header=0, sep="\t", index_col=0)
nhp_table = pd.read_csv(args.nhp_table,  header=0, sep="\t", index_col=0)
pop_table = pd.read_csv(args.pop_table,  header=0, sep="\t", index_col=0)

cell_table_df = pd.concat([clinical_table,nhp_table,pop_table])
cell_table_df = cell_table_df.dropna(how='all')
cell_table_df['BARCODE'].fillna(1, inplace=True)
cell_table_df = cell_table_df.dropna(subset=['CCS_JOB_ID']) # 2 samples

########################################################################################################
# # Outer left join
hifi_df['CCS Job ID'] = hifi_df['CCS Job ID'].astype(int).astype(str)
cell_table_df['CCS_JOB_ID'] = cell_table_df['CCS_JOB_ID'].astype(int).astype(str)

merged_df = (hifi_df.merge(cell_table_df, left_on=['Sample_Cell','Run_ID','CCS Job ID','Barcode'], right_on=['CELL','RUN_ID','CCS_JOB_ID','BARCODE'], how='left', indicator=True)
     .query('_merge == "left_only"')
     .drop('_merge', 1))

# hifi_df = hifi_df[['Sample Name','Project','Sample_Cell','Run_ID','CCS Job ID','Barcode']]
new_df = merged_df.drop(['CCS_JOB_ID','BARCODE','RUN_ID','CELL'],1)

non_pcbcc_df = new_df[new_df['Project'] != 'PCB-CC']

non_pcbcc_df = non_pcbcc_df[['Sample Name','Sample_Cell','Run_ID','CCS Job ID','Barcode']]
non_pcbcc_df.to_csv('not_transferred.tsv', sep='\t')


