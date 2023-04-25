# Filter >8 month PCB-CC Samples

import pandas as pd
import glob
from datetime import datetime, timedelta
import argparse

parser = argparse.ArgumentParser(description='Finds a list of external (PCB-CC) samples in HiFi_Libraries-Summary.tsv that are older than 8 months and have not been deleted.')
parser.add_argument('--hifi','-hifi', type=str,required=True,
                    help='Path to HiFi_Libraries-Summary table')

args = parser.parse_args()

#hifi table
df = pd.read_csv(args.hifi, sep='\t')
hifi_df = df[~df['Run_ID'].isna()]


pcbcc_df = hifi_df[hifi_df['Project'] == 'PCB-CC']

# PCB-CC older than 8 months
run_id_split = pcbcc_df['Run_ID'].str.split('_', n = 2, expand = True)
pcbcc_df['date'] = run_id_split[1].values.tolist()

pcbcc_df['date'] = pd.to_datetime(pcbcc_df['date'])
eight_months_ago = datetime.now() - timedelta(days=243)
older_dates = (pcbcc_df['date'] <= eight_months_ago)

pcbcc_old = pcbcc_df.loc[older_dates]

def doesExist(rec):
	SUBREADS_PATH = rec['Path_to_subreads']
	pat_glob = f'{SUBREADS_PATH}/*.subreads.bam'
	subread_list = glob.glob(pat_glob)
	if len(subread_list) == 0:
		return False
	else:
		return True

pcbcc_old['DOES_EXIST'] = pcbcc_old.apply(doesExist, axis=1)

external_df = pcbcc_old.loc[pcbcc_old['DOES_EXIST']==True]
external_df['CCS Job ID'] = external_df['CCS Job ID'].astype(int).astype(str)
external_df[['Sample Name','Sample_Cell','Run_ID','CCS Job ID','Barcode']].to_csv("external_to_delete.tsv", header=False, sep="\t", index=False)


