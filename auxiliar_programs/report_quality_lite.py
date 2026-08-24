#!/usr/bin/env python
# coding: utf-8
import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--psas', action='store_true')
args = parser.parse_args()

def load_csv_from_current_directory(filename):
    current_files = os.listdir(os.getcwd())
    if filename in current_files:
        return pd.read_csv(filename, comment='#')
    else:
        return pd.DataFrame()

str_frags = 'frags_mqc.csv'
str_peaks = 'peaks_mqc.csv'
str_fragle = 'fragle_mqc.csv'


df_frags = load_csv_from_current_directory(str_frags)
df_peaks = load_csv_from_current_directory(str_peaks)
df_fragle = load_csv_from_current_directory(str_fragle)

def load_enrichment_csvs():
    current_files = os.listdir(os.getcwd())
    
    enrichment_files = [
        f for f in current_files
        if f.startswith('enrichment') and f.endswith('.csv')]
    
    if enrichment_files:
        dataframes = [pd.read_csv(f) for f in enrichment_files]
        return pd.concat(dataframes, ignore_index=True)
    else:
        return pd.DataFrame()  # return empty if no enrichment files

df_enrichment = load_enrichment_csvs()

def load_psas_tsvs():
    psas_files = [
        f for f in os.listdir(os.getcwd())
        if f.endswith('.psas.tsv')]

    if psas_files:
        dataframes = [pd.read_csv(f, sep='\t', na_values=['NA']) for f in psas_files]
        return pd.concat(dataframes, ignore_index=True)
    return pd.DataFrame(columns=['sample_id', 'psas'])

df_psas = load_psas_tsvs()

def join_sample_dataframes(df_frags, df_peaks, df_fragle, df_enrichment, df_psas):
    merged = pd.merge(df_frags, df_peaks, on='SampleName', how='inner')

    if not df_fragle.empty:
        merged = pd.merge(merged,df_fragle,left_on='SampleName', right_on='Sample_ID')
        
    if not df_enrichment.empty:
        merged = pd.merge(merged, df_enrichment, on='SampleName', how='outer')
    if not df_psas.empty:
        merged = pd.merge(merged, df_psas, left_on='SampleName', right_on='sample_id', how='left')
    elif 'psas' not in merged.columns:
        merged['psas'] = ''
    return merged

dfJoin = join_sample_dataframes(df_frags, df_peaks, df_fragle, df_enrichment, df_psas)
dfJoin = dfJoin.drop(columns=['on_bp', 'off_bp', 'on_reads','Sample_ID', 'off_reads', 'sample_id'], errors='ignore')

dfJoin = dfJoin.rename(columns={
    'SampleName': 'Sample',
    'Fragments': 'TotalFragments',
    'Peaks': 'TotalPeaks',
    'ctDNA_Burden': 'ctDNA',
    'mark': 'Enrichment_Mark',
    'enrichment': 'Enrichment_Score'
})

if args.psas:
    psas_column_index = dfJoin.columns.get_loc('psas')
    psas_scores = pd.to_numeric(dfJoin['psas'], errors='coerce')
    dfJoin.insert(
        psas_column_index + 1,
        'pSAS_Prediction',
        psas_scores.gt(-0.0948).map({True: 'Yes', False: 'No'})
    )
    dfJoin = dfJoin.rename(columns={
        'psas': 'pSAS_Score'
    })
else:
    dfJoin = dfJoin.drop(columns=['psas'], errors='ignore')

dfJoin = dfJoin.where(pd.notnull(dfJoin), '')
dfJoin = dfJoin.sort_values(
    by=['Sample', 'Enrichment_Mark'],
    kind='mergesort'
)

filename = 'QualityMetrics.csv'
dfJoin.to_csv(filename, index=False, encoding='utf-8')
