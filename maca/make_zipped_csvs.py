
#!/usr/bin/env python
import os
from zipfile import ZipFile

import click
import pandas as pd



HTSEQ_ROWS_TO_DROP = ['EXP_ID',
 'WELL_MAPPING',
 '__alignment_not_unique',
 '__ambiguous',
 '__no_feature',
 '__not_aligned',
 '__too_low_aQual']

HTSEQ_METADATA_ROWS = ['__alignment_not_unique',
 '__ambiguous',
 '__no_feature',
 '__not_aligned',
 '__too_low_aQual']


def clean_htseq_mapping_stats(htseq, mapping_stats):
    """Remove metadata-type columns from htseq and combine with mapping_stats
    
    Parameters
    ----------
    htseq : pandas.DataFrame
        A genes-by-cells matrix output from HTSeq
    mapping_stats : pandas.DataFrame
        A features-by-cells concatenated output from STAR aligner
    """
    counts = htseq.loc[htseq.index.difference(HTSEQ_ROWS_TO_DROP)]
    
    mapping_stats = mapping_stats.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    metadata = pd.concat([mapping_stats, htseq.loc[HTSEQ_METADATA_ROWS]])
    metadata = metadata.dropna(how='all')
    
    # htseq outputs have double underscores
    metadata.index = ['htseq' + x if '__' in x else x.strip() for x in metadata.index]

    return counts, metadata


def write_counts_metadata(counts, metadata, zipped_folder, platename, cleaned_folder=None):
    counts_basename = f'{platename}.counts.csv'
    metadata_basename = f'{platename}.metadata.csv'
    
    if cleaned_folder is not None:
        counts.to_csv(os.path.join(cleaned_folder, counts_basename))
        metadata.to_csv(os.path.join(cleaned_folder, metadata_basename))
        print(f'\tWrote {counts_basename} and {metadata_basename} to {cleaned_folder}')
    
    zipname = os.path.join(zipped_folder, f'{platename}.zip')
    with ZipFile(zipname, 'w') as z:
        z.writestr(counts_basename, counts.to_csv())
        z.writestr(metadata_basename , metadata.to_csv())
    print(f'\tWrote cleaned counts and metadata to {zipname}')


@click.command()
@click.option('--input-folder', default='./')
@click.option('--cleaned-folder', default=None, 
              help='If provided, write the cleaned csvs to this folder')
@click.option('--zipped-folder', default='./zipped')
@click.option('--platename', default=None, help='if provided, only use this platename (good for debugging)')
def clean_and_zip(input_folder, cleaned_folder, zipped_folder, platename=None):
    
    if cleaned_folder is not None and not os.path.exists(cleaned_folder):
        os.mkdir(cleaned_folder)
        
    if not os.path.exists(zipped_folder):
        os.mkdir(zipped_folder)
    
    if platename is not None:
        print(f'Reading plate {platename}...')
        csv = os.path.join(input_folder, f'{platename}.htseq-count.csv')
        htseq = pd.read_csv(csv, index_col=0)
        
        csv = os.path.join(input_folder, f'{platename}.log.csv')
        mapping_stats = pd.read_csv(csv, index_col=0).dropna(how='all')
        
        counts, metadata = clean_htseq_mapping_stats(htseq, mapping_stats)
        write_counts_metadata(counts, metadata, zipped_folder, platename, cleaned_folder)
    
    else:
        for csv in glob.iglob(os.path.join(input_folder, '*.htseq-count.csv')):
            platename = os.path.basename(csv).split('.')[0]
            print(f'Reading plate {platename}...')
            
            htseq = pd.read_csv(csv, index_col=0)

            csv = os.path.join(input_folder, f'{platename}.log.csv')
            mapping_stats = pd.read_csv(csv, index_col=0).dropna(how='all')

            counts, metadata = clean_htseq_mapping_stats(htseq, mapping_stats)

            write_counts_metadata(counts, metadata, zipped_folder, cleaned_folder)

            
if __name__ == '__main__':
    clean_and_zip()