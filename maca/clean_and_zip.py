#!/usr/bin/env python

import glob
import os
import warnings
from zipfile import ZipFile

import click
import pandas as pd

# Suppress all warnings because pandas is annoying


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
    
    mapping_stats = mapping_stats.applymap(lambda x: x.strip() if
    isinstance(x, str) else x)
    metadata = pd.concat([mapping_stats, htseq.loc[HTSEQ_METADATA_ROWS]])
    metadata = metadata.dropna(how='all')
    
    # htseq outputs have double underscores
    metadata.index = ['htseq' + x if '__' in x else x.strip()
                      for x in metadata.index]
    return counts, metadata


def write(dataframe, prefix, folder='./', output_format='csv', rstats=False,
          **kwargs):
    """Output dataframe contents to a file
    
    Parameters
    ----------
    dataframe : pandas.DataFrame
        The data to write
    filename : str
        Basename or relative path to write the folder to, without the extension
    folder : str
        Absolute path folder location
    output_format : 'csv' | 'tsv' | 'tab' | 'table'
        Output format of the data. Default is "csv" and any of 'tsv', 'tab', 
        or 'table' will output tab-delimited files. 
    rstats : bool
        If True, removes the row label because R doesn't like that
    kwargs 
        Any other keyword arguments accepted by pandas.DataFrame.to_csv
    """
    if rstats:
        # Don't label the first column because R doesn't like that
        kwargs.update('index_label', False)
    if output_format == 'tsv' or output_format.startswith('tab'):
        kwargs.update('sep', '\t')

    filename = os.path.join(folder, prefix + '.' + output_format)
    dataframe.to_csv(filename)


def write_counts_metadata(counts, metadata, zipped_folder, platename,
                          cleaned_folder=None, output_format='csv',
                          rstats=False, **kwargs):
    counts_prefix = f'{platename}.counts'
    metadata_prefix = f'{platename}.metadata'
    
    if cleaned_folder is not None:
        write(counts, counts_prefix, cleaned_folder, output_format, rstats,
              **kwargs)
        write(metadata, metadata_prefix, cleaned_folder, output_format, rstats,
              **kwargs)
        print(f'\tWrote {counts_prefix} and {metadata_prefix} to '
              f'{cleaned_folder}')
    
    zipname = os.path.join(zipped_folder, f'{platename}.zip')
    with ZipFile(zipname, 'w') as z:
        z.writestr(counts_prefix, counts.to_csv())
        z.writestr(metadata_prefix , metadata.to_csv())
    print(f'\tWrote cleaned counts and metadata to {zipname}')


@click.command()
@click.option('--input-folder', default='./')
@click.option('--cleaned-folder', default=None, 
              help='If provided, write the cleaned csvs to this folder')
@click.option('--zipped-folder', default='./zipped')
@click.option('--platename', default=None,
              help='if provided, only use this platename (good for debugging)')
@click.option('--output-format', default='csv',
              help="Type of file to create. Accepted values: 'csv', 'tsv'")
@click.option('--rstats', '-r', is_flag=True,
              help='If added, then output to "R" statistical language format, '
                   'which does not have a column name in the first column')
def clean_and_zip(input_folder, cleaned_folder, zipped_folder, platename=None,
                  output_format='csv', rstats=False):
    
    if cleaned_folder is not None and not os.path.exists(cleaned_folder):
        os.mkdir(cleaned_folder)
        
    if not os.path.exists(zipped_folder):
        os.mkdir(zipped_folder)
    
    if platename is not None:
        print(f'Reading plate {platename}...')
        csv = os.path.join(input_folder, f'{platename}.htseq-count.csv')
        htseq = pd.read_csv(csv, index_col=0)
        
        csv = os.path.join(input_folder, f'{platename}.log.csv')
        with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mapping_stats = pd.read_csv(csv, index_col=0).dropna(how='all')
        
        counts, metadata = clean_htseq_mapping_stats(htseq, mapping_stats)
        write_counts_metadata(counts, metadata, zipped_folder, platename,
                              cleaned_folder, output_format, rstats)
    
    else:
        for csv in glob.iglob(os.path.join(input_folder, '*.htseq-count.csv')):
            platename = os.path.basename(csv).split('.')[0]
            print(f'Reading plate {platename}...')
            
            htseq = pd.read_csv(csv, index_col=0)

            csv = os.path.join(input_folder, f'{platename}.log.csv')
            with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    mapping_stats = pd.read_csv(csv, index_col=0).dropna(
                        how='all')

            counts, metadata = clean_htseq_mapping_stats(htseq, mapping_stats)

            write_counts_metadata(counts, metadata, zipped_folder,
                                  cleaned_folder, output_format, rstats)

            
if __name__ == '__main__':
    clean_and_zip()

