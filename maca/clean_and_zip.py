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


def make_basename(prefix, output_format):
    return prefix + '.' + output_format


def read_csv(csv):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return pd.read_csv(csv, index_col=0)


def write_counts_metadata(counts, metadata, folder, platename,
                          output_format='csv',
                          rstats=False, zipped=True,
                          **kwargs):
    """Output counts and metadata contents to a file
    
    Parameters
    ----------
    counts : pandas.DataFrame
        Matrix of gene or transcript counts
    metadata : pandas.DataFrame
        Matrix of cell metadata
    folder : str
        Absolute path folder location
    platename : str
        Prefix to add to the counts and metadata filenames, e.g. if the prefix 
        is 'MAA000901' then the files will be 'MAA000901.counts.csv' and 
        'MAA000901.metadata.csv' 
    output_format : 'csv' | 'tsv' | 'tab' | 'table'
        Output format of the data. Default is "csv" and any of 'tsv', 'tab', 
        or 'table' will output tab-delimited files. 
    rstats : bool
        If True, removes the row label because R doesn't like that
    zipped : bool
        If true, creates a single zip file containing both counts and metadata
    kwargs 
        Any other keyword arguments accepted by pandas.DataFrame.to_csv
    """
    if rstats:
        # Don't label the first column because R doesn't like that
        kwargs['index_label'] = False
    if output_format == 'tsv' or output_format.startswith('tab'):
        kwargs['sep'] = '\t'

    counts_basename = make_basename(f'{platename}.counts', output_format)
    metadata_basename = make_basename(f'{platename}.metadata', output_format)

    counts_filename = os.path.join(folder, counts_basename)
    metadata_filename = os.path.join(folder, metadata_basename)
    
    if not zipped:
        counts.to_csv(counts_filename, **kwargs)
        print('f\tWrote {counts_filename}')
        metadata.to_csv(metadata_filename, **kwargs)
        print('f\tWrote {metadata_filename}')

    if zipped:
        zipname = os.path.join(folder, f'{platename}.zip')
        with ZipFile(zipname, 'w') as z:
            z.writestr(counts_basename, counts.to_csv(**kwargs))
            z.writestr(metadata_basename, metadata.to_csv(**kwargs))
        print(f'\tWrote cleaned counts and metadata to {zipname}')


@click.command()
@click.option('--input-folder', default='./')
@click.option('--output-folder', default='./',
              help='Write the cleaned data to this folder')
@click.option('--platename', default=None,
              help='if provided, only use this platename (good for debugging)')
@click.option('--output-format', default='csv',
              help="Type of file to create. Accepted values: 'csv', 'tsv'")
@click.option('--rstats', '-r', is_flag=True,
              help='If added, then output to "R" statistical language format, '
                   'which does not have a column name in the first column')
@click.option('--zipped', is_flag=True,
              help="If added, make a zip file containing both the counts and "
                   "metadata")
def clean_and_zip(input_folder, output_folder, platename=None,
                  output_format='csv', rstats=False, zipped=False):
    if output_folder is not None and not os.path.exists(output_folder):
        os.mkdir(output_folder)

    if platename is not None:
        print(f'Reading plate {platename}...')
        csv = os.path.join(input_folder, f'{platename}.htseq-count.csv')
        htseq = read_csv(csv)

        csv = os.path.join(input_folder, f'{platename}.log.csv')
        mapping_stats = read_csv(csv).dropna(how='all')

        counts, metadata = clean_htseq_mapping_stats(htseq, mapping_stats)
        write_counts_metadata(counts, metadata, output_folder, platename,
                              output_format, rstats, zipped)
    
    else:
        for csv in glob.iglob(os.path.join(input_folder, '*.htseq-count.csv')):
            platename = os.path.basename(csv).split('.')[0]
            print(f'Reading plate {platename}...')
            
            htseq = read_csv(csv)

            csv = os.path.join(input_folder, f'{platename}.log.csv')
            mapping_stats = read_csv(csv).dropna(how='all')

            counts, metadata = clean_htseq_mapping_stats(htseq, mapping_stats)

            write_counts_metadata(counts, metadata, output_folder, platename,
                                  output_format, rstats, zipped)

            
if __name__ == '__main__':
    clean_and_zip()

