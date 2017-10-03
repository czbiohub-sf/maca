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

COUNTS_SUFFIX = '.htseq-count.csv'
METADATA_SUFFIX = '.log.csv'
CELL_DIMENSION = 'columns'


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
    mapping_stats.index = mapping_stats.index.map(lambda x: x.strip())

    # Remove "_S\d+" from the ends of the columns
    # mapping_stats.columns = mapping_stats.columns.str.replace('_S\d+', '')

    # metadata = pd.concat([mapping_stats, htseq.loc[HTSEQ_METADATA_ROWS]])
    # metadata = metadata.dropna(how='all')
    #
    # # htseq outputs have double underscores
    # metadata.index = ['htseq' + x if '__' in x else x.strip()
    #                   for x in metadata.index]
    return counts, mapping_stats


def make_basename(prefix, output_format):
    return prefix + '.' + output_format


def read_csv(csv, cell_dimension=CELL_DIMENSION):
    """Returns a feature-by-cell matrix
    
    By providing the cell dimension, this function will ALWAYS return a 
    feature-by-cell matrix 
    
    Parameters
    ----------
    csv : str
        Filename of the data to read
    cell_dimension : "row" | "col"
        Which dimension contains the cells
    """
    row = cell_dimension.startswith('row')
    col = cell_dimension.startswith('col')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = pd.read_csv(csv, index_col=0)

    # Transpose the table if cells are on the columns and not the rows
    if row and not col:
        table = table.T

    # Return a table with the cells as the columns
    return table


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
@click.argument('input_folder', nargs=1,
                type=click.Path(dir_okay=True, readable=True))
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
@click.option('--counts-suffix', default=COUNTS_SUFFIX,
              help="String at the end of the filename that indicates it is "
                      "a counts matrix, e.g. integers of read counts mapping "
                      "to genes for each cell. Must exactly match the end of "
                      "the file after a plate name, e.g. if there is a period "
                      "separating the plate name and this suffix, the suffix "
                      "should contian the period")
@click.option('--metadata-suffix', default=METADATA_SUFFIX,
              help="String at the end of the filename that indicates it is "
                      "a cell metadata file, e.g. number of reads per cell or "
                      "percent mapped reads. Must exactly match the end of the"
                      " file after a plate name, e.g. if there is a period "
                      "separating the plate name and this suffix, the suffix "
                      "should contian the period")
@click.option('--cell-dimension', default=CELL_DIMENSION,
              help="Specifies whether the cells are the rows or the columns "
                   "on the input files. Default is 'row', valid values are "
                   "'rows', 'cols', 'col', 'columns'")
def clean_and_zip(input_folder, output_folder, platename=None,
                  output_format='csv', rstats=False, zipped=False,
                  counts_suffix=COUNTS_SUFFIX, metadata_suffix=METADATA_SUFFIX,
                  cell_dimension=CELL_DIMENSION):
    """Combines counts and metadata files into single zipped plates files
    
    Always outputs features-by-cell matrices
    """
    import pdb; pdb.set_trace()
    if output_folder is not None and not os.path.exists(output_folder):
        os.mkdir(output_folder)

    if platename is not None:
        print(f'Reading plate {platename}...')
        csv = os.path.join(input_folder, f'{platename}{counts_suffix}')
        htseq = read_csv(csv, cell_dimension)

        csv = os.path.join(input_folder, f'{platename}{metadata_suffix}')
        mapping_stats = read_csv(csv, cell_dimension).dropna(how='all')

        counts, metadata = clean_htseq_mapping_stats(htseq, mapping_stats)
        write_counts_metadata(counts, metadata, output_folder, platename,
                              output_format, rstats, zipped)
    
    else:
        for csv in glob.iglob(os.path.join(input_folder, f'*{counts_suffix}')):
            platename = os.path.basename(csv).split('.')[0]
            print(f'Reading plate {platename}...')
            
            htseq = read_csv(csv, cell_dimension)

            csv = os.path.join(input_folder, f'{platename}{metadata_suffix}')
            mapping_stats = read_csv(csv, cell_dimension).dropna(how='all')

            counts, metadata = clean_htseq_mapping_stats(htseq, mapping_stats)

            write_counts_metadata(counts, metadata, output_folder, platename,
                                  output_format, rstats, zipped)

