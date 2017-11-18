import numpy as np


def clean_labels(x):
    try:
        return x.lower().strip().strip('()').replace('/', ' or ').replace(' ',
                                                                          '_').replace(
            '.', '_').replace('__', '_')
    except AttributeError:
        return x


def clean_annotation(df, tissue):
    # --- Bladder ---
    if tissue == 'Bladder':
        df['subannotation'] = df.annotation.str.extract(
            '(?P<subannotation>[AB]\d?$)')
        #         print(df.query('annotation == "Basal"').head())
        df['annotation'] = df.annotation.str.rstrip('AB12')
        #         print(df.query('annotation == "Basal"').head())
        df['annotation'] = df['annotation'] + ' cells'

    # --- Colon ---
    elif tissue == 'Colon':
        pattern = '(?P<annotation>[a-zA-Z -]+)(?P<subannotation>\d?)'
        df['annotation'] = df['annotation'].str.replace('Undiff.',
                                                        'undifferentiated')
        df = df.annotation.str.extract(pattern)
        df['subannotation'] = df.subannotation.replace('', np.nan)

    # --- Diaphragm ---
    elif tissue == 'Diaphragm':
        pattern = '(?P<annotation>[a-zA-Z /&]+)(?P<subannotation>\d?)'
        df['annotation'] = df['annotation'].replace('B cells & T-cells',
                                                    'immune cells')
        df = df.annotation.str.extract(pattern)
        df['subannotation'] = df.subannotation.replace('', np.nan)

    # --- Fat ---
    elif tissue == 'Fat':
        df['annotation'] = df['annotation'].str.replace('-', ' ')
        df['annotation'] = df['annotation'].str.replace('Mono/macro/DCs',
                                                        'Myeloid cells')
        df['annotation'] = df['annotation'].str.replace('NK cells',
                                                        'natural killer cells')
    # print(df.groupby('annotation').size().index)

    # --- Heart ---
    elif tissue == 'Heart':
        df['annotation'] = df.annotation.str.replace('Fb', 'fibroblasts')
        df['annotation'] = df.annotation.str.replace('Edc',
                                                     'endothelial_cells')
        df['annotation'] = df.annotation.str.replace('CMs', 'cardiomyocytes')
        df['annotation'] = df.annotation.str.replace('SMCs',
                                                     'smooth_muscle_cells')

        # Deal with Fb_1 and Immune_Cells_2
        rows = df.annotation.str.contains(r'\d$')
        pattern = '(?P<annotation>[a-zA-Z_]+)_(?P<subannotation>\d?)'
        df.loc[rows] = df.loc[rows, 'annotation'].str.extract(pattern)

        # Deal with edc_3_endocardial and edc_2_coronary_vascular
        rows = df.annotation.str.contains(r'_\d_')
        pattern = '(?P<annotation>[a-zA-Z_]+)_\d_(?P<subannotation>[a-zA-Z_]+)'
        df.loc[rows] = df.loc[rows, 'annotation'].str.extract(pattern)

    # --- Kidney ---
    elif tissue == "Kidney":
        #         df['annotation'] = df['annotation'].str.replace('tubules', 'tubule')
        #         rows = df.annotation.str.contains('(', regex=False)
        #         pattern = '(?P<annotation>[a-zA-Z ]+)(?P<subannotation> \([a-zA-Z ]+\)?)'
        #         df.loc[rows] = df.loc[rows].annotation.str.extract(pattern)
        df['subannotation'] = df['annotation'].str.extract(r'(\d)')
        rows = df.annotation == 'Proximal tubule cells'
        df.loc[rows, 'subannotation'] = '1'
        df['annotation'] = df.annotation.str.rstrip(' 1234')

    # --- Liver ---
    elif tissue == "Liver":
        # Remove newlines
        df['subannotation'] = df['subannotation'].str.replace('Female',
                                                              '').str.replace(
            'Male', '')

    # --- Lung ---
    elif tissue == "Lung":
        # Remove newlines
        df['annotation'] = df['annotation'].str.replace('\n', '')
        df['subannotation'] = df.annotation.str.extract(
            r'(Type [IV]+)').str.strip()
        df['annotation'] = df.annotation.str.replace('( Type [IV]+)',
                                                     '').str.strip().map(
            lambda x: x if x.endswith('s') else x + 's')
    # print('-- after cleaning 1 --\n', df.fillna('').groupby(['annotation']).size())

    #     # --- Marrow ---
    #     elif tissue == "Marrow":
    #         rows = df['annotation'].str.contains('B')
    #         df.loc[rows, 'subannotation'] = df.loc[rows, 'subannotation']

    # --- Pancreas ---
    elif tissue == "Pancreas":
        df['subannotation'] = df['subannotation'].str.replace('Alpha - ', '')

    # --- Skin ---
    elif tissue == 'Skin':
        df['annotation'] = df['annotation'].str.replace('IFE',
                                                    'interfollicular epidermis')

    # --- Spleen ---
    elif tissue == 'Spleen':
        df['annotation'] = df['annotation'].map(
            lambda x: x if x.endswith('s') else x + 's')
        df.annotation = df.annotation.str.replace('Follilular',
                                                  'Follicular').str.replace(
            'T1/T2/Follicular', 'follicular')
        rows = df.annotation.str.contains('[BT] cells')
        pattern = '(?P<subannotation>[a-zA-Z 48+]+) (?P<annotation>[BT] cells)'
        df.loc[rows] = df.annotation.str.extract(pattern)

    # --- Tongue ---
    elif tissue == "Tongue":
        df['annotation'] = df['annotation'].str.replace('Basal layer',
                                                        'basal_cells')
    # --- Trachea ---
    elif tissue == "Trachea":
        df['annotation'] = df['annotation'].str.replace('Immunue', 'Immune')

    # --- Thymus ---
    # SP: single positive CD4 or CD8
    # DP: double positive CD4 and CD8
    elif tissue == 'Thymus':
        rows = df['annotation'].str.startswith('thymocyte')
        df.loc[rows, 'subannotation'] = df.loc[rows, 'annotation'].str.extract(
            'thymocyte_\d_(.+)')
        df['subannotation'] = df['subannotation'].str.replace('DN',
                                                              'double_negative')
        df['subannotation'] = df['subannotation'].str.replace('DP',
                                                              'double_positive')
        df.loc[rows, 'annotation'] = 't_cell'

    df['annotation'] = df['annotation'].str.replace('&', 'and')
    df = df.applymap(clean_labels)
    return df