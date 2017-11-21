import numpy as np


ANNOTATION = 'annotation'
SUBANNOTATION = 'subannotation'


def clean_labels(x):
    try:
        return x.lower()\
            .strip('0123456789')\
            .strip()\
            .strip('()-_?')\
            .replace('/', ' or ')\
            .replace(' ', '_')\
            .replace('.', '_')\
            .replace('__', '_')
    except AttributeError:
        return x


def _fix_nk_cells(df):
    """Annotate natural killer (NK) cells as a subtype of T cells"""
    rows = df[ANNOTATION] == 'natural_killer_cells'
    df.loc[rows, ANNOTATION] = 't_cells'
    df.loc[rows, SUBANNOTATION] = 'natural_killer_cells'
    return df


def clean_annotation(df, tissue):
    # --- Aorta ---
    if tissue == "Aorta":
        df[ANNOTATION] = df[ANNOTATION].str.replace('hematopoetic',
                                                    'hematopoietic')
        rows = df[ANNOTATION] == 'heterogenous group of cells'
        df.loc[rows, ANNOTATION] = df.loc[rows, SUBANNOTATION]
        df.loc[rows, SUBANNOTATION] = np.nan

    # --- Bladder ---
    if tissue == 'Bladder':
        df[SUBANNOTATION] = df.annotation.str.extract(
            '(?P<subannotation>[AB]\d?$)')
        #         print(df.query('annotation == "Basal"').head())
        df[ANNOTATION] = df.annotation.str.rstrip('AB12')
        #         print(df.query('annotation == "Basal"').head())
        df[ANNOTATION] = df[ANNOTATION] + ' cells'

    # --- Brain_FACS_neurons ----
    if tissue == "Brain_FACS_neurons":
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'endothelial', 'endothelial_cells')
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'NPC', 'neural_progenitor_cell')
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'Berg.Glia', 'Bergmann_glia')
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'doublet', 'undetermined')

    # --- Colon ---
    elif tissue == 'Colon':
        pattern = '(?P<annotation>[a-zA-Z -]+)'
        df[ANNOTATION] = df[ANNOTATION].str.replace('Undiff.',
                                                        'undifferentiated')
        df[ANNOTATION] = df.annotation.str.extract(pattern)

        df[SUBANNOTATION] = np.nan
        df[SUBANNOTATION] = df.subannotation.replace('', np.nan)
        rows = df.annotation == 'Cycling undifferentiated Cell'
        df.loc[rows, SUBANNOTATION] = 'proliferating'
        rows = df.annotation == 'Non-Cycling undifferentiated Cell'
        df.loc[rows, SUBANNOTATION] = 'quiescent'

        rows = df.annotation.str.contains('undifferentiated')
        df.loc[rows, ANNOTATION] = 'undifferentiated_cell'

        df[ANNOTATION] = df[ANNOTATION].str.rstrip() + 's'

    # --- Diaphragm ---
    elif tissue == 'Diaphragm':
        pattern = '(?P<annotation>[a-zA-Z /&]+)(?P<subannotation>\d?)'
        df[ANNOTATION] = df[ANNOTATION].replace('B cells & T-cells',
                                                    'immune cells')
        df = df.annotation.str.extract(pattern)
        df[SUBANNOTATION] = df.subannotation.replace('', np.nan)
        df[SUBANNOTATION] = df.subannotation.replace('12', np.nan)

    # --- Fat ---
    elif tissue == 'Fat':
        df[ANNOTATION] = df[ANNOTATION].str.replace('-', ' ')
        df[ANNOTATION] = df[ANNOTATION].str.replace('Mono/macro/DCs',
                                                        'Myeloid cells')

        rows = df[ANNOTATION].str.contains('NK')
        df.loc[rows, ANNOTATION] = 't_cells'
        df.loc[rows, SUBANNOTATION] = 'natural killer cells'

        # Epithelial, endothelial --> endothelial_cells, epithelial_cells
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'thelial', 'thelial_cells')
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'Muscle?', 'muscle_cells')

    # --- Heart ---
    elif tissue == 'Heart':
        df[ANNOTATION] = df.annotation.str.replace('Fb', 'fibroblasts')
        df[ANNOTATION] = df.annotation.str.replace('Edc',
                                                     'endothelial_cells')
        df[ANNOTATION] = df.annotation.str.replace('CMs', 'cardiomyocytes')
        df[ANNOTATION] = df.annotation.str.replace('SMCs',
                                                     'smooth_muscle_cells')
        df[ANNOTATION] = df.annotation.replace('Myofibroblast',
                                                 'Myofibroblasts')

        # Deal with Fb_1 and Immune_Cells_2
        rows = df.annotation.str.contains(r'\d$')
        pattern = '(?P<annotation>[a-zA-Z_]+)_\d?'
        df.loc[rows, ANNOTATION] = df.loc[rows, ANNOTATION].str.extract(pattern)

        # Deal with edc_3_endocardial and edc_2_coronary_vascular
        rows = df.annotation.str.contains(r'_\d_')
        pattern = '(?P<annotation>[a-zA-Z_]+)_\d_(?P<subannotation>[a-zA-Z_]+)'
        df.loc[rows] = df.loc[rows, ANNOTATION].str.extract(pattern)

    # --- Kidney ---
    elif tissue == "Kidney":
#         df[ANNOTATION] = df[ANNOTATION].str.replace('tubules', 'tubule')
#         rows = df.annotation.str.contains('(', regex=False)
#         pattern = '(?P<annotation>[a-zA-Z ]+)(?P<subannotation> \([a-zA-Z ]+\)?)'
#         df.loc[rows] = df.loc[rows].annotation.str.extract(pattern)
#         df[SUBANNOTATION] = df[ANNOTATION].str.extract(r'(\d)')
#         rows = df.annotation == 'Proximal tubule cells'
#         df.loc[rows, SUBANNOTATION] = '1'
#         df[ANNOTATION] = df.annotation.str.rstrip(' 1234')
        rows = df[ANNOTATION].str.lower().str.startswith('proximal')
        df.loc[rows, SUBANNOTATION] = 'proximal'
        df.loc[rows, ANNOTATION] = 'tubule'
        rows = df[ANNOTATION].str.startswith('THICK')
        df.loc[rows, SUBANNOTATION] = 'thick_ascending'
        df.loc[rows, ANNOTATION] = 'tubule'

    # --- Liver ---
    elif tissue == "Liver":
        # Remove newlines
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'Female', '').str.replace('Male', '')
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'NPC', 'non-parenchymal cell')

    # --- Lung ---
    elif tissue == "Lung":
        # print(sorted(df[ANNOTATION].unique()))

        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'Aveolar', 'Alveolar')

        rows = df[ANNOTATION].str.contains('Alveolar Epithelial')
        df.loc[rows, ANNOTATION] = 'epithelial_cells'
        df.loc[rows, SUBANNOTATION] = 'alveolar'

        # Remove newlines
        df[ANNOTATION] = df[ANNOTATION].str.replace('\n', '')
        # df[SUBANNOTATION] = df.annotation.str.extract(
        #     r'(Type [IV]+)').str.strip()
        df[ANNOTATION] = df.annotation.str.replace(
            '( Type [IV]+)', '').str.strip().map(
            lambda x: x if x.endswith('s') else x + 's')

        df[ANNOTATION] = df[ANNOTATION].str.replace('Remaining ', '')

        rows = df[ANNOTATION] == 'Alveolar Macrophages'
        df.loc[rows, ANNOTATION] = 'macrophages'
        df.loc[rows, SUBANNOTATION] = 'alveolar'

        rows = df[ANNOTATION] == 'Interstital Macrophages'
        df.loc[rows, ANNOTATION] = 'macrophages'
        df.loc[rows, SUBANNOTATION] = 'interstitial'

        rows = df[ANNOTATION] == 'Unknown Immune Is'
        df.loc[rows, ANNOTATION] = 'immune_cells'
        df.loc[rows, SUBANNOTATION] = np.nan

        rows = df[ANNOTATION].str.contains('Natural Killer')
        df.loc[rows, ANNOTATION] = 't_cells'
        df.loc[rows, SUBANNOTATION] = 'natural_killer_cells'

    # -- Mammary ---
    elif tissue == "Mammary_Gland":
        rows = df[ANNOTATION].str.startswith('hormone')
        df.loc[rows, ANNOTATION] = 'luminal cells'
        df.loc[rows, SUBANNOTATION] = 'hormone responsive'

    # --- Marrow ---
    elif tissue == "Marrow":
        df = df.drop('plate.barcode', axis=1)

        rows = df[ANNOTATION] == 'Neutrophils'
        df.loc[rows, SUBANNOTATION] = 'quiescent'

        # Fix all B cell annotations (contain capital B)
        rows = df[ANNOTATION].str.contains('B')
        subset = df.loc[rows]
        pattern = '(.+)-B'
        df.loc[rows, SUBANNOTATION] = subset[ANNOTATION].str.extract(pattern)
        df.loc[rows, ANNOTATION] = 'b_cells'
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'Monocytes_Monocyte-Progenitors', 'monocytes')
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'Stem_Progenitors', 'hematopoietic stem cell')
        df[ANNOTATION] = df[ANNOTATION].str.replace('T_NK', 't_cells')

        # 'Immmature_Mature' --> "maturing"
        rows = df[SUBANNOTATION] == 'Immature_Mature'
        df.loc[rows, SUBANNOTATION] = 'maturing'

        # subannotation: MonoProgenitor --> progenitor
        # (since annotation says "monocyte" already)
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'MonoProgenitor', 'progenitor')

        # Split on dash for granulocytes and neutrophils
        rows = df[ANNOTATION].str.contains('-')
        pattern = '(?P<subannotation>.+)-(?P<annotation>.+)'
        df.loc[rows] = df.loc[rows, ANNOTATION].str.extract(pattern)[[ANNOTATION, SUBANNOTATION]]


        # Remove all numbers
        df[SUBANNOTATION] = df[SUBANNOTATION].str.rstrip('0123456789')
        df[SUBANNOTATION] = df[SUBANNOTATION].map(clean_labels)

        df[SUBANNOTATION] = df[SUBANNOTATION].replace(
            'monocyte', 'mature')
        df[SUBANNOTATION] = df[SUBANNOTATION].replace(
            't', np.nan)
        df[SUBANNOTATION] = df[SUBANNOTATION].replace(
            'resting', 'quiescent')
        df[SUBANNOTATION] = df[SUBANNOTATION].replace(
            'nk', 'natural_killer_cells')

    # --- Muscle ---
    elif tissue == "Muscle":
        df[ANNOTATION] = df[ANNOTATION].str.replace('-', '')

    # --- Pancreas ---
    elif tissue == "Pancreas":
        rows = df[SUBANNOTATION].str.contains('Alpha').fillna(False)
        df.loc[rows, SUBANNOTATION] = ''
        df[SUBANNOTATION] = df[SUBANNOTATION].replace('PP', 'pp_cells')

    # --- Skin ---
    elif tissue == 'Skin':
        rows = df[ANNOTATION].str.contains('IFE')
        df.loc[rows, SUBANNOTATION] = df.loc[rows, ANNOTATION]\
                                         .str.split().str[0] + '_cells'
        df.loc[rows, ANNOTATION] = 'interfollicular epidermis'
        df[ANNOTATION] = df[ANNOTATION].replace('Cell Cycle', 'proliferating_cells')

        rows = df[ANNOTATION] == 'Outer Bulge'
        df.loc[rows, ANNOTATION] == 'bulge_cells'
        df.loc[rows, SUBANNOTATION] == 'outer'
        rows = df[ANNOTATION] == 'Inner Bulge'
        df.loc[rows, ANNOTATION] == 'bulge_cells'
        df.loc[rows, SUBANNOTATION] == 'inner'

    # --- Spleen ---
    elif tissue == 'Spleen':
        df[ANNOTATION] = df[ANNOTATION].map(
            lambda x: x if x.endswith('s') else x + 's')
        # Spell check
        df.annotation = df.annotation.str.replace('Follilular',
                                                  'Follicular').str.replace(
            'T1/T2/Follicular', 'follicular')
        rows = df.annotation.str.contains('[BT] cells')
        pattern = '(?P<subannotation>[a-zA-Z 48+]+) (?P<annotation>[BT] cells)'
        df.loc[rows] = df.annotation.str.extract(pattern)

        rows = df[ANNOTATION].str.contains('Macrophages')
        df.loc[rows, ANNOTATION] = 'myeloid_cells'

        rows = df[ANNOTATION] == 'Natural Killer cells'
        df.loc[rows, ANNOTATION] = 't_cells'
        df.loc[rows, SUBANNOTATION] = 'natural_killer_cells'

    # --- Thymus ---
    elif tissue == 'Thymus':
        # Spellcheck
        df[ANNOTATION] = df[ANNOTATION].str.replace('differenation',
                                                    'differentiation')
        df[ANNOTATION] = df[ANNOTATION].str.replace('differentation',
                                                    'differentiation')

        # Deal with T cell group 1 separately since they have subannotations
        df[ANNOTATION] = df[ANNOTATION].replace(
            'thymocyte_1_mix of DN4_DP_immatureSPs', 't_cells')

        # Deal with thymocyte_2,3,4,5_subannotation
        rows = df[ANNOTATION].str.startswith('thymocyte')
        df.loc[rows, SUBANNOTATION] = df.loc[rows, ANNOTATION].str.extract(
            'thymocyte_\d_(.+)')

        # SP: single positive CD4 or CD8
        # DP: double positive CD4 and CD8
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'DN', 'double_negative')
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'DP', 'double_positive')
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'SP', 'single_positive')
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'SN', 'single_negative')
        df.loc[rows, ANNOTATION] = 't_cells'

    # --- Tongue ---
    elif tissue == "Tongue":
        df[ANNOTATION] = df[ANNOTATION].str.replace('Basal layer',
                                                    'basal_cells')
        df[SUBANNOTATION] = df[SUBANNOTATION].str.strip('0123456789')
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace('-', '_')

    # --- Trachea ---
    elif tissue == "Trachea":
        df[ANNOTATION] = df[ANNOTATION].str.replace('Immunue', 'Immune')


    df[ANNOTATION] = df[ANNOTATION].str.replace('&', 'and')
    df = _fix_nk_cells(df)
    df = df.applymap(clean_labels)
    df = df.replace('', np.nan)
    return df