import numpy as np
import re


ANNOTATION = 'annotation'
SUBANNOTATION = 'subannotation'

DONT_STRIP_NUMBERS = "Thymus", "Trachea"


def clean_labels(x, strip_numbers=True):
    try:
        if strip_numbers:
            x = x.strip('0123456789.')
        return x.lower()\
            .strip()\
            .strip('()-_?')\
            .replace(' ', '_')\
            .replace('.', '_')\
            .replace('__', '_')

            # .replace('/', ' or ')\
    except AttributeError:
        return x


def protect_cd8_cd4(x):
    # print(type(x))
    if isinstance(x, str):
        # print('replacing!')
        x = x.replace('cd8', 'cd8+')
        x = x.replace('cd4', 'cd4+')
        # print(x)
        return x
    else:
        return x


def print_sizes(df, message='\t--- within clean_annotation ---'):
    print(message)
    try:
        print(df.fillna('-').groupby([ANNOTATION, SUBANNOTATION]).size())
    except KeyError:
        print(df.fillna('-').groupby([ANNOTATION]).size())


def _fix_nk_cells(df):
    """Annotate natural killer (NK) cells as a subtype of T cells"""
    rows = df[ANNOTATION] == 'natural_killer_cells'
    df.loc[rows, ANNOTATION] = 't_cells'
    df.loc[rows, SUBANNOTATION] = 'natural_killer_cells'
    return df


def clean_annotation(df, tissue, debug=False):
    strip_numbers = tissue not in DONT_STRIP_NUMBERS
    df = df.applymap(lambda x: clean_labels(x, strip_numbers=strip_numbers))

    if debug:
        print_sizes(df, message='\t\t--- after cleaning labels ---')

    if tissue == "Thymus":
        df = df.applymap(protect_cd8_cd4)

    # if debug:
    #     print_sizes(df, message='\t--- after cleaning thymus ---')

    # --- Aorta ---
    if tissue == "Aorta":
        df[ANNOTATION] = df[ANNOTATION].str.replace('hematopoetic',
                                                    'hematopoietic')
        rows = df[ANNOTATION] == 'heterogenous_group_of_cells'
        df.loc[rows, ANNOTATION] = df.loc[rows, SUBANNOTATION]
        df.loc[rows, SUBANNOTATION] = np.nan

        df[ANNOTATION] = df[ANNOTATION].replace('erythroblasts_and_adipocytes',
                                                'unknown_cells')

    # --- Bladder ---
    elif tissue == "Bladder":
        df[SUBANNOTATION] = df[ANNOTATION].str.extract(
            '(?P<subannotation>[ab]$)', expand=False)
        #         print(df.query('annotation == "Basal"').head())
        df[ANNOTATION] = df[ANNOTATION].str.rstrip('ab12')
        #         print(df.query('annotation == "Basal"').head())
        df[ANNOTATION] = df[ANNOTATION].map(lambda x: x + ' cells' if not x.endswith('cells') else x)

    # --- Brain_FACS_neurons ----
    elif tissue == "Brain_Neurons":
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'endothelial', 'endothelial_cells')
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'npc', 'neural_progenitor_cell')
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'berg_glia', 'bergmann_glia')
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'doublet', 'undetermined')
        df[ANNOTATION] = df[ANNOTATION].replace(
            'opcs', 'oligodendrocyte_progenitor_cells')

        df[SUBANNOTATION] = df[SUBANNOTATION].replace(
            'vsmcs', 'vascular_smooth_muscle_cells')

        rows = df[SUBANNOTATION] == 'myelinating_oligodendrocytes'
        df.loc[rows, ANNOTATION] = 'oligodendrocytes'

    # --- Colon ---
    elif tissue == "Colon":
        pattern = '(?P<annotation>[a-zA-Z_-]+)'
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'undiff', 'undifferentiated')
        df[ANNOTATION] = df.annotation.str.extract(pattern, expand=False)

        df[SUBANNOTATION] = np.nan
        df[SUBANNOTATION] = df.subannotation.replace('', np.nan)
        rows = df.annotation == 'cycling_undifferentiated_cell'
        df.loc[rows, SUBANNOTATION] = 'proliferating'
        rows = df.annotation == 'non-cycling_undifferentiated_cell'
        df.loc[rows, SUBANNOTATION] = 'quiescent'

        rows = df.annotation.str.contains('undifferentiated')
        df.loc[rows, ANNOTATION] = 'undifferentiated_cell'

        df[ANNOTATION] = df[ANNOTATION].str.rstrip('_') + 's'

    # --- Diaphragm ---
    elif tissue == 'Diaphragm':
        pattern = '(?P<annotation>[a-zA-Z /&]+)(?P<subannotation>\d?)'
        df[ANNOTATION] = df[ANNOTATION].replace('b cells & t-cells',
                                                    'immune_cells')
        df = df.annotation.str.extract(pattern, expand=False)
        df[SUBANNOTATION] = df.subannotation.replace('', np.nan)
        df[SUBANNOTATION] = df.subannotation.replace('12', np.nan)
        df[ANNOTATION] = df[ANNOTATION].replace('fibro/adipogenic',
                                                'mesenchymal_stem_cells')

    # --- Fat ---
    elif tissue == 'Fat':
        df[ANNOTATION] = df[ANNOTATION].str.replace('-', ' ')
        df[ANNOTATION] = df[ANNOTATION].str.replace('mono/macro/dcs',
                                                        'myeloid_cells')

        rows = df[ANNOTATION].str.contains('nk')
        # df.loc[rows, ANNOTATION] = 't_cells'
        df.loc[rows, ANNOTATION] = 'natural killer cells'

        # Epithelial, endothelial --> endothelial_cells, epithelial_cells
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'thelial', 'thelial_cells')
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'muscle?', 'muscle_cells')

    # --- Heart ---
    elif tissue == 'Heart':
        df[ANNOTATION] = df.annotation.str.replace('fb', 'fibroblasts')
        df[ANNOTATION] = df.annotation.str.replace('edc',
                                                     'endothelial_cells')
        df[ANNOTATION] = df.annotation.str.replace('cm', 'cardiomyocyte')
        df[ANNOTATION] = df.annotation.str.replace('smc',
                                                     'smooth_muscle_cell')
        df[ANNOTATION] = df.annotation.replace('myofibroblast',
                                                 'myofibroblasts')

        # Extract endothelial cells' subannotations
        rows = df[ANNOTATION].str.endswith('endothelial_cells')
        pattern = '(?P<subannotation>[a-zA-Z_]+)_endothelial_cells'
        df.loc[rows, SUBANNOTATION] = df[ANNOTATION].str.extract(pattern, expand=False)
        df.loc[rows, ANNOTATION] = 'endothelial_cells'

        df[ANNOTATION] = df[ANNOTATION].replace('blood_cells',
                                                'erythrocytes')
        df[ANNOTATION] = df[ANNOTATION].replace('red_blood_cells',
                                                'erythrocytes')

        # Move subannotation to annotation for cluster of both fibroblasts and
        # smooth muscle cells
        rows = df[ANNOTATION] == 'fibroblasts_smooth_muscle_cell'
        # Add 's' for "fibroblast" --> "fibroblasts"
        df.loc[rows, ANNOTATION] = df.loc[rows, SUBANNOTATION].map(
            lambda x: x+'s' if not x.endswith('s') else x)
        df.loc[rows, SUBANNOTATION] = np.nan

        rows = df[SUBANNOTATION].str.lower().str.contains('jun').fillna(False)
        df.loc[rows, SUBANNOTATION] = np.nan

        # Deal with Fb_1 and Immune_Cells_2
        rows = df.annotation.str.contains(r'\d$')
        pattern = '(?P<annotation>[a-zA-Z_]+)_\d?'
        df.loc[rows, ANNOTATION] = df.loc[rows, ANNOTATION].str.extract(pattern, expand=False)

        # Deal with edc_3_endocardial and edc_2_coronary_vascular
        rows = df.annotation.str.contains(r'_\d_')
        pattern = '(?P<annotation>[a-zA-Z_]+)_\d_(?P<subannotation>[a-zA-Z_]+)'
        df.loc[rows] = df.loc[rows, ANNOTATION].str.extract(pattern, expand=False)

    # --- Kidney ---
    elif tissue == "Kidney":
        rows = df[ANNOTATION].str.lower().str.contains('tubule')
        df.loc[rows, SUBANNOTATION] = df.loc[rows, ANNOTATION].str.extract(
            '(?P<annotation>.+)(?=tubule)', expand=False)
        df.loc[rows, ANNOTATION] = 'tubule_cells'

        df[ANNOTATION] = df[ANNOTATION].str.replace('feneserated',
                                                    'fenestrated')
        df[ANNOTATION] = df[ANNOTATION].str.replace('capillaries',
                                                    'capillary_cells')

        df[ANNOTATION] = df[ANNOTATION].replace('other_immune', 'immune_cells')

    # --- Liver ---
    elif tissue == "Liver":

        # Plate CSV has no subannotations
        if df[SUBANNOTATION].notnull().any():
            # This is 10x annotations

            # Remove newlines
            df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
                'female', '').str.replace('male', '')
            df[ANNOTATION] = df[ANNOTATION].str.replace(
                'npc', 'non-parenchymal cell')

            df = df.applymap(lambda x: x.lstrip('fm-') if isinstance(x, str) else x)
            rows = df[ANNOTATION].str.contains('hep')
            df.loc[rows, SUBANNOTATION] = df[ANNOTATION].str.extract('-(.+)', expand=False)
            df.loc[rows, ANNOTATION] = 'hepatocytes'

            df[SUBANNOTATION] = df[SUBANNOTATION].replace('midlobule', 'midlobular')
        else:
            # This is plate annotations
            df[ANNOTATION] = df[ANNOTATION].replace('nk/nkt',
                                                    'natural_killer_cells')

    # --- Lung ---
    elif tissue == "Lung":
        # print(sorted(df[ANNOTATION].unique()))

        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'aveolar', 'alveolar')

        rows = df[ANNOTATION].str.contains('alveolar_epithelial')
        df.loc[rows, SUBANNOTATION] = 'alveolar_' + df[ANNOTATION].str.extract(
            r'(type_[iv]+)', expand=False).str.strip()
        df.loc[rows, ANNOTATION] = 'epithelial_cells'

        # Remove newlines
        df[ANNOTATION] = df[ANNOTATION].str.replace('\n', '')
        # df[SUBANNOTATION] = df.annotation.str.extract(
        #     r'(Type [IV]+)').str.strip()
        df[ANNOTATION] = df.annotation.str.replace(
            '(_type_[iv]+)', '').str.strip().map(
            lambda x: x if x.endswith('s') else x + 's')

        df[ANNOTATION] = df[ANNOTATION].str.replace('remaining ', '')

        rows = df[ANNOTATION] == 'alveolar_macrophages'
        df.loc[rows, ANNOTATION] = 'macrophages'
        df.loc[rows, SUBANNOTATION] = 'alveolar'

        rows = df[ANNOTATION] == 'interstital_macrophages'
        df.loc[rows, ANNOTATION] = 'macrophages'
        df.loc[rows, SUBANNOTATION] = 'interstitial'

        rows = df[ANNOTATION].str.startswith('unknown_immune')
        df.loc[rows, ANNOTATION] = 'immune_cells'
        df.loc[rows, SUBANNOTATION] = np.nan

        rows = df[ANNOTATION].str.contains('natural_killer')
        # df.loc[rows, ANNOTATION] = 't_cells'
        df.loc[rows, ANNOTATION] = 'natural_killer_cells'

        rows = df[ANNOTATION].str.contains('doublet')
        df.loc[rows, ANNOTATION] = 'unknown'

    # -- Mammary ---
    elif tissue == "Mammary":
        df[ANNOTATION] = df[ANNOTATION].str.lower().str.replace(
            'epithelial_', '')
        # print(df.fillna('-').groupby([ANNOTATION, SUBANNOTATION]).size())

        rows = df[ANNOTATION].str.startswith('hormone')
        df.loc[rows, ANNOTATION] = 'luminal_cells'
        df.loc[rows, SUBANNOTATION] = 'hormone_responsive'

        rows = df[ANNOTATION].str.contains('luminal_progenitors')
        df.loc[rows, ANNOTATION] = 'luminal_cells'
        df.loc[rows, SUBANNOTATION] = 'progenitors'

        rows = df[ANNOTATION].str.startswith('s100a4+/ccl5+')
        df.loc[rows, ANNOTATION] = 'stromal_cells'
        df.loc[rows, SUBANNOTATION] = 's100a4+/ccl5+'

        rows = df[ANNOTATION].str.contains('\+') & ~df[ANNOTATION].str.endswith('cells')
        pattern = '(?P<annotation>[a-z]+_cells)_(?P<subannotation>[a-z0-9\+\/_]+)'
        df.loc[rows] = df.loc[rows, ANNOTATION].str.extract(pattern, expand=True)

        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace('/', '_and_')

    # --- Marrow ---
    elif tissue == "Marrow":
        try:
            df = df.drop('plate.barcode', axis=1)
        except ValueError:
            pass

        rows = df[ANNOTATION] == 'neutrophils'
        df.loc[rows, SUBANNOTATION] = 'quiescent'

        # Fix all B cell annotations (contain capital B)
        rows = df[ANNOTATION].str.endswith('b')
        subset = df.loc[rows]
        pattern = '(.+)b'
        df.loc[rows, SUBANNOTATION] = subset[ANNOTATION].str.extract(pattern, expand=False)
        df.loc[rows, ANNOTATION] = 'b_cells'
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'monocytes_monocyte-progenitors', 'monocytes')
        df[ANNOTATION] = df[ANNOTATION].replace(
            'stem_progenitors', 'progenitors')

        # Replace T_NK-type populations with simply t_cells
        df[ANNOTATION] = df[ANNOTATION].str.replace('t_nkt', 't_cells')
        df[ANNOTATION] = df[ANNOTATION].str.replace('t_nk', 't_cells')

        # 'immmature_mature' --> "maturing"
        rows = df[SUBANNOTATION] == 'immature_mature'
        df.loc[rows, SUBANNOTATION] = 'maturing'

        # subannotation: MonoProgenitor --> progenitor
        # (since annotation says "monocyte" already)
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'monoprogenitor', 'progenitor')

        # Split on dash for granulocytes and neutrophils
        rows = df[ANNOTATION].str.contains('-')
        # erythroid-progenitors_erythroblasts -->  progenitors_erythroblasts
        pattern = '(?P<subannotation>.+)-(?P<annotation>.+)'
        df.loc[rows] = df.loc[rows, ANNOTATION].str.extract(pattern, expand=True)[[ANNOTATION, SUBANNOTATION]]

        # remove "erythroblasts" from "progenitors_erythroblasts"
        df[ANNOTATION] = df[ANNOTATION].replace('progenitors_erythroblasts',
                                                'progenitors')

        # Remove all numbers
        df[SUBANNOTATION] = df[SUBANNOTATION].str.rstrip('0123456789')
        df[SUBANNOTATION] = df[SUBANNOTATION].map(clean_labels)
        df[SUBANNOTATION] = df[SUBANNOTATION].replace('unknown', np.nan)

        df[SUBANNOTATION] = df[SUBANNOTATION].replace(
            'monocyte', 'mature')
        df[SUBANNOTATION] = df[SUBANNOTATION].replace(
            't', np.nan)
        df[SUBANNOTATION] = df[SUBANNOTATION].replace(
            'resting', 'quiescent')

        rows = df[SUBANNOTATION].str.startswith('nk').fillna(False)
        df.loc[rows, ANNOTATION] = 'natural_killer_cells'
        df.loc[rows, SUBANNOTATION] = np.nan
        # df[SUBANNOTATION] = df[SUBANNOTATION].replace(
        #     'nk', 'natural_killer_cells')

    # --- Muscle ---
    elif tissue == "Muscle":
        df[ANNOTATION] = df[ANNOTATION].fillna('unknown')
        df[ANNOTATION] = df[ANNOTATION].str.replace('-', '')

        # Only plate data has the erroneous jun annotations
        if df[SUBANNOTATION].notnull().any():
            # Remove Jun-negative, Jun-positive annotation from plate annotation
            rows = df[SUBANNOTATION].str.startswith("jun").fillna(False)
            df.loc[rows, SUBANNOTATION] = np.nan

        # Reduce ambiguity in annotation
        df[ANNOTATION] = df[ANNOTATION].replace('fibro/adipogenic_progenitors',
                                                'mesenchymal_stem_cells')

    # --- Pancreas ---
    elif tissue == "Pancreas":
        rows = df[SUBANNOTATION].str.contains('alpha').fillna(False)
        df.loc[rows, SUBANNOTATION] = np.nan

        df[ANNOTATION] = df[ANNOTATION].replace('duct_cells', 'ductal_cells')
        df[ANNOTATION] = df[ANNOTATION].replace('exocrine_cells',
                                                'acinar_cells')

        rows = df[ANNOTATION] == 'immune_cells'
        df.loc[rows, SUBANNOTATION] = np.nan

        rows = df[SUBANNOTATION] == 'pp'
        df.loc[rows, ANNOTATION] = 'pp_cells'
        df.loc[rows, SUBANNOTATION] = np.nan

    # --- Skin ---
    elif tissue == "Skin":
        rows = df[ANNOTATION].str.contains('ife')
        df.loc[rows, SUBANNOTATION] = df.loc[rows, ANNOTATION]\
                                         .str.split().str[0] + '_cells'
        df.loc[rows, ANNOTATION] = 'interfollicular_epidermis'

        rows = df[ANNOTATION] == 'cell_cycle'
        df.loc[rows, ANNOTATION] = 'interfollicular_epidermis'
        df.loc[rows, SUBANNOTATION] = 'proliferating'

        # "Bulge" cells are part of the hair follicle, so make hair follicle
        # the main annotation
        rows = df[ANNOTATION].str.endswith('bulge')
        df.loc[rows, SUBANNOTATION] = df.loc[rows, ANNOTATION]
        df.loc[rows, ANNOTATION] = 'hair_follicle_cells'

    # --- Spleen ---
    elif tissue == "Spleen":
        df[ANNOTATION] = df[ANNOTATION].str.strip('0123456789. ')
        df[ANNOTATION] = df[ANNOTATION].replace('macrophages/dendritic_cells/monocytes_and_others',
                                                'myeloid_cells')
        # print(df.head())

        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'ink', 'invariant_natural_killer_t')

        df[ANNOTATION] = df[ANNOTATION].map(
            lambda x: x if x.endswith('s') else x + 's')
        # Spell check
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            'follilular', 'follicular')
        df[ANNOTATION] = df[ANNOTATION].str.replace(
            't1/t2/follicular', 'follicular')
        rows = df.annotation.str.contains('_[bt]_cells$')
        pattern = '(?P<subannotation>[a-zA-Z_48+]+)_(?P<annotation>[bt]_cells)'
        df.loc[rows] = df.annotation.str.extract(pattern, expand=True)[
            [ANNOTATION, SUBANNOTATION]]

        rows = df[ANNOTATION].str.contains('.+macrophages.+').fillna(False)
        df.loc[rows, ANNOTATION] = 'myeloid_cells'

        rows = df[ANNOTATION] == 'natural_killer_cells'
        df.loc[rows, ANNOTATION] = 't_cells'
        df.loc[rows, SUBANNOTATION] = 'natural_killer_cells'

        rows = df[ANNOTATION] == 'plasmocytoid_dendritic_cells'
        df.loc[rows, ANNOTATION] = 'dendritic_cells'
        df.loc[rows, SUBANNOTATION] = 'plasmocytoid'

    # --- Thymus ---
    elif tissue == "Thymus":
        # Spellcheck
        df[ANNOTATION] = df[ANNOTATION].str.replace('differenation',
                                                    'differentiation')
        df[ANNOTATION] = df[ANNOTATION].str.replace('differentation',
                                                    'differentiation')

        # Deal with T cell group 1 separately since they have subannotations
        df[ANNOTATION] = df[ANNOTATION].replace(
            'thymocyte_1_mix_of_dn4_dp_immaturesps', 't_cells')

        # Make sure cells are plural
        df[ANNOTATION] = df[ANNOTATION].replace('stromal_mesenchymal_stem_cell',
                                                'stromal_stem_cells')
        # Make sure cells are plural
        df[ANNOTATION] = df[ANNOTATION].replace('stromal_mesenchymal_cell',
                                                'stromal_mesenchymal_cells')

        # Deal with thymocyte_2,3,4,5_subannotation
        rows = df[ANNOTATION].str.startswith('thymocyte') | ~df[ANNOTATION].str.contains('cell')
        df.loc[rows, SUBANNOTATION] = df.loc[rows, ANNOTATION].str.extract(
            '(thymocyte)?(?P<subannotation>[a-z0-9_-]+)', expand=False)[SUBANNOTATION]
        df[SUBANNOTATION] = df[SUBANNOTATION].map(lambda x: x + "+" if isinstance(x, str) and re.search('cd[48]$', x) is not None else x)
        df.loc[rows, SUBANNOTATION] = df.loc[rows, SUBANNOTATION].str.lstrip('_0123456789')
        # Deal with dn-stage4b, dn-stage dn-stage_4c

        # SP: single positive CD4 or CD8
        # DP: double positive CD4 and CD8
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'dn', 'double_negative')
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'dp', 'double_positive')
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'sp', 'single_positive')
        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace(
            'sn', 'single_negative')
        df.loc[rows, ANNOTATION] = 't_cells'

        # Replace "mix of cells" subannotation with annotation
        rows = (df[ANNOTATION] == 'mix_of_cells') & (df[SUBANNOTATION].str.startswith('stromal'))
        df.loc[rows, ANNOTATION] = df.loc[rows, SUBANNOTATION]
        df.loc[rows, SUBANNOTATION] = np.nan

        rows = df[SUBANNOTATION].str.startswith('single_positive').fillna(False)
        # df.loc[rows, SUBANNOTATION] = df.loc[rows, ANNOTATION]
        df.loc[rows, ANNOTATION] = 't_cells'

        df[SUBANNOTATION] = df[SUBANNOTATION].str.replace('stage_4c', 'stage4c')

    # --- Tongue ---
    elif tissue == "Tongue":

        rows = df[ANNOTATION] == 'basal/differentiating'
        df.loc[rows, ANNOTATION] = 'basal_cells'
        df.loc[rows, SUBANNOTATION] = np.nan

        rows = df[ANNOTATION] == 'filiform'
        df.loc[rows, ANNOTATION] = 'keratinocytes'
        df.loc[rows, SUBANNOTATION] = 'filiform_cells'

        rows = df[ANNOTATION] == 'stratified/differentiated_keratinocytes'
        df.loc[rows, ANNOTATION] = 'keratinocytes'
        df.loc[rows, SUBANNOTATION] = 'suprabasal_cells'

        rows = df[ANNOTATION] == 'maturing/nonkeratinized'
        df.loc[rows, ANNOTATION] = 'keratinocytes'
        df.loc[rows, SUBANNOTATION] = 'suprabasal'

        # Remaining cells are specific kinds of basal cells
        rows = df[SUBANNOTATION].isnull() & (df[ANNOTATION] != 'unknown_cell')
        df.loc[rows, SUBANNOTATION] = df.loc[rows, ANNOTATION].str.split('_').str[0]
        df.loc[rows, ANNOTATION] = 'basal_cells'

        df[ANNOTATION] = df[ANNOTATION].replace('unknown_cell', 'unknown_cells')

        # df[ANNOTATION] = df[ANNOTATION].str.replace('basal_layer',
        #                                             'basal_cells')
        # df[SUBANNOTATION] = df[SUBANNOTATION].str.strip('0123456789')
        # df[SUBANNOTATION] = df[SUBANNOTATION].str.replace('-', '_')

    # --- Trachea ---
    elif tissue == "Trachea":
        """
        Epcam are epithelial cells.  Within the epithelial cell types, there 
        are 1) Krt5 - basal cells, Scgb1a1 - secretory cells and 3) Foxj1 - 
        ciliated cells.  I removed Foxj1 from the most recent annotations 
        because they weren't well represented in both datasets and not 
        clustering properly in 10X dataset, but this is fine because they are 
        the arguably least interesting/important cell type in trachea.  
        Pdgfrb are stromal cells, apparent heterogeneity that is not well 
        understood.  Pecam1 are endothelial cells and Ptprc are immune cells.
        """
        gene_to_annotation = {'epcam': 'epithelial_cells',
                              # 'scgb1a1': 'epithelial_cells',
                              'pdgfrb_+_ptprc': 'neuroendocrine_cells',
                              'pecam': 'endothelial_cells',
                              'pecam1': 'endothelial_cells',
                              'ptprc': 'immune_cells',
                              'pdgfrb': 'stromal_cells'}
        gene_to_subannotation = {'krt5': 'basal_cells',
                                 'krt': 'basal_cells',
                                 'scgb1a1': 'secretory_cells',
                                 'scgb1a': 'secretory_cells',
                                 'foxj': 'ciliated_cells',
                                 'foxj1': 'ciliated_cells',
                                 }

        for gene, subannotation in gene_to_subannotation.items():
            rows = (df[ANNOTATION] == gene) | (df[SUBANNOTATION] == gene)
            df.loc[rows, ANNOTATION] = 'epithelial_cells'
            df.loc[rows, SUBANNOTATION] = subannotation

        for gene, annotation in gene_to_annotation.items():
            rows = (df[ANNOTATION] == gene) | (df[SUBANNOTATION] == gene)
            df.loc[rows, ANNOTATION] = annotation
            # df.loc[SUBANNOTATION] = subannotation

        # df[ANNOTATION] = df[ANNOTATION].str.replace('immunue', 'immune')

    if debug:
        print_sizes(df, '\t\t\t---- After tissue-specific cleaning ----')

    # Add "cells" if not already plural
    df[ANNOTATION] = df[ANNOTATION].map(
        lambda x: x + ' cells' if not x.endswith('s') else x)

    df[ANNOTATION] = df[ANNOTATION].str.replace('&', 'and')
    # df = _fix_nk_cells(df)
    df = df.applymap(lambda x: clean_labels(x, strip_numbers=strip_numbers))
    df = df.replace('', np.nan)

    # Remove _cells excdpt for b_cells and t_cells
    df[ANNOTATION] = df[ANNOTATION].map(lambda x: x.split('_cells')[0]
        if not x in ('t_cells', 'b_cells') else x)
    return df