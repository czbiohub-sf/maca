rm(list=ls())

setwd('~/Google Drive/MACA uploads')

library(Seurat)
library(cowplot)

objects = c(
            'Bladder/10x_Bladder_seurat_tiss.Robj',
            'Heart/10x_Heart_seurat_tiss.Robj',
            "Kidney/10x_Kidney_seurat_tiss.Robj",
            'Liver/20171117_plate_10x/10x_Liver_seurat_3m_201711017.Robj',
            "Lung/10x_Lung_seurat_tissue.10X.Robj",
            'Mammary_Gland/10x_Mammary_seurat_tiss.Robj',
            'Marrow/10x_Marrow_seurat_tiss.Robj',
            'Muscle/10x_Muscle_seurat_tiss.Robj',
            'Spleen/10x_Spleen_seurat_tiss.Robj',
            'Thymus/10x_Thymus_final_seurat_tiss.Robj',
            'Tongue/10x_Tongue_seurat_tiss.Robj',
            'Trachea/10x_Trachea_seurat_tiss.Robj'
)
print(length(objects()))

folder = '~/code/maca/metadata/number_of_cells_reads_genes_10x/'

extract_ngenes_ncells = function(tiss, object){
  tissue_of_interest = strsplit(object, '/')[[1]][1]
  print(tissue_of_interest)
  
  tissue_metadata = data.frame(
    c(dim(tiss@scale.data), dim(tiss@raw.data)[2]), 
    row.names=c('n_genes', 'n_cells_pass_qc', 'n_cells_sequenced'))
  colnames(tissue_metadata) = tissue_of_interest
  write.csv(tissue_metadata, paste0(folder, tissue_of_interest, 
                                    '_cell_numbers.csv'))
  
  write.csv(tiss@meta.data[c('nGene', 'nUMI', 'orig.ident')],
            paste0(folder, tissue_of_interest, 
                   '_nreads_ngenes.csv'))
}

cleaned_annotations = read.csv('~/code/maca/metadata/maca_3month_combined_cell_annotations_10x.csv', row.names=1)

figure_folder = '~/Google Drive/MACA_3mo_manuscript/Main figures/figure2/10x/'

plot_annotated_tsne = function(tiss, object_name, tissue_of_interest){
  TSNEPlot(object = tiss, do.label = TRUE, pt.size = 0.5, group.by='annotation_subannotation')
  ggsave(paste0(figure_folder, 'tsne_annotated_', tissue_of_interest, '.pdf'))
}

# Lung used a different variable name for their tissue
object_tissue = c("Lung")
skip_tissues = c("Tongue")


for (object_name in objects){
  load(object_name)
  print(ls())
  tissue_of_interest = strsplit(object_name, '/')[[1]][1]
  print(tissue_of_interest)
  tissue_annotations = cleaned_annotations[cleaned_annotations$tissue == tissue_of_interest, ]
  
  if (any(tissue_of_interest == skip_tissues)){
    next
  }

  if( any(tissue_of_interest == object_tissue)){
    extract_ngenes_ncells(tissue.10X, object_name)
    
    # Reassign metadata with cleaned annotations and plot TSNE
    tissue.10X@meta.data = tissue_annotations
    plot_annotated_tsne(tissue.10X, object_name, tissue_of_interest)
    rm(list=c('tissue.10X', 'tissue_of_interest'))
  } else{
    extract_ngenes_ncells(tiss, object_name)

    if (dim(tissue_annotations)[1] == 0){
      next
    }
    
    # Reassign metadata with cleaned annotations and plot TSNE
    tiss@meta.data = tissue_annotations
    plot_annotated_tsne(tiss, object_name, tissue_of_interest)
    rm(list=c('tiss', 'tissue_of_interest'))
  }
  
}