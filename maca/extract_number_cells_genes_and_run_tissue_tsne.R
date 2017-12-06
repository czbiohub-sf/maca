rm(list=ls())

setwd('~/Google Drive/MACA uploads')

library(Seurat)
library(cowplot)

objects = c('Aorta/Heart_seurat_tiss.Robj',
            'Bladder/Bladder_seurat_tiss.Robj',
            'Brain_Neurons/ALL_BRAIN_seurat.Robj',
            "Brain_Microglia/microglia_seurat_subtiss.Robj",
            'Colon/Colon_seurat_tiss.Robj',
            "Diaphragm/Muscle_Diaphragm_seurat_tiss.Robj",
            "Fat/Fat_seurat_tiss.Robj",
            'Heart/heart_seurat_tiss.Robj',
            "Kidney/Kidney_seurat_tiss.Robj",
            'Liver/Liver_seurat_tiss-1117.Robj',
            "Lung/SmartSeq2_Lung_seurat_tissue.Robj",
            'Mammary/Mammary_Gland_seurat_tiss.Robj',
            'Marrow/Marrow_seurat_tiss.Robj',
            'Muscle/Muscle_seurat_tiss.Robj',
            'Pancreas/Pancreas_seurat_tiss.Robj',
            'Skin/Skin_seurat_tiss.Robj',
            'Spleen/Spleen_seurat_tiss.Robj',
            'Thymus/Thymus_seurat_tiss.Robj',
            'Tongue/Tongue_seurat_tiss.Robj',
            'Trachea/Trachea_seurat_tiss.Robj'
            )
print(length(objects()))

folder = '~/code/maca/metadata/number_of_cells_reads_genes/'

extract_ngenes_ncells = function(tiss, object){
  tissue_of_interest = strsplit(object, '/')[[1]][1]
  print(tissue_of_interest)
  
  tissue_metadata = data.frame(
    c(dim(tiss@scale.data), dim(tiss@raw.data)[2]), 
    row.names=c('n_genes', 'n_cells_pass_qc', 'n_cells_sequenced'))
  colnames(tissue_metadata) = tissue_of_interest
  write.csv(tissue_metadata, paste0(folder, tissue_of_interest, 
                                    '_cell_numbers.csv'))
  
  write.csv(tiss@meta.data[c('nGene', 'nReads', 'orig.ident')],
            paste0(folder, tissue_of_interest, 
                   '_nreads_ngenes.csv'))
}

cleaned_annotations = read.csv('~/code/maca/metadata/maca_3month_annotations_plates.csv', row.names=1)

figure_folder = '~/Google Drive/MACA_3mo_manuscript/Main figures/figure2/plates/'

plot_annotated_tsne = function(tiss, object_name, tissue_of_interest) {
  p = TSNEPlot(
    object = tiss,
    do.label = FALSE,
    pt.size = 0.05,
    group.by = 'annotation',
    no.legend = TRUE,
    no.axes = TRUE,
    alpha = 0.5,
    do.return = TRUE
  ) #+ geom_point(alpha = 0.1)
  p + labs(title=tissue_of_interest)
  ggsave(
    paste0(
      figure_folder,
      'tsne_annotated_',
      tissue_of_interest,
      '.pdf'
    ),
    width = 2,
    height = 2
  )
  return(p)
}

object_tissue = c("Lung")
subset_tissues = c("Heart", "Muscle")


for (object_name in objects){
  load(object_name)
  print(ls())
  tissue_of_interest = strsplit(object_name, '/')[[1]][1]
  print(tissue_of_interest)
  tissue_annotations = cleaned_annotations[cleaned_annotations$tissue == tissue_of_interest, ]
  
  if( any(tissue_of_interest == object_tissue)){
    extract_ngenes_ncells(tissue, object_name)
    
    # Reassign metadata with cleaned annotations and plot TSNE
    tissue@meta.data = tissue_annotations
    p = plot_annotated_tsne(tissue, object_name, tissue_of_interest)
    rm(list=c('tissue', 'tissue_of_interest'))
  } else{
    extract_ngenes_ncells(tiss, object_name)
    
    # Heart and Aorta annotated some of the same cells but we use the Aorta annotations
    if (any(tissue_of_interest == subset_tissues)){
      tiss = SubsetData(object=tiss, cells.use=rownames(tissue_annotations))
    }

    # Reassign metadata with cleaned annotations and plot TSNE
    tiss@meta.data = tissue_annotations
    p = plot_annotated_tsne(tiss, object_name, tissue_of_interest)
    rm(list=c('tiss', 'tissue_of_interest'))
  }
  
}