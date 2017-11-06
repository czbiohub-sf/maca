library(R.utils)
library(dplyr)

# change this to the tissue you want to download
tissue_to_download <- 'Liver'
# change this to the location of your data (same as the MACA_Plate_Notebook location)
rootdir <- "~/src/maca"
tenx_data_dir = paste0(rootdir, '/data/10x')
dir.create(tenx_data_dir)

# don't change this one
download_url <- "https://s3.amazonaws.com/czbiohub-maca/10x_data/"

# read the metadata to get the plates we want
metadata <- read.csv(file = paste(rootdir, "metadata/MACA_10x.csv", sep='/'), sep=",", header = TRUE)

if(tissue_to_download == 'All'){
  tissue_channels = filter(metadata, mouse.age == 3)[, 'channel']
} else{
  tissue_channels = filter(metadata, tissue == tissue_to_download & mouse.age == 3)[, 'channel']
}


for (channel in tissue_channels) {
  channel_data_dir = paste0(rootdir, '/data/10x/',channel)
  dir.create(channel_data_dir)
  
  for(filetype in c('barcodes.tsv', 'genes.tsv', 'matrix.mtx')){
    url = paste0(download_url, channel, '/', filetype)
    destfile = file.path(channel_data_dir, filetype)
    
    res <- try(dwl.status <- download.file(url=url,
                                           destfile=destfile),
               silent=TRUE)
  }
}
