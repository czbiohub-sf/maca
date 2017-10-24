library(R.utils)
library(dplyr)

# change this to the tissue you want to download
tissue_to_download <- 'Pancreas'
# change this to the location of your data (same as the MACA_Plate_Notebook location)
rootdir <- "~/projects/maca"
plate_data_dir = paste0(rootdir, '/data/plates')
dir.create(plate_data_dir)

# don't change this one
download_url <- "https://s3.amazonaws.com/czbiohub-maca/gc_table_by_plates_processed_remux_redux_zipped"

# read the metadata to get the plates we want
metadata <- read.csv(file = paste(rootdir, "metadata/MACA_Metadata.csv", sep='/'), sep=",", header = TRUE)
colnames(metadata)[1] <- "plate.barcode"


if(tissue_to_download == 'All'){
  tissue_plates = filter(metadata, mouse.age == 3)[, 'plate.barcode']
} else{
  tissue_plates = filter(metadata, tissue == tissue_to_download & mouse.age == 3)[, 'plate.barcode']
}

# suppress warning messages for samples that havent been sequenced
oldw <- getOption("warn")
options(warn = -1)

for (plate in tissue_plates) {
  res <- try(dwl.status <- download.file(url=sprintf("%s/%s.zip", download_url, plate),
                                         destfile=file.path(plate_data_dir, sprintf("%s.zip", plate))),
             silent=TRUE)
  if (class(res) != "try-error") {
    if (dwl.status != 0) {
      message("\t-> download failed! :(")
    } else {
      unzip(zipfile=file.path(plate_data_dir, sprintf("%s.zip", plate)), exdir=file.path(plate_data_dir))
      file.remove(file.path(plate_data_dir, sprintf("%s.zip", plate)))
    }
  }
}

options(warn = oldw)

