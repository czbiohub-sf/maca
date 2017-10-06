library(R.utils)
library(dplyr)

tissue_to_download <- 'Heart'

download_url <- "https://s3.amazonaws.com/czbiohub-maca/gc_table_by_plates_zipped"
rootdir <- "~/projects/maca"

metadata <- read.csv(file = paste(rootdir, "metadata/MACA_Metadata.csv", sep='/'), sep=",", header = TRUE)
colnames(metadata)[1] <- "plate.barcode"

tissue_plates = filter(metadata, tissue == tissue_to_download & mouse.age == 3)[, 'plate.barcode']

for (plate in tissue_plates) {
  res <- try(dwl.status <- download.file(url=sprintf("%s/%s.zip", download_url, plate),
                                         destfile=file.path(rootdir, "data", sprintf("%s.zip", plate))),
             silent=TRUE)
  if (class(res) != "try-error") {
    if (dwl.status != 0) {
      message("\t-> download failed! :(")
    } else {
      unzip(zipfile=file.path(rootdir, 'data', sprintf("%s.zip", plate)), exdir=file.path(rootdir, 'data'))
      file.remove(file.path(rootdir, 'data', sprintf("%s.zip", plate)))
    }
  }
}
