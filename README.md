# maca
Scripts for analyzing and annotating data from the MACA project.

## Installation

### Installation via zip file

Click on the button labeled "Clone or download" and select "Download ZIP". Unzip the archive and you should be able to open  `MACA_Plate_Notebook.Rmd` in RStudio.

### Installation from GitHub

If you have git installed and you feel adventerous, you can install this repository using the command:

```
git clone https://github.com/czbiohub/maca
```

From then on you can download updates using `git pull origin master`. This might overwrite modifications you have made to your local copy, so we recommend making a copy of the scripts you are using.

The following commands are only necessary for installing the Python scripts in this repository&mdash;they are not needed for using Seurat.

```
cd maca
pip install -e .
```

## Requirements

Analysis requires [R 3.4](https://cran.cnr.berkeley.edu/), the [Java JDK](http://www.oracle.com/technetwork/java/javase/downloads/index.html), and the Seurat package. We also recommend installing [RStudio](https://www.rstudio.com/) for a nicer interface to R. Once you have installed R you can install Seurat using the command:

```
install.packages('Seurat')
```

## Usage

### Downloading data

Open `download_script.R` and modify the variable `tissue_to_download` to the tissue you are interested in. You can find the list of tissue names in `metadata/MACA_Metadata.csv`. Also change the value of `root_dir` to the location that you downloaded this repository. Run this script to download and unzip all of the available 3-month plate data for the tissue of interest.  

### Analyzing data

Open `MACA_Plate_Notebook.Rmd` in RStudio, modify `rootdir` and `tissue_of_interest` as before, and step through the instructions. If you encounter problems, you can email the MACA mailing list, [open an issue on GitHub](https://github.com/czbiohub/maca/issues/new), or ask the Slack channel.  
