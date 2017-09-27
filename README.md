# maca
Command line utilities for RNA-sequencing data

## Installation

```
git clone https://github.com/czbiohub/maca
cd maca
pip install -e .
```

## Usage
So far, there is only one command, `clean_and_zip`:

```
$ maca clean_and_zip --help
Usage: maca clean_and_zip [OPTIONS]

Options:
  --input-folder TEXT
  --cleaned-folder TEXT  If provided, write the cleaned csvs to this folder
  --zipped-folder TEXT
  --platename TEXT       if provided, only use this platename (good for
                         debugging)
  -h, --help             Show this message and exit.
```

Here's an example usage:

```
$ maca clean_and_zip --platename MAA000932 --input-folder /data1/maca/gc_table_by_plates/ --cleaned-folder ./cleaned
Reading plate MAA000932...
        Wrote cleaned counts and metadata to ./zipped/MAA000932.zip
```

### Example tab-delimited, zipped output to R

Here's an example that outputs tab-delimited files and zips them in a format
that R likes (no label in the first column)

```
$ maca clean_and_zip  --output-folder ../gc_table_by_plates_zipped/ --output-format tab --zipped --rstats
```