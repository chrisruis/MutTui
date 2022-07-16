# Combining multiple spectra into a single catalogue

Multiple spectra can be combined into a single catalogue using:
```
python3 combine_spectra_into_catalog.py -s spectrum1.csv spectrum2.csv spectrum3.csv -l S1 S2 S3 -o catalogue.csv
```

Provide the spectra to be combined with -s. Option -l specifies the names of the spectra in the combined catalogue file. Alternatively a file can be provided with -c containing the file name of the spectrum and the name it will be given in the catalogue. If -l and -c are not provided, the spectra will be named with their file names
