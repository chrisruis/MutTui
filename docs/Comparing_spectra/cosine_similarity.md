# Calculating cosine similarity between spectra

Cosine similarity is a useful measure to compare how similar spectra are to one another at the overall level

To calculate cosine similarity between a set of spectra run:
```
MutTui cluster -s spectrum1.csv spectrum2.csv spectrum3.csv -l S1 S2 S3 -o .
```

Provide the spectra to be compared with -s. The optional labels provided with -l will be the names of the spectra in the output file; provide one label for each input spectrum. If -l is not provided, the file names (including paths) will be used as spectrum names

The cosine similarities between all pairs of spectra are output as cosine_similarity.csv and the cosine distances (1 - cosine similarity) are output as sample_distances.csv

Alternatively, the spectra can first be combined into a catalogue (see [Combining spectra](https://chrisruis.github.io/MutTui/#/Combining_spectra)) which can be provided with -c without -s:
```
MutTui catalogue -c catalogue.csv -o .
```

In this case, the catalogue column names will be used as sample names

# Clustering spectra

**This is still being developed so will likely change soon**

The cluster_spectra.py script as run above will carry out clustering on the distances between spectra
