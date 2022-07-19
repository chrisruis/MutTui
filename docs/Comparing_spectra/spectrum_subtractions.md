# Subtracting spectra

One way to compare 2 spectra is to subtract the mutation proportions in one spectrum from those in the other spectrum

To do this run:
```
MutTui compare -s spectrum1.csv spectrum2.csv -o S1_S2
```

The 2nd spectrum provided with -s is subtracted from the first

2 output files are produced:

1) A spectrum plot of the subtraction, containing the difference in mutation proportion:

<img src="_figures/spectrum_subtraction.png" width = "500">

2) A dot plot comparing the proportions of all contextual mutations between the spectra:

<img src = "_figures/comparison_dot_plot.png" width = "500">

The dotted line in this plot is x=y
