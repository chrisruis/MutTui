#Plots mutation type proportions across a set of spectra within a catalogue
#Takes a catalogue as input and produces a single output file PDF with 6 or 12 panels for
#DNA and RNA spectra, respectively, showing mutation type proportions across the dataset
#Colours are setup to support up to 31 spectra being compared, colours are from https://g-thomson.github.io/Manu/
#Required packages: ggplot2, patchwork, argparse

library(ggplot2)
library(patchwork)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-c", help = "Catalogue of mutation type spectra to be plotted, in csv format")
parser$add_argument("--rna", dest = "rna", action = "store_true", default = FALSE, help = "Specify if using an RNA spectrum, otherwise a DNA spectrum is assumed")
parser$add_argument("-o", help = "Name of output PDF containing plots")
args <- parser$parse_args()

#Import catalogue
s <- read.csv(args$c, row.names = 1)
#Convert to proportions
for (eC in names(s)) {
  s[,eC] <- s[,eC]/sum(s[,eC])}

#Colours for plots
s_colours <- c("#CABEE9", "#7C7189", "#FAE093", "#D04E59", "#BC8E7D", "#2F3D70", "#325756", "#7d9fc2", "#C582B2", "#51806a", "#4d5f8e", "#A092B7", "#7ACCD7", "#115896", "#7C6C65", "#4C4C53", "#BA2F00", "#21282F", "#6C803A", "#7B5C34", "#AB7C47", "#CCAE42", "#D73202",  "#272318", "#D3CDBF","#7D9D33", "#CED38C", "#DCC949", "#BCA888", "#CD8862", "#775B24")

#Extract and plot mutation types
ca <- data.frame("Spectra" = names(s), "Proportion" = matrix(s["C>A",]))
caPlot <- ggplot(ca, aes(x = Spectra, y = Proportion, fill = Spectra)) + geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = s_colours) + guides(fill = "none") + ggtitle("C>A") + theme(axis.title = element_blank(), plot.title = element_text(size = 8))
cg <- data.frame("Spectra" = names(s), "Proportion" = matrix(s["C>G",]))
cgPlot <- ggplot(cg, aes(x = Spectra, y = Proportion, fill = Spectra)) + geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = s_colours) + guides(fill = "none") + ggtitle("C>G") + theme(axis.title = element_blank(), plot.title = element_text(size = 8))
ct <- data.frame("Spectra" = names(s), "Proportion" = matrix(s["C>T",]))
ctPlot <- ggplot(ct, aes(x = Spectra, y = Proportion, fill = Spectra)) + geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = s_colours) + guides(fill = "none") + ggtitle("C>T") + theme(axis.title = element_blank(), plot.title = element_text(size = 8))
ta <- data.frame("Spectra" = names(s), "Proportion" = matrix(s["T>A",]))
taPlot <- ggplot(ta, aes(x = Spectra, y = Proportion, fill = Spectra)) + geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = s_colours) + guides(fill = "none") + ggtitle("T>A") + theme(axis.title = element_blank(), plot.title = element_text(size = 8))
tc <- data.frame("Spectra" = names(s), "Proportion" = matrix(s["T>C",]))
tcPlot <- ggplot(tc, aes(x = Spectra, y = Proportion, fill = Spectra)) + geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = s_colours) + guides(fill = "none") + ggtitle("T>C") + theme(axis.title = element_blank(), plot.title = element_text(size = 8))
tg <- data.frame("Spectra" = names(s), "Proportion" = matrix(s["T>G",]))
tgPlot <- ggplot(tg, aes(x = Spectra, y = Proportion, fill = Spectra)) + geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = s_colours) + guides(fill = "none") + ggtitle("T>G") + theme(axis.title = element_blank(), plot.title = element_text(size = 8))

if (args$rna) {
  gt <- data.frame("Spectra" = names(s), "Proportion" = matrix(s["G>T",]))
  gtPlot <- ggplot(gt, aes(x = Spectra, y = Proportion, fill = Spectra)) + geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = s_colours) + guides(fill = "none") + ggtitle("G>T") + theme(axis.title = element_blank(), plot.title = element_text(size = 8))
  gc <- data.frame("Spectra" = names(s), "Proportion" = matrix(s["G>C",]))
  gcPlot <- ggplot(gc, aes(x = Spectra, y = Proportion, fill = Spectra)) + geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = s_colours) + guides(fill = "none") + ggtitle("G>C") + theme(axis.title = element_blank(), plot.title = element_text(size = 8))
  ga <- data.frame("Spectra" = names(s), "Proportion" = matrix(s["G>A",]))
  gaPlot <- ggplot(ga, aes(x = Spectra, y = Proportion, fill = Spectra)) + geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = s_colours) + guides(fill = "none") + ggtitle("G>A") + theme(axis.title = element_blank(), plot.title = element_text(size = 8))
  at <- data.frame("Spectra" = names(s), "Proportion" = matrix(s["A>T",]))
  atPlot <- ggplot(at, aes(x = Spectra, y = Proportion, fill = Spectra)) + geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = s_colours) + guides(fill = "none") + ggtitle("A>T") + theme(axis.title = element_blank(), plot.title = element_text(size = 8))
  ag <- data.frame("Spectra" = names(s), "Proportion" = matrix(s["A>G",]))
  agPlot <- ggplot(ag, aes(x = Spectra, y = Proportion, fill = Spectra)) + geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = s_colours) + guides(fill = "none") + ggtitle("A>G") + theme(axis.title = element_blank(), plot.title = element_text(size = 8))
  ac <- data.frame("Spectra" = names(s), "Proportion" = matrix(s["A>C",]))
  acPlot <- ggplot(ac, aes(x = Spectra, y = Proportion, fill = Spectra)) + geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = s_colours) + guides(fill = "none") + ggtitle("A>C") + theme(axis.title = element_blank(), plot.title = element_text(size = 8))}

if (args$rna) {
  cPlot <- caPlot + cgPlot + ctPlot + taPlot + tcPlot + tgPlot + gtPlot + gcPlot + gaPlot + atPlot + agPlot + acPlot + plot_layout(ncol = 2)} else {
  cPlot <- caPlot + cgPlot + ctPlot + taPlot + tcPlot + tgPlot + plot_layout(ncol = 1)}

pdf(args$o)
print(cPlot)
dev.off()