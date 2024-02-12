library(tidyverse)


setwd('/Users/elinorsterner/Documents/katzlab/allo/CladeGrabbing/Diff/')

clades <- read.csv('CladeSizesPerTaxon.csv')
clades <- clades[clades$sr_rh.txt > 2, ]



clade_sizes <- read.csv('clade_sizes.txt')
clade_sizes <- clade_sizes %>% filter(X9 > 2)


clades |>
  ggplot(aes(x = sr_rh.txt))+geom_histogram()+ggtitle("Shared: Number of Foram clades per tree (>2)")


clade_sizes |>
  ggplot(aes(x= X9))+
  geom_histogram()+
  ggtitle("Shared: number of Forams per clade, (>2)")

clade_sizes |>
  ggplot(aes(x= X9))+
  geom_boxplot()+
  ggtitle("Shared: number of Forams per clade (>2)")

