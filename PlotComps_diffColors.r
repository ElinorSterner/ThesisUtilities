#The purpose of this script is to make a simple jellyfish plot from a dataset like the one
#produced by the CUB.py script, also in the Github Utilities folder. This script was written
#by Auden Cote-L'Heureux and last updated in May 2023

library(tidyverse)

#Change this path
setwd('/Users/elinorsterner/Documents/katzlab/allo/for supplement/jelly/CUBOutputWholeR2G/SpreadSheets')

spp <- data.frame(read_tsv('Species_name.tsv'))

#You will need to change name of data frame below
#if you used the CUB.py script, it should be in the spreadsheets folder
#in the output and end in CompTrans.ENc.Raw.tsv
gc3 <- data.frame(read_tsv('CompTrans.ENc.Raw.tsv'))

gc3 <- data.frame(read_tsv('ENc.Raw.tsv'))|>
  mutate(taxon = paste(substr(SequenceID, 1, 5), substr(SequenceID,6,10), sep = '')) #this line reads in your 10-digit codes to a column in the data frame called taxon


gc3$GC3.Degen <- as.numeric(gc3$GC3.Degen)
gc3$ObsWrightENc_6Fold <- as.numeric(gc3$ObsWrightENc_6Fold)

#The data for the null expectation curve will be in the same folder as above
gc3$taxon_a <- gc3$taxon
for (i in seq_len(nrow(spp))) {
  gc3$taxon_a <- gsub(spp$ten_digit_code[i], spp$Species[i], gc3$taxon_a)
}

gc3 <- gc3 %>%
  group_by(taxon) %>%
  mutate(taxon_c = paste0(taxon_a,'\n',n()))
enc_null <- data.frame(read_tsv('ENc.Null.tsv'))


setwd('/Users/elinorsterner/Documents/katzlab/allo/for supplement/jelly/CUBOutputEpiClades/SpreadSheets')
spp <- data.frame(read_tsv('Species_name.tsv'))

#You will need to change name of data frame below
#if you used the CUB.py script, it should be in the spreadsheets folder
#in the output and end in CompTrans.ENc.Raw.tsv
gc4 <- data.frame(read_tsv('CompTrans.ENc.Raw.tsv'))

gc4 <- data.frame(read_tsv('ENc.Raw.tsv'))|>
  mutate(taxon = paste(substr(SequenceID, 1, 5), substr(SequenceID,6,10), sep = '')) #this line reads in your 10-digit codes to a column in the data frame called taxon



gc4$GC3.Degen <- as.numeric(gc4$GC3.Degen)
gc4$ObsWrightENc_6Fold <- as.numeric(gc4$ObsWrightENc_6Fold)

gc4$taxon_a <- gc4$taxon
for (i in seq_len(nrow(spp))) {
  gc4$taxon_a <- gsub(spp$ten_digit_code[i], spp$Species[i], gc4$taxon_a)
}

gc4 <- gc4 %>%
  group_by(taxon) %>%
  mutate(taxon_c = paste0(taxon_a,'\n',n()))

#The data for the null expectation curve will be in the same folder as above
enc_null4 <- data.frame(read_tsv('ENc.Null.tsv'))


# last set of seqs
setwd('/Users/elinorsterner/Documents/katzlab/allo/for supplement/jelly/CUBOutputDiffClades/SpreadSheets')
spp <- data.frame(read_tsv('Species_name.tsv'))

#You will need to change name of data frame below
#if you used the CUB.py script, it should be in the spreadsheets folder
#in the output and end in CompTrans.ENc.Raw.tsv
gc5 <- data.frame(read_tsv('CompTrans.ENc.Raw.tsv'))

gc5 <- data.frame(read_tsv('ENc.Raw.tsv'))|>
  mutate(taxon = paste(substr(SequenceID, 1, 5), substr(SequenceID,6,10), sep = '')) #this line reads in your 10-digit codes to a column in the data frame called taxon



gc5$GC3.Degen <- as.numeric(gc5$GC3.Degen)
gc5$ObsWrightENc_6Fold <- as.numeric(gc5$ObsWrightENc_6Fold)

gc5$taxon_a <- gc5$taxon
for (i in seq_len(nrow(spp))) {
  gc5$taxon_a <- gsub(spp$ten_digit_code[i], spp$Species[i], gc5$taxon_a)
}

gc5 <- gc5 %>%
  group_by(taxon) %>%
  mutate(taxon_c = paste0(taxon_a,'\n',n()))

#The data for the null expectation curve will be in the same folder as above
enc_null5 <- data.frame(read_tsv('ENc.Null.tsv'))




data_all <- data.frame(rbind(gc3,          # Combine both data frames
                             gc4,
                             gc5),
                       col = c(rep("#d9d9d9", nrow(gc3)),
                               rep("blue", nrow(gc4)),
                               rep("red", nrow(gc5))),
                       size = c(rep(1.5, nrow(gc3)),
                                rep(1.5, nrow(gc4)),
                                rep(10, nrow(gc5))))


gc3_plot <- ggplot(data = data_all, mapping = aes(as.numeric(GC3.Degen), as.numeric(ObsWrightENc_6Fold))) +
  geom_point(col = data_all$col, size = 0.25) +
  geom_line(data = enc_null, aes(GC3, ENc)) +
  theme_classic() +
  labs(x = '%GC at 3rd-pos 4-fold sites', y = 'Observed Wright ENc (6Fold)') +
  theme(
    legend.position = 'none'
  ) +
  facet_wrap(~taxon_a)

gc3_plot

