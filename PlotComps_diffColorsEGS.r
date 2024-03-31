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

gc3$coverage <- gsub(".*Cv(\\d+)_.*", "\\1", gc3$SequenceID) 


gc3$GC3.Degen <- as.numeric(gc3$GC3.Degen)
gc3$ObsWrightENc_6Fold <- as.numeric(gc3$ObsWrightENc_6Fold)

#The data for the null expectation curve will be in the same folder as above
gc3$taxon_a <- gc3$taxon
for (i in seq_len(nrow(spp))) {
  gc3$taxon_a <- gsub(spp$ten_digit_code[i], spp$Species[i], gc3$taxon_a)
}


gc3$Lifestage <- sapply(strsplit(gc3$taxon_a, " "), function(x) x[1])
gc3$Nucleus <- sapply(strsplit(gc3$taxon_a, " "), function(x) x[2])



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

gc4$coverage <- gsub(".*Cv(\\d+)_.*", "\\1", gc4$SequenceID) 

gc4$Lifestage <- sapply(strsplit(gc4$taxon_a, " "), function(x) x[1])
gc4$Nucleus <- sapply(strsplit(gc4$taxon_a, " "), function(x) x[2])

gc4$GF <- sapply(strsplit(gc4$SequenceID, "_"), function(x) x[length(x)])

gc4 <- gc4 %>%
  group_by(GF) %>%
  mutate(mean_cov_by_gf = round(mean(as.numeric(coverage))))

gc4 <- gc4 %>%
  group_by(GF, Lifestage) %>%
  mutate(mean_cov_by_gf_ls = round(mean(as.numeric(coverage))))


gc4$LsExpRank <- gc4|>
  group_by(mean_cov_by_gf_ls)|>
  group_indices()


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

gc5$coverage <- gsub(".*Cv(\\d+)_.*", "\\1", gc5$SequenceID) 


gc5$GC3.Degen <- as.numeric(gc5$GC3.Degen)
gc5$ObsWrightENc_6Fold <- as.numeric(gc5$ObsWrightENc_6Fold)

gc5$taxon_a <- gc5$taxon
for (i in seq_len(nrow(spp))) {
  gc5$taxon_a <- gsub(spp$ten_digit_code[i], spp$Species[i], gc5$taxon_a)
}

gc5$coverage <- gsub(".*Cv(\\d+)_.*", "\\1", gc5$SequenceID) 

gc5$Lifestage <- sapply(strsplit(gc5$taxon_a, " "), function(x) x[1])
gc5$Nucleus <- sapply(strsplit(gc5$taxon_a, " "), function(x) x[2])

gc5$GF <- sapply(strsplit(gc5$SequenceID, "_"), function(x) x[length(x)])

gc5 <- gc5 %>%
  group_by(GF) %>%
  mutate(mean_cov_by_gf = round(mean(as.numeric(coverage))))

gc5 <- gc5 %>%
  group_by(GF, Lifestage) %>%
  mutate(mean_cov_by_gf_ls = round(mean(as.numeric(coverage))))


gc5$LsExpRank <- gc5|>
  group_by(mean_cov_by_gf_ls)|>
  group_indices()


gc5 <- gc5 %>%
  group_by(taxon) %>%
  mutate(taxon_c = paste0(taxon_a,'\n',n()))

#The data for the null expectation curve will be in the same folder as above
enc_null5 <- data.frame(read_tsv('ENc.Null.tsv'))




data_all <- data.frame(rbind(gc3,          # Combine both data frames
                             gc4,
                             gc5),
                       col = c(rep("#555454", nrow(gc3)),
                               rep("blue", nrow(gc4)),
                               rep("red", nrow(gc5))),
                       size = c(rep(0.25, nrow(gc3)),
                                rep(0.25, nrow(gc4)),
                                rep(0.75, nrow(gc5))),
                       alpha = c(rep(0.1, nrow(gc3)),
                                rep(0.5, nrow(gc4)),
                                rep(0.5, nrow(gc5))),
                       data_set = c(rep("R2G", nrow(gc3)),
                                 rep("Epi", nrow(gc4)),
                                 rep("Other Trees", nrow(gc5))))


gc3_plot <- ggplot(data = data_all, mapping = aes(as.numeric(GC3.Degen), as.numeric(ObsWrightENc_6Fold))) +
  geom_point(col = data_all$col, size = data_all$size, alpha = data_all$alpha) +
  geom_line(data = enc_null, aes(GC3, ENc)) +
  theme_classic() +
  labs(x = '%GC at 3rd-pos 4-fold sites', y = 'Observed Wright ENc (6Fold)') +
  theme(
    legend.position = 'none'
  ) +
  facet_wrap(~taxon_a)

gc3_plot


# silly little plots about GC3 and cov
# gc3 = whole r2g

# gc4 = epi

# gc5 = diff 

b <- ggplot(data_all, aes(x=data_set, y=as.numeric(GC3.Degen))) + 
  geom_boxplot()
b

h <- ggplot(data_all, aes(x=as.numeric(GC3.Degen), fill=data_set)) +
  geom_histogram(alpha=0.5, position="identity")
h

e <- ggplot(gc4, aes(x=as.numeric(GC3.Degen))) +
  geom_histogram(alpha=0.5, position="identity")
e

o <- ggplot(gc5, aes(x=as.numeric(GC3.Degen))) +
  geom_histogram(alpha=0.5, position="identity")
o


cov_gc_epi <- ggplot(gc4, aes(x=as.numeric(GC3.Degen), y = as.numeric(coverage))) +
  geom_point() + ggtitle("Epi grabbed clades: sequence GC3 degen. and k-mer coverage")

cov_gc_epi

cov_gc_other <- ggplot(gc5, aes(x=as.numeric(GC3.Degen), y = as.numeric(coverage))) +
  geom_point() + ggtitle("diff OG list grabbed clades: sequence GC3 degen. and k-mer coverage")

cov_gc_other

cov_gc_r2g <- ggplot(gc3, aes(x=as.numeric(GC3.Degen), y = as.numeric(coverage))) +
  geom_point() + ggtitle("whole R2G: sequence GC3 degen. and k-mer coverage")

cov_gc_r2g


#  cov expression plots 

#starting with epi sequences

# ranking did not work
epi_rank <- ggplot(gc4, aes(x = reorder(LsExpRank, -LsExpRank), y = as.numeric(mean_cov_by_gf_ls), color=Lifestage)) +
  geom_jitter()
epi_rank


# EPI

ggplot(gc4, aes(x= reorder(GF, -mean_cov_by_gf), mean_cov_by_gf_ls, group = Lifestage, color = Lifestage)) +
  geom_line() +
  scale_colour_manual(values = c("#C52972", "#3F9DDA", "#7DAC36"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Gene family, ordered by overall mean k-mer") + ylab("Mean k-mer coverage")+
  ggtitle("Coverage by life stage of sequences in Epigenetic grabbed clades")

# OTHER

ggplot(gc5, aes(x= reorder(GF, -mean_cov_by_gf), mean_cov_by_gf_ls, group = Lifestage, color = Lifestage)) +
  geom_line() +
  scale_colour_manual(values = c("#C52972", "#3F9DDA", "#7DAC36"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Gene family, ordered by overall mean k-mer") + ylab("Mean k-mer coverage")+
  ggtitle("Coverage by life stage of sequences in other grabbed clades")


gc4$LsExpRank

ggplot(gc4, aes(x = GC3.Degen, y = as.numeric(coverage))) +
  geom_jitter()


