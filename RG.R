library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)

dir('CIRI','*.CIRI')
samples<-c('S33N', 'S33T', 'S33NT')


ciri<-lapply(samples, function(sample){
  ciri_file<-sprintf('CIRI/%s.CIRI',sample)
  data<-read.table(ciri_file, header=T, comment.char='', sep='\t', nrow=-1)
  data$junction_reads_ID<-NULL
  colnames(data)[5:8]<-paste(sample, colnames(data)[5:8], sep='_')
  data
  })

ciri_df<-Reduce(function(x,y) {
  merge(x,y,all=T, 
        by=c('circRNA_ID', 'chr', 
             'circRNA_start', 'circRNA_end', 
             'circRNA_type', 'gene_id'))},
  ciri)
head(ciri_df)
S33NT_RG<-read.table('CIRI/reads.S33NT.old.stat.gz',
                     comment.char = '', sep='\t', header=F, nrow=-1,
                     col.names=c('circRNA_ID', 'qname', 'mapq', 'cigar', 'pos', 'Sample')) %>%
  mutate(Sample = gsub('^.*:|.sort', '', Sample),
         Sample = paste0(Sample, '_junction_reads_byRG')) %>%
  select(circRNA_ID, qname, Sample) %>%
  unique %>%
  group_by(circRNA_ID, Sample) %>%
  summarise(n=n()) %>%
  spread(Sample, n)

ciri_df_RG<-merge(ciri_df, S33NT_RG, by='circRNA_ID') %>%
  mutate(S33N_junction_reads_byRG=ifelse(is.na(S33N_junction_reads_byRG),0, S33N_junction_reads_byRG),
         S33T_junction_reads_byRG=ifelse(is.na(S33T_junction_reads_byRG),0, S33T_junction_reads_byRG),
         S33N_X.junction_reads =ifelse(is.na(S33N_X.junction_reads),0, S33N_X.junction_reads),
         S33T_X.junction_reads =ifelse(is.na(S33T_X.junction_reads),0, S33T_X.junction_reads),
         NT_total = S33N_X.junction_reads + S33T_X.junction_reads,
         S33NT_X.junction_reads =ifelse(is.na(S33NT_X.junction_reads),0, S33NT_X.junction_reads),
         RG_total = S33N_junction_reads_byRG + S33T_junction_reads_byRG,
         isNTjunctionEqRGtotal = S33NT_X.junction_reads==RG_total,
         isNTtotalEqRGtotal = NT_total==RG_total,
         isNTjunctionEqNTtotal = S33NT_X.junction_reads==NT_total,
         NTjunctionDiffRGtotal = S33NT_X.junction_reads - RG_total,
         NTtotalDiffRGtotal = NT_total - RG_total,
         NTjunctionDiffNTtotal = S33NT_X.junction_reads - NT_total
         )

write.table(ciri_df_RG, file='ciri_RG.txt', row.names=F, sep='\t', quote=F)
