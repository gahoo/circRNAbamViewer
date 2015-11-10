library(UpSetR)
library(plyr)
library(dplyr)
library(tidyr)

ciri_files<-dir('CIRI')
ciri<-lapply(ciri_files, function(ciri_file){
  read.table(paste0('CIRI/', ciri_file), sep='\t',header=T) %>%
    select(circRNA_ID, p.values)
})

names(ciri)<-gsub('.CIRI.merged','',ciri_files)
ciri_df<-ldply(ciri,.id='Somatic') %>%
  mutate(p.values=1) %>%
  #rename(Identifier = circRNA_ID) %>%
  spread(Somatic, p.values, fill=0)
  

caclIntersect<-function(n){
  sum(sapply(1:n, function(i){choose(n, i)}))
}

caclIntersect2<-function(n){
  2^n-1
}

mutations <- read.csv( system.file("extdata", "mutations.csv", package = "UpSetR"), header=T, sep = ",")
upset(mutations, sets = c("PTEN", "TP53", "EGFR", "PIK3R1", "RB1"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")

upset(ciri_df, nsets=10, 
      #nintersects=caclIntersect2(length(ciri_df)-1),
      order.by='freq',
      decreasing=c(T,F)
      )


ciri_df$count<-rowSums(ciri_df[,2:11])

ciri_df$count %>%
  table %>%
  data.frame %>%
  ggplot(aes(x=., y=Freq)) +
  geom_bar(stat='identity')

