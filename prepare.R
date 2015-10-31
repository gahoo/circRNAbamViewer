# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gene_symbol<-as.data.frame(org.Hs.egSYMBOL)
data(genesymbol, package = "biovizBase")

splitQnames<-function(circRNA_IDs, qnames, ...){
  qnames %>%
    strsplit(split=',') %>%
    lapply(function(qname){
      data.frame(qname=qname, ..., stringsAsFactors = F)
    }) %>%
    "names<-"(circRNA_IDs) %>%
    plyr::ldply(.id='circRNA_ID') %>%
    mutate(circRNA_ID = as.character(circRNA_ID) )
}

getQnameMapq<-function(reads){
  reads %>%
    mcols %>%
    as.data.frame %>%
    select(qname, mapq)
}

loadCGCsymbols<-function(){
  cgc<-read.delim('cancer_gene_census_2012_0315', sep='\t', as.is=T)
  gene_alias<-as.data.frame(org.Hs.egALIAS2EG)
  gene_alias_symbol<-merge(gene_alias, gene_symbol)
  alias<-setdiff(cgc$Symbol, gene_symbol$symbol)
  alias_idx<-gene_alias_symbol$alias_symbol %in% alias
  #gene_idx<-gene_alias_symbol$gene_id %in% cgc$GeneID
  #cgc_symbols<-c(cgc$Symbol, gene_alias_symbol[gene_idx&alias_idx,]$symbol)
  c(cgc$Symbol, gene_alias_symbol[alias_idx,]$symbol)
}

cgc_symbols<-loadCGCsymbols()

is.empty.or<-function(x, v){
  if(length(x)==0){
    v
  }else if(length(x)==1){
    ifelse(is.na(x)|is.null(x), v, x)
  }else{
    x
  }
}

is.empty<-function(x){
  if(length(x)==0){
    TRUE
  }else if(length(x)==1){
    is.na(x)|is.null(x)
  }else{
    FALSE
  }
  
}

rmNAnullUniq<-function(x){
  unique(x[!is.na(x) & !is.null(x)])
}

selectedCircRNAReads<-function(reads, selected, circRNA_ID_qnames){
  reads %>%
    as.data.frame %>%
    mutate(overlap_start = start <= selected$circRNA_start &
             end >= selected$circRNA_start,
           overlap_end = start <= selected$circRNA_end &
             end >= selected$circRNA_end) %>%
    left_join(circRNA_ID_qnames, by='qname') %>%
    mutate(not_support_start = ifelse(overlap_start&is.na(circRNA_ID), T, F),
           not_support_end = ifelse(overlap_end&is.na(circRNA_ID), T, F),
           not_support = not_support_start|not_support_end)
}

countQnames<-function(filter_values, reads){
  column_names<-colnames(reads)
  idx_ranges<-grepl('\\.\\.\\.', filter_values)
  idx_bool<-grepl('true|false', filter_values)
  idx_string<-!(idx_ranges|idx_bool)
  
  filter_start<-gsub(' .*$','',filter_values[idx_ranges])
  filter_end<-gsub('^.* ','',filter_values[idx_ranges])
  sprintf("%s >= %s & %s <= %s",
          column_names[idx_ranges],
          filter_start,
          column_names[idx_ranges],
          filter_end) %>%
    paste(collapse = ' & ') ->
    range_filters
  
  bool_filters<-grepl('true', filter_values[idx_bool])
  sprintf("%s == %s", 
          column_names[idx_bool],
          bool_filters) %>%
    paste(collapse = ' & ') ->
    bool_filters
  
  if(range_filters == '' && bool_filters != ''){
    criteria <- sprintf("~ %s", bool_filters)
  }else if(range_filters != '' && bool_filters == ''){
    criteria <- sprintf("~ %s", range_filters)
  }else if(range_filters != '' && bool_filters != ''){
    criteria <- sprintf("~ %s & %s", range_filters, bool_filters)
  }else{
    criteria <- NULL
  }
  
  if(!is.null(criteria)){
    criteria<-as.formula(criteria)
    reads <- filter_(reads, criteria)
  }
  
  reads %>%
    "$"('qname') %>%
    unique %>%
    length
}

highlight_dt<-function(reads){
  reads %>%
    datatable(filter = 'top') %>%
    formatStyle(
      'mapq',
      background = styleColorBar(reads$mapq, 'steelblue'),
      backgroundSize = '100% 90%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center'
    )
}

plotReads<-function(reads, which){
  if(length(reads)==0){
    NULL
  }else{
    autoplot(reads, which=which)
  }
}

plotArc<-function(arc){
  if(length(arc)==0){
    message("empty")
    NULL
  }else{
    message(length(arc))
    ggbio() + 
      geom_arch(data=arc, 
                aes(color = type, 
                    height = junction_reads, 
                    size = mapq),
                alpha = 0.5)
  }
}

plotCoverage<-function(bam, which, cov_method){
  ggplot() + 
    stat_coverage(bam, 
                  which=which, 
                  method=cov_method) +
    ylab(cov_method)
}