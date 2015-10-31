library(shiny)
library(DT)
library(magrittr)
library(ggbio)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(dplyr)

source('prepare.R')

shinyServer(function(input, output, session) {
  
  norm_ciri<-reactive({
    message('load norm_ciri')
    read.table(input$norm_ciri, sep='\t', comment.char='', header=T, as.is=T)
  })
  
  norm_bam<-reactive({
    BamFile(input$norm_bam)
  })
  
  norm_bam_circ_region<-reactive({
    selected<-ciri_selected_row()
    what <- c("qname", "flag", "mapq")
    which <- GRanges(selected$chr, IRanges(selected$circRNA_start, selected$circRNA_end))
    # if(input$byGene){
    #   which<-genes(txdb, vals=list(gene_id=selected$gene_id))
    # }
    param <- ScanBamParam(which=which, what=what, tag=c('XA','SA', 'RG'))
    readGAlignments(norm_bam(), param=param)
  })
  
  norm_bam_circ_region_junction<-reactive({
    selected<-ciri_selected_row()
    bam_region<-norm_bam_circ_region()
    junction_reads_ID<-unlist(strsplit(selected$Normal.junction_reads_ID, split=','))
    idx<-mcols(bam_region)$qname %in% junction_reads_ID
    bam_region[idx]
  })
  
  output$norm_ciri<-DT::renderDataTable({
    norm_ciri() %>%
      datatable(filter = 'top')
  })
  
  norm_reads<-reactive({
    selected<-ciri_selected_row()
    if(input$norm_junction_only){
      reads<-norm_bam_circ_region_junction()
    }else{
      reads<-norm_bam_circ_region()
    }
    
    selectedCircRNAReads(reads, selected, circRNA_ID_qnames())
  })
  
  output$norm_reads_tb<-DT::renderDataTable({
    norm_reads() %>%
      highlight_dt
      #datatable(filter = 'top')
  })
  
  output$norm_reads_tb_uniq_qname<-renderText({
    filter_values<-input$norm_reads_tb_search_columns
    uniq_qname_cnt<-countQnames(filter_values, norm_reads())
    #str(filter_values)
    #str(input$norm_reads_tb_rows_all)
    #str(input$norm_reads_tb_rows_current)    
    paste0("Unique qname counts: ", uniq_qname_cnt)

  })
  
  tumor_ciri<-reactive({
    message('load tumor_ciri')
    read.table(input$tumor_ciri, sep='\t', comment.char='', header=T, as.is=T)
  })
  
  tumor_bam<-reactive({
    BamFile(input$tumor_bam)
  })
  
  tumor_bam_circ_region<-reactive({
    selected<-ciri_selected_row()
    what <- c("qname", "flag", "mapq")
    which <- GRanges(selected$chr, IRanges(selected$circRNA_start, selected$circRNA_end))
    # if(input$byGene){
    #   which<-genes(txdb, vals=list(gene_id=selected$gene_id))
    # }
    param <- ScanBamParam(which=which, what=what, tag=c('XA','SA', 'RG'))
    readGAlignments(tumor_bam(), param=param)
  })
  
  tumor_bam_circ_region_junction<-reactive({
    selected<-ciri_selected_row()
    bam_region<-tumor_bam_circ_region()
    junction_reads_ID<-unlist(strsplit(selected$Tumor.junction_reads_ID, split=','))
    idx<-mcols(bam_region)$qname %in% junction_reads_ID
    bam_region[idx]
  })
  
  output$tumor_ciri<-DT::renderDataTable({
    tumor_ciri() %>%
      datatable(filter = 'top')
  })
  
  tumor_reads<-reactive({
    selected<-ciri_selected_row()
    if(input$tumor_junction_only){
      reads<-tumor_bam_circ_region_junction()
    }else{
      reads<-tumor_bam_circ_region()
    }
    
    selectedCircRNAReads(reads, selected, circRNA_ID_qnames())
  })
  
  output$tumor_reads_tb<-DT::renderDataTable({
    tumor_reads() %>%
      highlight_dt
      #datatable(filter = 'top')
  })
  
  output$tumor_reads_tb_uniq_qname<-renderText({
    filter_values<-input$tumor_reads_tb_search_columns
    uniq_qname_cnt<-countQnames(filter_values, tumor_reads())
    paste0("Unique qname counts: ", uniq_qname_cnt)
    
  })
  
  norm_tumor_ciri<-reactive({
    ciri<-merge(norm_ciri(), tumor_ciri(),
                by=c('circRNA_ID', 'chr', 'circRNA_start', 'circRNA_end', 'gene_id', 'circRNA_type'), all=T)
    colnames(ciri)[7:16] <- paste(rep(c('Normal', 'Tumor'),each=5),
                                  c('junction_reads', 'SM_MS_SMS', 'non_junction_reads', 'junction_reads_ratio', 'junction_reads_ID'), sep='.')
    ciri %>%
      mutate(
        transcripts = gene_id,
        gene_id =  gsub(':.*$','',gene_id),
        gene_id =  gsub('^n/a$','',gene_id),
        length = circRNA_end - circRNA_start,
        log2RatioRatio = log2(Tumor.junction_reads_ratio/Normal.junction_reads_ratio),
        log2RatioRatio = ifelse(is.infinite(log2RatioRatio), NA, log2RatioRatio),
        inNormal = !is.na(Normal.junction_reads),
        inTumor = !is.na(Tumor.junction_reads),
        inBoth = inNormal & inTumor
      ) %>%
      left_join(gene_symbol, by='gene_id') %>%
      mutate(inCGC = symbol %in% cgc_symbols) %>%
      filter(inCGC == T)
      #filter(circRNA_start>=1623888 & circRNA_end<=3764177)
  })
  
  output$norm_tumor_ciri<-DT::renderDataTable({
    norm_tumor_ciri() %>%
      select(-Normal.junction_reads_ID, -Tumor.junction_reads_ID) %>%
      datatable(filter = 'top', #selection = 'single',
                extensions = 'ColVis', options = list(dom = 'C<"clear">lfrtip'))
  })
  
  circRNA_ID_qnames<-reactive({
    selected<-ciri_selected_row()
    
    norm_qnames<-with(
      selected, 
      splitQnames(circRNA_ID, Normal.junction_reads_ID, type='Normal')
    )
    
    tumor_qnames<-with(
      selected, 
      splitQnames(circRNA_ID, Tumor.junction_reads_ID, type='Tumor')
    )
    
    rbind(norm_qnames, tumor_qnames)
  })

  circ_arc<-reactive({
    selected<-ciri_selected_row()
    
    norm_mapq <- norm_bam_circ_region_junction() %>% getQnameMapq
    tumor_mapq <- tumor_bam_circ_region_junction() %>% getQnameMapq
    
    qname_mapq<-rbind(norm_mapq, tumor_mapq)
    circRNA_ID_mapq<-merge(circRNA_ID_qnames(), qname_mapq, by='qname') %>%
      group_by(circRNA_ID, type) %>%
      summarise(mapq=mean(mapq,na.rm=T))
    
    norm_tumor_rbind() %>%
      filter(circRNA_ID %in% selected$circRNA_ID) %>%
      left_join(circRNA_ID_mapq, by=c('circRNA_ID', 'type')) ->
      selected_norm_tumor
    
    with(selected_norm_tumor,
         GRanges(seqnames = chr,
            IRanges(start = circRNA_start,
                    end = circRNA_end),
            strand = "*",
            mapq = mapq,
            junction_reads = X.junction_reads,
            junction_ratio = junction_reads_ratio,
            type = type,
            element_type = circRNA_type)
         )
  })
  
  ciri_selected_row<-reactive({
    #row_id<-input$norm_tumor_ciri_row_last_clicked
    row_id<-input$norm_tumor_ciri_rows_selected
    if(is.null(row_id)){
      row_id<-1
    }
    norm_tumor_ciri()[row_id,] %>% as.list
  })
  
  output$helper<-renderText({
    message(input$norm_tumor_ciri_row_last_clicked)
    selected<-ciri_selected_row()
    str(selected)
    message(length(norm_bam_circ_region()))
    #str(circ_arc())
    str(as.data.frame(norm_bam_circ_region_junction()))
    input$norm_tumor_ciri_row_last_clicked
  })
  
  norm_tumor_rbind<-reactive({
    norm_ciri<-norm_ciri()
    norm_ciri$type<-"Normal"
    tumor_ciri<-tumor_ciri()
    tumor_ciri$type<-"Tumor"
    rbind(norm_ciri, tumor_ciri) %>%
      mutate(
        transcripts = gene_id,
        gene_id =  gsub(':.*$','',gene_id),
        length = circRNA_end - circRNA_start
      ) %>%
      left_join(gene_symbol, by='gene_id') %>%
      mutate(inCGC = symbol %in% cgc_symbols)
  })

  output$qc<-renderPlot({
    norm_tumor<-norm_tumor_rbind()
    data_type<-class(norm_tumor[[input$qc_column]])
    if(data_type == 'character' || data_type == 'factor'){
      norm_tumor %>%
        group_by_('type', input$qc_column) %>%
        summarise(count=n()) ->
        norm_tumor_summary
      ggplot(norm_tumor_summary) +
        aes_string(x=input$qc_column, y='count', fill='type') +
        geom_bar(stat='identity', position='dodge') +
        coord_flip()
    }else{
      p<-ggplot(norm_tumor) + aes_string(x=input$qc_column, fill='type')
      p_hist<- p + geom_histogram(position='dodge') + 
        labs(title=sprintf('%s histogram', input$qc_column),
             xlab=sprintf('%s number', input$qc_column))
      p_boxplot<- p + aes_string(x='type', y=input$qc_column) +
        geom_boxplot() + 
        labs(title=sprintf('%s boxplot', input$qc_column),
             xlab=sprintf('%s number', input$qc_column))
      
      if(input$qc_xlog){
        p_hist<-p_hist+scale_x_log10()
      }
      if(input$qc_ylog){
        p_hist<-p_hist+scale_y_log10()
        p_boxplot<-p_boxplot+scale_y_log10()
      }
      if(input$qc_facet){
        p_hist<-p_hist+facet_grid(type ~ .)
      }
      multiplot(p_hist, p_boxplot)
    }

  })

  output$qc_column_ctrl<-renderUI({
    list(selectInput("qc_column", 'columns to QC',
                     choice = names(norm_tumor_rbind()),
                     selected="X.junction_reads"),
         checkboxInput('qc_xlog', 'xlog', F),
         checkboxInput('qc_ylog', 'ylog', F),
         checkboxInput('qc_facet', 'facet', T)
      )
    
  })
  
  output$track<-renderPlot({
    selected<-ciri_selected_row()
    which <- GRanges(selected$chr, IRanges(selected$circRNA_start, selected$circRNA_end))
    #which_gene<-genes(txdb, vals=list(gene_id=as.numeric(selected$gene_id)))
    region_symbols<-rmNAnullUniq(selected$symbol)
    if(length(region_symbols)==0){
      which_gene<-which
    }else{
      which_gene<-genesymbol[region_symbols]
    }
    
    region_size<-max(end(which_gene))-min(start(which_gene))
    message(selected$gene_id, ":", selected$symbol,":",length(which_gene),":",region_size)
    
    if(region_size>500000){
      cov_method<-'estimate'
    }else{
      cov_method<-'raw'
    }
    
    norm_junction_reads<-plotReads(norm_bam_circ_region_junction(), which=which)
    norm_reads<-plotReads(norm_bam_circ_region(), which=which)
    
    tumor_junction_reads<-plotReads(tumor_bam_circ_region_junction(), which=which)
    tumor_reads<-plotReads(tumor_bam_circ_region(), which=which)
    
    norm_bam_gene_cov <- plotCoverage(norm_bam(), which=which_gene, cov_method=cov_method)
    tumor_bam_gene_cov <- plotCoverage(tumor_bam(), which=which_gene, cov_method=cov_method)
    # norm_bam_circ_region_junction_reads_cov <- ggplot() + stat_coverage(granges(norm_bam_circ_region_junction_reads))
    
    #gene_model <- ggbio() + geom_alignment(data=txdb, which = which_gene, stat = "reduce")
    transcripts <- ggbio() + geom_alignment(data=txdb, which = which_gene )
    
    #arc <- ggplot(circ_arc()) + geom_arch(aes(color = type, height = junction_reads, size = mapq), alpha = 0.5)
    arc<-plotArc(circ_arc())
    track_list<-list(
      arc=arc,
      #gene_model = gene_model,
      transcripts=transcripts,
      norm=norm_reads,
      norm_junction=norm_junction_reads,
      norm_coverage=norm_bam_gene_cov,
      tumor=tumor_reads,
      tumor_junction=tumor_junction_reads,
      tumor_coverage=norm_bam_gene_cov
    )
          
    track_list<-track_list[!sapply(track_list, is.null)]
    
    title<-with(selected, paste(circRNA_ID, gene_id, symbol, sep=':'))
    tracks(track_list, title=title)
  })
  
  nav_reads<-reactive({    
    goto <- parseGoTo(input$goto)
    #str(goto)
    
    norm_reads<-loadBamReads(norm_bam(), goto)
    tumor_reads<-loadBamReads(tumor_bam(), goto)
    
    mcols(norm_reads)$type<-"Normal"
    mcols(tumor_reads)$type<-"Tumor"
    all_reads<-c(norm_reads, tumor_reads)
    mcols(all_reads)$type<-as.factor(mcols(all_reads)$type)
    all_reads
    
  })
  
  output$nav_reads_tb<-DT::renderDataTable({
    selected<-ciri_selected_row()
    circRNA_ID_qnames<-circRNA_ID_qnames()
    circRNA_ID_qnames$type<-NULL
    
    nav_reads() %>%
    selectedCircRNAReads(selected, circRNA_ID_qnames) %>%
      as.data.frame %>%
      datatable(filter = 'top',
                extensions = 'TableTools', options = list(
                  dom = 'T<"clear">lfrtip',
                  tableTools = list(sSwfPath = copySWF()))
                )
  })
  
  observe({
    selected<-ciri_selected_row()
    value<-paste(selected$circRNA_ID, collapse = ',')
    updateTextInput(session, "goto", value = sprintf("%s,",value))
  })
  
}
)
