library(shiny)
library(DT)

collapsibleDiv<-function(id, ..., label='Show/Hide', .func=actionButton,
                         collapse = FALSE, class=NULL, icon=NULL, width=NULL){
  
  collapse_status<-ifelse(collapse, "on", "in")
  
  list(
    .func(
      sprintf("b_%s",id), label=label, icon=icon, class=class, width=width,
      "data-toggle"='collapse', "data-target"=sprintf('#%s',id)
    ),
    div(
      id = id, class = sprintf("collapse %s", collapse_status),
      ...
    )
  )
}

shinyUI(fluidPage(
  titlePanel("circRNA BAM Viewer"),
  collapsibleDiv(id='info', collapse = T,
                 label = 'Settings',
                 class = 'btn-info btn-xs',
                 icon = icon('info-sign',lib='glyphicon'),
                 
                 column(3,
                        textInput("norm_ciri", "Normal CIRI Path",
                                  value='CIRI/S33N.CIRI'),
                        textInput("tumor_ciri", "Tumor CIRI Path",
                                  value='CIRI/S33T.CIRI'),
                        textInput("norm_bam", "Normal BAM Path",
                                  value='bam/S33N/S33N.sort.cgc.bam'),
                        textInput("tumor_bam", "Tumor CIRI Path",
                                  value='bam/S33T/S33T.sort.cgc.bam'),
                        textOutput('helper')
                        ),
                 column(3
                        
                        )
                 ),
  collapsibleDiv(id='qc', collapse = T,
                 label = 'QC',
                 class = 'btn-info btn-xs',
                 uiOutput('qc_column_ctrl'),
                 plotOutput('qc')
  ),
  collapsibleDiv(id='normal_tumor_ciri', collapse = F,
                 label = 'Normal-Tumor CIRI',
                 class = 'btn-info btn-xs',
                 DT::dataTableOutput('norm_tumor_ciri')
  ),
  collapsibleDiv(id='track', collapse = T,
                 label = 'plot',
                 class = 'btn-info btn-xs',
                 plotOutput('track', height='800px')
  ),
  collapsibleDiv(id='normal_junction_reads', collapse = F,
                 label = 'Normal junction reads',
                 class = 'btn-info btn-xs',
                 checkboxInput('norm_junction_only', 'junction only', T),
                 DT::dataTableOutput('norm_reads_tb'),
                 textOutput('norm_reads_tb_uniq_qname')
  ),
  collapsibleDiv(id='tumor_junction_reads', collapse = T,
                 label = 'Tumor junction reads',
                 class = 'btn-info btn-xs',
                 checkboxInput('tumor_junction_only', 'junction only', T),
                 DT::dataTableOutput('tumor_reads_tb')
  ),
  collapsibleDiv(id='gene_reads', collapse = T,
                 label = 'Gene region reads',
                 class = 'btn-info btn-xs',
                 DT::dataTableOutput('gene_reads_tb')
  ),
  collapsibleDiv(id='reads_plot', collapse = T,
                 label = 'reads viz',
                 class = 'btn-info btn-xs',
                 plotOutput('reads_track', height='800px')
  ),
  collapsibleDiv(id='go_to_region', collapse = T,
                 label = 'Bam Navigation',
                 class = 'btn-info btn-xs',
                 DT::dataTableOutput('nav_reads_tb')
  )
))
