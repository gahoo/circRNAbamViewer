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
                 checkboxInput('byGene', 'By Gene', F),
                 plotOutput('track', height='800px')
  ),
  collapsibleDiv(id='normal_junction_reads', collapse = T,
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
                 DT::dataTableOutput('tumor_reads_tb'),
                 textOutput('tumor_reads_tb_uniq_qname')
  ),
  collapsibleDiv(id='reads_plot', collapse = T,
                 label = 'reads viz',
                 class = 'btn-info btn-xs',
                 collapsibleDiv(
                   id='reads_plot_settings', collapse = F,
                   label = 'settings',
                   class = 'btn-info btn-xs  pull-right',
                   checkboxInput('reads_track_follow_reads', 'follow reads', T),
                   checkboxInput('reads_track_show_arc', 'show arc', F),
                   checkboxInput('reads_track_show_transcript', 'show transcript', F),
                   checkboxInput('reads_track_show_repeats', 'show repeats', F),
                   conditionalPanel('input.reads_track_show_repeats == true',
                                    checkboxInput('reads_track_hide_uniq_repeats', 'hide unique and middle', F),
                                    textInput('reads_track_extend_bp',
                                              'Extend:', value=1000)
                                    ),
                   checkboxInput('reads_track_show_break_start', 'show break start', F),
                   checkboxInput('reads_track_show_break_end', 'show break end', F),
                   checkboxInput('reads_track_facet', 'facet by sample', F),
                   checkboxInput('reads_track_cur_page_reads', 'plot current page reads', F),
                   uiOutput('reads_track_ctrls')
                   ),
                 plotOutput('reads_track')
  ),
  collapsibleDiv(id='go_to_region', collapse = F,
                 label = 'Bam Navigation',
                 class = 'btn-info btn-xs',
                 textInput('goto', 'Go to:'),
                 DT::dataTableOutput('nav_reads_tb')
  )
))
