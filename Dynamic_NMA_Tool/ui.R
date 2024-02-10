shinyUI(fixedPage(
  #fluidPage(
  title = "Dynamic NMA Tool",
  
  shinyjs::useShinyjs(),
  
  tags$head(tags$style(
    "body { word-wrap: break-word; }")),
  
  tags$style(HTML("
                  .tabbable > .nav > li > a                  
                  {background-color: black;  color:white}"
  )),
  
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  ),
  
  
  
  
  fluidRow(
    column(width = 12, titlePanel(div(img(src="KI_Logo.png", width = 1200)))#,

           )),
  
  
  tabsetPanel(
    
    tabPanel("Outcome Selection",
             sidebarLayout(
               sidebarPanel(width=3,
                            style = "background-color: lightgrey",
                            
                            selectInput("Life_Stage",
                                        "Select life stage:",
                                        c("Pregnancy",
                                          "Exclusive Breastfeeding",
                                          "Complementary Feeding")),
                            
                            uiOutput("outcome_select")
               ),
               
               mainPanel(
                 fluidRow(
                   column(width=8,
                          fluidRow(
                            #column(1),
                            column(6,
                                   br(),
                                   uiOutput("show_properties"),
                                   uiOutput("model_props")
                            ),
                            column(3)
                          ),
                          fluidRow(
                            br(),
                            br(),
                            uiOutput("choose_layout"),
                            br(),
                            uiOutput("Network")
                          )
                   ),
                   column(width=4,
                          br(),
                          DT::dataTableOutput('mytable'),
                          br(),
                          br(),
                          uiOutput("show_weights"),
                          uiOutput("Links_Weights")
                   )
                   
                 )
               )
               
               
             )
    ),
    
    
    
    
    
    
    tabPanel("Evidence Base",
             mainPanel(
               br(),
               br(),
               fluidRow(
                 column(1),
                 column(8,
                        uiOutput("Study_tab_header"),
                        br(),
                        DT::dataTableOutput("contents")
                 )
               )
               
             )
    ),
    
    
    
    
    
    
    
    
    
    
    tabPanel("Forest Plots",
             sidebarLayout(
               sidebarPanel(width=3,
                            style = "background-color: lightgrey",
                            uiOutput("forest_cover"),
                            br(),
                            uiOutput("Forest_OR_RR"),
                            uiOutput("baseline_trt_select"),
                            uiOutput('all_vs_one_select'),
                            uiOutput("forest_trt_select")
               ),
               mainPanel(width = 9,
                 br(),
                 uiOutput("Forest_tab_header"),
                 br(),
                 uiOutput("Forest_Main")
               )
             )
    ),
    
    
    
    
    
    
    
    tabPanel("MCID",
             sidebarLayout(
               sidebarPanel(width=3,
                            style = "background-color: lightgrey",
                            uiOutput("MCID_select"),
                            br(),
                            uiOutput("MCID_highlight"),
                            br(),
                            tableOutput("trt_lst_renumb_MCID")
               ),
               mainPanel(
                 br(),
                 uiOutput("MCID_tab_header"),
                 br(),
                 uiOutput("display_table_select"),
                 fluidRow(
                   class = "rowhide",
                   column(width = 12,
                          formattable::formattableOutput("MCID_Table"),
                          tags$head(tags$style(type = "text/css", "#MCID_Table th {display:none;}"))
                   )
                 ),
                 fluidRow(
                   column(width=4,
                          br(),
                          uiOutput("trt_select_MCID")
                   ),
                   column(width=4,
                          br(),
                          uiOutput("MCID_Plot_Select")
                   ),
                   column(width=4,
                          br(),
                          uiOutput("MCID_base_vs_all_select")
                   )
                 ),
                 fluidRow(
                   column(width=12,
                          br(),
                          br(),
                          uiOutput("MCID_plot")
                   )
                 )
               )
             )
    ),
    
    
    
    
    
    
    tabPanel("Ranking",
             sidebarLayout(
               sidebarPanel(width=2,
                            style = "background-color: lightgrey",
                            uiOutput("rank_stats_choice")
               ),
               mainPanel(width = 10,
                 br(),
                 uiOutput("Ranking_tab_header"),
                 br(),
                 uiOutput("display_SUCRA_plot"),
                 fluidRow(class = "rowhide2",
                   br(),
                   column(width=12,
                          uiOutput("rank_plot")
                   )
                 ),
                 fluidRow(
                   br(),
                   uiOutput("sucra_bars_or_table"),

                   # column(width=5,
                   #        br(),
                   #        tableOutput("rank_tab_out")
                   # ),
                   # column(width=1),
                   # column(width=6,
                   #        br(),
                   #        uiOutput("SUCRA_barplot")
                   # )
                   column(12,
                     uiOutput("sucra_table_barplot")
                   )
                 )
               )
             )
    ),
    
    
    
    
    
    
    tabPanel("Cross Tables",
             sidebarLayout(
               sidebarPanel(width=3,
                            style = "background-color: lightgrey",
                            uiOutput("CrossTab_cover"),
                            br(),
                            uiOutput("CrossTabs_OR_RR"),
                            br(),
                            uiOutput("bottom_trt_select"),
                            br(),
                            uiOutput("trt_boxes")
               ),
               mainPanel(
                 fluidRow(
                   br(),
                   uiOutput("Crosstab_tab_header"),
                   br(),
                   column(12,
                          uiOutput("CrossTabs_Main"),
                          tags$head(tags$style(type = "text/css", "#CrossTabs_Main th {display:none;}"))
                   )
                 ),
                 fluidRow(
                   column(4,
                          br(),
                          br(),
                          DT::dataTableOutput("trt_lst_renumbered")
                   )
                 )
               )
             )
    ),
    
    
    
    
    
    
    tabPanel("Leverage Plot",
             fluidRow(width = 12,
               column(10,
                      fluidRow(
                        column(1),
                        column(10,
                               br(),
                               uiOutput("Leverage_tab_header"),
                               br(),
                               #plotOutput("Resids")
                               plotly::plotlyOutput("Resids_Plotly")
                        ),
                        column(1)
                      ),
                      fluidRow(
                        column(2),
                        column(7,
                               br(),
                               br(),
                               uiOutput("show_resids"),
                               uiOutput("RD_table")
                        ),
                        column(1)
                      )
               ),
               column(1)
               # ,
               # column(4,
               #        br(),
               #        br(),
               #        htmlOutput("extr_outlier_title"),
               #        br(),
               #        tableOutput("Extreme_Outliers"),
               #        br(),
               #        htmlOutput("outlier_title"),
               #        br(),
               #        tableOutput("Outliers")
               # )
             )
    ),
    
    
    
    
    
    
    
    tabPanel("About",
             mainPanel(
               fluidRow(
                 column(width=1),
                 column(width=10,
                        uiOutput("about")
                 )
               )
             )
    )
    
    
    
    
    
    
    
    
  )
  )
  )



