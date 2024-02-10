options(shiny.maxRequestSize = 5000*1024^2) 


source("NMA Functions Baseline Risk.R")



shinyServer(
  function(input, output, session){
    
    
    
    
    
    my_Life_Stage = reactive({
      input$Life_Stage
    })
    
    
    
    

    output$outcome_select = renderUI({
      if(my_Life_Stage() == "Pregnancy"){
        out_list = c("Low Birth Weight",
                     "Mean Birth Weight",
                     "Preterm Birth")
      } else if(my_Life_Stage() == "Exclusive Breastfeeding"){
        out_list = c("Change in LAZ", "Stunting")
      } else if(my_Life_Stage() == "Complementary Feeding"){
        out_list = c("Change in HAZ", "Stunting")
      }
      selectInput("Outcome", "Select outcome:", out_list)
    })
    
    
    
    
    
  
    my_input_outcome = reactive({
      input$Outcome
    })
    
    
    
    
    

    my_outcome = reactive({
      switch(my_input_outcome(),
             'Low Birth Weight' = 'LBW',
             'Mean Birth Weight' = 'MBW',
             'Preterm Birth' = 'Preterm',
             'Change in LAZ' = 'LAZ',
             'Change in HAZ' = 'HAZ',
             'Stunting' = ifelse(my_Life_Stage() == "Complementary Feeding", 
                                 'Stunting',
                                 '0-6_Stunting')
      )
    })
    
    
    
    
    
    my_filename = reactive({
      paste0("./RData Files/", my_outcome(), '_NMA_Object.RData')
    })
    
    
    
    
    
    my_NMA = reactive({
      env = new.env()
      load(my_filename(), env)
      env[[names(env)]]
    })
    
    
    
    
    
    
    treatList = reactive({
      trt_names = my_NMA()$trt_names
      treatList = as.data.frame(cbind(1:length(trt_names), 
                                      trt_names))
      names(treatList) = c("Number", "Treatment")
      
      treatList
    })
    
    
    
    
    
    
    output$mytable = DT::renderDataTable({
      data.frame(treatList())
    }, server = FALSE, escape = FALSE,  
    options = list( 
      dom = 't',
      ordering = FALSE,
      columnDefs = list(list(className = 'dt-center', targets = 0)),
      paging = FALSE,
      preDrawCallback = JS('function() { 
                           Shiny.unbindAll(this.api().table().node()); }'), 
      drawCallback = JS('function() { 
                        Shiny.bindAll(this.api().table().node()); } '),
      searching = FALSE
    ), 
    rownames = FALSE)
    
    
    
    
    
    
    output$choose_layout = renderUI({
        selectInput("layout", "Select network layout:", 
                    c("Kamada-Kawai",
                      "Circle", 
                      "Star", 
                      "Layered"),
                    width = '33%')
    })
    
    
    
    
    
    
    renumbered_t = reactive({
      temp = my_NMA()$BUGS_data
      t_inds = grep("t[[]", colnames(temp))
      t = temp[,t_inds]
      exist = sort(unique(c(t[!is.na(t)])))
      for(i in 1:length(exist)){
        for(j in 1:nrow(t)){
          for(k in 1:ncol(t)){
            if(!is.na(t[j,k]) & t[j,k]==exist[i]) t[j,k] = i 
          }
        }
      }
      
      temp[,t_inds] = t
      renumbered_t = temp
    })
    
    
    
    
    
    output$Network = renderPlot({
      if(input$layout == "Circle"){
        lay = "circ"
      } else if(input$layout == "Star"){
        lay = "star"
      } else if(input$layout == "Layered"){
        lay = "layered"
      } else if(input$layout == "Kamada-Kawai"){
        lay = "kk"
      }
      NetworkPlot(renumbered_t(), lay=lay)
    })
    
    
    
    
    
    my_Links = reactive({
      NetworkPlot(renumbered_t(), Plot = FALSE)
    })
    
    
    
    
    
    
    
    output$show_weights = renderUI({
      checkboxInput('Show_Weights',
                    "Display number of comparisons")
    })
    
    
    
    
    
    
    my_show_weights = reactive({
      input$Show_Weights
    })
    
    
    
    
    
    
    output$Links_Weights = renderUI({
      if(my_show_weights()){
        output$Links_Weights2 = DT::renderDataTable({
          my_Links()},
          server = FALSE,
          escape = FALSE,
          rownames = FALSE,
          selection = 'none',
          
          options = list( 
            dom = 't',
            paging = FALSE,
            searching = FALSE,
            ordering = FALSE,
            preDrawCallback = JS('function() { 
                                 Shiny.unbindAll(this.api().table().node()); }'), 
            drawCallback = JS('function() { 
                              Shiny.bindAll(this.api().table().node()); }'),
            columnDefs = list(list(className = 'dt-center', targets = 0:2))
          )
        )
        DT::dataTableOutput("Links_Weights2")
      }
    })
    
    
    
    
    
    my_height_network = reactive({
      n = length(unique(unlist(my_Links()[,1:2])))
      500 + (n - 16)*15*(n >= 16)
    })
    
    
    
    
    output$Network = renderUI({   
      output$Network2 = renderPlot({   
        if(input$layout == "Circle"){
          lay ="circ"
        } else if(input$layout == "Star"){
          lay ="star"
        } else if(input$layout == "Layered"){
          lay ="layered"
        } else if(input$layout == "Kamada-Kawai"){
          lay ="kk"
        } 
        
        NetworkPlot_new(BUGS_data = renumbered_t(), 
                        lay = lay, 
                        Plot = TRUE,
                        vertex.color = "orange")
      }, height = my_height_network())
      plotOutput("Network2")
    })
    
    
    
    
    
    my_properties = reactive({
      my_properties = properties_table(my_NMA())[1:6,]
    })
    
    
    
    
    
    output$show_properties = renderUI({
      checkboxInput('Show_Properties',
                    "Display statistical model attributes")
    })
    
    
    
    
    
    my_show_properties = reactive({
      input$Show_Properties
    })
    
    
    
    
    
    output$model_props = renderUI({
      if(my_show_properties()){
        output$model_props2 = DT::renderDataTable({
          my_properties()},
          server = FALSE,
          escape = FALSE,
          rownames = FALSE,
          selection = 'none',
          
          options = list( 
            dom = 't',
            paging = FALSE,
            searching = FALSE,
            ordering = FALSE,
            preDrawCallback = JS('function() { 
                                 Shiny.unbindAll(this.api().table().node()); }'), 
            drawCallback = JS('function() { 
                              Shiny.bindAll(this.api().table().node()); }'),
            columnDefs = list(list(className = 'dt-center', targets = 1))
            )
          )
        
        DT::dataTableOutput("model_props2")
      }
    })
    
    
    
    
    
    contentsrea = reactive({
        df = as.data.frame(my_NMA()$BUGS_data)
    })
    
    
    
    
    
    
    my_response = reactive({
      test1_1 = grep("se[[]", colnames(contentsrea()))
      test1_2 = grep("t[[]", colnames(contentsrea()))
      test1_3 = (length(test1_1) == length(test1_2))
      test2 = grep("E[[]", colnames(contentsrea()))
      if(length(test1_1) > 0 & test1_3) {
        my_response = "normal"
      } else if(length(test1_1) > 0){
        my_response = "normal_SMD"
      } else if(length(test2) > 0){
        my_response = "poisson"
      } else{
        my_response = "binomial"
      }
    })
    
    
    
    
    
    
    
    my_response_base_risk = reactive({
      if(!is.null(input$Base_adjusted)){
        if((my_response() == "normal") & (input$Base_adjusted)){
          my_response_base_risk ="normal_base_risk"
        } else if((my_response() == "binomial") & (input$Base_adjusted)){
          my_response_base_risk ="binomial_base_risk"
        } else if((my_response() == "poisson") & (input$Base_adjusted)){
          my_response_base_risk = "poisson_base_risk"
        } else{
          my_response_base_risk = my_response()
        }
      } else{
        my_response_base_risk = my_response()
      }
    })
    
    
    
    
    
    
    ColNames = reactive({
      resType = my_response()
      ColNames = names(contentsrea())
      comp.inds = c(grep("t[[]", ColNames), grep("ID", ColNames))
      if(resType %in% c("normal", "normal_SMD")){
        comp.inds = c(comp.inds, grep("y[[]", ColNames))
        comp.inds = c(comp.inds, grep("se[[]", ColNames))
      } else if(my_response() == "binomial"){
        comp.inds = c(comp.inds, grep("r[[]", ColNames))
        comp.inds = c(comp.inds, grep("n[[]", ColNames))
      } else if(my_response() == "poisson"){
        comp.inds = c(comp.inds, grep("E[[]", ColNames))
        comp.inds = c(comp.inds, grep("r[[]", ColNames))
      }
      comp.inds = c(comp.inds, grep("N[[]", ColNames))
      ColNames = c("None", ColNames[-c(1,comp.inds)])
    })
    
    
    
    
    
    
    my_events_are_bad_front = reactive({
      my_outcome() %in% c("LBW", 
                          "Preterm", 
                          "Stunting",
                          "0-6_Stunting")
    })
    
    
    
    
    
    
    
    #**********************************************************************************  
    #**********************************************************************************
    #*********                      Goodness of Fit Tab                       *********
    #**********************************************************************************  
    #**********************************************************************************
    my_RE = reactive({
      my_NMA()$RE
    })
    
    
    
    
    
    
    my_Leverage_Plot = reactive({
      residual_plot_new(my_NMA(), my_RE(), Plot = TRUE)
    })
    
    
    
    
    
    
    my_Leverage_Plotly = reactive({
      residual_plotly(my_NMA(), my_RE(), Plot = TRUE, resid_dev_table())
    })
    
    
    
    
    
    
    
    output$Resids = renderPlot({
      my_Leverage_Plot()
    })
    
    
    
    
    
    
    output$Resids_Plotly = plotly::renderPlotly({
      my_Leverage_Plotly()
    })
    
    
    
    
    
    
    my_Resids = reactive({
      residual_plot_new(my_NMA(), my_RE())
    })
    
    
    
    
    
    
    my_trt_names = reactive({
      my_trt_names = treatList()[,2:1]
    }) 
    
    
    
    
    
    output$Outliers = renderTable({
      outl = my_Resids()$outliers
      if(is.null(outl)) stop()
      else{
        O = nrow(outl)
        trts = character(O)
        temp = renumbered_t()
        inds = grep("t[[]", colnames(temp))
        ID_ind = grep("ID", colnames(temp))
        IDnms = c()
        
        if(length(ID_ind) > 0){
          IDs = as.matrix(temp[, ID_ind])
          IDnms = colnames(temp)[ID_ind]
        }
        
        t_table = temp[, inds]
        t_table = renumbering(t_table, my_trt_names())$t
        
        for(i in 1:O){
          trts[i] = my_NMA()$trt_names[t_table[outl[i, 1], outl[i,2]]]
        }
        
        if(length(ID_ind) == 0){
          out = data.frame(cbind(outl[,1], trts))
          names(out) = c("Study No.", "Treatment")
        } else{
          row_num = outl[,1]
          IDs = apply(IDs, 2, function(x){gsub("_", " ", x)})
          if(length(row_num) == 1){
            out = data.frame(t(c(row_num, IDs[row_num,], trts)))
          } else{
            out = data.frame(cbind(row_num, IDs[row_num,], trts))
          }
          names(out) = c("Study No.", IDnms, "Treatment")
        }
        
        out
      }
    })
    
    
    
    
    
    
    output$outlier_title = renderText({
      outl = my_Resids()$outliers
      if(is.null(outl)) return(NULL)
      else{
        HTML("<font size='4'><b> Potential outliers: </b></font>")
      }
    })
    
    
    
    
    
    
    output$extr_outlier_title = renderText({
      outl = my_Resids()$extreme_outliers
      if(is.null(outl)) return(NULL)
      else{
        HTML("<font size='4'><b> Outliers: </b></font>")
      }
    })
    
    
    
    
    
    output$Extreme_Outliers = renderTable({
      outl = my_Resids()$extreme_outliers
      if(is.null(outl)) return(NULL)
      else{
        O = nrow(outl)
        trts = character(O)
        temp = renumbered_t()
        inds = grep("t[[]", colnames(temp))
        ID_ind = grep("ID", colnames(temp))
        IDnms = c()
        
        if(length(ID_ind) > 0){
          IDs = as.matrix(temp[, ID_ind])
          IDnms = colnames(temp)[ID_ind]
        }
        
        t_table = temp[, inds]
        t_table = renumbering(t_table, my_trt_names())$t
        
        for(i in 1:O){
          trts[i] = my_NMA()$trt_names[t_table[outl[i, 1], outl[i,2]]]
        }
        
        if(length(IDnms) == 0){
          out = data.frame(cbind(outl[,1], trts))
          names(out) = c("Study No.", "Treatment")
        } else{
          row_num = outl[,1]
          if(length(row_num) == 1){
            out = data.frame(t(c(row_num, IDs[row_num,], trts)))
          } else{
            out = data.frame(cbind(row_num, IDs[row_num,], trts))
          }
          names(out) = c("Study No.", IDnms, "Treatment")
        }
        
        out
      }
    })
    
    
    
    
    
    
    resid_dev_table = reactive({
      RD_table(my_NMA())
    })
    
    
    
    
    
    
    output$show_resids = renderUI({
      checkboxInput(
        'Show_Resid',
        'Display residual deviance table')
    })
    
    
    
    
    
    
    
    my_show_resid = reactive({
      input$Show_Resid
    })
    
    
    
    
    
    
    output$RD_table = renderUI({
      if(my_show_resid()){
        output$RD_table2 = DT::renderDataTable({
          x = resid_dev_table()
          m = ncol(x)
          x
        }, server = FALSE, escape = FALSE,  
        options = list( 
          dom = 't',
          ordering = FALSE,
          columnDefs = list(list(className = 'dt-right', targets = 3:4)),
          paging = FALSE,
          preDrawCallback = JS('function() { 
                           Shiny.unbindAll(this.api().table().node()); }'), 
          drawCallback = JS('function() { 
                        Shiny.bindAll(this.api().table().node()); } '),
          searching = FALSE
        ), 
        rownames = FALSE)
        
        DT::dataTableOutput("RD_table2")
      }
    })
    
  
    
    
    
    
    
    
    #**********************************************************************************  
    #**********************************************************************************
    #*********                    Cross Table Tab                           *********
    #**********************************************************************************  
    #**********************************************************************************
    
    output$CrossTab_cover = renderUI({
      numericInput("Cover", "Select interval coverage:", 
                   min = .9, step = .005, value = .95, 
                   max = .995, width = '50%')
    })
    
   
    
    
    
    
    output$CrossTabs_OR_RR = renderUI({
      if(!(my_response() %in% c("normal", "normal_SMD", "poisson"))){
        radioButtons(inputId = "Select_OR_RR", 
                     label = "Select comparative measure:",
                     c("Odds Ratio", "Relative Risk"))
      }
    }) 
    
    
    
    
    
    
    my_OR_RR_select = reactive({
      input$Select_OR_RR
    })
    
    
    
    
    
    
    output$trt_boxes = renderUI({
      nms = my_NMA()$trt_names
      m = length(nms)
      checkboxGroupInput("trts_to_view", 
                         "Select treatments:", 
                         nms,
                         selected = nms)
    })
    
    
    
    
    
    trt_nms_to_view = reactive({
      input$trts_to_view
    })
    
    
    
    
    
    my_NMA_selected_CrossTab = reactive({
      my_NMA_selected_CrossTab = my_NMA()
      if(length(trt_nms_to_view()) != length(my_NMA()$trt_names)){
        trt_nums = which(!(my_NMA()$trt_names %in% trt_nms_to_view()))
        trt_remain = which(my_NMA()$trt_names %in% trt_nms_to_view())
        new_nums = c(1:length(trt_remain))
        mat = my_NMA()$meta.sim$sims.matrix
        
        inds = c()
        for(i in trt_nums){
          test1 = sapply(colnames(mat), function(x){grepl(paste0("[", i, ",", 
                                                                 collapse=""), 
                                                          x, fixed=T)})
          test2 = sapply(colnames(mat), function(x){grepl(paste0(",", i, "]",
                                                                 collapse=""), 
                                                          x, fixed=T)})
          inds = c(inds, which(test1), which(test2))
        }
        inds = unique(inds)
        mat = mat[,-inds]
        for(i in 1:length(trt_remain)){
          colnames(mat) = gsub(paste0("[[]", trt_remain[i], ","), 
                               paste0("[", new_nums[i], ","), 
                               colnames(mat))
          colnames(mat) = gsub(paste0(",", trt_remain[i], "[]]"), 
                               paste0(",", new_nums[i], "]"), 
                               colnames(mat))
        }
        trt_names = input$trts_to_view
        my_NMA_selected_CrossTab$trt_names = trt_names
        my_NMA_selected_CrossTab$meta.sim$sims.matrix = mat
        my_NMA_selected_CrossTab$nt = length(trt_remain)
      }
      
      my_NMA_selected_CrossTab
    })
    
    
    
    
    
    
    output$bottom_trt_select = renderUI({
      trts = my_NMA_selected_CrossTab()$trt_names
      selectInput('bottom_trt',
                  'Select bottom row:',
                  rev(trts))
    })
    
    
    
    
    
    my_bottom_row = reactive({
      input$bottom_trt
    })
    
    
    
    
    
    output$trt_lst_renumbered = DT::renderDataTable({
      nms = as.character(my_NMA()$trt_names)
      trts = nms[which(nms %in% input$trts_to_view)]
      trt_tab = as.data.frame(cbind(1:length(trts), trts))
      colnames(trt_tab) = c("Number", "Treatment")
      
      trt_tab
    },
    server = FALSE,
    escape = FALSE,
    rownames = FALSE,
    selection = 'none',
    
    options = list( 
      dom = 't',
      paging = FALSE,
      ordering = FALSE,
      searching = FALSE,
      preDrawCallback = JS('function() { 
                           Shiny.unbindAll(this.api().table().node()); }'), 
      drawCallback = JS('function() { 
                        Shiny.bindAll(this.api().table().node()); }'),
      columnDefs = list(list(className = 'dt-center', targets = 0))
      ))
    
    
    
    
    
    my_Mean_Diff_CrossTab = reactive({
      tab = Diff_CrossTab(my_NMA_selected_CrossTab(), 
                          coverage = input$Cover, 
                          my_bottom_row())
      tab
    })
    
    
    
    
    
    my_RR_CrossTab = reactive({
      tab = RR_CrossTab(my_NMA_selected_CrossTab(), 
                        coverage = input$Cover, 
                        my_bottom_row())
      tab
    })
    
    
    
    
    
    my_OR_CrossTab = reactive({
      tab = OR_CrossTab(my_NMA_selected_CrossTab(), 
                        coverage = input$Cover, 
                        my_bottom_row())
      tab
    })
    
    
    
    
    
    output$CrossTabs_Main = renderUI({
      output$CrossTabs_Main2 = formattable::renderFormattable({
        if(my_response() %in% c("normal", "normal_SMD")){
          M = my_Mean_Diff_CrossTab()
        } else if(my_response() == "poisson"){
          M = my_RR_CrossTab()
        } else if(my_OR_RR_select() == "Odds Ratio"){
          M = my_OR_CrossTab()
        } else{
          M = my_RR_CrossTab()
        }
        
        bottom_ind = which(my_NMA_selected_CrossTab()$trt_names == my_bottom_row())
        diag(M) = c(c(1:nrow(M))[-bottom_ind], bottom_ind)
        M = as.data.frame(M)
        
        
        if(my_response() %in% c("normal", "normal_SMD") & !my_events_are_bad_front()){
          signif = formatter("span", 
                             style = x ~ ifelse(grepl("\\*", x) & grepl("-", x), 
                                                formattable::style(color = "red", 
                                                                   font.weight = "bold",
                                                                   "text-align" = "center"), 
                                                ifelse(grepl("\\*", x), 
                                                       formattable::style(color = "blue", 
                                                                          font.weight = "bold",
                                                                          "text-align" = "center"),
                                                       formattable::style("text-align" = "center")))
          )
        } else if(my_response() %in% c("normal", "normal_SMD") & my_events_are_bad_front()){
          signif = formatter("span", 
                             style = x ~ ifelse(grepl("\\*", x) & grepl("-", x), 
                                                formattable::style(color = "blue", 
                                                                   font.weight = "bold",
                                                                   "text-align" = "center"), 
                                                ifelse(grepl("\\*", x), 
                                                       formattable::style(color = "red", 
                                                                          font.weight = "bold",
                                                                          "text-align" = "center"),
                                                       formattable::style("text-align" = "center")))
          )
        } else if(!my_events_are_bad_front()){
          signif = formatter("span", 
                             style = x ~ ifelse(grepl("\\*", x) & grepl("[[]0\\.|\\,0\\.", x), 
                                                formattable::style(color = "red", 
                                                                   font.weight = "bold",
                                                                   "text-align" = "center"), 
                                                ifelse(grepl("\\*", x), 
                                                       formattable::style(color = "blue", 
                                                                          font.weight = "bold",
                                                                          "text-align" = "center"),
                                                       formattable::style("text-align" = "center")))
          )
        } else{
          signif = formatter("span", 
                             style = x ~ ifelse(grepl("\\*", x) & grepl("[[]0\\.|\\,0\\.", x), 
                                                formattable::style(color = "blue", 
                                                                   font.weight = "bold",
                                                                   "text-align" = "center"), 
                                                ifelse(grepl("\\*", x), 
                                                       formattable::style(color = "red", 
                                                                          font.weight = "bold",
                                                                          "text-align" = "center"),
                                                       formattable::style("text-align" = "center")))
          )
        }
        
        M = formattable(M, list(area(col = colnames(M)) ~ signif),
                        align = rep("c", ncol(M)))
        M
      })
      
      formattable::formattableOutput("CrossTabs_Main2")
    })
    
    
    
    
    
    
    #**********************************************************************************  
    #**********************************************************************************
    #*********                     Forest Plot Tab                            *********
    #**********************************************************************************  
    #**********************************************************************************
    my_height = reactive({
      n = length(my_NMA()$trt_names)
      600 + (n - 16)*50*(n >= 16)
    })
    
    
    
    
    output$forest_trt_select =renderUI({
      checkboxGroupInput(inputId="forest_trt_Input",
                         label="Select treatments:", 
                         choices=my_NMA()$trt_names,
                         selected=my_NMA()$trt_names)
    })
    
    
    
    
    
    output$forest_cover = renderUI({
      numericInput("CoverForest", "Select interval coverage:", 
                   min = .9, step = .005, value = .95, 
                   max = .995, width = '50%')
    })
    
    
    
    
    
    my_forest_cover = reactive({
      input$CoverForest
    })
    
    
    
    
    
    output$baseline_trt_select = renderUI({
      selectInput("Baseline_Input",
                  "Select baseline treatment:", 
                  input$forest_trt_Input,
                  width = '75%')
    })
    
    
    
    
    
    my_base_input = reactive({
      input$Baseline_Input
    })
    
    
    
    
    
    output$Forest_OR_RR = renderUI({
      if(!(my_response() %in% c("normal", "normal_SMD", "poisson"))){
        radioButtons(inputId = "Select_OR_RR_Forest", 
                     label = "Select comparative measure:",
                     c("Odds Ratio", "Relative Risk"))
      }
    })
    
    
    
    
    
    my_OR_RR_forest = reactive({
      input$Select_OR_RR_Forest
    })
    
    
    
    
    
    output$all_vs_one_select = renderUI({
      radioButtons('all_vs_one_input',
                   ' ',
                   c('All vs. baseline',
                     'Baseline vs. all'))
    })
    
    
    
    
    
    my_all_vs_one = reactive({
      input$all_vs_one_input == 'All vs. baseline'
    })
    
    
    
    
    
    my_NMA_selected_new = reactive({
      my_NMA_selected_new = my_NMA()
      if(length(input$forest_trt_Input) != length(my_NMA()$trt_names)){
        trt_nums = which(!(my_NMA()$trt_names %in% input$forest_trt_Input))
        trt_remain = which(my_NMA()$trt_names %in% input$forest_trt_Input)
        new_nums = c(1:length(trt_remain))
        mat = my_NMA()$meta.sim$sims.matrix
        
        inds = c()
        for(i in trt_nums){
          test1 = sapply(colnames(mat), function(x){grepl(paste0("[", i, ",", 
                                                                 collapse=""), 
                                                          x, fixed=T)})
          test2 = sapply(colnames(mat), function(x){grepl(paste0(",", i, "]",
                                                                 collapse=""), 
                                                          x, fixed=T)})
          inds = c(inds, which(test1), which(test2))
        }
        inds = unique(inds)
        mat = mat[,-inds]
        for(i in 1:length(trt_remain)){
          colnames(mat) = gsub(paste0("[[]", trt_remain[i], ","), 
                               paste0("[", new_nums[i], ","), 
                               colnames(mat))
          colnames(mat) = gsub(paste0(",", trt_remain[i], "[]]"), 
                               paste0(",", new_nums[i], "]"), 
                               colnames(mat))
        }
        trt_names = input$forest_trt_Input
        my_NMA_selected_new$trt_names = trt_names
        my_NMA_selected_new$meta.sim$sims.matrix = mat
        my_NMA_selected_new$nt = length(trt_remain)
      }
      
      my_NMA_selected_new
    })
    
    
    
    
    
    
    
    output$Forest_Main = renderUI({
      output$Forest_Main2 = renderPlot({
        if(my_response() %in% c("normal", "normal_SMD")){
          Diff_forest(my_NMA_selected_new(), 
                      baseline = which(my_NMA_selected_new()$trt_names ==
                                         my_base_input()),
                      coverage = my_forest_cover(),
                      all_vs_one = my_all_vs_one())
        } else if(my_response() == "poisson"){
          RR_forest(my_NMA_selected_new(), 
                    baseline = which(my_NMA_selected_new()$trt_names ==
                                       my_base_input()), 
                    coverage = my_forest_cover(),
                    all_vs_one = my_all_vs_one())
        } else if(is.null(my_OR_RR_forest())) stop()
        else if(my_OR_RR_forest() == "Odds Ratio"){
          OR_forest(my_NMA_selected_new(), 
                    baseline = which(my_NMA_selected_new()$trt_names ==
                                       my_base_input()), 
                    coverage = my_forest_cover(),
                    all_vs_one = my_all_vs_one())
        } else{
          RR_forest(my_NMA_selected_new(), 
                    baseline = which(my_NMA_selected_new()$trt_names == 
                                       my_base_input()), 
                    coverage = my_forest_cover(),
                    all_vs_one = my_all_vs_one())
        }
      }, height = my_height()
      )
      
      plotOutput("Forest_Main2", width = "100%")
    })
    
    
    
    
    
    
    #**********************************************************************************  
    #**********************************************************************************
    #*********                         MCID Tab                               *********
    #**********************************************************************************  
    #**********************************************************************************
    output$MCID_select = renderUI({
      numericInput("MCID", 
                   ifelse(my_response() %in% c("binomial", "poisson"),
                          "Minimal Clinically Important Difference (RRR):",
                          "Minimal Clinically Important Difference:"),
                   value = 0.15, 
                   min = 0, max = 1,
                   step = 0.05,
                   width = '100%'
      )
    })
    
    
    
    
    
    my_MCID = reactive({
      input$MCID
    })
    
    
    
    
    
    
    my_MCID_table = reactive({
      pairwise_superiority_table(my_NMA(), 
                                 my_MCID(), 
                                 my_response(), 
                                 my_events_are_bad_front())
    })
    
    
    
    
    
    
    output$display_table_select = renderUI({
      if(!is.null(my_MCID_table())){
        checkboxInput("is_MCID_tab_disp",
                      "Display MCID table")
      }
    })
    
    
    
    
    
    
    my_MCID_tab_display = reactive({
      input$is_MCID_tab_disp
    })
    
    
    
    
    
    
    
    output$MCID_highlight = renderUI({
      if(my_MCID_tab_display()){
        numericInput("MCID_high", 
                     "Highlight probabilities greater than -",
                     value = .8, min = 0, 
                     max = 1, step = 0.05,
                     width = '100%')
      }
    })
    
    
    
    
    
    
    my_MCID_thold = reactive({
      input$MCID_high
    })
    
    
    
    
    
    output$trt_lst_renumb_MCID = renderTable({
      if(my_MCID_tab_display()){
        trts = as.character(my_NMA()$trt_names)
        trt_tab = as.data.frame(cbind(1:length(trts), trts))
        colnames(trt_tab) = c("Number", "Treatment")
        
        trt_tab
      }
    })
    
    
    
    
    
    
    output$MCID_Table = renderFormattable({
      if(my_MCID_tab_display()){
        M = my_MCID_table()
        if(anyNA(M)) stop()
        
        m = ncol(M)
        M = apply(M, 2, function(x){paste0(sprintf("%.1f", x*100), "%")}) 
        diag(M) = 1:m
        M = M %>% as.data.frame
        
        
        signif = formatter("span", 
                           style = x ~ ifelse(!grepl("%", x),
                                              formattable::style(color = "black", 
                                                                 font.weight = "bold"),
                                              ifelse(as.numeric(substr(x, 1, nchar(x)-1)) >= 100*input$MCID_high, 
                                                     formattable::style(color = "blue", 
                                                                        font.weight = "bold"),
                                                     NA))
        )
        
        
        M = formattable(M, list(area(col = colnames(M)) ~ signif),
                        align = rep("r", ncol(M)))
        M
      }
    })
    
    
    
    
    
    output$MCID_Plot_Select = renderUI({
      selectInput("MCID_Plot_Type",
                  "Select plot type:",
                  c("Bar plot", "Density plots"))
    })
    
    
    
    
    
    
    my_MCID_Plot_Select = reactive({
      input$MCID_Plot_Type
    })
    
    
    
    
    
    output$trt_select_MCID = renderUI({
      selectInput("MCID_baseline",
                  "Select baseline treatment:",
                  my_NMA()$trt_names)
    })
    
    
    
    
    
    
    my_MCID_baseline = reactive({
      input$MCID_baseline
    })
    
    
    
    
    
    
    output$MCID_base_vs_all_select = renderUI({
      radioButtons("MCID_base_vs_all",
                   '',
                   c("All vs. baseline",
                     "Baseline vs. all")
                   )
    })
    
    
    
    
    
    
    MCID_is_base_vs_all = reactive({
      input$MCID_base_vs_all == "Baseline vs. all"
    })
    
    
    
    
    
    
    my_MCID_plot = reactive({
      MCID_Probs_Plot(my_MCID_table(), 
                      my_MCID_baseline(),
                      my_MCID(),
                      MCID_is_base_vs_all())
    })
    
    
    
    
    
    output$MCID_Plot_Select = renderUI({
      selectInput("MCID_Plot_Type",
                  "Select plot type:",
                  c("Bar plot", "Density plots"))
    })
    
    
    
    
    
    
    my_MCID_plot_type = reactive({
      input$MCID_Plot_Type
    })
    
    
    
    
    
    post_RRR_samples = reactive({
      post_superiority_samples(my_NMA(), 
                               my_MCID_baseline(), 
                               my_events_are_bad_front(),
                               MCID_is_base_vs_all())
    })
    
    
    
    
    
    
    output$MCID_plot = renderUI({
      output$MCID_plot2 = renderPlot({
        if(my_MCID_plot_type() == "Bar plot"){
          MCID_Probs_Plot(my_MCID_table(), 
                          my_MCID_baseline(),
                          my_MCID(),
                          MCID_is_base_vs_all())
        } else{
          MCID_Density_Plots(my_NMA(), 
                             post_RRR_samples(), 
                             my_MCID(),
                             my_MCID_baseline(),
                             MCID_is_base_vs_all())
        }
      }, height = my_height())
      plotOutput("MCID_plot2")
    })
    
    
    
    
    
    
    #**********************************************************************************  
    #**********************************************************************************
    #*********                        Ranking Tab                             *********
    #**********************************************************************************  
    #**********************************************************************************   
    my_height_rank = reactive({
      n = length(unique(unlist(my_Links()[,1:2])))
      350 + (n - 16)*15*(n >= 16)
    })
    
    
    
    
    output$rank_stats_choice = renderUI({
      selectInput("rankSelect", "Statistics to display:", 
                  choices = c("SUCRA", "Rankogram"), 
                  width = '75%')
    })
    
    
    
    
    my_rank_mat_front = reactive({
      ranks(my_NMA(), my_events_are_bad_front())
    })
    
    
    
    
    
    my_rankogram_front = reactive({
      rankogram(my_rank_mat_front(), my_NMA())
    })
    
    
    
    
    
    my_is_SUCRA = reactive({
      my_is_SUCRA = input$rankSelect
    })
    
    
    
    
    
    output$rank_plot2 = renderPlot({
      if(my_is_SUCRA() == "SUCRA"){
        SUCRA_plots(my_rankogram_front())
      } else{
        rankogram_plot(my_rankogram_front())
      }
    })
    
    
    
    
    
    
    output$display_SUCRA_plot = renderUI({
      if(my_is_SUCRA() == "SUCRA" && !is.null(my_rankogram_front())){
        checkboxInput("SUCRA_disp_select",
                      "Display SUCRA plots")
      }
    })
    
    
    
    
    
    my_display_SUCRA_plots = reactive({
      (my_is_SUCRA() != "SUCRA") || input$SUCRA_disp_select
    })
    
    
    
    
    
    output$rank_plot = renderUI({
      plotOutput("rank_plot2", 
                 height = my_height_rank())
    })
    
    
    
    
    
    output$rank_tab_out = function(){
      if(my_is_SUCRA() == "SUCRA"){
        tab = SUCRA_table(my_rank_mat_front(), my_NMA()) 
        row.names(tab) = c()
        
        knitr::kable(tab, "html", align = "clr", escape = F) %>%
          kableExtra::kable_styling("striped", full_width = F)
      } else{
        tab = rankogram(my_rank_mat_front(), my_NMA())
        tab_new = t(tab)
        tab_new = tab_new[-1,]
        v = unname(names(tab)[-1])
        tab_new = as.data.frame(cbind(1:dim(tab_new)[1], 
                                      v, tab_new))
        names(tab_new) = c("Trt. No.", 
                           "Trt. Name", 
                           1:(ncol(tab_new) - 3), "Median")
        tab_new = tab_new[order(as.numeric(as.character(tab_new$Median))),]
        row.names(tab_new) = c()
        
        knitr::kable(tab_new, "html",
                     align = paste(c("c", 
                                     "l", rep("r", ncol(tab_new) - 3), "c"))
        ) %>%
          kableExtra::kable_styling("striped", full_width = F) %>%
          add_header_above(c("", "", "Rank" = ncol(tab_new) - 3, ""))
      }
    }
    
    
    
    
    
    
    my_SUCRA_Table = reactive({
      SUCRA_table(my_rank_mat_front(), my_NMA())
    })
    
    
    
    
    
    my_Rankogram = reactive({
      rankogram(my_rank_mat_front(), my_NMA())
    })
    
    
    
    
    
    my_SUCRA_Plots = reactive({
      SUCRA_plots(my_rankogram_front())
    })
    
    
    
    
    
    my_Rankogram_Plot = reactive({
      rankogram_plot(my_rankogram_front())
    })
    
    
    
    
    
    my_SUCRA_Bar_Plot = reactive({
      SUCRA_bars(my_rank_mat_front(), my_NMA())
    })
    
    
    
    
    
    output$sucra_bars_or_table = renderUI({
      if(my_is_SUCRA() == "SUCRA"){
        radioButtons("sucra_graph_or_table",
                     "",
                     c("Display table", 
                       "Display bar plot"))
      }
    })
    
    
    
    
    
    output$SUCRA_barplot = renderUI({
      output$SUCRA_barplot2 = renderPlot({
        if(my_is_SUCRA() == "SUCRA"){
          my_SUCRA_Bar_Plot()
        } else stop()
      }, height = my_height()
      )
      
      plotOutput("SUCRA_barplot2")
    })
    
    
    
    
    output$sucra_table_barplot = renderUI({
      if(my_is_SUCRA() == "SUCRA" && 
         input$sucra_graph_or_table == "Display bar plot"){
        uiOutput("SUCRA_barplot")
      } else{
        uiOutput("rank_tab_out")
      } 
    })
    
    
    
    
    
    
    
    #**********************************************************************************  
    #**********************************************************************************
    #*********                            Data Tab                            *********
    #**********************************************************************************  
    #**********************************************************************************    
    contents_construct = reactive({
      response = my_response()
      if(is.null(response)){ 
        return()
      } else{
        temp = renumbered_t()
        t_inds = grep("t[[]", colnames(temp))
        N_inds = grep("N[[]", colnames(temp))
        names(temp)[t_inds] = paste("Trt. arm", 1:length(t_inds))
        names(temp)[N_inds] = paste("Sample size arm", 1:length(t_inds))
        if(response == "binomial"){
          n_inds = grep("n[[]", colnames(temp))
          names(temp)[n_inds] = paste("Total counts arm", 1:length(t_inds))
          r_inds = grep("r[[]", colnames(temp))
          names(temp)[r_inds] = paste("No. of cases arm", 1:length(t_inds))
        } else if(response == "normal"){
          y_inds = grep("y[[]", colnames(temp))
          if(length(y_inds) == 0) return()
          names(temp)[y_inds] = paste("Mean arm", 1:length(t_inds))
          se_inds = grep("se[[]", colnames(temp))
          names(temp)[se_inds] = paste("Std. Err. arm", 1:length(t_inds))
        } else if(response == "normal_SMD"){
          y_inds = grep("y[[]", colnames(temp))
          if(length(y_inds) == 0) return()
          names(temp)[y_inds] = paste("SMD arm", 2:(length(t_inds)))
          se_inds = grep("se[[]", colnames(temp))
          names(temp)[se_inds] = paste("SE SMD arm", 2:(length(t_inds)))
        } else if(response == "poisson"){
          E_inds = grep("E[[]", colnames(temp))
          if(length(E_inds) == 0) return()
          names(temp)[E_inds] = paste("Exposure arm", 1:length(t_inds))
          r_inds = grep("r[[]", colnames(temp))
          names(temp)[r_inds] = paste("Events arm", 1:length(t_inds))
        }
        
        ID_inds = grep("ID|Title", colnames(temp))
        if(length(ID_inds) > 0){
          temp[,ID_inds] = apply(as.matrix(temp[,ID_inds]), 2, 
                                 function(x){gsub("_", " ", x)})
        }
        
        decimal_inds = which(apply(temp[,-ID_inds], 2,
                                   function(x){
                                     grepl("\\.", paste(as.character(na.omit(x)),
                                                        collapse=''))
                                   })
        )
        
        temp[,length(ID_inds) + decimal_inds] = apply(as.matrix(temp[,length(ID_inds) + decimal_inds]), 2, 
                                                      function(x){
                                                        sprintf("%.3f", as.numeric(as.character(x)))
                                                      })
        names(temp) = gsub("_", " ", names(temp))
        temp[is.na(temp)] = "-"
        temp[temp =="NA"] = "-"
        
        temp
      }
    })
    
    
    
    
    
    
      
    output$contents = DT::renderDataTable({
      M = contents_construct()
      inds = grep("StudyID|TrialID|Title", colnames(M))
      nms = colnames(M)[inds]
      M = as.matrix(M[,inds])
      colnames(M) = nms
      M[, colnames(M) == "Title"] = enc2utf8(gsub("_", " ", M[, colnames(M) == "Title"])) 
      
      df = data.frame(
        M,
        stringsAsFactors = FALSE,
        check.names = FALSE)
      
      df
    },
    
    server = FALSE,
    escape = FALSE,
    selection = 'none',
    
    options = list( 
      dom = 't',
      paging = FALSE,
      searching = FALSE,
      preDrawCallback = JS('function() { 
                           Shiny.unbindAll(this.api().table().node()); }'), 
      drawCallback = JS('function() { 
                        Shiny.bindAll(this.api().table().node()); }'),
      columnDefs = list(list(className = 'dt-left', targets = "_all"))
    )
    )
    
    
    
    
    
    
    my_header = reactive({
      paste0("<font size='5'><b>", 
             my_Life_Stage(), ": ", my_input_outcome(),
             "</b></font>")
    })
    
    
    
    
    
    output$Study_tab_header = renderUI({
      HTML(my_header())
    })
    
    
    
    
    output$Forest_tab_header = renderUI({
      HTML(my_header())
    })
    
    
    
    
    
    output$MCID_tab_header = renderUI({
      HTML(my_header())
    })
    
    
    
    
    
    output$about = renderUI({
      withMathJax(includeMarkdown('About.Rmd'))
    })
    
    
    
    
    
    output$Crosstab_tab_header = renderUI({
      HTML(my_header())
    })
    
    
    
    
    
    output$Leverage_tab_header = renderUI({
      HTML(my_header())
    })
    
    
    
    
    output$Ranking_tab_header = renderUI({
      HTML(my_header())
    })
    
    
    
    
    
    
    observe({
      if (is.null(input$is_MCID_tab_disp) || !input$is_MCID_tab_disp) {
        shinyjs::hide(selector = ".rowhide")
      } else {
        shinyjs::show(selector = ".rowhide")
      }
    })
    
    
    
    
    
    
    observe({
      if (is.null(input$SUCRA_disp_select) ||
          (!is.null(my_display_SUCRA_plots()) && 
          !my_display_SUCRA_plots())) {
        shinyjs::hide(selector = ".rowhide2")
      } else {
        shinyjs::show(selector = ".rowhide2")
      }
    })
    
    
    
    
    
    
      }
    )

