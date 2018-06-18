library(shiny)
library(scatterD3)

default_lines <- data.frame(slope = c(0, Inf), 
                            intercept = c(0, 0),
                            stroke = "#000",
                            stroke_width = 1,
                            stroke_dasharray = c(5, 5))
threshold_line <- data.frame(slope = 0, 
                             intercept = 30, 
                             stroke = "#F67E7D",
                             stroke_width = 2,
                             stroke_dasharray = "")

function(input, output, session) {
  
  data <- reactive({
    d[1:input$scatterD3_nb,]
  })
  
  lines <- reactive({
    if (input$scatterD3_threshold_line) {
      return(rbind(default_lines, threshold_line))
    }
    default_lines
  })
  
  output$scatterPlot <- renderScatterD3({
#     col_var <- if (input$scatterD3_col == "None") NULL else data()[,input$scatterD3_col]
#     symbol_var <- if (input$scatterD3_symbol == "None") NULL else data()[,input$scatterD3_symbol]
#     size_var <- if (input$scatterD3_size == "None") NULL else data()[,input$scatterD3_size]
    name_labels <- if (input$scatterD3_x_log == TRUE) treeNames else NULL
    if (input$dataSet == "veg") {
      displayData <- newwiwMDS
    }
    else {
      displayData <- newwiwRoetz
      name_labels <- if (input$scatterD3_x_log == TRUE) treeNamesRoetz else NULL
      Dtype <- DtypeRoetz
    }
    plotGrovesD3(displayData, groups=Dtype, 
                 point_size=input$scatterD3_point_size, 
                 colors=colours, col_lab="Method", 
                 point_opacity=input$scatterD3_opacity, 
                 labels_size = input$scatterD3_labsize,
                 ellipses=input$scatterD3_ellipses, 
                 ellipses_level=input$scatterD3_nb/100, 
                 treeNames=name_labels, # names to sanity check plot
                 xlab="MDS Axis 1",
                 ylab="MDS Axis 2",
                 axes_font_size="150%",
                 legend_font_size="150%",
                 transitions = TRUE,
                 click_callback = "function(id, index) {
                  if(id && typeof(Shiny) != 'undefined') {
                  Shiny.onInputChange('selected_point', index);}}")
    
  }) # endof scatterPlot
  
  observeEvent(input$reftop,{
    updateTabsetPanel(session, "mainPanel", selected = "References")
  })
  observeEvent(input$reftop2,{
    updateTabsetPanel(session, "mainPanel", selected = "References")
  })
  observeEvent(input$prevsum,{
    updateTabsetPanel(session, "mainPanel", selected = "Summary")
  })
  observeEvent(input$nextint,{
    updateTabsetPanel(session, "mainPanel", selected = "Introduction")
  })
  observeEvent(input$prevint,{
    updateTabsetPanel(session, "mainPanel", selected = "Introduction")
  })
  observeEvent(input$nextdis,{
    updateTabsetPanel(session, "mainPanel", selected = "Discussion")
  })
  observeEvent(input$prevdis,{
    updateTabsetPanel(session, "mainPanel", selected = "Discussion")
  })
  observeEvent(input$nextmed,{
    updateTabsetPanel(session, "mainPanel", selected = "Medians")
  })
  observeEvent(input$prevmed,{
    updateTabsetPanel(session, "mainPanel", selected = "Medians")
  })
  observeEvent(input$nextres,{
    updateTabsetPanel(session, "mainPanel", selected = "Results")
  })
  observeEvent(input$prevres,{
    updateTabsetPanel(session, "mainPanel", selected = "Results")
  })
  observeEvent(input$nextcon,{
    updateTabsetPanel(session, "mainPanel", selected = "Conclusions")
  })
  observeEvent(input$prevcon,{
    updateTabsetPanel(session, "mainPanel", selected = "Conclusions")
  })
  observeEvent(input$nextref,{
    updateTabsetPanel(session, "mainPanel", selected = "References")
  })


  #output$clicked_point <- renderText(input$selected_point)
  output$click_selected <- renderText(paste0(Dtype[input$selected_point], " ", treeNames[input$selected_point]))
  output$click_selected2 <- renderText(paste0(Dtype[input$selected_point], " ", treeNames[input$selected_point]))
  output$click_selected3 <- renderText(paste0(Dtype[input$selected_point], " ", treeNames[input$selected_point]))
  #output$pplot_dblclicked <- renderText(paste0("DoubleClick : ", input$pplot_dblclick))
  
#   observeEvent(input$dataSet, {
#     input$selected_point<-NULL
#   })
#   eventReactive(input$dataSet, {
#     input$selected_point<-NULL
#   })

#     observeEvent(input$pplot_dblclick, {
#       showModal(modalDialog(
#         renderPlot({plotCTree(cTreeListRoetz[[input$selected_point]])}), size="l"
#       )) })
      observeEvent(input$pplot_dblclick, {
      output$pplot_dblclicked <- 
        renderPlot({
          if (input$selected_point<=nTransPhylo){
            if (input$dataSet == "veg"){
              plotCTree(cTreeList[[input$selected_point]])}
            else{
              plotCTree(cTreeListRoetz[[input$selected_point]])}
          }
          if ((input$selected_point>nTransPhylo) & (input$selected_point<=nTransPhylo+nMcmc)){
            if (input$dataSet == "veg"){
              state <- myMCMCstate}
            else{
              state <- roetzMCMCstate}
            plot(state, plot.which = "sample", samplenr=(input$selected_point-nTransPhylo))}
          if (input$selected_point>(nTransPhylo+nMcmc) & (input$selected_point<=nTransPhylo+nMcmc+nScotti)){
            if (input$dataSet == "veg"){
              plot(scPTreeList[[(input$selected_point)-nTransPhylo-nMcmc]])}
            else{
              plot(scPTreeListRoetz[[(input$selected_point)-nTransPhylo-nMcmc]])}
          }
          if (input$selected_point>(nTransPhylo+nMcmc+nScotti)){
            if (input$dataSet == "veg"){
              plot(bPTreeList[[(input$selected_point)-nTransPhylo-nMcmc-nScotti]])}
            else{
              plot(bPTreeListRoetz[[(input$selected_point)-nTransPhylo-nMcmc-nScotti]])}
          }
        })
      updateTabsetPanel(session, "mainPanel", selected = "Phylogeny")
      })  

    observeEvent(input$tplot_dblclick, {
      output$tplot_dblclicked <- 
        renderPlot({
        if (input$dataSet == "veg") {
          plotTTree(cTTreeList[[input$selected_point]],w.shape=1.3, w.scale=1/0.3)}
      else {
        plotTTree(cTTreeListRoetz[[input$selected_point]],w.shape=1.3, w.scale=1/0.3)}  
        })
      updateTabsetPanel(session, "mainPanel", selected = "Transmission")
         
    })  

    observeEvent(input$itplot_dblclick, {
      output$itplot_dblclicked <- 
        renderPlot({
        if (input$dataSet == "veg") {
          plotTTree(cTTreeList[[input$selected_point]],w.shape=1.3, w.scale=1/0.3)}
      else {
        plotTTree(cTTreeListRoetz[[input$selected_point]],w.shape=1.3, w.scale=1/0.3)} 
        })
      updateTabsetPanel(session, "mainPanel", selected = "Trans Network")
         
    })  
  
  output$firstPlotTitle <- renderText({
    if(length(input$selected_point)>0){
      if (input$selected_point<=nTransPhylo){
        paste0("Phylogenetic tree with transmission events and coloured hosts")}
      if ((input$selected_point>nTransPhylo) & (input$selected_point<=nTransPhylo+nMcmc)){
        paste0("Phylogenetic tree with transmission events and coloured hosts")}
      if (input$selected_point>(nTransPhylo+nMcmc) & (input$selected_point<=nTransPhylo+nMcmc+nScotti)){
        paste0("Phylogenetic tree")
      }
      if (input$selected_point>(nTransPhylo+nMcmc+nScotti)){
        paste0("Phylogenetic tree")
      }
    }
  })

  output$chosenCTree <- renderPlot({
    if(length(input$selected_point)>0){
    if (input$selected_point<=nTransPhylo){
      if (input$dataSet == "veg"){
        plotCTree(cTreeList[[input$selected_point]])}
      else{
        plotCTree(cTreeListRoetz[[input$selected_point]])}
      }
    if ((input$selected_point>nTransPhylo) & (input$selected_point<=nTransPhylo+nMcmc)){
      if (input$dataSet == "veg"){
        state <- myMCMCstate}   # MAYBE THESE SHOULD BE FROM FILE
      else{
        state <- roetzMCMCstate}
      plot(state, plot.which = "sample", samplenr=(input$selected_point-nTransPhylo))}
    if (input$selected_point>(nTransPhylo+nMcmc) & (input$selected_point<=nTransPhylo+nMcmc+nScotti)){
      if (input$dataSet == "veg"){
        plot(scPTreeList[[(input$selected_point)-nTransPhylo-nMcmc]])}
      else{
        plot(scPTreeListRoetz[[(input$selected_point)-nTransPhylo-nMcmc]])}
    }
    if (input$selected_point>(nTransPhylo+nMcmc+nScotti)){
      if (input$dataSet == "veg"){
        plot(bPTreeList[[(input$selected_point)-nTransPhylo-nMcmc-nScotti]])}
      else{
        plot(bPTreeListRoetz[[(input$selected_point)-nTransPhylo-nMcmc-nScotti]])}
    }
    }
  })

  output$hairy <- renderPlot({
    plot(butterfly,ylab='Posterior probability',xlab='MCMC iterations',type='l',main='Figure 1. TransPhylo transmission tree posterior probabilites for Norwegian TB cluster.')
  })

  output$chosenTTree <- renderPlot({
    if(length(input$selected_point)>0){
    if (input$dataSet == "veg") {
      plotTTree(cTTreeList[[input$selected_point]],w.shape=1.3, w.scale=1/0.3)}
    else {
      plotTTree(cTTreeListRoetz[[input$selected_point]],w.shape=1.3, w.scale=1/0.3)}
    }
  })
  
  output$chosenNetworkPlot <- renderVisNetwork({
    if(length(input$selected_point)>0){
    if (input$dataSet == "veg"){
      networkTPlot(cTTreeList[[input$selected_point]])}
    else {
      networkTPlot(cTTreeListRoetz[[input$selected_point]])}
    }
  })

  output$MedianTree1 <- renderPlot({
    plotTTree(MedianTree1,w.shape=1.3, w.scale=1/0.3)
  })
  output$MedianTree2 <- renderPlot({
    plotTTree(MedianTree2,w.shape=1.3, w.scale=1/0.3)
  })
  output$MedianTree3 <- renderPlot({
    plotTTree(MedianTree3,w.shape=1.3, w.scale=1/0.3)
  })
  output$MedianTree4 <- renderPlot({
    plotTTree(MedianTree4,w.shape=1.3, w.scale=1/0.3)
  })

  output$MedianPlot1 <- renderVisNetwork({
    networkTPlot(MedianTree1)
  })
  output$MedianPlot2 <- renderVisNetwork({
    networkTPlot(MedianTree2)
  })
  output$MedianPlot3 <- renderVisNetwork({
    networkTPlot(MedianTree3)
  })
  output$MedianPlot4 <- renderVisNetwork({
    networkTPlot(MedianTree4)
  })

  url4 <- a("View Article.", href="https://academic.oup.com/mbe/article-lookup/doi/10.1093/molbev/msw275")
    output$ref4 <- renderUI({
    tagList("[4] Didelot X, Fraser C, Gardy J, Colijn C (2017) Genomic Infectious Disease Epidemiology in Partially Sampled and Ongoing Outbreaks. Mol Biol Evol. 34 (4): 997-1007.\n", url4)
  })
  url9 <- a("View Article.", href="https://academic.oup.com/mbe/article-lookup/doi/10.1093/molbev/msu121")
    output$ref9 <- renderUI({
    tagList("[9] Didelot X, Gardy J, Colijn C (2014) Bayesian Inference of Infectious Disease Transmission from Whole-Genome Sequence Data. Mol Biol Evol. 31 (7):1869–79.", url9)
  })
  url2 <- a("View Article.", href="http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005495")
    output$ref2 <- renderUI({
    tagList("[2] Klinkenberg D, Backer JA, Didelot X, Colijn C, Wallinga J (2017) Simultaneous inference of phylogenetic and transmission trees in infectious disease outbreaks. PLoS Comput Biol. 13(5).", url2)
  })
  url6 <- a("View Article.", href="https://arxiv.org/abs/1609.09051")
    output$ref6 <- renderUI({
    tagList("[6] Kendall M, Ayabina D, Colijn C (2016) Estimating transmission from genetic and epidemiological data: a metric to compare transmission trees. arXiv:1609.09051v1.", url6)
  })
  url3 <- a("View Article.", href="http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005130")
  output$ref3 <- renderUI({
    tagList("[3] De Maio N, Wu C-H, Wilson DJ (2016) SCOTTI: Efficient Reconstruction of Transmission within Outbreaks with the Structured Coalescent. PLoS Comput Biol. 12(9).", url3)
  })
  url1 <- a("View Article.", href="http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004613")
  output$ref1 <- renderUI({
    tagList("[1] Hall M, Woolhouse M, Rambaut A (2015) Epidemic Reconstruction in a Phylogenetics Framework: Transmission Trees as Partitions of the Node Set. PLoS Comput Biol. 11(12).", url1)
  })
  url8 <- a("View Article.", href="https://CRAN.R-project.org/package=shiny")
  output$ref8 <- renderUI({
    tagList("[8] Chang W, Cheng J, Allaire JJ, Xie Y, McPherson J (2017) Shiny: Web Application Framework for R. R package version 1.0.3.", url8)
  })
  url5 <- a("View Article.", href="http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1001387")
  output$ref5 <- renderUI({
    tagList("[5] Roetzer A, Diel R, Kohl TA, Rückert C, Nübel U, Blom J, Wirth T, Jaenicke S, Schuback S, Rüsch-Gerdes S, Supply P, Kalinowski J, Niemann S (2013) Whole Genome Sequencing versus Traditional Genotyping for Investigation of a Mycobacterium tuberculosis Outbreak: A Longitudinal Molecular Epidemiological Study.", url5)
  })
  url7 <- a("View Article.", href="https://link.springer.com/chapter/10.1007%2F978-3-540-33037-0_14")
  output$ref7 <- renderUI({
    tagList("[7] Cox TF, Cox MAA (2000) Multidimensional Scaling. CRC Press.")#, url7)
  })
  


}