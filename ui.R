library(shiny)
library(shinyBS)
library(shinythemes)
library(scatterD3)
library(markdown) # for the navbar code
library(TransPhylo)
library(visNetwork)

navbarPage("",
           theme = shinytheme("cerulean"),
           id = "mainPanel",
           
  tabPanel("Summary",
           HTML("<a name='sumtop'></a>"),      
    titlePanel("How much can we learn about disease transmission from genetic data?"),
    
    hr(),
    fluidRow(
      column(1
      ),
      column(10,
             titlePanel("Lay Summary"),
             hr(),
             includeHTML("LaySummary.html")
      ),
      column(1
      )
    ),
    fluidRow(
      column(1
      ),
      column(10,
             hr(),
             titlePanel("Abstract"),
             hr(),
             includeHTML("abstract.html")#,
             #actionLink("reftop", "RefPage")
      ),
      column(1
      )
    ),
    fluidRow(
      column(10,
             hr(),
             HTML("<strong>James Stimson (June 2017)</strong>")
       ),
      column(2, hr(), HTML("<a href='#sumtop'>Top</a>"),actionLink("nextint", "Next>"))
      
     
    )
    ),
  tabPanel("Introduction",HTML("<a name='inttop'></a>"),
           titlePanel("Approach to the problem, and some details of methods"),
           
           hr(),
           fluidRow(
             column(1),
             column(10,
                    HTML("<h4>Background and Motivation for Approach</h4>"),
                    HTML("<p>The problem of discovering exactly how infectious diseases spread has in general been difficult, time-consuming and expensive to solve using traditional epidemiological techniques. 
                         </p>"),
                    HTML("<p>In recent years the use of genetic data from disease causing pathogens has given rise to the hope that much more rapid and accurate analysis of disease transmission mechanisms will be made possible.
                         There are several stochastic processes (see <a href=\"http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005130\">[3]</a>) which need to be modelled to capture 
                         the underlying dynamics fully, including <ul>
                        <li>Genetic mutation rate of the pathogen</li>
                        <li>Time from infection of host to transmission to another host</li>
                        <li>Time from infection to sampling</li>
                        <li>Time from coalescence in the pathogen's phylogeny to sampling or transmission time</li></ul>                          
                         To this end, Bayesian approaches have been widely employed to tackle the complexity of modelling required.</p>"),
                    HTML("<p>Identifying the pattern of disease transmission is also made more difficult where diseases have long latent periods in a host: 
                         in other words, where there is a long time between infection and presentation of the disease, as is the case with Tuberculosis (TB) which can lie dormant for many years.</p>"),
                    HTML("<p>Here, the <strong>objective</strong> is to show that the use of a particular metric <a href=\"https://arxiv.org/abs/1609.09051\">[6]</a> on the space of transmission trees 
                         provides a suitable way to compare the results of these approaches. </p>")
             ),
             column(1
             )
           ),
          fluidRow(
             column(1),
             column(10,
                    hr(),
                    HTML("<h4>Introduction to the Methods</h4>"),
                    HTML("<p>The four methods we are comparing - <strong>Beastlier</strong>, <strong>PhyBreak</strong>, <strong>SCOTTI</strong> and <strong>TransPhylo</strong> - have been developed with the goal of inferring transmission trees given the same essential input information, comprising genetic data of host pathogens together with the dates at which the samples were taken. </p>"),
                    HTML("<p>Unlike popular methods which preceded them, they model the population dynamics of the pathogen within hosts. 
                         This is important given the long latency periods of many diseases combined with rapid mutation rates. 
                          When transmissions between hosts are simulated, it assumed that this happens with only one host strain making the jump between hosts:
                         the transmission events form 'complete' bottlenecks.</p>"),
                    HTML("<p>However the model set up, underlying assumptions, and implementation details all differ.</p>")
             ),
             column(1
             )
           ),
           fluidRow(
             column(1
             ),
             column(10,
                    hr(),
                    HTML("<h4>Common Features</h4>"),
                    includeHTML("Introduction.html")
             ),
             column(1
             )
           ),
           fluidRow(
             column(1
             ),
             column(10,
                    hr(),
                    HTML("<h4>Differences</h4>"),
                    HTML("Here some of the main differences between the methods are highlighted, with particular reference to the inputs needed for the models. 
                         This is not a comprehensive exposition of all the differences in the models, for which the reader is referred to the original papers."),
                    actionLink("reftop", "[References 1-4]"),
                    includeHTML("Differences.html")
             ),
             column(1
             )
           ),
          fluidRow(
            column(10,
                   hr(),
                   HTML("<strong>James Stimson (June 2017)</strong>")
            ),
            column(2, hr(),actionLink("prevsum", "<Prev"), HTML("<a href='#inttop'>Top</a>"),actionLink("nextdis", "Next>"))
          )
  ),  
  tabPanel("Discussion",HTML("<a name='distop'></a>"),
           titlePanel("Methods and interpretation of results"),
           
           hr(),
           fluidRow(
             column(1),
             column(10,
           
                    HTML("<h4>Posterior Probabilities</h4>"),
                    p("When choosing the number of simulations to run for each algorithm and data set, the posterior probabilities of 
                         the MCMC are examined, in order to check that the process has converged to a stable solution. Plots of simulation 
                      against posterior probability can be drawn for this purpose. Figure 1, below, shows the plot for the simulations 
                      run for the Norwegian Cluster using TransPhylo. This was run for 100,000 simulations, with every tenth one of these use 
                      to draw the plot. The scale on the y axis is the natural logarithm of the transmission tree posterior probability.")              
             ),
             column(1)
           ),
           fluidRow(
             column(1),
             column(10,
                    plotOutput("hairy")
                    
             ),
             column(1)
           ),
           fluidRow(
             column(1),
             column(10,
                    p("Similar plots can be drawn for each of the other cases."),
                    HTML("SCOTTI was also run for 100,000 simulations, whereas Beastlier and PhyBreak were run with 20,000. Since the models are simlauting different processes, these numbers are not directly comparable."),
                    hr(),
                    includeHTML("random_notes.html"),
                    HTML("The metric used defines a Euclidean distance between each of the transmission trees. 
                         These distances are in a high dimensional space, so in order to visualise them we use multidimensional scaling (MDS)"),
                    actionLink("reftop2", "[7]"),
                    HTML("techniques to project them down to a two-dimensional representation: these are MDS Axis 1 and MDS Axis 2 in the main MDS plot on the 'Results' page.
                         The MDS algorithm works by using an eigenvalue decomposition on the multidimensional matrix of Euclidean distance between points. The two principal eigenvalues 
                         are then used to produce a two-dimensional representation for the scatter plot in which distances are preserved as well as possible, though not perfectly:
                         it is fundamentally impossible to preserve all aspects of high-dimensional structure when projecting down to two dimensions."),
                    hr(),
                    includeHTML("results_preamble.html")
                    
             ),
             column(1)
           ),
           fluidRow(
             column(10,
                    hr(),
                    HTML("<strong>James Stimson (June 2017)</strong>")
             ),
             column(2, hr(),actionLink("prevint", "<Prev"),HTML("<a href='#distop'>Top</a>"),actionLink("nextmed", "Next>"))
           )
  ),
  tabPanel("Medians",HTML("<a name='medtop'></a>"),
           titlePanel("A look at some representative data"),
           
           hr(),
           fluidRow(
             #column(1),
             column(12,
                 
                    HTML("<h4>Median Trees</h4>"),
                    HTML("<p>The transmission tree metric gives us a way of extracting representative trees from our result sets. 
                      Because the trees have well-defined distances between them, we can use an algorithm to select the tree which is most 'central' for any given set of trees.
                      These can then be used to compare visually different clusters against each other.</p>
                      <p>Below, we illustate this by showing side-by-side the median timed transmission trees and (further below) the transmission 
                      network for each the four methods for the <strong>Norwegian</strong> TB cluster.</p>"),
                    includeHTML("TransTrees.html")
             )#,
             #column(1)
           ),
           fluidRow(
             column(3,
                    hr(),
                    HTML("<h4>TransPhylo</h4>"),
                    HTML("<u>Median transmission tree against date</u>"),
                    plotOutput("MedianTree1")  
                    
             ),
             column(3,
                    hr(),
                    HTML("<h4>PhyBreak</h4>"),
                    HTML("<u>Median transmission tree against date</u>"),
                    plotOutput("MedianTree2")
             ),
             column(3,
                    hr(),
                    HTML("<h4>Scotti</h4>"),
                    HTML("<u>Median transmission tree against date</u>"),
                    plotOutput("MedianTree3")
             ),
             column(3,
                    hr(),
                    HTML("<h4>Beastlier</h4>"),
                    HTML("<u>Median transmission tree against date</u>"),
                    plotOutput("MedianTree4")          
             )
             
           ),
           
           fluidRow(
             column(3,
                    hr(),
                    HTML("<h4>TransPhylo</h4>"),
                    HTML("<u>Undated median transmission network</u>"),
                    visNetworkOutput("MedianPlot1")
             ),
             column(3,
                    hr(),
                    HTML("<h4>PhyBreak</h4>"),
                    HTML("<u>Undated median transmission network</u>"),
                    visNetworkOutput("MedianPlot2")             
             ),
             column(3,
                    hr(),
                    HTML("<h4>Scotti</h4>"),
                    HTML("<u>Undated median transmission network</u>"),
                    visNetworkOutput("MedianPlot3")
             ),
             column(3,
                    hr(),
                    HTML("<h4>Beastlier</h4>"),
                    HTML("<u>Undated median transmission network</u>"),
                    visNetworkOutput("MedianPlot4")
             )
             
           ),
           fluidRow(
             column(10,
                    hr(),
                    HTML("<strong>James Stimson (June 2017)</strong>")
             ),
             column(2, hr(),actionLink("prevdis", "<Prev"),HTML("<a href='#medtop'>Top</a>"),actionLink("nextres", "Next>"))
           )
  ),
  tabPanel("Results",HTML("<a name='restop'></a>"),
           tags$head(tags$style(HTML('

                        .modal-lg {
                       width: 100%;
                        height: 1000pt;
                        }
                      '))),
    fluidPage(
      
       
      titlePanel("Complete result set for transmission tree inference"),
      
      hr(),
      div(class="row",
          div(class="col-md-12",
              div(class="alert alert-warning alert-dismissible",
                  HTML('<button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>'),
                  HTML("<ul><li>Use the control panel on the left to customise the main display, which is <strong>zoomable</strong>.</li> 
                      <li>Hover over a dot to display tooltips. While the tooltip is visible, click to bring up the <strong>three</strong> individual plots below.</li>
                       <li>If in any doubt about current selection, check the 'Selection' displayed on the left hand side below the control panel.</li>
                       <li>Use the 'Data Set' toggle to switch between the two available data sets. Try points from each of the four methods: they have differing displays.</li>
                       <li>The three individual plots are discussed in the 'Individual Plots' section on the 'Discussion' page and on the 'Medians' page.</li></ul>")))),
      HTML("<h4><u>Scatter Plot: MDS projection of distances between transmission trees sampled from the posterior</u></h4>"),
      sidebarLayout(
        sidebarPanel(
          
        radioButtons("dataSet", "Data Set",
                       c("Norwegian TB Cluster (quick)"="veg", "Hamburg TB Outbreak (slight delay)"="roe")),
        checkboxInput("scatterD3_ellipses", "Confidence ellipses", value = TRUE),
        sliderInput("scatterD3_nb", "Confidence level (%) :",
                     min = 50, max = 100, step = 1, value = 95),
        checkboxInput("scatterD3_x_log", "Show point labels", value = FALSE),
        sliderInput("scatterD3_labsize", "Labels size :",
                    min = 6, max = 26, value = 10),
    
        sliderInput("scatterD3_point_size", "Points size :", min = 16, max = 96, value = 32, step = 1),
        sliderInput("scatterD3_opacity", "Points opacity :", min = 0, max = 1, value = 1, step = 0.05),
    #    checkboxInput("scatterD3_transitions", "Use transitions", value = TRUE),
        tags$p(actionButton("scatterD3-reset-zoom", HTML("<span class='glyphicon glyphicon-search' aria-hidden='true'></span> Reset Zoom")))
        
      ),
      mainPanel(scatterD3Output("scatterPlot"))    
   
    ),
    bsModal("modal_detail", "Your plot", "go", size = "large",plotOutput("modal_detail")),
    fluidRow(
      column(width = 1,
             hr(),
             HTML("<h4>Selection</h4>"),
             textOutput("click_selected3")          
      ),
      column(width = 10,
             #titlePanel("Phylogenetic Tree (double click to enlarge)"),
             hr(),
             HTML("<h4><u>Transmission Network (Movable and Zoomable)</u></h4>"),
             visNetworkOutput("chosenNetworkPlot", height="800px")     ## NEED TO ADD DOUBLE CLICK SOMEHOW (OR CLICK) OR ADD A BUTTON
      ),
      column(width = 1)
    ),   
    fluidRow(
      column(width = 1,
             hr(),
             HTML("<h4>Selection</h4>"),
             textOutput("click_selected2")          
      ),
      column(width = 10,
             #titlePanel("Transmission Tree (double click to enlarge)"),
             hr(),
             HTML("<h4><u>Timed Transmission Tree (Static)</u></h4>"),
             plotOutput("chosenTTree", dblclick = "tplot_dblclick", height="800px") #, height="500px", width="400px")
      ),
      column(width = 1)
    ),

    fluidRow(
      column(width = 1,
             hr(),
             HTML("<h4>Selection</h4>"),
             textOutput("click_selected")          
      ),
      column(width = 10,
            
             hr(),
             HTML("<h4><u>Timed Phylogenetic Tree (Static)</u></h4>"),
             #textOutput("firstPlotTitle"),
             plotOutput("chosenCTree", dblclick = "pplot_dblclick", height="800px") 
      ),
      column(width = 1)
    ),
    fluidRow(
      column(10,
             hr(),
             HTML("<strong>James Stimson (June 2017)</strong>")
      ),
      column(2, hr(),actionLink("prevmed", "<Prev"),HTML("<a href='#restop'>Top</a>"),actionLink("nextcon", "Next>"))
    ))
    
  ),
#   tabPanel("Phylogeny",
#      fluidPage(
#        
#        fluidRow(
#          #column(width=12, titlePanel(paste0("Phylogenetic Tree for ","click_selected")),
#          column(width=12, titlePanel(textOutput("click_selected2")),
#                 plotOutput("pplot_dblclicked", height="1000px", width="1000px")  )
#          ),
#        fluidRow(
#          column(12,
#                 hr(),
#                 HTML("<strong>James Stimson (June 2017)</strong>")
#          )
#        )
#      )
#   ),
#   tabPanel("Transmission",
#            fluidPage(
#              
#              fluidRow(
#                column(width=12, titlePanel("Transmission Tree"),
#                       plotOutput("tplot_dblclicked", height="1000px", width="1000px")  )
#              ),
#              fluidRow(
#                column(12,
#                       hr(),
#                       HTML("<strong>James Stimson (June 2017)</strong>")
#                )
#              )   
#            )
#   ),
#   tabPanel("Trans Network",
#            fluidPage(
#              
#              fluidRow(
#                column(width=12, titlePanel("Interactive Transmission Network"),
#                       plotOutput("itplot_dblclicked", height="1000px", width="1000px")  )
#              ),
#              fluidRow(
#                column(12,
#                       hr(),
#                       HTML("<strong>James Stimson (June 2017)</strong>")
#                )
#              )
#            )
#   ),
  

  tabPanel("Conclusions",HTML("<a name='contop'></a>"),
           titlePanel("Concluding remarks"),
           
           hr(),
         fluidRow(
           column(1),
           column(10,
                  includeHTML("conclusions.html"),
                  hr(),
                  HTML("<i>This website was built using <strong>shiny</strong>.</i><a href=\"https://CRAN.R-project.org/package=shiny\">[8]</a>")
                  
           ),
           column(1)
         ),
         fluidRow(
           column(10,
                  hr(),
                  HTML("<strong>James Stimson (June 2017)</strong>")
           ),
           column(2, hr(),actionLink("prevres", "<Prev"),HTML("<a href='#contop'>Top</a>"),actionLink("nextref", "Next>"))
         )
  ),
  tabPanel("References",
           #id = 'reftop',
    
    
    
    fluidRow(
      column(1),
      column(10,
             hr(),
             HTML("<a name='reftop'></a>"),
             HTML("<a name='reftop2'></a>"),
             HTML("<a name='reftop3'></a>"),
             HTML("<i>Thanks to Caroline Colijn for help and advice during the project.</i>"),
             hr(),
             HTML("<h4>References</h4>"),       
      
      #a(strong("[1] ")),
      br(),
      uiOutput("ref1"),
      br(),
      uiOutput("ref2"),
      br(),
      uiOutput("ref3"),
      br(),
      uiOutput("ref4"),
      br(),
      uiOutput("ref5"),
      br(),
      uiOutput("ref6"),
      br(),
      uiOutput("ref7"),
      br(),
      uiOutput("ref8"),
      #br(),
      #uiOutput("ref9"),
      br()
      ),
    column(1)
    ),
    
    fluidRow(
      column(10,
             hr(),
             HTML("<strong>James Stimson (June 2017)</strong>")
      )
      ,
      column(2, hr(),actionLink("prevcon", "<Prev"),HTML("<a href='#reftop3'>Top</a>"))
    )

  )
)



