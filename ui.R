#
# Define the user interface of the app
#
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

fluidPage(
  # allows HTML element to span a few columns
  # (used for checkboxes)
  tags$head(
    tags$style(HTML("
    .multicol {
    -webkit-column-count: 4; /* Chrome, Safari, Opera */
    -moz-column-count: 4;    /* Firefox */
    column-count: 4;         /* generic */
    }"))
  ),

  # content
  fluidRow(column(12,

    # selector for classification dates
    # = applies to all elements below
    # = allows to monitor progress through time
    fluidRow(column(10, offset=1,
      sliderInput("dates", label="Date range of classifications to consider", min=dates[1], max=dates[2], value=dates, width="100%")
    )),


    tabsetPanel(

      # inspect classifications through time
      tabPanel("Timeline",
        checkboxInput("log", label="Log-transform nb of classifications", value=FALSE),
        plotOutput("p_time"),
        verbatimTextOutput("t_stats"),
        div()
      ),

      # inspect number of classifs per user
      tabPanel("Users",
        sliderInput("minn", label="Display only users with more than this number of classifications", min=100, max=10000, step=100, value=1000),
        plotOutput("p_who"),
        div()
      ),

      # inspect effort in terms of number of users per frame
      tabPanel("Effort",
        p("This is a bit long to compute. Please wait."),
        plotOutput("p_effort"),
        div()
      ),

      # inspect classification results, in their real life setting
      tabPanel("Results",
        fluidRow(
          column(2, checkboxGroupInput("transects", "Plot transects", choices=transects, selected=transects)),
          column(8, class="multicol", checkboxGroupInput("species", "Plot species", choices=species, selected="copepod")),
          column(2, selectInput("geom", "Plot as", choices=c("Points", "Density", "Both"), selected="Both"))
        ),
        fluidRow(column(12, plotOutput("p_transects", height="auto"))),
        fluidRow(
          column(6, h4("Stats per transect"), tableOutput("t_transects_t")),
          column(6, h4("Stats per species"), tableOutput("t_transects_s"))
        )
      )
    )
  ))
)
