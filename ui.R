# Author: Matt Watts
# Date: 10 Dec 2014
# Purpose: RunMarxanTas web app ui.R

require(shiny)

shinyUI(pageWithSidebar(

    headerPanel("Tas Activity: Run Marxan"),

    sidebarPanel(
        actionButton("mrun","Run"), 
        br(),
        br(),
        textOutput("textfeedback"),
        br(),
        numericInput("blm", "Boundary Length Modifier:",1,min=0),
        br(),
        br(),
        selectInput("feature", "Choose a species to edit:",
                    choices = c("All features","Acaciaforests","Acaciashrublands","Eucshrublands","Eucalyptuswoodlands","Freshwater","Grasslands",
                                "Lowforest","MalleeShrublands","Mulgawoodlands","Openshrublands","Saltbush","Saltlagoon","Sedgelands",
                                "TropicalRainforest","Tussockgrasslands","Azurekingfsh","Blindvelvetworm","Bartailedgodwit",
                                "Claspleathheath","Fairywren","Hoarysunray","MaskedOwl","OrangeParrot","SeaEagle","SwiftParrot")),
        numericInput("prop", "Target:",0.1,min=0,max=1,step=0.1),
        numericInput("spf", "SPF:",10,min=0),
        br(),
        actionButton("savetargetspf","Save Target and SPF"),
        conditionalPanel(condition = "input.tabs == 'Map'",
            br(),
            br(),
            radioButtons("map", "Map to display:",
                         list("Best solution" = "bestmap",
                              "Solution M" = "runMmap",
                              "Selection frequency" = "ssolnNmap"))
                    ),
        conditionalPanel(condition = "input.tabs == 'Map' & (input.map == 'bestmap' | input.map == 'runMmap' | input.map == 'ilpmap')",
                         HTML("<img src='http://marxan.net/images/green.png' /></a>"),
                         HTML("Selected"),
                         br(),
                         HTML("<img src='http://marxan.net/images/turquoise.png' /></a>"),
                         HTML("Existing Reserves"),
                         br(),
                         HTML("<img src='http://marxan.net/images/grey.png' /></a>"),
                         HTML("Excluded")
                        ),
        conditionalPanel(condition = "input.tabs == 'Map' & input.map == 'ssolnNmap'",
                         HTML("Selection frequency"),
                         br(),
                         HTML("<img src='http://marxan.net/images/blue10.png' /></a>"),
                         HTML("91-100"),
                         br(),
                         HTML("<img src='http://marxan.net/images/blue9.png' /></a>"),
                         HTML("81-90"),
                         br(),
                         HTML("<img src='http://marxan.net/images/blue8.png' /></a>"),
                         HTML("71-80"),
                         br(),
                         HTML("<img src='http://marxan.net/images/blue7.png' /></a>"),
                         HTML("61-70"),
                         br(),
                         HTML("<img src='http://marxan.net/images/blue6.png' /></a>"),
                         HTML("51-60"),
                         br(),
                         HTML("<img src='http://marxan.net/images/blue5.png' /></a>"),
                         HTML("41-50"),
                         br(),
                         HTML("<img src='http://marxan.net/images/blue4.png' /></a>"),
                         HTML("31-40"),
                         br(),
                         HTML("<img src='http://marxan.net/images/blue3.png' /></a>"),
                         HTML("21-30"),
                         br(),
                         HTML("<img src='http://marxan.net/images/blue2.png' /></a>"),
                         HTML("11-20"),
                         br(),
                         HTML("<img src='http://marxan.net/images/blue1.png' /></a>"),
                         HTML("1-10"),
                         br(),
                         HTML("<img src='http://marxan.net/images/turquoise.png' /></a>"),
                         HTML("Existing Reserves"),
                         br(),
                         HTML("<img src='http://marxan.net/images/grey.png' /></a>"),
                         HTML("Excluded")
                        ),
        conditionalPanel(condition = "input.tabs == 'Table'",
            br(),
            br(),
            radioButtons("table", "Table to display:",
                         list("Summary" = "sumtable",
                              "Conservation Features" = "spec",
                              "Best solution Missing values" = "mvbesttable",
                              "Solution M Missing values" = "mvNtable"))
        ),
        conditionalPanel(condition = "(input.tabs == 'Map' & input.map == 'runMmap') | (input.tabs == 'Table' & input.table == 'mvNtable')",
            br(),
            br(),
            sliderInput("m", "Solution M:",
                        value = 1,
                        min = 1,
                        max = 100, step = 1)

        ),
        conditionalPanel(condition = "input.prop == -1",
                         numericInput("refreshinput", "Refresh Input", 0))
    ),

    mainPanel(
        tabsetPanel(id="tabs",
            tabPanel("Map", plotOutput('marxanmap')),
            tabPanel("Table", tableOutput('marxantable')),
            tabPanel("Cluster", plotOutput("plot2ds"),
                                plotOutput("plotdendogram"))
                  )
             )
))
