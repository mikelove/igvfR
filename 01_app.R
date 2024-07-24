library(shiny)
library(tibble)
library(dplyr)
library(tidyr)

# UI page layout: left-side input, right-side output
ui <- fluidPage(
  titlePanel("IGVF shiny demo"),
  sidebarLayout(
    sidebarPanel(
      textInput("query_symbol", "Gene symbol", value=""),
      numericInput("score_thres", "Minimum score", value=.85, 
                   min=0, max=1, step = .01),
      actionButton("go", "Go!")
    ),
    mainPanel(
      plotOutput("score_histogram"),
      tableOutput("regulatory_regions")
    )
  )
)

# server defines the logic of the page, builds the table, etc.
server <- function(input, output) {
  
  # this happens when you press 'Go'
  observeEvent(input$go, {
    
    # lookup ENSG
    gene_id <- sym2gene |> 
      dplyr::filter(symbol == input$query_symbol) |>
      dplyr::pull(gene_id)
    
    # below we use the Arango Query Language to query the db,
    # alternatively could use the HTTP API for simple queries
    
    # send query to db, retrieve list of results (regions)
    cursor <- db$aql$execute(
      "FOR l in regulatory_regions_genes \
       FILTER l._to == @geneid \
       FILTER l.`score:long` > @thres \
       return l",
      bind_vars=list(geneid=paste0("genes/",gene_id), 
                     thres=input$score_thres)
    )
    
    # make a list of hits until the cursor is empty
    hits <- list()
    while (!cursor$empty()) {
      hits <- c(hits, list(cursor$pop()))
    }
    
    # region to gene table
    tab0 <- as_tibble(
      do.call(
        rbind, lapply(hits, as.data.frame)
        )
      )

    # new column names and sorted by score
    tab <- tab0 |>
      dplyr::select(from=X_from,
                    score=score.long, 
                    source=source,
                    context=biological_context) |>
      dplyr::mutate(
        from = sub("regulatory_regions/", "", from),
        context = sub("ontology_terms/", "", context)) |>
      tidyr::separate(from, 
                      into=c("type","chrom","start","end","ref")) |>
      dplyr::arrange(desc(score))
  
    # histogram of scores
    output$score_histogram <- renderPlot({
      hist(tab$score, xlab="score", main="")
    })
    
    # render a table of the region data
    output$regulatory_regions <- renderTable({
      tab
    })
    
  })
      
}

# run the app
shinyApp(ui = ui, server = server)
