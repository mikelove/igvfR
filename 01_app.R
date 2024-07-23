library(shiny)
library(tibble)
library(dplyr)
library(DT)

# page layout: left-side input, right-side output
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
      tableOutput("regulatory_regions")
    )
  )
)

server <- function(input, output) {
  makeTable <- eventReactive(input$go, {
    
    # lookup ENSG
    gene_id <- sym2gene |> 
      dplyr::filter(symbol == input$query_symbol) |>
      dplyr::pull(gene_id)
    
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
    
    # e2g table
    tab <- as_tibble(do.call(rbind, lapply(hits, as.data.frame)))

    # return with new column names and sorted by score
    tab |>
      dplyr::select(from=X_from, to=X_to, 
                    score=score.long, 
                    context=biological_context) |>
      dplyr::mutate(
        from = sub("regulatory_regions/", "", from),
        to = sub("genes/", "", to)
      ) |>
      dplyr::arrange(desc(score))
    
  })
  output$regulatory_regions <- renderTable({
    tab <- makeTable()
    tab
  })
}

shinyApp(ui = ui, server = server)
