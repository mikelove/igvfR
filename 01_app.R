library(shiny)
library(tibble)
library(dplyr)
library(tidyr)
# library(plyranges) # just using as_granges()
library(plotgardener)

# UI page layout: left-side input, right-side output
ui <- fluidPage(
  titlePanel("IGVF shiny demo"),
  sidebarLayout(
    sidebarPanel(
      textInput("query_symbol", "Gene symbol", value="GCK"),
      numericInput("score_thres", "Minimum score", value=.85, 
                   min=0, max=1, step = .01),
      actionButton("go", "Go!")
    ),
    mainPanel(
      plotOutput("genomic_plot"),
      tableOutput("regulatory_regions_table")
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
      dplyr::select(-ref) |>
      dplyr::mutate_at(c("start","end"), as.integer) |>
      dplyr::arrange(desc(score))
  
    # make a GRanges object of regions
    regions_gr <- tab |>
      dplyr::rename(seqnames = chrom) |>
      plyranges::as_granges()
    
    mid <- round(mean(tab$start))
    halfwindow <- 1e5
    par <- pgParams(
      chrom = tab$chrom[1], 
      chromstart = mid - halfwindow,
      chromend = mid + halfwindow,
      assembly = "hg38", 
      just = c("left", "bottom")
    )
    
    pal <- colorRampPalette(c("blue3", "purple"))
    
    # draw the regulatory regions in a plotgardener plot
    output$genomic_plot <- renderPlot({
      pageCreate(width = 10, height = 4, 
                 showGuides = FALSE)
      plotRanges(
        regions_gr,
        fill = colorby("score", palette=pal),
        params = par, 
        x = .5, y = 2.5,
        width = 9, height = 2
      )
      plotGenes(
        params = par,
        geneHighlights = data.frame(
          "gene" = c(input$query_symbol),
          "color" = c("magenta")
        ),
        x = .5, y = 3.5,
        width = 9, height = 1
      )
      plotGenomeLabel(
        params = par,
        x = .5, y = 3.75,
        length = 9
      )
    })
    
    # render a table of the region data
    output$regulatory_regions_table <- renderTable({
      tab
    })
    
  })
      
}

# run the app
shinyApp(ui = ui, server = server)
