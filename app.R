## app.R

version <- "0.7"
update  <- "Jan 13th, 2026"
changes <- paste(
  "Removed BSgenome dependency (Ensembl sequence/region).", "\t", 
  "Forced protospacer mode requires SNP within editing window; enzyme assigned by A/C in window.","\t",
  "Progress bar + no-guides modal.","\t",
  "Improved reactivity: auto-run (debounced) for single rsID, and reliable single-click execution."
)

library(shiny)
library(DT)
library(dplyr)
library(purrr)
library(shinyjs)

source("global.R")

ui <- fluidPage(
  shinyjs::useShinyjs(),
  pageWithSidebar(
    headerPanel(paste0("snpGuide v", version)),

    sidebarPanel(
      p(""),
      em(paste0("Updated ", update, ".")),
      p(""),
      em(paste0("Latest Change: ", changes)),
      p(""),
      em("Please email yuriy.baglaenko@cchmc.org with any questions, comments, or concerns."),
      p(""),
      tags$a(href = "https://www.baglaenkolab.com", em("Check us out at baglaenkolab.com")),
      p(""),

      textInput("pam", "PAM sequence", value = "NGN"),
      p(""),

      sliderInput(
        inputId = "editing_window",
        label = "Choose base editing window (protospacer positions)",
        min = 1, max = 20,
        value = c(4, 8),
        step = 1
      ),
      p(""),

      textInput("rsID", "Input rsID to target: example rs61839660", value = ""),
      p(""),

      fileInput("upload", "Upload batch rsIDs as .txt file", accept = ".txt"),
      downloadButton("download_example", "Download example input"),
      p(""),
      tags$hr(),

      actionButton("run_normal", "Design guides (requires editable transition)"),
      p(""),
      actionButton("run_forced", "Design protospacers (SNP in window; enzyme by A/C in window)"),
      p(""),

      checkboxGroupInput(
        "force_bases",
        "Forced mode: bases to consider editable in window",
        choices = c("A", "C", "G", "T"),
        selected = c("A", "C")
      ),
      checkboxInput("both_strands", "Forced mode: search both strands", value = TRUE),
      checkboxInput(
        "keep_only_editable",
        "Forced mode: keep only protospacers with editable bases in window",
        value = TRUE
      ),

      p(""),
      checkboxInput(
        "auto_run",
        "Auto-run for single rsID after typing stops (recommended)",
        value = TRUE
      )
    ),

    mainPanel(
      h2("Predicted guides"),
      br(),
      DTOutput("df"),
      br(),
      downloadButton("download_guide", "Download all guides")
    )
  )
)

server <- function(input, output, session) {

  # -------------------------
  # Modal must be defined BEFORE it can be called reliably
  # -------------------------
  show_no_guides_modal <- function() {
    showModal(modalDialog(
      title = "No guides found",
      HTML(
        "<b>No valid protospacers were designed.</b><br><br>",
        "Possible reasons:<br>",
        "<ul>",
        "<li>No PAM sites in this region</li>",
        "<li>SNP not positioned within the editing window</li>",
        "<li>No editable bases (e.g., A/C) in the window given your settings</li>",
        "<li>Variant type incompatible with selected mode (normal mode only)</li>",
        "</ul>",
        "Try adjusting the editing window, PAM, or using forced mode."
      ),
      easyClose = TRUE,
      footer = modalButton("OK")
    ))
  }

  # -------------------------
  # State
  # -------------------------
  results <- reactiveVal(NULL)
  running <- reactiveVal(FALSE)
  mode <- reactiveVal("normal")  # "normal" or "forced"

  # Debounced snapshot of inputs (prevents spamming while typing)
  inputs_debounced <- reactive({
    list(
      rsID = trimws(input$rsID),
      pam = input$pam,
      editing_window = input$editing_window,
      upload = input$upload,
      force_bases = input$force_bases,
      both_strands = input$both_strands,
      keep_only_editable = input$keep_only_editable
    )
  }) %>% debounce(700)

  # -------------------------
  # Compute helper with progress + button disable
  # -------------------------
  run_with_progress <- function(rsids, force_mode) {
    rsids <- unique(trimws(rsids))
    rsids <- rsids[nzchar(rsids)]
    n <- length(rsids)
    if (n == 0) return(NULL)

    withProgress(message = "Designing guides...", value = 0, {
      incProgress(0.05, detail = "Initializing")

      out_list <- vector("list", n)
      for (i in seq_along(rsids)) {
        incProgress(0.9 / n, detail = paste0("Processing ", rsids[i], " (", i, "/", n, ")"))

        out_list[[i]] <- purrr::possibly(snpGuide, otherwise = NULL)(
          rsid = rsids[i],
          pam = input$pam,
          guide_length = 20,
          min_editing_window = input$editing_window[1],
          max_editing_window = input$editing_window[2],
          flank = 20,
          force_protospacers = force_mode,
          force_edit_bases = input$force_bases,
          search_both_strands = input$both_strands,
          keep_only_window_has_editable = input$keep_only_editable
        )
      }

      incProgress(0.05, detail = "Finalizing results")
      bind_rows(out_list)
    })
  }

  run_analysis <- function(force_mode) {
    if (isTRUE(running())) return(invisible(NULL))

    # Validate input early (no silent req() failures)
    if (!nzchar(trimws(input$rsID)) && is.null(input$upload)) {
      showNotification("Enter an rsID or upload a batch .txt file.", type = "warning", duration = 5)
      return(invisible(NULL))
    }

    running(TRUE)
    shinyjs::disable("run_normal")
    shinyjs::disable("run_forced")

    on.exit({
      running(FALSE)
      shinyjs::enable("run_normal")
      shinyjs::enable("run_forced")
    }, add = TRUE)

    out <- if (is.null(input$upload)) {
      run_with_progress(rsids = trimws(input$rsID), force_mode = force_mode)
    } else {
      rsids <- readLines(input$upload$datapath)
      run_with_progress(rsids = rsids, force_mode = force_mode)
    }

    results(out)

    if (is.null(out) || nrow(out) == 0) show_no_guides_modal()

    invisible(NULL)
  }

  # -------------------------
  # Buttons: always run on first click
  # -------------------------
  observeEvent(input$run_normal, {
    mode("normal")
    run_analysis(force_mode = FALSE)
  }, ignoreInit = TRUE)

  observeEvent(input$run_forced, {
    mode("forced")
    run_analysis(force_mode = TRUE)
  }, ignoreInit = TRUE)

# Auto-run (single rsID only) after typing stops
observeEvent(inputs_debounced(), {
  if (!isTRUE(input$auto_run)) return()

  x <- inputs_debounced()

  # Don't auto-run in batch mode
  if (!is.null(x$upload)) return()

  # Must have rsID
  if (!nzchar(x$rsID)) return()

  # Avoid re-entrancy
  if (isTRUE(running())) return()

  # Run immediately inside reactive context
  run_analysis(force_mode = (mode() == "forced"))
}, ignoreInit = TRUE)

  # -------------------------
  # Outputs
  # -------------------------
  output$df <- renderDT({
    req(results())
    DT::datatable(
      results(),
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$download_guide <- downloadHandler(
    filename = function() paste(Sys.Date(), "snpGuides", "csv", sep = "_"),
    content = function(file) {
      g <- results()
      if (is.null(g) || nrow(g) == 0) {
        write.csv(data.frame(), file, row.names = FALSE)
      } else {
        write.csv(g, file, row.names = FALSE)
      }
    }
  )

  output$download_example <- downloadHandler(
    filename = function() "example_input.txt",
    content = function(file) {
      writeLines(
        c("rs12722517", "rs7920946", "rs791593", "rs11597237", "rs12722497"),
        file
      )
    }
  )
}

shinyApp(ui, server)
