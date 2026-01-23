library(shiny)
library(leaflet)
library(sf)
library(dplyr)
library(RColorBrewer)
library(htmltools)
library(here)
library(bslib)
library(shinyBS)
library(stringr)
library(purrr)


# ----------------------------------
# 0) Set directory for ML logo
# ----------------------------------

register_assets <- function() {
  assets_dir <- here::here("www")
  if (dir.exists(assets_dir)) {
    shiny::addResourcePath("assets", assets_dir)
  } else {
    warning("Assets directory not found: ", assets_dir)
  }
}
register_assets()




# -------------------------------
# 1) Load & prepare spatial data
# -------------------------------
# LSOA boundaries
map_data <- st_read(here("Data","map_data.shp")) 

#check geometry exists - should return FALSE

any(is.na(sf::st_geometry(map_data)))

# check CRS

sf::st_crs(map_data)



# --------------------------
# 2) Load ICB patient data
# --------------------------

ICB_data <- read.csv(here("Outputs","ICB_data.csv"))


# --------------------------------------
# 3A) Distinct mapping of PCN to Sub-ICB
# --------------------------------------
pcn_map <- map_data |>
  st_drop_geometry() |>
  distinct(Sub_ICB, PCN_NAME)


# -------------------------------------------
# 3B) Load fixed PCN colour palette from CSV
# -------------------------------------------

`%||%` <- function(a,b) if (is.null(a)) b else a

# Validate hex
is_hex <- function(x) {
  grepl("^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", x %||% "", perl = TRUE)
}

# Build the distinct PCN list from the map 
pcn_map <- map_data |>
  sf::st_drop_geometry() |>
  dplyr::distinct(Sub_ICB, PCN_NAME)

# Load & normalise your colour file
pcn_colors_path <- here::here("Data", "pcn_colours.csv")
if (!file.exists(pcn_colors_path)) {
  stop(sprintf("Cannot find file: %s", pcn_colors_path))
}

pcn_colors_raw <- read.csv(pcn_colors_path, stringsAsFactors = FALSE, check.names = FALSE)

# Normalise column names to our canonical: Sub_ICB, PCN_NAME, color
# Accepts 'PCN_Name' or 'PCN_NAME', and 'Colour' (UK) or 'Color' (US)
pcn_colors <- pcn_colors_raw |>
  dplyr::rename(
    Sub_ICB  = dplyr::any_of(c("Sub_ICB")),
    PCN_NAME = dplyr::any_of(c("PCN_NAME","PCN_Name")),
    color    = dplyr::any_of(c("Color","Colour"))
  ) |>
  dplyr::mutate(
    Sub_ICB  = trimws(Sub_ICB),
    PCN_NAME = trimws(PCN_NAME),
    color    = toupper(trimws(color))
  )

# Check required columns are present after rename
missing_cols <- setdiff(c("Sub_ICB","PCN_NAME","color"), names(pcn_colors))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "pcn_colours.csv is missing required columns (after normalisation): %s",
    paste(missing_cols, collapse = ", ")
  ))
}

# Validate hex & drop bad rows
bad_hex <- which(!is_hex(pcn_colors$color))
if (length(bad_hex) > 0) {
  warning(sprintf(
    "Dropping %d rows with invalid hex codes in pcn_colors.csv: %s",
    length(bad_hex),
    paste(pcn_colors$PCN_NAME[bad_hex], collapse = "; ")
  ))
  pcn_colors <- pcn_colors[-bad_hex, , drop = FALSE]
}

# Warn on duplicates in the colour file for the same (Sub_ICB, PCN_NAME)
#dups <- pcn_colors |>
#  dplyr::count(Sub_ICB, PCN_NAME) |>
#  dplyr::filter(n > 1)
#if (nrow(dups) > 0) {
#  warning(
#    "Duplicate colour rows found in pcn_colors.csv for:\n",
#    paste(apply(dups, 1, function(r) paste0("  - ", r[["Sub_ICB"]], " / ", r[["PCN_NAME"]], " (", r[["n"]], ")")), collapse = "\n"),
#    "\nKeeping the first occurrence for each."
#  )
#  pcn_colors <- pcn_colors |>
#    dplyr::group_by(Sub_ICB, PCN_NAME) |>
#    dplyr::slice(1) |>
#    dplyr::ungroup()
#}

# ---- Join fixed colours by BOTH keys, keep only needed fields ----
pcn_palette_df <- pcn_map |>
  dplyr::left_join(
    pcn_colors[, c("Sub_ICB", "PCN_NAME", "color")],
    by = c("Sub_ICB", "PCN_NAME")
  )




# --- Helpers for grouped UI ---
safe_id <- function(x) gsub("[^A-Za-z0-9]", "_", x)

# Split PCNs by SubICB for building UI controls
subicb_groups <- split(pcn_map$PCN_NAME, pcn_map$Sub_ICB)

# Stable input IDs per group
group_ids <- setNames(safe_id(names(subicb_groups)), names(subicb_groups))



# ----------------------------------------------
# 4) Functions: union polygons + grouped legend
# ----------------------------------------------
# Union selected LSOAs into one polygon per PCN
build_union <- function(pcns, threshold) {
  if (length(pcns) == 0) return(NULL)
  
  # Filter once then group by PCN and union
  sub <- map_data |> 
    filter(PCN_NAME %in% pcns,
           !is.na(Cum_Perc),
           Cum_Perc <= threshold)
  
  if (nrow(sub) == 0) return(NULL)
  
  # Robust union by PCN (keeps sf class and avoids named vectors)
  
  polys <- sub |>
    group_by(PCN_NAME, Sub_ICB) |>
    summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
    sf::st_transform(4326)  # ensure coordinates are truly WGS84 lon/lat
  
  
  if (nrow(polys) == 0) return(NULL)
  polys
}



library(scales)  # for comma formatting

# Sum patients for the same filter used in build_union()
compute_pcn_patients <- function(pcns, threshold, patient_col = "Patients") {
  # Validate column exists
  if (!patient_col %in% names(map_data)) {
    warning(sprintf("Column '%s' not found in map_data; returning NA counts.", patient_col))
    return(
      map_data %>%
        st_drop_geometry() %>%
        filter(PCN_NAME %in% pcns,
               !is.na(Cum_Perc),
               Cum_Perc <= threshold) %>%
        distinct(PCN_NAME, Sub_ICB) %>%
        mutate(patients = NA_integer_)
    )
  }
  
  # Aggregate patient totals by PCN (respecting threshold and current filter)
  map_data %>%
    st_drop_geometry() %>%
    filter(PCN_NAME %in% pcns,
           !is.na(Cum_Perc),
           Cum_Perc <= threshold) %>%
    group_by(PCN_NAME, Sub_ICB) %>%
    summarise(patients = sum(.data[[patient_col]], na.rm = TRUE), .groups = "drop")
}



ensure_patients_col <- function(sf_obj) {
  if (!"patients" %in% names(sf_obj)) sf_obj$patients <- NA_integer_
  sf_obj
}



# Create grouped HTML legend 


make_grouped_legend <- function(pal_df) {
  groups <- split(pal_df, pal_df$Sub_ICB)
  
  container <- tags$div(
    style = "background:white;padding:10px 12px;border-radius:4px;
             box-shadow:0 1px 4px rgba(0,0,0,0.2);max-height:320px;overflow:auto;
             z-index:1000; position:relative;",   # position:relative for the handle
    tags$div(style="font-weight:600;margin-bottom:6px;", "PCNs by Sub‑ICB"),
    lapply(groups, function(df) {
      tags$div(
        style = "margin-bottom:8px;",
        tags$div(style="font-weight:600;margin-bottom:4px;", df$Sub_ICB[1]),
        tags$div(
          lapply(seq_len(nrow(df)), function(i) {
            tags$div(
              style = "display:flex;align-items:center;margin:2px 0;",
              tags$span(
                style = paste0("background:", df$color[i],
                               ";width:12px;height:12px;display:inline-block;",
                               "margin-right:6px;border:1px solid #555;")
              ),
              tags$span(df$PCN_NAME[i])
            )
          })
        )
      )
    }),
    # --- resize handle (bottom-right) ---
    tags$div(id = "pcn-legend-resize", style = "
        position:absolute; right:4px; bottom:4px; width:14px; height:14px;
        cursor:se-resize; opacity:0.75; border-right:2px solid #999; border-bottom:2px solid #999;
        border-radius:0;")  # diagonal corner look
  )
  
  HTML(as.character(container))
}



# --- West Midlands bbox and view helpers (EPSG:4326) ---

WEST_MIDLANDS_BBOX <- list(
  xmin = -3.20,  # west
  ymin = 51.95,  # south
  xmax = -1.30,  # east
  ymax = 53.10   # north
)

# Expand a bbox asymmetrically to create more space top-right (legend side)
expand_bbox_asymmetric <- function(bb, pad_bl = c(0.02, 0.02), pad_tr = c(0.25, 0.28)) {
  w <- bb[["xmax"]] - bb[["xmin"]]
  h <- bb[["ymax"]] - bb[["ymin"]]
  bb[["xmin"]] <- bb[["xmin"]] - w * pad_bl[1]
  bb[["ymin"]] <- bb[["ymin"]] - h * pad_bl[2]
  bb[["xmax"]] <- bb[["xmax"]] + w * pad_tr[1]
  bb[["ymax"]] <- bb[["ymax"]] + h * pad_tr[2]
  bb
}

# Fit to a bbox with asymmetric padding
fit_bounds_asymmetric <- function(map_id = "map", bb,
                                  pad_bl = c(0.02, 0.02),
                                  pad_tr = c(0.25, 0.28)) {
  bb2 <- expand_bbox_asymmetric(bb, pad_bl = pad_bl, pad_tr = pad_tr)
  leaflet::leafletProxy(map_id) %>%
    leaflet::fitBounds(
      lng1 = bb2[["xmin"]], lat1 = bb2[["ymin"]],
      lng2 = bb2[["xmax"]], lat2 = bb2[["ymax"]]
    )
  invisible(NULL)
}

# Fit to West Midlands with asymmetric padding so polygons sit bottom-left
fit_to_west_midlands <- function(map_id = "map",
                                 pad_bl = c(0.02, 0.02),
                                 pad_tr = c(0.25, 0.28)) {
  fit_bounds_asymmetric(map_id, WEST_MIDLANDS_BBOX, pad_bl = pad_bl, pad_tr = pad_tr)
}

# Fit to current layer bounds (snug), else fallback to WM
fit_to_layer_bounds <- function(map_id = "map", sf_layer,
                                pad_bl = c(0.02, 0.02),
                                pad_tr = c(0.25, 0.28)) {
  if (is.null(sf_layer) || nrow(sf_layer) == 0) {
    return(fit_to_west_midlands(map_id, pad_bl = pad_bl, pad_tr = pad_tr))
  }
  bb <- as.list(sf::st_bbox(sf::st_transform(sf_layer, 4326)))
  fit_bounds_asymmetric(map_id, bb, pad_bl = pad_bl, pad_tr = pad_tr)
}






# --- Helpers to improve visibility ---


# Compute luminance (0 = black, 1 = white) for hex colours
luminance_hex <- function(hex) {
  rgb <- grDevices::col2rgb(hex) / 255
  srgb_to_linear <- function(c) ifelse(c <= 0.03928, c/12.92, ((c + 0.055)/1.055)^2.4)
  lin <- srgb_to_linear(rgb)
  as.numeric(0.2126*lin[1,] + 0.7152*lin[2,] + 0.0722*lin[3,])
}

# Style sf polygons by luminance: stronger borders for light fills
style_polys <- function(sf_obj, color_col = "color") {
  stopifnot(color_col %in% names(sf_obj))
  sf_obj |>
    dplyr::mutate(
      lum = luminance_hex(.data[[color_col]]),
      stroke_weight = dplyr::case_when(
        lum >= 0.80 ~ 3.2,
        lum >= 0.70 ~ 2.4,
        TRUE        ~ 1.6
      ),
      stroke_col = dplyr::case_when(
        lum >= 0.80 ~ "#1F1F1F",
        lum >= 0.70 ~ "#2B2B2B",
        TRUE        ~ "#3A3A3A"
      ),
      fill_opacity = dplyr::case_when(
        lum < 0.25 ~ 0.35,
        lum < 0.45 ~ 0.40,
        lum < 0.65 ~ 0.45,
        TRUE       ~ 0.55
      )
    )
}

# Optional: add a subtle glow underlay to separate polygons from basemap
#add_glow_underlay <- function(map_id, sf_obj) {
#  leaflet::leafletProxy(map_id) %>%
#    leaflet::addPolygons(
#      data = sf_obj,
#      fillColor   = "transparent",
#      fillOpacity = 0,
#      color       = "#000000",
#      weight      = 3.5,       # slightly thicker than main stroke
#      opacity     = 0.25,
#      smoothFactor = 0
#    )
#}

# Swap dark glow underlay for a white halo instead

add_white_halo <- function(map_id, sf_obj) {
  leaflet::leafletProxy(map_id) %>%
    leaflet::addPolygons(
      data        = sf_obj,
      fill        = FALSE,       # ensure no fill is drawn
      fillColor   = "transparent",
      fillOpacity = 0,
      color       = "#FFFFFF",
      weight      = 7.0,         # thicker than your main stroke so it peeks out
      opacity     = 0.95,
      smoothFactor = 0
    )
}



# -------------------------------
# 5) UI
# -------------------------------

ui <- fluidPage(
  theme = bslib::bs_theme(version = 5), # Bootstrap 5
  
  titlePanel(
    div(
      style = "display:flex; align-items:center; gap:40px;",
      # Bigger, crisp logo; use SVG if available
      tags$img(src = "assets/logo.png",
               alt = "NHS Midlands & Lancashire CSU",
               style = "height:48px; width:auto; display:block;"
      ),
      tags$span(
        "SSOT and STW PCN Boundaries - Jan 2026",
        style = "font-weight:600;"
      )
    )
  ),
  
  
  # DEBUG BANNER — place right after titlePanel
  #tags$div(
  #  style = "background:#fff3cd;color:#664d03;border:1px solid #ffecb5;padding:6px;margin-bottom:8px;",
  #  paste("Working dir:", getwd()),
  #  tags$br(),
  #  paste("www/ contents:", paste(list.files("www"), collapse = ", ")),
  #  tags$br(),
  #  paste("logo.png exists:", file.exists("www/logo.png"))
  #),
  
  
  # Add Information box
  
  actionLink("info_btn", label = NULL, icon = icon("circle-info"),
             class = "btn btn-link", 
             style = "position:absolute; right:18px; top:8px; font-size:18px;"),
  
  
  
  
  tags$style(HTML("

/* =========================
   Sub‑ICB PCN selector blocks
   ========================= */

.subicb-block {
  margin: 6px 0 14px 0;
  padding-bottom: 10px;
  border-bottom: 1px solid #e6e6e6;
}

/* Sub‑ICB header checkbox */
.subicb-block > .form-group.shiny-input-container.checkbox {
  margin-bottom: 6px;
}
.subicb-block > .form-group.shiny-input-container.checkbox label,
.subicb-block > .form-group.shiny-input-container .form-check-label {
  font-weight: 700 !important;     /* forces bold */
  font-size: 10px !important;      /* header slightly larger than the PCNs */
  color: #1f1f1f;
}

/* PCN list indented under Sub‑ICB */
.subicb-block .shiny-input-checkboxgroup {
  margin-top: 0px;
  padding-left: 18px;
}

/* Tighter PCN spacing */
.subicb-block .shiny-input-checkboxgroup .checkbox {
  margin-top: 2px;
  margin-bottom: 2px;
}

/* =========================
   Sidebar base typography (replaces .sidebarPanel section)
   ========================= */

.app-sidebar {
  font-size: 14px;      /* change to 8px if you want smaller */
  line-height: 1.4;
}

.app-sidebar label,
.app-sidebar .control-label,
.app-sidebar .form-label,
.app-sidebar .form-check-label,
.app-sidebar .help-block,
.app-sidebar .shiny-text-output,
.app-sidebar .shiny-html-output {
  font-size: inherit !important;
  line-height: inherit !important;
}

.app-sidebar input[type='checkbox']{
  transform: scale(0.9);
  transform-origin: left center;
}

/* Scrollable sidebar content */
.sidebar-scroll {
  max-height: 75vh;
  overflow-y: auto;
  padding-right: 8px;
}

/* Sticky footer button */
.sidebar-sticky-footer {
  position: sticky;
  bottom: 0;
  background: #fff;
  padding: 8px 0;
  box-shadow: 0 -4px 8px rgba(0,0,0,0.05);
}

")),
  


  
  tags$style(HTML("
/* =========================
   Floating legend panel
   ========================= */

.legend-panel {
  background: white;
  padding: 10px 12px;
  border-radius: 4px;
  box-shadow: 0 1px 4px rgba(0,0,0,0.2);
  z-index: 1000;
  max-height: 320px;
  overflow: auto;
  resize: both;
  font-size: 12px;
  line-height: 1.25;
}
")),
  
  
  
  # Draggable + resizable legend overlay
  absolutePanel(
    id = "legend_panel",
    class = "legend-panel",
    draggable = TRUE,               # Shiny makes it draggable
    top = 80, right = 20,           # initial position (like top-right control)
    width = 340,                    # initial width
    # height can be auto; we'll cap/max-height via CSS
    uiOutput("legend_ui")           # we'll fill this from server
  ),
  
  
  
  
  sidebarLayout(
    sidebarPanel(
      class = "app-sidebar",
      
      # Scrollable content (everything except the update button)
      div(class = "sidebar-scroll",
          # Global on/off
          checkboxInput("all_toggle", "Select all PCNs", TRUE),
          
          # Grouped PCN selectors built dynamically per SubICB
          uiOutput("grouped_pcn_ui"),
          
          tags$hr(),
          sliderInput("threshold", "Proportion of Patients (%):",
                      min = 50, max = 98, value = 98, step = 1),
          # checkboxInput("show_debug", "Show diagnostics in console", value = FALSE),
          
          tags$hr(),
          # Display note requested
          helpText(
            "Note: PCN catchments can be selected to show where at least 50% and up to 98% of their patients reside, based on data at Oct 2025.",
            "In the legend, ABC PCN and Rugeley & Great Haywood PCN are assigned to the Sub-ICB with the majority of their practices, although their catchments are defined on their total patient list."
          )
      ),
      # Sticky footer with the Update Map button (always visible)
      div(class = "sidebar-sticky-footer",
          actionButton("update", "Update Map", class = "btn btn-primary", width = "100%")
      )
    ),
    mainPanel(
      leafletOutput("map", height = 1000)
    )
  )
)


# -------------------------------
# 6) Server (cleaned)
# -------------------------------
server <- function(input, output, session) {
  
  # 6a) Base map
  
  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(worldCopyJump = FALSE)) %>%
      addProviderTiles("CartoDB.Positron", options = providerTileOptions(opacity = 1)) %>%
      setView(lng = -2.18, lat = 53.00, zoom = 6)  # Birmingham centre-ish
    # Optionally restrict panning to England:
    # setMaxBounds(lng1 = -6.0, lat1 = 49.9, lng2 = 2.1, lat2 = 55.9)
    
    # Optional: restrict panning to West Midlands bounds
    # setMaxBounds(lng1 = WEST_MIDLANDS_BBOX$xmin, lat1 = WEST_MIDLANDS_BBOX$ymin,
    #              lng2 = WEST_MIDLANDS_BBOX$xmax, lat2 = WEST_MIDLANDS_BBOX$ymax)
    
    
  })
  
  
  # 6b) Build grouped PCN UI (one bsCollapsePanel per Sub-ICB)
  
  output$grouped_pcn_ui <- renderUI({
    tagList(
      lapply(names(subicb_groups), function(g) {
        gid <- group_ids[[g]]
        
        tags$div(
          class = "subicb-block",
          
          # This becomes the "heading"
          checkboxInput(
            inputId = paste0("toggle_", gid),
            label   = paste0("Select all PCNs in ", g),
            value   = TRUE
          ),
          
          # PCNs underneath
          checkboxGroupInput(
            inputId  = paste0("pcn_", gid),
            label    = NULL,
            choices  = unname(sort(subicb_groups[[g]])),
            selected = unname(sort(subicb_groups[[g]]))
          )
        )
      })
    )
  })
  
  
  # 6c) Reactively gather all selected PCNs across groups
  selected_pcns <- reactive({
    unlist(lapply(names(subicb_groups), function(g) {
      gid <- group_ids[[g]]
      input[[paste0("pcn_", gid)]]
    }), use.names = FALSE)
  })
  
  
  
  # 6c1) - added to keep each Sub-ICB select-all checkbox in sync with manual PCN selections
  lapply(names(subicb_groups), function(g) {
    gid <- group_ids[[g]]
    
    observeEvent(input[[paste0("pcn_", gid)]], {
      chosen <- input[[paste0("pcn_", gid)]] %||% character(0)
      all_in_group <- subicb_groups[[g]]
      
      updateCheckboxInput(
        session,
        paste0("toggle_", gid),
        value = length(chosen) == length(all_in_group)
      )
    }, ignoreInit = TRUE)
  })
  
  
  
  # 6d) Per-SubICB toggle logic (binary select/clear per group)
  lapply(names(subicb_groups), function(g) {
    gid <- group_ids[[g]]
    
    observeEvent(input[[paste0("toggle_", gid)]], {
      sel <- if (isTRUE(input[[paste0("toggle_", gid)]])) subicb_groups[[g]] else character(0)
      updateCheckboxGroupInput(session, paste0("pcn_", gid), selected = sort(sel))
    }, ignoreInit = TRUE)
  })
  
  # 6e) Global select-all toggle: affects every group
  observeEvent(input$all_toggle, {
    lapply(names(subicb_groups), function(g) {
      gid <- group_ids[[g]]
      sel <- if (isTRUE(input$all_toggle)) subicb_groups[[g]] else character(0)
      updateCheckboxGroupInput(session, paste0("pcn_", gid), selected = sort(sel))
      updateCheckboxInput(session, paste0("toggle_", gid), value = isTRUE(input$all_toggle))
    })
  }, ignoreInit = TRUE)
  
  
  # 6f) Information Box
  
  observeEvent(input$info_btn, {
    showModal(modalDialog(
      title = "About this app",
      easyClose = TRUE,
      footer = modalButton("Close"),
      size = "l",
      tagList(
        tags$p("This application allows users to visualise PCN boundaries based on the proportion of patients registered at January 2026. For each Lower Layer Super Output Area (LSOA), the total number of patients registered at a PCN are summed.  The LSOAs are then sorted in descending order of patients, so LSOAs with the highest density of patients are displayed first. Cumulative patient counts are calculated as well as the percentage of the total PCN population those represent."),  
        tags$p("The cumulative percentage is used in the threshold slider.  Users can select the proportion of patients they wish to view, and the underlying data at LSOA is merged to produce a single boundary per PCN. The number of patients within the boundary is displayed in the map labels."),
        tags$ul(
          tags$li("Coverage threshold: adjustable to show where at least 50% and up to 98% of patients reside."),
          tags$li("Data source: ",
                  tags$a(
                    href = "https://digital.nhs.uk/data-and-information/publications/statistical/patients-registered-at-a-gp-practice",
                    target = "_blank",
                    "NHS England GP Registered Patients"
                  ),
                  "(Jan 2026)."
          ))
      ),
      tags$hr(),
      tags$p(
        strong("Contact: "), "Sarah Glover - NHS Midlands & Lancashire CSU — BI Analytics Hub.",
        tags$br(),
        "Email: ", tags$a(href="mailto:analytics.mlcsu@nhs.net", "analytics.mlcsu@nhs.uk")
      )
    )
    )
  })
  
  
  # 6g) Initial render: all PCNs at 98% + grouped legend
  observeEvent(TRUE, {
    pcns0  <- sort(unique(map_data$PCN_NAME))
    polys0 <- build_union(pcns0, 98)
    
    leafletProxy("map") %>%
      clearControls()
    
    if (!is.null(polys0)) {
      
      
      
      polys0 <- polys0 %>%
        dplyr::left_join(pcn_palette_df, by = c("Sub_ICB", "PCN_NAME"))
      
      
      # Join patient totals (safe even if column missing)
      totals0 <- compute_pcn_patients(pcns0, threshold = 98, patient_col = "PCN_Pat")
      polys0 <- polys0 %>%
        dplyr::left_join(totals0, by = c("Sub_ICB", "PCN_NAME")) %>%
        ensure_patients_col() %>%
        style_polys(color_col = "color")
      
      
      # Legend based on the coloured layer (now into the overlay panel)
      output$legend_ui <- renderUI({
        make_grouped_legend(
          polys0 %>% sf::st_drop_geometry() %>% dplyr::select(Sub_ICB, PCN_NAME, color)
        )
      })
      
      
      
      # Optional: white halo under outlines for extra separation
      add_white_halo("map", polys0)
      
      leafletProxy("map") %>%
        addPolygons(
          data        = polys0,
          fillColor   = ~color,
          fillOpacity = 0.6,
          color       = ~stroke_col,
          weight      = ~stroke_weight,
          opacity     = 1,
          label = ~sprintf(
            "%s\nPatients: %s",
            PCN_NAME, ifelse(is.na(patients), "n/a", scales::comma(patients))
          ),
          labelOptions = labelOptions(
            noHide = FALSE,
            direction = "auto",
            style = list(
              "font-weight" = "600",
              "color" = "#1f1f1f",
              "background" = "rgba(255,255,255,0.9)",
              "border" = "1px solid #ccc",
              "padding" = "4px 6px"
            )
          ),
          popup = ~paste0(
            "<b>PCN:</b> ", PCN_NAME,
            "<br><b>Sub‑ICB:</b> ", Sub_ICB,
            ifelse(is.na(patients), "", paste0("<br><b>Patients:</b> ", scales::comma(patients)))
          ),
          highlightOptions = highlightOptions(weight = 4, color = "#000000", opacity = 1, bringToFront = TRUE)
        )
      
      fit_to_layer_bounds("map", polys0)
      
    } else {
      fit_to_west_midlands("map")
    }
  }, once = TRUE)
  
  
  # 6h) Update map on button click
  
  observeEvent(input$update, {
    req(length(selected_pcns()) > 0)
    
    # Build union for selected PCNs at chosen threshold
    polys <- build_union(selected_pcns(), input$threshold)
    
    leafletProxy("map") %>%
      clearShapes() %>%
      clearControls()
    
    if (!is.null(polys)) {
      
      
      polys <- polys %>%
        dplyr::left_join(pcn_palette_df, by = c("Sub_ICB", "PCN_NAME"))
      
      
      # Join patient totals and style
      totals <- compute_pcn_patients(
        selected_pcns(),
        threshold   = input$threshold,
        patient_col = "PCN_Pat"     # change here if your column name differs
      )
      
      polys <- polys %>%
        dplyr::left_join(totals, by = c("Sub_ICB", "PCN_NAME")) %>%
        ensure_patients_col() %>%
        style_polys(color_col = "color")
      
      
      # Update the overlay legend when selection/threshold changes
      output$legend_ui <- renderUI({
        make_grouped_legend(
          polys %>% sf::st_drop_geometry() %>% dplyr::select(Sub_ICB, PCN_NAME, color)
        )
      })
      
      
      # Optional: white halo
      add_white_halo("map", polys)
      
      leafletProxy("map") %>%
        addPolygons(
          data        = polys,
          fillColor   = ~color,
          fillOpacity = 0.6,
          color       = ~stroke_col,
          weight      = ~stroke_weight,
          opacity     = 1,
          label = ~sprintf(
            "%s\nPatients: %s",
            PCN_NAME, ifelse(is.na(patients), "n/a", scales::comma(patients))
          ),
          labelOptions = labelOptions(
            noHide = FALSE,
            direction = "auto",
            style = list(
              "font-weight" = "600",
              "color" = "#1f1f1f",
              "background" = "rgba(255,255,255,0.9)",
              "border" = "1px solid #ccc",
              "padding" = "4px 6px"
            )
          ),
          popup = ~paste0(
            "<b>PCN:</b> ", PCN_NAME,
            "<br><b>Sub‑ICB:</b> ", Sub_ICB,
            ifelse(is.na(patients), "", paste0("<br><b>Patients:</b> ", scales::comma(patients)))
          ),
          highlightOptions = highlightOptions(weight = 4, color = "#000000", opacity = 1, bringToFront = TRUE)
        )
      
      fit_to_layer_bounds("map", polys)
      
    } else {
      showNotification("No polygons for the selected PCNs/threshold. Try 98% or check data.", type = "message")
      fit_to_west_midlands("map")
    }
  })
}


# -------------------------------
# 7) Launch
# -------------------------------
options(shiny.launch.browser = TRUE)
shinyApp(ui, server)



