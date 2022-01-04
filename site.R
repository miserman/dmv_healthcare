library(community)
page_head(
  title = "Capital Region Healthcare Facility Access",
  description = "Comparison of floating catchment area ratios of health care facilities in the DMV area."
)
page_navbar(
  title = "Health Care Facility Access",
  input_button("Reset", "reset_selection", "reset.selection", class = "btn-link"),
  list(
    name = "Settings",
    backdrop = "false",
    items = list(
      input_switch("Dark Theme", id = "settings.theme_dark"),
      input_select(
        "Color Palette", options = "palettes", default = "rdylbu7", id = "settings.palette",
        floating_label = FALSE
      ),
      input_switch(
        "Color by Order", id = "settings.color_by_order",
        title = paste(
          "Switch from coloring by value to coloring by sorted index.",
          "This may help differentiate regions with similar values."
        )
      ),
      input_switch("Hide URL Settings", id = "settings.hide_url_parameters"),
      input_number("Digits", "settings.digits", min = 0, max = 6, floating_label = FALSE),
      input_select(
        "Summary Level", options = c("dataset", "all"), default = "dataset",
        display = c("All Regions", "Selected Region"), id = "settings.summary_selection",
        floating_label = FALSE,
        title = paste(
          "Determins which regions are included in summaries for box-plots and color scaling;",
          "All-Regions are county-wide, Selected Region Types are filtered by the Region Types input, and",
          "Selected Region are filtered by region selection."
        )
      ),
      '<p class="section-heading">Map Options</p>',
      input_switch("Show Background Shapes", id = "settings.background_shapes"),
      input_button("Clear Settings", "reset_storage", "clear_storage", class = "btn-danger footer")
    )
  ),
  list(
    name = "About",
    items = list(
      page_text(c(
        paste0(
          "This site was made by the [Social and Decision Analytics Division]",
          "(https://biocomplexity.virginia.edu/institute/divisions/social-and-decision-analytics)",
          " of the [Biocomplexity Institute](https://biocomplexity.virginia.edu)."
        ),
        "View its source on [GitHub](https://github.com/uva-bi-sdad/dmv_healthcare).",
        "Credits",
        paste(
          "Built in [R](https://www.r-project.org) with the",
          "[community](https://uva-bi-sdad.github.io/community) package, using these resources:"
        )
      ), class = c("", "", "h5")),
      output_credits()
    )
  )
)
page_menu(
  page_section(
    type = "col",
    page_section(
      type = "row",
      wraps = "col",
      input_select(
        "County", options = "ids", dataset = "counties", dataview = "primary_view",
        id = "selected_county", reset_button = TRUE
      ),
      input_select(
        "Census Tract", options = "ids", dataset = "tracts", dataview = "primary_view",
        id = "selected_tract", reset_button = TRUE
      ),
      conditions = c("", "selected_county")
    )
  ),
  page_section(
    type = "col",
    page_section(
      type = "row",
      wraps = "col",
      input_select(
        "Y Variable (mapped)", options = "variables",
        default = "access_3sfca", depends = "shapes",  id = "selected_y"
      ),
      input_select(
        "X Variable", options = "variables",
        default = "population", depends = "shapes", id = "selected_x"
      )
    )
  ),
  position = "top",
  default_open = TRUE
)
input_variable("shapes", list(
  "selected_county && !selected_tract" = "tracts",
  "selected_tract" = "blockgroups"
), "counties")

input_variable("region_select", list(
  "shapes == tracts" = "selected_tract"
), "selected_county")

input_variable("selected_region", list(
  "selected_tract" = "selected_tract"
), "selected_county")

input_dataview(
  "primary_view",
  y = "selected_y",
  dataset = "shapes",
  ids = "selected_region"
)
page_section(
  type = "col",
  output_text(c(
    "(National Capital Region)[r selected_county]",
    "? > {selected_county}[r selected_tract]",
    "? > {selected_tract}"
  )),
  output_text(list(
    "default" = "National Capital Region Counties",
    "selected_county" = "{selected_county} Census Tracts",
    "selected_tract" = "{selected_tract} Census Block Groups"
  ), tag = "h1", class = "text-center"),
  page_section(
    type = "row",
    wraps = "col",
    sizes = c(NA, 4),
    output_map(
      lapply(c("counties", "tracts", "blockgroups"), function(s) list(
        name = s,
        url = paste0("../dmv_healthcare/docs/data/", s, ".geojson")
      )),
      dataview = "primary_view",
      click = "region_select",
      id = "main_map",
      subto = "main_plot",
      options = list(
        attributionControl = FALSE,
        scrollWheelZoom = FALSE,
        center = c(38.938, -77.315),
        zoom = 7,
        height = "430px"
      ),
      background_shapes = "tracts",
      tiles = list(
        light = list(url = "https://stamen-tiles-{s}.a.ssl.fastly.net/toner-lite/{z}/{x}/{y}{r}.png"),
        dark = list(url = "https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png")
      ),
      attribution = list(
        list(
          name = "Stamen toner-light",
          url = "https://stamen.com",
          description = "Light-theme map tiles by Stamen Design"
        ),
        list(
          name = "CARTO Dark Matter",
          url = "https://carto.com/attributions",
          description = "Dark-theme map tiles by CARTO"
        ),
        list(
          name = "OpenStreetMap",
          url = "https://www.openstreetmap.org/copyright"
        )
      )
    ),
    page_section(
      type = "d-flex flex-column col align-items-end compact",
      output_info(
        title = "variables.short_name",
        body = "variables.source",
        dataview = "primary_view",
        id = "variable_info_pane",
      ),
      page_section(
        wraps = "row",
        output_info(
          title = "features.name",
          default = c(title = "National Capital Region", body = "Hover over or select a region for more information."),
          dataview = "primary_view",
          subto = c("main_map", "main_plot")
        ),
        output_info(
          body = c(
            "variables.long_name" = "selected_y",
            "variables.statement"
          ),
          row_style = c("stack", "table"),
          dataview = "primary_view",
          subto = c("main_map", "main_plot"),
          variable_info = FALSE
        )
      ),
      output_legend("settings.palette", "Below", "Region Median", "Above"),
      wraps = c("row", "row mb-auto", "row")
    )
  ),
  page_section(
    type = "row",
    wraps = "col",
    sizes = c(7, 5),
    page_tabgroup(
      list(
        name = "Plot",
        output_plot(
          x = "selected_x", y = "selected_y", dataview = "primary_view",
          click = "region_select", subto = "main_map", id = "main_plot",
          options = list(
            layout = list(
              showlegend = FALSE,
              xaxis = list(fixedrange = TRUE),
              yaxis = list(fixedrange = TRUE, zeroline = FALSE)
            ),
            config = list(modeBarButtonsToRemove = c("select2d", "lasso2d", "sendDataToCloud"))
          )
        )
      ),
      list(
        name = "Data",
        output_table(
          dataview = "primary_view", wide = FALSE,
          features = c(ID = "id", Name = "name"),
          options = list(
            scrollY = 400,
            rowGroup = list(dataSrc = "features.name"),
            columnDefs = list(list(targets = "features.name", visible = FALSE)),
            buttons = c('copy', 'csv', 'excel', 'print'),
            dom = "<'row't><'row'<'col-sm'B><'col'f>>"
          )
        )
      )
    ),
    output_table("selected_y", dataview = "primary_view", options = list(
      info = FALSE,
      searching = FALSE
    ))
  )
)
site_build('../dmv_healthcare')
