library(community)
page_head(
  title = "Capital Region Healthcare Facility Access",
  icon = "https://raw.githubusercontent.com/uva-bi-sdad/catchment/main/logo.svg",
  description = "Comparison of floating catchment area ratios of health care facilities in the DMV area."
)
page_navbar(
  title = "Healthcare Facility Access",
  logo = "https://raw.githubusercontent.com/uva-bi-sdad/catchment/main/logo.svg",
  input_button("Reset", "reset_selection", "reset.selection", class = "btn-link"),
  list(
    name = "Settings",
    backdrop = "false",
    items = list(
      input_switch("Dark Theme", id = "settings.theme_dark"),
      input_select("Color Palette", options = "palettes", id = "settings.palette", floating_label = FALSE),
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
        paste(
          "The [Catchment](https://uva-bi-sdad.github.io/catchment) package was used to calculate",
          "floating catchment area ratios."
        ),
        paste(
          "The [case study](https://uva-bi-sdad.github.io/catchment/articles/casestudy-dmv.html) article walkes",
          "through the data collection and calculation of floating catchment area ratios."
        ),
        "Credits",
        paste(
          "Built in [R](https://www.r-project.org) with the",
          "[community](https://uva-bi-sdad.github.io/community) package, using these resources:"
        )
      ), class = c(character(4), "h5")),
      output_credits()
    )
  )
)

output_text("National Capital Region", tag = "h1", class = "text-center")

input_variable("shapes_a", list(
  "county_a && !tract_a" = "tract",
  "tract_a" = "block_group"
), "county")
input_variable("region_select_a", list(
  "shapes_a == tract" = "tract_a"
), "county_a")
input_variable("region_a", list(
  "tract_a" = "tract_a"
), "county_a")
input_dataview(
  "view_a",
  y = "variable_a",
  dataset = "shapes_a",
  ids = "region_a"
)

input_variable("shapes_b", list(
  "county_b && !tract_b" = "tract",
  "tract_b" = "block_group"
), "county")
input_variable("region_select_b", list(
  "shapes_b == tract" = "tract_b"
), "county_b")
input_variable("region_b", list(
  "tract_b" = "tract_b"
), "county_b")
input_dataview(
  "view_b",
  y = "variable_b",
  dataset = "shapes_b",
  ids = "region_b"
)


output_info(
  title = "features.name",
  body = c(
    "variables.long_name" = "variable_a",
    "variables.statement"
  ),
  default = c(title = ""),
  row_style = c("stack", "table"),
  dataview = "view_a",
  subto = c("map_a", "plot_a"),
  variable_info = FALSE,
  floating = TRUE
)

output_info(
  title = "features.name",
  body = c(
    "variables.long_name" = "variable_b",
    "variables.statement"
  ),
  default = c(title = ""),
  row_style = c("stack", "table"),
  dataview = "view_b",
  subto = c("map_b", "plot_b"),
  variable_info = FALSE,
  floating = TRUE
)

page_section(
  wraps = "col",
  page_section(
    output_text(list(
      "default" = "All Counties",
      "county_a" = "{county_a} Census tract",
      "tract_a" = "{tract_a} Block Groups"
    ), tag = "h2", class = "text-center"),
    page_section(
      type = "row",
      wraps = "col",
      page_section(
        page_section(
          type = "row form-row",
          input_select(
            "County", options = "ids", dataset = "county", dataview = "view_a",
            id = "county_a", reset_button = TRUE
          ),
          input_select(
            "Census Tract", options = "ids", dataset = "tract", dataview = "view_a",
            id = "tract_a", reset_button = TRUE
          )
        ),
        input_select(
          "Variable A", options = "variables",
          default = "doctors_3sfca", depends = "shapes_a",  id = "variable_a"
        ),
        output_info(
          title = "variables.short_name",
          body = "variables.source",
          dataview = "view_a"
        )
      ),
      output_map(
        lapply(
          list(c("counties", "county"), c("tracts", "tract"), c("blockgroups", "block_group")),
          function(s) list(
            name = s[2],
            url = paste0("../dmv_healthcare/docs/data/", s[1], ".geojson")
          )
        ),
        dataview = "view_a",
        click = "region_select_a",
        id = "map_a",
        subto = "plot_a",
        options = list(
          attributionControl = FALSE,
          scrollWheelZoom = FALSE,
          center = c(38.938, -77.315),
          zoom = 7,
          height = "400px"
        ),
        background_shapes = "tract",
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
      )
    ),
    output_plot(
      x = "selected_x", y = "variable_a", dataview = "view_a",
      click = "region_select_a", subto = "map_a", id = "plot_a",
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
  page_section(
    output_text(list(
      "default" = "All Counties",
      "county_b" = "{county_b} Census tract",
      "tract_b" = "{tract_b} Block Groups"
    ), tag = "h2", class = "text-center"),
    page_section(
      wraps = "col",
      page_section(
        page_section(
          type = "row form-row",
          input_select(
            "County", options = "ids", dataset = "county", dataview = "view_b",
            id = "county_b", reset_button = TRUE
          ),
          input_select(
            "Census Tract", options = "ids", dataset = "tract", dataview = "view_b",
            id = "tract_b", reset_button = TRUE
          )
        ),
        input_select(
          "Variable B", options = "variables",
          default = "doctors_3sfca_step", depends = "shapes_b",  id = "variable_b"
        ),
        output_info(
          title = "variables.short_name",
          body = "variables.source",
          dataview = "view_b"
        )
      ),
      output_map(
        lapply(
          list(c("counties", "county"), c("tracts", "tract"), c("blockgroups", "block_group")),
          function(s) list(
            name = s[2],
            url = paste0("../dmv_healthcare/docs/data/", s[1], ".geojson")
          )
        ),
        dataview = "view_b",
        click = "region_select_b",
        id = "map_b",
        subto = "plot_b",
        options = list(
          attributionControl = FALSE,
          scrollWheelZoom = FALSE,
          center = c(38.938, -77.315),
          zoom = 7,
          height = "400px"
        ),
        background_shapes = "tract",
        tiles = list(
          light = list(url = "https://stamen-tiles-{s}.a.ssl.fastly.net/toner-lite/{z}/{x}/{y}{r}.png"),
          dark = list(url = "https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png")
        )
      )
    ),
    output_plot(
      x = "selected_x", y = "variable_b", dataview = "view_b",
      click = "region_select_b", subto = "map_b", id = "plot_b",
      options = list(
        layout = list(
          showlegend = FALSE,
          xaxis = list(fixedrange = TRUE),
          yaxis = list(fixedrange = TRUE, zeroline = FALSE)
        ),
        config = list(modeBarButtonsToRemove = c("select2d", "lasso2d", "sendDataToCloud"))
      )
    )
  )
)

input_select(
  "X Variable", options = "variables",
  default = "population", depends = "shapes", id = "selected_x"
)
output_legend("settings.palette", "Below", "Region Median", "Above")

site_build('../dmv_healthcare')
