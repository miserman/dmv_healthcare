# specify the main directory
maindir <- "../dmv_healthcare/docs/data/"

# make a directory to store working data files in
oridir <- paste0(maindir, "original/")
dir.create(oridir, FALSE, TRUE)

#
# provider information
#

# download primary care provider information from the Centers for Medicare & Medicaid Services
# https://data.cms.gov

## set up the url to download from the API
providers_url <- paste0(
  "https://data.cms.gov/data-api/v1/dataset/", # base url
  "a399e5c1-1cd1-4cbe-957f-d2cc8fe5d897/data/", # dataset id, which specifies year
  "?filter[region][condition][path]=Rndrng_Prvdr_State_Abrvtn", # filter: state %in% c("DC", "MD", "VA")
  "&filter[region][condition][operator]=IN",
  "&filter[region][condition][value][1]=DC",
  "&filter[region][condition][value][2]=MD",
  "&filter[region][condition][value][3]=VA",
  "&offset=" # pagination, as only 1000 rows are returned at a time
)

## retrieve the data in steps
providers <- list()
offset <- 0
while(length(raw <- jsonlite::read_json(paste0(providers_url, offset)))){
  cat("\rdownloading file", offset / 1000 + 1, "     ")
  providers[[length(providers) + 1]] <- do.call(rbind, lapply(raw, unlist))
  offset <- offset + 1000
}
providers <- as.data.frame(do.call(rbind, providers))

## focus only on medical doctors who aren't also dentists
providers <- providers[
  grepl("\\bm\\.?d(?:\\W|$)", providers$Rndrng_Prvdr_Crdntls, TRUE) &
    !grepl("d\\.?d\\.?s", providers$Rndrng_Prvdr_Crdntls, TRUE),
]

## format address strings, and get the coordinates of each unique entry
address_parts <- c("Rndrng_Prvdr_St1", "Rndrng_Prvdr_City", "Rndrng_Prvdr_State_Abrvtn", "Rndrng_Prvdr_Zip5")

providers$Rndrng_Prvdr_St1 <- sub("([nesw])\\.([nesw])\\.*", "\\1\\2", providers$Rndrng_Prvdr_St1, TRUE)
providers$Rndrng_Prvdr_St1 <- sub("\\s*#\\s.*$", "", providers$Rndrng_Prvdr_St1)
providers$Rndrng_Prvdr_St1 <- sub("\\s*([nesw]{2})\\s.*$", " \\1", providers$Rndrng_Prvdr_St1, TRUE)
providers$address <- do.call(paste, c(unname(as.list(providers[, address_parts])), sep = ", "))

## collapse to locations based on address
vars <- c("address", "Rndrng_Prvdr_Gndr", grep("^(tot|drug|med|bene)_", colnames(providers), TRUE, value = TRUE))
vars <- vars[!vars %in% c("Drug_Sprsn_Ind", "Med_Sprsn_Ind")]
provider_locations <- do.call(rbind, lapply(unique(providers$addresses), function(a){
  d <- providers[providers$addresses == a, vars]
  d[d == ""] <- NA
  as.data.frame(list(
    address = a,
    doctors = nrow(d),
    prop_women = mean(d$Rndrng_Prvdr_Gndr == "F"),
    as.list(colMeans(matrix(
      as.numeric(as.matrix(d[, -(1:2)])), nrow(d),
      dimnames = list(NULL, vars[-(1:2)])
    ), na.rm = TRUE))
  ))
}))
provider_locations[is.na(provider_locations)] <- NA

## geocode addresses; takes a while
library(parallel)
cl <- makeCluster(detectCores() - 2)
address_coords <- do.call(rbind, parLapply(cl, provider_locations$address, function(a){
  coords <- tidygeocoder::geo(a, progress_bar = FALSE, quiet = TRUE, method = "arcgis")
  if(is.na(coords$long)) coords <- tidygeocoder::geo(a, progress_bar = FALSE, quiet = TRUE)
  coords
}))
stopCluster(cl)

### add to provider locations
provider_locations[, c("lat", "long")] <- address_coords[, -1]

## remove provider locations that could not be geocoded
provider_locations <- provider_locations[!is.na(provider_locations$long),]

## make unique IDs for each provider location
provider_locations$id <- paste0("l", seq_len(nrow(provider_locations)))

## save provider locations dataset
write.csv(provider_locations, paste0(maindir, "providers.csv"), row.names = FALSE)

#
# population information
#

# define counties that are part of the capital region
dmv_counties <- list(
  dc = "District of Columbia",
  md = c("Prince George's", "Montgomery"),
  va = c("Arlington", "Loudoun", "Fairfax", "Alexandria")
)
data <- list()
shapes <- list()

# download / load
for(state in c("dc", "md", "va")){
  # shapes
  counties <- download_census_shapes(oridir, state, "county", paste0(state, "_counties"))
  tracts <- download_census_shapes(oridir, state, "tract", paste0(state, "_tracts"))
  blockgroups <- download_census_shapes(oridir, state, "bg", paste0(state, "_blockgroups"))
  
  ## store subsets to combine later
  counties <- counties[counties$NAME %in% dmv_counties[[state]],]
  counties[counties$NAME == "Fairfax", "NAME"] <- c("Fairfax City", "Fairfax")
  shapes[[state]] <- list(
    counties = counties,
    tracts = tracts[substr(tracts$GEOID, 1, 5) %in% counties$GEOID,],
    blockgroups = blockgroups[substr(blockgroups$GEOID, 1, 5) %in% counties$GEOID,]
  )
  
  # population data
  data[[state]] <- download_census_population(
    oridir, state, 2019, include_margins = TRUE, include_commutes = TRUE,
    counties = counties$GEOID, verbose = TRUE
  )
}

## create and save combined shapes
library(sf)
library(rmapshaper)

for(level in names(shapes$dc)){
  st_write(
    ms_simplify(do.call(rbind, lapply(shapes, "[[", level)), keep_shapes = TRUE),
    paste0(maindir, paste0(level, ".geojson"))
  )
}

## create and save final files
data_combined <- do.call(rbind, lapply(names(data), function(state){
  d <- data[[state]]$estimates
  s <- shapes[[state]]$blockgroups
  rownames(s) <- s$GEOID
  total <- d$TOTAL.POPULATION_Total
  total[total == 0] <- 1
  data.frame(
    GEOID = d$GEOID,
    population = d$TOTAL.POPULATION_Total,
    percent_female = d$SEX.BY.AGE_Female_Female / total * 100,
    percent_white = d$RACE_Total_White.alone / total * 100,
    percent_over_49 = rowSums(d[, grep("[5-8][05]", colnames(d))]) / total * 100,
    st_coordinates(st_centroid(st_geometry(s[d$GEOID,])))
  )
}))
write.csv(data_combined, paste0(maindir, "data.csv"), row.names = FALSE)

### get travel times between each included block group and provider location
### would need to be split if server is default
library(osrm)
traveltimes <- osrmTable(
  src = data_combined[, c("GEOID", "X", "Y")],
  dst = provider_locations[, c("id", "long", "lat")],
  osrm.server = Sys.getenv("OSRM_SERVER")
)$duration
write.csv(
  cbind(GEOID = rownames(traveltimes), as.data.frame(as.matrix(traveltimes))),
  paste0(maindir, "traveltimes.csv"), row.names = FALSE
)
system2("bzip2", shQuote(paste0(maindir, "traveltimes.csv")))

library(Matrix)
commutes <- sparseMatrix(
  {}, {}, x = 0,
  dims = rowSums(vapply(data, function(d) dim(d$commutes), numeric(2))),
  dimnames = rep(list(do.call(c, unname(lapply(data, function(d) colnames(d$commutes))))), 2)
)
for(d in data) commutes[rownames(d$commutes), colnames(d$commutes)] <- d$commutes
write.csv(
  cbind(GEOID = rownames(commutes), as.data.frame(as.matrix(unname(commutes)))),
  paste0(maindir, "commutes.csv"), row.names = FALSE
)
system2("bzip2", shQuote(paste0(maindir, "commutes.csv")))

#
# Floating Catchment Area calculations
#

# reload datasets if needed (can start from here)
if(!exists("provider_locations")) provider_locations <- read.csv(paste0(maindir, "providers.csv"))
if(!exists("data_combined")) data_combined <- read.csv(paste0(maindir, "data.csv"))
if(!exists("traveltimes")){
  con <- bzfile(paste0(maindir, "traveltimes.csv.bz2"))
  traveltimes <- read.csv(con, row.names = 1)
  traveltimes <- as(as.matrix(traveltimes), "dgCMatrix")
}
if(!exists("commutes")){
  con <- bzfile(paste0(maindir, "commutes.csv.bz2"))
  commutes <- read.csv(con, row.names = 1)
  colnames(commutes) <- rownames(commutes)
  commutes <- as(as.matrix(commutes), "dgCMatrix")
}

# baseline -- 3-step with step weights up to 60 minutes
weight <- list(c(60, .042), c(30, .377), c(20, .704), c(10, .962))

data_combined$access_3sfca <- catchment_ratio(
  data_combined, provider_locations, traveltimes, weight, normalize_weight = TRUE,
  consumers_value = "population", providers_id = "id", providers_value = "doctors",
  return_type = "region", verbose = TRUE
)

# 2-step version
data_combined$access_2sfca <- catchment_ratio(
  data_combined, provider_locations, traveltimes, weight,
  consumers_value = "population", providers_id = "id", providers_value = "doctors",
  return_type = "region", verbose = TRUE
)

# shorter range version
data_combined$access_3sfca_30 <- catchment_ratio(
  data_combined, provider_locations, traveltimes, weight[-1], normalize_weight = TRUE,
  consumers_value = "population", providers_id = "id", providers_value = "doctors",
  return_type = "region", verbose = TRUE
)

# bounded continuous version
data_combined$access_3sfca_gauss <- catchment_ratio(
  data_combined, provider_locations, traveltimes, "gaussian", max_cost = 60, normalize_weight = TRUE,
  consumers_value = "population", providers_id = "id", providers_value = "doctors",
  return_type = "region", verbose = TRUE
)

# a commuter-adjusted version
data_combined$access_3sfca_commute <- catchment_ratio(
  data_combined, provider_locations, traveltimes, weight, normalize_weight = TRUE,
  consumers_value = "population", providers_id = "id", providers_value = "doctors",
  consumers_commutes = commutes, return_type = "region", verbose = TRUE
)

# make final datasets at each geography level
vars = colnames(data_combined)[-c(1, 6:7)]
write.csv(data_combined[, -(6:7)], paste0(maindir, "blockgroups.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(split(data_combined, substr(data_combined$GEOID, 1, 11)), function(d){
  if(is.null(dim(d))) d <- as.data.frame(as.list(d))
  data.frame(
    GEOID = substr(d[1, "GEOID"], 1, 11),
    as.list(colSums(d[, vars], na.rm = TRUE) / c(1, rep(c(nrow(d), 1), c(3, 5)))
    ))
})), paste0(maindir, "tracts.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(split(data_combined, substr(data_combined$GEOID, 1, 5)), function(d) data.frame(
  GEOID = substr(d[1, "GEOID"], 1, 5),
  as.list(colSums(d[, vars], na.rm = TRUE) / c(1, rep(c(nrow(d), 1), c(3, 5)))
  )))), paste0(maindir, "counties.csv"), row.names = FALSE)

#
# add data to site, along with metadata
#

data_add(
  c(
    counties = "counties.csv",
    tracts = "tracts.csv",
    blockgroups = "blockgroups.csv"
  ),
  rep(list(list(
    ids = list(variable = "GEOID"),
    variables = list(
      population = list(
        long_name = "Total Population",
        description = "Total population as estimated by the U.S. Census Bureau.",
        statement = "There are {value} people in {features.name}.",
        type = "count",
        source = list(
          list(
            name = "American Community Survey",
            date_accessed = 2021,
            url = "https://www.census.gov/programs-surveys/acs.html"
          )
        )
      ),
      percent_female = list(
        long_name = "Female Percent of the Population",
        description = "Percent of the population identified as female.",
        statement = "{value} of the population in {features.name} identified as Female.",
        type = "percent",
        source = list(
          list(
            name = "American Community Survey",
            date_accessed = 2021,
            url = "https://www.census.gov/programs-surveys/acs.html"
          )
        )
      ),
      percent_White = list(
        long_name = "White Percent of the Population",
        description = "Percent of the population identified as White.",
        statement = "{value} of the population in {features.name} identified as White",
        type = "percent",
        source = list(
          list(
            name = "American Community Survey",
            date_accessed = 2021,
            url = "https://www.census.gov/programs-surveys/acs.html"
          )
        )
      ),
      percent_over_49 = list(
        long_name = "Older Percent of the Population",
        description = "Percent of the population who are 50 years old or older.",
        statement = "{value} of the population in {features.name} reported being over 49 years old.",
        type = "percent",
        source = list(
          list(
            name = "American Community Survey",
            date_accessed = 2021,
            url = "https://www.census.gov/programs-surveys/acs.html"
          )
        )
      ),
      access_3sfca = list(
        long_name = "Doctors (3-Step Floating Catchment Area)",
        short_name = "Doctors (3SFCA)",
        description = paste(
          "Number of doctors available in the region, as calculated within",
          "floating catchment areas (60-minute radius) with normalized weights."
        ),
        statement = "There are {value} doctors available in {features.name}.",
        type = "sum",
        citations = "wan12",
        source = list(
          list(
            name = "American Community Survey",
            date_accessed = 2021,
            url = "https://www.census.gov/programs-surveys/acs.html"
          ),
          list(
            name = "Centers for Medicare & Medicaid Services",
            date_accessed = 2021,
            url = "https://data.cms.gov"
          )
        )
      ),
      access_2sfca = list(
        long_name = "Doctors (2-Step Floating Catchment Area)",
        short_name = "Doctors (2SFCA)",
        description = paste(
          "Number of doctors available in the region, as calculated within",
          "floating catchment areas (60-minute radius) without normalized weights."
        ),
        statement = "There are {value} doctors available in {features.name}.",
        type = "sum",
        citations = "wan12",
        source = list(
          list(
            name = "American Community Survey",
            date_accessed = 2021,
            url = "https://www.census.gov/programs-surveys/acs.html"
          ),
          list(
            name = "Centers for Medicare & Medicaid Services",
            date_accessed = 2021,
            url = "https://data.cms.gov"
          )
        )
      ),
      access_3sfca_30 = list(
        long_name = "Doctors (3-Step Floating Catchment Area Ratio, 30 minutes)",
        short_name = "Doctors (3SFCA30)",
        description = paste(
          "Number of doctors available in the region, as calculated within",
          "smaller floating catchment areas (30-minute radius) with normalized weights."
        ),
        statement = "There are {value} doctors available in {features.name}.",
        type = "sum",
        citations = "wan12",
        source = list(
          list(
            name = "American Community Survey",
            date_accessed = 2021,
            url = "https://www.census.gov/programs-surveys/acs.html"
          ),
          list(
            name = "Centers for Medicare & Medicaid Services",
            date_accessed = 2021,
            url = "https://data.cms.gov"
          )
        )
      ),
      access_3sfca_gauss = list(
        long_name = "Doctors (3-Step Floating Catchment Area Ratio, Gaussian)",
        short_name = "Doctors (KD3SFCA)",
        description = paste(
          "Number of doctors available in the region, as calculated within",
          "floating catchment areas (60-minute radius) with normalized, continuous weights."
        ),
        statement = "There are {value} doctors available in {features.name}.",
        type = "sum",
        citations = "wan12",
        source = list(
          list(
            name = "American Community Survey",
            date_accessed = 2021,
            url = "https://www.census.gov/programs-surveys/acs.html"
          ),
          list(
            name = "Centers for Medicare & Medicaid Services",
            date_accessed = 2021,
            url = "https://data.cms.gov"
          )
        )
      ),
      access_3sfca_commute = list(
        long_name = "Doctors (3-Step Floating Catchment Area Ratio, commuter-based)",
        short_name = "Doctors (CB3SFCA)",
        description = paste(
          "Number of doctors available in the region, as calculated within",
          "floating catchment areas (60-minute radius) with normalized weights",
          "adjusted for multiple origins (home and work)."
        ),
        statement = "There are {value} doctors available in {features.name}.",
        type = "sum",
        citations = "wan12",
        source = list(
          list(
            name = "American Community Survey",
            date_accessed = 2021,
            url = "https://www.census.gov/programs-surveys/acs.html"
          ),
          list(
            name = "Centers for Medicare & Medicaid Services",
            date_accessed = 2021,
            url = "https://data.cms.gov"
          ),
          list(
            name = "Longitudinal Employer-Household Dynamics",
            date_accessed = 2019,
            url = "https://lehd.ces.census.gov"
          )
        )
      ),
      "_references" = list(
        lou09 = list(
          author = list(
            list(given = "Wei", family = "Lou"),
            list(given = "Yi", family = "Qi")
          ),
          year = 2009,
          title = paste(
            "An enhanced two-step floating catchment area (e2sfca) method for measuring spatial",
            "accessibility to primary care physicians"
          ),
          journal = "Health & Place",
          volume = 15,
          page = "1100-1107",
          doi = "10.1016/j.healthplace.2009.06.002"
        ),
        wan12 = list(
          author = list(
            list(given = "Neng", family = "Wan"),
            list(given = "Bin", family = "Zou"),
            list(given = "Troy", family = "Sternberg")
          ),
          year = 2012,
          title = "A three-step floating catchment area method for analyzing spatial access to health services",
          journal = "International Journal of Geographical Information Science",
          volume = 26,
          page = "1073-1089",
          doi = "10.1080/13658816.2011.624987"
        )
      )
    )
  )), 3),
  dir = maindir,
  refresh = TRUE
)
