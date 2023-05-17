# specify the main directory
maindir <- "../dmv_healthcare/docs/data/"

# make a directory to store working data files in
oridir <- paste0(maindir, "original/")
dir.create(oridir, FALSE, TRUE)

#
# population information
#

library(community)
library(catchment)

# define counties that are part of the capital region
dmv_counties <- list(
  dc = "District of Columbia",
  md = c("Charles", "Frederick", "Montgomery", "Prince George's"),
  va = c("Alexandria", "Arlington", "Fairfax", "Falls Church", "Loudoun", "Manassas", "Manassas Park", "Prince William")
)
data <- list()
shapes <- list()

# download / load
for (state in c("dc", "md", "va")) {
  # shapes
  counties <- download_census_shapes(oridir, state, "county", paste0(state, "_counties"))
  tracts <- download_census_shapes(oridir, state, "tract", paste0(state, "_tracts"))
  blockgroups <- download_census_shapes(oridir, state, "bg", paste0(state, "_blockgroups"))

  ## store subsets to combine later
  counties <- counties[counties$NAME %in% dmv_counties[[state]], ]
  counties[counties$NAME == "Fairfax", "NAME"] <- c("Fairfax City", "Fairfax")
  shapes[[state]] <- list(
    counties = counties,
    tracts = tracts[substr(tracts$GEOID, 1, 5) %in% counties$GEOID, ],
    blockgroups = blockgroups[substr(blockgroups$GEOID, 1, 5) %in% counties$GEOID, ]
  )

  # population data
  data[[state]] <- download_census_population(
    oridir, state, 2019,
    include_margins = TRUE, include_commutes = TRUE,
    counties = counties$GEOID, verbose = TRUE
  )
}

## create and save combined shapes
library(sf)
library(rmapshaper)

for (level in names(shapes$dc)) {
  st_write(
    ms_simplify(do.call(rbind, lapply(shapes, "[[", level)), keep_shapes = TRUE),
    paste0(maindir, paste0(level, ".geojson"))
  )
}

## create and save final files
data_combined <- do.call(rbind, lapply(names(data), function(state) {
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
    st_coordinates(st_centroid(st_geometry(s[as.character(d$GEOID), ])))
  )
}))
write.csv(data_combined, paste0(maindir, "data.csv"), row.names = FALSE)

library(Matrix)
commutes <- sparseMatrix(
  {},
  {},
  x = 0,
  dims = rowSums(vapply(data, function(d) dim(d$commutes), numeric(2))),
  dimnames = rep(list(do.call(c, unname(lapply(data, function(d) colnames(d$commutes))))), 2)
)
for (d in data) commutes[rownames(d$commutes), colnames(d$commutes)] <- d$commutes
write.csv(
  cbind(GEOID = rownames(commutes), as.data.frame(as.matrix(unname(commutes)))),
  paste0(maindir, "commutes.csv"),
  row.names = FALSE
)
system2("bzip2", shQuote(paste0(maindir, "commutes.csv")))


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
while (length(raw <- jsonlite::read_json(paste0(providers_url, offset)))) {
  cat("\rdownloading file", offset / 1000 + 1, "     ")
  providers[[length(providers) + 1]] <- do.call(rbind, lapply(raw, unlist))
  offset <- offset + 1000
}
providers <- as.data.frame(do.call(rbind, providers))

## get a set of ZIP codes within the focal counties
county_shapes <- read_sf(paste0(maindir, "counties.geojson"), as_tibble = FALSE)
geography_ref <- read.csv("https://www2.census.gov/geo/docs/maps-data/data/rel/zcta_county_rel_10.txt")
zips <- unique(unlist(lapply(names(dmv_counties), function(state) {
  GEOIDs <- county_shapes[county_shapes$NAME %in% dmv_counties[[state]], "GEOID", drop = TRUE]
  formatC(geography_ref[geography_ref$GEOID %in% GEOIDs, "ZCTA5"], width = 5, flag = 0)
}), use.names = FALSE))

## focus only on medical doctors who aren't also dentists, within the selected counties
providers <- providers[
  grepl("\\bm\\.?d(?:\\W|$)", providers$Rndrng_Prvdr_Crdntls, TRUE) &
    !grepl("d\\.?d\\.?s", providers$Rndrng_Prvdr_Crdntls, TRUE) &
    providers$Rndrng_Prvdr_Zip5 %in% zips,
]

## format address strings, and get the coordinates of each unique entry
address_parts <- c(
  "Rndrng_Prvdr_St1", "Rndrng_Prvdr_City", "Rndrng_Prvdr_State_Abrvtn", "Rndrng_Prvdr_Zip5"
)

providers$Rndrng_Prvdr_St1 <- sub("([nesw])\\.([nesw])\\.*", "\\1\\2", providers$Rndrng_Prvdr_St1, TRUE)
providers$Rndrng_Prvdr_St1 <- sub("\\s*#\\s.*$", "", providers$Rndrng_Prvdr_St1)
providers$Rndrng_Prvdr_St1 <- sub("\\s*([nesw]{2})\\s.*$", " \\1", providers$Rndrng_Prvdr_St1, TRUE)
providers$address <- do.call(paste, c(unname(as.list(providers[, address_parts])), sep = ", "))

## collapse to locations based on address
vars <- c(
  "address", "X", "Y", "Rndrng_Prvdr_Gndr",
  grep("^(tot|drug|med|bene)_", colnames(providers), TRUE, value = TRUE)
)
vars <- vars[!vars %in% c("Drug_Sprsn_Ind", "Med_Sprsn_Ind")]
addresses <- unique(providers$address)

## geocode addresses; takes a while
library(parallel)

cl <- makeCluster(detectCores() - 2)
address_coords <- as.data.frame(do.call(rbind, parLapply(cl, addresses, function(a) {
  coords <- tidygeocoder::geo(a, progress_bar = FALSE, quiet = TRUE, method = "arcgis")
  if (is.na(coords$long)) coords <- tidygeocoder::geo(a, progress_bar = FALSE, quiet = TRUE)
  coords
})))
rownames(address_coords) <- address_coords$address
stopCluster(cl)

## add coordinates to providers data
providers[, c("Y", "X")] <- address_coords[providers$address, c("lat", "long")]
providers <- providers[!is.na(providers$X), ]
providers$locid <- paste0(providers$X, ",", providers$Y)

provider_locations <- do.call(rbind, lapply(unique(providers$locid), function(l) {
  d <- providers[providers$locid == l, vars]
  d[d == ""] <- NA
  as.data.frame(list(
    address = d[1, "address"],
    X = d[1, "X"],
    Y = d[1, "Y"],
    doctors = nrow(d),
    prop_women = mean(d$Rndrng_Prvdr_Gndr == "F"),
    as.list(colMeans(matrix(
      as.numeric(as.matrix(d[, -(1:4)])), nrow(d),
      dimnames = list(NULL, vars[-(1:4)])
    ), na.rm = TRUE))
  ))
}))
provider_locations[is.na(provider_locations)] <- NA

## identify zip codes that cross counties
zip_cross <- substr(unique(do.call(
  paste0,
  geography_ref[geography_ref$ZCTA5 %in% zips, c("ZCTA5", "GEOID")]
)), 1, 5)
zip_cross <- zip_cross[duplicated(zip_cross)]

### check that those are actually in the focal counties
potential_ex <- provider_locations[
  grepl(paste0("(?:", paste(zip_cross, collapse = "|"), ")$"), provider_locations$address) &
    !grepl(paste0("(?:", paste(unlist(dmv_counties), collapse = "|"), "),"), provider_locations$address),
]
potential_ex$county <- tidygeocoder::reverse_geo(potential_ex$Y, potential_ex$X, full_results = TRUE)$county
provider_locations <- provider_locations[
  !provider_locations$address %in% potential_ex[
    !is.na(potential_ex$county) &
      !grepl(paste0("(?:", paste(unlist(dmv_counties), collapse = "|"), ")"), potential_ex$county),
    "address"
  ],
]

## make unique IDs for each provider location
provider_locations$ID <- paste0("l", seq_len(nrow(provider_locations)))

## save provider locations dataset
write.csv(provider_locations, paste0(maindir, "providers.csv"), row.names = FALSE)

## get travel times between each included block group and provider location
library(osrm)
options(osrm.server = Sys.getenv("OSRM_SERVER"))
traveltimes <- osrmTable(
  src = data_combined[, c("GEOID", "X", "Y")],
  dst = provider_locations[, c("ID", "X", "Y")]
)$duration
write.csv(
  cbind(GEOID = rownames(traveltimes), as.data.frame(as.matrix(traveltimes))),
  paste0(maindir, "traveltimes.csv"),
  row.names = FALSE
)
system2("bzip2", shQuote(paste0(maindir, "traveltimes.csv")))

#
# Floating Catchment Area calculations
#

library(Matrix)
library(catchment)

# reload datasets if needed (can start from here)
if (!exists("maindir")) maindir <- "../dmv_healthcare/docs/data/"
if (!exists("provider_locations")) provider_locations <- read.csv(paste0(maindir, "providers.csv"))
if (!exists("data_combined")) data_combined <- read.csv(paste0(maindir, "data.csv"))
if (!exists("traveltimes")) {
  con <- bzfile(paste0(maindir, "traveltimes.csv.bz2"))
  traveltimes <- read.csv(con, row.names = 1)
  traveltimes <- as(as.matrix(traveltimes), "dgCMatrix")
}
if (!exists("commutes")) {
  con <- bzfile(paste0(maindir, "commutes.csv.bz2"))
  commutes <- read.csv(con, row.names = 1)
  colnames(commutes) <- rownames(commutes)
  commutes <- as(as.matrix(commutes), "dgCMatrix")
}

# baseline -- 3-step with Gaussian weights
data_combined$doctors_3sfca <- catchment_ratio(
  data_combined, provider_locations, traveltimes, "gaussian",
  normalize_weight = TRUE,
  consumers_value = "population", providers_id = "ID", providers_value = "doctors",
  scale = 18, return_type = 1000
)

# 2-step version
data_combined$doctors_kd2sfca <- catchment_ratio(
  data_combined, provider_locations, traveltimes, "gaussian",
  consumers_value = "population", providers_id = "ID", providers_value = "doctors",
  scale = 18, return_type = 1000
)

## versions with euclidean distances
data_combined$doctors_3sfca_euclidean <- catchment_ratio(
  data_combined, provider_locations,
  weight = "gaussian", normalize_weight = TRUE,
  consumers_value = "population", providers_id = "ID", providers_value = "doctors",
  scale = 18, return_type = 1000
)
data_combined$doctors_kd2sfca_euclidean <- catchment_ratio(
  data_combined, provider_locations,
  weight = "gaussian",
  consumers_value = "population", providers_id = "ID", providers_value = "doctors",
  scale = 18, return_type = 1000
)

# shorter range version
data_combined$doctors_3sfca_30 <- catchment_ratio(
  data_combined, provider_locations, traveltimes, "gaussian",
  max_cost = 30, normalize_weight = TRUE,
  consumers_value = "population", providers_id = "ID", providers_value = "doctors",
  scale = 18, return_type = 1000
)

# step version
weight <- list(c(60, .042), c(30, .377), c(20, .704), c(10, .962))
data_combined$doctors_3sfca_step <- catchment_ratio(
  data_combined, provider_locations, traveltimes, weight,
  normalize_weight = TRUE, return_type = 1000,
  consumers_value = "population", providers_id = "ID", providers_value = "doctors"
)

# step 2-step version
data_combined$doctors_e2sfca <- catchment_ratio(
  data_combined, provider_locations, traveltimes, weight,
  return_type = 1000,
  consumers_value = "population", providers_id = "ID", providers_value = "doctors"
)

# a commuter-adjusted version
data_combined$doctors_3sfca_commute <- catchment_ratio(
  data_combined, provider_locations, traveltimes, "gaussian",
  normalize_weight = TRUE, return_type = 1000,
  consumers_value = "population", providers_id = "ID", providers_value = "doctors",
  scale = 18, consumers_commutes = commutes
)

doctors_vars <- grep("doctors_", colnames(data_combined), fixed = TRUE, value = TRUE)

# save complete dataset
write.csv(data_combined[, -(6:7)], paste0(maindir, "blockgroups.csv"), row.names = FALSE)

# make county-level aggregates

## make a new set with just population variables initially
counties <- do.call(rbind, lapply(split(data_combined, substr(data_combined$GEOID, 1, 5)), function(d) {
  data.frame(
    GEOID = substr(d[1, "GEOID"], 1, 5),
    population = sum(d$population, na.rm = TRUE),
    as.list(colMeans(d[, grepl("percent_", colnames(d), fixed = TRUE)], na.rm = TRUE))
  )
}))

## add aggregated catchment ratios, matched by GEOID substring
for (variable in doctors_vars) {
  counties[, variable] <- catchment_aggregate(data_combined, 5, value = variable)
}

## save aggregated dataset
write.csv(counties, paste0(maindir, "counties.csv"), row.names = FALSE)

# make tract-level aggregates
tracts <- do.call(rbind, lapply(split(data_combined, substr(data_combined$GEOID, 1, 11)), function(d) {
  if (is.null(dim(d))) d <- as.data.frame(as.list(d))
  data.frame(
    GEOID = substr(d[1, "GEOID"], 1, 11),
    population = sum(d$population, na.rm = TRUE),
    as.list(colMeans(d[, grepl("percent_", colnames(d), fixed = TRUE)], na.rm = TRUE))
  )
}))
for (variable in doctors_vars) {
  tracts[, variable] <- catchment_aggregate(data_combined, 11, value = variable)
}
write.csv(tracts, paste0(maindir, "tracts.csv"), row.names = FALSE)

#
# add data to site, along with metadata
#

data_add(
  c(
    county = "counties.csv",
    tract = "tracts.csv",
    block_group = "blockgroups.csv"
  ),
  list(
    ids = list(
      variable = "GEOID",
      map = "https://raw.githubusercontent.com/uva-bi-sdad/capital_region/main/docs/data/entity_info.json"
    ),
    variables = list(
      population = list(
        long_name = "Total Population",
        description = "Total population as estimated by the U.S. Census Bureau.",
        statement = "There are {value} people in {features.name}.",
        type = "count",
        source = list(
          list(
            name = "American Community Survey",
            date_accessed = 2019,
            url = "https://www.census.gov/programs-surveys/acs.html"
          )
        )
      ),
      percent_female = list(
        long_name = "Female Percent of the Population",
        description = "Percent of the population identified as female.",
        statement = "{value} of the population in {features.name} identified as Female.",
        type = "percent",
        sources = list(
          list(
            name = "American Community Survey",
            date_accessed = 2019,
            url = "https://www.census.gov/programs-surveys/acs.html"
          )
        )
      ),
      percent_White = list(
        long_name = "White Percent of the Population",
        description = "Percent of the population identified as White.",
        statement = "{value} of the population in {features.name} identified as White",
        type = "percent",
        sources = list(
          list(
            name = "American Community Survey",
            date_accessed = 2019,
            url = "https://www.census.gov/programs-surveys/acs.html"
          )
        )
      ),
      percent_over_49 = list(
        long_name = "Older Percent of the Population",
        description = "Percent of the population who are 50 years old or older.",
        statement = "{value} of the population in {features.name} reported being over 49 years old.",
        type = "percent",
        sources = list(
          list(
            name = "American Community Survey",
            date_accessed = 2019,
            url = "https://www.census.gov/programs-surveys/acs.html"
          )
        )
      ),
      doctors_3sfca = list(
        long_name = "Doctors (3-Step Floating Catchment Area)",
        short_name = "Doctors (3SFCA)",
        description = paste(
          "Number of doctors available per 1,000 people, as calculated within",
          "floating catchment areas with normalized Gaussian weights."
        ),
        statement = "There are {value} doctors available per 1,000 people in {features.name}.",
        type = "sum",
        citations = "wan12",
        sources = list(
          list(
            name = "American Community Survey",
            date_accessed = 2019,
            url = "https://www.census.gov/programs-surveys/acs.html"
          ),
          list(
            name = "Centers for Medicare & Medicaid Services",
            date_accessed = 2019,
            url = "https://data.cms.gov"
          )
        )
      ),
      doctors_kd2sfca = list(
        long_name = "Doctors (Kernel Density 2-Step Floating Catchment Area)",
        short_name = "Doctors (KD2SFCA)",
        description = paste(
          "Number of doctors available per 1,000 people, as calculated within",
          "floating catchment areas with Gaussian weights."
        ),
        statement = "There are {value} doctors available per 1,000 people in {features.name}.",
        type = "sum",
        citations = "dai10",
        sources = list(
          list(
            name = "American Community Survey",
            date_accessed = 2019,
            url = "https://www.census.gov/programs-surveys/acs.html"
          ),
          list(
            name = "Centers for Medicare & Medicaid Services",
            date_accessed = 2019,
            url = "https://data.cms.gov"
          )
        )
      ),
      doctors_3sfca_euclidean = list(
        long_name = "Doctors (3-Step Floating Catchment Area, Euclidean)",
        short_name = "Doctors (3SFCA, Euclidean)",
        description = paste(
          "Number of doctors available per 1,000 people, as calculated within",
          "floating catchment areas with normalized Gaussian weights based on Euclidean distance."
        ),
        statement = "There are {value} doctors available per 1,000 people in {features.name}.",
        type = "sum",
        citations = "wan12",
        sources = list(
          list(
            name = "American Community Survey",
            date_accessed = 2019,
            url = "https://www.census.gov/programs-surveys/acs.html"
          ),
          list(
            name = "Centers for Medicare & Medicaid Services",
            date_accessed = 2019,
            url = "https://data.cms.gov"
          )
        )
      ),
      doctors_kd2sfca_euclidean = list(
        long_name = "Doctors (Kernel Density 2-Step Floating Catchment Area, Euclidean)",
        short_name = "Doctors (KD2SFCA, Euclidean)",
        description = paste(
          "Number of doctors available per 1,000 people, as calculated within",
          "floating catchment areas with Gaussian weights based on Euclidean distance."
        ),
        statement = "There are {value} doctors available per 1,000 people in {features.name}.",
        type = "sum",
        citations = "dai10",
        sources = list(
          list(
            name = "American Community Survey",
            date_accessed = 2019,
            url = "https://www.census.gov/programs-surveys/acs.html"
          ),
          list(
            name = "Centers for Medicare & Medicaid Services",
            date_accessed = 2019,
            url = "https://data.cms.gov"
          )
        )
      ),
      doctors_3sfca_30 = list(
        long_name = "Doctors (3-Step Floating Catchment Area Ratio, 30 minutes)",
        short_name = "Doctors (3SFCA30)",
        description = paste(
          "Number of doctors available per 1,000 people, as calculated within",
          "smaller floating catchment areas with normalized Gaussian weights within 30 minute buffers."
        ),
        statement = "There are {value} doctors available per 1,000 people in {features.name}.",
        type = "sum",
        citations = "wan12",
        sources = list(
          list(
            name = "American Community Survey",
            date_accessed = 2019,
            url = "https://www.census.gov/programs-surveys/acs.html"
          ),
          list(
            name = "Centers for Medicare & Medicaid Services",
            date_accessed = 2019,
            url = "https://data.cms.gov"
          )
        )
      ),
      doctors_3sfca_step = list(
        long_name = "Doctors (3-Step Floating Catchment Area Ratio)",
        short_name = "Doctors (3SFCA, Step)",
        description = paste(
          "Number of doctors available per 1,000 people, as calculated within",
          "floating catchment areas (60-minute radius) with normalized step weights."
        ),
        statement = "There are {value} doctors available per 1,000 people in {features.name}.",
        type = "sum",
        citations = "wan12",
        sources = list(
          list(
            name = "American Community Survey",
            date_accessed = 2019,
            url = "https://www.census.gov/programs-surveys/acs.html"
          ),
          list(
            name = "Centers for Medicare & Medicaid Services",
            date_accessed = 2019,
            url = "https://data.cms.gov"
          )
        )
      ),
      doctors_e2sfca = list(
        long_name = "Doctors (Enhanced 2-Step Floating Catchment Area Ratio)",
        short_name = "Doctors (E2SFCA)",
        description = paste(
          "Number of doctors available per 1,000 people, as calculated within",
          "floating catchment areas (60-minute radius) with step weights."
        ),
        statement = "There are {value} doctors available per 1,000 people in {features.name}.",
        type = "sum",
        citations = "luo09",
        sources = list(
          list(
            name = "American Community Survey",
            date_accessed = 2019,
            url = "https://www.census.gov/programs-surveys/acs.html"
          ),
          list(
            name = "Centers for Medicare & Medicaid Services",
            date_accessed = 2019,
            url = "https://data.cms.gov"
          )
        )
      ),
      doctors_3sfca_commute = list(
        long_name = "Doctors (3-Step Floating Catchment Area Ratio, commuter-based)",
        short_name = "Doctors (CB3SFCA)",
        description = paste(
          "Number of doctors available per 1,000 people, as calculated within",
          "floating catchment areas (60-minute radius) with normalized weights",
          "adjusted for multiple origins (home and work)."
        ),
        statement = "There are {value} doctors available per 1,000 people in {features.name}.",
        type = "sum",
        citations = c("wan12", "fransen15"),
        sources = list(
          list(
            name = "American Community Survey",
            date_accessed = 2019,
            url = "https://www.census.gov/programs-surveys/acs.html"
          ),
          list(
            name = "Centers for Medicare & Medicaid Services",
            date_accessed = 2019,
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
        dai10 = list(
          title = "Black residential segregation, disparities in spatial access to health care facilities, and late-stage breast cancer diagnosis in metropolitan Detroit",
          author = list(list(given = "Dajun", family = "Dai")),
          journal = "Health & place",
          volume = "16",
          pages = "1038--1052",
          year = "2010",
          publisher = "Elsevier",
          doi = "10.1016/j.healthplace.2010.06.012"
        ),
        luo09 = list(
          title = "An enhanced two-step floating catchment area (E2SFCA) method for measuring spatial accessibility to primary care physicians",
          author = list(list(given = "Wei", family = "Lou"), list(given = "Yi", family = "Qi")),
          journal = "Health & place",
          volume = "15",
          pages = "1100--1107",
          year = "2009",
          publisher = "Elsevier",
          doi = "10.1016/j.healthplace.2009.06.002"
        ),
        wan12 = list(
          title = "A three-step floating catchment area method for analyzing spatial access to health services",
          author = list(list(given = "Neng", family = "Wan"), list(given = "Bin", family = "Zou"), list(given = "Troy", family = "Sternberg")),
          journal = "International Journal of Geographical Information Science",
          volume = "26",
          pages = "1073--1089",
          year = "2012",
          publisher = "Taylor & Francis",
          doi = "10.1080/13658816.2011.624987"
        ),
        fransen15 = list(
          title = "A commuter-based two-step floating catchment area method for measuring spatial accessibility of daycare centers",
          author = list(
            list(given = "Koos", family = "Fransen"), list(given = "Tijs", family = "Neutens"),
            list(given = "Philippe", family = "De Maeyer"), list(given = "Greet", family = "Deruyter")
          ),
          journal = "Health & place",
          volume = "32",
          pages = "65--73",
          year = "2015",
          publisher = "Elsevier",
          doi = "10.1016/j.healthplace.2015.01.002"
        )
      )
    )
  ),
  dir = maindir
)
