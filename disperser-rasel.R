#rm(list=ls())

#installing packages
install.packages(c("Rcpp", "downloader", "rgdal"))


#installing SplitR package
devtools::install_github("lhenneman/SplitR", force = TRUE)

# devtools::install_github("lhenneman/hyspdisp", force = TRUE)


install.packages("USAboundariesData", repos = "http://packages.ropensci.org", type = "source")

####### problem installing package (vignette building failed) ######


# devtools::install_github("munshimdrasel/disperseR", force = TRUE, build_vignettes = TRUE)

# package installation without vignette

devtools::install_github("munshimdrasel/disperseR", force = TRUE, build_vignettes = FALSE)

#setting R directory
setwd("/Users/munshirasel/hello-world/munshimdrasel/mello/mello2")

#loading libraries
library(disperseR) # our package
library(ncdf4)
library(data.table)
library(tidyverse)
library(parallel)
library(sf)
library(viridis)
library(ggplot2)
library(scales)
library(ggsn)
library(gridExtra)
library(ggmap)
library(ggrepel)
library(fst)

#creating directory for disperseR
disperseR::create_dirs(location = "/Users/munshirasel/hello-world/munshimdrasel/mello/mello2")

# download data
disperseR::get_data(data = "all",
                    start.year = "2005",
                    start.month = "01",
                    end.year = "2006",
                    end.month = "05")


# view units data
view(disperseR::units-rasel)
# pick out units to run, top two SOx emitters in 2006 & 2005
unitsrun2005 <- disperseR::units %>%
  dplyr::filter(year == 2005) %>% # only get data for 2005
  dplyr::top_n(2, SOx)  # sort and take the two rows with the biggest value for SOx
unitsrun2006 <- disperseR::units %>%
  dplyr::filter(year == 2006) %>%  # only get data for 2006
  dplyr::top_n(2, SOx)  # sort and take the two rows with the biggest value for SOx

# append together and transform to data table
unitsrun<-data.table::data.table(rbind(unitsrun2005, unitsrun2006))

# find unique combos of Latitude, Longitude, and Height
unitslatlonh <- unique( unitsrun[ ,.( Latitude, Longitude, Height, year)] )
unitslatlonh[, unit_run_ref:=1:nrow( unitslatlonh)]
unitsrun_trim <- merge( unitsrun, unitslatlonh)[ !duplicated( unit_run_ref)]

# define data.table with all emissions events
input_refs <- disperseR::define_inputs(units = unitsrun,
                                       startday = '2005-11-01',
                                       endday = '2006-02-28',
                                       start.hours =  c(0, 6, 12, 18),
                                       duration = 120)

head(input_refs, 10)
# subset the input refs
input_refs_subset <- input_refs[format(as.Date(input_refs$start_day,
                                               format = "%Y-%m-%d"),
                                       format = "%d") == "01" & start_hour == 0]

head (input_refs_subset, 10)
# run disperser
hysp_raw <- disperseR::run_disperser_parallel(input.refs = input_refs_subset,
                                              pbl.height = pblheight,
                                              species = 'so2',
                                              proc_dir = proc_dir,
                                              overwrite = FALSE, ## FALSE BY DEFAULT
                                              npart = 100,
                                              keep.hysplit.files = FALSE, ## FALSE BY DEFAULT
                                              mc.cores = parallel::detectCores())

# Link results to spatial domains


yearmons <- disperseR::get_yearmon(start.year = "2005",
                                   start.month = "07",
                                   end.year = "2006",
                                   end.month = "06")

unitsrun

linked_zips <- disperseR::link_all_units(
  units.run = unitsrun,
  link.to = 'zips',
  mc.cores = parallel::detectCores(),
  year.mons = yearmons,
  pbl.height = pblheight,
  crosswalk. = crosswalk,
  duration.run.hours = 240,
  res.link = 12000,
  overwrite = FALSE)
#> processed unit 3136-1
#> processed unit 3149-1
#> processed unit 3136-2

# link all units to counties
linked_counties <- disperseR::link_all_units(
  units.run=unitsrun,
  link.to = 'counties',
  mc.cores = parallel::detectCores(),
  year.mons = yearmons,
  pbl.height = pblheight,
  counties. = USAboundaries::us_counties( ),
  crosswalk. = NULL,
  duration.run.hours = 240,
  overwrite = FALSE)
#> processed unit 3136-1
#> processed unit 3149-1
#> processed unit 3136-2

####problem######
# link all units to grids
linked_grids <- disperseR::link_all_units(
  units.run=unitsrun,
  link.to = 'grids',
  mc.cores = parallel::detectCores(),
  year.mons = yearmons,
  pbl.height = pblheight,
  crosswalk. = NULL,
  duration.run.hours = 240,
  overwrite = FALSE)
#> processed unit 3136-1
#> processed unit 3149-1
#> processed unit 3136-2


head(linked_zips)
head(linked_counties)
head(linked_grids)


unique(linked_zips$comb)


# Visualization of the results.

impact_table_zip_single <- disperseR::create_impact_table_single(
  data.linked=linked_zips,
  link.to = 'zips',
  data.units = unitsrun,
  zcta.dataset = zcta_dataset,
  map.unitID = "3136-1",
  map.month = "200511",
  metric = 'N')

impact_table_county_single <- disperseR::create_impact_table_single(
  data.linked=linked_counties,
  link.to = 'counties',
  data.units = unitsrun,
  counties. = USAboundaries::us_counties( ),
  map.unitID = "3136-1",
  map.month = "200511",
  metric = 'N')
impact_table_grid_single <- disperseR::create_impact_table_single(
  data.linked=linked_grids,
  link.to = 'grids',
  data.units = unitsrun,
  map.unitID = "3136-1",
  map.month = "200511",
  metric = 'N')

head(impact_table_zip_single)

link_plot_zips <- disperseR::plot_impact_single(
  data.linked = linked_zips,
  link.to = 'zips',
  map.unitID = "3136-1",
  map.month = "20061",
  data.units = unitsrun,
  zcta.dataset = zcta_dataset,
  metric = 'N',
  graph.dir = graph_dir,
  zoom = T, # TRUE by default
  legend.name = 'HyADS raw exposure',
  # other parameters passed to ggplot2::theme()
  axis.text = element_blank(),
  legend.position = c( .75, .15))

link_plot_grids <- disperseR::plot_impact_single(
  data.linked = linked_grids,
  link.to = 'grids',
  map.unitID = "3136-1",
  map.month = "20061",
  data.units = unitsrun,
  metric = 'N',
  graph.dir = graph_dir,
  zoom = F, # TRUE by default (false meaning to show the whole country)
  legend.name = 'HyADS raw exposure',
  # other parameters passed to ggplot2::theme()
  axis.text = element_blank(),
  legend.position = c( .75, .15))

link_plot_counties <- disperseR::plot_impact_single(
  data.linked = linked_counties,
  link.to = 'counties',
  map.unitID = "3136-1",
  map.month = "20061",
  counties. = USAboundaries::us_counties( ),
  data.units = unitsrun,
  metric = 'N',
  graph.dir = graph_dir,
  zoom = T, # TRUE by default (true means to show that location area only)
  legend.name = 'HyADS raw exposure',
  # other parameters passed to ggplot2::theme()
  axis.text = element_blank(),
  legend.position = c( .75, .15))

link_plot_zips
link_plot_grids
link_plot_counties

# Combine all results into RData file.
combined_ziplinks <- disperseR::combine_monthly_links(
  month_YYYYMMs = yearmons,
  link.to = 'zips',
  filename = 'hyads_vig_unwgted_zips.RData')

combined_countylinks <- disperseR::combine_monthly_links(
  month_YYYYMMs = yearmons,
  link.to = 'counties',
  filename = 'hyads_vig_unwgted_counties.RData')

combined_gridlinks <- disperseR::combine_monthly_links(
  month_YYYYMMs = yearmons,
  link.to = 'grids',
  filename = 'hyads_vig_unwgted_grids.RData')

names(combined_ziplinks)


# Calculate and extract useful information from the results

# Weight the results by emissions


exp_ann_unit_zip <- disperseR::calculate_exposure(
  year.E = 2005,
  year.D = 2005,
  link.to = 'zips',
  pollutant = 'SO2.tons',
  rda_file = file.path(rdata_dir, "hyads_vig_unwgted_zips.RData"),
  exp_dir = exp_dir,
  units.mo = PP.units.monthly1995_2017,
  source.agg = 'unit',
  time.agg = 'month',
  return.monthly.data = T)

exp_ann_unit_grids <- disperseR::calculate_exposure(
  year.E = 2005,
  year.D = 2005,
  link.to = 'grids',
  pollutant = 'SO2.tons',
  rda_file = file.path(rdata_dir, "hyads_vig_unwgted_grids.RData"),
  exp_dir = exp_dir,
  units.mo = PP.units.monthly1995_2017,
  source.agg = 'unit',
  time.agg = 'month',
  return.monthly.data = T)

exp_ann_unit_counties <- disperseR::calculate_exposure(
  year.E = 2005,
  year.D = 2005,
  link.to = 'counties',
  pollutant = 'SO2.tons',
  rda_file = file.path(rdata_dir, "hyads_vig_unwgted_counties.RData"),
  exp_dir = exp_dir,
  units.mo = PP.units.monthly1995_2017,
  source.agg = 'unit',
  time.agg = 'month',
  return.monthly.data = T)

zip_exp_ann_plot <- disperseR::plot_impact_weighted(
  data.linked = exp_ann_unit_zip,
  data.units = unitsrun,
  link.to = 'zips',
  zcta.dataset = zcta_dataset,
  time.agg = 'year',
  metric = 'hyads',
  legend.name = 'Aggregate HyADS exposure',
  zoom = T, # TRUE by default
  graph.dir = graph_dir,
  map.month = NULL, # NULL by default change if time.agg = 'month'
  # other parameters passed to ggplot2::theme()
  axis.text = element_blank(),
  legend.position = c( .75, .15)) # 0 by default

counties_exp_ann_plot <- disperseR::plot_impact_weighted(
  data.linked = exp_ann_unit_counties,
  data.units = unitsrun,
  link.to = 'counties',
  counties. = USAboundaries::us_counties( ),
  time.agg = 'year',
  metric = 'hyads',
  legend.name = 'Aggregate HyADS exposure',
  zoom = T, # TRUE by default
  graph.dir = graph_dir,
  map.month = NULL, # NULL by default change if time.agg = 'month'
  # other parameters passed to ggplot2::theme()
  axis.text = element_blank(),
  legend.position = c( .75, .15)) # 0 by default

grids_exp_ann_plot <- disperseR::plot_impact_weighted(
  data.linked = exp_ann_unit_grids,
  data.units = unitsrun,
  link.to = 'grids',
  time.agg = 'year',
  metric = 'hyads',
  legend.name = 'Aggregate HyADS exposure',
  zoom = T, # TRUE by default
  graph.dir = graph_dir,
  map.month = NULL, # NULL by default change if time.agg = 'month'
  # other parameters passed to ggplot2::theme()
  axis.text = element_blank(),
  legend.position = c( .75, .15)) # 0 by default

zip_exp_ann_plot

counties_exp_ann_plot

grids_exp_ann_plot

#plotting monthly exposure


exp_mon_unit_zip <- disperseR::calculate_exposure(
  year.E = 2005,
  year.D = 2005,
  link.to = 'zips',
  pollutant = 'SO2.tons',
  rda_file = file.path(rdata_dir, "hyads_vig_unwgted_zips.RData"),
  exp_dir = exp_dir,
  units.mo = PP.units.monthly1995_2017,
  source.agg = 'unit',
  time.agg = 'month',
  return.monthly.data = T)

exp_mon_unit_grids <- disperseR::calculate_exposure(
  year.E = 2005,
  year.D = 2005,
  link.to = 'grids',
  pollutant = 'SO2.tons',
  rda_file = file.path(rdata_dir, "hyads_vig_unwgted_grids.RData"),
  exp_dir = exp_dir,
  units.mo = PP.units.monthly1995_2017,
  source.agg = 'unit',
  time.agg = 'month',
  return.monthly.data = T)

exp_mon_unit_counties <- disperseR::calculate_exposure(
  year.E = 2005,
  year.D = 2005,
  link.to = 'counties',
  pollutant = 'SO2.tons',
  rda_file = file.path(rdata_dir, "hyads_vig_unwgted_counties.RData"),
  exp_dir = exp_dir,
  units.mo = PP.units.monthly1995_2017,
  source.agg = 'unit',
  time.agg = 'month',
  return.monthly.data = T)

zip_exp_mon_plot <- disperseR::plot_impact_weighted(
  data.linked = exp_mon_unit_zip,
  data.units = unitsrun,
  zcta.dataset = zcta_dataset,
  time.agg = 'month',
  map.month = "200511",
  metric = 'hyads',
  legend.name = 'Montly HyADS exposure',
  zoom = T, # TRUE by default
  graph.dir = graph_dir)

zip_exp_mon_plot

# Plot unit-specific impacts over time

zip_exp_unit_mon2005 <- disperseR::calculate_exposure(
  rda_file = file.path(rdata_dir, "hyads_vig_unwgted_zips.RData"),
  units.mo = PP.units.monthly1995_2017,
  link.to = 'zips',
  year.E = 2005,
  year.D = 2005,
  pollutant = 'SO2.tons',
  source.agg = 'unit', # note!
  time.agg = 'month',
  exp_dir = exp_dir,
  return.monthly.data = T)

zip_exp_unit_mon2006 <- disperseR::calculate_exposure(
  rda_file = file.path(rdata_dir, "hyads_vig_unwgted_zips.RData"),
  units.mo = PP.units.monthly1995_2017,
  link.to = 'zips',
  year.E = 2006,
  year.D = 2006,
  pollutant = 'SO2.tons',
  source.agg = 'unit', # note!
  time.agg = 'month',
  exp_dir = exp_dir,
  return.monthly.data = T)

zip_exp_unit_mon <- rbind(zip_exp_unit_mon2005, zip_exp_unit_mon2006)

zipcodes <- c("13039","21798", "03804")


###????
zip_exp_unit <- disperseR::plot_impact_unit(
  data.linked = zip_exp_unit_mon,
  zip.codes = zipcodes,
  graph.dir = graph_dir)
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?


# Rank facilities.

zip_exp_ann_unit <- disperseR::calculate_exposure(
  year.E = 2005,
  year.D = 2005,
  link.to = 'zips',
  pollutant = 'SO2.tons',
  rda_file = file.path(rdata_dir, "hyads_vig_unwgted_zips.RData"),
  exp_dir = exp_dir,
  units.mo = PP.units.monthly1995_2017,
  source.agg = 'unit',
  time.agg = 'year')

zip_exp_ann_unit[, year := 2005]

unitRanks2005 <- disperseR::rankfacs_by_popwgt_location(
  data.linked = zip_exp_ann_unit,
  crosswalk. = crosswalk,
  rank.by = c('hyads'),
  state.value = 'PA',
  year = 2005)

unitRanks2005


# Plot ranked facilities.

plotUnitsRanked <- disperseR::plot_units_ranked(
  data.units = unitsrun,
  data.ranked = unitRanks2005,
  year = 2005,
  graph.dir = graph_dir)

plotUnitsRanked
#> $ggbar

#>
#> $ggmap


##End of code

# not uploading main folders due to large size

