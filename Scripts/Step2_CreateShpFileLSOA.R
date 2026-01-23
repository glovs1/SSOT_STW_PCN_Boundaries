library(sf)

-----------------------------------------------------------
  # 1) Load ICB patient data 
-------------------------------------------------------------

ICB_data <- read.csv(here("Outputs","ICB_data.csv"))

# -------------------------------
# 2) Load & prepare spatial data
# -------------------------------

lsoa_shapes <- st_read("C:/Users/Sarah.Glover/OneDrive - NHS/BI Analytics Hub - Specialism Resources/Location Intelligence/Data/LSOAs/SHAPE/LSOA_(Dec_2021)_Boundaries_Full_Clipped_EW_(BFC).shp")

#check geometry exists - should return FALSE

any(is.na(sf::st_geometry(lsoa_shapes)))

# check CRS

sf::st_crs(lsoa_shapes)


# Then transform to WGS84 for Leaflet
lsoa_shapes <- st_transform(lsoa_shapes, 4326)

# Check CRS is now 4326
st_crs(lsoa_shapes)

#check geometry exists now transformed to WGS84
any(is.na(sf::st_geometry(lsoa_shapes)))


# -----------------------------------------------------------------------------------------------------------------------------------------
# 3) Join ICB patient data to the LSOA boundaries and filter LSOAs up to 98% of patients 
# st_make_valid() function from sf package that repairs broken spatial geometries e.g. slivers, polygons with holes touching a boundary etc
# -----------------------------------------------------------------------------------------------------------------------------------------

map_data <- lsoa_shapes  |> 
  left_join(ICB_data, by = c("LSOA21CD" = "LSOA_CODE"))  |>    
  filter(!is.na(PCN_NAME))  |>                                  # keep only rows mapped to a PCN
  filter(substr(LSOA21CD, 1, 1) == "E") |>
  filter(CumulativePercent < 98) |> 
  st_make_valid()


# -----------------------------------------------------------------------------------------------------
# 3A) A basic check if you have LSOAs in one file which have not joined to the LSOA spatial file
# this not needed as I filtered for LSOAs beginning with E only
# -----------------------------------------------------------------------------------------------------    


# Drop geometry from the LSOA layer so we can do a simple attribute join
##lsoa_attrs <- lsoa_shapes |> 
##  st_drop_geometry() |> 
##  select(LSOA21CD, everything())   # <-- replace with your actual LSOA code column

# Safe join that preserves all SSOT rows
##SSOT_join <- SSOT_data |> 
##  left_join(lsoa_attrs, by = c("LSOA_CODE" = "LSOA21CD"))  # <-- change key name if needed

##nrow(SSOT_join)  

# Find which SSOT rows did NOT match any LSOA
##ssot_unmatched <- anti_join(SSOT_data, lsoa_attrs, by = c("LSOA_CODE" = "LSOA21CD"))
##nrow(ssot_unmatched)  

# Inspect a few and tally by PCN/ICB if useful
##head(ssot_unmatched, 10)
##ssot_unmatched |> count(PCN_NAME, sort = TRUE)


# Writing out map_data to see how many patients there are <98%

##library(tidyverse)

##map_data_csv <- map_data |> 
##  st_drop_geometry() 
  
##write_csv(map_data_csv, file = here("Outputs","map_data.csv"))




# ------------------------------------------------------------------------------------------------------
# 4) Inspect names in SHAPE file as writing to shape files truncates column names to only 10 characters
--------------------------------------------------------------------------------------------------------

names(map_data)  
  
name_map <- c(
  "PCN_Patients" = "PCN_Pat",
  "Proportion_PCN" ="PCN_Prop",
  "LSOA_Population" = "LSOA_Pop",
  "Proportion_LSOA" = "LSOA_Prop",
  "CumulativePatients" = "Cum_Pat",
  "CumulativePercent" = "Cum_Perc",
  "SUB_ICB_LOCATION_NAME" = "Sub_ICB"
)


# Apply the mapping only to columns that exist
to_rename <- intersect(names(map_data), names(name_map))
names(map_data)[match(to_rename, names(map_data))] <- name_map[to_rename]

# Double-check no name exceeds 10 characters (for shapefile safety)
stopifnot(all(nchar(names(map_data)) <= 10))


# ---------------------------
# 5) Write to Shapefile
# -----------------------------


write_sf(map_data,"C:/Users/Sarah.Glover/OneDrive - NHS/BI Analytics Hub - Work Requests - Work Requests/A0198/A0198_PCN_Boundary_SSOT_STW/Data/map_data.shp", delete_layer = TRUE)


# ----------------------------------------------------------------------------------------------------------------------------
# 6) Check to see if the objects written to shapefile match the map_data objects ie no lsoa's are missing from the shape file
# ----------------------------------------------------------------------------------------------------------------------------

# out_path <- "C:/Users/Sarah.Glover/OneDrive - NHS/BI Analytics Hub - Work Requests - Work Requests/A0198/A0198_PCN_Boundary_SSOt_STW/Data/map_data.shp"

# tmp <- sf::st_read(out_path, quiet = TRUE)
# nrow(tmp)                         # should match the data file rows
# any(sf::st_is_empty(tmp))         # should be FALSE
# summary(tmp$Shape__Are)           # likely has some NAs now


