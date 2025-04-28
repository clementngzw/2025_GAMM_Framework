#### Chapter 3: A non-linear statistical framework to investigate changes in life history patterns within and among populations ####

## Set the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

## Libraries
## Packages used for wrangling 
library(plyr)
library(tidyverse)
library(data.table)
library(readxl)

# Packages used to check model assumptions
library(zoo)
library(naniar)
library(DHARMa)

# Additive mixed models and model selection
library(mgcv)
library(AICcmodavg)
library(marginaleffects)
library(itsadug)

# Maps
library(ggmap)
library(ggforce)
library(rnaturalearth)
library(cowplot)

cbPalette <- c("#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2", "#F0E442", "#CC79A7", "#999999") 

c. <- function (x) scale(x, scale = FALSE) 

## `Mypairs` function by Alain Zuur and Elena Ieno
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

Mypairs <- function(Z) {
  MyVarx <- colnames(Z)
  pairs(Z, labels = MyVarx,
        cex.labels =  2,
        lower.panel = function(x, y, digits=2, prefix="", cex.cor = 6) {
          panel.cor(x, y, digits, prefix, cex.cor)}, 
        upper.panel =  function(x, y) points(x, y, 
                                             pch = 16, cex = 0.8, 
                                             col = gray(0.1)))
  #print(P)
}

set.seed(100)
par(mfrow = c(2, 2))

#load("2025_otolith_chemistry_v2.RData") 

#################################################################################################################################
#################################################################################################################################

### Data coding ####
## Convert LAICPMS files (.d) to .csv

## Move to the folder containing the compiled data
parent.folder <- "../compiled_data"
sub.folders <- list.dirs(parent.folder, recursive=TRUE)[-1]

files <- list.files(sub.folders, full.names = TRUE, pattern = "\\.d$") 

# Loop through each .d file and convert it to .csv
for (file in files) {
  # Read the file using readr package's read() function
  data <- read.delim(file)
  # Extract the file name without extension
  file_name <- tools::file_path_sans_ext(file)
  # Convert to CSV file by appending the .csv extension
  csv_file <- paste0(file_name, ".csv")
  # Write the data to a CSV file using write.csv function
  write.csv(data, file = csv_file, row.names = FALSE)
  # Print the conversion status
  cat("Converted", file, "to", csv_file, "\n")
}

data

data <- data.frame() 
csv_files <- list.files(sub.folders, full.names = TRUE, pattern = "\\.csv$")

# Loop through each .csv file
for (file in csv_files) {
  # Read the .csv files into a data frame 
  sheet <- read.csv(file, header = TRUE)
  # Add FishID column using file names
  sheet$FishID <- basename(file)
  # Add BatchID column using folder names
  sheet$BatchID <- dirname(file)
  # Coerce time elapsed column to numeric
  sheet$ElapsedTime_s <- as.numeric(sheet$ElapsedTime_s)
  # Create a new column "distance" by multiplying laser speed with time elapsed
  sheet$distance <- 5 * sheet$ElapsedTime_s
  # Stack all files into single dataframe
  data <- rbind(data, sheet)
}
# Note: Amend speed as necessary based on your settings.  

# Remove folder and file names 
data <- data %>% 
  tidyr::separate(col="FishID", sep=".csv", into=c("FishID", "delim1"), remove=TRUE) %>%
  mutate(FishID = gsub(' ', '', FishID)) %>%
  subset(select = -c(delim1))

data <- data %>% 
  tidyr::separate(col="BatchID", sep=c(26, 32), into=c("delim1","Batch", "delim2"), remove=TRUE) %>% 
  subset(select = -c(delim1, delim2))

# Rename columns to E_ppm
data <- data %>% 
  dplyr::rename(Li_ppm = Li_ppm_m7, 
                B_ppm = B_ppm_m11,
                Na_ppm = Na_ppm_m23,
                Mg_ppm = Mg_ppm_m24,
                P_ppm = P_ppm_m31,
                Ca_ppm = Ca_ppm_m44,
                Mn_ppm = Mn_ppm_m55,
                Cu_ppm = Cu_ppm_m63,
                Zn_ppm = Zn_ppm_m66,
                Sr_ppm = Sr_ppm_m88,
                Ba_ppm = Ba_ppm_m137, 
                Pb_ppm = Pb_ppm_m208)

data.table::setDT(data) # Convert to `data.table` format for quicker processing 
data[data < 0] = 0 # Convert all negative values to zero; it will be converted to an infinitesimally small value later. 
data

## Create new data columns, converting PPM to Element:Ca conversion 
data$'Li_Ca' <- ((data$Li_ppm / 6.94)   / (388000 / 40.08) * 1000000) 
data$'B_Ca'  <- ((data$B_ppm  / 10.81)  / (388000 / 40.08) * 1000000) 
data$'Na_Ca' <- ((data$Na_ppm / 22.99)  / (388000 / 40.08) * 1000000) 
data$'Mg_Ca' <- ((data$Mg_ppm / 24.31)  / (388000 / 40.08) * 1000000)
data$'P_Ca'  <- ((data$P_ppm  / 30.97)  / (388000 / 40.08) * 1000000) 
data$'Mn_Ca' <- ((data$Mn_ppm / 54.94)  / (388000 / 40.08) * 1000000)
data$'Cu_Ca' <- ((data$Cu_ppm / 63.55)  / (388000 / 40.08) * 1000000)
data$'Zn_Ca' <- ((data$Zn_ppm / 65.38)  / (388000 / 40.08) * 1000000)
data$'Sr_Ca' <- ((data$Sr_ppm / 87.62)  / (388000 / 40.08) * 100000 )
data$'Ba_Ca' <- ((data$Ba_ppm / 137.33) / (388000 / 40.08) * 1000000)
data$'Pb_Ca' <- ((data$Pb_ppm / 207.20) / (388000 / 40.08) * 1000000)

# Remove data observations that are not from either Lutjanus johnii or L. malabaricus. 
# These were mistakenly identified by our collaborator. 
data <- data %>% filter(FishID != "AB885-01")
data <- data %>% filter(FishID != "AB888-01")
data <- data %>% filter(FishID != "AB889-01")

## Remove the valterite sample
data %>% 
  filter(FishID == "AB575-01") %>%
  ggplot(aes(x=distance, y=Mg_Ca))+ 
  geom_point(size=0.1)+
  geom_line(linewidth= 0.1, alpha= 0.1)+
  theme(legend.position = "none")

data <- data %>% filter(FishID != "AB575-01")

## Format file names 
data <- data %>% mutate(FishID = gsub("-01", "", FishID))
setcolorder(data, c('FishID','Batch','distance','ElapsedTime_s','Ca43_CPS','IntStdWv'))
data

### Internal and external standards ####
standards <- data %>% filter(!grepl('AA|AB', FishID)) 
standards <- standards %>%
  group_by(FishID, Batch) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)
standards

write.csv(standards, 'standards.csv')

data <- data %>% filter(grepl('AA|AB', FishID)) 

transect_length <- data %>% 
  subset(select = c(FishID, distance)) %>%
  group_by(FishID) %>% 
  mutate(distance = max(distance)) %>%
  distinct(FishID, distance)
transect_length

write.csv(transect_length, 'transect_length.csv')

#################################################################################################################################
#################################################################################################################################

## Merge the sample metadata with the element data ####
# There is a pivot table in the Excel document in columns 15-18, which will be removed. 
fish_data <- readxl::read_xlsx("../sample_summary.xlsx", sheet = 4) 
fish_data

colSums(is.na(fish_data)) 
fish_data <- fish_data %>% 
  subset(select = -c(1, 9, 12:18)) %>%
  dplyr::rename(Site = Location)

unique(fish_data$Site)
fish_data <- fish_data %>% 
  mutate(Type = case_when(Site == "Malacca Straits" | Site == "Riau Archipelago" | Site == "Eastern Makassar" | Site == "Sorong" ~ "Market Sampling",
                          TRUE ~ "Research Collaborations"))
fish_data <- fish_data %>% 
  mutate(Region = case_when(Site == "Malacca Straits" | Site == "Kuala Selangor" | Site == "Hutan Melintang" ~ "Malacca Straits",
                            Site == "Eastern Makassar" ~ "Makassar",
                            Site == "Sorong, West Papua" ~ "Sorong", 
                            Site == "Chia Yi Dong / Tung Shi" ~ "Dongshi",
                            Site == "Yi Lan Da Xi" ~ "Yilan",
                            Site == "Long Feng Yu Gang" ~ "Miaoli",
                            Site == "Ping Dong / Tung Heng Chun" ~ "Pingtung",
                            TRUE ~ Site))

fish_data <- fish_data %>% 
  dplyr::rename(FishID = SampleID)

data <- merge(data, fish_data, by = "FishID")
data <- data %>% dplyr::rename('capture_year' = 'Year')
data <- data %>% dplyr::rename('capture_month' = 'Month')
data 

## Check the relationship between left and right otoliths ####
otolith_data <- readxl::read_xlsx("../otolith_data.xlsx", sheet = 1) 
otolith_data <- otolith_data %>% dplyr::rename(TL = 4, SL = 5, WW = 6, L_State = 7, R_State = 8, L_Wt = 9, R_Wt = 10)
otolith_data

otolith_data <- otolith_data %>% 
  mutate(L_Wt = case_when(L_State == "Chipped" | L_State == "Broken" ~ NA, TRUE ~ L_Wt)) %>% 
  mutate(R_Wt = case_when(R_State == "Chipped" | R_State == "Broken" ~ NA, TRUE ~ R_Wt))

## Assess the relationship between left and right otoliths
checkTable <- otolith_data %>% drop_na(L_Wt, R_Wt)
t.test(checkTable$L_Wt, checkTable$R_Wt)

## Assess the correlation between fish length and otolith weight
checkTable <- otolith_data %>% drop_na(SL, R_Wt)
cor.test(checkTable$SL, checkTable$R_Wt, method = "pearson") 

#fwrite(data, file = "example_data.csv")

#################################################################################################################################
#################################################################################################################################

### Check physical samples for damage ####
unique(data$FishID) # List sample IDs

# AA231 -- Broken from around 3970 µm to 4190 µm. 
# AB016 -- The transect path passed through the sulcus, from 197.452 µm to 445.86 µm. 
# AB187 -- The transect path passed through the sulcus, from 164.975 µm to 418.782 µm. 
# AB357 -- The transect path passed through the sulcus. To remove the first 121.202 µm. 
# AB679 -- The transect path passed through the sulcus. To remove the first 323.33 µm.

## Example from AA231
data %>%
  filter(FishID == "AA231") %>%
  ggplot(aes(x = distance, y = Mg_Ca)) + 
  geom_line(linewidth = 0.5, alpha = 0.3) + 
  geom_point(size = 0.5, alpha = 0.3) + 
  theme(legend.position = "none")
## There are a few options on what to do with missing values. You may replace with NAs or interpolate the values. 
## We will address this in a later line. 

## Extract growth information from .txt files 
## All growth increments are measured from the edge (smallest value) to the core and beyond (highest value), noting that the ablation path does not usually start directly on the core. 
## We will amend the age measurements in code. 
growth_data <- list.files(path = "../../Chapter 3_Growth/Images", recursive = TRUE,
                          pattern = "\\.txt$", 
                          full.names = TRUE)
growth_data <- rbindlist(sapply(growth_data, fread, simplify = FALSE), use.names = TRUE, fill = TRUE, idcol = "TxtFileName")

## Rename columns 
growth_data <- growth_data %>% 
  subset(select = -c(ypos, xpos, V1)) %>%
  dplyr::rename("Increment" = "IncNum v1.3c") %>%
  dplyr::rename("FishID" = "TxtFileName") %>%
  dplyr::rename("Width" = "Thickness (µm)")

growth_data <- growth_data %>% 
  tidyr::separate(col="FishID", sep=c(30, 35), into=c("delim1", "FishID", "delim2"), remove=TRUE) %>% 
  subset(select = -c(delim1, delim2))

growth_data$Width[growth_data$Width == 0.000] <- NA
growth_data <- na.omit(growth_data) # remove points on the edge
growth_data

## Maximum increment value represents the distance to remove from the ablation path. 
growth_data <- growth_data %>% 
  group_by(FishID) %>% 
  mutate(Core = case_when(Increment == max(Increment) ~ Width)) %>%
  fill(Core, .direction = "downup") %>%
  ungroup() 

## Remove the additional points beyond the ablation path
growth_data <- growth_data %>% 
  group_by(FishID) %>% 
  mutate(Core_Delete = case_when(Increment == (max(Increment)) ~ "Delete")) %>%
  ungroup()
growth_data <- growth_data[- grep("Delete", growth_data$Core_Delete),]
growth_data <- subset(growth_data, select = -c(Core_Delete))

## Extract edge increment
growth_data <- growth_data %>% 
  group_by(FishID) %>% 
  mutate(Edge = case_when(Increment == (min(Increment)) ~ Width)) %>%
  fill(Edge, .direction = "downup") %>%
  ungroup()

## Convert Age to the correct order.
# This is largely dependent on your growth measurement methodology. 
growth_data <- growth_data %>% 
  group_by(FishID) %>%
  arrange(desc(Increment), .by_group = TRUE) %>% 
  mutate(Age = 1:n()) %>%
  ungroup()
growth_data

### Merge the growth information with sample metadata
growth_data <- merge(growth_data, fish_data[, c('FishID','Species','Year','Month','Region')], by = "FishID")
growth_data <- growth_data %>% rename('capture_year' = 'Year')
growth_data <- growth_data %>% rename('capture_month' = 'Month')
growth_data <- growth_data %>%
  mutate(formation_month = case_when(Species == "Lutjanus johnii" ~ "9", 
                                     Species == "Lutjanus malabaricus" ~ "10"))

## Assign a corrected Year variable based on the capture year and formation period
growth_data$Increment <- (growth_data$Increment - 1) * - 1 
growth_data <- growth_data %>% 
  group_by(FishID) %>% 
  mutate(Year = case_when((as.numeric(capture_month) <= as.numeric(formation_month)) ~ (capture_year + Increment), 
                          (as.numeric(capture_month) > as.numeric(formation_month)) & Age != max(Age) ~ (capture_year + Increment + 1), 
                          (as.numeric(capture_month) > as.numeric(formation_month)) & Age == max(Age) ~ (capture_year + Increment))) %>%
  ungroup()
growth_data <- subset(growth_data, select = -c(Increment, capture_year, capture_month)) # Year measurements have now been corrected. 

## Cohort variable 
growth_data <- growth_data %>% 
  group_by(FishID) %>% 
  mutate(Cohort = min(Year) - min(Age))

colnames(growth_data)
setcolorder(growth_data, c('FishID','Age','Year','Species','Cohort','Region','Width','Core','Edge'))
growth_data

## Check number of growth bands per Age class by Region
growth_data %>% 
  filter(Species == "Lutjanus johnii" & Region == "Riau Archipelago" | Region == "Sorong") %>%
  group_by(Region, Age) %>% # replace with Year as required
  summarise(n = n()) %>% 
  print(n = Inf)

growth_data %>% 
  filter(Species == "Lutjanus malabaricus" & Region == "Riau Archipelago" | Region == "Sorong") %>%
  group_by(Region, Age) %>% # replace with Year as required
  summarise(n = n()) %>%
  print(n = Inf)

#################################################################################################################################
#################################################################################################################################

### Remove measurements outside the ablation transect 
core_data <- growth_data %>% 
  subset(select = c("FishID","Core")) %>% 
  distinct(FishID, Core)
data1 <- merge(data, y = core_data[, c("FishID","Core")], by = "FishID", all.x=TRUE)
data1 <- data1 %>% 
  group_by(FishID) %>% 
  filter(distance > Core) %>%
  ungroup()
data1

## Remove rows where there are multiple NA observations 
colSums(is.na(data1))
data1 <- data1 %>% filter(rowSums(across(18:28, ~ is.na(.x))) < 10) 

## Adjust columns to start from zero within each group
## c_distance represents the corrected distance. 
data1 <- data1 %>%
  group_by(FishID) %>%
  mutate(across(distance, ~ . - first(.), .names = "c_distance")) %>%
  ungroup()

## Create a cumulative sum distance column for the growth data
growth_data <- growth_data %>% 
  group_by(FishID) %>% 
  arrange(FishID, Age) %>% 
  mutate(c_distance = cumsum(Width))
growth_data$diff <- growth_data$c_distance # We will need to duplicate this for rolling join alignment
data1$c_distance <- round(data1$c_distance, 3)
data1

#write.csv(growth_data, "growth_data.csv")

## Replace areas with cracks or damage with NAs
data1 <- data1 %>% mutate(across(7:29, ~ if_else(distance >= 3970.000 & distance <= 4190.000 & FishID == "AA231", NA_real_, .)))
data1 <- data1 %>% mutate(across(7:29, ~ if_else(c_distance >= 197.452 & c_distance <= 445.86 & FishID == "AB016", NA_real_, .)))
data1 <- data1 %>% mutate(across(7:29, ~ if_else(distance >= 164.975 & distance <= 418.782 & FishID == "AB187", NA_real_, .)))
data1 <- data1 %>% mutate(across(7:29, ~ if_else(c_distance <= 121.202 & FishID == "AB357", NA_real_, .)))
data1 <- data1 %>% mutate(across(7:29, ~ if_else(c_distance <= 323.33 & FishID == "AB679", NA_real_, .)))

## Create nearest rolling join variable in the data frame
setDT(data1) ; setDT(growth_data)
data1 <- growth_data[,c("FishID","c_distance","Age","Year","Cohort","diff")][data1, roll = "nearest", on = .(FishID, c_distance)]
data1$difference <- abs(data1$c_distance - data1$diff) # The difference variable represents the distance from the first increment. 

## Minimise the difference between the location of each Age and the ablation spot
data1 <- data1 %>% 
  group_by(FishID, Age) %>% 
  dplyr::mutate(Age = case_when(difference == min(difference) ~ Age))
data1 <- data1 %>% mutate(Year = case_when(Age != 'NA' ~ Year, TRUE ~ NA)) # We will interpolate Age and Year in the next chunk.
data1 <- as.data.frame(data1)
data1

### Create rolling median and mean to remove outliers and potential data fuzziness ####
median_size <- 5
data1 <- data1 %>%
  group_by(FishID) %>%
  mutate(sLi_Ca = rollmedian(Li_Ca, k = median_size, align = "center", fill = NA),
         sB_Ca  = rollmedian(B_Ca,  k = median_size, align = "center", fill = NA),
         sNa_Ca = rollmedian(Na_Ca, k = median_size, align = "center", fill = NA),
         sMg_Ca = rollmedian(Mg_Ca, k = median_size, align = "center", fill = NA),
         sP_Ca  = rollmedian(P_Ca,  k = median_size, align = "center", fill = NA),
         sMn_Ca = rollmedian(Mn_Ca, k = median_size, align = "center", fill = NA),
         sCu_Ca = rollmedian(Cu_Ca, k = median_size, align = "center", fill = NA),
         sZn_Ca = rollmedian(Zn_Ca, k = median_size, align = "center", fill = NA),
         sSr_Ca = rollmedian(Sr_Ca, k = median_size, align = "center", fill = NA),
         sBa_Ca = rollmedian(Ba_Ca, k = median_size, align = "center", fill = NA),
         sPb_Ca = rollmedian(Pb_Ca, k = median_size, align = "center", fill = NA))

mean_size <- 10
data1 <- data1 %>%
  group_by(FishID) %>%
  mutate(sLi_Ca = rollmean(sLi_Ca, k = mean_size, align = "center", fill = NA),
         sB_Ca  = rollmean(sB_Ca,  k = mean_size, align = "center", fill = NA),
         sNa_Ca = rollmean(sNa_Ca, k = mean_size, align = "center", fill = NA),
         sMg_Ca = rollmean(sMg_Ca, k = mean_size, align = "center", fill = NA),
         sP_Ca  = rollmean(sP_Ca,  k = mean_size, align = "center", fill = NA),
         sMn_Ca = rollmean(sMn_Ca, k = mean_size, align = "center", fill = NA),
         sCu_Ca = rollmean(sCu_Ca, k = mean_size, align = "center", fill = NA),
         sZn_Ca = rollmean(sZn_Ca, k = mean_size, align = "center", fill = NA),
         sSr_Ca = rollmean(sSr_Ca, k = mean_size, align = "center", fill = NA),
         sBa_Ca = rollmean(sBa_Ca, k = mean_size, align = "center", fill = NA),
         sPb_Ca = rollmean(sPb_Ca, k = mean_size, align = "center", fill = NA))

data1 <- data1 %>% subset(select = -c(12:34)) 
data1 <- as.data.frame(data1)

setcolorder(data1, c('FishID','Batch','ElapsedTime_s','c_distance','Age','Year','diff','difference','distance',
                     'Species','Country','Region','Site','Cohort','capture_year','capture_month','TL','SL','WW','Type',
                     'Ca43_CPS','IntStdWv','Core'))
data1 <- data1 %>% filter(rowSums(across(24:34, ~ is.na(.x))) < 10) 
data1

### Interpolating Age based on the missing NA values present ####
## Define the number of decimal places in both columns
data1$Age <- round(data1$Age, 10)
data1$Year <- round(data1$Year, 10)

## Replace the minimum Age and Cohort value
data1 <- data1 %>%
  group_by(FishID) %>%
  mutate(Age = ifelse(c_distance == min(c_distance), 0, Age)) %>%
  mutate(Year = ifelse(c_distance == min(c_distance), min(Cohort), Year)) %>%
  ungroup()

## Age interpolation
data1 <- data1 %>%
  group_by(FishID) %>%
  mutate(c_Age = zoo::na.approx(Age, na.rm = FALSE)) %>%
  ungroup()

## Year interpolation
data1 <- data1 %>%
  group_by(FishID) %>%
  mutate(c_Year = zoo::na.approx(Year, na.rm = FALSE)) %>%
  ungroup()

## Year fill for each increment
data1 <- data1 %>%
  group_by(FishID) %>%
  fill(Year, .direction = "downup") %>%
  ungroup() 

colSums(is.na(data1))
data1 <- data1 %>% 
  drop_na(c_Age, c_Year, sLi_Ca, sB_Ca, sNa_Ca, sMg_Ca, sP_Ca, sMn_Ca, sCu_Ca, sZn_Ca, sSr_Ca, sBa_Ca, sPb_Ca)

data1 <- data1 %>% 
  subset(select = -c(Age, diff, difference, distance)) %>%
  dplyr::rename(Age = c_Age)

setcolorder(data1, c('FishID','Batch','ElapsedTime_s','c_distance','Age','Year','c_Year','Cohort',
                     'Species','Country','Region','Site','TL','SL','WW','Type',
                     'capture_year','capture_month',
                     'Ca43_CPS','IntStdWv','Core'))
data1

#################################################################################################################################
#################################################################################################################################

### Data preparation 
## For this case example, we will be focusing solely on samples from the following Regions, Malacca Straits, Riau Archipelago, and Sorong. 
## The full data frame is data1. 
element_data <- data1
element_data <- element_data %>% 
  filter(Region == "Malacca Straits" | Region ==  "Riau Archipelago" | Region ==  "Sorong")
element_data[c('FishID','Species','Country','Region','Site')] <- lapply(element_data[c('FishID','Species','Country','Region','Site')], as.character)
element_data[c('FishID','Species','Country','Region','Site')] <- lapply(element_data[c('FishID','Species','Country','Region','Site')], as.factor) 
## Categorical variables must be coded `as.factor()` in the mgcv package. 

## Replace zeros; this indicates that it is not a true absence of E:Ca but that measurements below the detection limit.
element_data[22:32][element_data[22:32] == 0] <- 0.0000000001 
element_data

## Visualise Strontium data 
element_data %>% 
  ggplot(aes(x=c_distance, y=sSr_Ca, colour=FishID, group=FishID))+ 
  geom_point(size=0.1, alpha= 0.1, pch='.')+
  geom_line(linewidth= 0.1, alpha= 0.1)+
  facet_wrap(Species~Region)+
  theme(legend.position = "none")

### Create separate data frames for each species
Lm_data <- subset(element_data, element_data$Species == "Lutjanus malabaricus") 
Lj_data <- subset(element_data, element_data$Species == "Lutjanus johnii") 

### Lutjanus malabaricus
setDT(Lm_data) 
colSums(is.na(Lm_data)) ; naniar::gg_miss_var(Lm_data)  

## Create new variables
Lm_data$fAge <- as.factor(Lm_data$Age)
Lm_data$fYear <- as.factor(Lm_data$Year)
Lm_data$region_year <- as.factor(paste(Lm_data$Region, Lm_data$Year, sep = "_")) 

## Truncate data based on the sample depth per Year and Region class bins
growth_data %>% 
  filter(Species == "Lutjanus malabaricus" & Region == "Riau Archipelago" | Region == "Sorong") %>%
  group_by(Region) %>% 
  summarise(n = n()) %>%
  print(n = Inf)
Lm_data <- Lm_data %>% filter(!((Region == "Riau Archipelago" & Year <= "2016") | (Region == "Sorong" & Year <= "2006")))
Lm_data <- Lm_data %>% filter(Age <= 12) 

## Number of measurements per individual
Lm_data %>%
  group_by(FishID) %>%
  dplyr::summarise(count = n()) %>%
  print(n = Inf)

## Effect of FishID 
Lm_data %>%
  ggplot(aes(y = sSr_Ca, x = FishID)) +
  geom_boxplot(alpha = 0.5) + 
  geom_hline(yintercept = mean(Lm_data$sSr_Ca), alpha = 0.2) +
  theme_bw() # There is an individual effect

## Plot the relationship between E:Ca response variable and the continuous covariate of interest
## Continuous covariates: Age, Year
Mypairs(Lm_data[,c("Age", "Year")]) 

Lm_data %>% 
  ggplot(aes(x=c_distance, y=sBa_Ca, colour=FishID, group=FishID))+ 
  geom_point(size=0.5, alpha= 0.1)+
  geom_line(linewidth= 0.1, alpha= 0.1)+
  facet_wrap(.~Region)+
  theme(legend.position = "none")

Lm_data %>%
  ggplot(aes(y = sSr_Ca, x = Age, colour = Region)) +
  geom_point(pch = '.') + 
  geom_smooth(aes(group = Region), method = 'gam', colour = 'black') + 
  theme_bw()

Lm_data %>%
  ggplot(aes(y = sSr_Ca, x = Year, colour = Region)) +
  geom_point(pch = '.') + 
  geom_smooth(aes(group = Region), method = 'gam', colour = 'black') + 
  theme_bw()

## Categorical covariates: Region (two levels: Sorong, Riau Archipelago)
table(Lm_data$Region)

colnames(Lm_data)
setcolorder(Lm_data, c('FishID','Batch','ElapsedTime_s','c_distance','Age','fAge','Year','fYear','c_Year','Cohort',
                       'Species','Country','Region','region_year','Site','TL','SL','WW','Type',
                       'capture_year','capture_month',
                       'Ca43_CPS','IntStdWv','Core'))

## Summary statistics 
Lm_AAC <- growth_data %>% 
  filter((Species == "Lutjanus malabaricus" & Region == "Riau Archipelago") | (Species == "Lutjanus malabaricus" & Region == "Sorong")) %>% 
  group_by(FishID) %>% 
  mutate(AAC = max(Age), Cap_Year = (Cohort + AAC))

Lm_AAC %>%
  distinct(FishID, Region, AAC) %>%
  group_by(Region) %>%
  summarise(mean = mean(AAC), sd = sd(AAC), min = min(AAC), max = max(AAC), n = n_distinct(FishID))

Lm_AAC %>% 
  filter(Region != "NA") %>%
  group_by(AAC, Region) %>% 
  summarise(n=n_distinct(FishID)) %>% 
  ggplot(., aes(x = AAC, y = n))+ 
  geom_bar(stat = 'identity',color = 'grey5',fill = 'grey45') +
  scale_y_continuous(limits = c(0, 14)) + 
  scale_x_continuous(limits = c(0, 23)) +
  facet_wrap(.~Region, ncol = 2) +
  labs(y = 'Count', x = 'Age-at-Capture') + 
  theme_bw() +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 

Lm_AAC <- as.data.frame(Lm_AAC)
Lm_AAC <- Lm_AAC %>% filter(Age == "1") %>% arrange(Region, AAC)

s_plot1 <- Lm_AAC %>% 
  filter(Region == "Riau Archipelago") %>%
  mutate(OrderFactor = factor(FishID, levels = unique(.$FishID[order(.$Cap_Year, desc(.$Cap_Year - .$Cohort))]))) %>% 
  ggplot(aes(y = Year, x = OrderFactor)) +
  geom_linerange(aes(ymin = Cohort, ymax = Cap_Year)) + 
  scale_y_continuous(limits = c(2010, 2023), breaks = seq(2010, 2023, by = 4), labels = seq(2010, 2023, by = 4)) +
  facet_wrap(.~Region) +
  coord_flip() + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        strip.text = element_text(size = 16))
s_plot1

ggsave("s_plot1.png", plot = s_plot1, dpi = 300, width = 210, height = 297, units = "mm")

s_plot2 <- Lm_AAC %>% 
  filter(Region == "Sorong") %>%
  mutate(OrderFactor = factor(FishID, levels = unique(.$FishID[order(.$Cap_Year, desc(.$Cap_Year - .$Cohort))]))) %>% 
  ggplot(aes(y = Year, x = OrderFactor)) +
  geom_linerange(aes(ymin = Cohort, ymax = Cap_Year)) + 
  scale_y_continuous(limits = c(2000, 2023), breaks = seq(2000, 2023, by = 5), labels = seq(2000, 2023, by = 5)) +
  facet_wrap(.~Region) +
  coord_flip() + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 16))
s_plot2

ggsave("s_plot2.png", plot = s_plot2, dpi = 300, width = 210, height = 297, units = "mm")

#################################################################################################################################
#################################################################################################################################

### Aim: Investigate drivers of E:Ca changes in tropical snappers across the Indo-Pacific Region
# Additive mixed models partition variance in Region and Age.
# Smoothers are typically specified using the default bs = 'tp' function. 
# k represents the maximum allowable degrees of freedom. While it is possible to set a high k to explain the wiggliness of the data, it comes at a computational cost. 
# Note: Zuur also mentioned that it is hard to explain such complex, wiggly patterns above k = 10. 

# As mentioned by Gavin Simpson in a Stack Overview post, 
# "The general strategy should be to set k to be as large as needed to create a basis rich enough that the 
# true (but unknown) function or a close approximation to that true function is representable by the basis."

# We are using the Gamma distribution with a log-link function because we have continuous positive values in our data. 
# Zero values have also been replaced with 0.0000000001, which indicates that it is not a true absence of E:Ca but that measurements are below the detection limit.

# Random effects are specified using the bs = 're' function. 
# Due to the size of our data frame, we will also be using the bam() function rather than gam() to lower the computational time. 
# We are using fREML to speed up comparisons. 
# select = TRUE penalises the linear components of each smoother term.

# `discrete = ` discretises covariates into a number of bins to speed up computational time. 
# Our models with ~150,000 observations do not take a long time to run, and I have encountered issues with model prediction when `discrete = TRUE`.

# You may also set a penalty, m-1, for the squared first derivative of the function, rather than the second derivative. 
# This reduces collinearity between the global smoother and the group-specific terms which may lead to high uncertainty around the global smoother.

### Model construction ####
## Important considerations during model construction, from Pedersen et al., (2019)
# Should each group have its own smoother, or will a common (global) smoother suffice?
# Do all of the group-specific smoothers have the same wiggliness, or should each group have its own smoothing parameter?
# Will the smoothers for each group have a similar shape to one another a shared global smoother?

# See Pedersen, E.J., Miller, D.L., Simpson, G.L. and Ross, N., 2019. Hierarchical generalized additive models in ecology: an introduction with mgcv. PeerJ, 7, p.e6876.
# https://doi.org/10.7717/peerj.6876 

# Also see useful discussion on Stack Exchange by Dr. Gavin Simpson,
# https://stats.stackexchange.com/questions/637423/conceptual-interpretation-of-bs-fs-and-by-term-in-gam. 

### Strontium ####

## Base model 1
# In our base model, we allowing for dependency in E:Ca response among each individual across Ages (i.e., random Age slope by FishID intercept).
# The base model includes a global Age smoother (G model) without Region-specific patterns (Pedersen et al., 2019). 
Lm_Sr_M1a <- bam(sSr_Ca ~ 
                   s(Age, k = 24, bs = "tp", m = 2) +
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

# We are fitting a thin plate regression spline (other splines can also be used, like cubic regression splines).
# k = 10, which allows for a maximum estimated degrees of freedom of 10, as a high k value may complicate our ability to interpret patterns in the data. 
# k should be defined based on a biological understanding of the system. 

# bs = 're' creates a random effect.
# Random intercepts adjust the height of other model terms with a constant value: s(fac, bs = "re")
# Random slopes adjust the slope of the trend of a numeric predictor: s(fac, x0, bs = "re")

# select = TRUE penalises the linear components of each smoother term.
# method = 'fREML' produces fast Restricted Maximum Likelihood Estimations. 
# family = 'Gamma(link = "log")' calls the gamma distibution in `mgcv` with a log-link function. 

## Check the numerical output of the data
summary(Lm_Sr_M1a)
anova(Lm_Sr_M1a)
plot(fitted(Lm_Sr_M1a), residuals(Lm_Sr_M1a)) # check residual plots 

## Check diagnostic information for our baseline model
# According to ?gam.check(), the p-value is computed by simulation: the residuals are randomly re-shuffled `k.rep` times to obtain the null distribution of the differencing variance estimator, if there is no pattern in the residuals. 
set.seed(100) ; k.check(Lm_Sr_M1a)
set.seed(100) ; gam.check(Lm_Sr_M1a) 

## Generate preliminary plots
plot.gam(Lm_Sr_M1a, all.terms = T)
p <- marginaleffects::plot_predictions(Lm_Sr_M1a, condition = c("Age"), type = 'response') 
p

## Model 2 
# In our second model, we will construct a G model (Pedersen et al., 2019). 
# Here, we allowing for dependency in E:Ca response among each individual across Ages. 
# We are also interested in differences in E:Ca among Regions, and this is fitted as a random effect. It can be fitted as a factor. 
# The base model includes a global Age smoother (G model) and Region random effect term to allow for differences among sampling locations (Pedersen et al., 2019). 
Lm_Sr_M1b <- bam(sSr_Ca ~ 
                   s(Age, k = 24, bs = "tp", m = 2) +
                   s(Region, bs = 're') + 
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model 3
# In our third model, we will construct a GS model (Pedersen et al., 2019). 
# There is a global (common) Age smoother, plus group-level Region Age smoothers (factor smooth) with the same wiggliness. 
# Essentially, it tests the amount of variation in E:Ca explained by the common Age smoother, plus fully penalised group-level Region smooth functions of Age, 
# against a null of no effect. 
Lm_Sr_M1c <- bam(sSr_Ca ~ 
                   s(Age, k = 24, bs = "tp", m = 2) +
                   s(Age, Region, k = 24, bs = "fs", m = 2) +
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

# The `bs = "fs"` term in s(Age, Region) creates a smooth for each Region level using k (n = 10) basis functions. 
# This factor smooth is fully penalized, meaning that it has no penalty null space (no functions for which the penalty has no effect). 
# As part of the penalty construction, each group (factor level) will have its own intercept. 

# See the comment by Gavin Simpson in https://stats.stackexchange.com/questions/637423/conceptual-interpretation-of-bs-fs-and-by-term-in-gam . 
# See also https://stats.stackexchange.com/questions/482509/in-mgcvgam-what-is-the-difference-between-sx-by-cat-and-sx-cat-bs-fs. 

## Model 4
# In our fourth model, we will construct a GI model structure for Age (Pedersen et al., 2019).
# A Region random intercept (or factor) must be included as group-specific intercepts (means) are not incorporated into factor `by` variable smoothers. 
# Here, we are interested to know whether there is a global Age smoother plus group-level Region smoothers with the different wiggliness. 
Lm_Sr_M1d <- bam(sSr_Ca ~ 
                   s(Age, k = 24, bs = "tp", m = 2) +
                   s(Age, by = Region, k = 24, bs = "tp", m = 2) +
                   s(Region, bs = 're') + 
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

# Unlike `bs = "fs"`, the `by` smooth variant creates a separate smoothing parameter for each level of Region. 
# This means that one smoothing parameter for Region may be wiggly, while the other may produce a smooth pattern. 

## Important note:
# Pedersen et al., (2019) penalises the squared first derivative of the s(Age, by = Region) term to reduce collinearity between the global smoother and 
# group-specific terms, which has been reported to produce high uncertainty around the global smoother. The factor `by` smooth has a zero penalty null space 
# once the constant term is removed by the sum to zero constraint.

# Similarly, we noted that Model GI produced high uncertainty in the estimates when we penalised the squared first derivatives for group-level smoothers using m = 1. 
# Erroneous predictions occur when the data is overfitted to a sparse number of data points. Importantly, note that the linear function is in the span of the basis 
# but it is unpenalized, and it cannot be shrunk back to a zero function. 
# See https://stats.stackexchange.com/questions/625419/choice-of-m-order-of-derivative-for-mgcv-splines. 

## Model 5
# In our fifth model, we will construct an S model (Pedersen et al., 2019).
# This model describes a model with group-specific Age smoothers with identical wiggliness.
# Essentially, it tests the amount of variation in E:Ca explained by fully penalised group-level Region smooth functions of Age against a null of no effect. 
Lm_Sr_M1e <- bam(sSr_Ca ~ 
                   s(Age, Region, k = 24, bs = "fs", m = 2) +
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

# One may not be particularly interested in the effects of individual Regions, which is assessed using the `by` smooth. Instead, we think of this as a single 
# "smooth" as a single Age plot with each Region smooth superimposed to get some idea of the variation within and between individual functions.

## Model 6
# Our sixth model is model I (Pedersen et al., 2019).
# This model describes a model with group-specific smoothers with different wiggliness.
# Unlike `bs = "fs"`, the `by` smooth variant creates a separate smoothing parameter for each level of Region. We are also assuming there is no global effect present. 
Lm_Sr_M1f <- bam(sSr_Ca ~ 
                   s(Age, by = Region, k = 24, bs = "tp", m = 2) +
                   s(Region, bs = 're') + 
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model selection using Akaike's Information Criterion ####
bbmle::AICtab( Lm_Sr_M1a, Lm_Sr_M1b, Lm_Sr_M1c, Lm_Sr_M1d, Lm_Sr_M1e, Lm_Sr_M1f, 
               base=T,logLik=T,weights=T ) 
## Optimal model structure is Lm_Sr_M1d. 
summary(Lm_Sr_M1d)

# As per Gavin Simpson, AICc should not be used for GAMMs as there is an additional penalty for non parametric methods. 
# See https://stats.stackexchange.com/questions/552880/by-group-random-effect-gam. 

## Assess for temporal patterns in E:Ca 
## Fit a discrete Year slope and Region random intercept. 
Lm_Sr_M2 <- bam(sSr_Ca ~ 
                  s(Age, k = 24, bs = "tp", m = 2) +
                  s(Age, by = Region, k = 24, bs = "tp", m = 2) +
                  s(Region, bs = 're') + 
                  s(Region, Year, bs = 're') + 
                  s(FishID, bs = 're') + 
                  s(FishID, Age, bs = 're'), 
                select = TRUE, 
                method = "fREML", 
                family = Gamma(link = "log"),
                data = Lm_data) 

bbmle::AICtab( Lm_Sr_M1d, Lm_Sr_M2, base=T,logLik=T,weights=T ) 
## Model selection indicates that the null model without temporal effects provides a marginally better fit. 
## Complex models are penalised by AIC model selection criteria. 

## Construct the optimal model for Strontium for Lutjanus malabaricus ####
Lm_Sr <- Lm_Sr_M1d

## Check the numerical output of the data
summary(Lm_Sr)
anova(Lm_Sr)
plot(fitted(Lm_Sr), residuals(Lm_Sr)) # check residual plots 

## Check diagnostic information for our baseline model
# According to ?gam.check(), the p-value is computed by simulation: the residuals are randomly re-shuffled `k.rep` times to obtain the null distribution of the differencing variance estimator, if there is no pattern in the residuals. 
set.seed(100) ; k.check(Lm_Sr)
set.seed(100) ; gam.check(Lm_Sr) 

## Generate preliminary plots
plot.gam(Lm_Sr, all.terms = T)

## Plot conditional effects of Age for each Region, while setting the group-specific random effects or offsets to zero. 
## This contrasts with marginal effect plots, where one is interested in the Age effect for each Region, while integrating or holding the random effect values to a constant. 
p <- marginaleffects::plot_predictions(Lm_Sr, condition = c("Age","Region"), type = 'response') + theme_bw()
p

#################################################################################################################################
#################################################################################################################################

## Assess the numerical results from the optimal GAMM ####
## R code from Highland Statistics Ltd

## Extract the r parameter from the fitted gamma GAMM model
summary(Lm_Sr) # See scale estimated from `summary()` 
r <- 1 / 0.0094024
r 
#  sSr_Ca ~ Gamma(mu_ijk, r =  106.3558) 
#  E(sSr_Ca)   = mu_ijk
#  var(sSr_Ca) = mu_ijk ^2 / 106.3558

## Attain the sigma of the random effect of FishID
gam.vcomp(Lm_Sr) 

## Perform model validation for the Gamma GAMM 
Lm_data$E  <- resid(Lm_Sr, type = "pearson")  # Pearson residual values
Lm_data$mu <- fitted(Lm_Sr)                   # Gamma fitted values

## Plot residuals versus fitted values
Lm_data %>% 
  filter(Age < 12) %>%
  ggplot(aes(y = E, x = mu)) + 
  geom_point(size = 0.5, alpha = 0.3) +
  xlab("Fitted values") + ylab("Residuals") +
  theme_bw()

# Plot observed versus fitted values
Lm_data %>% 
  ggplot(aes(y = sSr_Ca, x = mu)) +
  geom_point(size = 0.5, alpha = 0.5) +
  xlab("Fitted values") + ylab("sSr_Ca") + 
  theme_bw()

## Plot residuals against Age data 
Lm_data %>%
  ggplot(aes(y = E, x = Age, colour = Region)) + 
  geom_point(size = 0.5, alpha = 0.5) + 
  geom_smooth(colour = 'black') +
  facet_wrap(.~Region) +
  xlab("Age") + ylab("Residuals") + 
  theme_bw()

### Simulate Residuals using the DHARMa package 
## This allows us to check if there are any problems with the data. 
Lm_Sr_Residuals <- simulateResiduals(fittedModel = Lm_Sr, plot = FALSE)

### Assess for uniformity and check for outliers ####
## According to the package vignette for DHARMa, Florian Hartig states that random effect estimates are excluded in the predictions when plotting 
## residual versus predicted data. He adds that DHARMa residuals often show a slight pattern in the residuals even if the model is correctly specified, 
## and tests for this can get significant for large sample sizes. 

## See section, 'General remarks on interperting residual patterns and tests'
## https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#general-remarks-on-interperting-residual-patterns-and-tests 

par(mfrow = c(1,1))
plotQQunif(Lm_Sr_Residuals, 
           testUniformity = TRUE, 
           testOutliers = TRUE, 
           testDispersion = TRUE) 
## Overall, our QQ plot looks okay. DHARMa flags errors for uniformity, outliers, and dispersion. 
## However, this is likely due to the large sample size that is present in our data. 

DHARMa::plotResiduals(Lm_Sr_Residuals, quantreg = TRUE, smoothScatter = TRUE) 
## There are significant quantile deviations and presence of outliers. But this is not a cause for concern given our large sample size. 

plotResiduals(Lm_Sr_Residuals, form = Lm_data$Age, smoothScatter = TRUE) 
# Looks good.

#################################################################################################################################
#################################################################################################################################

### Check for dependency in the random effects ####
## R code from Highland Statistics Ltd

## We need to extract the random effects and check whether they are independent and normally distributed.
coef(Lm_Sr)

# The idea is to specify a data frame with only unique values for s(FishID), and everything else is set to something constant. 
# We will then use the predict() function to predict the values of the random effects and SEs. 
data2 <- data.frame(FishID = levels(Lm_data$FishID))

## Covariates can be set to something irrelevant. 
data2$Age <- mean(Lm_data$Age)
data2$Region <- "Riau Archipelago"

## We then use the predict function with type = "terms"
p.Lm_Sr <- predict(Lm_Sr, data2, 
                   type = "terms", se.fit = TRUE)
head(p.Lm_Sr$fit)

## Extract the estimated random effects and their standard errors,
data2$ai <- p.Lm_Sr[["fit"]][ , "s(FishID)"]
data2$ai.se <- p.Lm_Sr[["se.fit"]][ , "s(FishID)"]

## Plot in a dotchart
data2$rowID <- 1:nrow(data2)

p <- data2 %>% 
  ggplot(aes(x = ai, y =rowID)) +
  geom_point(size = 1, col = grey(0.5)) +
  geom_errorbarh(aes(y = rowID, xmax = ai + 1.96 * ai.se, xmin = ai - 1.96 * ai.se), col = "blue", height = 0.05, alpha = 0.2) +
  geom_vline(xintercept = 0, lty = 2, col = "red") +
  xlab("Estimated random effect") + ylab("FishID") +
  scale_x_continuous(limits = c(-0.75, 0.75), breaks = c(-0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75), labels = c(-0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75)) +
  theme_bw()
p

#################################################################################################################################
#################################################################################################################################

### Preliminary plot using the `itsadug` package
par(mfrow = c(1, 1))
itsadug::plot_smooth(Lm_Sr, view = "Age", plot_all = c("Region"), rm.ranef = FALSE) # without extracting random effects
itsadug::plot_smooth(Lm_Sr, view = "Age", plot_all = c("Region"), transform = exp, rm.ranef = FALSE) # transform to original scale
# Note that this is equivalent to plot_predictions(..., type = 'response'). 

### Model visualisation for the optimal GAMM model ####
## Plot changes in Sr across Age for each Region for the average Year. 
Lm_Sr_Age_Data <- ddply(Lm_data, .(Region), summarize, 
                        Age = seq(min(Age), max(Age), length = 1000)) 
Lm_Sr_Age_Data$FishID <- "NA"
Lm_Sr_Age_Data$FishID <- factor(Lm_Sr_Age_Data$FishID)

## Predict how Sr varies with Age for each Region after removing random effects
sapply(Lm_Sr$smooth,"[[","label") 
Lm_Sr_PredictAge <- predict.bam(Lm_Sr, 
                                newdata = Lm_Sr_Age_Data, 
                                exclude = c("s(Region)",
                                            "s(FishID)",
                                            "s(FishID,Age)"),
                                se.fit = TRUE)
## Important note: 
# Each smooth term should be excluded in the same way it was labelled in summary(), not during model construction. 

Lm_Sr_Age_Data$mu <- exp(Lm_Sr_PredictAge$fit)
Lm_Sr_Age_Data$ul <- exp(Lm_Sr_PredictAge$fit + 1.96 * Lm_Sr_PredictAge$se.fit )
Lm_Sr_Age_Data$ll <- exp(Lm_Sr_PredictAge$fit - 1.96 * Lm_Sr_PredictAge$se.fit ) 

Lm_Sr_Age_Data <- Lm_Sr_Age_Data %>% filter(Age <= 12)
head(Lm_Sr_Age_Data, 5)

### Plot marginal effect of Age:Region on Sr:Ca, setting all random effect values to a constant value. #### 
plot1 <- Lm_Sr_Age_Data %>%
  ggplot(aes(x = Age, y = mu)) + 
  geom_line(aes(colour = Region)) + 
  geom_ribbon(aes(ymax = ul, ymin = ll, fill = Region), alpha = 0.5) + 
  scale_color_manual(values=c("#E69F00","#0072B2")) + 
  scale_fill_manual(values=c("#E69F00","#0072B2")) +
  labs(y = 'Predicted Sr:Ca (mmol/mol)') + 
  scale_x_continuous(limits = c(0, 12), breaks = c(0,3,6,9,12), labels = c("0","3","6","9","12")) +
  scale_y_continuous(limits = c(180, 800), breaks = c(200, 400, 600, 800), labels = c("200", "400", "600", "800")) +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14),
        legend.background = element_rect(fill='white',linewidth=0.5,linetype="solid",colour ="grey15"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1.0, 1.0, 1.0, 1.0, "cm"))
plot1 

ggsave("plot1.png", plot = plot1, dpi = 300, width = 297, height = 210, units = "mm")

#################################################################################################################################
#################################################################################################################################

### Barium ####
## Construct different models to explore the relationship between Barium with Age

## Model 1- Global Ba pattern 
Lm_Ba_M1a <- bam(sBa_Ca ~ 
                   s(Age, k = 24, bs = "tp", m = 2) +
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model 2- Global Ba pattern with Region random effect
Lm_Ba_M1b <- bam(sBa_Ca ~ 
                   s(Age, k = 24, bs = "tp", m = 2) +
                   s(Region, bs = 're') + 
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model 3- Global Age pattern with group-level Region smoothers with identical wiggliness
Lm_Ba_M1c <- bam(sBa_Ca ~ 
                   s(Age, k = 24, bs = "tp", m = 2) +
                   s(Age, Region, k = 24, bs = "fs", m = 2) +
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model 4- Global Age pattern with group-level Region smoothers with different wiggliness
Lm_Ba_M1d <- bam(sBa_Ca ~ 
                   s(Age, k = 24, bs = "tp", m = 2) +
                   s(Age, by = Region, k = 24, bs = "tp", m = 2) +
                   s(Region, bs = 're') + 
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model 5- Group-level Age smoothers for each Region with identical wiggliness
Lm_Ba_M1e <- bam(sBa_Ca ~ 
                   s(Age, Region, k = 24, bs = "fs", m = 2) +
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model 6- Group-level Age smoothers for each Region with identical wiggliness
Lm_Ba_M1f <- bam(sBa_Ca ~ 
                   s(Age, by = Region, k = 24, bs = "tp", m = 2) +
                   s(Region, bs = 're') + 
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model selection using Akaike's Information Criterion ####
bbmle::AICtab( Lm_Ba_M1a, Lm_Ba_M1b, Lm_Ba_M1c, Lm_Ba_M1d, Lm_Ba_M1e, Lm_Ba_M1f, 
               base=T,logLik=T,weights=T ) 
summary(Lm_Ba_M1e)
## Optimal model is M1e. 
## The variation in E:Ca is best explained by fully penalised group-level Region smooth functions of Age against a null of no effect. 

#################################################################################################################################
#################################################################################################################################

## Assess for temporal patterns in E:Ca 
## Fit a discrete Year slope and Region random intercept. 
Lm_Ba_M2 <- bam(sBa_Ca ~ 
                  s(Age, Region, k = 24, bs = "fs", m = 2) +
                  s(Region, bs = 're') + 
                  s(Region, Year, bs = 're') + 
                  s(FishID, bs = 're') + 
                  s(FishID, Age, bs = 're'), 
                select = TRUE, 
                method = "fREML", 
                family = Gamma(link = "log"),
                data = Lm_data) 

bbmle::AICtab( Lm_Ba_M1e, Lm_Ba_M2, 
               base=T,logLik=T,weights=T ) 
## Model selection indicates that the null model provides a better fit. 
## Complex models are penalised by AIC model selection criteria. 

## Construct the optimal model for Barium for Lutjanus malabaricus ####
Lm_Ba <- Lm_Ba_M2

## Check the numerical output of the data
summary(Lm_Ba)
anova(Lm_Ba)
plot(fitted(Lm_Ba), residuals(Lm_Ba)) # check residual plots 

## Check diagnostic information for our baseline model
# According to ?gam.check(), the p-value is computed by simulation: the residuals are randomly re-shuffled `k.rep` times to obtain the null distribution of the differencing variance estimator, if there is no pattern in the residuals. 
set.seed(100) ; k.check(Lm_Ba)
set.seed(100) ; gam.check(Lm_Ba) 

## Generate preliminary plots
plot.gam(Lm_Ba, all.terms = T)

## Plot conditional effects of Age for each Region, while setting the individual-specific random effects or offsets to zero. 
## This contrasts with marginal effect plots, where one is interested in the Age effect for each Region, while integrating or holding the random effect values to a constant. 
p <- marginaleffects::plot_predictions(Lm_Ba, condition = c("Age","Region"), type = 'response') + theme_bw()
p

### Simulate Residuals using the DHARMa package 
Lm_Ba_Residuals <- simulateResiduals(fittedModel = Lm_Ba, plot = FALSE)

## Plot residuals against Age data 
plotResiduals(Lm_Ba_Residuals, form = Lm_data$Age, smoothScatter = TRUE) 

#################################################################################################################################
#################################################################################################################################

### Preliminary plot using the `itsadug` package
par(mfrow = c(1, 1))
itsadug::plot_smooth(Lm_Ba, view = "Age", plot_all = c("Region"), rm.ranef = FALSE) # without extracting random effects
itsadug::plot_smooth(Lm_Ba, view = "Age", plot_all = c("Region"), transform = exp, rm.ranef = FALSE) # transform to original scale
# Note that this is equivalent to plot_predictions(..., type = 'response'). 

### Model visualisation for the optimal GAMM model ####
## Plot changes in Ba across Age for each Region for the average Year. 
Lm_Ba_Age_Data <- ddply(Lm_data, .(Region), summarize, 
                        Age = seq(min(Age), max(Age), length = 1000)) 
Lm_Ba_Age_Data$FishID <- "NA"
Lm_Ba_Age_Data$FishID <- factor(Lm_Ba_Age_Data$FishID)
Lm_Ba_Age_Data$Year <- mean(Lm_data$Year)

## Predict how Ba varies with Age assuming a constant Year and Cohort effect.
sapply(Lm_Ba$smooth,"[[","label")
Lm_Ba_PredictAge <- predict.bam(Lm_Ba, 
                                newdata = Lm_Ba_Age_Data, 
                                exclude = c("s(Region)",
                                            "s(Region,Year)",
                                            "s(FishID)",
                                            "s(FishID,Age)"),
                                se.fit = TRUE)

Lm_Ba_Age_Data$mu <- exp(Lm_Ba_PredictAge$fit)
Lm_Ba_Age_Data$ul <- exp(Lm_Ba_PredictAge$fit + 1.96 * Lm_Ba_PredictAge$se.fit )
Lm_Ba_Age_Data$ll <- exp(Lm_Ba_PredictAge$fit - 1.96 * Lm_Ba_PredictAge$se.fit ) 

Lm_Ba_Age_Data <- Lm_Ba_Age_Data %>% filter(Age <= 12)
head(Lm_Ba_Age_Data, 5)

### Plot marginal effect of Age:Region on Ba:Ca, setting all random effect values to a constant value. #### 
plot2 <- Lm_Ba_Age_Data %>%
  ggplot(aes(x = Age, y = mu)) + 
  geom_line(aes(colour = Region)) + 
  geom_ribbon(aes(ymax = ul, ymin = ll, fill = Region), alpha = 0.5) + 
  scale_color_manual(values=c("#E69F00","#0072B2")) + 
  scale_fill_manual(values=c("#E69F00","#0072B2")) +
  labs(y = 'Predicted Ba:Ca (µmol/mol)') + 
  scale_x_continuous(limits = c(0, 12), breaks = c(0,3,6,9,12), labels = c("0","3","6","9","12")) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0.0, 2.5, 5.0, 7.5, 10.0), labels = c("0.0", "2.5", "5.0", "7.5", "10.0")) +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14),
        legend.background = element_rect(fill='white',linewidth=0.5,linetype="solid",colour ="grey15"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1.0, 1.0, 1.0, 1.0, "cm"))
plot2

ggsave("plot2.png", plot = plot2, dpi = 300, width = 297, height = 210, units = "mm")

#################################################################################################################################
#################################################################################################################################

### Magnesium ####
## Construct different models to explore the relationship between Magnesium with Age

## Model 1- Global Mg pattern 
Lm_Mg_M1a <- bam(sMg_Ca ~ 
                   s(Age, k = 24, bs = "tp", m = 2) +
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model 2- Global Mg pattern with Region random effect
Lm_Mg_M1b <- bam(sMg_Ca ~ 
                   s(Age, k = 24, bs = "tp", m = 2) +
                   s(Region, bs = 're') + 
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model 3- Global Age pattern with group-level Region smoothers with identical wiggliness
Lm_Mg_M1c <- bam(sMg_Ca ~ 
                   s(Age, k = 24, bs = "tp", m = 2) +
                   s(Age, Region, k = 24, bs = "fs", m = 2) +
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model 4- Global Age pattern with group-level Region smoothers with different wiggliness
Lm_Mg_M1d <- bam(sMg_Ca ~ 
                   s(Age, k = 24, bs = "tp", m = 2) +
                   s(Age, by = Region, k = 24, bs = "tp", m = 2) +
                   s(Region, bs = 're') + 
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model 5- Group-level Age smoothers for each Region with identical wiggliness
Lm_Mg_M1e <- bam(sMg_Ca ~ 
                   s(Age, Region, k = 24, bs = "fs", m = 2) +
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model 6- Group-level Age smoothers for each Region with identical wiggliness
Lm_Mg_M1f <- bam(sMg_Ca ~ 
                   s(Age, by = Region, k = 24, bs = "tp", m = 2) +
                   s(Region, bs = 're') + 
                   s(FishID, bs = 're') + 
                   s(FishID, Age, bs = 're'), 
                 select = TRUE, 
                 method = "fREML", 
                 family = Gamma(link = "log"),
                 data = Lm_data) 

## Model selection using Akaike's Information Criterion ####
bbmle::AICtab( Lm_Mg_M1a, Lm_Mg_M1b, Lm_Mg_M1c, Lm_Mg_M1d, Lm_Mg_M1e, Lm_Mg_M1f, 
               base=T,logLik=T,weights=T ) 
summary(Lm_Mg_M1e)
## Optimal model is M1e. 
## The variation in E:Ca is best explained by fully penalised group-level Region smooth functions of Age against a null of no effect. 

#################################################################################################################################
#################################################################################################################################

## Assess for temporal patterns in E:Ca 
## Fit a discrete Year slope and Region random intercept. 
Lm_Mg_M2 <- bam(sMg_Ca ~ 
                  s(Age, Region, k = 24, bs = "fs", m = 2) +
                  s(Region, bs = 're') + 
                  s(Region, Year, bs = 're') + 
                  s(FishID, bs = 're') + 
                  s(FishID, Age, bs = 're'), 
                select = TRUE, 
                method = "fREML", 
                family = Gamma(link = "log"),
                data = Lm_data) 

bbmle::AICtab( Lm_Mg_M1e, Lm_Mg_M2, base=T,logLik=T,weights=T ) 
## Model selection indicates that the null model provides a better fit. 
## Complex models are penalised by AIC model selection criteria. 
Lm_Mg <- Lm_Mg_M2

## Check the numerical output of the data
summary(Lm_Mg)
anova(Lm_Mg)
plot(fitted(Lm_Mg), residuals(Lm_Mg)) # check residual plots 

## Check diagnostic information for our baseline model
# According to ?gam.check(), the p-value is computed by simulation: the residuals are randomly re-shuffled `k.rep` times to obtain the null distribution of the differencing variance estimator, if there is no pattern in the residuals. 
set.seed(100) ; k.check(Lm_Mg)
set.seed(100) ; gam.check(Lm_Mg) 

## Generate preliminary plots
plot.gam(Lm_Mg, all.terms = T)

## Plot conditional effects of Age for each Region, while setting the individual-specific random effects or offsets to zero. 
## This contrasts with marginal effect plots, where one is interested in the Age effect for each Region, while integrating or holding the random effect values to a constant. 
p <- marginaleffects::plot_predictions(Lm_Mg, condition = c("Age","Region"), type = 'response') + theme_bw()
p

### Simulate Residuals using the DHARMa package 
Lm_Mg_Residuals <- simulateResiduals(fittedModel = Lm_Mg, plot = FALSE)

## Plot residuals against Age data 
plotResiduals(Lm_Mg_Residuals, form = Lm_data$Age, smoothScatter = TRUE) 

#################################################################################################################################
#################################################################################################################################

### Preliminary plot using the `itsadug` package
par(mfrow = c(1, 1))
itsadug::plot_smooth(Lm_Mg, view = "Age", plot_all = c("Region"), rm.ranef = FALSE) # without extracting random effects
itsadug::plot_smooth(Lm_Mg, view = "Age", plot_all = c("Region"), transform = exp, rm.ranef = FALSE) # transform to original scale
# Note that this is equivalent to plot_predictions(..., type = 'response'). 

### Model visualisation for the optimal GAMM model ####
## Plot changes in Mg across Age for each Region for the average Year. 
Lm_Mg_Age_Data <- ddply(Lm_data, .(Region), summarize, 
                        Age = seq(min(Age), max(Age), length = 1000)) 
Lm_Mg_Age_Data$FishID <- "NA"
Lm_Mg_Age_Data$FishID <- factor(Lm_Mg_Age_Data$FishID)
Lm_Mg_Age_Data$Year <- mean(Lm_data$Year)

## Predict how Mg varies with Age assuming a constant Year and Cohort effect.
sapply(Lm_Mg$smooth,"[[","label")
Lm_Mg_PredictAge <- predict.bam(Lm_Mg, 
                                newdata = Lm_Mg_Age_Data, 
                                exclude = c("s(Region)",
                                            "s(Region,Year)",
                                            "s(FishID)",
                                            "s(FishID,Age)"),
                                se.fit = TRUE)

Lm_Mg_Age_Data$mu <- exp(Lm_Mg_PredictAge$fit)
Lm_Mg_Age_Data$ul <- exp(Lm_Mg_PredictAge$fit + 1.96 * Lm_Mg_PredictAge$se.fit )
Lm_Mg_Age_Data$ll <- exp(Lm_Mg_PredictAge$fit - 1.96 * Lm_Mg_PredictAge$se.fit ) 

Lm_Mg_Age_Data <- Lm_Mg_Age_Data %>% filter(Age <= 12)
head(Lm_Mg_Age_Data, 5)

### Plot marginal effect of Age:Region on Mg:Ca, setting all random effect values to a constant value. #### 
plot3 <- Lm_Mg_Age_Data %>%
  ggplot(aes(x = Age, y = mu)) + 
  geom_line(aes(colour = Region)) + 
  geom_ribbon(aes(ymax = ul, ymin = ll, fill = Region), alpha = 0.5) + 
  scale_color_manual(values=c("#E69F00","#0072B2")) + 
  scale_fill_manual(values=c("#E69F00","#0072B2")) +
  labs(y = 'Predicted Mg:Ca (µmol/mol)') + 
  scale_x_continuous(limits = c(0, 12), breaks = c(0,3,6,9,12), labels = c("0","3","6","9","12")) +
  scale_y_continuous(limits = c(0, 400), breaks = c(0, 100, 200, 300, 400), labels = c("0", "100", "200", "300", "400")) +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14),
        legend.background = element_rect(fill='white',linewidth=0.5,linetype="solid",colour ="grey15"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1.0, 1.0, 1.0, 1.0, "cm"))
plot3

ggsave("plot3.png", plot = plot3, dpi = 300, width = 297, height = 210, units = "mm")

#################################################################################################################################
#################################################################################################################################

## Plot geographic map ####
register_stadiamaps(key = "4eef9e3f-e6c1-4911-bdb3-8c243f7ef226")
p <- get_stadiamap(bbox = c(top = 12, bottom = -15, left = 95, right = 140), maptype = c("stamen_terrain_background"), crop = FALSE, zoom = 7) 
ggmap(p) 

## Base map
base_map <- ggmap(p) + 
  ylab("Latitude") + xlab("Longitude") + 
  geom_circle(aes(x0 = 106.8, y0 = 1.3, r = 1.4), colour = 'black', fill = 'firebrick3', alpha = 0.2) +
  geom_circle(aes(x0 = 130.0, y0 = 0.4, r = 1.4), colour = 'black', fill = 'firebrick3', alpha = 0.2) +
  
  annotate('text', x = 106.8, y = 1.3, color = 'white', size = 4, label = 'RA', fontface = 2) +
  annotate('text', x = 130.0, y = 0.4, color = 'white', size = 4, label = 'SR', fontface = 2) +
  
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        plot.margin = grid::unit(c(1.0, 1.0, 1.0, 1.0), "cm"))
base_map

## World map
world <- ne_countries(scale = "medium", returnclass = "sf")
world_map <- ggplot(world)+ 
  geom_sf() + 
  geom_rect(aes(xmin = 95, xmax = 140, ymin = -7, ymax = 9), colour = 'red', fill = NA, linewidth = 0.7) + 
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(50, 200), ylim = c(-40, 40), expand = FALSE) +
  theme(panel.border = element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.6),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = grid::unit(c(0, 0, -1, -1), "mm"))
world_map

## Ovelay maps
plot4 <- ggdraw() +
  draw_plot(base_map) +
  draw_plot(world_map, height = 0.180, x = -0.283, y = 0.137) 
plot4

ggsave("plot4.png", plot = plot4, dpi = 300, width = 297, height = 210, units = "mm")

#################################################################################################################################
#################################################################################################################################

## Filter for samples only, ignoring the NIST and MACS standards 
## Check the data for potential patterns
data %>%
  filter(grepl('AA|AB', FishID)) %>%
  ggplot(aes(x = distance, y = Mn_Ca, colour = FishID, group = FishID)) + 
  geom_line(alpha = 0.3) + 
  geom_point(alpha = 0.3) + 
  theme(legend.position = "none")

## There is an outlier data point present, which is likely made up of vaterite instead of aragonite. 
data %>%
  filter(grepl('AA|AB', FishID)) %>%
  ggplot(aes(x = distance, y = Mg_Ca, colour = FishID, group = FishID)) + 
  geom_line(alpha = 0.3) + 
  geom_point(alpha = 0.3) + 
  theme(legend.position = "none")

data %>%
  #filter(Mg_Ca < 1000) %>%
  filter(grepl('AA|AB', FishID)) %>%
  ggplot(aes(x = distance, y = Sr_Ca, colour = FishID, group = FishID)) + 
  geom_line(alpha = 0.3) + 
  geom_point(alpha = 0.3) + 
  theme(legend.position = "none")

data %>%
  filter(grepl('AA|AB', FishID)) %>%
  arrange(desc(Mg_Ca))

#################################################################################################################################
#################################################################################################################################

#save.image("otolith_chemistry_v2.RData")
#save.image("2025_otolith_chemistry_v2.RData")

#################################################################################################################################
#################################################################################################################################
