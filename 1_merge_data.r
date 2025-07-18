rm(list = ls())

library(tidyverse)
library(stringr)
library(lubridate)

getwd()
setwd("C:/Users/yekim/OneDrive - Desert Research Institute/Documents/research/phase1_interpolation")

#####################
## 1. OpenET
#####################
OpenET <- read.csv("./OpenET_PhaseII_model_ET_dataset/daily_data.csv")

# site name
colnames(OpenET)[1] <- "site"
site <- unique(OpenET$site)

# Date as r date
OpenET$DATE <- as.Date(OpenET$DATE,format = "%m/%d/%Y")

#####################
## 2. Flux data
#####################
getwd()
setwd("./flux_ET_dataset/daily_data_files")
getwd()

for(i in c(1:length(site))){
	file_name <- paste0(site[i],"_daily_data.csv")
	df <- read.csv(file_name)
	
	# Date as r Date
	colnames(df)[1] <- "DATE"
	df$DATE <- as.Date(df$DATE)
	
	# select variables
	if(sum(colnames(df) == 'ET') == 0){df$ET <- NA}
	if(sum(colnames(df) == 'ET_corr') == 0){df$ET_corr <- NA}
	if(sum(colnames(df) == 'ebr') == 0){df$ebr <- NA}
	if(sum(colnames(df) == 'LE_subday_gaps') == 0){df$LE_subday_gaps <- NA}
	
	# add site column
	df$site <- site[i]
	
	df <- df %>% select(site,DATE, ET, ET_corr, ebr, LE_subday_gaps)
	df <- df %>% filter(!is.na(ET_corr))
	
	#merge OpenET data
	merged <- left_join(df,OpenET,by = c("site","DATE"))

	#merge all sites
	if(i==1){
		Flux_OpenET <- merged
	} else {
		Flux_OpenET <- full_join(Flux_OpenET,merged)
	}
}

setwd("..")
setwd("..")

#####################
## 3. gridmet data
#####################
getwd()
setwd("./gridMET")

for(i in c(1:length(site))){
	file_name <- paste0(site[i],".csv")
	df <- read.csv(file_name)
	
	# Date as r Date
	colnames(df)[1] <- "DATE"
	df$DATE <- as.Date(df$DATE)
	
	# add site column
	df$site <- site[i]
	
	#merge all sites
	if(i==1){
		gridmet <- df %>% select(site,DATE,gridMET_ETo)
	} else {
		gridmet <- full_join(gridmet,df %>% select(site,DATE,gridMET_ETo))
	}
}

Flux_OpenET_gridMET <- left_join(Flux_OpenET,gridmet)

setwd("..")

#####################
## 4. ETo bias correction
#####################
ratio <- read.csv("monthly_eto_correction_ratio.csv")

Flux_OpenET_gridMET <- Flux_OpenET_gridMET %>% 
	mutate(YEAR = year(DATE),
				 MON = month(DATE))

Flux_OpenET_gridMET <- Flux_OpenET_gridMET %>% left_join(ratio)
Flux_OpenET_gridMET <- Flux_OpenET_gridMET %>%
	mutate(gridMET_ETo_corr = gridMET_ETo * ratio)

Flux_OpenET_gridMET <- Flux_OpenET_gridMET %>% 
	select(!YEAR) %>% select(!MON) %>% select(!ratio)

head(Flux_OpenET_gridMET)

### save
setwd("C:/Users/yekim/OneDrive - Desert Research Institute/Documents/research/EToF_interpolation_uncertainty_analysis")
write.csv(Flux_OpenET_gridMET,"daily_data.csv",row.names = F)


