rm(list = ls())

library(tidyverse)
library(stringr)
library(lubridate)
library(bigleaf)

constants <- bigleaf.constants()

getwd()
setwd("C:/Users/yekim/OneDrive - Desert Research Institute/Documents/research/phase1_interpolation")


### meta data
meta <- read.csv("./flux_ET_dataset/station_metadata.csv",skip = 1)
colnames(meta)[1] <- "site"
head(meta)
meta_use <- meta %>% select(site,Latitude,Longitude,Elevation..m.)


#####################
## 1. gridmet data
#####################
getwd()
setwd("./gridMET")
ls <- list.files()

for(i in c(1:length(ls))){
	file_name <- ls[i]
	df <- read.csv(file_name)
	
	# Date as r Date
	colnames(df)[1] <- "DATE"
	df$DATE <- as.Date(df$DATE)
	
	# add site column
	df$site <- str_remove(ls[i],".csv")
	
	#merge all sites
	if(i==1){
		gridmet <- df %>% select(site,DATE,gridMET_q,gridMET_srad,gridMET_u10,gridMET_tmax,gridMET_tmin)
	} else {
		gridmet <- full_join(gridmet,df %>% select(site,DATE,gridMET_q,gridMET_srad,gridMET_u10,gridMET_tmax,gridMET_tmin))
	}
}

setwd("..")

#####################
## 2. ra
#####################
df <- gridmet %>% left_join(meta_use)

#Ra
calc_Ra <- function(date, lat, lon) {
	# Ensure the date is in Date format
	date <- as.Date(date)
	
	# Calculate the day of year (J)
	J <- as.numeric(format(date, "%j"))
	
	# Convert latitude from degrees to radians
	phi <- lat * pi / 180
	
	# Solar constant in MJ m-2 min-1
	Gsc <- 0.0820
	
	# Relative distance between earth and sun (dr)
	dr <- 1 + 0.033 * cos(2 * pi / 365 * J)
	
	# Solar declination (delta) in radians
	delta <- 0.409 * sin(2 * pi / 365 * J - 1.39)
	
	# Sunset hour angle (omega_s) in radians
	# Note: Ensure that the argument to acos is within [-1, 1] to avoid NaNs.
	omega_s <- acos(-tan(phi) * tan(delta))
	
	# Calculate extraterrestrial radiation (Ra) in MJ m-2 d-1
	Ra <- (24 * 60 / pi) * Gsc * dr * 
		(omega_s * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(omega_s))
	
	# MJ m-2 d-1 to W m-2
	Ra_Wm2 <- Ra/(60 * 60 * 24) * (10^6)
	
	# Return results as a data frame
	return(Ra_Wm2)
}

## apply
df <- df %>%
	mutate(Ra_Wm2 = calc_Ra(DATE,Latitude,Longitude))

#####################
## 3. SFE 
#####################

# air pressure and gamma from elevation
PA <- meta_use %>%
	mutate(PA = 101.3 * ((293 - 0.0065*Elevation..m.)/293)^5.26,
				 gamma = 0.000665 * PA
	) %>%
	select(site,PA,gamma)

df <- df %>% left_join(PA)

# calculate saturation humidity variables
df <- df %>%
	mutate(gridMET_tmax_K = gridMET_tmax,
				 gridMET_tmin_K = gridMET_tmin,
				 gridMET_tmax = gridMET_tmax_K - 273.16,
				 gridMET_tmin = gridMET_tmin_K - 273.16,
				 gridMET_tavg = (gridMET_tmax +gridMET_tmin)/2,
				 
				 gridMET_Delta = 2503 * exp(17.27*gridMET_tavg/(gridMET_tavg + 237.3))/(gridMET_tavg + 237.3)^2,
				 
				 gridMET_esat_tmax = 0.6108 * exp(17.27*gridMET_tmax/(gridMET_tmax + 237.3)),
				 gridMET_esat_tmin = 0.6108 * exp(17.27*gridMET_tmin/(gridMET_tmin + 237.3)),
				 gridMET_esat_tavg = 0.6108 * exp(17.27*gridMET_tavg/(gridMET_tavg + 237.3)),
				 gridMET_es = (gridMET_esat_tmax + gridMET_esat_tmin)/2,
				 
				 density = air.density(gridMET_tavg,PA)
				 
	)

# calculate humidity
df <- df %>%
	mutate(gridMET_ea = gridMET_q * PA/((1 - constants$eps) * gridMET_q + constants$eps),
				 
				 gridMET_RH= gridMET_ea/gridMET_es,
				 gridMET_RH= ifelse(gridMET_RH > 0.99, 0.99, gridMET_RH),
				 gridMET_ea = gridMET_RH*gridMET_es,
	)

# caclulate u2
df <- df %>%
	mutate(gridMET_u2 = gridMET_u10 * 4.87 /log(67.8*10 - 5.42),
				 gridMET_ga_H = gridMET_u2/208,
				 gridMET_ra_H = 1/gridMET_ga_H
	)

# radiation 
df <- df %>% 
	mutate(rso = (0.75+2e-5 * Elevation..m.)*Ra_Wm2,
				 gridMET_srad_corr = ifelse(gridMET_srad < 20, 20, gridMET_srad),
				 Rs_ratio = ifelse(gridMET_srad /rso > 1, 1, ifelse(gridMET_srad /rso < 0.3, 0.3, gridMET_srad/rso)),
				 gridMET_Rns = gridMET_srad_corr * (1 - 0.23), #ASCE standard reference ET 0.23
				 
				 fcd = 1.35 * Rs_ratio - 0.35,
				 gridMET_Rnl = fcd * 5.672454e-08 * (0.34 - 0.14 * sqrt(gridMET_ea)) * 
				 	(gridMET_tmax_K^4 + gridMET_tmin_K^4)/2,
				 
				 gridMET_Rn = gridMET_Rns - gridMET_Rnl,
				 gridMET_Rn = ifelse(gridMET_Rn < 5, 5, gridMET_Rn),
				 
				 gridMET_Rn_MJ_day = gridMET_Rn * 60 * 60 * 24 / (10^6),
				 
				 # SFE
				 SFE = (0.408*gridMET_RH*gridMET_Delta*gridMET_Rn_MJ_day)/
				 	(gridMET_RH*gridMET_Delta+gamma),
				 EA_SFE = (gamma*2.6*(1+0.54*gridMET_u2)*(gridMET_es - gridMET_ea))/
				 	(gridMET_RH*gridMET_Delta+gamma)/2,
	)

df_use <- df %>% select(site, DATE, SFE,EA_SFE)

#####################
## 4. Merge data 
#####################
setwd("C:/Users/yekim/OneDrive - Desert Research Institute/Documents/research/EToF_interpolation_uncertainty_analysis")

Flux_OpenET_gridMET <- read.csv("daily_data.csv")
Flux_OpenET_gridMET <- Flux_OpenET_gridMET %>% 
	mutate(DATE=as.Date(DATE))
Flux_OpenET_gridMET <- Flux_OpenET_gridMET %>% left_join(df_use)
write.csv(Flux_OpenET_gridMET,"daily_data.csv",row.names = F)
