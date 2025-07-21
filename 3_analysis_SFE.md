interpolation uncertainty analysis_SFE method
================
Yeonuk Kim

## Data preparation

``` r
rm(list = ls())

library("tidyverse")
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library("lubridate")
library(Metrics)

# interpolation max days
maxdays <- 32

### daily EC, OpenET, gridMET ETo
df <- read.csv("daily_data.csv")
df <- df %>% 
    mutate(DATE=as.Date(DATE),
                 YEAR = year(DATE),
                 MON = month(DATE)) 


# filter data to OpenET data available period
df <- df %>%
    mutate(RS_OBS = ifelse(!is.na(eeMETRIC),T,F),
                 RS_OBS_SIMS = ifelse(!is.na(SIMS),T,F),
                 ET_OBS = ifelse(!is.na(ET_corr),T,F)
    )

vaild_date <- df %>%
    filter(RS_OBS) %>%
    group_by(site) %>%
    summarise(start = min(DATE) - maxdays,
                        end = max(DATE) + maxdays)


df <- df %>% left_join(vaild_date)
```

    ## Joining with `by = join_by(site)`

``` r
df <- df %>% filter(DATE >= start & DATE <= end) %>%
    select(!start) %>% select(!end)


###
meta <- read.csv("station_metadata.csv",skip = 1)
colnames(meta)[1] <- "site"

meta_use <- meta %>% 
    select(site,General.classification,Land.cover.type) %>%
    rename(site=site, 
                 class1 = General.classification,
                 class2 = Land.cover.type) %>%
    mutate(is_crop = ifelse(class1 == "Croplands",TRUE,FALSE),
                 class1 = factor(class1,
                                                levels = c("Croplands","Mixed Forests","Evergreen Forests",
                                                                     "Shrublands","Grasslands","Wetland/Riparian",""))
                 )

df <- df %>% left_join(meta_use)
```

    ## Joining with `by = join_by(site)`

## Interpolation preparation

``` r
# nearest temporal distance
df <- df %>%
    group_by(site) %>%        
    arrange(DATE) %>%   
    mutate(nearest_distance = sapply(DATE, function(d) {
        if (any(RS_OBS)) {
            min(abs(as.numeric(difftime(DATE[RS_OBS], d, units = "days"))))
        } else {
            Inf  # No RS_OBS in the group at all
        }
    }))

# grouping for interpolation (if gap is longer than 32, set as different group)
cutting_long_gap <- function(x, maxgap){
    # Find positions where x equals maxgap
    pos <- which(x == maxgap+1)
    
    # Start with all observations in group 1
    group <- rep(1, length(x))
    
    if(length(pos) >= 1){
        # For the first occurrence, increment the group for observations after it
        if(pos[1] < length(x)){
            group[(pos[1] + 1):length(x)] <- group[(pos[1] + 1):length(x)] + 1
        }
        # For subsequent occurrences, increment the group starting from the occurrence itself
        if(length(pos) >= 2){
            for(i in 2:length(pos)){
                group[pos[i]:length(x)] <- group[pos[i]:length(x)] + 1
            }
        }
    }
    return(group)
}

df <- df %>%
    group_by(site) %>%
    arrange(DATE) %>%
    mutate(gap_group = cutting_long_gap(nearest_distance,maxdays)) %>%
    ungroup()

# EAF
df <- df %>%
    mutate(EAF_geeSEBAL = (geeSEBAL-SFE)/EA_SFE,
                 EAF_PT.JPL = (PT.JPL-SFE)/EA_SFE,
                 EAF_SSEBop = (SSEBop-SFE)/EA_SFE,
                 EAF_SIMS = (SIMS-SFE)/EA_SFE,
                 EAF_eeMETRIC = (eeMETRIC-SFE)/EA_SFE,
                 EAF_DisALEXI = (DisALEXI-SFE)/EA_SFE,
                 EAF_Ensemble = (Ensemble-SFE)/EA_SFE,
                 
                 EAF_geeSEBAL = ifelse(EAF_geeSEBAL > 1, 1,ifelse(EAF_geeSEBAL < -1, -1, EAF_geeSEBAL)),
                 EAF_PT.JPL = ifelse(EAF_PT.JPL > 1, 1,ifelse(EAF_PT.JPL < -1, -1, EAF_PT.JPL)),
                 EAF_SSEBop = ifelse(EAF_SSEBop > 1, 1,ifelse(EAF_SSEBop < -1, -1, EAF_SSEBop)),
                 EAF_SIMS = ifelse(EAF_SIMS > 1, 1,ifelse(EAF_SIMS < -1, -1, EAF_SIMS)),
                 EAF_eeMETRIC = ifelse(EAF_eeMETRIC > 1, 1,ifelse(EAF_eeMETRIC < -1, -1, EAF_eeMETRIC)),
                 EAF_DisALEXI = ifelse(EAF_DisALEXI > 1, 1,ifelse(EAF_DisALEXI < -1, -1, EAF_DisALEXI)),
                 EAF_Ensemble = ifelse(EAF_Ensemble > 1, 1,ifelse(EAF_Ensemble < -1, -1, EAF_Ensemble)),
                 )
```

## Interpolation

``` r
### site filter
selected_site <- df %>% group_by(site) %>%
    summarise(n = sum(RS_OBS)) %>%
    filter(n > 1) %>%
    pull(site)

selected_site_SIMS <- df %>% group_by(site) %>%
    summarise(n = sum(RS_OBS_SIMS)) %>%
    filter(n > 1) %>%
    pull(site)

df <- df %>% filter(site %in% selected_site)
df_SIMS <- df %>% filter(site %in% selected_site_SIMS)


# linear interpolated
# Define a safe interpolation function
safe_approx <- function(x, y, new_x, rule = 2) {
    # Only interpolate if there are at least 2 non-NA y values
    if(sum(!is.na(y)) >= 2) {
        # Use only the non-NA indices for interpolation
        approx(x = x[!is.na(y)], y = y[!is.na(y)], 
                     xout = new_x, method = "linear", rule = rule)$y
    } else if(sum(!is.na(y)) == 1){
        approx(x = x[!is.na(y)], y = y[!is.na(y)], 
                     xout = new_x, method = "constant", rule = rule)$y
    } else {
        rep(NA, length(new_x))
    }
}

# interpolation
df <- df %>%
    group_by(site,gap_group) %>%        
    arrange(DATE) %>%   
    mutate(interpolated_EAF_geeSEBAL = safe_approx(seq_along(EAF_geeSEBAL), EAF_geeSEBAL, seq_along(EAF_geeSEBAL)),
                 interpolated_EAF_PT.JPL = safe_approx(seq_along(EAF_PT.JPL), EAF_PT.JPL, seq_along(EAF_PT.JPL)),
                 interpolated_EAF_SSEBop = safe_approx(seq_along(EAF_SSEBop), EAF_SSEBop, seq_along(EAF_SSEBop)),
                 interpolated_EAF_eeMETRIC = safe_approx(seq_along(EAF_eeMETRIC), EAF_eeMETRIC, seq_along(EAF_eeMETRIC)),
                 interpolated_EAF_DisALEXI = safe_approx(seq_along(EAF_DisALEXI), EAF_DisALEXI, seq_along(EAF_DisALEXI)),
                 interpolated_EAF_Ensemble = safe_approx(seq_along(EAF_Ensemble), EAF_Ensemble, seq_along(EAF_Ensemble)),
                 )

df_SIMS <- df_SIMS %>%
    group_by(site,gap_group) %>%        
    arrange(DATE) %>%   
    mutate(interpolated_EAF_SIMS = safe_approx(seq_along(EAF_SIMS), EAF_SIMS, seq_along(EAF_SIMS)))

df <- df %>% left_join(df_SIMS %>% select(site,DATE,interpolated_EAF_SIMS))
```

    ## Adding missing grouping variables: `gap_group`
    ## Joining with `by = join_by(site, DATE, gap_group)`

``` r
# Gap-filled (GF) ET for each model
df <- df %>%
    mutate(GF_geeSEBAL = interpolated_EAF_geeSEBAL * EA_SFE + SFE,
                 GF_PT.JPL = interpolated_EAF_PT.JPL * EA_SFE + SFE,
                 GF_SSEBop = interpolated_EAF_SSEBop * EA_SFE + SFE,
                 GF_eeMETRIC = interpolated_EAF_eeMETRIC * EA_SFE + SFE,
                 GF_DisALEXI = interpolated_EAF_DisALEXI * EA_SFE + SFE,
                 GF_SIMS = interpolated_EAF_SIMS * EA_SFE + SFE,
                 GF_Ensemble = interpolated_EAF_Ensemble * EA_SFE + SFE,
                 
                 
                 GF_geeSEBAL = ifelse(GF_geeSEBAL < 0, 0,GF_geeSEBAL),
                 GF_PT.JPL = ifelse(GF_PT.JPL < 0, 0,GF_PT.JPL),
                 GF_SSEBop = ifelse(GF_SSEBop < 0, 0,GF_SSEBop),
                 GF_SIMS = ifelse(GF_SIMS < 0, 0,GF_SIMS),
                 GF_eeMETRIC = ifelse(GF_eeMETRIC < 0, 0,GF_eeMETRIC),
                 GF_DisALEXI = ifelse(GF_DisALEXI  < 0, 0,GF_DisALEXI),
                 GF_Ensemble = ifelse(GF_Ensemble < 0, 0,GF_Ensemble),
                 
                 )
```

## number of data

``` r
for(i in c(0:32)){

    rmse_df <- df %>%
        ungroup() %>%
        group_by(class1) %>%
        filter(nearest_distance == i & ET_OBS == T & !is.na(GF_geeSEBAL)) %>%
        summarise(ETm = mean(ET_corr,na.rm=T),
                            distance = i,
                            n = n()
        ) 
    
    if(i == 0){
        result <- rmse_df
    } else{
        result <- rbind(result,rmse_df)
    }
}

result %>%
    ggplot(aes(distance,n)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(
        y = "number of data", 
        x = "Temporal distance from satellite (days)"
    ) +
    scale_y_log10()+
    facet_wrap(~class1, scales = "free_y")
```

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
result %>%
    ggplot(aes(distance, ETm)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(
        y = "mean ET (mm/day)", 
        x = "Temporal distance from satellite (days)"
    ) +
    facet_wrap(~class1, scales = "free_y")
```

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

## correlation

``` r
for(i in c(0:32)){

    rmse_df <- df %>%
        ungroup() %>%
        group_by(class1) %>%
        filter(nearest_distance == i & ET_OBS == T & !is.na(GF_geeSEBAL)) %>%
        summarise(ETm = mean(ET_corr,na.rm=T),
                            
                            geeSEBAL = cor(ET_corr,GF_geeSEBAL),
                            PT.JPL = cor(ET_corr,GF_PT.JPL),
                            SSEBop = cor(ET_corr,GF_SSEBop),
                            eeMETRIC = cor(ET_corr,GF_eeMETRIC),
                            DisALEXI = cor(ET_corr,GF_DisALEXI),
                            SIMS = cor(ET_corr,GF_SIMS),
                            Ensemble = cor(ET_corr,GF_Ensemble),
                            
                            distance = i,
                            n = n()
        ) 
    
    if(i == 0){
        result <- rmse_df
    } else{
        result <- rbind(result,rmse_df)
    }
}

result %>%
    select(!n) %>%
    select(!ETm) %>%
    gather(key = "model",value = "cor", -c(class1,distance)) %>%
    mutate(model = factor(model,
                                                levels = c("Ensemble","geeSEBAL","PT.JPL","SSEBop",
                                                                     "eeMETRIC","DisALEXI","SIMS"))
                 ) %>%
    ggplot(aes(distance,cor,color = model, shape=model)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(
        y = "Pearson correlation", 
        x = "Temporal distance from satellite (days)"
    ) +
    facet_wrap(~class1, scales = "free_y")
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values because more
    ## than 6 becomes difficult to discriminate
    ## ℹ you have requested 7 values. Consider specifying shapes manually if you need
    ##   that many have them.

    ## Warning: Removed 165 rows containing missing values or values outside the scale range
    ## (`geom_line()`).

    ## Warning: Removed 198 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
result %>%
    select(!n) %>%
    select(!ETm) %>%
    gather(key = "model",value = "cor", -c(class1,distance)) %>%
    mutate(model = factor(model,
                                                levels = c("Ensemble","geeSEBAL","PT.JPL","SSEBop",
                                                                     "eeMETRIC","DisALEXI","SIMS"))
                 ) %>%
    ggplot(aes(distance,cor,color = model, shape=model)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(
        y = "Pearson correlation", 
        x = "Temporal distance from satellite (days)"
    ) +
    facet_wrap(~class1)
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values because more
    ## than 6 becomes difficult to discriminate
    ## ℹ you have requested 7 values. Consider specifying shapes manually if you need
    ##   that many have them.

    ## Warning: Removed 165 rows containing missing values or values outside the scale range
    ## (`geom_line()`).

    ## Warning: Removed 198 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
result %>%
    filter(class1 == "Croplands") %>%
    select(!n) %>%
    select(!ETm) %>%
    gather(key = "model",value = "cor", -c(class1,distance)) %>%
    mutate(model = factor(model,
                                                levels = c("Ensemble","geeSEBAL","PT.JPL","SSEBop",
                                                                     "eeMETRIC","DisALEXI","SIMS"))
                 ) %>%
    ggplot(aes(distance,cor,color = model, shape=model)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(title = "crop",
        y = "Pearson correlation", 
        x = "Temporal distance from satellite (days)"
    )
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values because more
    ## than 6 becomes difficult to discriminate
    ## ℹ you have requested 7 values. Consider specifying shapes manually if you need
    ##   that many have them.

    ## Warning: Removed 33 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

## Normalized RMSE

``` r
for(i in c(0:32)){

    rmse_df <- df %>%
        ungroup() %>%
        group_by(class1) %>%
        filter(nearest_distance == i & ET_OBS == T & !is.na(GF_geeSEBAL)) %>%
        summarise(ETm = mean(ET_corr,na.rm=T),
                            
                            geeSEBAL = rmse(ET_corr,GF_geeSEBAL)/ETm,
                            PT.JPL = rmse(ET_corr,GF_PT.JPL)/ETm,
                            SSEBop = rmse(ET_corr,GF_SSEBop)/ETm,
                            eeMETRIC = rmse(ET_corr,GF_eeMETRIC)/ETm,
                            DisALEXI = rmse(ET_corr,GF_DisALEXI)/ETm,
                            SIMS = rmse(ET_corr,GF_SIMS)/ETm,
                            Ensemble = rmse(ET_corr,GF_Ensemble)/ETm,
                            
                            distance = i,
                            n = n()
        ) 
    
    if(i == 0){
        result <- rmse_df
    } else{
        result <- rbind(result,rmse_df)
    }
}

result %>%
    select(!n) %>%
    select(!ETm) %>%
    gather(key = "model",value = "nRMSE", -c(class1,distance)) %>%
    mutate(model = factor(model,
                                                levels = c("Ensemble","geeSEBAL","PT.JPL","SSEBop",
                                                                     "eeMETRIC","DisALEXI","SIMS"))
                 ) %>%
    ggplot(aes(distance,nRMSE*100,color = model, shape=model)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(
        y = "Normalized RMSE (%)", 
        x = "Temporal distance from satellite (days)"
    ) +
    facet_wrap(~class1, scales = "free_y")
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values because more
    ## than 6 becomes difficult to discriminate
    ## ℹ you have requested 7 values. Consider specifying shapes manually if you need
    ##   that many have them.

    ## Warning: Removed 165 rows containing missing values or values outside the scale range
    ## (`geom_line()`).

    ## Warning: Removed 198 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
result %>%
    select(!n) %>%
    select(!ETm) %>%
    gather(key = "model",value = "nRMSE", -c(class1,distance)) %>%
    mutate(model = factor(model,
                                                levels = c("Ensemble","geeSEBAL","PT.JPL","SSEBop",
                                                                     "eeMETRIC","DisALEXI","SIMS"))
                 ) %>%
    ggplot(aes(distance,nRMSE*100,color = model, shape=model)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(
        y = "Normalized RMSE (%)", 
        x = "Temporal distance from satellite (days)"
    ) +
    facet_wrap(~class1)
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values because more
    ## than 6 becomes difficult to discriminate
    ## ℹ you have requested 7 values. Consider specifying shapes manually if you need
    ##   that many have them.

    ## Warning: Removed 165 rows containing missing values or values outside the scale range
    ## (`geom_line()`).

    ## Warning: Removed 198 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
result %>%
    filter(class1 == "Croplands") %>%
    select(!n) %>%
    select(!ETm) %>%
    gather(key = "model",value = "nRMSE", -c(class1,distance)) %>%
    mutate(model = factor(model,
                                                levels = c("Ensemble","geeSEBAL","PT.JPL","SSEBop",
                                                                     "eeMETRIC","DisALEXI","SIMS"))
                 ) %>%
    ggplot(aes(distance,nRMSE*100,color = model, shape=model)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(title = "crop",
        y = "Normalized RMSE (%)", 
        x = "Temporal distance from satellite (days)"
    )
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values because more
    ## than 6 becomes difficult to discriminate
    ## ℹ you have requested 7 values. Consider specifying shapes manually if you need
    ##   that many have them.

    ## Warning: Removed 33 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

## Normalized MAE

``` r
for(i in c(0:32)){

    rmse_df <- df %>%
        ungroup() %>%
        group_by(class1) %>%
        filter(nearest_distance == i & ET_OBS == T & !is.na(GF_geeSEBAL)) %>%
        summarise(ETm = mean(ET_corr,na.rm=T),
                            
                            geeSEBAL = mae(ET_corr,GF_geeSEBAL)/ETm,
                            PT.JPL = mae(ET_corr,GF_PT.JPL)/ETm,
                            SSEBop = mae(ET_corr,GF_SSEBop)/ETm,
                            eeMETRIC = mae(ET_corr,GF_eeMETRIC)/ETm,
                            DisALEXI = mae(ET_corr,GF_DisALEXI)/ETm,
                            SIMS = mae(ET_corr,GF_SIMS)/ETm,
                            Ensemble = mae(ET_corr,GF_Ensemble)/ETm,
                            
                            distance = i,
                            n = n()
        ) 
    
    if(i == 0){
        result <- rmse_df
    } else{
        result <- rbind(result,rmse_df)
    }
}

result %>%
    select(!n) %>%
    select(!ETm) %>%
    gather(key = "model",value = "nMAE", -c(class1,distance)) %>%
    mutate(model = factor(model,
                                                levels = c("Ensemble","geeSEBAL","PT.JPL","SSEBop",
                                                                     "eeMETRIC","DisALEXI","SIMS"))
                 ) %>%
    ggplot(aes(distance,nMAE*100,color = model, shape=model)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(
        y = "Normalized MAE (%)", 
        x = "Temporal distance from satellite (days)"
    ) +
    facet_wrap(~class1, scales = "free_y")
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values because more
    ## than 6 becomes difficult to discriminate
    ## ℹ you have requested 7 values. Consider specifying shapes manually if you need
    ##   that many have them.

    ## Warning: Removed 165 rows containing missing values or values outside the scale range
    ## (`geom_line()`).

    ## Warning: Removed 198 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
result %>%
    select(!n) %>%
    select(!ETm) %>%
    gather(key = "model",value = "nMAE", -c(class1,distance)) %>%
    mutate(model = factor(model,
                                                levels = c("Ensemble","geeSEBAL","PT.JPL","SSEBop",
                                                                     "eeMETRIC","DisALEXI","SIMS"))
                 ) %>%
    ggplot(aes(distance,nMAE*100,color = model, shape=model)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(
        y = "Normalized MAE (%)", 
        x = "Temporal distance from satellite (days)"
    ) +
    facet_wrap(~class1)
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values because more
    ## than 6 becomes difficult to discriminate
    ## ℹ you have requested 7 values. Consider specifying shapes manually if you need
    ##   that many have them.

    ## Warning: Removed 165 rows containing missing values or values outside the scale range
    ## (`geom_line()`).

    ## Warning: Removed 198 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
result %>%
    filter(class1 == "Croplands") %>%
    select(!n) %>%
    select(!ETm) %>%
    gather(key = "model",value = "nMAE", -c(class1,distance)) %>%
    mutate(model = factor(model,
                                                levels = c("Ensemble","geeSEBAL","PT.JPL","SSEBop",
                                                                     "eeMETRIC","DisALEXI","SIMS"))
                 ) %>%
    ggplot(aes(distance,nMAE*100,color = model, shape=model)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(title = "crop",
        y = "Normalized MAE (%)", 
        x = "Temporal distance from satellite (days)"
    )
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values because more
    ## than 6 becomes difficult to discriminate
    ## ℹ you have requested 7 values. Consider specifying shapes manually if you need
    ##   that many have them.

    ## Warning: Removed 33 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->

## Normalized MBE

``` r
mbe <- function (actual, predicted) {return(mean(predicted - actual,na.rm=T))}

for(i in c(0:32)){

    rmse_df <- df %>%
        ungroup() %>%
        group_by(class1) %>%
        filter(nearest_distance == i & ET_OBS == T & !is.na(GF_geeSEBAL)) %>%
        summarise(ETm = mean(ET_corr,na.rm=T),
                            
                            geeSEBAL = mbe(ET_corr,GF_geeSEBAL)/ETm,
                            PT.JPL = mbe(ET_corr,GF_PT.JPL)/ETm,
                            SSEBop = mbe(ET_corr,GF_SSEBop)/ETm,
                            eeMETRIC = mbe(ET_corr,GF_eeMETRIC)/ETm,
                            DisALEXI = mbe(ET_corr,GF_DisALEXI)/ETm,
                            SIMS = mbe(ET_corr,GF_SIMS)/ETm,
                            Ensemble = mbe(ET_corr,GF_Ensemble)/ETm,
                            
                            distance = i,
                            n = n()
        ) 
    
    if(i == 0){
        result <- rmse_df
    } else{
        result <- rbind(result,rmse_df)
    }
}

result %>%
    select(!n) %>%
    select(!ETm) %>%
    gather(key = "model",value = "nMBE", -c(class1,distance)) %>%
    mutate(model = factor(model,
                                                levels = c("Ensemble","geeSEBAL","PT.JPL","SSEBop",
                                                                     "eeMETRIC","DisALEXI","SIMS"))
                 ) %>%
    ggplot(aes(distance,nMBE*100,color = model, shape=model)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(
        y = "Normalized MBE (%)", 
        x = "Temporal distance from satellite (days)"
    ) +
    facet_wrap(~class1, scales = "free_y")
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values because more
    ## than 6 becomes difficult to discriminate
    ## ℹ you have requested 7 values. Consider specifying shapes manually if you need
    ##   that many have them.

    ## Warning: Removed 165 rows containing missing values or values outside the scale range
    ## (`geom_line()`).

    ## Warning: Removed 198 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
result %>%
    select(!n) %>%
    select(!ETm) %>%
    gather(key = "model",value = "nMBE", -c(class1,distance)) %>%
    mutate(model = factor(model,
                                                levels = c("Ensemble","geeSEBAL","PT.JPL","SSEBop",
                                                                     "eeMETRIC","DisALEXI","SIMS"))
                 ) %>%
    ggplot(aes(distance,nMBE*100,color = model, shape=model)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(
        y = "Normalized MBE (%)", 
        x = "Temporal distance from satellite (days)"
    ) +
    facet_wrap(~class1)
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values because more
    ## than 6 becomes difficult to discriminate
    ## ℹ you have requested 7 values. Consider specifying shapes manually if you need
    ##   that many have them.

    ## Warning: Removed 165 rows containing missing values or values outside the scale range
    ## (`geom_line()`).

    ## Warning: Removed 198 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
result %>%
    filter(class1 == "Croplands") %>%
    select(!n) %>%
    select(!ETm) %>%
    gather(key = "model",value = "nMBE", -c(class1,distance)) %>%
    mutate(model = factor(model,
                                                levels = c("Ensemble","geeSEBAL","PT.JPL","SSEBop",
                                                                     "eeMETRIC","DisALEXI","SIMS"))
                 ) %>%
    ggplot(aes(distance,nMBE*100,color = model, shape=model)) +
    geom_line() + 
    geom_point() +
    theme_bw() +
    labs(title = "crop",
        y = "Normalized MBE (%)", 
        x = "Temporal distance from satellite (days)"
    )
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values because more
    ## than 6 becomes difficult to discriminate
    ## ℹ you have requested 7 values. Consider specifying shapes manually if you need
    ##   that many have them.

    ## Warning: Removed 33 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](3_analysis_SFE_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->
