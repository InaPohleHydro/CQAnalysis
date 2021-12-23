##Required packages
#library(trend) 
library(Kendall) #for Kendall test used in trend_detection()
library(hydroTSM) #for flow duration curve used in quantiles_from_fdc() and 
library(minpack.lm) #for fitting concentration discharge models in 
library(hydrostats) #for baseflow separation used in rising_peak_falling()

#Analyses cQ type and chemostatic vs. chemodynamic behaviour - integrates all other functions
analyse_type_cQ_and_chemostatic_chemodynamic <- function(cQ_input_data){ # data frame containing Date, Discharge, concentration
  minimum_valid_datasets <- 10 #the minimum of paired concentration and discharge data sets. if below minimum, the computation will not be carried out
  names(cQ_input_data) <- c("Date","Discharge","Concentration")
  cQ_data_with_rising_peak_falling <- rising_peak_falling(cQ_data = cQ_input_data)
  print("FINISHED STEP 1: rising, peak, falling")
  rising_cQ_parameters <- concentration_discharge_hydrograph_segmentation_all_quantiles(concentration = cQ_data_with_rising_peak_falling$Concentration[cQ_data_with_rising_peak_falling$rising_peak_falling %in% c("rising","peak")], discharge = cQ_data_with_rising_peak_falling$Discharge[cQ_data_with_rising_peak_falling$rising_peak_falling %in% c("rising","peak")])
  rising_cQ_parameters_best <- rising_cQ_parameters[rising_cQ_parameters$r_both == max(rising_cQ_parameters$r_both),][1,]
  
  falling_cQ_parameters <- concentration_discharge_hydrograph_segmentation_all_quantiles(concentration = cQ_data_with_rising_peak_falling$Concentration[cQ_data_with_rising_peak_falling$rising_peak_falling %in% c("falling")], discharge = cQ_data_with_rising_peak_falling$Discharge[cQ_data_with_rising_peak_falling$rising_peak_falling %in% c("falling")])
  falling_cQ_parameters_best <- falling_cQ_parameters[falling_cQ_parameters$r_both == max(falling_cQ_parameters$r_both),][1,]
  if(rising_cQ_parameters_best$r_both > 0 && falling_cQ_parameters$r_both > 0){
  print("FINISHED STEP 2: cQ Parameters")
  
  rising_cQ <- modelled_concentration_from_discharge(discharge = cQ_data_with_rising_peak_falling$Discharge[cQ_data_with_rising_peak_falling$rising_peak_falling %in% c("rising","peak")], a_lower = rising_cQ_parameters_best$a_lower, b_lower = rising_cQ_parameters_best$b_lower, a_higher = rising_cQ_parameters_best$a_higher, b_higher = rising_cQ_parameters_best$b_higher, segmentation_quantile = rising_cQ_parameters_best$quantiles_considered)
  falling_cQ <- modelled_concentration_from_discharge(discharge = cQ_data_with_rising_peak_falling$Discharge[cQ_data_with_rising_peak_falling$rising_peak_falling %in% c("falling")], a_lower = falling_cQ_parameters_best$a_lower, b_lower = falling_cQ_parameters_best$b_lower, a_higher = falling_cQ_parameters_best$a_higher, b_higher = falling_cQ_parameters_best$b_higher, segmentation_quantile = falling_cQ_parameters_best$quantiles_considered)
  print("FINISHED STEP 3: Modelled concentrations")
  
  cQ_type <- type_concentration_discharge(rising_modelled_concentration = rising_cQ$Concentration,
                                         falling_modelled_concentration = falling_cQ$Concentration)
  print("FINISHED STEP 4: cQ Type")
  
  cd_cs <- chemodynamic_chemostatic(concentration = cQ_input_data$Concentration, discharge = cQ_input_data$Discharge)
  print("FINISHED STEP 5: Chemostatic / chemodynamic behaviour")
  results <- c(cQ_type, cd_cs)
  print("ANALYSIS COMPLETED")
  
  print(paste("RESULT: cQ Type is",cQ_type,"and behaviour is",cd_cs))
  results} else {
    print("too few paired concentration and discharge data for analysis")}
}

##Chemostatic versus chemodynamic behaviour based on the definition by Musolff et al. 2017 (https://doi.org/10.1002/2017GL072630)
chemodynamic_chemostatic <- function(concentration, discharge ){
  discharge <- discharge[!is.na(discharge)]
  concentration <- concentration[!is.na(concentration)]
  cv_discharge <- sd(discharge)/mean(discharge)
  cv_concentration <- sd(concentration)/mean(concentration)
  chemostatic.chemodynamic <- NA
  chemostatic.chemodynamic[cv_concentration/cv_discharge > 0.5] <- "chemodynamic"
  chemostatic.chemodynamic[cv_concentration/cv_discharge <= 0.5] <- "chemostatic"
  chemostatic.chemodynamic  
}
#Calculation of concentration discharge models with hydrograph segmentation according to Moatar et al. 2017 https://doi.org/10.1002/2016WR019635
concentration_discharge_hydrograph_segmentation <- function(concentration, discharge){
  data <- data.frame(concentration,discharge)
  data <- data[!is.na(data$Concentration),]  
  data <- data[!is.na(data$Discharge),]  
  Q50 <- quantile_from_fdc(data$Discharge, 0.5)
  low <- data[data$Discharge < Q50,]
  high <- data[data$Discharge > Q50,]
  curve_nlslrc <- nlsLM(concentration ~ a * discharge**b,
                        start=list(a=(mean(data$Concentration/data$Discharge)),  b=-0.5),data = low) #fitting curve of the shape c =  a * Q **b
  powerlawfunction_low <- nls(concentration~a*discharge**b, 
                              start = list(a = coef(curve_nlslrc)[1], b = coef(curve_nlslrc)[2]),
                              data = low, control = list(maxiter = 1000, warnOnly = TRUE))
  powerlawfunction_low_summary <-  summary(powerlawfunction_low)  
  a_low <- powerlawfunction_low_summary$parameters[1,1]
  b_low <- powerlawfunction_low_summary$parameters[2,1]
  #significant slope
  significant_slope_low <- powerlawfunction_low_summary$parameters[2,4] < 0.05
  low$estimate <-   powerlawfunction_low_summary$coefficients[1,1] * low$Discharge**powerlawfunction_low_summary$coefficients[2,1]
  r_low <- cor(low$Concentration, low$estimate)
  curve_nlslrc <- nlsLM(concentration ~ a * discharge**b,
                  start=list(a=(mean(data$Concentration/data$Discharge)), b=-0.5),data = high) #fitting curve of the shape c =  a * Q **b
  powerlawfunction_high <- nls(concentration~a*discharge**b, 
                               start = list(a = coef(curve_nlslrc)[1], b = coef(curve_nlslrc)[2]),
                               data = high, control = list(maxiter = 1000, warnOnly = TRUE))
  powerlawfunction_high_summary <-  summary(powerlawfunction_high)  
  a_high <- powerlawfunction_high_summary$parameters[1,1]
  b_high <- powerlawfunction_high_summary$parameters[2,1]
  #significant slope
  significant_slope_high <- powerlawfunction_high_summary$parameters[2,4] < 0.05
  high$estimate <-   powerlawfunction_high_summary$coefficients[1,1] * high$Discharge**powerlawfunction_high_summary$coefficients[2,1]
  r_high <- cor(high$Concentration, high$estimate)
  low_and_high <- rbind(low,high)
  r_low_and_high <- cor(low_and_high$Concentration, low_and_high$estimate)
  cQ50_low <-   a_low * Q50**b_low
  cQ50_high <-   a_high * Q50**b_high
  cQ50 <- mean(c(cQ50_low,cQ50_high))
  #classify by types
  low_type <- "flat"
  if (b_low > 0 & significant_slope_low){
    low_type <- "up"
  }
  if (b_low < 0 & significant_slope_low){
    low_type <- "down"
  }
  high_type <- "flat"
  if (b_high > 0 & significant_slope_high){
    high_type <- "up"
  }
  if (b_high < 0 & significant_slope_high){
    high_type <- "down"
  }
  type <- paste(low_type,".",high_type,sep="")
  results <- data.frame(a_low,b_low,as.character(significant_slope_low),
                        r_low,a_high,b_high,as.character(significant_slope_high),
                        r_high,r_low_and_high,cQ50,low_type,high_type,type)
  return(results)
}

#Calculation of concentration discharge models with a modified hydrograph segmentation by separating into very high flow, high flow, medium flow, low flow, very low flow
concentration_discharge_hydrograph_segmentation_five_segments <- function(concentration, discharge ){
  data <- data.frame(concentration,discharge)
  data <- data[!is.na(data$Concentration),]  
  data <- data[!is.na(data$Discharge),]  
  q10 <- quantile(data$Discharge, probs = .10) 
  q25 <- quantile(data$Discharge, probs = .25) 
  q75 <- quantile(data$Discharge, probs = .75) 
  q90 <- quantile(data$Discharge, probs = .90) 
  very_low <- data[data$Discharge < q10,]
  medium_low <- data[data$Discharge >= q10 & data$Discharge < q25,]
  low <- data[data$Discharge < q25,]
  medium <- data[data$Discharge >= q25 & data$Discharge < q75,]
  high <- data[data$Discharge >= q75,]
  medium_high <- data[data$Discharge >= q75 & data$Discharge < q90,]
  very_high <- data[data$Discharge >= q90,]
  if (length(very_low[,1]) > minimum_valid_datasets & length(medium_low[,1]) > minimum_valid_datasets & length(medium_high[,1]) > minimum_valid_datasets & length(very_high[,1]) > minimum_valid_datasets){
    curve_nlslrc <- nlsLM(concentration ~ a * discharge**b, 
                         start=list(a=(mean(data$Concentration/data$Discharge)),
                                    b=-0.5),data = low) #fitting curve of the shape c =  a * Q **b
    cQ_very_low  <- nls(concentration~a*discharge**b, 
                        start = list(a = coef(curve_nlslrc)[1], b = coef(curve_nlslrc)[2]), data =  very_low , control = list(maxiter = 1000, warnOnly = TRUE))
    cQ_very_low_summary <-  summary(cQ_very_low)  
    a_very_low <- cQ_very_low_summary$parameters[1,1]
    b_very_low <- cQ_very_low_summary$parameters[2,1]
    significant_b_very_low <- cQ_very_low_summary$parameters[2,4] < 0.05
    very_low$estimate <-   a_very_low * very_low$Discharge**b_very_low
    r_very_low <- cor(very_low$Concentration, very_low$estimate)
    very_low_type <- "flat"
    if (b_very_low > 0 & significant_b_very_low){
      very_low_type <- "up"
    }
    if (b_very_low < 0 & significant_b_very_low){
      very_low_type <- "down"
    }
    curve_nlslrc <- nlsLM(concentration ~  a * discharge**b,
                         start=list(a=(mean(data$Concentration/data$Discharge)),
                                    b=0),data = medium_low) #fitting curve of the shape c =  a * Q **b
    cQ_medium_low  <- nls(concentration~a*discharge**b, start = list(a = coef(curve_nlslrc)[1], b = coef(curve_nlslrc)[2]), data =  medium_low , control = list(maxiter = 1000, warnOnly = TRUE))
    cQ_medium_low_summary <-  summary(cQ_medium_low)  
    a_medium_low <- cQ_medium_low_summary$parameters[1,1]
    b_medium_low <- cQ_medium_low_summary$parameters[2,1]
    significant_b_medium_low <- cQ_medium_low_summary$parameters[2,4] < 0.05
    medium_low$estimate <-   a_medium_low * medium_low$Discharge**b_medium_low
    r_medium_low <- cor(medium_low$Concentration, medium_low$estimate)
    medium_low_type <- "flat"
    if (b_medium_low > 0 & significant_b_medium_low){
      medium_low_type <- "up"}
    if (b_medium_low < 0 & significant_b_medium_low){
      medium_low_type <- "down"}
    conc_medium_very_low <- mean(c(medium_low$estimate[medium_low$Discharge == min(medium_low$Discharge)],very_low$estimate[very_low$Discharge == max(very_low$Discharge)]))
    curve_nlslrc <- nlsLM(concentration ~ a * discharge**b,
                         start=list(a=(mean(data$Concentration/data$Discharge)),
                                    b=-0.5),data = low) #fitting curve of the shape c =  a * Q **b
    cQ_low  <- nls(concentration~a*discharge**b, start = list(a = coef(curve_nlslrc)[1], b = coef(curve_nlslrc)[2]), data =  low , control = list(maxiter = 1000, warnOnly = TRUE))
    cQ_low_summary <-  summary(cQ_low)  
    a_low <- cQ_low_summary$parameters[1,1]
    b_low <- cQ_low_summary$parameters[2,1]
    significant_b_low <- cQ_low_summary$parameters[2,4] < 0.05
    low$estimate <-   a_low * low$Discharge**b_low
    r_low <- cor(low$Concentration, low$estimate)
    low_type <- "flat"
    if (b_low > 0 & significant_b_low){
      low_type <- "up"}
    if (b_low < 0 & significant_b_low){
      low_type <- "down"}
    curve_nlslrc <- nlsLM(concentration ~ a * discharge**b,
                         start=list(a=(mean(data$Concentration/data$Discharge)),
                                    b=-0.5),data = medium) #fitting curve of the shape c =  a * Q **b
    cQ_medium <- nls(concentration~a*discharge**b, start = list(a = coef(curve_nlslrc)[1], b = coef(curve_nlslrc)[2]), data =  medium , control = list(maxiter = 1000, warnOnly = TRUE))
    cQ_medium_summary <-  summary(cQ_medium)  
    a_medium <- cQ_medium_summary$parameters[1,1]
    b_medium <- cQ_medium_summary$parameters[2,1]
    significant_b_medium <- cQ_medium_summary$parameters[2,4] < 0.05
    medium$estimate <-   a_medium * medium$Discharge**b_medium
    r_medium <- cor(medium$Concentration, medium$estimate)
    medium_type <- "flat"
    if (b_medium > 0 & significant_b_medium){
      medium_type <- "up"}
    if (b_medium < 0 & significant_b_medium){
      medium_type <- "down"}
    conc_low_medium <- mean(c(medium$estimate[medium$Discharge == min(medium$Discharge)],low$estimate[low$Discharge == max(low$Discharge)]))
    curve_nlslrc <- nlsLM(concentration ~ a * discharge**b,
                         start=list(a=(mean(data$Concentration/data$Discharge)),
                                    b=-0.5),data = high) #fitting curve of the shape c =  a * Q **b
    cQ_high  <- nls(concentration~a*discharge**b, start = list(a = coef(curve_nlslrc)[1], b = coef(curve_nlslrc)[2]), data =  high , control = list(maxiter = 1000, warnOnly = TRUE))
    cQ_high_summary <-  summary(cQ_high)  
    a_high <- cQ_high_summary$parameters[1,1]
    b_high <- cQ_high_summary$parameters[2,1]
    significant_b_high <- cQ_high_summary$parameters[2,4] < 0.05
    high$estimate <-   a_high * high$Discharge**b_high
    r_high <- cor(high$Concentration, high$estimate)
    high_type <- "flat"
    if (b_high > 0 & significant_b_high){
      high_type <- "up"
    }
    if (b_high < 0 & significant_b_high){
      high_type <- "down"
    }
    conc_medium_high <- mean(c(medium$estimate[medium$Discharge == max(medium$Discharge)],high$estimate[high$Discharge == min(high$Discharge)]))
    curve_nlslrc <- nlsLM(concentration ~ a * discharge**b,
                         start=list(a=(mean(data$Concentration/data$Discharge)),
                                    b=-0.5),data = medium_high) #fitting curve of the shape c =  a * Q **b
    cQ_medium_high  <- nls(concentration~a*discharge**b, start = list(a = coef(curve_nlslrc)[1], b = coef(curve_nlslrc)[2]), data =  medium_high , control = list(maxiter = 1000, warnOnly = TRUE))
    cQ_medium_high_summary <-  summary(cQ_medium_high)  
    a_medium_high <- cQ_medium_high_summary$parameters[1,1]
    b_medium_high <- cQ_medium_high_summary$parameters[2,1]
    significant_b_medium_high <- cQ_medium_high_summary$parameters[2,4] < 0.05
    medium_high$estimate <-   a_medium_high * medium_high$Discharge**b_medium_high
    r_medium_high <- cor(medium_high$Concentration, medium_high$estimate)
    medium_high_type <- "flat"
    if (b_medium_high > 0 & significant_b_medium_high){
      medium_high_type <- "up"
    }
    if (b_medium_high < 0 & significant_b_medium_high){
      medium_high_type <- "down"
    }
    curve_nlslrc <- nlsLM(concentration ~ a * discharge**b,
                         start=list(a=(mean(data$Concentration/data$Discharge)),
                                    b=-0.5),data = very_high) #fitting curve of the shape c =  a * Q **b
    cQ_very_high  <- nls(concentration~a*discharge**b, start = list(a = coef(curve_nlslrc)[1], b = coef(curve_nlslrc)[2]), data =  very_high , control = list(maxiter = 1000, warnOnly = TRUE))
    cQ_very_high_summary <-  summary(cQ_very_high)  
    a_very_high <- cQ_very_high_summary$parameters[1,1]
    b_very_high <- cQ_very_high_summary$parameters[2,1]
    significant_b_very_high <- cQ_very_high_summary$parameters[2,4] < 0.05
    very_high$estimate <-   a_very_high * very_high$Discharge**b_very_high
    r_very_high <- cor(very_high$Concentration, very_high$estimate)
    very_high_type <- "flat"
    if (b_very_high > 0 & significant_b_very_high){
      very_high_type <- "up"
    }
    if (b_very_high < 0 & significant_b_very_high){
      very_high_type <- "down"
    }
conc_medium_very_high <- mean(c(medium$estimate[medium_high$Discharge == max(medium_high$Discharge)],very_high$estimate[very_high$Discharge == min(very_high$Discharge)]))
five_segments <- rbind(very_low,medium_low,medium,medium_high,very_high)
correlation_five_segments <- cor(five_segments$Concentration,five_segments$estimate)
three_segments <- rbind(low,medium,high)
correlation_three_segments <- cor(three_segments$Concentration,three_segments$estimate)
type_five_segments <- paste(very_low_type,medium_low_type,medium_type,medium_high_type,very_high_type,sep = ".")
type_three_segments <- paste(low_type,medium_type,high_type,sep = ".")
    
results <- data.frame(a_very_low,b_very_low,as.character(significant_b_very_low),r_very_low,very_low_type,
                          a_medium_low,b_medium_low,as.character(significant_b_medium_low),r_medium_low,medium_low_type,
                          a_low,b_low,as.character(significant_b_low),r_low,low_type,
                          a_medium,b_medium,as.character(significant_b_medium),r_medium,medium_type,
                          a_high,b_high,as.character(significant_b_high),r_high,high_type,
                          a_medium_high,b_medium_high,as.character(significant_b_medium_high),r_medium_high,medium_high_type,
                          a_very_high,b_very_high,as.character(significant_b_very_high),r_very_high,very_high_type,
                          conc_medium_very_low,conc_low_medium,conc_medium_high,        conc_medium_very_high,type_five_segments,correlation_five_segments,type_three_segments,correlation_three_segments)}   else {
                            results <- data.frame(a_very_low,b_very_low,significant_b_very_low,r_very_low,very_low_type,
                                                  a_medium_low,b_medium_low,significant_b_medium_low,r_medium_low,medium_low_type,
                                                  a_low,b_low,significant_b_low,r_low,low_type,
                                                  a_medium,b_medium,significant_b_medium,r_medium,medium_type,
                                                  a_high,b_high,significant_b_high,r_high,high_type,
                                                  a_medium_high,b_medium_high,significant_b_medium_high,r_medium_high,medium_high_type,
                                                  a_very_high,b_very_high,significant_b_very_high,r_very_high,very_high_type,
                                                  conc_medium_very_low,conc_low_medium,conc_medium_high,        conc_medium_very_high,type_five_segments,correlation_five_segments,type_three_segments,correlation_three_segments) * NA
              }
  return(results)
}

#Calculation of concentration discharge models with hydrograph segmentation for all discharge quantiles
concentration_discharge_hydrograph_segmentation_all_quantiles <- function(concentration, discharge){
  options(warn=-1)
  data <- data.frame(discharge,concentration)  
  names(data) <- c("Discharge","Concentration")
  data <- data[!is.na(data$Concentration),]  
  data <- data[!is.na(data$Discharge),]  
    quantile_variation <- data.frame()  
    for (quantiles_considered in c(1:99)){
      quantile_discharge <- quantile_from_fdc(data$Discharge, quantiles_considered/100)
    if (quantile_discharge < max(data$Discharge) & quantile_discharge > min(data$Discharge)){
      lower <- data[data$Discharge < quantile_discharge,]
      higher <- data[data$Discharge >= quantile_discharge,]
      #there is a problem when concentrations are always similar (at the detection limit)
      #add some  noise
      if (length(unique(lower$Concentration)) > minimum_valid_datasets & length(unique(higher$Concentration)) > minimum_valid_datasets 
          #& length(unique(lower$Concentration)) > length((lower$Concentration))/5  & length(unique(higher$Concentration)) > length((higher$Concentration))/5  
      ){
        for (i in c(1:length(lower$Concentration))){
          lower$Concentration[i] <- lower$Concentration[i] + min(lower$Concentration)*0.1*rnorm(1)
          lower$Discharge[i] <- abs(lower$Discharge[i] + min(lower$Discharge)*0.1*rnorm(1))
        }
        for (i in c(1:length(higher$Concentration))){
          higher$Concentration[i] <- higher$Concentration[i] + min(higher$Concentration)*0.1*rnorm(1)
          higher$Discharge[i] <- higher$Discharge[i] + min(lower$Discharge)*0.1*rnorm(1)
        }
        curve_nlslrc <- nlsLM(concentration ~ a * abs(discharge)**b,
                             start=list(a=(mean(lower$Concentration/abs(lower$Discharge))),
                                        b=-0.5),data = lower) #fitting curve of the shape c =  a * Q **b
        cQ_lower  <- nls(concentration~a*discharge**b, start = list(a = coef(curve_nlslrc)[1], b = coef(curve_nlslrc)[2]), data =  lower , control = list(maxiter = 10000, warnOnly = TRUE))
        cQ_lower_summary <-  summary(cQ_lower)  
        a_lower <- cQ_lower_summary$parameters[1,1]
        b_lower <- cQ_lower_summary$parameters[2,1]
        significant_b_lower <- cQ_lower_summary$parameters[2,4] < 0.05
        lower$estimate <-   a_lower * lower$Discharge**b_lower
        r_lower <- cor(lower$Concentration, lower$estimate)
        lower_type <- "flat"
        if (b_lower > 0 & significant_b_lower){
          lower_type <- "up"
        }
        if (b_lower < 0 & significant_b_lower){
          lower_type <- "down"
        }
        curve_nlslrc <- nlsLM(concentration ~  a * discharge**b,
                             start=list(a=(mean(higher$Concentration/higher$Discharge)),
                                        b=-0.5),data = higher) #fitting curve of the shape c =  a * Q **b
        cQ_higher  <- nls(concentration~a*discharge**b, start = list(a = coef(curve_nlslrc)[1], b = coef(curve_nlslrc)[2]), data =  higher , control = list(maxiter = 10000, warnOnly = TRUE))
        cQ_higher_summary <-  summary(cQ_higher)  
        a_higher <- cQ_higher_summary$parameters[1,1]
        b_higher <- cQ_higher_summary$parameters[2,1]
        significant_b_higher <- cQ_higher_summary$parameters[2,4] < 0.05
        higher$estimate <-   a_higher * higher$Discharge**b_higher
        r_higher <- cor(higher$Concentration, higher$estimate)
        higher_type <- "flat"
        if (b_higher > 0 & significant_b_higher){
          higher_type <- "up"
        }
        if (b_higher < 0 & significant_b_higher){
          higher_type <- "down"
        }
        conc_low_high <- mean(c(lower$estimate[lower$Discharge == max(lower$Discharge)],higher$estimate[higher$Discharge == min(higher$Discharge)]))
        both <- rbind(higher,lower)
        r_both <- cor(both$Concentration,both$estimate)
        both_type <- paste(lower_type,higher_type,sep=".")
        curr_quantile <- data.frame(a_lower,b_lower,significant_b_lower,r_lower,
                                    a_higher,b_higher,significant_b_higher,r_higher,r_both,lower_type,higher_type,both_type,quantiles_considered,quantile_discharge)
        quantile_variation <- rbind(quantile_variation,curr_quantile)
      }
    }
  }
    if (length(quantile_variation) == 0){
    quantile_variation <- data.frame(message = "too few paired concentration and discharge data for analysis", r_both = 0)
  }
results <- quantile_variation
options(warn=0)

  return(results)
}

###Modelled concentration from discharge using C-Q relationship based on hydrograph segmentation into two segments
modelled_concentration_from_discharge <- function(discharge, a_lower, a_higher, b_lower, b_higher, segmentation_quantile){
  modelled_concentrations <- data.frame(1:100)
  names(modelled_concentrations)[1] <- "Q_quantile"
for (i in c(1:100)){
    modelled_concentrations$Q[i] <- quantile_from_fdc(discharge = discharge, specified_quantile =  modelled_concentrations$Q_quantile[i]/100)
}
  modelled_concentrations$Concentration <- NA
  modelled_concentrations$Concentration[modelled_concentrations$Q_quantile < segmentation_quantile] <- a_lower * modelled_concentrations$Q[modelled_concentrations$Q_quantile < segmentation_quantile] ** b_lower 
  modelled_concentrations$Concentration[modelled_concentrations$Q_quantile >= segmentation_quantile] <- a_higher * modelled_concentrations$Q[modelled_concentrations$Q_quantile >= segmentation_quantile] ** b_higher
  modelled_concentrations
}


#Determine whether discharge is on rising or falling limb or peak discharge based on baseflow separation. 
#The function derives flow events. Flow events are consecutive intervals where baseflow exceeded (minus a tolerance as the function would otherwise not be able to distinguish between two small flow events in a row) and which have a peak that is higher than a defined minimim peak quantile (otherwhise we would have a number of very small events.) 
#The periods in between the flow events are assigned to the previous flow event (when before the "bottom") or to the next flow event (after the "bottom")
#For each event the peak is derived and the timeseries split into time before and after the peak - rising & falling limb

#maybe this is a thing for later & for assessing the critical dry period duration 
rising_peak_falling <- function(cQ_data, #data.frame containing Date, Discharge, concentration
                                tolerance_quantile = 0.90){
  cQ_data <- cQ_data[!is.na(cQ_data$Discharge),]
  cQ_data <- cQ_data[!duplicated(cQ_data),]
  cQ_data$Date <- as.POSIXct(cQ_data$Date)
  date_Q <- data.frame(Date = cQ_data$Date, Q = cQ_data$Discharge) 
  date_Q$Date <- as.POSIXct(date_Q$Date)
  date_Q <- date_Q[!is.na(date_Q$Q),]
  date_Q <- date_Q[!duplicated(date_Q),]
  q_base <- baseflows(date_Q, a = 0.975, n.reflected = 30, ts = "daily") # separate baseflow using baseflows function in hydrostats
 #events are times when discharge is above baseflow
    q_base$greater_than_baseflow <- (q_base$Q > (q_base$bf + quantile(q_base$Q,(1-tolerance_quantile),na.rm = TRUE)))
    #part behind the "+" solves issue (b) with adding some tolerance when discharge almost goes down to baseflow based on the tolerance quantile
  #each event starts the day after the lowest point when discharge is below baseflow and ends at the next lowest point (arbitrary, but need to make some )
  #for each event: max flow & rising and falling limb
  q_base$event_number <- NA
  q_base$between_event_number <- NA
  event_number <- 0
  between_event_number <- 1
  #those where q > q_base are simple to allocate to events 
  #the times in between have to be allocated in the next iteration
  first_q_base <- which(!q_base$greater_than_baseflow)[1] # we need to start at the first encounter of q = q_base
##define events and periods in between
  for (time in c((1+first_q_base):nrow(q_base))){
    if (!q_base$greater_than_baseflow[time]){
      q_base$between_event_number[time] <- between_event_number
      if (q_base$greater_than_baseflow[time-1]){
        event_number <- event_number + 1
      }}
    if (q_base$greater_than_baseflow[time]){
      q_base$event_number[time] <- event_number
      if (!q_base$greater_than_baseflow[time-1]){
        between_event_number <- between_event_number + 1
      }
    }
  }
#assign time between events to events (before bottom: previous events, after bottom: next event)
  for (between in c(1:(max(between_event_number) - 1))){
    event_before <- as.numeric(rownames(q_base[!is.na(q_base$between_event_number) & q_base$between_event_number == between,]))[1] - 1 #last event interval before between period
    if (!is.na(event_before)){
      event_after<- tail(as.numeric(rownames(q_base[!is.na(q_base$between_event_number) & q_base$between_event_number == between,])),1) + 1 #first event interval after between period
      between_period <- q_base[!is.na(q_base$between_event_number) & q_base$between_event_number == between,]
      bottom <- as.numeric(rownames(between_period[between_period$Q == min(between_period$Q),])[1]) #chose the first point in the between period with lowest discharge
      q_base[!is.na(q_base$between_event_number) & q_base$between_event_number == between & as.numeric(rownames(q_base)) <= bottom,]$event_number <- q_base[event_before,]$event_number #assign number of last event
      if (nrow(between_period[as.numeric(rownames(between_period))>bottom,]) > 0){
        q_base[!is.na(q_base$between_event_number) & q_base$between_event_number == between & as.numeric(rownames(q_base)) > bottom,]$event_number <- q_base[event_after,]$event_number
      }
    }
  }
#for each event determine maximum and assign rising and falling limb
 q_base$rising_peak_falling <- NA
  for (event_number in c(1:(max(event_number) - 1))){
    event <- q_base[!is.na(q_base$event_number) & q_base$event_number == event_number,] 
    peak <- as.numeric(rownames(event[event$Q == max(event$Q),])[1])
 if (nrow(event[as.numeric(rownames(event))<peak,]) > 0){
      q_base[!is.na(q_base$event_number) & q_base$event_number == event_number & as.numeric(rownames(q_base)) < peak,]$rising_peak_falling <- "rising"}
    q_base[!is.na(q_base$event_number) & q_base$event_number == event_number & as.numeric(rownames(q_base)) == peak,]$rising_peak_falling <- "peak"
    if (nrow(event[as.numeric(rownames(event))>peak,]) > 0){
      q_base[!is.na(q_base$event_number) & q_base$event_number == event_number & as.numeric(rownames(q_base)) > peak,]$rising_peak_falling <- "falling"}
  }
  cQ_data_with_rising_peak_falling <- cbind(cQ_data,q_base)
  cQ_data_with_rising_peak_falling <- cQ_data_with_rising_peak_falling[,!names(cQ_data_with_rising_peak_falling) %in%  c("Q","Date1")]
  return(cQ_data_with_rising_peak_falling)  
}

##quantile from FDC: Determines discharge quantile for discharge values using the fdc function in hydrotsm. Used to fit concentration discharge models
quantile_from_fdc <- function(discharge, specified_quantile ){
  discharge <- na.omit(discharge)
  discharge_fdc <- data.frame(discharge,fdc(discharge,plot = FALSE))
  names(discharge_fdc) <- c("Discharge","fdc_discharge")
  Q_quantile <- max(discharge_fdc[discharge_fdc$fdc_discharge >= specified_quantile,1])  
  return(Q_quantile)
}

##Trend detection: used to assess export regime based on modelled concentrations
trend_detection <- function(data,significance_level = 0.05){
  trend <- MannKendall(data)
  p.value <- as.numeric(trend$sl)
  significant <- p.value < significance_level
  return(significant)
}

##Trend magnitude: used to assess export regime based on modelled concentrationss
trend_magnitude <- function(date,data){
  if (length(na.omit(data))>2){
    lin_mod <- lm(data~date,na.action = na.exclude)
    slope <- as.numeric(lin_mod$coefficients[2])
  } else {
    slope <- NA}
  return(slope)
}

### Type of concentration discharge relationship: Assessed based on modelled concentrations for rising and falling limb (each ordered from high to low discharge)
type_concentration_discharge <- function(rising_modelled_concentration, 
                                         falling_modelled_concentration){
  highflow_lowflow <- 2 #set default values
  rising_falling <- 0 #set default values
  cQ_type <- c("enrichment,clockwise",
              "enrichment,nohysteresis",
              "enrichment,anticlockwise",
              "constancy,clockwise",
              "constancy,nohysteresis",
              "constancy,anticlockwise",
              "dilution,clockwise",
              "dilution,nohysteresis",
              "dilution,anticlockwise")
  #enrichment,clockwise: rising > falling, Highflow > lowflow
  #enrichment,nohysteresis: rising = falling, Highflow > lowflow
  #enrichment,anticlockwise: rising < falling, Highflow > lowflow
  #constancy,clockwise: rising > falling, Highflow = lowflow
  #constancy,nohysteresis: rising = falling, Highflow = lowflow
  #constancy,anticlockwise: rising < falling, Highflow = lowflow
  #dilution,clockwise: rising > falling, Highflow < lowflow
  #dilution,nohysteresis: rising = falling, Highflow < lowflow
  #dilution,anticlockwise: rising < falling, Highflow < lowflow
  ###highflow > lowflow? 
  trend_detection_rising <- trend_detection(rising_modelled_concentration)
  trend_magnitude_rising <- trend_magnitude(date = c(100:1), data = rising_modelled_concentration)
  trend_detection_falling <- trend_detection(falling_modelled_concentration)
  trend_magnitude_falling <- trend_magnitude(date = c(100:1), data = falling_modelled_concentration)
    if (trend_detection_rising & trend_detection_falling & 
      trend_magnitude_rising > 0 & trend_magnitude_falling > 0){
    highflow_lowflow <- 1}
  if (trend_detection_rising & trend_detection_falling & 
      trend_magnitude_rising < 0 & trend_magnitude_falling < 0){
    highflow_lowflow <- 3}
    ## rising > falling ? 
 categ <- as.factor(c(rep("r",100),rep("f",100)))
  ktest <- kruskal.test(c(rising_modelled_concentration),
                        falling_modelled_concentration,categ)
  ktest$p.value
  if (ktest$p.value < 0.05 & 
      mean(as.numeric(rising_modelled_concentration)) > mean(as.numeric(falling_modelled_concentration))){
    rising_falling <- -1}
  if (ktest$p.value < 0.05 & 
      mean(as.numeric(rising_modelled_concentration)) < mean(as.numeric(falling_modelled_concentration))){
    rising_falling <- 1 }
  group <- rising_falling + (3 * highflow_lowflow) -1
  cQ_type[group]
}



