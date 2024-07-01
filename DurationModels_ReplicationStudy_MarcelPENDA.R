#' ---
#' title: "Transition & Duration Models - Homework: Replication Study of Lalive et al. (2006)"
#' author: "Marcel Richard Penda"
#' date: "2024-04-10"
#' output:
#'   html_document: default
#'   pdf_document: default
#' ---
#' 
#' # Setup
#' 
#' ## Libraries
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Loading packages 
library(pacman)
p_load(foreign, tidyverse, survival, KernSmooth, 
       tidyverse, data.table, broom, parallel, here, plotly, ggplot2, stargazer, magrittr,skimr,janitor,  tidymodels, ADAPTS, caret, yardstick, rlang, parsnip, sandwich, lmtest, haven, tinytex, rdrobust,dplyr, plotrix, plyr,readxl, usmap, stringr, finalfit, scales,tidyr, gridExtra, patchwork, EventStudy, fixest,kableExtra,wesanderson, gtsummary, maps, cowplot, corrplot, ggcorrplot, ggthemes, wesanderson, mgcv, lmtest)

# Deactivate scientific notation
options(scipen = 999) 

#' 
#' ## Read Data
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
unemp <- read.csv("EconomicDataAustria.csv")  
unemp <- unemp[,1:134] ## get rid of some superfluous variables
unemp <- as_tibble(unemp)
head(unemp)

#' 
#' # Q1: Summary of the Paper
#' The research paper "How Changes in Financial Incentives Affect the Duration of Unemployment" by Lalive, van Ours, and Zweimuller (2006) examines the impact of different changes in unemployment policies on the duration of unemployment periods for different groups of workers. Therefore the authors exploit the Austrian unemployment policy change of 1989 which introduced extensions of unemployment benefits according to individual's characteristics. More precisely, the Austrian government increased the unemployment replacement rate (RR) as well as the potential duration of unemployment benefits (PDB) from 30 to 32 or 52 weeks depending on the individual's wage, age, and work experience. Given the structural policy and the stability of the Austrian economy at that time, the setting can be considered as natural experiment. In addition, the new system is implemented during a phase of employment growth, ensuring that any extension of the duration is not influenced by fluctuations in the labor market.
#' 
#' To conduct their study, two data sources are consulted, providing longitudinal individual data. The authors combine data on the individual's employment/unemployment records and earnings from the Austrian Social Security Database with data on the individual's socioeconomic information from the Austrian unemployment register. From this data, the authors extract 225,821 individuals who entered the state of unemployment between 1987 (two years before policy treatment) and 1991 (two years after policy treatment). About 1% of the data is right-censored, since whether an exit into employment (85%) or a non-job exit destination (14%) was recorded. The median unemployment spell lasts 12 months.
#' 
#' Following the policy regulation, the authors categorize the individuals based on the determinants (i.e. age, income, and experience) of the RR and PBD eligibility, resulting in four different groups of unemployed workers:
#' 
#' (1)  Eligible for PDB: Older individuals (threshold at age 40) earning more than 12,610 ATS and with high experience.
#' 
#' (2)  Eligible for increase in RR: Young (and older) unemployed individuals earning less than 12,610 ATS with low and high experience (with low experience).
#' 
#' (3)  Eligible for PDB and increase in RR: Older individuals earning less than than 12,610 ATS with high experience.
#' 
#' (4)  Control group: Young (and older) unemployed individuals earning more than 12,610 ATS with low and high experience (with low experience).
#' 
#' Consequently, we have heterogeneity in treatment and in the treatment magnitude since the RR and PDB increases are not comparable. Due to additional unemployment benefits for regions dominated by the steel industry, the study concentrates on non-steel industry regions only. Based on this data, the authors first conduct a simple difference-in-difference analysis and estimate the survivor functions as well as exit hazards for the different groups of workers. The key results of their study come from the estimation of a piecewise exponential model (PWE) measuring the treatment effects for different time intervals.
#' 
#' The authors find that an increase in the RR effects the rates of exiting the state of unemployment at the beginning, whereas the effect of extended PBD occurs especially arround the expiration dates of the PBD. Moreover, older individuals are more effected by the PBD and the combined disincentives. Furthermore, the PWE shows that the treatment effects are different for different intervals. 
#' 
#' One of the study's main contributions is the understanding of unemployment incentives. More precisely, the authors conclude that for policy makers to influence/reduce unemployment rates, the PBD represents a more effective tool compared to the RR. Moreover, with their study the authors provide empirical evidence for exiting theoretical labor economic models. While these models provide predictions only about the sign/directions of unemployment incentives, these results provide indication about how large these effects are.
#' 
#' However, before 1989 the Austrian unemployment benefits were rather low compared to other European countries, making the results difficult to transfer to other countries.
#' 
#' # Q2: Difference-in-Differences Analysis 
#' 
#' To replicate the authors' difference-in-differences results, we prepare our dependent variable, i.e. durations in weeks, such that our data is right-censored after 104 weeks. Put differently, all duration values exceeding 104 weeks are set to the duration of 104 weeks such that we obtain a right-censored subsample.
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Computation of average spells when durations are truncated at 104 weeks
unemp %>%
  mutate(dur104 = dur,
         dur104 = ifelse(dur104 > 104, 104, dur104)) -> unemp

#' 
#' In a next step, we calculate the difference-in-differences estimator: 
#' $$
#' \hat{\Delta}_{DD} = (Y^T_A - Y^T_B) - (Y_A^C - Y_B^C) \\ \; \text{where} \; Y^T_A \; \text{and} \; Y^T_B \; \text{denote the average unemployment duration of the treated group} \\ \text{and} \; Y_A^C \; Y_B^C \; \text{denote the average unemployment duration of the control group}
#' $$
#' 
#' To do so, we first compute the average unemployment duration for all treated groups and the control sample. Based on the latter the difference-in-differences estimator is estimated
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create objects with group names
group <- c("PBD", "RR",  "PBD and RR")
group_se <- c("PBD SE", "RR SE", "PBD and RR SE")
group_names <- c("PBD", "PBD SE", "RR", "RR SE", "PBD and RR", "PBD and RR SE", "control", "control SE")

# Initialize DiD_table with appropriate dimensions
DiD_table <- data.frame(matrix(nrow = length(group_names), ncol = 5))
colnames(DiD_table) <- c("Group", "Before 1989", "After 1989", "Change", "Diff-in-Diff")
DiD_table$Group <- group_names

# Average duration for control
Y_CB <- mean(unemp$dur104[unemp$type == "control" & unemp$after == 0])
Y_CA <- mean(unemp$dur104[unemp$type == "control" & unemp$after == 1])
change_control <- Y_CA - Y_CB

# standard errors for control group
se_Y_CB <- sd(unemp$dur104[unemp$type == "control" & unemp$after == 0])/sqrt(length(unemp$dur104[unemp$type == "control" & unemp$after == 0]))
se_Y_CA <- sd(unemp$dur104[unemp$type == "control" & unemp$after == 1])/sqrt(length(unemp$dur104[unemp$type == "control" & unemp$after == 1]))
se_change_control <- sqrt(se_Y_CB^2 + se_Y_CA^2)

# values for control group
values_control <- c(Y_CB, Y_CA, change_control, NA)
se_control <- c(se_Y_CB, se_Y_CA, se_change_control, NA) 

# Assign control group values to the corresponding row in DiD_table
DiD_table[DiD_table$Group == "control", -1] <- values_control
DiD_table[DiD_table$Group == "control SE", -1] <- se_control

# Calculate values for different groups
for (i in group) {
  # Calculate values for different groups
  Y_TB <- mean(unemp$dur104[unemp$type == i & unemp$after == 0])
  Y_TA <- mean(unemp$dur104[unemp$type == i & unemp$after == 1])
  change_var <- Y_TA - Y_TB
  did_var <- change_var - change_control
  
  # calculate standard errors
  se_Y_TB <- sd(unemp$dur104[unemp$type == i & unemp$after == 0])/sqrt(length(unemp$dur104[unemp$type == i & unemp$after == 0]))
  se_Y_TA <- sd(unemp$dur104[unemp$type == i & unemp$after == 1])/sqrt(length(unemp$dur104[unemp$type == i & unemp$after == 1]))
  se_change_var <- sqrt(se_Y_TB^2 + se_Y_TA^2)
  se_did <- sqrt(se_change_var^2 + se_change_control^2)
  
  # Values for variables and standard erros
  values_var <- c(Y_TB,  Y_TA,  change_var,  did_var)
  se_var <- c(se_Y_TB, se_Y_TA, se_change_var, se_did)
  
  # Assign values to corresponding group in DiD_table
  DiD_table[DiD_table$Group == i, -1] <- values_var
  DiD_table[DiD_table$Group == paste(i, "SE", sep = " "), -1] <- se_var
}

print(DiD_table)

#' 
#' # Q3: Survival Function
#' 
#' To reproduce the survival functions for each group, we apply the `survfit` function to our right-censored duration variable `dur104`. We plot the survival functions, i.e. the probability that unemployment spells last longer than x weeks, for the period before (solid lines) and after the policy implementation  (dashed lines) for each group. The function generates the erratic survival probability functions based on the Kaplan-Meier estimator:
#' 
#' $$
#' \hat{S}(t) = \prod_{r=1}^{n} \frac{N_r-E_r}{N_r} \\ 
#' \text{where} \; N_r \; \text{denotes the number of observations in the risk set for the interval} \; r \; \\
#' \text{and} \; E_r \; \text{denotes the number of observations exiting the the risk set in the interval}.
#' $$
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
options(repr.plot.width = 6, repr.plot.height =15)
group <- c("PBD", "RR", "PBD and RR", "control")
for (i in group) {
  # survival function after 1989
  plot(survfit(Surv(dur104, uncc) ~ 1, # do not include any covariates
                     conf.int = 0, # do not show confidence intervals
                     data = unemp[unemp$type == i & unemp$after == 1, ]),
       main = i,
       xlab = "Unemployment duration (weeks)", 
       ylab = "Survival Probability", 
       lty = "solid")
  
  # survival function before 1989
  lines(survfit(Surv(dur104, uncc) ~ 1, # do not include any covariats
                conf.int = 0, # do not show confidence intervals
                data = unemp[unemp$type == i & unemp$after == 0, ]), 
        lty = "dashed")
  
  # horizontal lines
  abline(h = seq(0, 1, by = 0.1), col = "lightgray", lty = "solid")
  
  # Change ticks on axis
  axis(side = 2, at = seq(0, max(1), by = 0.1))
  
  
  # Include legend
  legend("topright", legend = c("After", "Before"), lty = c("solid", "dashed"))
}

#' 
#' For all three groups, we can observe the general trend that the probabilities of workers to survive in the risk set is higher before the change. Put differently, generally, unemployment spells last longer before the policy introduction in 1989 for all four groups. 
#' 
#' A more detailed analysis yields that after 20 weeks, the survival function of the PBD group after the policy change shows higher values, i.e. less exits from unemployment, compared to the survival function of the PBD group before 1989. The difference in the duration of unemployment spells slowly diminishes and disappears after 65 weeks. Hence, we can conclude that the effects of the policy change, i.e. longer unemployment spells, can be mainly observed about 15 weeks after a worker became unemployed and last until the 65th week. This is in line with the difference-in-differences estimations, which indicate that the PBD policy change lead to an average increase of 1.07% in unemployment duration. We can observe similar effects for the RR and PBD+RR groups also supporting the previous diff-in-diff findings.
#' 
#' The survival functions of the control group before and after the policy implementation show almost no differences.This suggests that workers who were not disincentivized by any increase in the potential duration of unemployment benefits (PBD), higher replacement rates (RR), or both (PBD+RR) did not change their job search behavior.
#' 
#' 
#' # Q4: KM estimates of the unemployment exit hazard
#' 
#' $$
#' \lambda = -\frac{d\ln(1-F(t))}{dt}
#' $$
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Initiate objects containing groups and indicator for policy implementation
groups <- c("PBD", "RR", "PBD and RR", "control")
policy = c(0,1)

par(mfrow = c(1,1))

for (group in groups){
  for (state in policy){
    km_fit <- survfit(Surv(dur104, uncc) ~ 1, # do not include any covariates
                       conf.int = 0, # do not show confidence intervals
                       data = unemp[unemp$type == group & unemp$after == state, ])
    km_summary <- summary(km_fit)
    
    # To make discrete function less erratic, pick smaller number of points, here one for each week
    ## Calculate the number of weeks
    num_weeks <- ceiling(max(km_summary$time))
    
    ## Initialize vectors to store selected survival probabilities and time points
    surv_prob <- numeric(num_weeks)
    time_pts <- numeric(num_weeks)
    
    ## Select one value per week
      for (week in 1:num_weeks) {
      # Find the closest time point to the current week
      time_point <- km_summary$time[which.min(abs(km_summary$time - week))]
    
      # Find the index of the closest time point
      index <- which(km_summary$time == time_point)
    
      # Store the selected survival probability and time point
      surv_prob[week] <- km_summary$surv[index]
      time_pts[week] <- time_point
      assign(paste("time_", state, sep = ""), time_pts)
    
      # Calculate the hazards
      hazard <- -diff(log(surv_prob)) / diff(time_pts)
      assign(paste("hazard_", state, sep = ""), hazard)
      }
  }
  
  # Plotting the hazard functions
  ## Replicating Figure 4 in Lalive et al. (erratic graphs)
  plot(time_1[-1], hazard_1, type = "l", 
       main = paste("Hazard Estimates for", sep = " ", group),
       xlab = "Unemplyoment Duration (weeks)", 
       ylab = "Hazard Estimate",
       ylim = c(0, 0.25))
  lines(time_0[-1], hazard_0, lty = "dashed")
  ### Include legend
  legend("topright", legend = c("After", "Before"), lty = c("solid", "dashed"))
  
  
  ## Replicating grpahs from [Christian Schluter’s] M2 project (smoothed graphs)
  ### KernSmooth::locpoly cannot handle NAs, replace NAs with hazard mean
  #### To avoid "Inf" values at the end of the hazard vector
  hazard_1 <- hazard_1[-(length(hazard_1))] # exclude last hazard value which is "Inf"
  hazard_0 <- hazard_0[-(length(hazard_0))] # exclude last hazard value which is "Inf"
  
  
  ### Calculate the mean of non-NA values
  mean_hazard_1 <- mean(hazard_1, na.rm = TRUE) 
  mean_hazard_0 <- mean(hazard_0, na.rm = TRUE)
  ## Replace NA values with the mean
  hazard_1[is.na(hazard_1)] <- mean_hazard_1
  hazard_0[is.na(hazard_0)] <- mean_hazard_0
  
  ### Define polynomial models of degree 3
  poly_model_1 = locpoly(time_1[-c(1, length(time_1))], # exclude last time value, respectively
               hazard_1, 
               degree = 3,
               bandwidth=3)
  
  poly_model_0 = locpoly(time_0[-c(1, length(time_0))], # exclude last time value, respectively
               hazard_0,
               degree = 3,
               bandwidth=3)
  
  ### Plot the polynomial functions
  plot(poly_model_1$x, poly_model_1$y, type = "l",
       main = paste("Smooth Hazard Estimates for", sep = " ", group),
       xlab = "Unemplyoment Duration (weeks)", 
       ylab = "Hazard Estimate",
       ylim = c(0, 0.15))
  lines(poly_model_0$x, poly_model_0$y, lty = "dashed")
  legend("topright", legend = c("After", "Before"), lty = c("solid", "dashed"))
}

#' 
#' 
#' # Q5: Estimte the causal treatment effect in a PH model
#' 
#' The Proportional Hazard (PH) model incoroporates time-invariant covariates such as sex or ethnicity by conditioning on them. One approach to account for these covariate effect, i.e. individual characteristics that impact the hazard but do not change over time, is the Cox's PH model:
#' 
#' $$
#' \lambda(t;x) = K(x)\;\lambda_0(t) \;  \\ \text{where} \; \lambda_0(t) \; \text{denotes the baseline hazard common to all observations}, \\
#' x \; \text{denotes a vetor of covariates}, \\
#' \text{and K(.) describes the effect of} \; x \; \text{on the baseline hazard}.
#' $$
#' In Lalive et al.(2006) the baseline duration dependence is of central interest since it captures the exit rate for an homogeneous group of workers. The authors allow the duration dependence to shift in every four-week interval. To do so, we split the data into four-weeks intervals and consider the intervals as stand-alone explanatory as well as interction variable.
#' 
#' Splitting data into intervals
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create variable indicating if any policy change (PBD, RR, or PBD+RR) apply to observation
unemp %>%
  mutate(all = tr*(t39+t52)) -> unemp

# Define the breaks for the four-week interval
breaks <- seq(from=3,to=59, by=4)

# Create labels for the intervals 
labels <- paste("(", c(0,breaks), ",", c(breaks,104), "]",sep="")

# Create split data set at specified times
usplit <- survSplit(Surv(dur104, uncc)~., # uncc is the censoring indicator
                    data = unemp,
                    cut = breaks, # cut data at predfined breaks
                    end = "time", # defined name of event time variable
                    start = "start", # define name of start time variable
                    event = "nonjob_exit", # name of censoring indicator
                    episode = "interval" # define name of episode variable
                    )

# Create new variable exposure indicating the duration and interval
usplit %>%
  mutate(duration = time - start, # create duration variable
        interval = factor(interval+1, # given intervals are complete intervals, i.e. add one to obtain the started interval
                          labels = labels)) ->  usplit

#' 
#' Estimating Piecewise Proportional Hazard Model (Treatment Effects) using the above-mentioned Cox's PH model.
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
pwe <- coxph(Surv(duration, nonjob_exit)
             ~ age 
             + single
             + married
             + divorced
             + sfrau
             + f_marr 
             + f_single 
             + f_divor
             + lehre
             + med_educ 
             + hi_educ
             + bc 
             + lwage 
             + ten72 
             + pnon_10 
             + seasonal 
             + manuf 
             + y1988 
             + y1989 
             + y1990 
             + y1991 
             + q2 
             + q3 
             + q4
             + interval # duration dependence
             + interval * tr 
             + interval * t39 
             + interval * t52 
             + interval * all 
             + interval * after
             # Treatment effects (alternatively use the seperate variables)
             + interval * tr * after
             + interval * t39 * after
             + interval * t52 * after 
             + interval * t39 * tr * after
             + interval * t52 * tr * after,
             data = usplit
              )

#' 
#' 
#' Results
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
stargazer(pwe, 
          dep.var.caption = "", dep.var.labels = "",
          omit.table.layout = "n", star.cutoffs = NA,
          keep.stat = c("n", "ll"), no.space = TRUE,
          header = FALSE,
          title = "PWE Model Results (Treatment Effects)", type = "text"
)

## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

