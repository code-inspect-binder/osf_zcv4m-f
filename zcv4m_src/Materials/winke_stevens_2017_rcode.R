###################################################
### winke_stevens_rcode.R
### Created by Tim Winke on 7 Feb 2011 (tim.winke@gmail.com)
### Finalized by Jeffrey Stevens on 26 Sep 2017 (jeffrey.r.stevens@gmail.com)
### Summary: This script calculates descriptive statistics and generates figures
###   for the analysis of the cooperative memory experiment for the following article:
###   Winke, T. & Stevens, J.R. (2017). Is cooperative memory special? 
###    The role of costly errors, context, and social network size when remembering 
###    cooperative actions. Frontiers in Robotics and AI. doi:10.3389/frobt.2017.00052.
### Instructions: Place this file and the data files (winke_stevens_2017_data.csv) 
###   in the same directory.  Create a folder called "figures". In R set the 
###	  working directory to this directory.  Type 
### 	> source("winke_stevens_2017_rcode.R")
### 	This will run the script, adding all of the calculated variables to the
### 	workspace and saving PDF versions of the figures in the figures directory.
### Description of the data columns for stevens_past_interacts.txt:
### experiment - experiment number (1=original, 2=social network size replication, 0=Stevens et al. 2011)
### subject - name of subject
###	date - date of session
###	age - age of subject
###	gender - gender of participant
### payoff_scheme - payoff scheme (Standard or Costly)
### context - memory context (Cooperation or Newspaper)
###	trial - trial number
### intervening - number of intervening trials between interactions with a partner
###	accuracy - binary value for whether participant was correct on this trial (0 = incorrect, 1 = correct)
###	partner_action - partner's action chosen in previous trial
###	total_payment - total amount of money received (in euros)
### contacts - number of social contacts for each participant
### License: This script is released under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
###   International license (CC BY-NC-SA 4.0). You may share and adapt this content with attribution, for 
###   non-commercial purposes if you ShareAlike (distribute any contributions under the same license).
###################################################

##############################
## Load libraries, clear variables, and define functions
##############################

# Load libraries
library(BayesFactor) # required to calculate Bayes factors
library(cowplot)  # required for plot_grid
library(ggplot2)  # required for plotting
library(lme4) # required for GLMMs
library(plyr)  # required for ddply
library(dplyr)  # required for ddply
library(papaja) # required for formatting
library(tidyr)  # required for spread

# Clear variables
rm(list=ls())  # clear all variables

################
# Convert BIC values to Bayes factor
################

bic_bf <- function(null, alternative) {
  new_bf <- exp((null - alternative) / 2) # convert BICs to Bayes factor
  names(new_bf) <- NULL   # remove BIC label
  return(new_bf)  # return Bayes factor
}

################
# Calculate statistical leverage to determine outliers
################

leverage <- function(x){
  1/length(x) + (x - mean(x)) ^ 2 / sum((x-mean(x)) ^ 2)
} 

################
# Calculate standard error of the mean
################

sem <- function(data) {
  sd(data, na.rm = TRUE) / sqrt(length(data))
}

##############################
## Load data
##############################

all_data <- read.csv("winke_stevens_2017_data.csv") # input all data
all_data$accuracy <- as.numeric(as.character(all_data$accuracy))  # convert to numeric

expt1_data <- subset(all_data, experiment == 1)   # subset experiment 1 data
expt1_data$subject <- as.factor(expt1_data$subject)		# convert subject to factor
expt1_data$partner_coop <- factor(expt1_data$partner_action, labels = c(1, 0))  # create partner_coop from partner_action
expt1_data$partner_coop <- as.numeric(as.character(expt1_data$partner_coop))  # convert to numeric

expt2_data <- subset(all_data, experiment == 2)   # subset experiment 2 data

stevens_data <- subset(all_data, experiment == 0)   # subset Stevens et al. (2011) data

##############################
## Analyze distribution of partner's action
##############################

# Aggregate partner cooperation by distinguishing b/w cooperation and newspaper context
cooperating <- subset(expt1_data, context == "Cooperation", drop = TRUE)  # subset cooperate/defect condition
cooperating_pc <- cooperating %>% group_by(subject, context) %>% summarize(cooplc_means = mean(partner_coop))  # calculate mean partner cooperation proportion for each subject

reading <- subset (expt1_data, context=="Neutral", drop=TRUE) # subset read/not read condition
reading_pc <- reading %>% group_by(subject, context) %>% summarize(readlc_means = mean(partner_coop))  # calculate mean partner cooperation proportion for each subject

all_pc <- expt1_data %>% group_by(subject, context) %>% summarize(alllc_means = mean(partner_coop))  # calculate mean partner cooperation proportion for each subject and context
all_pc$included <- ifelse(all_pc$alllc_means > 0.4 & all_pc$alllc_means < 0.6, "Included", "Not included")  # label subjects as included or not included depending on the proportion of cooperation for their partners

# Add partner cooperation means to expt1_data
expt1_data <- merge(expt1_data, all_pc, by = c("subject", "context"), all.x = TRUE)
expt1_data <- merge(expt1_data, cooperating_pc, by = "subject", all.x = TRUE, suffixes = c("", "y"))
expt1_data <- merge(expt1_data, reading_pc, by = "subject", all.x = TRUE, suffixes = c("", "y2"))
expt1_data <- expt1_data[-c(18, 20)]  # remove extra columns

# Subset data that includes only partner cooperation means within 0.4 - 0.6
bound_data <- subset(expt1_data, expt1_data$cooplc_means >= 0.4 & expt1_data$cooplc_means <= 0.6 & expt1_data$readlc_means >= 0.4 & expt1_data$readlc_means <= 0.6)

# Remove participants with uncorrected txt-files to have balanced groups
bound_data <- subset(bound_data, bound_data$subject != 12 & bound_data$subject != 32)
bound_data$subject <- bound_data$subject[, drop = TRUE]   # drop unused levels of factor subject

num_subjects <- length(unique(expt1_data$subject))  # calculate total number of subjects
num_subjects_clean <- length(unique(bound_data$subject))  # calculate number of subjects in bounded data set
num_subjects_removed <- num_subjects - num_subjects_clean   # calculate number of subjects removed

# Aggregate partner cooperation for bounded data set
bound_pc <- aggregate(partner_coop ~ subject * context, data = bound_data, mean)
names(bound_pc) = c("subject", "context",  "lastchoice_means")

# Plot histogram of chances that partner's choice was positive for all data
coop_hist_ggplot <- ggplot(all_pc, aes(x = alllc_means * 100)) +
  geom_histogram(aes(fill = included), bins = 50) +  # plot histogram
  scale_fill_manual(values = c("black", "grey50"), name="", label=c("Included", "Not included")) +  # color values inside and outside of 0.4-0.6 differently
  labs(x = "Percent partner positive actions", y = "Number of participants") + # label axes
  theme_classic() + # use classic theme
  theme(axis.title=element_text(size=45), axis.text=element_text(size=30), legend.text=element_text(size=30), legend.position = c(0.25, 0.9), legend.key.size = unit(2.5, 'lines'))
png(file = "figures/partner_action_histogram.png", width = 1200, height = 750) # open device
plot(coop_hist_ggplot) # plot figure
dev.off() # close device

##############################
# Demographic info
##############################

##############
# Experiment 1
##############

expt1_data_subj <- bound_data %>% group_by(subject, age, gender, total_payment) %>% summarize(accuracy = mean(accuracy * 100), contacts = mean(contacts)) # calculate mean accuracy by subject, age, gender, and total_payment
gender1 <- table(expt1_data_subj$gender) # calculate gender frequency
age_mean1 <- mean(expt1_data_subj$age)   # calculate age mean
age_sd1 <- sd(expt1_data_subj$age)       # calculate age sd
age_min1 <- min(expt1_data_subj$age)     # calculate age min
age_max1 <- max(expt1_data_subj$age)     # calculate age max

##############
# Experiment 2
##############

expt2_data_subj <- expt2_data %>% group_by(subject, age, gender, total_payment) %>% summarize(accuracy = mean(accuracy), contacts = mean(contacts))   # calculate mean accuracy by subject, age, gender, and total_payment
num_subjects2 <- length(unique(expt2_data_subj$subject))  # calculate number of subjects
gender2 <- table(expt2_data_subj$gender) # calculate gender frequency
age_mean2 <- mean(expt2_data_subj$age)   # calculate age mean
age_sd2 <- sd(expt2_data_subj$age)       # calculate age sd
age_min2 <- min(expt2_data_subj$age)     # calculate age min
age_max2 <- max(expt2_data_subj$age)     # calculate age max

##############################
# Compare current data to Stevens et al. (2011)
##############################
## Prepare data
stevens_data_subj <- stevens_data %>% group_by(subject) %>% summarize(accuracy = mean(accuracy * 100, na.rm = TRUE))   # calculate mean accuracy per subject
stevens_mean <- mean(stevens_data_subj$accuracy) # calculate overall mean accuracy
stevens_sd <- sd(stevens_data_subj$accuracy) # calculate overall sd of accuracy

expt1_compare_data <- subset(bound_data, payoff_scheme == "Standard" & context == "Cooperation") # subset standard payoff scheme and cooperation context to compare to other data sets
expt1_compare_data_subj <- expt1_compare_data %>% group_by(subject) %>% summarize(accuracy = mean(accuracy * 100, na.rm = TRUE)) # calculate mean accuracy per subject
expt1_compare_mean <- mean(expt1_compare_data_subj$accuracy) # calculate overall mean accuracy
expt1_compare_sd <- sd(expt1_compare_data_subj$accuracy) # calculate overall sd of accuracy

expt2_compare_data_subj <- expt2_data %>% group_by(subject) %>% summarize(accuracy = mean(accuracy * 100, na.rm = TRUE)) # calculate mean accuracy per subject
expt2_compare_mean <- mean(expt2_compare_data_subj$accuracy) # calculate overall mean accuracy
expt2_compare_sd <- sd(expt2_compare_data_subj$accuracy) # calculate overall sd of accuracy

## Compare mean accuracy
# Experiment 1 and Stevens et al. (2011)
compare1_ttest <- t.test(expt1_compare_data_subj$accuracy, stevens_data_subj$accuracy) # conduct t-test
compare1_diff <- compare1_ttest$estimate[1] - compare1_ttest$estimate[2] # calculate mean difference
compare1_diff_ci <- compare1_ttest$conf.int[2] - compare1_diff  # calculate 95% CI
compare1_diff_d <- compare1_diff / max(sd(expt1_compare_data_subj$accuracy), sd(stevens_data_subj$accuracy))  # calculate Cohen's d
compare1_ttestbf <- ttestBF(stevens_data_subj$accuracy, expt1_compare_data_subj$accuracy)  # calculate Bayes factor
compare1_ttestbf_bf <- extractBF(compare1_ttestbf)$bf # extract Bayes factor

# Experiment 2 and Stevens et al. (2011)
compare2_ttest <- t.test(expt2_compare_data_subj$accuracy, stevens_data_subj$accuracy) # conduct t-test
compare2_diff <- compare2_ttest$estimate[2] - compare2_ttest$estimate[1] # calculate mean difference
compare2_diff_ci <- compare2_diff - compare2_ttest$conf.int[2]  # calculate 95% CI
compare2_diff_d <- compare2_diff / max(sd(expt2_compare_data_subj$accuracy), sd(stevens_data_subj$accuracy))  # calculate Cohen's d
compare2_ttestbf <- ttestBF(stevens_data_subj$accuracy, expt2_compare_data_subj$accuracy)  # calculate Bayes factor
compare2_ttestbf_bf <- extractBF(compare2_ttestbf)$bf # extract Bayes factor

# Experiments 1 and 2
compare2b_ttest <- t.test(expt2_compare_data_subj$accuracy, expt1_compare_data_subj$accuracy) # conduct t-test
compare2b_diff <- compare2b_ttest$estimate[2] - compare2b_ttest$estimate[1] # calculate mean difference
compare2b_diff_ci <- compare2b_diff - compare2b_ttest$conf.int[2]  # calculate 95% CI
compare2b_diff_d <- compare2b_diff / max(sd(expt2_compare_data_subj$accuracy), sd(expt1_compare_data_subj$accuracy))  # calculate Cohen's d
compare2b_ttestbf <- ttestBF(expt1_compare_data_subj$accuracy, expt2_compare_data_subj$accuracy)  # calculate Bayes factor
compare2b_ttestbf_bf <- extractBF(compare2b_ttestbf)$bf # extract Bayes factor

# Compare intervening events
stevens_intervene_subj <- stevens_data %>% group_by(subject, intervening) %>% summarize(accuracy = mean(accuracy, na.rm = TRUE))  # calculate mean per subject and number of intervening events
expt1_intervene_subj <- subset(expt1_compare_data, !is.na(intervening)) %>% group_by(subject, intervening) %>% summarize(accuracy = mean(accuracy, na.rm = TRUE), sem = sd(accuracy, na.rm = TRUE) / sqrt(length(accuracy)))  # calculate mean per subject and number of intervening events
expt2_intervene_subj <- subset(expt2_data, !is.na(intervening)) %>% group_by(subject, intervening) %>% summarize(accuracy = mean(accuracy, na.rm = TRUE), sem = sd(accuracy, na.rm = TRUE) / sqrt(length(accuracy)))  # calculate mean per subject and number of intervening events

stevens_intervene <- stevens_intervene_subj %>% group_by(intervening) %>% summarize(accuracy = mean(accuracy * 100))  # calculate mean per number of intervening events
stevens_intervene_sem <- stevens_intervene_subj %>% group_by(intervening) %>% summarize(sem = sem(accuracy))  # calculate sem per number of intervening events
stevens_intervene$sem <- stevens_intervene_sem$sem  # add sem column

expt1_intervene <- expt1_intervene_subj %>% group_by(intervening) %>% summarize(accuracy = mean(accuracy * 100))  # calculate mean per number of intervening events
expt1_intervene_sem <- expt1_intervene_subj %>% group_by(intervening) %>% summarize(sem = sem(accuracy))  # calculate sem per number of intervening events
expt1_intervene$sem <- expt1_intervene_sem$sem  # add sem column

expt2_intervene <- expt2_intervene_subj %>% group_by(intervening) %>% summarize(accuracy = mean(accuracy * 100))  # calculate mean per number of intervening events
expt2_intervene_sem <- expt2_intervene_subj %>% group_by(intervening) %>% summarize(sem = sem(accuracy))  # calculate sem per number of intervening events
expt2_intervene$sem <- expt2_intervene_sem$sem  # add sem column

intervene_data <- data.frame(experiment = c(rep("Stevens", length(stevens_intervene$sem)), rep(1, length(expt1_intervene$sem)), rep(2, length(expt2_intervene$sem))), rbind.fill(stevens_intervene, expt1_intervene, expt2_intervene)) # create date frame of both data sets
intervene_data$lower <- intervene_data$accuracy - intervene_data$sem  # calculate lower sem
intervene_data$upper <- intervene_data$accuracy + intervene_data$sem  # calculate upper sem

# Create intervening events plot
intervene_ggplot <- ggplot(subset(intervene_data, intervening < 17), aes(x = intervening, y = accuracy, color = experiment)) +
  geom_line(size = 1, aes(linetype = experiment)) + # add lines for experiments
  geom_point(size = 5, aes(shape = experiment)) + # add points for experiments
  labs(x = "Number of intervening events",y = "Memory accuracy (%)") + # label axes
  scale_linetype_manual(values=c(1, 2, 3), name="", labels=c("Experiment 1", "Experiment 2", "Stevens et al. (2011)")) +  # define linetypes
  scale_color_manual(name="", values = c("black", "grey30", "grey60"), labels = c("Experiment 1", "Experiment 2", "Stevens et al. (2011)")) +  # define line colors
  scale_shape_discrete(name="", label=c("Experiment 1", "Experiment 2", "Stevens et al. (2011)")) +  # define point shapes
  scale_fill_discrete(name="", label=c("Experiment 1", "Experiment 2", "Stevens et al. (2011)")) +  # define point colors
  theme_classic() + # use classic theme
  theme(axis.title=element_text(size=45), axis.text=element_text(size=30), legend.text=element_text(size=30), legend.position = c(0.75, 0.9), legend.key.size = unit(2.75, 'lines'))
png(file = "figures/intervene.png", width = 750, height = 750) # open device
plot(intervene_ggplot) # plot figure
dev.off() # close device

##############################
# Analyze accuracy as a function of payoff scheme (Standard or Costly), context (Cooperation or Newspaper), and partner action (Cooperate or Defect)
##############################

# Conduct binomial GLMM of payoff scheme * partner action + context for memory accuracy
accuracy_glmer_full <- glmer(accuracy ~ payoff_scheme * partner_action * context + (1 | subject), bound_data, family = binomial(link = "logit"))   # calculate GLMM of full model
# accuracy_glmer_full_ci <- confint(accuracy_glmer_full, method = "profile")    # calculate profile likelihood CIs ****takes a long time to run, so values hard coded below
accuracy_glmer_full_ci <- data.frame(fact = c("payoff_scheme", "partner_action", "context", "payoff_action", "payoff_context", "action_context", "payoff_action_context"), lower = c(-0.36, 0.28, 0.04, -0.47, -0.45, -0.58, -0.57), upper = c(0.52, 0.73, 0.47, 0.17, 0.15, 0.07, 0.33))    # profile likelihood CIs

# Conduct binomial GLMMs for memory accuracy to calculate BIC values to transform to Bayes factors
accuracy_glmer_null <- summary(glmer(accuracy ~ (1 | subject), bound_data, family = binomial(link = "logit")))    # calculate GLMM for null model
accuracy_glmer_payoff <- summary(glmer(accuracy ~ payoff_scheme  + (1 | subject), bound_data, family = binomial(link = "logit")))    # calculate GLMM for payoff scheme
accuracy_glmer_context <- summary(glmer(accuracy ~ context + (1 | subject), bound_data, family = binomial(link = "logit")))    # calculate GLMM for context
accuracy_glmer_action <- summary(glmer(accuracy ~ partner_action + (1 | subject), bound_data, family = binomial(link = "logit")))    # calculate GLMM for partner action
accuracy_glmer_payoff_action <- summary(glmer(accuracy ~ payoff_scheme + partner_action + (1 | subject), bound_data, family = binomial(link = "logit")))    # calculate GLMM for payoff_scheme + partner_action
accuracy_glmer_payoff_action_inter <- summary(glmer(accuracy ~ payoff_scheme * partner_action + (1 | subject), bound_data, family = binomial(link = "logit")))    # calculate GLMM for payoff_scheme * partner_action
accuracy_glmer_action_context <- summary(glmer(accuracy ~ partner_action + context + (1 | subject), bound_data, family = binomial(link = "logit")))    # calculate GLMM for partner_action + context
accuracy_glmer_action_context_inter <- summary(glmer(accuracy ~ partner_action * context + (1 | subject), bound_data, family = binomial(link = "logit")))    # calculate GLMM for partner_action * context
accuracy_glmer_payoff_action_context <- summary(glmer(accuracy ~ payoff_scheme + partner_action + context + (1 | subject), bound_data, family = binomial(link = "logit")))   # calculate GLMM for payoff_scheme * partner_action * context

# Transform BICs to Bayes factors
accuracy_glmer_payoff_bf <- bic_bf(accuracy_glmer_null$AICtab[2], accuracy_glmer_payoff$AICtab[2])    # calculate Bayes factor for payoff scheme
accuracy_glmer_context_bf <- bic_bf(accuracy_glmer_null$AICtab[2], accuracy_glmer_context$AICtab[2])    # calculate Bayes factor for context
accuracy_glmer_action_bf <- bic_bf(accuracy_glmer_null$AICtab[2], accuracy_glmer_action$AICtab[2])    # calculate Bayes factor for partner action
accuracy_glmer_payoff_action_bf <- bic_bf(accuracy_glmer_payoff_action$AICtab[2], accuracy_glmer_payoff_action_inter$AICtab[2])    # calculate Bayes factor for payoff scheme * partner_action interaction
accuracy_glmer_action_context_bf <- bic_bf(accuracy_glmer_action_context$AICtab[2], accuracy_glmer_action_context_inter$AICtab[2])    # calculate Bayes factor for payoff scheme * partner action interaction
accuracy_glmer_payoff_action_context_bf <- bic_bf(accuracy_glmer_payoff_action_context$AICtab[2], summary(accuracy_glmer_full)$AICtab[2])    # calculate Bayes factor for payoff scheme * partner action interaction

##############################
# Analyze accuracy as a function of payoff scheme (Standard or Costly) and partner action (Cooperate or Defect) for last 40 trials
##############################

# Subset data from last 40 trials
last_data <- subset(bound_data, context == "Cooperation" & ((trial > 59 & trial < 101) | (trial > 159)))

# Conduct binomial GLMM of payoff scheme * partner action + context for memory accuracy
last_accuracy_glmer_full <- glmer(accuracy ~ payoff_scheme * partner_action + (1 | subject), last_data, family = binomial(link = "logit"))   # calculate GLMM of full model
# last_accuracy_glmer_full_ci <- confint(last_accuracy_glmer_full, method = "profile")    # calculate profile likelihood CIs ****takes a long time to run, so values hard coded below
last_accuracy_glmer_full_ci <- data.frame(fact = c("payoff_scheme", "partner_action", "payoff_action"), lower = c(-0.6, 0.32, -0.81), upper = c(0.67, 1.07, 0.22))    # profile likelihood CIs

# Conduct binomial GLMMs for memory accuracy to calculate BIC values to transform to Bayes factors
last_accuracy_glmer_null <- summary(glmer(accuracy ~ (1 | subject), last_data, family = binomial(link = "logit")))    # calculate GLMM for null model
last_accuracy_glmer_payoff <- summary(glmer(accuracy ~ payoff_scheme  + (1 | subject), last_data, family = binomial(link = "logit")))    # calculate GLMM for payoff scheme
last_accuracy_glmer_action <- summary(glmer(accuracy ~ partner_action + (1 | subject), last_data, family = binomial(link = "logit")))    # calculate GLMM for partner action
last_accuracy_glmer_payoff_action <- summary(glmer(accuracy ~ payoff_scheme + partner_action + (1 | subject), last_data, family = binomial(link = "logit")))    # calculate GLMM for payoff_scheme + partner_action
last_accuracy_glmer_payoff_action_inter <- summary(glmer(accuracy ~ payoff_scheme * partner_action + (1 | subject), last_data, family = binomial(link = "logit")))    # calculate GLMM for payoff_scheme * partner_action

# Transform BICs to Bayes factors
last_accuracy_glmer_payoff_bf <- bic_bf(last_accuracy_glmer_null$AICtab[2], last_accuracy_glmer_payoff$AICtab[2])    # calculate Bayes factor for payoff scheme
last_accuracy_glmer_action_bf <- bic_bf(last_accuracy_glmer_null$AICtab[2], last_accuracy_glmer_action$AICtab[2])    # calculate Bayes factor for partner action
last_accuracy_glmer_payoff_action_bf <- bic_bf(last_accuracy_glmer_payoff_action$AICtab[2], last_accuracy_glmer_payoff_action_inter$AICtab[2])    # calculate Bayes factor for payoff scheme * partner_action interaction

##############################
# Plot memory accuracy as a function of payoff scheme, context, and partner action 
##############################

## Payoff scheme
# Prepare data
payoff_scheme_subj <- bound_data %>% group_by(payoff_scheme, subject) %>% summarize(accuracy = mean(accuracy * 100)) # calculate mean accuracy per payoff scheme and subject
payoff_scheme_subj$payoff_scheme <- factor(payoff_scheme_subj$payoff_scheme, levels = c("Standard", "Costly"))  # reorder payoff matrices
payoff_scheme <- payoff_scheme_subj %>% group_by(payoff_scheme) %>% summarize(accuracy = mean(accuracy)) # calculate mean accuracy per payoff scheme
payoff_scheme$ci <- c(ci(subset(payoff_scheme_subj, payoff_scheme == "Costly")$accuracy), ci(subset(payoff_scheme_subj, payoff_scheme == "Standard")$accuracy)) # calculate 95% CIs for accuracy per payoff scheme
payoff_scheme$lower <- payoff_scheme$accuracy - payoff_scheme$ci  # create lower CI
payoff_scheme$upper <- payoff_scheme$accuracy + payoff_scheme$ci  # create upper CI

# Plot data
payoff_scheme_ggplot <- ggplot(payoff_scheme_subj, aes(payoff_scheme, accuracy)) +
  geom_boxplot(color = "grey50", coef = 5) + # create boxplot with whiskers as full range
  geom_point(data = payoff_scheme, aes(payoff_scheme, accuracy), size = 6) +  # create points for means
  geom_linerange(data = payoff_scheme, aes(ymin = lower, ymax = upper), size = 2) + # create lines for 95% CIs
  labs(x = "Payoff scheme",y = "Memory accuracy (%)") + # label axes
  lims(y = c(47, 100.5)) +  # scale y-axis
  theme_classic() + # use classic theme
  theme(axis.title=element_text(size=45), axis.text=element_text(size=30), plot.margin = grid::unit(c(0, 0, 0.5, 1.5), "lines"))  # adjust label font sizes

## Context
# Prepare data
context_subj <- bound_data %>% group_by(context, subject) %>% summarize(accuracy = mean(accuracy * 100)) # calculate mean accuracy per payoff scheme and subject
context <- context_subj %>% group_by(context) %>% summarize(accuracy = mean(accuracy)) # calculate mean accuracy per payoff scheme
context$ci <- wsci(context_subj, dv = "accuracy", id = "subject", factors = "context")$accuracy # calculate within-subjects 95% CIs for accuracy per payoff scheme
context$lower <- context$accuracy - context$ci  # create lower CI
context$upper <- context$accuracy + context$ci  # create upper CI

# Plot data
context_ggplot <- ggplot(context_subj, aes(context, accuracy)) +
  geom_boxplot(color = "grey50", coef = 5) + # create boxplot with whiskers as full range
  geom_point(data = context, aes(context, accuracy), size = 6) +  # create points for means
  geom_linerange(data = context, aes(ymin = lower, ymax = upper), size = 2) + # create lines for 95% CIs
  labs(x = "Context", y = "Memory accuracy (%)") + # label axes
  lims(y = c(47, 100.5)) +  # scale y-axis
  theme_classic() + # use classic theme
  theme(axis.title=element_text(size=45), axis.text=element_text(size=30), plot.margin = grid::unit(c(0, 0, 0.5, 2.5), "lines"))  # adjust label font sizes

# Generate figure file
png(file = "figures/payoff_context.png", width = 1200, height = 750) # open device
plot(plot_grid(payoff_scheme_ggplot, context_ggplot, ncol = 2, labels = c("(a)", "(b)"), label_size = 40, hjust = 0.05, vjust = 1))
dev.off() # close device

## Partner action
# Prepare data
partner_action_subj <- bound_data %>% group_by(partner_action, subject) %>% summarize(accuracy = mean(accuracy * 100)) # calculate mean accuracy per payoff scheme
partner_action <- partner_action_subj %>% group_by(partner_action) %>% summarize(accuracy = mean(accuracy), sd = sd(accuracy, na.rm=TRUE)) # calculate mean accuracy per payoff scheme
partner_action$ci <- wsci(partner_action_subj, dv = "accuracy", id = "subject", factors = c("partner_action"))$accuracy # calculate within-subjects 95% CIs for accuracy per payoff scheme
partner_action$lower <- partner_action$accuracy - partner_action$ci  # create lower CI
partner_action$upper <- partner_action$accuracy + partner_action$ci  # create upper CI

partner_payoff_subj <- bound_data %>% group_by(payoff_scheme, partner_action, subject) %>% summarize(accuracy = mean(accuracy * 100)) # calculate mean accuracy per payoff scheme and subject
partner_payoff_subj$payoff_scheme <- factor(partner_payoff_subj$payoff_scheme, levels = c("Standard", "Costly"))  # reorder payoff matrices
partner_payoff <- partner_payoff_subj %>% group_by(payoff_scheme, partner_action) %>% summarize(accuracy = mean(accuracy)) # calculate mean accuracy per payoff scheme
partner_payoff$ci <- wsci(partner_payoff_subj, dv = "accuracy", id = "subject", factors = c("payoff_scheme", "partner_action"))$accuracy # calculate within-subjects 95% CIs for accuracy per payoff scheme
partner_payoff$lower <- partner_payoff$accuracy - partner_payoff$ci  # create lower CI
partner_payoff$upper <- partner_payoff$accuracy + partner_payoff$ci  # create upper CI

# Plot data
partner_payoff_ggplot <- ggplot(partner_payoff_subj, aes(x = partner_action, y = accuracy)) +
  facet_wrap(~payoff_scheme) +
  geom_boxplot(color = "grey50", coef = 5) + # create boxplot with whiskers as full range
  geom_point(data = partner_payoff, aes(partner_action, accuracy), size = 6) +  # create points for means
  geom_linerange(data = partner_payoff, aes(ymin = lower, ymax = upper), size = 2) + # create lines for 95% CIs
  labs(x = "Partner action",y = "Memory accuracy (%)") + # label axes
  theme_classic() + # use classic theme
  theme(axis.title=element_text(size=45), axis.text=element_text(size=30), strip.text.x = element_text(size = 30), strip.background = element_rect(fill="grey80"))  # adjust label font sizes
png(file = "figures/partner_payoff.png", width = 750, height = 750) # open device
plot(partner_payoff_ggplot)  # plot figure
dev.off() # close device

##############################
# Social contact and memory
##############################

##############
# Experiment 1
##############

# Calculate correlation between mean social network size and mean memory accuracy
coop_contacts_cor <- cor.test(expt1_data_subj$accuracy, expt1_data_subj$contacts) # calculate network size/accuracy correlation
coop_contacts_bfdf <- data.frame(accuracy = expt1_data_subj$accuracy, contacts = expt1_data_subj$contacts)  # create new data frame for Bayesian analysis
coop_contacts_lmbf <- lmBF(accuracy ~ contacts, data = coop_contacts_bfdf)  # calculate Bayes regression
coop_contacts_bf <- extractBF(coop_contacts_lmbf)$bf  # extract Bayes factor

# Plot data
coop_contact_ggplot <- ggplot(expt1_data_subj, aes(contacts, accuracy)) +
  geom_point(size = 5) +  # set point size
  geom_smooth(method = "lm", color = "black") +  # create linear regression line
  labs(x = "Number of contacts",y = "Memory accuracy (%)") +  # create axis labels
  lims(y = c(40, 105)) +  # scale y-axis
  annotate("text", 1250, 104, label = paste("r =", round(coop_contacts_cor$estimate, 2)), size = 12) +  # add r value
  annotate("text", 1250, 100, label = paste("BF =", round(coop_contacts_bf, 1)), size = 12) +  # add Bayes factor
  theme(axis.title=element_text(size=45), axis.text=element_text(size=30), plot.margin = grid::unit(c(0, 2, 0, 0), "lines"))    # set font sizes and buffers

# Test for outliers in social contact data
zmod <- (0.6745*(expt1_data_subj$contacts - median(expt1_data_subj$contacts))) / median(abs(expt1_data_subj$contacts - median(expt1_data_subj$contacts))) # calculate Iglewicz & Hoaglin's (1993) modified Z score
outliers <- which(zmod > 3.5) # find index number of outliers
contacts_lev <- leverage(expt1_data_subj$contacts)  # calculate leverage of contacts
influence_h <- 2 * 2 / length(expt1_data_subj$contacts) # calculate influence of contacts
high_influence <- which(contacts_lev > influence_h) # find high influence data  ***SAME AS Z-SCORE OUTLIERS***
expt1_data_subj_nooutliers <- expt1_data_subj[-c(outliers), ] # remove high influence data points

# Calculate correlation between mean Number of contacts and mean memory accuracy
coop_contacts_cor_nooutliers <- cor.test(expt1_data_subj_nooutliers$accuracy, expt1_data_subj_nooutliers$contacts) # calculate network size/accuracy correlation
coop_contacts_nooutliers_bfdf <- data.frame(accuracy = expt1_data_subj_nooutliers$accuracy, contacts = expt1_data_subj_nooutliers$contacts)  # create new data frame for Bayesian analysis
coop_contacts_nooutliers_lmbf <- lmBF(accuracy ~ contacts, data = coop_contacts_nooutliers_bfdf)  # calculate Bayes regression
coop_contacts_nooutliers_bf <- extractBF(coop_contacts_nooutliers_lmbf)$bf  # extract Bayes factor

# Plot data
coop_contact_nooutliers_ggplot <- ggplot(expt1_data_subj_nooutliers, aes(contacts, accuracy)) +
  geom_point(size = 5) +  # set point size
  geom_smooth(method = "lm", color = "black") +  # create linear regression line
  labs(x = "Number of contacts",y = "Memory accuracy (%)") + # create axis labels
  lims(y = c(40, 105)) +  # scale y-axis
  annotate("text", 300, 104, label = paste("r =", round(coop_contacts_cor_nooutliers$estimate, 2)), size = 12) +  # add r value
  annotate("text", 300, 100, label = paste("BF =", round(coop_contacts_nooutliers_bf, 1)), size = 12) +  # add Bayes factor
  theme(axis.title=element_text(size=45), axis.text=element_text(size=30), plot.margin = grid::unit(c(0, 2, 0, 1), "lines"))  # set font sizes and buffers

##############
# Experiment 2
##############

# Test for outliers in social contact data
zmod2 <- (0.6745*(expt2_data_subj$contacts - median(expt2_data_subj$contacts))) / median(abs(expt2_data_subj$contacts - median(expt2_data_subj$contacts))) # calculate Iglewicz & Hoaglin's (1993) modified Z score (from Rosenthal & Rosnow 2008)
outliers2 <- which(zmod2 > 3.5) # find index number of outliers ***NONE***

## Calculate correlation between mean number of contacts and mean memory accuracy
coop_contacts_cor2 <- cor.test(expt2_data_subj$accuracy, expt2_data_subj$contacts)  # calculate network size/accuracy correlation
coop_contacts2_bfdf <- data.frame(accuracy = expt2_data_subj$accuracy, contacts = expt2_data_subj$contacts) # create new data frame for Bayesian analysis
coop_contacts2_lmbf <- lmBF(accuracy ~ contacts, data = coop_contacts2_bfdf)  # calculate Bayes regression
coop_contacts2_bf <- extractBF(coop_contacts2_lmbf)$bf  # extract Bayes factor

# Plot data
coop_contact2_ggplot <- ggplot(expt2_data_subj, aes(contacts, accuracy * 100)) +
  geom_point(size = 5) +  # set point size
  geom_smooth(method = "lm", color = "black") +  # create linear regression line
  labs(x = "Number of contacts",y = "Memory accuracy (%)") + # create axis labels
  lims(y = c(40, 105)) +  # scale y-axis
  annotate("text", 12.5, 104, label = paste("r =", round(coop_contacts_cor2$estimate, 2)), size = 12) +  # add r value
  annotate("text", 12.5, 100, label = paste("BF =", round(coop_contacts2_bf, 1)), size = 12) +  # add Bayes factor
  theme(axis.title=element_text(size=45), axis.text=element_text(size=30), plot.margin = grid::unit(c(0, 1, 0, 0), "lines"))   # set font sizes and buffers 

png(file = "figures/coop_contact_all.png", width = 2550, height = 750) # open device
plot(plot_grid(coop_contact_ggplot, coop_contact_nooutliers_ggplot,  coop_contact2_ggplot, ncol = 3, labels = c("(a)", "(b)", "(c)"), label_size = 40, hjust = 0, vjust = 1))  # plot figure
dev.off() # close device

