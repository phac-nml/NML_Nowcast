##########################################################
# Load in packages & source functions
##########################################################
#Package list
packages <- c("plyr", "tidyverse", "nnet", "emmeans", "HelpersMG", "reshape2", "ggplot2", "ggrepel", "ggeffects", "ggpubr", "rjson","dplyr",
              "rstudioapi", "plotly", "ggthemes", "ggrepel", "scales", "plyr", "rpivotTable", "data.table", "DT", "htmltools")

#Load packages
tmp<-lapply(packages, require, character.only=T); rm(tmp)

#Set dir to current file folder - replaces here package & project dependency 
setwd(dirname(getSourceEditorContext()$path) %>%
        sub(., pattern="scripts", replacement=""))

##########################################################
# Load functions
##########################################################
source("scripts/fcn_NML_nowcast_core_webversion.r")

##########################################################
# Generate/update lineage table to stay up-to-date, saved in "data" folder 
##########################################################
# create output data to save the lineage hierarchy associations files generated in the next line
dir.create(file.path("data"), showWarnings = FALSE) # where "data" is the folder where the input data & lineage hierarchy is sourced

linhierfile <- genLinHierTable(dir2save="data") # return the location of lineage hierarchy file    

##########################################################
# Setting day 0
# If resetting to some other date, make sure it's a Sunday
##########################################################
week0day1 <- as.Date("2020-01-05") 

##########################################################
# Set data path
##########################################################
# Set file name of the lineage calls
# Each row is a sequence entry,
# columns are "date" for collection day, and "pango_lineage" for the lineage call from Pangolin
# A toy dataset is currently included in data/ - change as necessary
timeFile <- "data/exampleInput.csv"

##########################################################
# Load in data & pre-process
##########################################################
timeDat.s <- getSurveilData(timeFile, week0day1=week0day1)


### Choosing the end time to model
# No matter which ever day of the week specified here, the 'actual_time_end' will be the Saturday of the epiweek this time_end is in.
# Backtrack 2 epiweeks from recorded data & set it as model end date 
time_end <- as.Date((max(timeDat.s$week, na.rm=T) - 2)*7 + week0day1)
#time_end <- as.Date("2022-12-11") #Example to set a definite end date

##########################################################
# Run Nowcast model to generate figures and tables
##########################################################

# create output folder to save modelling results
dir.create(file.path("output"), showWarnings = FALSE) # where "output" is the name of the folder where modelling results are placed into


### Option 1: Run Nowcast with automatically selected lineages to model
### toestim="standard50": with automatically selected lineages with at least 50 sequences within any of the 8 modelled weeks + the 2 projected weeks
### ref2use="BA.5.2": Using the "BA.5.2" category as reference
standard <- runNowcast(timeDat.s = timeDat.s, 
                       week0day1 = week0day1, 
                       time_end = time_end, 
                       toestim="standard50",
                       ref2use="BA.5.2",
                       model_weeks=8, 
                       nweekPred=2, 
                       loi=NULL, 
                       output_dir="output",
                       linhierfile = linhierfile)


### Option 2: Run Nowcast with specified lineages with BA.5 lineage as reference
### toestim="loi": using specified lineages in the loi param
### ref2use="BA.5": find and use the lineage with the highest proportion as reference
loi1 <- runNowcast(timeDat.s = timeDat.s, 
                   week0day1 = week0day1, 
                   time_end = time_end, 
                   toestim="loi",
                   loi=c("BA.5", "BA.4", "BA.4.6","BF.7","BF.10","BQ.1","BQ.1.1"), 
                   ref2use="BA.5",
                   model_weeks=8, 
                   nweekPred=2, 
                   output_dir="output",
                   linhierfile = linhierfile)


### Option 2.1: Run Nowcast with specified lineages using specified reference 
### toestim="loi": using specified lineages in the loi param
### ref2use="BQ.1": Choosing own reference (this has to be one of the lineages modelled)
loi2 <- runNowcast(timeDat.s = timeDat.s, 
                   week0day1 = week0day1, 
                   time_end = time_end, 
                   toestim="loi",
                   loi=c("BA.5", "BA.4", "BA.4.6","BF.7","BF.10","BQ.1","BQ.1.1"), 
                   ref2use="BQ.1",
                   model_weeks=8, 
                   nweekPred=2, 
                   output_dir="output",
                   linhierfile = linhierfile)


### Option 3: Run Nowcast with automatically selected lineages & re-group these within specified parental supergroups to model (default in web version)
### toestim="webteam50": Uses 2 filters, in the first filter is "standard50" but it uses a table of full counts instead of partial to determine the top 5 lineages (via top counts) within any lineage group of filter#2.
# In the second filter, the selected lineages are re-assessed to be grouped under the umbrella of a parental supergroup, in which several sub-lineages share direct ancestry. If a lineage does not fall within any supergroup then it is designated its own group (these are shown under the “Other” supergroup in the website). Under any given supergroup, an assessment is made to collapse any lineages that are too small and similar to their related sister lineages in order to obtain a total of five subgroups (i.e. only ≤ 5 sub-lineages with the highest counts over the modelling period remain under a parental supergroup, the excised are re-grouped into sub-lineages sharing direct ancestry),
### ref2use="top": find and use the lineage with the highest proportion as reference
### supergroups=c("x","y","x"): these are the primary parental groups which descendant sublineages will fall under if they contain a count of ≥ 50 in any of the 8 modelling (model_weeks) + 2 forecasted (nweekPred) weeks specfied below in runNowcast()
webteam <- runNowcast(timeDat.s = timeDat.s, 
                      week0day1 = week0day1, 
                      time_end = time_end, 
                      toestim="webteam50",
                      ref2use="top",
                      supergroups=c("BQ.1", "BQ.1.1", "BQ.*", "BA.5","BA.2.75", "BA.4.6", "BA.5.2", "XBB"),
                      model_weeks=8, 
                      nweekPred=2, 
                      output_dir="output",
                      linhierfile = linhierfile)

##########################################################
# Render basic Nowcast report
# Note - users should customize and select visualizations to tailor the reports to suit own needs
  # The user must specify a variable of saved outputs from the modelling above (i.e. standard, loi1, loi2, webteam)
    # This input is placed in the 'modelRMD <- VARIABLE' below, where variable is the saved model outputs
##########################################################

modelRMD <- standard #Set stored model variable here to generate a report (i.e. standard, loi1, loi2, webteam)

#Render markdown report
rmarkdown::render(input="scripts/Reporting_script.Rmd",
                  output_format="html_document",
                  output_dir="output",
                  output_file=paste0("NMLNowcast_BarebonesReport_",modelRMD$actual_time_end,"_",modelRMD$params$toestim,"_",modelRMD$params$model_weeks,"_ref",modelRMD$params$ref2use,".html"), #Set your own output file name here
                  runtime="static",
                  params= list(output_dir="output",
                               data=modelRMD
                               ))

