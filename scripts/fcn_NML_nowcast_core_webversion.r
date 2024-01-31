########################################################
## Functions to be called
## Written by Julie Chih-yu Chen
## Modified and cleaned up on 2024-01-31
########################################################

## This script is only dependent on the two upstream pangolin files: alias_key.json and (the first column of) lineage_notes.txt
## To create the parental lineage table for each lineage using the alias json file and the full lineage list from 
## output saved in .csv file: a list with each row being the lineage, along with its ancestors
## return file directory to the lineage table.
genLinHierTable <- function(dir2save="data"){
  
  alias_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json"
  lineages_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt"
  
  ### load alias json
  aliasData <- fromJSON(file=alias_url)
  ## skipping recombination which has multiple parents here
  aliasData.s <- sapply(aliasData,function(y){ret=ifelse(length(y)==1,y,NA);ret})
  aliasData.s <- na.omit(aliasData.s)
  
  ### load lineage list
  ## The following lines dealt w extra \t in some lines 2023-08-07
  linlist <- readLines(lineages_url)
  linlist <- data.frame(t(sapply(strsplit(linlist, "\t"), function(y) y[1:2]))) #ignoring things after the 2nd tab where it appears
  colnames(linlist)=linlist[1,]; linlist <- linlist[-1,]
  
  ## remove withdrawn ones with the *
  linlist <- linlist[-grep("^\\*",linlist[,1]),]
  ## create a new full lineage name without the aliases
  linlist[,3]=linlist$Lineage
  for(y in 3:length(aliasData.s)){ ## skipping A and B
    id <- grep(paste0("^",names(aliasData.s)[y],"\\."),linlist$Lineage) 
    linlist[id,3]=gsub(names(aliasData.s)[y],aliasData.s[y],linlist[id,3])
  }
  
  ### get ancestors in list
  fullHier2 <- lapply(linlist[,3],function(y){
    ## split the levels
    tmp <- strsplit(y,"\\.")[[1]]
    
    ## reconstruct all ancestors, including itself (no alias). 
    ancs <- sapply(1:length(tmp), function(x){paste(tmp[1:x],collapse=".")})
    
    ## output the ancestral lineage names with corresponding aliases, so we can match it up after. 
    ## * this step is dependent on the sorting of lineage_note file, from ancestor to descendants, which has been the case.
    ## In the order of child (self) to parent and ancestries
    rev(linlist[,1][linlist[,3]%in%ancs])
  } )
  names(fullHier2)=linlist[,1]
  
  ### generate the table of hierarchy file
  #####################
  ### long data format to match what Adrian had, saving the csv file and the date of extraction in .txt file
  #####################
  ancestry <- ldply(fullHier2, data.frame)
  colnames(ancestry)=c("children","lineage")
  
  directParent <- do.call(rbind,lapply(fullHier2, function(y)y[c(1:2)]))
  colnames(directParent)=c("lineage","parent")
  ancestry2 <- merge(data.frame(order=1:nrow(ancestry),ancestry),directParent,by="lineage",all.x=T, sort=F)
  ancestry2 <- ancestry2[order(ancestry2$order),c("lineage","parent","children")]
  
  write.csv(ancestry2, file = file.path(dir2save, "lineage-hierarchy.csv"),quote=F, row.names=F)
  write.table(format(Sys.time(), "%b %e %Y"), file = file.path(dir2save, "lineage-hierarchy.txt"),quote=F, row.names=F,col.names=F)
  print("Lineage Hierarchy file saved.")
  file.path(dir2save, "lineage-hierarchy.csv")
}

####################################
### read data and data cleaning
####################################
getSurveilData <- function(timeFile, week0day1){
  # Load in table
  timeDat.s <- read.table(file=timeFile, sep=",", header=T, quote="",fill=T, as.is=T)
  
  #Format dates
  timeDat.s$date <- as.Date(timeDat.s$date)
  
  #Standardize weeks to cdc NowCast "epiweek" based upon week0day1 -> Sunday standardization
  timeDat.s$week <- as.numeric(timeDat.s$date  - week0day1) %/% 7
  
  #Convert weeks to a date relative to week0day1
  timeDat.s$weekDate <- timeDat.s$week*7 + week0day1 # the corresponding Sunday
  
  currentweek <- as.numeric(Sys.Date()  - week0day1) %/% 7
  
  #Remove sequencing dates with a negative week or current week (these are errors)
  timeDat.s <- timeDat.s %>% filter(week > 0 & week < currentweek)
  #Remove all pango_lineage's which are "Unassigned, None, none, or Unknown"
  timeDat.s <- timeDat.s %>% filter(!pango_lineage %in% c("Unassigned","None","none","Unknown"))
  #Remove all "NA" pango lineage
  timeDat.s <- timeDat.s %>% filter(!is.na(pango_lineage))
  timeDat.s
}



####################################
# Function to grab sublineages for each lineage lin
####################################
getSublineages <- function(lin, linhier){
  
  p2c <- paste(linhier$lineage, linhier$children, sep="_")
  
  # get descendants
  str <- paste0(lin,"_")
  tmp <- grep(str,p2c)
  c(unique(gsub(str,"",p2c[tmp]))) #return itself and the descendants names
}

####################################
# Function to make sure sublineages did not overlap between entries of a list
# ie. Given a list with BA.1 and BA.1.1 and their descendants, make sure that BA.1.1 as well as the descendants of BA.1.1 aren't in BA.1's vector of descendants
# input is a list with the list names being the lineage names of topn or loi
####################################
lineageExclusive <- function(linlist){
  linname <- names(linlist)
  
  # !Need to make sure each of the topn (and its descendants) is not repeated in another, so removing from descendant list
  for( i in 1:length(linlist)){
    for(j in 1:length(linlist)){
      # if these are different entries
      if(i!=j){
        # check if the lineage of interest j is already in the descendants of i
        if(sum(linlist[[i]]%in%names(linlist)[j])>0){
          # if so, remove lineage j and its descendants from lineage i descendants
          linlist[[i]]=setdiff(linlist[[i]],linlist[[j]])
        }
      }
    }#end j
  }#end i
  linlist
}


####################################
# Function to generate the trend plot
# input: 
#	  dat: Weekly proportion table with actual proportion (prop) and predicted proportion (predicted) information for variants
#	  actual_time_end: the actual modelling end date
####################################
trendPlot <- function(dat, actualTimeEnd){
  dat2plot = dat[!is.na(dat$week.x),] #plotting just for the weeks and not in between
  dat2plot[,c("prop","predicted", "std.error","conf.low","conf.high")] = dat2plot[,c("prop","predicted", "std.error","conf.low","conf.high")]*100
  p <- ggplot(dat2plot, aes(x = weekDate, y = predicted, group = variant, col= variant, fill=variant)) +
    geom_line() + #fitted line
    ylab("Proportion (%)")+
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, linetype=0) + #CI
    geom_vline(xintercept=actualTimeEnd, linetype="dashed", color="darkgrey", linewidth=1)+
    geom_point(aes(x = weekDate, y =prop, group=variant)) # the actual weekly proportions
  suppressWarnings(print(p)) # when forecasting a week without data, there's a warning on NA for no actual data
}

####################################
# Function to generate the week info plot with predicted proportions (predicted) vs selection coefficients for the variants (week.coef)
# input: dat: Table with the information for a specific week
#         byday: If TRUE, divide the selection coefficients by 7 to get the daily coefficients
####################################
weekInfoPlot <- function(dat, byday = T){
  if(byday){
    lab4y = "Selection coefficient (Daily)"
    dat[,c("week.coef", "se", "CI.week.2.5..", "CI.week.97.5..")]=dat[,c("week.coef", "se", "CI.week.2.5..", "CI.week.97.5..")]/7
  }else{
    lab4y = "Selection coefficient (Weekly)"
  }
  
  ## converting to predicted %
  dat[,c("predicted", "std.error","conf.low","conf.high")]= dat[,c("predicted", "std.error","conf.low","conf.high")]*100
  
  # generate plot
  weekp <- ggplot(dat, aes(x=predicted, y=week.coef))+
    geom_point()+
    geom_text_repel(data=dat, aes(x=predicted, y=week.coef,
                                  label=variant, size=5, col=variant))+ #, max.overlaps = Inf
    geom_errorbarh(aes(xmin=conf.low, xmax=conf.high, colour=variant))+
    geom_errorbar(aes(ymin=CI.week.2.5.., ymax=CI.week.97.5.., colour=variant))+
    xlab(paste("Predicted proportion (%) -", dat$weekDate[1]))+
    ylab(lab4y) +
    scale_x_log10()+
    theme(legend.position="none")
  print(weekp)
}



####################################
## Run Nowcast modeling
## Returns a list
####################################
runNowcast <- function(timeDat.s, #data table with 
                       week0day1, 
                       time_end, 
                       toestim="standard50", ## defaulted to choosing anything above 50 sequences in the model and projecting weeks; "loi" for own choice of lineages to model; 
                       loi=NULL, ## specify the loi here if setting toestim to be "loi"
                       supergroups=NULL, ## for 2tier setup - specify the supergroups here if setting toestim to be "webteamX"
                       ref2use="Other", ## "Other":using Other as reference; "top": using top lineage; Or specify the lineage directly
                       model_weeks=8, ## modelling 8 weeks until end the actual model end date
                       nweekPred=2, ## projecting 2 weeks onward
                       output_dir="output", ## directory to save the figures and tables
                       linhierfile, ## the lineage file
                       modelw="Freq", ## using "Freq" for model weight or "prop"
                       weight.adj=TRUE, ## if extreme weights are to be adjusted
                       inclusive=FALSE ## if including the last nweekPred weeks (incomplete) in the modelling
){
  
  ## create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE)
  
  ## save params for return
  paramsUsed <- list(
    week0day1 = week0day1,
    time_end = time_end,
    toestim = toestim, 
    loi = loi,
    supergroups = supergroups,
    ref2use = ref2use, 
    model_weeks = model_weeks, 
    nweekPred = nweekPred, 
    output_dir = output_dir, 
    modelw = modelw, 
    weight.adj = weight.adj, 
    inclusive =  inclusive
  )
  
  
  ##########################################################
  # Model dates & setting training+validation weeks
  ##########################################################
  # Maximum week included in the model -> standardized to a Sunday based upon static start week0day1 (start of epiweek)
  model_week_max <- as.numeric(as.Date(time_end) - week0day1) %/% 7
  #Convert maximum modelled week to a date (+6 to specify end of that epiweek--a Saturday)
  actual_time_end <- model_week_max*7 + week0day1 + 6
  print(paste("The actual model end date is",actual_time_end))
  
  # Setting the training weeks and validation weeks
  
  #Training weeks
  tweeks <- (model_week_max-model_weeks+1):model_week_max
  
  if(inclusive){ ## if modelling including the last two incomplete weeks
    #Validation weeks: 
    vweeks <- (model_week_max-nweekPred+1):(model_week_max)
  }else{ ## if not modelling the last two incomplete weeks
    #Validation weeks
    vweeks <- (model_week_max+1):(model_week_max + nweekPred)
  }
  
  #####################################################
  ###### sorting out dates/weeks
  #####################################################
  # 1 week before first week included in the model
  model_week_min <- model_week_max - model_weeks
  # a midpoint week that will be used to center week values
  model_week_mid <- round(model_weeks/2) # center it around 0 instead of using 1:model_weeks
  # add the new model_week info to the dataframe
  timeDat.s$model_week <- timeDat.s$week - model_week_min - model_week_mid
  
  
  #####################################################
  ###### creating the variant column, making groups
  #####################################################
  ## if the choice is to estimate for specific lineages specify here
  ## can make toestim to be a scalar or vector later on
  if(toestim=="loi"){
    if(is.null(loi)){ # if loi isn't specified
      stop("The loi parameter is needed to specify the lineages to model in loi mode.")
    }
    linFocus=loi; names(linFocus)=loi
  }
  
  ## if using top5 lineages of the last two training weeks
  if(length(grep("^top",toestim))==1){
    ## get top n specified
    topFilter=as.numeric(gsub("top","",toestim))
    
    ## getting the top 5 lineages and the counts within the last 2 weeks of training period, otherwise it's too far back to use all weeks of training period
    topnWeek=2
    toestim=paste0(toestim,"last",topnWeek)
    topn <- head(sort(table(timeDat.s[timeDat.s$week>=(model_week_max-(topnWeek-1)) & timeDat.s$week <= model_week_max,"pango_lineage"]),decreasing=T), topFilter)
    print(topn)
    linFocus=topn
  }
  
  ## if the choice is to estimate lineages with at least x% within modelling perio, specify percx
  if(length(grep("^perc",toestim))==1){
    ## get percentage specified
    percFilter=as.numeric(gsub("perc","",toestim))
    
    ## aggregate weekly data into long format
    WeekCntsLin <- data.frame(table(timeDat.s$pango_lineage,timeDat.s$week))
    colnames(WeekCntsLin)=c("variant","week", "Freq")
    WeekCntsLin$week=as.numeric(as.character(WeekCntsLin$week)) # make sure the week column is numeric
    
    propModelWeeks <- do.call(rbind,lapply(tweeks,function(y){
      tmp <- WeekCntsLin[WeekCntsLin$week==y,]
      data.frame(tmp,prop=tmp$Freq/sum(tmp$Freq))
    }))
    percLin <- as.character(unique(propModelWeeks[propModelWeeks$prop >= percFilter/100,"variant"]))
    print(percLin)
    linFocus=percLin; names(linFocus)=percLin
  }
  
  ## From Nelson: if the choice is to estimate lineages with at least top counts within modelling period, specify countsx
  if(length(grep("^counts",toestim))==1){
    ## get counts specified
    countsFilter <- as.numeric(gsub("counts","",toestim))
    WeekCntsLin <- timeDat.s %>% group_by(week) %>% filter(week %in% (model_week_max-model_weeks+1):model_week_max) 
    WeekCntsLin <- as.data.frame(table(c(WeekCntsLin$pango_lineage)))
    WeekCntsLin <- WeekCntsLin %>% filter(Freq >= countsFilter) %>% dplyr::rename("variant" = "Var1")
    linFocus <- as.vector(WeekCntsLin$variant); names(linFocus) <- linFocus
  }
  
  ## If the choice is to estimate lineages with at least x counts per week within modelling period + 2 prediction weeks, specify standardx
  if(length(grep("^standard",toestim))==1){
    ## get counts specified
    countsFilter <- as.numeric(gsub("standard","",toestim))
    WeekCntsLin <- timeDat.s%>% filter(week %in% c(tweeks,vweeks)) 
    WeekCntsLin <- as.data.frame(table(WeekCntsLin$pango_lineage, WeekCntsLin$week))
    WeekCntsLin <- WeekCntsLin %>% filter(Freq >= countsFilter) %>% dplyr::rename("variant" = "Var1")
    linFocus <- as.vector(unique(WeekCntsLin$variant)); names(linFocus) <- linFocus
  }
  
  ## If the choice is webteamX (2filter) which uses standardX for modelling but classifies these organisms under an LOI supergroup
  if(length(grep("^webteam",toestim))==1){
    if(grep(toestim, pattern="^webteam") == 1 & is.null(supergroups)){ #1st - set supergroups if left at NULL (note - this may fail over time due to the natural dissapearance of these supergroups)
      # Parental supergroups
      # If no default parents/supergroups are set when using the "webteam" filter
      stop("The supergroups parameter is needed to specify the lineages to model in webteamX mode.")
    }
    #Filter1 - standard50: this checks for sublineages with >= 50 counts within the specified x modelling (model_weeks) + y forecasted (nweekPred) weeks
    ## get counts specified
    countsFilter <- as.numeric(gsub("webteam","",toestim)) #set filter criteria to counts over 50 due to prior assessment
    WeekCntsLin <- timeDat.s%>% filter(week %in% c(tweeks,vweeks)) #filter for all lineages within x modelling + y forecasted weeks
    WeekCntsLin <- as.data.frame(table(WeekCntsLin$pango_lineage, WeekCntsLin$week)) #generate a table of counts for these lineages
    countsFilterPassed <- WeekCntsLin %>% filter(Freq >= countsFilter) %>% pull(Var1) #filter for lineages with counts >= 50 over this x+y period
    WeekCntsLin <- WeekCntsLin %>% filter(Var1 %in% countsFilterPassed) %>% dplyr::rename("variant" = "Var1") #filter for these passed lineages from the tweek/vweek filtered table to get a complete counts table used to determine top5 per supergroup in filter#2    
    #Filter2 - Group the passed lineages back into executively decided supergroups
    ## Get parental supergroups & re-group linFocus into these supergroups
    # get parental groupings
    supergroups <- supergroups
    names(supergroups) <- supergroups
    linhier <- read.csv(linhierfile) ##load lineage hierarchy file
    linDescWebteam<-lapply(names(supergroups), getSublineages, linhier) #create a comprehensive list of all sublineages falling under each of these specified supergroups
    names(linDescWebteam)=names(supergroups)
    linDescWebteam<-lineageExclusive(linDescWebteam) # making sure the descendants weren't redundant for the targeted lineages
    # Re-group lineages passing filter#1 back into these supergroups
    #Assign supergroups to variants in a new column in propDF2.1a
    for(i in 1:nrow(WeekCntsLin)){
      #Note - some lineages (Like XBM) does not follow into the strict parental groups which messes up the search and placement code below. Therefore, set them into the "Others" group
      #1st check if variant can be found in parental groupings
      ifelse(any(!is.na(lapply(linDescWebteam, function(x) which(x %in% WeekCntsLin$variant[i])) >= 1)) !=F,
             #If yes, then assign them
             yes=WeekCntsLin$parent[i] <- names(linDescWebteam)[!is.na(lapply(linDescWebteam, function(x) which(x %in% WeekCntsLin$variant[i])) >= 1)],
             #If not, then put them into their own supergroup -> put them into the "Other" supergroup instead of forming their own
             no=WeekCntsLin$parent[i] <- as.character("Other")) #This assigns unspecified lineages into the "Other" supergroup if they don't fall within specified supergroups
    }
    # Check if any supergroup contains >5 variants - for accessibility purposes, we re-group sublineages with their closest ancestor if >5 sublineages fall under any given supergroup
    #Summarize the groupings 
    summarizedgroups0 <- WeekCntsLin %>% group_by(parent) %>% dplyr::summarise(n_distinct(variant), .groups="keep") #note - during automation, using summarise(unique(variant)) will shoot an error, solved this using n_distinct 
    #get counts of #variants mapped to a parent
    summarizedgroups <- as.vector(summarizedgroups0$`n_distinct(variant)`); names(summarizedgroups) <- summarizedgroups0$parent
    # For supergroups with >5 variants, find which have the lowest counts over the x+y model weeks
    #Subset the weeklyextra50counts (filter#1) table for supergroups w/ >5 variants & summarize variant counts for all weeks (each supergroup)
    groupsgreaterthan5 <- WeekCntsLin %>% filter(parent %in% names(which(summarizedgroups>5))) %>% group_by(variant, parent) %>% dplyr::summarise(total=sum(Freq)) %>% arrange(total)
    #Filter for the lowest count variants to ensure only 5 variants exist in each parent
    lowestcountvariants <- groupsgreaterthan5 %>% filter(!variant %in% 
                                                           (groupsgreaterthan5 %>% group_by(parent) %>% dplyr::top_n(total, n=5) %>% pull(variant)))
    # Remove these low count variants from weeklyextra50counts (filter#1) table 
    WeekCntsLin <- WeekCntsLin %>% filter(!variant %in% lowestcountvariants$variant)  
    # Set these remaining sublineages for downstream modelling
    linFocus <- as.vector(unique(WeekCntsLin$variant)); names(linFocus) <- linFocus
    #Check if the specified supergroups are in linFocus, if not then add them back in
    linFocus <- c(linFocus, setdiff(supergroups, linFocus)[!grepl(setdiff(supergroups, linFocus), pattern="\\.\\*")]) ; names(linFocus) <- linFocus #Exclude any supergroups with ".*")
  }
  
  ###############
  ### get the list of lineages to focus on through user criteria and add the variant column
  ###############
  # leveraging lineage hierarchy file 
  linhier <- read.csv(linhierfile)
  
  # Returning the list of lineages with their sublineages
  linDesc <- lapply(names(linFocus), getSublineages, linhier=linhier)
  names(linDesc)=names(linFocus)
  
  # making sure the descendants weren't redundant for the targeted lineages
  linDesc <- lineageExclusive(linDesc)
  
  #230112 Return final lineage groupings to be reported
  lineageGroupings <- linDesc
  
  # add a column that notes lineage of interest and its descendants to model, leave the rest as Other
  timeDat.s$variant=timeDat.s$pango_lineage
  
  # label all corresponding descendants to the target parental lineage
  for(lin in names(linFocus)){
    timeDat.s$variant[timeDat.s$variant%in%linDesc[[lin]]]=lin
  }
  
  # make the rest that havent been renamed as Other
  timeDat.s$variant[!timeDat.s$variant%in% names(linFocus)]="Other"
  
  variantlevels <- c("Other",names(linFocus)) # setting Other as reference
  timeDat.s$variant <- factor(timeDat.s$variant,levels=variantlevels)
  
  ## relevel the reference to the variant with highest prevalence in the modelling period - 221108 Nelson added by request
  if(ref2use=="top"){ 
    WeekCntsLin2 <- timeDat.s %>% filter(week %in% c(tweeks))   
    WeekCntsLin2$variant <- as.character(WeekCntsLin2$variant) #remove factor levels
    WeekCntsLin2 <- as.data.frame(table(c(WeekCntsLin2$variant))) #Get groups to select ref (b/c in loi the top count lineage is likely grouped within parent so relevel doesn't work)
    ref2use = as.character(WeekCntsLin2[which.max(WeekCntsLin2$Freq),"Var1"])
    
  }else{ ## pre-specified to one lineage
    #make sure specified lineage is within the ones to model
    if(!ref2use %in% levels(timeDat.s$variant)){
      warning("The reference lineage supplied is not listed to be modelled!")
    }
  }
  
  ## setting reference accordingly
  print(paste("Lineage reference used",ref2use))
  timeDat.s$variant <- relevel(timeDat.s$variant, ref= ref2use)
  
  ## aggregate daily data into long format
  dateCntsLong <- data.frame(table(timeDat.s$variant,timeDat.s$date))
  colnames(dateCntsLong)=c("variant","date","Freq")
  
  ## getting the total number of sequencing data by day
  dateCntsTotal <- aggregate(data.frame(total=dateCntsLong$Freq), list(dateCntsLong$date),FUN = sum)
  rownames(dateCntsTotal)=dateCntsTotal[,1]
  dateCntsLong <- data.frame(dateCntsLong,dateCntsTotal[dateCntsLong$date,] )
  dateCntsLong$prop <- dateCntsLong$Freq/dateCntsLong$total
  dateCntsLong$date=as.Date(as.character(dateCntsLong$date))
  dateCntsLong$weeknum=(dateCntsLong$date-week0day1)/7
  
  #########################################################
  ### Model nnet: aggregated input with weight
  #########################################################
  ## aggregate weekly data into long format
  weekCntsLong <- data.frame(table(timeDat.s$variant,timeDat.s$week))
  colnames(weekCntsLong)=c("variant","week", "Freq")
  weekCntsLong$week=as.numeric(as.character(weekCntsLong$week)) # make sure the week column is numeric
  
  ## getting the total number of sequencing data by week
  totalseqByweek=aggregate(data.frame(weekCntsLong$Freq), list(weekCntsLong$week),FUN = sum)
  rownames(totalseqByweek)=totalseqByweek[,1]
  
  ## adding in the corresponding total of the week as a column
  weekCntsLong <- data.frame(weekCntsLong, total=totalseqByweek[as.character(weekCntsLong[,"week"]),"weekCntsLong.Freq"])
  
  weekCntsLong <- data.frame( weekCntsLong, prop=weekCntsLong$Freq/weekCntsLong$total, weekDate=weekCntsLong$week*7+week0day1)
  
  ## setting the training weeks and validation weeks
  ## Subsetting data for model training and validating/predicting
  tandvWeeks <- unique(c(tweeks, vweeks))
  dataTV <- weekCntsLong[weekCntsLong$week%in%tandvWeeks, ]
  
  ## creating weights according to the weekly counts or proportion
  ## modelw="Freq": Modelling using counts
  weight=dataTV[,modelw]

  ## Weight adjustment for extreme weights, currently adjusting only the minimum weight to be nonzero
  if(weight.adj){
    ## Taking extremes only within the modelled period
    ## maximum weight; trim weights > 99th percentile
    #max_weight <- quantile(x = weight[dataTV$week%in%tweeks], probs = .99, na.rm=T) 
    ## get the minimum weight to replace weights of 0
    min_weight <- min(weight[dataTV$week%in%tweeks][weight[dataTV$week%in%tweeks] > 0], na.rm=T) 
    
    weight[weight==0]=min_weight
    #weight[weight>max_weight]=max_weight
  }
  
  dataTV <- data.frame(dataTV, weight=weight)
  
  ## centering the weeks for modelling
  weekStd = floor(mean(dataTV$week))
  dataTV$week = dataTV$week - weekStd
  
  ## fitting a model with the weights
  mnModel3 <- nnet::multinom(formula = formula("variant  ~ week"),
                             data      = dataTV[dataTV$week%in%(tweeks- weekStd), ],
                             weights   = dataTV[dataTV$week%in%(tweeks- weekStd),"weight"],
                             maxit     = 1000,
                             Hess      = TRUE,
                             trace     = FALSE)
  
  ## get model coefficients, se and the confidence intervals
  mnMOdel3.coef <- coef(mnModel3)
  mnModel3.coefDF <- data.frame(week.coef=coef(mnModel3)[,2], 
                                se=summary(mnModel3)$standard.errors[,2], 
                                CI=data.frame(apply(confint(mnModel3, level=0.95),c(1,2),function(y)y))[,c(2,4)])
  
  ## prediction dates
  prediction_weeks=data.frame(week=seq(min(dataTV$week), max(vweeks)-weekStd,1/7))
  
  
  ########################################
  ## getting PI through Hessian matrix
  ########################################
  ## for getting the predicted proportions and se, could rewrite to go aroundthis
  predemm <- ggemmeans(mnModel3, terms = paste0("week [",paste(c(min(dataTV$week), max(vweeks)-weekStd),collapse=":"),", by=0.25]"), ci.lvl = 0.95)
  predemmdf <- as.data.frame(predemm)
  colnames(predemmdf)[1] = "week"
  colnames(predemmdf)[colnames(predemmdf)=="response.level"] = "variant"
  
  ## using Hessian matrix for sampling coefficients to get PIs
  ## using Hessian matrix
  sampTimes <- 1000
  paramFormat <- reshape2::melt(coef(mnModel3))
  paramFormat2 <- paramFormat[,3]
  names(paramFormat2) = paste(paramFormat[,1],paramFormat[,2],sep=":")
  
  ## Sampling with Hessian
  dfHessian <- RandomFromHessianOrMCMC(Hessian=mnModel3$Hessian, fitted.parameters=paramFormat2,
                                       method="Hessian", replicates=sampTimes, silent=T)$random  
  
  ## get coefficient CI from RandomFromHessian, similar to the confint from model
  coefCI <- t(apply(dfHessian,2,quantile,probs=c(.025,.975)))
  coefCI <- coefCI[seq(2,nrow(coefCI),by=2),]; rownames(coefCI)=gsub(":week","",rownames(coefCI ))
  
  ## function to get 95 prediction intervals from sampled coefficients
  ## usage: hessianPI <- getPI_HelpersMG(dfHessian, weeki=4)
  getPI_HelpersMG <- function(sampledCoef, weeki){
    epredictMat <- do.call(rbind,lapply(1:nrow(sampledCoef), function(i){
      #i=1 # the iteration
      interCoef <- t(matrix(as.numeric(sampledCoef[i,]),nrow=2))
      tmpc <- exp(c(0,weeki*interCoef[,2]+interCoef[,1]))
      pprop <- tmpc/sum(tmpc) # the proportion
      pprop
    }))
    colnames(epredictMat) = c("ref",gsub(":week","",names(sampledCoef)[seq(2,ncol(sampledCoef), by=2)]))
    
    t(apply(epredictMat, 2, function(y){c(mean=mean(y), quantile(y, probs=c(.025, .975)))}))
  }
  
  ## getting PI for the period
  hPIs <- do.call(rbind,lapply(sort(unique(predemmdf$week)),function(y){
    getPI_HelpersMG(dfHessian, weeki=y)
  }))
  
  predemmdf$conf.low=hPIs[,2]
  predemmdf$conf.high=hPIs[,3]
  
  ## to transform the week back from standardization
  predemmdf$week = predemmdf$week + weekStd
  dataTV$week = dataTV$week + weekStd
  
  predemmdf <- data.frame(predemmdf, weekDate=predemmdf$week*7+week0day1)
  propDF <- merge(dataTV, predemmdf, by=c("weekDate","variant"),all.y=T)
  
  ## In case the prediction weeks (vweeks) extend beyond data availability
  propDF[propDF[,"week.y"]%in%vweeks,"week.x"] = propDF[propDF[,"week.y"]%in%vweeks,"week.y"]
  
  dateCntsLong.s <- dateCntsLong[dateCntsLong$weeknum>=min(tandvWeeks)&dateCntsLong$weeknum<=max(tandvWeeks),]
  
  ########################################
  ## N: Due to popular demand, append asterisks to lineages that are grouped during modeling & add this column to propDF
  ## This is needed for the webteam version
  ########################################                                                  
  ## Find which focused lineage groups contain sublineages
  lineageGroupingsWsublineages <- lapply(linDesc, function(x){length(x) > 1}) %>% unlist()      
  ## Append asterisk to new column in propDF
  ## Initialize empty column
  propDF$sublineage_grouped <- NA
  ## Append "*" if grouped by sublineages
  for(i in 1:nrow(propDF)){
    #For each variant row (propDF$variant[i]), isolate said variant from lineageGroupingsWsublineages & check if contains sublineages
    if(isTRUE(lineageGroupingsWsublineages[names(lineageGroupingsWsublineages) %in% propDF$variant[i]])){
      propDF$sublineage_grouped[i] <- paste0(propDF$variant[i], "*")
    }else{
      #If no sublineage grouped, then don't append asterisk
      propDF$sublineage_grouped[i] <- paste0(propDF$variant[i]) #paste0() b/c in factor format, otherwise it pulls the index#
    }
  }
  ## Set as factor
  propDF$sublineage_grouped <- factor(propDF$sublineage_grouped, levels=unique(propDF$sublineage_grouped))
  
  ## Also append "*" to dateCntsLong.s
  for(i in 1:nrow(dateCntsLong.s)){
    #For each variant row (propDFnoDaily$variant[i]), isolate said variant from lineageGroupingsWsublineages & check if contains sublineages
    if(isTRUE(lineageGroupingsWsublineages[names(lineageGroupingsWsublineages) %in% dateCntsLong.s$variant[i]])){
      dateCntsLong.s$sublineage_grouped[i] <- paste0(dateCntsLong.s$variant[i], "*")
    }else{
      #If no sublineage grouped, then don't append asterisk
      dateCntsLong.s$sublineage_grouped[i] <- paste0(dateCntsLong.s$variant[i]) #paste0() b/c in factor format, otherwise it pulls the index#
    }
  }
  ## Set as factor
  dateCntsLong.s$sublineage_grouped <- factor(dateCntsLong.s$sublineage_grouped, levels=unique(dateCntsLong.s$sublineage_grouped))     

  ########################################   
  ## Generating plots and results                                                 
  ########################################   
  ## Save trend plot
  pdf(file.path(output_dir,paste0("prop_trend_",actual_time_end,"_",toestim,"_",model_weeks,"_ref",ref2use,"_",modelw,".pdf")),height=6, width=8.5)
  trendPlot(propDF, actual_time_end)
  dev.off()
  ## PropDF with coefficients for table generation
    #Format coefficient table for join
    cjt <- rbind(rep(0,ncol(mnModel3.coefDF)), mnModel3.coefDF) %>% #Add 0 values for the reference organism into the coefficient table                
      rownames_to_column("variant")
    cjt[1,1] <- ref2use #change the name of the variant for newly added row pertaining to the reference from 1 -> ref 
    #Add these coef. values into the proportions table
    propDF.r <- propDF %>% left_join(cjt, by="variant") 
    
  ## Week info plot
  propDF.s <- propDF[which(propDF$week.x==max(vweeks)),]
  propDF.s <- data.frame(propDF.s, rbind(rep(0,ncol(mnModel3.coefDF)), mnModel3.coefDF))
    
  pdf(file.path(output_dir,paste0("prop_coef_",actual_time_end,"_",toestim,"_",model_weeks,"_ref",ref2use,"_",modelw,".pdf")),height=6, width=6)
  ## generate plot
  weekInfoPlot(propDF.s)
  dev.off()
  
  ### using propDF
  propDF$prop.diff <- abs(propDF$prop-propDF$predicted)
  propDF$withinInd <- propDF$week.x %in% tweeks
  propDF$projInd <- propDF$week.x %in% vweeks
  
  ### save output
  write.table(propDF %>% mutate_if(is.numeric, round, digits = 6),
              sep="\t", file=file.path(output_dir,paste0("prop_",actual_time_end,"_",toestim,"_",model_weeks,"_ref",ref2use,"_",modelw,".txt")), quote=F, row.names=F, col.names=T)
  ### saving the very last week, and including week coefficients as well as it's CI from the model
  write.table(propDF.s %>% mutate_if(is.numeric, round, digits = 6),
              sep="\t", file=file.path(output_dir,paste0("prop_",actual_time_end,"_",toestim,"_",model_weeks,"_ref",ref2use,"_",modelw,"_lastweek.txt")), quote=F, row.names=F, col.names=T)
  
  
  ## mse (per week) for modelled weeks
  WithinMmsq <- sum(propDF[propDF$withinInd,"prop.diff"]^2)/length(unique(propDF[propDF$withinInd,"week.x"]))/(length(linFocus)+1)
  print(paste("Within model mse",round(WithinMmsq,6)))
  
  ## mse (per week) for projected weeks
  ProjMmsq <- sum(propDF[propDF$projInd,"prop.diff"]^2)/length(unique(propDF[propDF$projInd,"week.x"]))/(length(linFocus)+1)
  print(paste("Prediction mse",round(ProjMmsq,6)))

  ######################################## 
  ## "webteamX" model - additional processing of output file
  ## This is needed for the webteam version
  ######################################## 
  ## Appending the supergroup umbrella lineages to each variant in the modelled output(propDF)                                                   
    if(isTRUE(grep(toestim, pattern = "^webteam") == 1)){ # Check if "webteamX" filter is specified  
      propDF$supergroup <- NA #If yes, then append a "supergroup" column to the output table
      for(i in 1:nrow(propDF)){ #Map the specified supergroups to the filtered variant lineages
        ifelse(any(WeekCntsLin$variant %in% propDF$variant[i]), #check if variant[i] is mapped to supergroup (weekcntslin) as set by filter above
               yes=propDF$supergroup[i] <- unique(WeekCntsLin$parent[which(WeekCntsLin$variant %in% propDF$variant[i])]), #If running webteam (2filter) then append a supergroups column mapping each focused lineage
               no=propDF$supergroup[i] <- "Other" #If not in supergroups then placed into "Other" as supergroup
        )
      }
      #One additional check for supergroups lacking the matching sublineage  i.e) the specified supergroup lacks itself as a variant because it did not pass the webteamX filter (i.e. XBB.1.5 supergroup lacks XBB.1.5 as major lineage) - these need to be set within their own supergroup to extend beyond the top5 limit
      propDF$supergroup[which(propDF$variant %in% linFocus[which(!linFocus %in% unique(WeekCntsLin$variant))])] <- linFocus[which(!linFocus %in% unique(WeekCntsLin$variant))]
      
      #Set as factor
      propDF$supergroup <- factor(propDF$supergroup, levels=unique(propDF$supergroup))
      propDF <- propDF %>% relocate(supergroup, .after=sublineage_grouped) #move the supergroup column
      
      #Filter propDF output table down
      propDFweb <- propDF %>% filter(!is.na(week.x)) #remove daily data
      ifelse(length(propDFweb %>% filter(projInd == T) %>% distinct(weekDate) %>% pull()) >= 3, #Remove 3rd prediction week (only want to show 2 on the webside)
             yes=propDFweb <- propDFweb[!(propDFweb$week.x %in% (propDFweb %>% group_by(week.x) %>% tally() %>% pull(week.x) %>% tail(., n=1))),], #Filter away 3rd predicted week
             no= propDFweb <- propDFweb
      )
    
    ## Reformat output CSV tables to a different format - can be used for a more targeted way of reporting
      #CSV file#1 - actual proportions (8 model weeks) + predicted proportions (n predicted weeks) - IN PERCENTAGES (*100)
      propDF2.1a <- propDFweb %>% dplyr::select(c("weekDate", "variant", "predicted", "prop", "week.x", "sublineage_grouped", "supergroup", "conf.low", "conf.high", "projInd"))
      #Want first 8 model weeks as ACTUAL proportions and the last 2 weeks as PREDICTED proportions
      propDF2.1a$proportions <- c(propDF2.1a %>% filter(projInd == F) %>% pull(prop)*100, #Actual proportions for modelled 8 weeks
                                  propDF2.1a %>% filter(projInd == T) %>% pull(predicted)*100) #Predicted proportions for projected 2 weeks
      #Remove CI for the 8 modelled weeks (actual proportions) - only want CI for the forecasted weeks
      propDF2.1a$conf.h <- c(rep(NA, times=nrow(propDF2.1a %>% filter(projInd == F))),
                             propDF2.1a %>% filter(projInd == T) %>% pull(conf.high)*100)
      propDF2.1a$conf.l <- c(rep(NA, times=nrow(propDF2.1a %>% filter(projInd == F))),
                             propDF2.1a %>% filter(projInd == T) %>% pull(conf.low)*100)
      #Remove unnecessary columns
      propDF2.1a <- propDF2.1a %>% dplyr::select(-c("predicted", "prop", "week.x", "projInd", "conf.high", "conf.low"))
      #Set weekDate as date format
      propDF2.1a$weekDate <- as.Date(propDF2.1a$weekDate)
      #Reorder and re-name
      propDF2.1a <- propDF2.1a %>%
        dplyr::rename("week of collection" = "weekDate") %>% 
        dplyr::rename("variant_aggregated" = "sublineage_grouped") %>% 
        group_by(supergroup) %>%
        arrange(`week of collection`)
      
      #CSV file#2 - table of total "actual" counts per epiweek (excluding forecast weeks)
      propDF2.1b <- propDFweb %>% dplyr::select(c("weekDate", "total")) %>% distinct()
      propDF2.1b.init <- propDF2.1b #save df to return before replacing last 2 weeks
      #Remove counts for the 2 forecasted weeks
      propDF2.1b$total <- replace(propDF2.1b$total, propDF2.1b$total==tail(propDF2.1b$total, n=2), "-")
      #Set weekDate as date format
      propDF2.1b$weekDate <- as.Date(propDF2.1b$weekDate)
      propDF2.1b <- propDF2.1b %>%
        dplyr::rename("week of collection" = "weekDate") %>% 
        arrange(`week of collection`)

    ## Write CSV tables
    # write.csv(propDF2.1a, file.path(output_dir,paste0(Sys.Date(), "_NMLNowcast_WebteamVerC_proportions.csv")), row.names = F)
    # write.csv(propDF2.1b, file.path(output_dir,paste0(Sys.Date(), "_NMLNowcast_WebteamVerC_actualcounts.csv")), row.names = F)
    write.table(propDFweb %>% mutate_if(is.numeric, round, digits = 6),
                sep="\t", file=file.path(output_dir,paste0("propWeb_",actual_time_end,"_",toestim,"_",model_weeks,"_ref",ref2use,"_",modelw,".txt")), quote=F, row.names=F, col.names=T)
    
  }else{
    propDF <- propDF #No supergroup column if running any filter other than "webteamX"
    propDF2.1a <- NULL #Blank return element - only evaluated for webteam filter but not other filters
    propDF2.1b.init <- NULL #Blank return element - only evaluated for webteam filter but not other filters
  }
  
  ## return the list
  return(
    list(propDF = propDF, 
         propDF.s = propDF.s, 
         propDF.r = propDF.r,
         WithinMmsq = WithinMmsq,
         ProjMmsq = ProjMmsq,
         actual_time_end = actual_time_end,
         theModel = mnModel3,
         mnMOdel3.coef = mnMOdel3.coef,
         tweeks = tweeks,
         vweeks = vweeks,
         mnModel3.coefDF = mnModel3.coefDF,
         dateCntsLong.s = dateCntsLong.s,
         ref2use = ref2use,
         propDF2.1a = propDF2.1a,
         propDF2.1b.init = propDF2.1b.init,       
         params = paramsUsed
    )
  )
  
}

####################################
### generate logit plots from the Nowcast model
### automated logit plot for all variants against the reference
### The input is the runNowcast output (list)
### The output is a list of logit plots
### tosave: whether to print & save the figures into a pdf file within specified output folder
####################################
logitplot <- function(ncOutput, tosave=FALSE){
  propDF0 <- ncOutput$propDF
  preds <- dcast(propDF0[,c("weekDate","variant","predicted")], weekDate~variant, value.var="predicted") 
  
  propDF0 <- na.omit(propDF0) # to get weekly counts only
  freqs <- dcast(propDF0[,c("weekDate","variant","Freq")], weekDate~variant, value.var="Freq") 
  
  # get names of all variants modeled
  varis <- levels(ncOutput$propDF$variant)
  
  loplot <- list()
  
  #loop through all variants & compare to ref, the first level
  yrange=c(-3,3)
  for(i in 2:length(varis)){
    ## log10 can be -Inf for newly emerging lineages that had 0 counts in earlier weeks, the points are plotted at the bottom of the figure
    freqs$logit = log10(freqs[,varis[i]]/freqs[,varis[1]])
    preds$logit = log10(preds[,varis[i]]/preds[,varis[1]])
    
    ## saving the ranges for all variants to specify ylim later on
    allys <- c(yrange,preds$logit, freqs$logit)
    yrange = range(allys[allys != -Inf])
    
    loplot[[i-1]] <- ggplot(freqs, aes(x=weekDate, y=logit))+
      geom_point()+
      geom_line(data=preds, aes(x=weekDate, y=logit, group=1))+
      ylab(paste("log10(",varis[i],"/",varis[1],")"))+
      theme(legend.position="none")
  }
  ## add the range as ylim
  loplot <- lapply(loplot, function(y) y+ylim(yrange[1],yrange[2]))
  
  ## for saving the plots into a pdf file
  if(tosave){
    pdf(file.path(ncOutput$params$output_dir,paste0("logitplot_",ncOutput$actual_time_end,"_",ncOutput$params$toestim,"_",ncOutput$params$model_weeks,"_ref",ncOutput$params$ref2use,"_",ncOutput$params$modelw,".pdf")), height=min(8,2*ceiling(length(loplot)/3)), width=8)
    print(ggarrange(plotlist=loplot, ncol=3, nrow=min(4,ceiling(length(loplot)/3))))
    dev.off()
  }
  
  names(loplot) <- varis[2:length(varis)]
  loplot 
}
