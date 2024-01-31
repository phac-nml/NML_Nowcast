########################################################
## Plotting functions to be called 
## Written by Nelson Mok
## Modified and cleaned up on 2022-04-03
########################################################

#Trend plot
#Where
#dataTP = propDF (preferably propDF.original)
#npredweeks = number of predicted weeks to be visualized (params$nweekPred)
#nweeksmodel = number of weeks used for modeling (ie) trainng weeks - length(tweeks))
#dateformat = x label date format (as a character - i.e) "%b-%d")
#cutoffline = vline between model end & forecast dates
trendplot <- function(dataTP, npredweeks, nweeksmodel, dateformat, cutoffline){
  dat2plot = dataTP[!is.na(dataTP$week.x),]
  dat2plot[,c("prop","predicted", "std.error","conf.low","conf.high")] = dat2plot[,c("prop","predicted", "std.error","conf.low","conf.high")]*100
  TP <- ggplot(dat2plot %>% filter(weekDate <= actual_time_end+(npredweeks)*7), 
               aes(x = weekDate, y = round(predicted,2), 
                   group = sublineage_grouped, 
                   col= sublineage_grouped, 
                   fill=sublineage_grouped)) +
    geom_line() + #fitted line = predicted data
    ylab("Proportion (%)")+
    xlab(paste0(sub(tail(dataTP$weekDate, n=1), pattern="-.*", replacement=""),"\n Week of sample collection")) +
    theme(axis.title = element_text(size=15, face="bold"),
          axis.text.x=element_text(angle=45, vjust=1, hjust=1),
          plot.title=element_text(hjust=0.5),
          legend.title=element_blank()) +
    scale_x_date(breaks=dataTP %>% filter(!is.na(week.x)) %>% pull(weekDate), #This sets the weeks visualized on the plot
                 limits=c(actual_time_end-(nweeksmodel+1)*7,actual_time_end+(npredweeks)*7), #min must include 1 week before the number of model weeks (params$model_weeks) | max is equal to the number of forecasted weeks (params$nweekPred)
                 date_labels=dateformat, expand=c(0,-2.6)) +
    geom_ribbon(aes(ymin = round(conf.low,2), ymax = round(conf.high,2)), alpha = 0.1, linetype=0) + #CI
    geom_vline(xintercept=as.numeric(cutoffline), linetype="dashed", color="darkgrey", size=1)+
    geom_point(aes(x = weekDate, y = round(prop,2), group=sublineage_grouped))+ # the actual weekly proportions as points
    ggtitle(paste0("Nowcast - ", toestim, " Model"))
  return(TP)
}

#Stacked barplot
#Where 
#dataSBP = propDF (preferably propDF.originalnoDaily)
#forecastbarend = #weeks predicted in model (ie) params$nweekPred)
#nweekmodeled = number of weeks modeled (used for training)
#acutalorpred = show actual proportions or predicted proportions (prop = actual, predicted = forecasted) - set in character format
#datefmt = x label date format (as a character - i.e) "%b-%d")
#yinc = increments on the y-axis
stackedbarplot <- function(dataSBP, forecastbarend, nweekmodeled, actualorpred, datefmt, yinc){
  dataSBP1 <- dataSBP %>% filter(!is.na(week.x))
  if(actualorpred=="prop"){
    SBP <- ggplot(dataSBP1, aes(x=weekDate, y=round(prop*100, 2), 
                               fill=sublineage_grouped,
                               text=paste0("</br> Variant: ", sublineage_grouped,
                                           "</br> Week: ", weekDate,
                                           "</br> Counts Actual: ", Freq,
                                           "</br> Proportion Actual: ", round(prop*100, 2), "%"))) +
      geom_rect(aes(xmin=actual_time_end-2.5, xmax=actual_time_end+(forecastbarend)*7, #geom_rect = set where the nowcast forecast bar hovers above & how many weeks (black hover bar)
                    ymin=105, ymax=110), fill="black", alpha=0.3, inherit.aes = F) +
      geom_col(colour="black") + #geom_col = actual stacked bars per variant
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
            panel.grid = element_blank(),
            panel.background = element_rect(fill="white"),
            axis.title = element_text(size=15, face="bold"),
            legend.title = element_blank()) +
      scale_x_date(breaks=dataSBP1$weekDate, #This sets the number of weeks visualized on the plot
                   limits=c(actual_time_end-(nweekmodeled+1)*7,actual_time_end+(forecastbarend)*7), #min must include 1 week before the number of model weeks (params$model_weeks) | max is equal to the number of forecasted weeks (params$nweekPred)
                   date_labels=datefmt, expand=c(0,-2.6)) + #Set "expand" to limit xlimits. To view n weeks, start from model end date and subtract (n*7) weeks
      scale_y_continuous(breaks = seq(0,100,by=yinc)) +
      xlab(paste0(sub(tail(dataSBP1$weekDate, n=1), pattern="-.*", replacement=""),"\n Week of sample collection")) + ylab(element_blank())
  }else{
    SBP <- ggplot(dataSBP1, aes(x=weekDate, y=round(predicted*100, 2), 
                               fill=sublineage_grouped,
                               text=paste0("</br> Variant: ", sublineage_grouped,
                                           "</br> Week: ", weekDate,
                                           "</br> Proportion Predicted: ", round(predicted*100, 2), "%"))) +
      geom_rect(aes(xmin=actual_time_end-2.5, xmax=actual_time_end+(forecastbarend)*7, #geom_rect = set where the nowcast forecast bar hovers above & how many weeks
                    ymin=105, ymax=110), fill="black", alpha=0.3, inherit.aes = F) +
      geom_col(colour="black") + #geom_col = actual stacked bars per variant
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
            panel.grid = element_blank(),
            panel.background = element_rect(fill="white"),
            axis.title = element_text(size=15, face="bold"),
            legend.title = element_blank()) +
      scale_x_date(breaks=dataSBP1$weekDate, #This sets the number of weeks visualized on the plot
                   limits=c(actual_time_end-(nweekmodeled+1)*7,actual_time_end+(forecastbarend)*7), #min must include 1 week before the number of model weeks (params$model_weeks) | max is equal to the number of forecasted weeks (params$nweekPred)
                   date_labels=datefmt, expand=c(0,-2.6)) + #Set "expand" to limit xlimits. To view n weeks, start from model end date and subtract (n*7) weeks
      scale_y_continuous(breaks = seq(0,100,by=yinc)) +
      xlab(paste0(sub(tail(dataSBP1$weekDate, n=1), pattern="-.*", replacement=""),"\n Week of sample collection")) + ylab(element_blank())
  }
  return(SBP)
}


#Pivot table
#Where datapvt = propDF (recommended propDF.original which has the appended asterisk to grouped lineages)
plotpivottable <- function(datapvt){
  propDFpivot <- rpivotTable(datapvt %>% 
                               filter(!is.na(week.x)) %>% 
                               dplyr::select(weekDate, sublineage_grouped, prop, predicted, Freq) %>% #sublineage_grouped only found in propDF.original (so won't work for just propDF)
                               mutate_at("prop", ~prop*100) %>% 
                               mutate_at("predicted", ~predicted*100) %>% 
                               dplyr::rename("Actual Proportion (%)" = "prop") %>% 
                               dplyr::rename("Predicted Proportion (%)" = "predicted") %>% 
                               dplyr::rename("Week" = "weekDate") %>% 
                               dplyr::rename("Variant" = "sublineage_grouped") %>% 
                               dplyr::rename("Counts" = "Freq"),
                             rows="Variant", cols="Week", aggregatorName="Sum", vals="Predicted Proportion (%)")
  return(propDFpivot)
}

#WOW plot classic (week over week plot non-interactive)
#Where dataT = propDF (preferred propDF.s, which filters for data pertaining to just the 2nd last forecast week)
#To change forecast week visualized -> change propDF.s filter week
#Note - this parameter is the same as WOWplotly()
WOWplot <- function(dataT, byday=T){
  if(byday){
    lab4y = "Selection coefficient (Daily)"
    dataT[,c("week.coef", "se", "CI.week.2.5..", "CI.week.97.5..")]=dataT[,c("week.coef", "se", "CI.week.2.5..", "CI.week.97.5..")]/7
  }else{
    lab4y = "Selection coefficient (Weekly)"
  }
  
  ## converting to predicted %
  dataT[,c("predicted", "std.error","conf.low","conf.high")]= dataT[,c("predicted", "std.error","conf.low","conf.high")]*100
  
  # generate plot
  weekp<-ggplot(dataT, 
                aes(x=predicted, y=week.coef, color=variant))+
    geom_point()+
    geom_text_repel(data=dataT, aes(x=predicted, y=week.coef,
                                    label=variant, size=5),
                    show.legend = F,
                    force=0.5)+ #, max.overlaps = Inf
    geom_errorbarh(aes(xmin=conf.low, xmax=conf.high, colour=variant))+
    geom_errorbar(aes(ymin=CI.week.2.5.., ymax=CI.week.97.5.., colour=variant))+
    xlab(paste("Predicted proportion (%) -", dataT$weekDate[1]))+
    ylab(lab4y) +
    scale_x_log10(labels = label_comma())+
    theme(legend.position="none",
          panel.grid = element_line(colour="grey"),
          panel.background = element_rect(fill="white"),
          axis.title = element_text(size=15, face="bold"),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(paste0(toestim, " filter - ", length(unique(dataT$variant)), " lineages"))
  return(weekp)
}

#Table of coefficients
#Where ctable = table of coefficients (i.e) ctalbe = CoefTable)
CoeffDT <- function(ctable){
  return(
    DT::datatable(data.frame(ctable),
                  caption=HTML("<b><center><font size='5'> \n The selected reference is: ", 
                               ref2use,
                               "</b></center></font>"),
                  filter='top',
                  extensions="Buttons",
                  options=list(pageLength=100,
                               columnDefs = list(list(className = 'dt-center', targets = "_all")),
                               buttons=c("copy","csv"),
                               dom="Blfrtip")) %>%
      formatStyle(columns=colnames(.), fontSize="10%")
  )
}


#logit plots from the Nowcast model
logitplot <- function(ncOutput, tosave=FALSE){ #where ncOutput = variable saved for runNowcast()
  propDF0 <- ncOutput
  preds <- dcast(propDF0[,c("weekDate","variant","predicted")], weekDate~variant, value.var="predicted") 
  
  propDF0 <- na.omit(propDF0) # to get weekly counts only
  freqs <- dcast(propDF0[,c("weekDate","variant","Freq")], weekDate~variant, value.var="Freq") 
  
  # get names of all variants modeled
  varis <- levels(ncOutput$variant)
  
  #Initialize empty list
  loplot <- vector("list", length=length(varis)-1); names(loplot) <- varis[2:length(varis)] #start at #2 to exclude ref
  
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
      ## may need to adjust the lim for extremes but this is already 1000 times
      #ylim(-3,3)+
      ylab(paste("log10(",varis[i],"/",varis[1],")"))+
      theme(legend.position="none",
            axis.text.x=element_text(angle=45, vjust=1, hjust=1))
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
  return(loplot)
}


#######################################################

#Plotly & Datatable plotting functions (interactive)

#Trend plotly
#Where trendplot = whatever trendplot is saved as (ie) p1 <- trendplot())
trendplotly <- function(trendplot){
  defaultW <- getOption("warn") 
  options(warn = -1) #hide warnings - need this specialized method for plotly objects
  p1x <- ggplotly(trendplot, tooltip=c("group", "date", "weekDate", "y", "ymax", "ymin")) %>% 
    layout(dragmode=F, 
           legend=list(title=list(text="<b>Lineage</b>"))
    )
  #Changing labels
  #B/c multiple ggplot objects, can't modify text labels - need to manually change
  for(i in 1:length(p1x$x$data)){
    p1x$x$data[[i]]$text <- gsub(p1x$x$data[[i]]$text, pattern="weekDate", replacement="Week")
    p1x$x$data[[i]]$text <- gsub(p1x$x$data[[i]]$text, pattern="date", replacement="Date")
    p1x$x$data[[i]]$text <- gsub(p1x$x$data[[i]]$text, pattern="round\\(predicted, 2\\)", replacement="Proportion Predicted")
    p1x$x$data[[i]]$text <- gsub(p1x$x$data[[i]]$text, pattern="round\\(prop, 2\\)", replacement="Proportion Actual")
    p1x$x$data[[i]]$text <- gsub(p1x$x$data[[i]]$text, pattern="round\\(conf\\.low, 2\\)", replacement="Lower CI")
    p1x$x$data[[i]]$text <- gsub(p1x$x$data[[i]]$text, pattern="round\\(conf\\.high, 2\\)", replacement="Upper CI")
    p1x$x$data[[i]]$text <- gsub(p1x$x$data[[i]]$text, pattern="sublineage_grouped", replacement="Variant")
  }
  return(p1x)
  options(warn = defaultW)#reset to allow warnings
}

#Stackedbar plotly
#Where ggstackedactual = ggplot of stackedbarplot(actualorpred=="prop") [ie) ggstackedactual = stackedbar_actual]
#Where ggstackedpredicted = ggplot of stackedbarplot(actualorpred=="predicted") [ie) ggstackedpredicted = stackedbar_predicted]
stackedbarplotly <- function(ggstackedactual, ggstackedpredicted){
  defaultW <- getOption("warn") 
  stackedptly <- subplot(ggplotly(ggstackedactual, tooltip="text"), style(ggplotly(ggstackedpredicted, tooltip="text"), showlegend=F), 
                         nrows=2, shareX=T, titleY=T, shareY=F) %>% 
    layout(dragmode=F, 
           plot_bgcolor='rgb(255,255,255)',
           legend=list(title=list(text="<b>Lineage</b>")),
           annotations = list(list(y=1.0, text="Actual", showarrow=F, xref="paper", yref="paper"),
                              list(y=0.45, text="Nowcast", showarrow=F, xref="paper", yref="paper"),
                              list(x=-0.07, y=0.5, text="<b>Proportion (%)</b>", font=list(size=20),
                                   textangle=270, showarrow=F, xref="paper", yref="paper"))
    )
  return(stackedptly)
  options(warn = defaultW)
}

