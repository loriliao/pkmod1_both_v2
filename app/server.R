# Server inputs: Dose, dosing regimen, dosing frequency,
# dosing cycle definition, number of dosing cycles

# read objects from "rx_shiny_data.rda" in the  AppDir folder,
# objects include: pkMod1, pkMod1_params, inits

load("rx_shiny_data.rda")

######################Simulate with variability######################
# Perform simulation
sim_params<-function(nsub,ppkparams){ 
  sigma<-matrix(c(0.0798,0.0628,0.0628,0.0661),2,2)
  mv<-mvrnorm(n=nsub,rep(0,2),sigma)
  CL<-ppkparams["CL"]*exp(mv[,1])
  V1<-ppkparams["V1"]*exp(mv[,2])
  
  sigma1<-matrix(c(0.000179),1,1)
  mv1<-mvrnorm(n=nsub,rep(0,1),sigma1)
  V2<-ppkparams["V2"]*exp(mv1[,1])
  Q<-ppkparams["Q"]
  
  sigma2<-matrix(c(0.0629,0.0505,0.0505,0.0521),2,2)
  mv2<-mvrnorm(n=nsub,rep(0,2),sigma2)
  alpha<-ppkparams["alpha"]*exp(mv2[,1])
  beta <-ppkparams["beta"]*exp(mv2[,2])+1   
  
  params.all<- data.frame(ID=1:200,CL=CL,V1=V1,Q=Q,V2=V2,alpha=alpha,beta=beta, row.names = NULL)
  return(params.all)
}

get.PK <- function(ds1a, ds1b, ds2, ds3, lday, dfreq1, dfreq2, bidfreq, params.mod) {    
  time<-sort(c(0:((lday+5)*24))) 
  res <- NULL #Create an empty matrix for storing results
  et_all<-NULL
  
  dose1a<-ds1a*1000  #convert mg to mcg
  if (is.na(ds1b)) ds1b=0 
  dose1b<-ds1b*1000  #convert mg to mcg
  dose2<-ds2*1000  #convert mg to mcg
  dose3<-ds3*1000  #convert mg to mcg
  bidfreq<-as.numeric(bidfreq)
  
  #define time-varying parameters
  params.mod<-params.mod %>%
    mutate(CL_1a=CL*(ds1a/300)^0.401,
           CL_1b=CL*(ds1b/300)^0.401,
           CL_2=CL*(ds2/300)^0.401,
           CL_3=CL*(ds3/300)^0.401,
           V1_1a=V1*(ds1a/300)^0.252,
           V1_1b=V1*(ds1b/300)^0.252,
           V1_2=V1*(ds2/300)^0.252,
           V1_3=V1*(ds3/300)^0.252)%>%
    mutate(CL_1b=ifelse(CL_1b==0,CL_1a,CL_1b),
           V1_1b=ifelse(V1_1b==0,V1_1a,V1_1b))

  for (i in 1:nrow(params.mod)){
    ## loading doses
    dosing <- et(timeUnits="hr") %>%
      add.dosing(dose=dose1a,start.time=0,nbr.doses=1)
    if(dose1b>0) dosing <- dosing %>%
      add.dosing(dose=dose1b,start.time=bidfreq,nbr.doses=1)
    
    dosing <- dosing %>%
      add.dosing(dose=dose2,start.time=24,nbr.doses=1)
    if(dfreq1==2) dosing <- dosing %>%
      add.dosing(dose=dose2,start.time=24+bidfreq,nbr.doses=1)
    
    ## maintenance dose
    if(dfreq2==1) dosing <- dosing %>%
      add.dosing(dose=dose3,start.time=48,dosing.interval=24,nbr.doses=(lday-2)) else 
        if(dfreq2==2) dosing <- dosing %>%
      add.dosing(dose=dose3,start.time=48,dosing.interval=24,nbr.doses=(lday-2)*2)%>%
      add.dosing(dose=dose3,start.time=48+bidfreq,dosing.interval=24,nbr.doses=(lday-2)*2)

    ##Event table
    et <- dosing  %>%
      add.sampling(set_units(time,hours))
    
    #add time-varying parameters
    et<- et %>% 
      mutate(CL=ifelse(time < set_units(12, h)|(time==set_units(12, h)&evid==0), params.mod[i,"CL_1a"], 
                       ifelse(time < set_units(24, h)|(time==set_units(24, h)&evid==0), params.mod[i,"CL_1b"],
                              ifelse(time < set_units(48, h)|(time==set_units(48, h)&evid==0),params.mod[i,"CL_2"],
                                     params.mod[i,"CL_3"]))),
             V1=ifelse(time < set_units(12, h)|(time==set_units(12, h)&evid==0), params.mod[i,"V1_1a"], 
                       ifelse(time < set_units(24, h)|(time==set_units(24, h)&evid==0), params.mod[i,"V1_1b"],
                              ifelse(time < set_units(48, h)|(time==set_units(48, h)&evid==0),params.mod[i,"V1_2"],
                                     params.mod[i,"V1_3"])))
      )
    
    x<-rxSolve(pkMod1,params.mod[i,] %>% dplyr::select(-c(CL,V1,CL_1a,CL_1b,CL_2,CL_3,V1_1a,V1_1b,V1_2,V1_3)),
               et,inits,method = "lsoda")
    res<-cbind(res,x[,"C1"])
    et_all<-rbind(et_all,data.frame(pat=i,et,params.mod[i,] %>%
                                      dplyr::select(-c(CL,V1,CL_1a,CL_1b,CL_2,CL_3,V1_1a,V1_1b,V1_2,V1_3)),
                                    row.names=NULL))
  }
  res<- data.frame(time=time+24,res)  #Correct time to start from Day 1
  write.csv(et_all,"et_all.csv",row.names=F)
  
  return(res)
}

sum.PK<-function(ds){
  res.q.t <- apply(ds, 1, quantile, prob = c(.05, .5, .95))
  res.pk<-data.frame(time=ds$time,t(res.q.t))   #convert hr to day
  names(res.pk)<-c("Time, d","5th Pctl","Median","95th Pctl")
  return(res.pk)
}

calc_auc<-function(data)
{
  auccalc<-data
  auccalc$h=auccalc$time*auccalc$PK
  auccalc$lagH = Lag(auccalc$h,shift=1)
  auccalc$lagTime = Lag(auccalc$time,shift=1)
  auccalc$lagConc = Lag(auccalc$PK,shift=1)
  auccalc$trapAuc=(auccalc$time-auccalc$lagTime)*(auccalc$PK+auccalc$lagConc)/2
  return(sum(auccalc$trapAuc,na.rm=T))
}

get.PKparams<-function(ds,start,end){
  doseint<-ds %>% filter(time>=start,time<=end) %>%
    gather(ID,PK,-time) %>%
    mutate(ID=as.numeric(str_remove(ID,"X")))
  auc1<-NULL
  for (i in unique(doseint$ID)){
    auc1<-rbind(auc1,data.frame(ID=i,
                                AUC=calc_auc(doseint %>%
                                               filter(ID==i))))
  }
  cmax_min<-inner_join(doseint %>%
                         group_by(ID) %>%
                         filter(time==end) %>%
                         dplyr::select(-time) %>%
                         rename(Cmin=PK),
                       doseint %>%
                         group_by(ID) %>%
                         reframe(Cmax=max(PK)))
  pkparams<-inner_join(cmax_min,auc1) %>%
    mutate(intStart=start/24, intEnd=end/24)
  return(pkparams)
}

sum.PKparams<-function(ds,start,end){
  sumPKparams<-data.frame(STAT=c("5th Pctl","Median","95th Pctl"),
                          ds %>% filter(intStart*24==start) %>%
                            group_by(NULL) %>% 
                            reframe(Cmin=quantile(Cmin,p=c(.05,.5,.95)),
                                    Cmax=quantile(Cmax,p=c(.05,.5,.95)),
                                    AUC=quantile(AUC,p=c(.05,.5,.95)))) %>%
    mutate(intStart=start/24, intEnd=end/24)
  names(sumPKparams)<-c("Statistic","Cmin, ng/mL","Cmax, ng/mL",
                           "AUC(24), hr*ng/mL","Interval start, d", "Interval end, d")
  return(sumPKparams)
}

# Define server logic 
shinyServer(function(input, output) {
  # disable the downdload button on page load
  shinyjs::disable("downloadCSV")
  shinyjs::disable("downloadPlot")
  
  ds1a <- reactive({input$dose1a}) 
  ds1b <- reactive({input$dose1b}) 
  ds2 <- reactive({input$dose2})
  ds3 <- reactive({input$dose3}) 
  dfreq1 <- reactive({input$dfreq1}) 
  dfreq2 <- reactive({input$dfreq2})
  bidfreq <- reactive({input$bidfreq})
  lday <- eventReactive(input$submit_mod, {input$ldoseDay})   #Maintenance dose up thru this day
  day1 <- reactive({input$day1})  #derive PK params for this day
  day2 <- reactive({input$day2})  #derive PK params for this day 
  semilog<-reactive({switch(input$semilog, "linear"=F, "semi-log"=T)})
  
  sim.params <- eventReactive(input$submit_mod, {
    pkparams<-sim_params(200,pkMod1_params) 
    return(pkparams)
  })
  
  get.cp <- eventReactive(input$submit_mod, {   
    get.PK(ds1a(), ds1b(), ds2(), ds3(), lday(), dfreq1(), dfreq2(), bidfreq(), sim.params())
  })
  
  plotInput <- reactive({
    res<-get.cp()
    res.q.t <- apply(res, 1, quantile, prob = c(.05, .5, .95))
    res.pk<-data.frame(time=res$time/24,t(res.q.t))   #convert hr to day
    
    # write.csv(res.q.t,"simPK.csv",row.names=F)
    p1<-ggplot(res.pk, aes(time,X50.))+
      labs(x="Time, d",y="Conc, ng/mL") +
      geom_line() + 
      geom_ribbon(aes(ymin=X5.,ymax=X95.),alpha=.2)+
      scale_x_continuous(breaks=c(1:(lday()+6)))+
      theme_bw()+
      theme(panel.grid.minor = element_blank())
    
    if (semilog()) {
      p1<-p1+
        scale_y_log10()
    }
    return (p1)
  })
  
  output$CpPlot <- renderPlot({
    print(plotInput())
  })
  
  datasetInput <- reactive({
    df_PK<-get.cp()
    df_sampPK<-df_PK[,1:11] %>% mutate(time=time/24)
    names(df_sampPK)<-c("Time, d",paste("Subj",1:10))
    df_sumPK<-sum.PK(df_PK %>% mutate(time=time/24))
    
    start1<-day1()*24
    end1<-(day1()+1)*24
    start2<-day2()*24
    end2<-(day2()+1)*24
    if (!is.na(start1)&!is.na(end1)) {
      df1<-get.PKparams(get.cp(),start1,end1)
      df3<-sum.PKparams(df1,start1,end1)
    } else {
      df1<-NULL
      df3<-NULL
    }
    if (!is.na(start2)&!is.na(end2)) {
      df2<-get.PKparams(get.cp(),start2,end2)
      df4<-sum.PKparams(df2,start2,end2)
    } else {
      df2<-NULL
      df4<-NULL
    }
    df_indPKparams<-rbind(df1,
                          df2)
    df_sumPKparams<-rbind(df3,
                            df4)
      
    switch(input$dataset,
           "PK (Sample of 10 Subjects)" = df_sampPK,
           "PK parameters" = df_indPKparams,
           "Summary of PK" = df_sumPK,
           "Summary of PK parameters" = df_sumPKparams)
  })
  
  # Table of selected dataset ----
  output$table <- renderTable({
    if (!is.null(datasetInput())) datasetInput() else NULL
  })
  
  #Save all output to zip folder
  observe({
    if (input$submit_mod > 0) {
      # enable the download button
      shinyjs::enable("downloadCSV")
      shinyjs::enable("downloadPlot")
    }
  })
  
  to_download <- reactiveValues()
  observe({
    start1<-day1()*24
    end1<-(day1()+1)*24
    start2<-day2()*24
    end2<-(day2()+1)*24
    
    to_download$params<-sim.params()
    to_download$statPK<-sum.PK(get.cp() %>% mutate(time=time/24))
    to_download$simPK<-get.cp() %>% mutate(time=time/24)
    
    if (!is.na(start1)&!is.na(end1)) {
      df1<-get.PKparams(get.cp(),start1,end1)
      df3<-sum.PKparams(df1,start1,end1)
    } else {
      df1<-NULL
      df3<-NULL
    }
    if (!is.na(start2)&!is.na(end2)) {
      df2<-get.PKparams(get.cp(),start2,end2)
      df4<-sum.PKparams(df2,start2,end2)
    } else {
      df2<-NULL
      df4<-NULL
    }
    to_download$PKparams<-rbind(df1,
                          df2)
    to_download$statPKparams<-rbind(df3,
                          df4)
    
  })
  output$downloadCSV <- downloadHandler(
    filename = function(){
      paste("simData", ".zip", sep = "")
    },
    content = function(file){
      temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
      dir.create(temp_directory)
      
      reactiveValuesToList(to_download) %>%
        imap(function(x,y){
          if(!is.null(x)){
            file_name <- glue("{y}.csv")
            readr::write_csv(x, file.path(temp_directory, file_name))
          }
        })
      
      zip::zip(
        zipfile = file,
        files = dir(temp_directory),
        root = temp_directory
      )
    },
    contentType = "application/zip"
  )
  
  output$downloadPlot <- downloadHandler(
    filename = function() { "simPK plot.png" },
    content = function(file) {
      ggsave(file, plot = plotInput(), device = "png")
    }
  )
})

