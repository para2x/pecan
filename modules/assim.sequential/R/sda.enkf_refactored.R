sda.enkf.refactored <- function(settings,
                                obs.mean,
                                obs.cov,
                                IC = NULL,
                                Q = NULL,
                                adjustment = TRUE,
                                control=list(trace=T,interactivePlot=T)) {
  #My personal notes -----------------------
  # Three analysis function was initially developed into this:
  # 1-EnKF
  # 2-Generalized Ensubmle Filter: -tobit / -wish
  #------------- Some important variables  
  # muf/Pf - forcast mean and covariance
  # Y/R    - Observed data and covariance
  # mu.a/Pa  - afetr analysis - new mean and covariance
  # nt is the length of observed  
  # When processvar == FALSE it means we are doin EnKF and it's TRUE Generlized Ensumble Filter
  # Generlized Ensumble Filter NEEDS process variance to avoid filter divergence and it does not
  # have analytical solution - needs MCMC
  # X stores IC of state variables and then collects state variables in each loop
  # Y stores the observed mean
  #------------------------  
  ymd_hms <- lubridate::ymd_hms
  hms     <- lubridate::hms
  second  <- lubridate::second
  ###-------------------------------------------------------------------###
  ### read settings                                                     ###
  ###-------------------------------------------------------------------### 
  model      <- settings$model$type
  write      <- settings$database$bety$write
  defaults   <- settings$pfts
  outdir     <- settings$modeloutdir # currently model runs locally, this will change if remote is enabled
  rundir     <- settings$host$rundir
  host       <- settings$host
  forecast.time.step <- settings$state.data.assimilation$forecast.time.step  #idea for later generalizing
  nens       <- as.numeric(settings$state.data.assimilation$n.ensemble)
  processvar <- settings$state.data.assimilation$process.variance
  sample_parameters <- settings$state.data.assimilation$sample.parameters
  var.names <- unlist(sapply(settings$state.data.assimilation$state.variable, 
                             function(x) {
                               x$variable.name
                             }, 
                             USE.NAMES = FALSE), 
                      use.names = FALSE)
  names(var.names) <- NULL
  #filtering obs data based on years specifited in setting > state.data.assimilation
  assimyears<- year(settings$state.data.assimilation$start.date):year(settings$state.data.assimilation$end.date) # years that assimilations will be done for - obs will be subsetted based on this
  obs.mean<-obs.mean[sapply(year(names(obs.mean)),function(obs.year) obs.year%in%(assimyears))]
  obs.cov<-obs.cov[sapply(year(names(obs.cov)),function(obs.year) obs.year%in%(assimyears))]
  # dir address based on the end date
  if(!dir.exists("SDA")) dir.create("SDA",showWarnings = F)
  ###-------------------------------------------------------------------###
  ### tests before data assimilation                                    ###
  ###-------------------------------------------------------------------###  
  
  # at some point add a lot of error checking 
  # read time from data if data is missing you still need
  # to have NAs or NULL with date name vector to read the correct netcdfs by read_restart
  
  obs.times <- names(obs.mean)
  obs.times.POSIX <- ymd_hms(obs.times)

  for (i in seq_along(obs.times)) {
    if (is.na(obs.times.POSIX[i])) {
      if (is.na(lubridate::ymd(obs.times[i]))) {
        print("Error: no dates associated with observations")
      } else {
        ### Data does not have time associated with dates 
        ### Adding 12:59:59PM assuming next time step starts one second later
        print("Pumpkin Warning: adding one minute before midnight time assumption to dates associated with data")
        obs.times.POSIX[i] <- ymd_hms(paste(obs.times[i], "23:59:59"))
        # if(nchar(year(obs.times.POSIX[i]))==3){
        #   #TODO: BROKEN: need to add leading zeros to years with less than 4 digits
        #   obs.times.POSIX[i] <- paste0('0',ymd_hms(paste(obs.times[i], "23:59:59")))
        # } 
      }
    }
  }
  obs.times <- obs.times.POSIX
  
  # need explicit forecast length variable in settings start time, stop time, restart time if
  # restart time is not provided restart in stop time
  
  ###-------------------------------------------------------------------###
  ### set up for data assimilation                                      ###-----
  ###-------------------------------------------------------------------###  
  
  nt          <- length(obs.times)
  FORECAST    <- ANALYSIS <- list()
  enkf.params <- list()
  aqq         <- NULL
  bqq         <- numeric(nt + 1)
  CI.X1       <- matrix(0, 3, nt)
  CI.X2       <- CI.X1
  q.bar        <- NULL #default process covariance matrix
  
  ##### Creating matrices that describe the bounds of the state variables
  ##### interval is remade everytime depending on the data at time t
  ##### state.interval stays constant and converts new.analysis to be within the correct bounds
  interval    <- NULL
  state.interval <- cbind(as.numeric(lapply(settings$state.data.assimilation$state.variables,'[[','min_value')),
                          as.numeric(lapply(settings$state.data.assimilation$state.variables,'[[','max_value')))
  rownames(state.interval) <- var.names
  
  wish.df <- function(Om, X, i, j, col) {
    (Om[i, j]^2 + Om[i, i] * Om[j, j]) / var(X[, col])
  }
  
  if(var.names=="Fcomp"){
    y_star_create<-y_star_create_Fcomp
  }

  # weight matrix
  wt.mat <- matrix(NA, nrow = nens, ncol = nt)
  
  ###-------------------------------------------------------------------###
  ### loop over time                                                    ###----
  ###-------------------------------------------------------------------### 
  
  for(t in seq_len(nt)){#) {
   
    # do we have obs for this time - what year is it ?
    obs <- which(!is.na(obs.mean[[t]]))
    obs.year<-year(names(obs.mean)[t])
    cat("Start ",obs.year,"====================================================")
    #- genereating the ensumbles
    if (t>1){
      #print(run.id)
      restart.arg<-list(runid = run.id, 
                                start.time = strptime(obs.times[t-1],format="%Y-%m-%d %H:%M:%S"),
                                stop.time = strptime(obs.times[t],format="%Y-%m-%d %H:%M:%S"), 
                                settings = settings,
                                new.state = new.state, 
                                new.params = new.params, 
                                inputs = inputs, 
                                RENAME = TRUE)
    }else{
      restart.arg<-NULL
    }
    #-- main function for making ensembles
    ensemble.gen(settings,
                 settings$state.data.assimilation$n.ensemble,
                 restart=restart.arg,
                 #years=year(settings$state.data.assimilation$start.date):year(settings$state.data.assimilation$end.date)
                 )->ens.outs
    #collectting the basics of how the simulations were run: runid/inputs and simulation outputs
    ens.outputs<-ens.outs[[1]] # this is just the variables what would be left out is the params
    params<-ens.outs[[2]]
    run.id<-ens.outs[[3]]
    inputs<-ens.outs[[4]]
    #print(run.id)

    #--- this could be expanded to find the exact date in ens.outputs
    X<-lapply(ens.outputs,function(lis){
      #take out the year of observation
      (lis[[obs.year%>%as.character()]])[[1]]%>%unlist()
      })

    X <- do.call(rbind, X)

  
    FORECAST[[t]] <- X
    mu.f <- as.numeric(apply(X, 2, mean, na.rm = TRUE))
    Pf <- cov(X)
    pmiss <- which(diag(Pf) == 0)
    diag(Pf)[pmiss] <- 0.1 ## hack for zero variance

  
    ###-------------------------------------------------------------------###
    ###  preparing OBS                                                   ###----
    ###-------------------------------------------------------------------###  
    if (any(obs)) {

      choose <- na.omit(charmatch(colnames(X),names(obs.mean[[t]])))
      
      Y <- unlist(obs.mean[[t]][choose])
      Y[is.na(Y)] <- 0 
      
      R <- as.matrix(obs.cov[[t]][choose,choose])
      R[is.na(R)]<-0
      
      if (length(obs.mean[[t]]) > 1) {
        diag(R)[which(diag(R)==0)] <- min(diag(R)[which(diag(R) != 0)])/2
        diag(Pf)[which(diag(Pf)==0)] <- min(diag(Pf)[which(diag(Pf) != 0)])/5
      }
      
      ### TO DO: plotting not going to work because of observation operator i.e. y and x are on different scales
      ###-------------------------------------------------------------------###
      ### Analysis                                                          ###----
      ###-------------------------------------------------------------------###
      #method=EnKF/GNF

      if(processvar == FALSE){an.method<-"EnKF"  }else{    an.method<-"GEF"   }  
      #-analysis function
        enkf.params[[t]] <-Analysis.sda(settings,method=an.method, Forcast=list(Pf=Pf,mu.f=mu.f,Q=Q,X=X),
                                                                   Observed=list(R=R,Y=Y),
                                                                   choose)
    

      Pa<- enkf.params[[t]]$Pa
      mu.a<- enkf.params[[t]]$mu.a
      #-- writing Trace--------------------
      if(control$trace) {
        cat ("\n Obs data in this iteration *************** \n")
        print(Y)
        print(R)
        cat ("\n Forcast and analysis output *************** \n")
        print(enkf.params)
      }
      
      } else {
      ###-------------------------------------------------------------------###
      ### No Observations -- Starts Here                                    ###----
      ###-------------------------------------------------------------------### 
      ### no process variance -- forecast is the same as the analysis ###
      if (processvar==FALSE) {
        mu.a <- mu.f
        Pa   <- Pf + Q
        ### yes process variance -- no data
      } else {
        mu.a <- mu.f
        if(is.null(q.bar)){
          q.bar <- diag(ncol(X))
          print('Process variance not estimated. Analysis has been given uninformative process variance')
        } 
        Pa   <- Pf + solve(q.bar)
      }
      enkf.params[[t]] <- list(mu.f = mu.f, Pf = Pf, mu.a = mu.a, Pa = Pa)
    }
    ###-------------------------------------------------------------------###
    ### adjustement/update state matrix                                   ###----
    ###-------------------------------------------------------------------### 
    if(adjustment == TRUE){
      S_f  <- svd(Pf)
      L_f  <- S_f$d
      V_f  <- S_f$v
      
      ## normalize
      Z <- X*0
      
      for(i in seq_len(nrow(X))){
        if(processvar == TRUE) {
          Z[i,] <- 1/sqrt(L_f) * t(V_f)%*%(X.new[i,]-mu.f)
        }else{
          Z[i,] <- 1/sqrt(L_f) * t(V_f)%*%(X[i,]-mu.f)
        }
      }
      Z[is.na(Z)]<-0
      
      ## analysis
      S_a  <- svd(Pa)

      L_a  <- S_a$d
      V_a  <- S_a$v

      ## analysis ensemble
      X_a <- X*0
      for(i in seq_len(nrow(X))){

        X_a[i,] <- V_a %*%diag(sqrt(L_a))%*%Z[i,] + mu.a
      }
      
      # # calculate likelihoods
#      for(i in seq_len(nens)){
#        wt.mat[i,t]<-dmnorm_chol(FORECAST[[t]][i,], mu.a, solve(Pa), log = TRUE)
#      }
      
      if(sum(mu.a - colMeans(X_a)) > 1 | sum(mu.a - colMeans(X_a)) < -1) logger.warn('Problem with ensemble adjustment (1)')
      if(sum(diag(Pa) - diag(cov(X_a))) > 5 | sum(diag(Pa) - diag(cov(X_a))) < -5) logger.warn('Problem with ensemble adjustment (2)')
      
      analysis <- as.data.frame(X_a)
    }else{
 
      analysis <- as.data.frame(rmvnorm(as.numeric(nrow(X)), mu.a, Pa, method = "svd"))
    }
    
    colnames(analysis) <- colnames(X)

    ##### Mapping analysis vectors to be in bounds of state variables
    if(processvar==TRUE){
      for(i in 1:ncol(analysis)){
        int.save <- state.interval[which(startsWith(colnames(analysis)[i],
                                                    var.names)),]
        analysis[analysis[,i] < int.save[1],i] <- int.save[1]
        analysis[analysis[,i] > int.save[2],i] <- int.save[2]
      }
    }
    
    ## in the future will have to be separated from analysis
    new.state  <- analysis
    new.params <- params
    ANALYSIS[[t]] <- analysis
    ### Interactive plotting ------------------------------------------------------   
    if (t > 1 & control$interactivePlot) { #
      interactive.plotting.sda(settings,t,obs.times,obs.mean,obs.cov,obs,X,FORECAST,ANALYSIS)
    }
    ###-------------------------------------------------------------------###
    ### save outputs                                                      ###----
    ###-------------------------------------------------------------------### 
    save(t, FORECAST, ANALYSIS, enkf.params, file = file.path(settings$outdir,"SDA", "sda.output.Rdata"))
  } ### end loop over time
  ###-------------------------------------------
  ###-------------------------------------------------------------------###
  ### create diagnostics                                                ###-----
  ###-------------------------------------------------------------------### 
  
  ### LOAD CLIMATE ### HACK ### LINKAGES SPECIFIC----
  if (model == "LINKAGES") {
    climate_file <- settings$run$inputs$met$path
    load(climate_file)
    temp.mat     <- temp.mat[year(obs.times) - 853, ]
    precip.mat   <- precip.mat[year(obs.times) - 853, ]
  } else {
    print("climate diagnostics under development")
  }
  
  ###-------------------------------------------------------------------###
  ### time series plots                                                 ###-----
  ###-------------------------------------------------------------------### 
  postana.timeser.plotting.sda(settings,t,obs.times,obs.mean,obs.cov,obs,X,FORECAST,ANALYSIS)
  ###-------------------------------------------------------------------###
  ### bias diagnostics                                                  ###----
  ###-------------------------------------------------------------------###
  postana.bias.plotting.sda(settings,t,obs.times,obs.mean,obs.cov,obs,X,FORECAST,ANALYSIS)
  ###-------------------------------------------------------------------###
  ### process variance plots                                            ###-----
  ###-------------------------------------------------------------------### 
  if (processvar) postana.bias.plotting.sda(t,obs.times,X,aqq,bqq)
    
  

  
} # sda.enkf

