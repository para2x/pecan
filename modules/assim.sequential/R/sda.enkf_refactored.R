sda.enkf.refactored <- function(settings, obs.mean, obs.cov, IC = NULL, Q = NULL, adjustment = TRUE, restart=NULL) {
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

  t1         <- 1
  pink       <- col2rgb("deeppink")
  alphapink  <- rgb(pink[1], pink[2], pink[3], 180, max = 255)
  green      <- col2rgb("green")
  alphagreen <- rgb(green[1], green[2], green[3], 75, max = 255)
  blue       <- col2rgb("blue")
  alphablue  <- rgb(blue[1], blue[2], blue[3], 75, max = 255)
  purple       <- col2rgb("purple")
  alphapurple <- rgb(purple[1], purple[2], purple[3], 75, max = 255)
  brown       <- col2rgb("brown")
  alphabrown <- rgb(brown[1], brown[2], brown[3], 75, max = 255)
  
  # weight matrix
  wt.mat <- matrix(NA, nrow = nens, ncol = nt)
  
  ###-------------------------------------------------------------------###
  ### loop over time                                                    ###----
  ###-------------------------------------------------------------------### 
  
  for(t in seq_len(nt)) {
    PEcAn.assim.sequential:::ensemble.gen(settings,settings$state.data.assimilation$n.ensemble)->ens.outputs
    
    
    X <- do.call(rbind, X)

    FORECAST[[t]] <- X
    
    obs <- which(!is.na(obs.mean[[t]]))
    
    mu.f <- as.numeric(apply(X, 2, mean, na.rm = TRUE))
    Pf <- cov(X)
    pmiss <- which(diag(Pf) == 0)
    diag(Pf)[pmiss] <- 0.1 ## hack for zero variance
    
    ###-------------------------------------------------------------------###
    ### analysis                                                          ###----
    ###-------------------------------------------------------------------###  
    if (any(obs)) {
      # if no observations skip analysis
      # choose <- na.omit(charmatch(
      #   na.omit(unlist(lapply(strsplit(colnames(X),
      #                                  split = var.names),
      #                         function(x) x[2]))),
      #   na.omit(unlist(lapply(strsplit(names(obs.mean[[t]]),
      #                                  split = var.names), #TO DO don't hardcode this
      #                         function(x) x[2]))))) #matches y to model
      # 
      
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
      enkf.params[[t]] <- ifelse(processvar == FALSE,
                                  Analysis.sda(method="EnKF"),
                                  Analysis.sda(method="GEF")
                                  ) 
      
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
    
    ANALYSIS[[t]] <- analysis
    ### Interactive plotting ------------------------------------------------------   
    if (interactive() & t > 1) { #
     
      t1 <- 1
      names.y <- unique(unlist(lapply(obs.mean[t1:t], function(x) { names(x) })))
      Ybar <- t(sapply(obs.mean[t1:t], function(x) {
        tmp <- rep(NA, length(names.y))
        names(tmp) <- names.y
        mch <- match(names(x), names.y)
        tmp[mch] <- x[mch]
        tmp
      }))
      
      if(any(obs)){
        Y.order <- na.omit(pmatch(colnames(X), colnames(Ybar)))
        Ybar <- Ybar[,Y.order]
        Ybar[is.na(Ybar)] <- 0
        YCI <- t(as.matrix(sapply(obs.cov[t1:t], function(x) {
          if (length(x)<2) {
            rep(NA, length(names.y))
          }
          sqrt(diag(x))
        })))
        
        YCI <- YCI[,Y.order]
        YCI[is.na(YCI)] <- 0

      }else{
        YCI <- matrix(NA,nrow=length(t1:t), ncol=max(length(names.y),1))
      }
      
      par(mfrow = c(2, 1))
      for (i in 1:ncol(FORECAST[[t]])) { #
        
        Xbar <- plyr::laply(FORECAST[t1:t], function(x) { mean(x[, i]/rowSums(x[,1:9]), na.rm = TRUE) })
        Xci  <- plyr::laply(FORECAST[t1:t], function(x) { quantile(x[, i]/rowSums(x[,1:9]), c(0.025, 0.975), na.rm = TRUE) })
        
        Xa <- plyr::laply(ANALYSIS[t1:t], function(x) { mean(x[, i]/rowSums(x[,1:9]), na.rm = TRUE) })
        XaCI <- plyr::laply(ANALYSIS[t1:t], function(x) { quantile(x[, i]/rowSums(x[,1:9]), c(0.025, 0.975), na.rm = TRUE) })
        
        ylab.names <- unlist(sapply(settings$state.data.assimilation$state.variable, 
                                    function(x) { x })[2, ], use.names = FALSE)
        
        # observation / data
        if (i <= ncol(Ybar) & any(obs)) {
          plot(as.Date(obs.times[t1:t]), 
               Xbar, 
               ylim = range(c(XaCI, Xci, Ybar[,i]), na.rm = TRUE), 
               type = "n", 
               xlab = "Year", 
               ylab = ylab.names[grep(colnames(X)[i], var.names)], 
               main = colnames(X)[i])
          ciEnvelope(as.Date(obs.times[t1:t]),
                     as.numeric(Ybar[, i]) - as.numeric(YCI[, i]) * 1.96, 
                     as.numeric(Ybar[, i]) + as.numeric(YCI[, i]) * 1.96, 
                     col = alphagreen)
          lines(as.Date(obs.times[t1:t]), 
                as.numeric(Ybar[, i]), 
                type = "l", 
                col = "darkgreen", 
                lwd = 2)
        }else{
          plot(as.Date(obs.times[t1:t]), 
               Xbar, 
               ylim = range(c(XaCI, Xci), na.rm = TRUE), 
               type = "n", 
               xlab = "Year", 
               ylab = ylab.names[grep(colnames(X)[i], var.names)], 
               main = colnames(X)[i])
        }
        
        # forecast
        ciEnvelope(as.Date(obs.times[t1:t]), Xci[, 1], Xci[, 2], col = alphablue)  #col='lightblue')
        lines(as.Date(obs.times[t1:t]), Xbar, col = "darkblue", type = "l", lwd = 2)
        
        # analysis
        ciEnvelope(as.Date(obs.times[t1:t]), XaCI[, 1], XaCI[, 2], col = alphapink)
        lines(as.Date(obs.times[t1:t]), Xa, col = "black", lty = 2, lwd = 2)
        #legend('topright',c('Forecast','Data','Analysis'),col=c(alphablue,alphagreen,alphapink),lty=1,lwd=5)
      }
    }
    #dev.off()
    ###-------------------------------------------------------------------###
    ### forecast step                                                     ###----
    ###-------------------------------------------------------------------### 
    if (t < nt) {

      ###-------------------------------------------------------------------###
      ### Run model                                                         ###
      ###-------------------------------------------------------------------### 
    
    }
    
    ###-------------------------------------------------------------------###
    ### save outputs                                                      ###
    ###-------------------------------------------------------------------### 
    save(t, FORECAST, ANALYSIS, enkf.params, file = file.path(settings$outdir, "sda.output.Rdata"))

    
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
  
  pdf(file.path(settings$outdir, "sda.enkf.time-series.pdf"))
  
  names.y <- unique(unlist(lapply(obs.mean[t1:t], function(x) { names(x) })))
  Ybar <- t(sapply(obs.mean[t1:t], function(x) {
    tmp <- rep(NA, length(names.y))
    names(tmp) <- names.y
    mch <- match(names(x), names.y)
    tmp[mch] <- x[mch]
    tmp
  }))
  Y.order <- na.omit(pmatch(colnames(FORECAST[[t]]), colnames(Ybar)))
  Ybar <- Ybar[,Y.order]
  YCI <- t(as.matrix(sapply(obs.cov[t1:t], function(x) {
    if (is.null(x)) {
      rep(NA, length(names.y))
    }
    sqrt(diag(x))
  })))
  
  Ybar[is.na(Ybar)]<-0
  YCI[is.na(YCI)]<-0
  
  YCI <- YCI[,Y.order]
  Xsum <- plyr::laply(FORECAST, function(x) { mean(rowSums(x[,1:length(names.y)], na.rm = TRUE)) })[t1:t]
  Xasum <- plyr::laply(ANALYSIS, function(x) { mean(rowSums(x[,1:length(names.y)], na.rm = TRUE)) })[t1:t]
  
  for (i in seq_len(ncol(X))) {
    Xbar <- plyr::laply(FORECAST[t1:t], function(x) {
      mean(x[, i], na.rm = TRUE) }) #/rowSums(x[,1:9],na.rm = T)
    Xci <- plyr::laply(FORECAST[t1:t], function(x) { 
      quantile(x[, i], c(0.025, 0.975),na.rm = T) })
    
    Xci[is.na(Xci)]<-0
    
    Xbar <- Xbar
    Xci <- Xci
    
    Xa <- plyr::laply(ANALYSIS[t1:t], function(x) { 
      
      mean(x[, i],na.rm = T) })
    XaCI <- plyr::laply(ANALYSIS[t1:t], function(x) { 
      quantile(x[, i], c(0.025, 0.975),na.rm = T )})
    
    Xa <- Xa
    XaCI <- XaCI
    
    plot(as.Date(obs.times[t1:t]),
         Xbar, 
         ylim = range(c(XaCI, Xci), na.rm = TRUE),
         type = "n", 
         xlab = "Year", 
         ylab = ylab.names[grep(colnames(X)[i], var.names)],
         main = colnames(X)[i])
    
    # observation / data
    if (i<ncol(X)) { #
      ciEnvelope(as.Date(obs.times[t1:t]), 
                 as.numeric(Ybar[, i]) - as.numeric(YCI[, i]) * 1.96, 
                 as.numeric(Ybar[, i]) + as.numeric(YCI[, i]) * 1.96, 
                 col = alphagreen)
      lines(as.Date(obs.times[t1:t]), 
            as.numeric(Ybar[, i]), 
            type = "l", col = "darkgreen", lwd = 2)
    }
    
    # forecast
    ciEnvelope(as.Date(obs.times[t1:t]), Xci[, 1], Xci[, 2], col = alphablue)  #col='lightblue') #alphablue
    lines(as.Date(obs.times[t1:t]), Xbar, col = "darkblue", type = "l", lwd = 2) #"darkblue"
    
    # analysis
    ciEnvelope(as.Date(obs.times[t1:t]), XaCI[, 1], XaCI[, 2], col = alphapink) #alphapink
    lines(as.Date(obs.times[t1:t]), Xa, col = "black", lty = 2, lwd = 2) #"black"
    
    legend('topright',c('Forecast','Data','Analysis'),col=c(alphablue,alphagreen,alphapink),lty=1,lwd=5)
    
  }
  
  dev.off()
  ###-------------------------------------------------------------------###
  ### bias diagnostics                                                  ###----
  ###-------------------------------------------------------------------###
  pdf(file.path(settings$outdir, "bias.diagnostic.pdf"))
  for (i in seq_along(obs.mean[[1]])) {
    Xbar <- plyr::laply(FORECAST[t1:t], function(x) { mean(x[, i], na.rm = TRUE) })
    Xci <- plyr::laply(FORECAST[t1:t], function(x) { quantile(x[, i], c(0.025, 0.975)) })
    
    Xa <- plyr::laply(ANALYSIS[t1:t], function(x) { mean(x[, i], na.rm = TRUE) })
    XaCI <- plyr::laply(ANALYSIS[t1:t], function(x) { quantile(x[, i], c(0.025, 0.975)) })
    
    if(length(which(is.na(Ybar[,i])))>=length(t1:t)) next()
    reg <- lm(Xbar[t1:t] - unlist(Ybar[, i]) ~ c(t1:t))
    plot(t1:t, 
         Xbar - unlist(Ybar[, i]),
         pch = 16, cex = 1, 
         ylim = c(min(Xci[, 1] - unlist(Ybar[, i])), max(Xci[,2] - unlist(Ybar[, i]))), 
         xlab = "Time", 
         ylab = "Error", 
         main = paste(colnames(X)[i], " Error = Forecast - Data"))
    ciEnvelope(rev(t1:t), 
               rev(Xci[, 1] - unlist(Ybar[, i])), 
               rev(Xci[, 2] - unlist(Ybar[, i])),
               col = alphabrown)
    abline(h = 0, lty = 2, lwd = 2)
    abline(reg)
    mtext(paste("slope =", signif(summary(reg)$coefficients[2], digits = 3), 
                "intercept =", signif(summary(reg)$coefficients[1], digits = 3)))
    # d<-density(c(Xbar[t1:t] - unlist(Ybar[t1:t,i]))) lines(d$y+1,d$x)
    
    # forecast minus analysis = update
    reg1 <- lm(Xbar - Xa ~ c(t1:t))
    plot(t1:t, 
         Xbar - Xa, 
         pch = 16, cex = 1, 
         ylim = c(min(Xbar - XaCI[, 2]), max(Xbar - XaCI[, 1])), 
         xlab = "Time", ylab = "Update", 
         main = paste(colnames(X)[i], 
                      "Update = Forecast - Analysis"))
    ciEnvelope(rev(t1:t), 
               rev(Xbar - XaCI[, 1]), 
               rev(Xbar - XaCI[, 2]), 
               col = alphapurple)
    abline(h = 0, lty = 2, lwd = 2)
    abline(reg1)
    mtext(paste("slope =", signif(summary(reg1)$coefficients[2], digits = 3),
                "intercept =", signif(summary(reg1)$coefficients[1], 
                                      digits = 3)))
    # d<-density(c(Xbar[t1:t] - Xa[t1:t])) lines(d$y+1,d$x)
  }
  dev.off()
  
  ###-------------------------------------------------------------------###
  ### process variance plots                                            ###-----
  ###-------------------------------------------------------------------### 
  if (processvar) {
    
    library(corrplot)
    pdf('process.var.plots.pdf')
    
    cor.mat <- cov2cor(aqq[t,,] / bqq[t])
    colnames(cor.mat) <- colnames(X)
    rownames(cor.mat) <- colnames(X)
    par(mfrow = c(1, 1), mai = c(1, 1, 4, 1))
    corrplot(cor.mat, type = "upper", tl.srt = 45,order='FPC')
    
    par(mfrow=c(1,1))   
    plot(as.Date(obs.times[t1:t]), bqq[t1:t],
         pch = 16, cex = 1,
         ylab = "Degrees of Freedom", xlab = "Time")
    
    dev.off()
    
  }
  
  ###-------------------------------------------------------------------###
  ### climate plots                                                     ###
  ###-------------------------------------------------------------------### 
  
  # plot(rowMeans(temp.mat[5:t,]),
  #      Xbar[5:t] -  unlist(Ybar[5:t,i]),
  #      xlim=range(rowMeans(temp.mat[5:t,])),
  #      ylim = range(Xbar[5:t] -  unlist(Ybar[5:t,i])),pch=16,cex=1,
  #      xlab="Average Monthly Temp",
  #      ylab="Error",
  #      main=colnames(Ybar)[i])
  # 
  # plot(rowSums(precip.mat[5:t,]),
  #      Xbar[5:t] - unlist(Ybar[5:t,i]),
  #      xlim=range(rowSums(precip.mat[5:t,])),
  #      ylim = range(Xbar [5:t]- unlist(Ybar[5:t,i])),
  #      pch=16,cex=1,xlab="Total Yearly Precip",
  #      ylab="Error",main=colnames(Ybar)[i])
  # 
  # plot(rowMeans(temp.mat[5:t,]),Xbar[5:t] - Xa[5:t],pch=16,
  #      cex=1,xlab="Average Monthly Temp",
  #      ylab="Update",main=colnames(Ybar)[i])
  # plot(rowSums(precip.mat[5:t,]),Xbar[5:t] - Xa[5:t],pch=16,
  #      cex=1, xlab="Total Yearly Precip",
  #      ylab="Update",main=colnames(Ybar)[i])
  
} # sda.enkf

