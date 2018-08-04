#-------------------------------------------------------------------------------
# Copyright (c) 2012 University of Illinois, NCSA.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the 
# University of Illinois/NCSA Open Source License
# which accompanies this distribution, and is available at
# http://opensource.ncsa.illinois.edu/license.html
#-------------------------------------------------------------------------------

##' @title Read restart function for SDA with SIPNET
##' 
##' @author Ann Raiho \email{araiho@@nd.edu}, Hamze Dokoohaki \email{hamzed@bu.edu}
##' 
##' @inheritParams PEcAn.ModelName::read_restart.ModelName
##' 
##' @description Read Restart for SIPNET
##' 
##' @return X.vec      vector of forecasts
##' @export
read_restart.SIPNET <- function(outdir, runid, stop.time, settings, var.names, params,
                                var.map=list(c("leaf_carbon_content","litter_carbon_content","TotSoilCarb","SoilMoistFrac","SWE","GWBI",   "AbvGrndWood","AGB"),
                                             c("LeafC",              "Litter",               "TotSoilCarb","SoilMoistFrac","SWE","GWBI",   "AbvGrndWood",'AGB'), # what we want
                                             c('','','','','',                                                                   "kg/m^2/s","kg/m^2",  ""), # unit in
                                             c('','','','','',                                                                   "Mg/ha/yr","Mg/ha",   ""), #unit out
                                             c('','','','','',                                                                    mean,     "",        "") #preprocess function. Sending the function explicitly                                
                                             ), 
                                timez="UTC",
                                When=NULL,
                                control=list(trace=F),...
                                ) {
#see if there is something else coming
  dots<-list(...)
  if (length(dots)>0) lapply(names(dots),function(name){assign(name,dots[[name]], pos=1 )})
  
  prior.sla <- params[[which(!names(params) %in% c("soil", "soil_SDA", "restart"))[1]]]$SLA
  forecast <- list()
  # additional varnames, because we need these deterministic relationships
  #  var.names <- c(var.names, "fine_root_carbon_content", "coarse_root_carbon_content")
  # Read ensemble output

  if("AbvGrndWood" %in% var.names) var.names<-c(var.names,'fine_root_carbon_content','coarse_root_carbon_content')
  
  ens <- read.output(runid = runid, 
                     outdir = file.path(outdir, runid), 
                     start.year = lubridate::year(stop.time), 
                     end.year = lubridate::year(stop.time),
                     variables = var.names)
  if(all(is.na(unlist(ens)))) PEcAn.logger::logger.error("No output has been generated. Either asking for the wrong year or error in simulation run.")
  forecast <- list()
  ##I'm kaing the timefrmae for the simulations given the year and the length of the output. leap years needs to be check for 30 or 31 in date
  where.index<-c()
  #When<-as.POSIXct(When,tz=timez) # making the timezone right for When
  Timeframe<-seq(as.POSIXct(paste0(year(stop.time),"/01/01"),tz=timez),
                 as.POSIXct(paste0(year(stop.time),"/12/30"),tz=timez),
                 length.out = length(ens[[1]])
                 )
  trunc(Timeframe,"hours")->Timeframe
  if(!is.null(When)){
    #-- find the closest simulation to what is it asked
    lapply(When,function(x){
      which(abs(Timeframe-x) == min(abs(Timeframe - x)))[1]->ww
      ww
    })%>%unlist()%>%na.omit()%>%as.numeric()->where.index
    
  }
  # removing the first and the last - when you ask for a different year that's not in time frime it would send out either the first or the last one
  where.index<-where.index[where.index!=1 & where.index!=length(Timeframe)]
  if(is.null(When) | length(where.index)==0)    where.index <- length(ens[[1]])

  #The first vector in the arg is the model's title names for outputs and second is what we wanted it to be
  # the third and fourth are the units in and out
  purrr::pwalk(var.map,
              function(modelT,ForcastT,unit.in,unit.out,FUN){
      if(modelT%in%names(ens)){
        # doing the proprocess - if we need to take mean or etc
        if(typeof(FUN)!='character'){preproc.val<-FUN(unlist(ens[[modelT]]))}else{preproc.val<-ens[[modelT]][where.index]}
        #do the unit conversion
        if(unit.in!=''){value<-udunits2::ud.convert(preproc.val,unit.in,unit.out)}else{value<-preproc.val}
        # print(str(value))
        # store it back to forecast
        forecast<<-c(forecast,setNames(list(value),ForcastT))
      }          
      
  })
  
  
  # Finding some params-------------
  if ("AbvGrndWood" %in% var.names) {
    # calculate fractions, store in params, will use in write_restart
    wood_total_C    <- ens$AbvGrndWood[where.index] + ens$fine_root_carbon_content[where.index] + ens$coarse_root_carbon_content[where.index]
    abvGrndWoodFrac <- ens$AbvGrndWood[where.index]  / wood_total_C
    coarseRootFrac  <- ens$coarse_root_carbon_content[where.index] / wood_total_C
    fineRootFrac    <- ens$fine_root_carbon_content[where.index]   / wood_total_C
    params$restart <- c(abvGrndWoodFrac, coarseRootFrac, fineRootFrac)
    names(params$restart) <- c("abvGrndWoodFrac", "coarseRootFrac", "fineRootFrac")
  }
  

  #print(runid)
  X_tmp <- list(X =unlist(forecast), params = params, When=Timeframe[where.index])
                
  return(X_tmp)
} # read_restart.SIPNET
