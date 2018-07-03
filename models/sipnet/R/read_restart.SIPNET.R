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
                                var.map=list(c("leaf_carbon_content","litter_carbon_content","TotSoilCarb","SoilMoistFrac","SWE","GWBI","AbvGrndWood","AGB"),
                                             c("LeafC","Litter","TotSoilCarb","SoilMoistFrac","SWE","GWBI","AbvGrndWood",'AGB'), # what we want
                                             c('','','','','',"kg/m^2/s","kg/m^2",""), # unit in
                                             c('','','','','',"Mg/ha/yr","Mg/ha",""), #unit out
                                             c('','','','','',mean,"",mean) #preprocess function. Sending the function explicitly                                
                                             )
                                ) {
  
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
  last <- length(ens[[1]])
  forecast <- list()

  #if (length(var.names[!var.names%in%names(ens))>0) PEcAn.logger::logger.warn('There are varibales that they are not maped and so they will not be in the ensemble output.')
    
  #The first vector in the arg is the model's title names for outputs and second is what we wanted it to be
  # the third and fourth are the units in and out
  purrr::pwalk(var.map,
              function(modelT,ForcastT,unit.in,unit.out,FUN){
      if(modelT%in%names(ens)){
        # doing the proprocess - if we need to take mean or etc
        preproc.val<-ifelse(typeof(FUN)!='character',FUN(unlist(ens[[modelT]])),ens[[modelT]][last])
        #do the unit conversion
        value<-ifelse(unit.in!='',udunits2::ud.convert(preproc.val,unit.in,unit.out),preproc.val)
        # store it back to forecast
        forecast<<-c(forecast,setNames(value,ForcastT))
      }          
      
  })

  # Finding some params-------------
  if ("AbvGrndWood" %in% var.names) {
    # calculate fractions, store in params, will use in write_restart
    wood_total_C    <- ens$AbvGrndWood[last] + ens$fine_root_carbon_content[last] + ens$coarse_root_carbon_content[last]
    abvGrndWoodFrac <- ens$AbvGrndWood[last]  / wood_total_C
    coarseRootFrac  <- ens$coarse_root_carbon_content[last] / wood_total_C
    fineRootFrac    <- ens$fine_root_carbon_content[last]   / wood_total_C
    params$restart <- c(abvGrndWoodFrac, coarseRootFrac, fineRootFrac)
    names(params$restart) <- c("abvGrndWoodFrac", "coarseRootFrac", "fineRootFrac")
  }
  print(runid)
  
  X_tmp <- list(X = unlist(forecast), params = params)
                
  return(X_tmp)
} # read_restart.SIPNET
