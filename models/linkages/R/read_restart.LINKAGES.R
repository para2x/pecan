#-------------------------------------------------------------------------------
# Copyright (c) 2012 University of Illinois, NCSA.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the 
# University of Illinois/NCSA Open Source License
# which accompanies this distribution, and is available at
# http://opensource.ncsa.illinois.edu/license.html
#-------------------------------------------------------------------------------

##' @title read_restart.LINKAGES
##' @name  read_restart.LINKAGES
##' @author Ann Raiho \email{araiho@@nd.edu} & Hamze Dokoohaki \email{hamzed@bu.edu}
##' 
##' @param outdir      output directory
##' @param runid       run ID
##' @param stop.time   year that is being read
##' @param multi.settings    PEcAn settings object
##' @param var.names   var.names to be extracted
##' 
##' @description Read Restart for LINKAGES
##' 
##' @return X.vec      vector of forecasts
##' @export
##' 
read_restart.LINKAGES <- function(outdir, runid, stop.time, settings, var.names = NULL, params = NULL,
                                  var.map=list(c('AGB.pft','TotSoilCarb'),
                                               c('AGB.pft','TotSoilCarb'), # what we want
                                               c('',''), # unit in
                                               c('',''), #unit out
                                               c('','') #preprocess function. Sending the function explicitly                                
                                  ),
                                  timez="UTC",
                                  When=NULL) {
  
  # Read ensemble output
  ens <- read.output(runid = runid, 
                     outdir = file.path(outdir, runid), 
                     start.year = lubridate::year(stop.time), 
                     end.year = lubridate::year(stop.time), 
                     variables = var.names)  # change to just 'AGB' for plot level biomass


  # Add PFT name to variable if applicable
  pft.names <- numeric(length(settings$pfts))
  for (i in seq_along(settings$pfts)) {
    pft.names[i] <- settings$pfts[i]$pft$name
  }
  #ens.pft.names <- grep("pft", names(ens))
  #names(ens[[grep("pft", names(ens))]]) <- pft.names

  forecast <- list()
  
  #The first vector in the arg is the model's title names for outputs and second is what we wanted it to be
  # the third and fourth are the units in and out
  purrr::pwalk(var.map,
               function(modelT,ForcastT,unit.in,unit.out,FUN){
                 if(modelT%in%names(ens)){
                   # doing the proprocess - if we need to take mean or etc
                   if(typeof(FUN)!='character'){
                     preproc.val<-FUN(unlist(ens[[modelT]]))
                   } else{
                     preproc.val<-unlist(ens[[modelT]])
                   }
                
                   #do the unit conversion
                   if(unit.in!=''){
                      value<-udunits2::ud.convert(preproc.val,unit.in,unit.out)
                   }else{
                     value<-preproc.val
                   }
                   #see how many ouputs are genrated. as many as the pfts ?
                   if(length(value)==length(pft.names)){nns<-paste0(ForcastT,pft.names)}else{nns<-ForcastT}
                   # store it back to forecast
                   forecast<<-c(forecast,setNames(value,nns))
                 }          
                 
               })
 
  # Put forecast into vector
  print(runid)
  return(list(unlist(forecast),NULL))
}
