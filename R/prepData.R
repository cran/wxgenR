#' Read in data to train simulator
#'
#' Read training data and setup variables to facilitate simulation.
#'
#' @param trainingData A path to .csv file (or matrix, data frame in your environment)
#' with the following variables
#'  is required: year, month, day, prcp (daily precipitation),
#'   temp (daily temperature),
#'    and season (1, 2, ..., N, for N seasons - up to 20 seasons will work).
#'    Any units will work for precipitation and temperature, as long as they are
#'    consistent. Can be station data, basin averages, grid cells, etc.
#'
#' @param sdate Start date of training data (yyyymmdd).
#' If empty, the start date will be the beginning of your time series.
#' @param edate End date of training data (yyyymmdd).
#' If empty, the end date will be the end of your time series.
#'
#' @examples
#'
#' prepData(trainingData = "./MetData.csv",  sdate = 20000101, edate = 20201231)
#'
# @import lubridate
#'
#' @noRd


"prepData" <- function(trainingData, syr, smo, eyr, emo, unitSystem, traceThreshold){

  #starting and ending date of simulation
  if(is.null(syr) == T){
    syr = trainingData$year[1]
  }else{
    trainingData = subset(trainingData, year >= syr)
  }

  if(is.null(smo)){
    smo = trainingData$month[1] #if no start month, use beginning of record (already subset by years)
  }else if(is.character(smo)){
    smo = match_month(smo) #otherwise convert character month to numeric
  }

  sdy = trainingData$day[1]

  # trainingData = subset(trainingData, year >= syr & month >= smo & day >= 1)

  sdate = syr*10^4+smo*10^2+sdy

  if(is.null(eyr) == T){
    eyr = tail(trainingData$year,1)
  }else{
    trainingData = subset(trainingData, year <= eyr)
  }

  if(is.null(emo)){
    emo = tail(trainingData$month, 1) #if no end month, use end of record (already subset by years)
    edy = tail(trainingData$day,1)
  }else if(is.character(emo)){
    emo = match_month(emo) #otherwise convert character month to numeric
    edy = tail(subset(trainingData, month == emo)$day,1)
  }else if(is.numeric(emo)){
    edy = tail(subset(trainingData, month == emo)$day,1)
  }

  # trainingData = subset(trainingData, year <= eyr & month <= emo & day <= 30)

  edate = eyr*10^4 + emo*10^2 + edy

  #convert units to U.S. Customary if necessary
  if(unitSystem == "metric" | unitSystem == "Metric"){
    trainingData$prcp = trainingData$prcp/25.4     #mm to inches
    trainingData$temp = trainingData$temp*1.8 + 32 #deg C to deg F
    if(traceThreshold != 0.005) traceThreshold = traceThreshold/25.4           #mm to inches
  }

  dat.d = trainingData

  dat.d$date1 = dat.d$year*10000 + dat.d$month*100 + dat.d$day

  i1 = which(dat.d$date1 == sdate)
  i2 = which(dat.d$date1 == edate)

  dat.d = dat.d[i1:i2,]
  yr.d = dat.d$year
  mo.d = dat.d$month
  da.d = dat.d$day
  wk.d = week(ymd(paste(yr.d, mo.d, da.d, sep="-")))
  wk.d[wk.d == 53] = 52
  dat.d$week = wk.d
  dat.d$date = ymd(paste(yr.d, mo.d, da.d, sep="-"))
  dat.d$states = dat.d$season
  #
  #default
  return(list(dat.d = dat.d, syr = syr, eyr = eyr, smo = smo, emo = emo,
              sdate = sdate, edate = edate, traceThreshold = traceThreshold))
} #end function


