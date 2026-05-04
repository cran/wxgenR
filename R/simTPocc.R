#' Simulate precipitation occurrence and temperature magnitude
#'
#' Simulate daily precipitation occurrence (wet or dry)
#' and daily temperature for the desired length of time.
#' One of the covariates in the temperature model is
#' precipitation occurrence.
#'
#' @param aseed Random number seed.
#' @param dat.d Training data processed from prepData wrapper function.
#' @param nsim Number of simulation years.
#' @param nrealz Number of realizations.
#' @param coeftmp Coefficients from linear model for predicting daily temperature.
#' @param tpm.sd Standard deviation of temperature model residuals by month.
#' @param tpm.y2 Transition probability matrix for each year in training data.
#' @param tpm.y Aggregate transition probability matrix for all years in training data.
#'
# @import lubridate
# @rawNamespace import(stats, except = filter)
#'
#' @noRd
#'

"simTPocc" <- function(aseed, dat.d, smo, emo, nsim, nrealz, coeftmp, tmp.sd, tpm.y2, tpm.y, tempPerturb, pcpOccFlag){
  #simulate precipitation occurrence and temperature magnitude
  #
  #initialize arrays
  X <- matrix(NA, nrow = nsim*366, ncol = nrealz)     #saves prcp occurrence
  simyr1 <- matrix(NA, nrow = nsim, ncol = nrealz)    #saves ALL selected simulation years
  Xjday <- matrix(NA, nrow = nsim*366, ncol = nrealz) #saves julian day when prcp occurs
  Xtemp <- matrix(NA, nrow = nsim*366, ncol = nrealz) #saves temperature
  Xseas <- matrix(NA, nrow = nsim*366, ncol = nrealz) #saves season
  Xweek <- matrix(NA, nrow = nsim*366, ncol = nrealz) #saves week
  yr.d = dat.d$year
  uyr = unique(yr.d)
  nyr = length(uyr)
  set.seed(aseed) #set seed
  #

  #get month for a given julian day
  it1 = get_reference_366_start(dat.d, smo, "precipitation occurrence and temperature simulation")
  it2 = it1+366-1
  jdaymth <- dat.d$month[it1:it2]
  zz <- rep(jdaymth, nsim) #simulation month
  yy <- rep(1:nsim, each=366) #simulation year

  #
  #realization loop
  irealz = 1
  for (irealz in 1:nrealz){
    message(irealz)
    prcpocc <- temp <- matrix(NA, nrow=366, ncol=nsim) #366-day
    pseas <- pweek <- matrix(NA, nrow=366, ncol=nsim)  #366-day
    simyr <- rep(NA, nsim)

    #loop through simulation years
    isim=1
    for (isim in 1:nsim){

      if(smo == 1){ #leap year or not for calendar year

        iyr = sample(1:nyr, 1)          #randomly select a year index
        simyr[isim] = uyr[iyr]
        simyr1[isim,irealz] = simyr[isim]

        leapflag = 0
        if (leap_year(uyr[iyr])) leapflag = 1
        nt = 365 + leapflag
        startdate <- ymd(paste(uyr[iyr], smo, 1, sep="-")) #1st date of year uyr[iyr]

        # startdate <- ymd(paste(uyr[iyr], 1, 1, sep="-")) #jan 1 of year uyr[iyr]
        it1 = which(dat.d$date == startdate) #starting index of data
        it2 = it1+nt-1
        # dframe = subset(dat.d, )
        dframe <- dat.d[it1:it2,] #data subset for simulation

      }else if(smo > 1){ #leap year or not for water year

        iyr = sample(1:(nyr-1), 1)          #randomly select a year index
        smpl_yr = uyr[-1][iyr]
        simyr[isim] = smpl_yr          #remove year 1 since not a water year included in data
        simyr1[isim,irealz] = simyr[isim]

        leapflag = 0
        if (leap_year(smpl_yr)) leapflag = 1
        nt = 365 + leapflag
        startdate <- ymd(paste(smpl_yr-1, smo, 1, sep="-")) #1st date of year uyr[iyr]

        it1 = which(dat.d$date == startdate) #starting index of data
        it2 = it1+nt-1
        # dframe = subset(dat.d, )
        dframe <- dat.d[it1:it2,] #data subset for simulation

      }
      prcpocc[1,isim] = dframe$oc[1]

      if(is.na(dframe$temp[1]) == F){
        temp[1,isim] = dframe$temp[1]      #temperature for day=1 of simulation
      }else if(is.na(dframe$temp[1]) == T & is.na(dframe$tavgm[1]) == F){ #if day 1 temp is missing, use monthly average
        temp[1,isim] = dframe$tavgm[1]      #temperature for day=1 of simulation
      }else{ #if monthly average is NA for that year, randomly sample other years until a monthly average is found
        while(is.na(temp[1,isim])){
          if(smo == 1){
            iyr.t = sample(1:nyr, 1)          #randomly select a year index
            simyr.t = uyr[iyr.t]
            startdate.t <- ymd(paste(simyr.t, smo, 1, sep="-"))
          }else{
            iyr.t = sample(1:(nyr-1), 1)
            simyr.t = uyr[-1][iyr.t]
            startdate.t <- ymd(paste(simyr.t-1, smo, 1, sep="-"))
          }
          simyr1[isim,irealz] = simyr[isim]

          it1.t = which(dat.d$date == startdate.t) #starting index of data
          it2.t = it1.t+nt-1
          dframe.t <- dat.d[it1.t:it2.t,] #data subset for simulation

          temp[1,isim] = dframe.t$tavgm[1]

        }
      }

      if(is.na(dframe$oc[1])){
        pstate=1 #if occurence is NA, set to dry
      }else if(dframe$oc[1]==0){
        pstate=1      #dry, corresponds to row number 1
      }else if(dframe$oc[1]==1){
        pstate=2      #wet, corresponds to row number 2
      }

      aseas <- dframe$season[1]          #selected season for sim day 1
      aweek <- dframe$week[1]            #week number for sim day 1
      ptvec <- tpm.y2[pstate,,aseas,iyr] #prcp, prob transition vector

      #Use 30-year average probabilities if transition is not found within the season that year
      if(is.na(ptvec[1]) == TRUE | is.na(ptvec[2]) == TRUE) ptvec=tpm.y[pstate,,aseas]
      pseas[1,isim] = aseas
      pweek[1,isim] = aweek

      #loop through days in simulated year
      it = 2
      for (it in 2:nt){
        #precipitation occurrence
        u = runif(1)
        pstate = selectState(u,cumsum(ptvec))
        prcpocc[it,isim] = pstate - 1     #0 - dry; 1 - wet
        aseas <- dframe$season[it]
        aweek <- dframe$week[it]
        pseas[it,isim] = aseas
        pweek[it,isim] = aweek
        ptvec <- tpm.y2[pstate,,aseas,iyr]
        if(is.na(ptvec[1]) == TRUE | is.na(ptvec[2]) == TRUE) ptvec = tpm.y[pstate,,aseas]

        #temperature simulation using linear model coefficients
        if(tempPerturb == T & pcpOccFlag == T){
          temp[it,isim] = sum(coeftmp*c(1,temp[(it-1),isim], dframe$ct[it],
                                       dframe$st[it], prcpocc[it,isim],
                                       dframe$tavgm[it]
                                       )
                              ) + rnorm(n=1,mean=0,sd=tmp.sd[jdaymth[it]])
        } else if(tempPerturb == F & pcpOccFlag == T){
          temp[it,isim] = sum(coeftmp*c(1,temp[(it-1),isim], dframe$ct[it],
                                        dframe$st[it], prcpocc[it,isim],
                                        dframe$tavgm[it]
                                        )
                              )
        } else if(tempPerturb == T & pcpOccFlag == F){
          temp[it,isim] = sum(coeftmp*c(1,temp[(it-1),isim], dframe$ct[it],
                                        dframe$st[it],
                                        dframe$tavgm[it]
                                        )
                              ) + rnorm(n=1,mean=0,sd=tmp.sd[jdaymth[it]])
        } else if(tempPerturb == F & pcpOccFlag == F){
          temp[it,isim] = sum(coeftmp*c(1,temp[(it-1),isim], dframe$ct[it],
                                        dframe$st[it],
                                        dframe$tavgm[it])
                              )
        }

      } #it
     } #isim
     #make arrays
     #precipitation occurrence
     x1 <- prcpocc                    #dimension of prcpocc is 366 x nsim
     x2 <- vector()
     x2 <- append(x2,x1[,1:nsim])
     X[,irealz] = x2                  #precipitation occurrence
     #temperature
     x3 <- temp                       #dimension of temp is 366 x nsim
     x4 <- vector()
     x4 <- append(x4, x3[,1:nsim])
     Xtemp[,irealz] = x4              #temperature
     #seasonality
     x5 <- pseas
     x6 <- vector()
     x6 <- append(x6,x5[,1:nsim])
     Xseas[,irealz] = x6              #season
     #week
     x7 <- pweek
     x8 <- vector()
     x8 <- append(x8, x7[,1:nsim])
     Xweek[,irealz] = x8              #week
  } #irealz
  #
  #get julian day of precipitation occurrence
  it1 = seq(1, dim(X)[1], 366)
  it2 = it1+366-1
  for (irealz in 1:nrealz){
    for(isim in 1:nsim){
      i1 = it1[isim]
      i2 = it2[isim]
      xx <- X[i1:i2, irealz]
      pday <- which(xx == 1)
      idxlist = i1+pday-1
      Xjday[idxlist, irealz] = pday
    } #isim
  } #irealz
  #
  X1 <- cbind(yy,zz,X) #simulation year, simulation month, prcp occurrence
  #
  #default
  olist=list("simyr1" = simyr1, "X" = X, "X1" = X1, "Xjday" = Xjday,
             "Xseas" = Xseas, "Xweek" = Xweek, "Xtemp" = Xtemp
            )
  return(olist)
} #end function

