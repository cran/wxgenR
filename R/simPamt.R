#' Precipitation simulator
#'
#' Simulate daily precipitation depth.
#'
#' @param dat.d Training data processed from prepData function.
#' @param syr Start year.
#' @param eyr End year.
#' @param wwidth Window set for finding surrounding days: +/- wwidth.
#' @param nsim Number of simulation years.
#' @param nrealz Number of realizations.
#' @param Xjday Simulation Julian day when prcp occurs, from simTPocc function
#' @param Xseas Simulation season index from simTPocc function
#' @param ekflag Simulate outside historical envelope?
#' @param numbCores Enable parallel computing for precipitation simulation, set number of cores to enable (must be a positive integer greater than or equal to 2). Turned off by default; if set to 0 or 1 it will run as single thread. Use function 'detectCores()' from 'parallel' package to show the number of available cores on your machine.
#' @param aseed Specify a seed for reproducibility.
#'
# @import lubridate
# @import parallel
# @import doParallel
# @import foreach
# @import sm
# @rawNamespace import(stats, except = filter)
#'
#' @noRd
#'

"simPamt" <- function(dat.d,syr,eyr,smo,emo,sdate,edate,wwidth,nsim,nrealz,Xjday,Xseas,ekflag,awinFlag,numbCores,traceThreshold,aseed){

  if (is.null(numbCores) || !is.numeric(numbCores) || numbCores < 2) {
    numbCores = 1
  }

  # require("sm")
  #simulate precipitation amounts and select dates
  #
  #get month for a given julian day
  it1 = get_reference_366_start(dat.d, smo, "precipitation amount simulation")
  it2 = it1+366-1
  jdaymth <- dat.d$month[it1:it2] #calculations are based on a 366-day year
  #
  #get dates in window for each julian day 1-366
  lw = length(wwidth)
  if(lw == 1){
    Xdates = getDatesInWindow(syr,eyr,smo,emo,sdate,edate,wwidth,leapflag=T)
  } else if(lw > 1){
    Xdates.vw = list(length = lw)
    for(i in 1:lw){
      # print(wwidth[i])
      Xdates.vw[[i]] = getDatesInWindow(syr,eyr,smo,emo,sdate,edate,wwidth[i],leapflag=T)
    }
    names(Xdates.vw) = wwidth
  }

  #Run routine for Epanchnikov kernal if enabled
  if (ekflag){

    #get the bandwidth for each julian day
    bSJ <- if(lw > 1) matrix(NA, nrow = 366, ncol = lw) else vector()
    diwprcp <- vector()              #days in window with precipitation

    jday = 1
    for (jday in 1:366){

      # print(jday)

      seas.bandwidth <- if(lw > 1) 1:lw else 1

      for(seas.j in seas.bandwidth){
        #if variable window width, calculate the bandwidth for each season-specific window
        if(lw > 1){
          Xdates = Xdates.vw[[seas.j]]
          wwidth.adapt = wwidth[seas.j]
        } else{
          wwidth.adapt = wwidth
        }

        diw <- na.omit(Xdates[jday,])  #dates in window for a given julian day

        idxlist <- vector()

        iday = 1
        for (iday in 1:length(diw)){
          idxlist[iday] = which(dat.d$date1==diw[iday])
        } #iday

        baprcp <- dat.d$prcp[idxlist]       #basin average precipitation
                                            #also includes 0 prcp amount
                                            #for days within the window

        pamt <- baprcp[which(baprcp>=traceThreshold)] #precipitation amount vector

        #### Begin adaptive window width process if less than 2 prcp amounts exist
        while(length(pamt) < 2 | all(diff(pamt) == 0)){

          wwidth.adapt = wwidth.adapt + 1

          #get dates in window for each julian day 1-366
          Xdates.adapt=getDatesInWindow(syr,eyr,smo,emo,sdate,edate,wwidth.adapt,leapflag=T)
          diw <- na.omit(Xdates.adapt[jday,])  #dates in window for a given julian day
          idxlist <- vector()

          for (iday in 1:length(diw)){
            idxlist[iday]=which(dat.d$date1==diw[iday])
          } #iday
          baprcp <- dat.d$prcp[idxlist]

          pamt <- baprcp[which(baprcp>=traceThreshold)] #precipitation amount vector

        }

        diwprcp[jday] = length(pamt)
        logpamt <- log(pamt)                #log-transformed precipitation amount vector
        if(lw > 1){
          bSJ[jday, seas.j] = sm::hsj(logpamt) #Sheather-Jones plug-in bandwidth
        }else{
          bSJ[jday] = sm::hsj(logpamt)
        }
      }

    } #jday
  } #ekflag
  #
  #simulate precipitation amount and select prcp date
  # Xpdate <- Xpamt <- matrix(NA, nrow=nsim*366, ncol=nrealz)

  if(numbCores >= 2){

    Xpdate <- Xpamt <- matrix(NA, nrow=nsim*366, ncol=1)

    # library(foreach)
    # library(doParallel)

    if(numbCores > detectCores()){
      message(paste0("numbCores is set to more than the ", detectCores(), " available cores on your machine. Setting numbCores to ", detectCores()-1, " cores (leaving one core free)."))
      message(paste0("If you wish, interupt code to set a new value for numbCores that is between 2 and ", detectCores(), ", otherwise the simulation will continue."))
    }

    cl = makePSOCKcluster(numbCores)
    registerDoParallel(cl)

    # startTime <- Sys.time() #benchmark run time

    irealz = 1
    # result <- foreach(irealz=1:nrealz, .export=c('repan', 'getDatesInWindow')) %dopar% {
    result <- foreach(irealz=1:nrealz, .export=c('repan', 'getDatesInWindow', 'get_sim_season'),
                      .options.RNG = if (exists("aseed")) aseed else NULL) %dorng% {

    # for (irelz in 1:nrealz){

      message(paste0("-- Starting trace number ", irealz, " --"))
      xp = Xjday[,irealz]
      nxp = length(xp)

      ixp = 1
      # foreach(ixp=1:nxp, .combine = c) %do% {
      for (ixp in 1:nxp){

        jd = xp[ixp]

        if (is.na(jd)){
          Xpamt[ixp] = 0.0
        }else{

          #if variable window width, find the season in which current jday exists
          if(lw > 1){
            seas.j = get_sim_season(Xseas, ixp, irealz)
            if(is.na(seas.j) | seas.j < 1 | seas.j > lw) stop("Unable to match simulated day to a valid season-specific window width.", call. = FALSE)
            Xdates = Xdates.vw[[seas.j]]
            #set adaptive window width to window width for current season
            wwidth.adapt = wwidth[seas.j]
          } else{
            wwidth.adapt = wwidth
          }

          diw <- na.omit(Xdates[jd,])  #dates in window for a given julian day
          idxlist <- vector()

          iday = 1
          for (iday in 1:length(diw)){
            idxlist[iday] = which(dat.d$date1 == diw[iday])
          } #iday

          baprcp <- dat.d$prcp[idxlist]       #e.g., basin average precipitation
          #also includes 0 prcp amount
          #days within the window

          np = length(which(baprcp>=traceThreshold))      #number of prcp days in window
          pdate <- diw[which(baprcp>=traceThreshold)]   #dates in window where prcp occurred
          pamt <- baprcp[which(baprcp>=traceThreshold)] #precipitation amount vector

          #Adaptive window width routine
          while(np < 2){

            wwidth.adapt = wwidth.adapt + 1

            #get dates in window for each julian day 1-366
            Xdates.adapt = getDatesInWindow(syr,eyr,smo,emo,sdate,edate,wwidth.adapt,leapflag=T)

            diw <- na.omit(Xdates.adapt[jd,])  #dates in window for a given julian day
            idxlist <- vector()

            iday = 1
            for (iday in 1:length(diw)){
              idxlist[iday]=which(dat.d$date1==diw[iday])
            } #iday

            baprcp <- dat.d$prcp[idxlist]       #e.g., basin average precipitation
            #also includes 0 prcp amount
            #days within the window

            np = length(which(baprcp>=traceThreshold))      #number of prcp days in window
            pdate <- diw[which(baprcp>=traceThreshold)]   #dates in window where prcp occurred
            pamt <- baprcp[which(baprcp>=traceThreshold)] #precipitation amount vector

            if(awinFlag == T){
              message(paste0("\n Window width too small on Julian day ", jd,", increasing window to ", wwidth.adapt*2+1, " days"))
            }
          }

          logpamt <- log(pamt)                #log-transformed prcp amount vector
          aindex = sample(1:np, 1)               #randomly pick a prcp day
          Xpdate[ixp] = pdate[aindex]
          ybar = logpamt[aindex]
          Xpamt[ixp] = exp(ybar)
          if (ekflag){
            rek = repan(1)                      #simulate a random number from the EKD
            bandwidth = if(lw > 1) bSJ[jd, seas.j] else bSJ[jd]
            Xpamt[ixp] = exp(ybar+rek*bandwidth)
          } #ekflag
        }
      } #ixp

      message("\n")

      list(Xpamt, Xpdate)

    } #irealz

    stopCluster(cl)

    # endTime = Sys.time()

    # timeP = difftime(endTime, startTime, units='mins')

    Xpamt = as.data.frame(do.call(cbind,lapply(result,function(x){x[[1]]})))
    Xpdate = as.data.frame(do.call(cbind,lapply(result,function(x){x[[2]]})))

    #
    #default
    olist=list("Xpamt"=Xpamt,"Xpdate"=Xpdate)
    if (ekflag) olist=list("Xpamt"=Xpamt,"Xpdate"=Xpdate,"bSJ"=bSJ)
    return(olist)

  }else{ #non-parallel loop

      Xpdate <- Xpamt <- matrix(NA, nrow=nsim*366, ncol=nrealz)

      # startTime <- Sys.time() #benchmark run time

      irealz = 1
      for (irealz in 1:nrealz){

      message(paste0("-- Starting trace number ", irealz, " --"))
      xp = Xjday[,irealz]
      nxp = length(xp)

      ixp = 1
      for (ixp in 1:nxp){

        jd = xp[ixp]

        if (is.na(jd)){
          Xpamt[ixp,irealz] = 0.0
        }else{

          #if variable window width, find the season in which current jday exists
          if(lw > 1){
            seas.j = get_sim_season(Xseas, ixp, irealz)
            if(is.na(seas.j) | seas.j < 1 | seas.j > lw) stop("Unable to match simulated day to a valid season-specific window width.", call. = FALSE)
            Xdates = Xdates.vw[[seas.j]]
            #set adaptive window width to window width for current season
            wwidth.adapt = wwidth[seas.j]
          } else{
            wwidth.adapt = wwidth
          }

          diw <- na.omit(Xdates[jd,])  #dates in window for a given julian day
          idxlist <- vector()

          iday = 1
          for (iday in 1:length(diw)){
            idxlist[iday] = which(dat.d$date1 == diw[iday])
          } #iday

          baprcp <- dat.d$prcp[idxlist]       #e.g., basin average precipitation
                                              #also includes 0 prcp amount
                                              #days within the window

          np = length(which(baprcp>=traceThreshold))      #number of prcp days in window
          pdate <- diw[which(baprcp>=traceThreshold)]   #dates in window where prcp occurred
          pamt <- baprcp[which(baprcp>=traceThreshold)] #precipitation amount vector

          #adaptive window width routine
          while(np < 2){

            wwidth.adapt = wwidth.adapt + 1

            #get dates in window for each julian day 1-366
            Xdates.adapt = getDatesInWindow(syr,eyr,smo,emo,sdate,edate,wwidth.adapt,leapflag=T)

            diw <- na.omit(Xdates.adapt[jd,])  #dates in window for a given julian day
            idxlist <- vector()

            iday = 1
            for (iday in 1:length(diw)){
              idxlist[iday]=which(dat.d$date1==diw[iday])
            } #iday

            baprcp <- dat.d$prcp[idxlist]       #e.g., basin average precipitation
            #also includes 0 prcp amount
            #days within the window

            np = length(which(baprcp>=traceThreshold))      #number of prcp days in window
            pdate <- diw[which(baprcp>=traceThreshold)]   #dates in window where prcp occurred
            pamt <- baprcp[which(baprcp>=traceThreshold)] #precipitation amount vector

            if(awinFlag == T){
              message(paste0("\n Window width too small on Julian day ", jd,", increasing window to ", wwidth.adapt*2+1, " days"))
            }

          }

          logpamt <- log(pamt)                #log-transformed prcp amount vector
          aindex = sample(1:np, 1)               #randomly pick a prcp day
          Xpdate[ixp,irealz] = pdate[aindex]
          ybar = logpamt[aindex]
          Xpamt[ixp,irealz] = exp(ybar)
          if (ekflag){
            rek = repan(1)                      #simulate a random number from the EKD
            bandwidth = if(lw > 1) bSJ[jd, seas.j] else bSJ[jd]
            Xpamt[ixp,irealz] = exp(ybar+rek*bandwidth)
          } #ekflag
        }
      } #ixp

      message("\n")
    } #irealz

    # endTime = Sys.time()

    # timeNP = difftime(endTime, startTime, units='mins')

    #
    #default
    olist=list("Xpamt"=Xpamt,"Xpdate"=Xpdate)
    if (ekflag) olist=list("Xpamt"=Xpamt,"Xpdate"=Xpdate,"bSJ"=bSJ)
    return(olist)

  } #end non-parallel loop

} #end function

