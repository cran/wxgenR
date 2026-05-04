#' Runs weather generator
#'
#' Runs the weather generator based on user inputs.\cr
#' \cr
#' Your input/training data MUST have the following variables,
#'  in this order: year, month, day, prcp, temp, season. These variables are case sensitive
#'  and must be spelled as specified here.\cr
#'  \cr
#'  Your training data must be a temporally complete time series (i.e., the time series includes all expected timestamps, even if data is missing).
#'  The training data is expected by default to start at the beginning of the calendar year (January 1) but custom year definitions (e.g., water years)
#'  can be set using the 'smo' and 'emo' arguments to define start and end months, respectively.\cr
#'  \cr
#'  Use starting- and ending- years to subset your input data if desired;
#'  otherwise starting and ending dates will default to the beginning and end of your dataset.\cr
#'  \cr
#' Using 'ekflag = T' will generate simulations outside of the historical envelope
#' via an Epanechnikov kernel. For more details on the Epanechnikov kernel and its use
#' in a weather generator, see Rajagopalan et al. (1996).\cr
#'  \cr
#'  Leap years may be included in the simulated weather if they are included in your training data,
#'  so non-leap years include a row of 'NA' values at the end of the calendar year as a book-keeping
#'  measure so that the total number of rows in each trace is the same.\cr
#'  \cr
#'  The weather generator can handle missing precipitation and temperature data if it is
#'  marked as `NA` in your training data. It will set `NA` precipitation values to 0 and pass along `NA` temperature values
#'  if that date is sampled for the simulations. Consider replacing any missing data with monthly or
#'  daily averages to avoid `NA` values in your simulated weather.
#'
#'
#' @param trainingData Either a matrix, dataframe, or path to a .csv file with the following variables
#'  is required: year, month, day, prcp (daily precipitation),
#'   temp (daily temperature),
#'    and season (1, 2, ..., N, for N seasons - up to 26 seasons will work but seasons need to be defined in a meaningful way).
#'    Units must be either U.S. Customary (inches, degrees F) or metric (mm, degrees C) and must be specified with
#'    the `unitSystem` input variable. Input data can be station-based, basin averages, grid cells, etc.
#'        Input data MUST have these variables: year, month, day, prcp, temp, season.
#' @param syr Optional: subset training data to specific start year (defaults to beginning of training data). Subset will begin on the first day available in `syr`.
#' @param smo Training data start month (you can also use to subset your training data).
# @param sdy Training data start day (you can also use to subset your training data).
#' @param eyr Optional: subset training data to specific end year (defaults to end of training data). Subset will end on the last day available in `eyr`.
#' @param emo Training data end month (you can also use to subset your training data).
# @param edy Training data end day (you can also use to subset your training data).
#' @param nsim Number of simulation years.
#' @param nrealz Number of realizations or traces (i.e., ensemble size).
#' @param aseed Specify a seed for reproducibility.
#' @param wwidth Set the sampling window for each day of year, a lower value for `wwidth` will sample fewer surrounding days (lower variability) and a higher value will sample more days (higher variability). Typical setting of `wwidth` is between 2 and 15, resulting in a daily sampling window of 5 days and 31 days, respectively.
#' Can either be a single number for a uniform window width through the year, or a vector of window widths specific to each season in the training data. In the case of variable window widths, the number of window widths should be equal to the number of seasons.
#' @param unitSystem Specify the unit system of your training data. Input a string that is either "U.S. Customary" or "Metric". U.S. Customary corresponds to inches and degrees Fahrenheit, while Metric corresponds to millimeter and degrees Celsius.
#' If Metric is specified, units will automatically be converted to U.S. Customary for weather simulation, then re-converted to Metric for results output.
#' @param ekflag Simulate outside historical envelope using an Epanechnikov kernel? (T/F)
#' @param awinFlag Set to T or TRUE if you would like to see the results of the adaptive window width.
#'  If only one or zero precipitation values (>0.01 inches) are found within the initial window width you set from a day where precipitation occurred,
#'  it will be iteratively increased until two or more precipitation values are found. By default, the results are not shown.
#' @param tempPerturb Set to T or TRUE if you would like to add random noise to the
#' temperature simulations based on a normal distribution fit on the training data.
#' @param pcpOccFlag Set to TRUE to use precipitation occurrence as a variable in the temperature simulation model or set to FALSE to omit precipitation occurrence as a variable.
#' Simulated daily temperature uses concurrent daily precipitation occurrence as a variable if enabled. By default, this is turned off.
#' @param traceThreshold Threshold for determining whether precipitation depth is considered a trace amount or not. Precipitation depths below this value will be considered trace amounts and will not be used for simulation. A default value of 0.005-inches is used based on National Weather Service guidance. If using a custom trace depth, ensure that it is in the same unit system as your training data and specified by the `unitSystem` flag.
#' @param numbCores Enable parallel computing for precipitation simulation, set number of cores to enable (must be a positive integer greater than or equal to 2). Turned off by default; if set to 0 or 1 it will run as single thread. Use function 'detectCores()' from 'parallel' package to show the number of available cores on your machine.
#' @param returnTempModel Optional flag to return the fitted linear model for daily temperature simulation along with simulation results. Enable by setting TRUE (FALSE by default).
#'
#' @return Returns a list containing both inputs to the weather generator as well as outputs.
#' \itemize{
#'   \item dat.d - User inputs to weather generator, saved for future use.
#'   \item simyr1 - The years sampled for each trace.
#'   \item X - The simulated daily dry/wet sequences for each trace (0 = dry, 1 = wet).
#'   \item Xseas - The simulated season by day for each trace.
#'   \item Xpdate - If precipitation was simulated to occur on a given day, this is the date from which historical precipitation is sampled.
#'   \item Xpamt - The simulated daily precipitation depth.
#'   \item Xtemp - The simulated daily mean temperature.
#' }
#'
#' @examples
#'
#' \donttest{
#'
#' data(LowerSantaCruzRiverBasinAZ)
#'
#' head(LowerSantaCruzRiverBasinAZ)
#'
#' #No input for `syr` because we want the training period to begin at the beginning of the data
#' #record (1970), but set `eyr` = 1990 because we want to subset training period to end in 1990.
#'
#' wx(trainingData = LowerSantaCruzRiverBasinAZ,
#'  eyr = 1990, nsim = 3, nrealz = 3, aseed = 23,
#'   wwidth = 3, unitSystem = "U.S. Customary",
#'    ekflag = TRUE, awinFlag = TRUE, tempPerturb = TRUE,
#'     pcpOccFlag = FALSE, numbCores = NULL)
#'
#'}
#'
#' @export
#'
# @importFrom plyr ddply
# @importFrom dplyr group_by summarise left_join glimpse mutate relocate if_else filter
# @import lubridate
# @import msm
# @import sm
# @import doParallel
# @import parallel
# @import foreach
# @import utils
# @import magrittr
#'
#'
"wx" <- function(trainingData, syr = NULL, eyr = NULL, smo = NULL, emo = NULL,
                 nsim, nrealz, aseed, wwidth, unitSystem,
                 ekflag, awinFlag = FALSE, tempPerturb, pcpOccFlag = FALSE,
                 traceThreshold = 0.005, numbCores = NULL, returnTempModel = F
                ){
  #weather generator
  #
  # require("lubridate")
  if(typeof(trainingData) == "character"){
    trainingData = read.table(trainingData, header=T, sep = ",")
  }

  #
  #read data and setup variables to facilitate simulation
  preppedData <- prepData(trainingData, syr, smo, eyr, emo, unitSystem, traceThreshold)
  dat.d = preppedData$dat.d
  syr = preppedData$syr
  eyr = preppedData$eyr
  smo = preppedData$smo
  emo = preppedData$emo
  sdate = preppedData$sdate
  edate = preppedData$edate
  traceThreshold = preppedData$traceThreshold

  #
  #calculate seasonal precipitation transition prob matrix for each year
  ptpm <- getPtpm(dat.d, traceThreshold)
  tpm.y2 <- ptpm$tpm.y2
  tpm.y <- ptpm$tpm.y
  #
  #calculate parameters for temperature simulation
  z.t <- getTpars(dat.d, pcpOccFlag, traceThreshold, smo, emo, returnTempModel)
  dat.d=z.t$dat.d     #updated with tavgm, sine and cosine terms
  coeftmp=z.t$coeftmp
  tmp.sd=z.t$tmp.sd
  #
  #simulate precipitation occurrence and temperature
  message("...Simulate precipitation occurrence and temperature...")
  z <- simTPocc(aseed,dat.d,smo,emo,nsim,nrealz,coeftmp,tmp.sd,tpm.y2,tpm.y,tempPerturb,pcpOccFlag)
  simyr1=z$simyr1
  X=z$X
  Xjday=z$Xjday
  Xseas=z$Xseas
  Xtemp=z$Xtemp
  #
  #simulate precipitation amount
  message("...Simulate precipitation amount...")
  wwidth = as.numeric(wwidth)
  z <- simPamt(dat.d,syr,eyr,smo,emo,sdate,edate,wwidth,nsim,nrealz,Xjday,Xseas,ekflag,awinFlag,numbCores,traceThreshold,aseed)
  Xpamt <- z$Xpamt
  Xpdate <- z$Xpdate
  if (ekflag) bSJ <- z$bSJ
  #

  #re-convert units back to metric if necessary
  if(unitSystem == "metric" | unitSystem == "Metric"){
    #simulations
    Xpamt = Xpamt*25.4        #inches to mm
    Xtemp = (Xtemp-32)*(5/9)  #deg F to deg C
    #observed
    dat.d$prcp = dat.d$prcp*25.4 #inches to mm
    dat.d$temp = (dat.d$temp-32)*(5/9) #deg F to deg C
    dat.d$tavgm = (dat.d$tavgm-32)*(5/9) #deg F to deg C
  }

  #default
  olist=list("dat.d"=dat.d,"simyr1"=simyr1,"X"=X,"Xseas"=Xseas,
             "Xpdate"=Xpdate,"Xpamt"=Xpamt,"Xtemp"=Xtemp
            )

  # if (ekflag) olist=c(olist)

  if(returnTempModel == F){
    olist=list("dat.d"=dat.d,"simyr1"=simyr1,"X"=X,"Xseas"=Xseas,
               "Xpdate"=Xpdate,"Xpamt"=Xpamt,"Xtemp"=Xtemp
    )
  }else{
    olist=list("dat.d"=dat.d,"simyr1"=simyr1,"X"=X,"Xseas"=Xseas,
               "Xpdate"=Xpdate,"Xpamt"=Xpamt,"Xtemp"=Xtemp,"tmp.mod" = z.t$tmp.mod
    )
  }

  wxgenR_fun_message("wx")

  return(olist)
} #end function

