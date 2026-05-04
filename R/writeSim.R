#' Write simulations to file
#'
#' Write simulation results to .csv files (one .csv file is generated for each trace).
#' Inputs include the weather simulations stored in the list object output from the `wx()` function as well as the `nsim` and `nrealz`
#' variables that were inputs to the `wx()` function.\cr
#' \cr
#' A debug flag allows for more detailed reports (debug = TRUE), but setting 'debug = FALSE' is generally
#' recommended for more concise output. Keeping 'debug = FALSE' will also include a simulation
#' time stamp (year, month, day) beginning in year 1.\cr
#' \cr
#' This function will write the .csv files to your working directory.\cr
#' \cr
#' Leap years may be included in the simulated weather if they are included in your training data,
#' so non-leap years include a row of 'NA' values at the end of the calendar year as a book-keeping
#' measure so that the total number of rows in each trace is the same.
#'
#' @param wxOutput Weather simulations output from `wx()` function.
#' @param nsim Number of simulation years.
#' @param nrealz Number of realizations (ensemble size).
#' @param path Specified path to where simulation output shall be written. Defaults to current working directory (path = NULL).
#' Specified path should be a character string of the folder location ending with '/'.
#' @param debug Option to include additional variables in the .csv file
#' outputs for debugging and advanced analysis. Includes sampling date, etc. Default = FALSE (off).
#'     If debug is off, the weather simulations will have a simulation year time stamp
#'     (beginning in year 1) as well as month and day time stamps.
#'
#' @return No return value, called to write simulation results to file.
#'
#' @examples
#'
#' \donttest{
#' z = wx(trainingData = LowerSantaCruzRiverBasinAZ,
#'  eyr = 1990, nsim = 5, nrealz = 5, aseed = 23,
#'   wwidth = 3, unitSystem = "U.S. Customary",
#'    ekflag = TRUE, awinFlag = TRUE, tempPerturb = TRUE,
#'     pcpOccFlag = FALSE, numbCores = NULL)
#'
#'
#' writeSim(wxOutput = z, nsim = 5, nrealz = 5, path = paste0(tempdir(), "/"), debug = FALSE)
#'}
#'
#' @export
#'
# @import utils
#'
#'
"writeSim" <- function(wxOutput, nsim, nrealz, path = NULL, debug = FALSE){

  #parse variables from wx() output
  dat.d = wxOutput$dat.d
  simyr1 = wxOutput$simyr1
  X = wxOutput$X
  Xseas = wxOutput$Xseas
  Xpdate = wxOutput$Xpdate
  Xpamt = wxOutput$Xpamt
  Xtemp = wxOutput$Xtemp

  #write simulation output
  #
  it1 <- seq(1, length(X[,1]), 366)
  it2 = it1+366-1

  #loop through realization
  for (irealz in 1:nrealz){
  outmat <- vector()
  if(is.null(path) == TRUE){
    fout = paste("Realization_", sprintf("%03d", irealz), ".csv", sep="")
  }else if(is.null(path) == FALSE){
    fout = paste(path, "Realization_", sprintf("%03d", irealz), ".csv", sep="")
  }

  if(debug == TRUE){

    #loop through simulation years
    for (isim in 1:nsim){
      leapflag = FALSE
      ayr = simyr1[isim, irealz]
      if (leap_year(ayr)) leapflag = TRUE
      col1 = rep(isim, 366)        #column 1, simulation year
      d1 = ayr*10^4+01*10^2+01; d2 = ayr*10^4+12*10^2+31
      i1 = which(dat.d$date1 == d1)
      i2 = which(dat.d$date1 == d2)
      col2 = dat.d$date1[i1:i2]    #column 2, simulation date
      if (leapflag == FALSE) col2 = c(col2,NA)
      i1 = it1[isim]
      i2 = it2[isim]
      col3 = Xseas[i1:i2, irealz]  #column 3, simulation season
      col4 = X[i1:i2, irealz]      #column 4, precipitation occurrence
      col5 = Xpdate[i1:i2, irealz] #column 5, precipation resampling date
      col6 = Xpamt[i1:i2, irealz]  #column 6, resampled precipitation amount
      col7 = Xtemp[i1:i2, irealz]  #column 7, simulated temperature
      outmat = rbind(outmat, cbind(col1, col2, col3, col4, col5, col6, col7))
    } #isim

    narows = (which(is.na(outmat[,2]) == TRUE))*-1
    outmat = outmat[narows,]

    colnames(outmat) = c("Year", "Date", "Season", "PrcpOccurence", "PrcpResamplingDate", "ResampledPrcpAmt", "SimulatedTemp")

  } else{

    #loop through simulation years
    for (isim in 1:nsim){
      leapflag = FALSE
      ayr = simyr1[isim, irealz]
      if (leap_year(ayr)) leapflag = TRUE
      col1 = rep(isim, 366)        #column 1, simulation year
      d1 = ayr*10^4+01*10^2+01; d2 = ayr*10^4+12*10^2+31
      i1 = which(dat.d$date1 == d1)
      i2 = which(dat.d$date1 == d2)
      col2 = dat.d$date1[i1:i2]   #column 2, simulation date
      if (leapflag == FALSE) col2 = c(col2,NA)
      i1 = it1[isim]
      i2 = it2[isim]
      col3 = Xseas[i1:i2, irealz]  #column 3, simulation season
      col4 = X[i1:i2, irealz]      #column 4, precipitation occurrence
      col5 = Xpdate[i1:i2, irealz] #column 5, precipation resampling date
      col6 = Xpamt[i1:i2, irealz]  #column 6, resampled precipitation amount
      col7 = Xtemp[i1:i2, irealz]  #column7, simulated temperature

      #create time series of 'simulation day'
      sim.yr = rep(isim, length(col2))
      sim.month = month(ymd(col2))
      sim.day = day(ymd(col2))

      outmat = rbind(outmat, cbind(sim.yr, sim.month, sim.day, col6, col7, col3))
    } #isim

    colnames(outmat) = c("simulation year", "month", "day", "prcp", "temp", "season")

  }

  write.table(outmat, fout, row.names = FALSE, col.names = TRUE, sep=",")

} #irealz

  wxgenR_fun_message("writeSim")

} #end function


