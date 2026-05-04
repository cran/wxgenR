#' N states probability
#'
#' Returns an integer vector corresponding to n states broken by equal
#' probability or equal distance.
#'
#' @noRd
#'
# @rawNamespace import(stats, except = filter)
#'
#'
nStatProb <-   function(x, n, limit.type = 'prob', limits = NULL, tie = 1, altobs = NULL ){
  # returns an integer vector corresponding to n states broken by equal
  # probability or equal distance
  #
  limit <-
    if(limit.type == 'prob')
      quantile(x,seq(0,1,1/n))
  else if(limit.type == 'equal')
    seq(min(x),max(x),by=diff(range(x))/n)
  else if(limit.type == 'manual')
    limits

  if(!is.null(altobs)) limit <- quantile(altobs,seq(0,1,1/n))

  b <- integer(length(x))

  for(i in 1:(n+1)){
    filter <-
      if(tie == 1)
        x >= limit[i] & x <= limit[i+1]
    else
      x > limit[i] & x <= limit[i+1]

    #only need to set the 1's because b is already 0's
    b[filter] <- as.integer(i-1)
  }

  # if(class(x) == 'ts')
  if(inherits(x, 'ts')){
    return(ts(b,start=start(x),end=end(x)))
  }else{
    return(b)
  }
} #end function

#' TPM
#'
#' Checks transition probability matrix.
#'
# @import msm
#'
#' @noRd
#'
#'
transProbMatrix <- function(x,ns=NULL,limits=NULL,tie=0){

  # require(msm)

  if(is.null(ns)){
    ns <- max(x)
    states <- x
    if(length(unique(states)) > 26) stop('Too many states, specify a smaller number.')
  }
  # else{
  #   states <- ntile.ts(x,n=ns,limit.type='manual',limits=limits,tie=tie)
  # }

  st <- statetable.msm(state,data=list(state=states))
  st/apply(st,1,sum)

} #end function

#Convert month strings to numeric
# Standardize input to match month names or abbreviations
match_month <- function(month) {
  month <- tolower(month)
  match <- match(tolower(substr(month, 1, 3)), tolower(month.abb))
  return(match)
}

#get days in month of any start and end month sequence
days_in_months <- function(sd, ed) {
  # Define days in each month for a leap year
  days_in_month_leap <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

  # Ensure sd and ed are numeric
  if(is.character(sd)) sd <- match_month(sd)
  if(is.character(ed)) ed <- match_month(ed)

  # Handle cases where ed < sd (e.g., water year October to September)
  if(ed < sd) ed <- ed + 12

  # Create a sequence of months from sd to ed
  months_seq <- (sd:ed) %% 12
  months_seq[months_seq == 0] <- 12

  # Extract the relevant days from the leap year vector
  days_in_month <- days_in_month_leap[months_seq]

  return(days_in_month)
}

# Find a complete 366-day period starting on the simulation start month.
get_reference_366_start <- function(dat.d, smo, context = "weather generator") {
  candidates <- which(dat.d$month == smo & dat.d$day == 1)

  for (idx in candidates) {
    idx2 <- idx + 365
    if (idx2 > nrow(dat.d)) next

    dates <- as.Date(dat.d$date[idx:idx2])
    if (any(is.na(dates))) next
    if (!all(as.numeric(diff(dates)) == 1)) next
    if (!any(dat.d$month[idx:idx2] == 2 & dat.d$day[idx:idx2] == 29)) next

    return(idx)
  }

  stop(paste0(
    "wxgenR requires at least one complete 366-day period starting in month ",
    smo, " for ", context,
    ". Include a leap year in the training period, or adjust `syr`, `eyr`, `smo`, and `emo`."
  ), call. = FALSE)
}

get_sim_season <- function(Xseas, row, col) {
  rows <- c(row, row - 1, row + 1)
  rows <- rows[rows >= 1 & rows <= nrow(Xseas)]
  seasons <- Xseas[rows, col]
  seasons <- seasons[!is.na(seasons)]

  if (length(seasons) == 0) return(NA_integer_)
  as.integer(seasons[1])
}

wxgenR_fun_messages_enabled <- function() {
  isTRUE(getOption("wxgenR.funMessages", interactive()))
}

wxgenR_fun_message <- function(stage) {
  if (!wxgenR_fun_messages_enabled()) return(invisible(NULL))

  messages <- list(
    wx = c(
      "wxgenR: weather generated. The atmosphere has signed off on the paperwork.",
      "wxgenR: simulation complete. Clouds have been persuaded into matrix form.",
      "wxgenR: done. The stochastic weather machine is cooling down.",
      "wxgenR: weather generated. A cosmic gumbo of precipitation, temperature, and seasonality."
    ),
    writeSim = c(
      "wxgenR: simulations written. The files are wearing tiny hard hats.",
      "wxgenR: output saved. Your traces have left the building.",
      "wxgenR: files complete. The commas behaved themselves."
    ),
    multisite_shuffle = c(
      "wxgenR: multisite shuffle complete. Station ranks have changed seats politely.",
      "wxgenR: shuffle done. Spatial correlation has entered the chat.",
      "wxgenR: multisite results ready. The stations are now in coordinated formation.",
      "wxgenR: multisite shuffle complete. The station network is kind of a cosmic gumbo."
    ),
    generate_TmaxTmin = c(
      "wxgenR: Tmax and Tmin generated. The diurnal range has found its lane.",
      "wxgenR: temperature post-processing complete. Maximum and minimum are on speaking terms.",
      "wxgenR: Tmax/Tmin complete. The thermometer has been briefed."
    )
  )

  stage_messages <- messages[[stage]]
  if (length(stage_messages) == 0) return(invisible(NULL))

  tick <- floor((as.numeric(Sys.time()) %% 86400) * 1000)
  msg <- stage_messages[(tick %% length(stage_messages)) + 1]
  message(msg)

  invisible(NULL)
}
