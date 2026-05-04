#' Generate maximum and minimum daily temperature from daily average temperature
#'
#' Estimate maximum and minimum temperature from `wx` simulation results of daily average temperature. This function uses quantile mapping
#' to develop a relationship between daily average temperature and daily maximum temperature. Average and maximum temperature can be contained in your training
#' data for example. Then, the relationship between average and maximum temperature is applied to the simulated daily average temperature output by `wx` to estimate
#' a simulated daily maximum temperature. Daily minimum temperature is then estimated using the standard equation Tavg = 0.5(Tmax + Tmin).
#' \cr\cr
#' Quantile mapping is trained and applied separately by calendar month.
#' The quantile mapping function comes from the `qmap` R package, specifically using empirical quantiles.
#' \cr\cr
#'  You will likely need to post-process the `wx` function output in order
#'  to process each simulation trace since the `wx` output is in wide format by trace
#'   and this function is best applied in long format. Code examples for processing both the training data and simulated data
#'   to put into suitable  format for the `generate_TmaxTmin` function are available here:\cr\cr
#'   https://github.com/dwoodson-usbr/wxgenR/tree/master/vignettes
#'
#' @param df.train A dataframe containing the daily training data used to develop the relationship
#'  between average and maximum temperature. At a minimum should contain variables named `tmax`, `temp`, and `month`
#'  which represent maximum temperature, average temperature, and numeric month in which the daily observation exists, respectively.
#' @param df.sim A dataframe containing the daily simulation results from the `wx` function in long format. At a minimum should contain
#'  variables named `sim_temp` and `month` which respectively represent the simulated daily average temperature from `wx`
#'  and the numeric month in which each data point resides.
#'
#'
#'
#' @return Returns a list containing results and metadata from the multisite shuffling in 'long' format for easy analysis and visualization.
#' \itemize{
#'   \item df.sim - Dataframe containing simulated maximum and minimum temperature generated from daily average temperature.
#'   \item qmap.monthly - List of dataframes containing fitted quantile mapping models for each calendar month.
#' }
#'
#' @examples
#' \donttest{
#' # Example with toy data
#' df.train = data.frame(
#'   temp = runif(10, 50, 70),
#'   tmax = runif(10, 60, 80),
#'   month = rep(1:2, each = 5)
#' )
#'
#' df.sim = data.frame(
#'   sim_temp = runif(10, 50, 70),
#'   month = rep(1:2, each = 5)
#' )
#'
#' result = generate_TmaxTmin(df.train, df.sim)
#' head(result$df.sim)
#' }

#'
#' @references
#'
#' Gudmundsson L (2025). qmap: Statistical transformations for post-processing climate model output. R package version 1.0-6.
#'
#' Gudmundsson L, Bremnes JB, Haugen JE, Engen-Skaugen T (2012). “Technical Note: Downscaling RCM precipitation to the station scale using statistical transformations - a comparison of methods.” Hydrology and Earth System Sciences, 16(9), 3383–3390. doi:10.5194/hess-16-3383-2012.
#'
#' @export
#'
#' @importFrom dplyr mutate filter
#' @importFrom qmap fitQmapQUANT doQmapQUANT
#' @import magrittr
#'
#'

"generate_TmaxTmin" = function(df.train, df.sim){

  df.sim = df.sim %>%
    mutate(sim_tmax = NA, sim_tmin = NA)

  qmap.monthly = list()

  mo = 1
  for(mo in 1:12){

    df.mo = df.train %>%
      filter(month == mo)

    qmap.fit = fitQmapQUANT(obs = df.mo$tmax, mod = df.mo$temp, wet.day	= F)

    qmap.monthly[[mo]] = qmap.fit

    df.sim = df.sim %>%
      mutate(
        sim_tmax = ifelse(month == mo,
                          doQmapQUANT(x = sim_temp, fobj = qmap.fit, type = "linear"), sim_tmax)
      )

  }

  df.sim = df.sim %>%
    mutate(sim_tmin = 2*sim_temp-sim_tmax)

  # Identify rows where sim_temp is not between sim_tmin and sim_tmax
  invalid_rows = df.sim %>%
    filter(!(sim_tmax > sim_temp & sim_temp > sim_tmin))

  # Check if there are any violations and print appropriate message
  if (nrow(invalid_rows) == 0) {
    message("tmin and tmax simulation complete: No temperature violations found.")
  } else {
    message("Temperature violations found: Displaying affected rows.")
    print(invalid_rows)  # Adjust 'n' as needed to display more rows
  }

  wxgenR_fun_message("generate_TmaxTmin")

  return(list(df.sim = df.sim, qmap.monthly = qmap.monthly))

}
