#' Multisite shuffling of wxgenR simulation results
#'
#' Runs a postprocessor on `wx` simulation results to shuffle
#' multiple stations' simulations such that wxgenR can be used in multisite applications.
#' Specifically, the `multisite_shuffle` function uses an approach developed by Iman and Conover (1982) and later applied by Clark et al. (2004)
#' to capture the rank correlation among observed station data and introduce it to those stations' simulations. The Iman and Conover approach is implemented
#' using the `cornode` function from the `mc2d` package.
#' Note that `mc2d::cornode()` preserves the first variable passed to the
#' shuffling routine as an anchor. In `multisite_shuffle()`, variables are
#' ordered by station, with precipitation before temperature, so the first
#' station's precipitation series may remain unchanged between the unshuffled
#' and shuffled outputs. This behavior is expected and does not indicate a
#' failed shuffle.
#'
#' @param list.sta.wx A list containing multiple stations' observed and simulated data to be used in multisite shuffling.
#' Each list element should be named for the station's data it holds and should contain the dataframe output from the `wx` for that station.
#' @param numbCores Enable parallel computing for multisite shuffling,
#'  set number of cores to enable (must be a positive integer greater than or equal to 2).
#'   Turned off by default; if set to 0 or 1 it will run as single thread.
#'    Use function 'detectCores()' from 'parallel' package to show the number of available cores on your machine.
#' @param aseed Specify a seed for reproducibility.
#'
#' @return Returns a list containing results and metadata from the multisite shuffling in 'long' format for easy analysis and visualization.
#' \itemize{
#'   \item shuffledResultsOnly - Dataframe containing shuffled simulations for all traces and stations ('sim_prcp' and 'sim_temp').
#'   \item shuffledResultsAndObs - Dataframe containing shuffled simulations as well as corresponding observations/training data ('prcp' and 'temp' are the observed data).
#'   \item shuffledAndUnshuffled - Dataframe containing shuffled simulations, unshuffled simulations, observations. Unshuffled vs and shuffled are indicated by the 'Tag' variable.
#'   \item nameMap - Dataframe mapping station IDs and variables to internal names used during shuffling.
#' }
#'
#' @examples
#' \donttest{
#' # Simulated example with two stations
#'
#' data(BOCO_sims)
#'
#' ms = multisite_shuffle(BOCO_sims, numbCores = 2, aseed = 123)
#'
#' print(ms)
#' }
#'
#' @references
#' Iman, Ronald & Conover, William. (1982). A Distribution-Free Approach to Inducing Rank Correlation Among Input Variates. Communications in Statistics-simulation and Computation - COMMUN STATIST-SIMULAT COMPUT. 11. 311-334. 10.1080/03610918208812265.
#'
#' Clark, M., Gangopadhyay, S., Hay, L., Rajagopalan, B., & Wilby, R. (2004). The Schaake Shuffle: A Method for Reconstructing Space-Time Variability in Forecasted Precipitation and Temperature Fields. Journal of Hydrometeorology, 5(1), 243-262. https://doi.org/10.1175/1525-7541(2004)005<0243:TSSAMF>2.0.CO;2
#'
#' R. Pouillot, M.-L. Delignette-Muller (2010), Evaluating variability and uncertainty in microbial quantitative risk assessment using two R packages. International Journal of Food Microbiology. 142(3):330-40
#'
#' @export
#'
#' @importFrom dplyr mutate filter bind_rows bind_cols arrange left_join group_by summarise select
#' @importFrom tidyr pivot_longer pivot_wider starts_with
#' @importFrom lubridate ymd leap_year
#' @importFrom magrittr %>%
#' @importFrom mc2d cornode
#' @importFrom foreach foreach %dopar%
#' @importFrom doRNG registerDoRNG
#' @import parallel
#' @import doParallel
#'
#'

"multisite_shuffle" = function(list.sta.wx, numbCores = NULL, aseed = NULL){

  stations = names(list.sta.wx)

  if (is.null(stations) || any(is.na(stations)) || any(stations == "")) {
    stop("list.sta.wx must be a named list. Each list element name is treated as the station ID.", call. = FALSE)
  }

  if (anyDuplicated(stations) > 0) {
    stop("Station names in list.sta.wx must be unique.", call. = FALSE)
  }

  if (!is.null(aseed)) set.seed(aseed)

  make_name_map = function(stations,
                           hist_vars = c("prcp", "temp"),
                           sim_vars = c("sim_prcp", "sim_temp")) {

    if (length(hist_vars) != length(sim_vars)) {
      stop("hist_vars and sim_vars must have the same length and matching order.", call. = FALSE)
    }

    out = do.call(rbind, lapply(seq_along(stations), function(i) {
      data.frame(
        station = stations[i],
        var_index = seq_along(hist_vars),
        hist_variable = hist_vars,
        sim_variable = sim_vars,
        internal_name = paste0("V", sprintf("%05d", ((i - 1) * length(hist_vars)) + seq_along(hist_vars))),
        stringsAsFactors = FALSE
      )
    }))

    rownames(out) = NULL
    out
  }

  make_positive_definite_cor = function(x, eps = 1e-06) {
    if (is.null(x)) {
      return(matrix(numeric(0), nrow = 0, ncol = 0))
    }

    x = as.matrix(x)
    dn = dimnames(x)
    n = ncol(x)

    if (n == 0) return(x)
    if (n == 1) {
      x[1, 1] = 1
      return(x)
    }

    x[!is.finite(x)] = 0
    x = (x + t(x)) / 2
    diag(x) = 1

    e = eigen(x, symmetric = TRUE)
    values = pmax(e$values, eps)
    x = e$vectors %*% diag(values, nrow = length(values)) %*% t(e$vectors)
    x = (x + t(x)) / 2
    x = stats::cov2cor(x)
    x[!is.finite(x)] = 0
    diag(x) = 1

    min_eigen = min(eigen(x, symmetric = TRUE, only.values = TRUE)$values)
    if (min_eigen <= eps / 10) {
      x = x + diag(abs(min_eigen) + eps, nrow = n)
      x = stats::cov2cor(x)
      diag(x) = 1
    }

    dimnames(x) = dn
    x
  }

  safe_spearman_cor = function(x) {
    x = as.data.frame(x)
    n = ncol(x)
    nms = colnames(x)

    if (n == 0) {
      out = matrix(numeric(0), nrow = 0, ncol = 0)
      return(out)
    }

    if (n == 1 || nrow(x) < 2) {
      out = diag(1, nrow = n)
      dimnames(out) = list(nms, nms)
      return(out)
    }

    out = suppressWarnings(stats::cor(x, method = "spearman", use = "pairwise.complete.obs"))

    if (is.null(out)) {
      out = diag(1, nrow = n)
      dimnames(out) = list(nms, nms)
      return(out)
    }

    out = as.matrix(out)
    out[!is.finite(out)] = 0
    diag(out) = 1
    dimnames(out) = list(nms, nms)
    make_positive_definite_cor(out)
  }

  shuffle_with_target = function(sim_data, target_cor, aseed) {
    sim_data = as.data.frame(sim_data)
    nms = colnames(sim_data)
    out = sim_data

    if (is.null(target_cor) || nrow(target_cor) == 0 || ncol(target_cor) == 0) {
      target_cor = diag(1, nrow = length(nms))
      dimnames(target_cor) = list(nms, nms)
    } else {
      target_cor = target_cor[nms, nms, drop = FALSE]
    }

    varying = vapply(sim_data, function(x) {
      x = x[is.finite(x)]
      length(unique(x)) > 1
    }, logical(1))

    if (sum(varying) >= 2 && nrow(sim_data) >= 2) {
      target_active = make_positive_definite_cor(target_cor[varying, varying, drop = FALSE])

      shuffled = tryCatch(
        {
          invisible(capture.output({
            shuffled_result = mc2d::cornode(as.matrix(sim_data[, varying, drop = FALSE]),
                                            target = target_active, result = TRUE,
                                            outrank = FALSE, seed = aseed)
          }))
          shuffled_result
        },
        error = function(e) {
          warning(paste0(
            "cornode failed with sampled target correlation matrix; retrying with an identity target. Original error: ",
            conditionMessage(e)
          ), call. = FALSE)

          identity_target = diag(1, nrow = sum(varying))
          dimnames(identity_target) = dimnames(target_active)

          tryCatch(
            {
              invisible(capture.output({
                shuffled_result = mc2d::cornode(as.matrix(sim_data[, varying, drop = FALSE]),
                                                target = identity_target, result = TRUE,
                                                outrank = FALSE, seed = aseed)
              }))
              shuffled_result
            },
            error = function(e2) {
              warning(paste0(
                "cornode failed with an identity target; leaving this monthly block unshuffled. Error: ",
                conditionMessage(e2)
              ), call. = FALSE)
              as.matrix(sim_data[, varying, drop = FALSE])
            }
          )
        }
      )

      shuffled = as.data.frame(shuffled)
      colnames(shuffled) = nms[varying]
      out[, varying] = shuffled
    }

    list(data = out, cor = safe_spearman_cor(out))
  }

  shuffle_one_trace = function(sim, df.combo, stations, years, list.hist.cor, name_map, aseed) {
    message(sim)

    list.shuff.cor = list()
    df.shuff = NULL

    df.sim = df.combo %>%
      filter(simulation == sim)

    df.combo.sim = NULL

    for(sta in stations){
      sim_map = name_map[name_map$station == sta, ]
      sim_map = sim_map[order(sim_map$var_index), ]

      df.sta = df.sim %>%
        filter(station == sta)

      df.sta = df.sta[, sim_map$sim_variable, drop = FALSE]
      colnames(df.sta) = sim_map$internal_name

      df.combo.sim = bind_cols(df.combo.sim, df.sta)
    }

    metaData = df.sim %>%
      filter(station == stations[1]) %>%
      dplyr::select(date, year, month)

    df.combo.sim = df.combo.sim %>%
      bind_cols(metaData)

    for(yr in years){
      list.shuff.cor[[as.character(yr)]] = list()

      for(mo in 1:12){
        df.combo.sim.mo = df.combo.sim %>%
          filter(year == yr & month == mo) %>%
          dplyr::select(-c(date, year, month))

        if (nrow(df.combo.sim.mo) == 0) next

        df.dates.mo = df.combo.sim %>%
          filter(year == yr & month == mo) %>%
          dplyr::select(date)

        target.cor.yr = sample(years, 1)
        target.cor = list.hist.cor[[as.character(target.cor.yr)]][[as.character(mo)]]

        shuffled = shuffle_with_target(df.combo.sim.mo, target.cor, aseed)

        df.shuffle = data.frame(date = df.dates.mo$date, shuffled$data, check.names = FALSE)

        df.shuffle.long = df.shuffle %>%
          pivot_longer(cols = -date, names_to = "internal_name", values_to = "value") %>%
          left_join(name_map %>% dplyr::select(internal_name, variable = sim_variable, station),
                    by = "internal_name") %>%
          dplyr::select(-internal_name) %>%
          pivot_wider(names_from = "variable", values_from = "value") %>%
          mutate(simulation = sim, Tag = "Shuffled") %>%
          arrange(date, station)

        if (any(is.na(df.shuffle.long$station))) {
          stop("Internal multisite shuffle name mapping failed.", call. = FALSE)
        }

        df.shuff = bind_rows(df.shuff, df.shuffle.long)
        list.shuff.cor[[as.character(yr)]][[as.character(mo)]] = shuffled$cor
      }
    }

    list(df.shuff = df.shuff, list.shuff.cor = list.shuff.cor)
  }

  worker_env = new.env(parent = parent.env(environment()))
  worker_env$make_positive_definite_cor = make_positive_definite_cor
  worker_env$safe_spearman_cor = safe_spearman_cor
  worker_env$shuffle_with_target = shuffle_with_target
  worker_env$shuffle_one_trace = shuffle_one_trace

  environment(worker_env$make_positive_definite_cor) = worker_env
  environment(worker_env$safe_spearman_cor) = worker_env
  environment(worker_env$shuffle_with_target) = worker_env
  environment(worker_env$shuffle_one_trace) = worker_env

  make_positive_definite_cor = worker_env$make_positive_definite_cor
  safe_spearman_cor = worker_env$safe_spearman_cor
  shuffle_with_target = worker_env$shuffle_with_target
  shuffle_one_trace = worker_env$shuffle_one_trace

  name_map = make_name_map(stations)

  df.combo = NULL
  df.hist = NULL

  for(sta in stations){
    df = list.sta.wx[[sta]]

    data = df$dat.d
    psim = as.data.frame(df$Xpamt)
    tsim = as.data.frame(df$Xtemp)

    if (is.null(data) || is.null(psim) || is.null(tsim)) {
      stop(paste0("Station ", sta, " is missing dat.d, Xpamt, or Xtemp from wx() output."), call. = FALSE)
    }

    monthly_climatology = data %>%
      group_by(month) %>%
      summarise(monthly_temp = mean(temp, na.rm = TRUE), .groups = "drop")

    years = unique(data$year)
    non_leap_years = years[!leap_year(years)]

    dec31_season = data %>%
      filter(month == 12, day == 31) %>%
      dplyr::select(year, season)

    extra_na_rows = data.frame(
      year = non_leap_years,
      month = rep(12, length(non_leap_years)),
      day = rep(32, length(non_leap_years)),
      prcp = rep(NA_real_, length(non_leap_years)),
      temp = rep(NA_real_, length(non_leap_years)),
      stringsAsFactors = FALSE
    ) %>%
      left_join(dec31_season, by = "year")

    data_fixed = data[, c("year", "month", "day", "prcp", "temp", "season")] %>%
      bind_rows(extra_na_rows) %>%
      arrange(year, month, day)

    if (nrow(data_fixed) != nrow(psim) || nrow(data_fixed) != nrow(tsim)) {
      stop(paste0(
        "Row mismatch for station ", sta, ": padded observed data has ", nrow(data_fixed),
        " rows, Xpamt has ", nrow(psim), " rows, and Xtemp has ", nrow(tsim),
        " rows. nsim likely does not match the padded observed/training period used by multisite_shuffle()."
      ), call. = FALSE)
    }

    colnames(psim) = paste0("sim_prcp_", seq_len(ncol(psim)))
    colnames(tsim) = paste0("sim_temp_", seq_len(ncol(tsim)))

    df.c.p = cbind(data_fixed, psim)
    df.c.t = cbind(data_fixed, tsim)

    df.prcp.long = df.c.p %>%
      pivot_longer(cols = starts_with("sim_prcp_"),
                   names_to = "simulation",
                   values_to = "sim_prcp") %>%
      mutate(simulation = gsub("sim_prcp_", "sim_", simulation))

    df.temp.long = df.c.t %>%
      pivot_longer(cols = starts_with("sim_temp_"),
                   names_to = "simulation",
                   values_to = "sim_temp") %>%
      mutate(simulation = gsub("sim_temp_", "sim_", simulation))

    df.long = left_join(df.prcp.long, df.temp.long,
                        by = c("year", "month", "day", "prcp", "temp", "season", "simulation")) %>%
      left_join(monthly_climatology, by = "month") %>%
      mutate(temp = ifelse(is.na(temp) | is.nan(temp), monthly_temp, temp),
             sim_temp = ifelse(is.na(sim_temp) | is.nan(sim_temp), monthly_temp, sim_temp),
             prcp = ifelse(is.na(prcp) | is.nan(prcp), 0, prcp),
             sim_prcp = ifelse(is.na(sim_prcp) | is.nan(sim_prcp), 0, sim_prcp)
      ) %>%
      dplyr::select(-monthly_temp) %>%
      filter(!(month == 12 & day == 32)) %>%
      mutate(date = ymd(paste(year, month, day, sep = "-")),
             station = sta)

    missing_temp = sum(is.na(df.long$temp) | is.nan(df.long$temp))
    missing_sim_temp = sum(is.na(df.long$sim_temp) | is.nan(df.long$sim_temp))
    missing_prcp = sum(is.na(df.long$prcp) | is.nan(df.long$prcp))
    missing_sim_prcp = sum(is.na(df.long$sim_prcp) | is.nan(df.long$sim_prcp))

    if (missing_temp > 0 | missing_sim_temp > 0 | missing_prcp > 0 | missing_sim_prcp > 0) {
      warning("Missing values found:\n",
              if (missing_temp > 0) paste("temp:", missing_temp, "\n") else "",
              if (missing_sim_temp > 0) paste("sim_temp:", missing_sim_temp, "\n") else "",
              if (missing_prcp > 0) paste("prcp:", missing_prcp, "\n") else "",
              if (missing_sim_prcp > 0) paste("sim_prcp:", missing_sim_prcp, "\n") else "",
              call. = FALSE)
    } else {
      message(paste0(sta, " - All checks passed: No missing values in temp, sim_temp, prcp, or sim_prcp."))
    }

    df.combo = bind_rows(df.combo, df.long)

    hist_map = name_map[name_map$station == sta, ]
    hist_map = hist_map[order(hist_map$var_index), ]

    df.hist.sta = filter(df.long, simulation == "sim_1")
    df.hist.sta = df.hist.sta[, hist_map$hist_variable, drop = FALSE]
    colnames(df.hist.sta) = hist_map$internal_name

    metaData = df.long %>%
      filter(simulation == "sim_1") %>%
      dplyr::select(date, year, month)

    if (is.null(df.hist)) {
      df.hist = bind_cols(metaData, df.hist.sta)
    } else {
      if (!identical(df.hist$date, metaData$date)) {
        stop("All stations must have matching date sequences for multisite shuffling.", call. = FALSE)
      }
      df.hist = bind_cols(df.hist, df.hist.sta)
    }
  }

  df.combo$Tag = "Unshuffled"
  years = sort(unique(df.hist$year))

  list.hist.cor = list()

  for(yr in years){
    list.hist.cor[[as.character(yr)]] = list()

    for(mo in 1:12){
      df.hist.mo = df.hist %>%
        filter(year == yr & month == mo) %>%
        dplyr::select(-c(date, year, month))

      list.hist.cor[[as.character(yr)]][[as.character(mo)]] = safe_spearman_cor(df.hist.mo)
    }
  }

  if (is.null(numbCores) || !is.numeric(numbCores) || numbCores < 2) {
    numbCores = 1
  }

  simz = unique(df.combo$simulation)

  if(numbCores == 1){
    df.shuff = NULL
    list.shuff.cor = list()

    for(sim in simz){
      res = shuffle_one_trace(sim = sim, df.combo = df.combo, stations = stations,
                              years = years, list.hist.cor = list.hist.cor,
                              name_map = name_map, aseed = aseed)

      df.shuff = bind_rows(df.shuff, res$df.shuff)
      list.shuff.cor[[sim]] = res$list.shuff.cor
    }

  } else {
    available_cores = parallel::detectCores()

    if(!is.na(available_cores) && numbCores > available_cores){
      numbCores = max(1, available_cores - 1)
      message(paste0("numbCores is set above the available cores on your machine. Setting numbCores to ", numbCores, "."))
    }

    if(numbCores < 2){
      df.shuff = NULL
      list.shuff.cor = list()

      for(sim in simz){
        res = shuffle_one_trace(sim = sim, df.combo = df.combo, stations = stations,
                                years = years, list.hist.cor = list.hist.cor,
                                name_map = name_map, aseed = aseed)

        df.shuff = bind_rows(df.shuff, res$df.shuff)
        list.shuff.cor[[sim]] = res$list.shuff.cor
      }
    } else {
      cl = parallel::makeCluster(numbCores)
      cluster_open = TRUE
      on.exit(if (cluster_open) parallel::stopCluster(cl), add = TRUE)
      doParallel::registerDoParallel(cl)
      doRNG::registerDoRNG(aseed)

      results = foreach::foreach(sim = simz,
                                 .packages = c("dplyr", "tidyr", "mc2d", "magrittr")) %dorng% {
        shuffle_one_trace(sim = sim, df.combo = df.combo, stations = stations,
                          years = years, list.hist.cor = list.hist.cor,
                          name_map = name_map, aseed = aseed)
      }

      parallel::stopCluster(cl)
      cluster_open = FALSE

      df.shuff = bind_rows(lapply(results, function(x) x$df.shuff))
      list.shuff.cor = stats::setNames(lapply(results, function(x) x$list.shuff.cor), simz)
    }
  }

  df.shuff = data.frame(df.shuff)

  df.shuff.full = df.shuff %>%
    left_join(
      df.combo %>% dplyr::select(year, month, day, prcp, temp, season, date, station, simulation),
      by = c("date", "station", "simulation")
    ) %>%
    arrange(date, station)

  df.merge = bind_rows(df.combo, df.shuff.full)

  out = list(shuffledResultsOnly = df.shuff,
             shuffledResultsAndObs = df.shuff.full,
             shuffledAndUnshuffled = df.merge,
             shuffledCorrelations = list.shuff.cor,
             nameMap = name_map)

  wxgenR_fun_message("multisite_shuffle")

  return(out)
}
