# Peter Boyd
# CTMISD FUNCTIONS
# March 29, 2022


# Continuous Triggering Model Independent Stochastic Declustering ---------

#' Model Independent Stochastic Declustering
#'
#' This function uses nonparametric procedures to analyze a Hawkes
#' process in temporal or spatio-temporal domain, with or without marks
#' through the Model independent stochastic declustering algorithm.
#'
#' @param dates a vector of dates as "yyyy-mm-dd".
#' @param ref_date a date that marks the beginning of the process, taken to be minimum date.
#' @param lat a vector of latitudes, omit if not using spatial data
#' @param lon a vector of longitudes, omit if not using spatial data
#' @param lat_lim vector of minimum, maximum, and height of pixels for nonstationary background rate
#' @param lon_lim vector of minimum, maximum, and width of pixels for nonstationary background rate
#' @param marks a vecotr of marks, or magnitudes, omit if not using marked data
#' @param time_breaks a vector of cutoff values for temporal bins of time differences
#' @param space_breaks a vector of cutoff values for spatial bins of distance differences
#' @param mark_breaks a vector of cutoff values for magnitude bins
#' @param ref_date a date to serve as time 0, defaults to earliest observation
#' @param time_of_day character string that lists the time of day of events, as hour:minute:second
#' @param just_times TRUE or FALSE object. TRUE if \code{dates} object is a vector of only times,
#' indicating time elapsed since beginning of catalog.
#' @param time_unit character string that specifies the desired unit of time
#' @param n_ints integer defining the number of temporal and spatial intervals to be used
#' to estimate triggering functions
#' @param g_span smoothing parameter (span) of LOESS used in estimating g(t)
#' @param h_span smoothing parameter (span) of LOESS used in estimating h(s)
#' @param k_span smoothing parameter (span) of LOESS used in estimating k(m)
#' @param nonstat_br TRUE if using a nonstationary background rate
#' @param dist vector specifying the sampling distribution of intervals and required parameters
#' @param np integer defining the number of neighboring points to use for variable bandwidth
#' of kernel density smoothing for nonstationary background rate
#' @param time_len double between 0 and 1 defining the maximum temporal interval length as a proportion
#' of the temporal domain
#' @param space_len double between 0 and 1 defining the maximum spatial interval length as a proportion
#' of the spatial domain
#' @param dist_unit character string that specifies the desired unit of distance: meter, kilometer, or mile
#' @param stop_when scalar that serves as conversion criterion, 1e-3 as default
#' @param show_progress when TRUE, algorithm will print the iteration number and maximum
#' pairwise change in the probability matrix.
#'
#' @return Probability matrix \code{p0} containing the probabilities that event
#' \code{i} is an offspring of event \code{j}, \code{i > j}. Diagonal elements
#' represent the probability that event \code{i} is a background event.
#' @return \code{br} is the estimated background rate of the process
#' @return \code{perc_diag} is the proportion of mass lying on the diagonal of matrix \code{p0}
#' @return \code{perc_br} is the proportion of events in which the maximum probabilistic assignment
#' is as a background event
#' @return \code{n_iterations} is the number of iterations executed until convergence
#' @return \code{locs} is a data frame listing midpoint latitude, midpoint longitude,
#' x and y index, and background rate of the pixel each event lies in, for nonstationary background rate
#' @return \code{br_grid} is a data frame with the background rate and midpoint of each pixel
#' @return \code{g_integral} is the Reimann sum used to scale g(t) as a density
#' @return \code{h_integral} is the Reimann sum used to scale h(s) as a density
#' @return \code{fitg} is a list containing the output of the LOESS fit of g(t)
#' @return \code{fith} is a list containing the output of the LOESS fit of h(s)
#' @return \code{fitk} is a list containing the output of the LOESS fit of k(m)
#' @return \code{input} is a list of all inputs
#'
#' @example
#' data("hm.csv")
#' out = ctmisd(dates = hm$t,
#'              ref_date = "1999-10-16",
#'              lat = hm$lat,
#'              lon = hm$lon,
#'              marks = hm$m,
#'              just_times = TRUE,
#'              time_unit = "day",
#'              dist_unit = "mile",
#'              n_ints = 3000,
#'              g_span = 0.05,
#'              h_span = 0.20,
#'              k_span = 0.10,
#'              nonstat_br = TRUE,
#'              dist = c("beta", 1/2, 2),
#'              time_len = 0.005,
#'              space_len = 0.01,
#'              )

#' @export
ctmisd <- function(dates,
                   ref_date = 0,
                   lat = NA,
                   lon = NA,
                   lon_lim = NA,
                   lat_lim = NA,
                   marks = NA,
                   stopwhen = 1e-3,
                   time_of_day = NA,
                   just_times = FALSE,
                   time_unit = "day",
                   dist_unit = "mile",
                   n_ints = 1000,
                   g_span = 0.1,
                   h_span = 0.1,
                   k_span = 0.1,
                   nonstat_br = TRUE,
                   dist = c("beta", 0.5, 1),
                   np = 20,
                   time_len = 0.01,
                   space_len = 0.01,
                   show_progress = FALSE) {
  # ONE: data maintenance
  if(is.na(lat) == TRUE) {
    lat = rep(0, length(dates))
  }

  if(is.na(lon) == TRUE) {
    lon = rep(0, length(dates))
  }

  if(is.na(marks) == TRUE) {
    marks = rep(0, length(dates))
  }

  df = data.frame(dates, lat, lon, marks)
  # put dates in correct format
  if (just_times == FALSE) {
    if (class(dates) == "Date") {
      dates_clean = dates
    } else {
      dates_clean = lubridate::as_date(dates)
    }

    #ref_date = lubridate::mdy(ref_date)
    # create continuous time variable since beginning
    times = lubridate::time_length(lubridate::interval(ref_date, dates_clean), time_unit)
    if (is.na(time_of_day[1]) == TRUE) {
      # add random time of day on date if not given
      times = times + runif(length(times), 0, 1)
    } else {
      times_clean = lubridate::hms(time_of_day)
      h = lubridate::hour(times_clean) / 24
      m = lubridate::minute(times_clean) / 1440
      s = lubridate::second(times_clean) / 86400
      times = times + h + m + s
    }
  } else {
    #ref_date = lubridate::mdy(ref_date)
    times = dates
  }
  # ensure all events are correctly ordered
  df$times = times
  df = df[order(df$times), ]
  df = df[which(df$times > 0), ]

  # output ordered events
  lat = df$lat
  lon = df$lon
  marks = df$marks
  times = df$times
  dates = df$dates

  # TWO: establish time, space differences for each event
  time_mat = get_time(times)

  n = length(lat)
  dist_mat = matrix(data = NA,
                    nrow = n,
                    ncol = n)
  R = 6371e3
  # use Haversine formula for distances
  for (i in 1:n) {
    for (j in 1:n) {
      phi1 = lat[i] * pi / 180
      phi2 = lat[j] * pi / 180
      latdiff = (lat[j] - lat[i]) * pi / 180
      londiff = (lon[j] - lon[i]) * pi / 180

      a = sin(latdiff / 2) * sin(latdiff / 2) +
        cos(phi1) * cos(phi2) * sin(londiff / 2) * sin(londiff / 2)
      b = 2 * atan2(sqrt(a), sqrt(1 - a))
      d = R * b

      dist_mat[i, j] = d
      #print(c(i,j))
    }
  }

  # convert distances to correct units
  if (dist_unit == "mile") {
    dist_mat = dist_mat * 0.000621371
  } else if (dist_unit == "kilometer") {
    dist_mat = dist_mat * 0.001
  } else {
    dist_mat = dist_mat
  }


  # sort distances. easier computation of number of
  # points for nonstationary background rate
  dist_mat2 = matrix(data = NA,
                     nrow = n,
                     ncol = n)
  for (i in 1:n) {
    dist_mat2[i, ] = sort(dist_mat[i, ])
  }
  # FIX THIS LATER
  if (sum(lat) != 0) {
    dist_mat = ifelse(dist_mat == 0, 0.0001, dist_mat)
  }
  # THREE: nonstationary background rate
  # use xg, xh as vector of differences for regression
  time_mat[upper.tri(time_mat, diag = FALSE)] = NA
  dist_mat[upper.tri(dist_mat, diag = FALSE)] = NA
  diag(dist_mat) = 0
  xg = time_mat
  xg = as.vector(na.omit(as.vector(t(xg))))
  xh = dist_mat
  xh = as.vector(na.omit(as.vector(t(xh))))
  xk = marks

  # create spatial grid to contain events
  # ADD IF STATEMENT HERE
  if(is.na(lat_lim) == T) {
    lat_lim = c(min(lat) - 0.01,
                max(lat) + 0.01,
                (max(lat) + 0.01 - min(lat) - 0.01) / 10)
  }

  if(is.nat(lon_lim) == T) {
    lon_lim = c(min(lon)-0.01,
                max(lon)+0.01,
                (max(lon) +0.01 - min(lon) - 0.01) / 10)
  }
  x_grid = seq(lon_lim[1], lon_lim[2], lon_lim[3])
  y_grid = seq(lat_lim[1], lat_lim[2], lat_lim[3])
  x_grid[length(x_grid)] = lon_lim[2]
  y_grid[length(y_grid)] = lat_lim[2]
  # grid is entire window, pix is midpoints of grid
  x_pix = seq(lon_lim[1] + lon_lim[3] * 0.5,
              lon_lim[2] - lon_lim[3] * 0.5,
              lon_lim[3])
  y_pix = seq(lat_lim[1] + lat_lim[3] * 0.5,
              lat_lim[2] - lat_lim[3] * 0.5,
              lat_lim[3])
  # calculate radii such that np number of events are within
  # distance d_i
  di = calc_d(dist_mat2, np = np)# was doing 24

  # put all events into a pixel
  # select proper br in prob calcs
  pix = get_pix(x_grid, y_grid, lat, lon, x_pix, y_pix)
  pix = data.frame(pix)
  names(pix) = c("lon", "lat", "x", "y")

  # FOUR: random interval creation
  # create time intervals
  time_ints = matrix(nrow = n_ints, ncol = 4)
  dist_ints = matrix(nrow = n_ints, ncol = 4)

  for (i in 1:n_ints) {
    if (dist[1] == "beta") {
      x_time = as.numeric(rbeta(1, as.numeric(dist[2]), as.numeric(dist[3])) *
                            ceiling(max(times)))
      x_space = rbeta(1, as.numeric(dist[2]), as.numeric(dist[3])) *
        ceiling(max(dist_mat, na.rm = T))
    } else {
      x_time = runif(1, 0, 1) * ceiling(max(times))
      x_space = runif(1, 0, 1) * ceiling(max(dist_mat, na.rm = T))
    }
    t_diff = time_len * ceiling(max(times))
    y_time = x_time + runif(1,-1 * t_diff, t_diff)
    y_time = ifelse(y_time < 0, 0, y_time)
    time_ints[i, 1] = min(x_time, y_time)
    time_ints[i, 2] = max(x_time, y_time)
    time_ints[i, 4] = (time_ints[i, 2] + time_ints[i, 1]) / 2

    s_diff = space_len * ceiling(max(dist_mat, na.rm = T))
    y_space = x_space + runif(1,-1 * s_diff, s_diff)
    y_space = ifelse(y_space < 0, 0, y_space)
    dist_ints[i, 1] = min(x_space, y_space)
    dist_ints[i, 2] = max(x_space, y_space)
    dist_ints[i, 4] = (dist_ints[i, 2] + dist_ints[i, 1]) / 2
  }

  # create mark intervals
  mark_ints = matrix(nrow = n_ints, ncol = 4)
  # mark_quants = quantile(xk, probs = seq(0.1, 1, 0.1))
  # mark_ints = get_ints(min(xk), mark_quants, n_ints)
  # mark_ints[,4] = (mark_ints[,1] + mark_ints[,2]) / 2
  for (i in 1:n_ints) {
    x = rbeta(1, as.numeric(dist[2]), as.numeric(dist[3])) * (max(xk, na.rm = T))
    y = x + runif(1,-0.05, 0.05)
    y = ifelse(y < 3, 3, y)
    #x = runif(2, 0, 1)*ceiling(max(dist_mat, na.rm = T))
    mark_ints[i, 1] = min(x, y)
    mark_ints[i, 2] = max(x, y)
    mark_ints[i, 4] = (mark_ints[i, 2] + mark_ints[i, 1]) / 2
  }

  # FIVE: initialize prior to iteration
  p0 = init_p0(times)
  max_diff = 1
  n_iterations = 0
  md = c()

  # SIX: iterate
  #while( max_diff > stopwhen){
  while (max_diff > stopwhen) {
    # decide br method to use
    if (nonstat_br == TRUE) {
      tau = calc_tau(x_pix, y_pix, p0, di,
                     lon, lat, times)
      # z = calc_z(times, x_pix, y_pix, tau) # this function crashed I think
      z = sum(tau[, 1] * diff(x_pix)[1] * diff(y_pix)[1] * ceiling(max(times)))
      br = calc_br_nonstat(p0, times, z, tau)

      br_grid = data.frame(br)
      names(br_grid) = c("br", "lon", "lat")

      locs = dplyr::left_join(pix, br_grid, by = c("lon", "lat"))
      locs$br = ifelse(is.na(locs$br) == TRUE, 0, locs$br)
      br = as.vector(locs$br)
    } else {
      br = calc_br(p0, times)
      br = rep(br, n)
      br_grid = NA
      locs = NA
    }

    # SEVEN: calculate triggering components
    h = get_h(p0, dist_mat, dist_ints, n_ints)
    # get_g crashed, exporting division worked.
    g_vals = get_g(p0, time_mat, time_ints, n_ints)
    g = time_ints
    g[, 3] = g_vals[, 1] / g_vals[, 2]
    k = get_k(p0, marks, mark_ints, n_ints)

    # remove intervals containing 0 events
    g = na.omit(g)
    h = na.omit(h)
    k = na.omit(k)
    attributes(g)$na.action = NULL
    attributes(h)$na.action = NULL
    attributes(k)$na.action = NULL

    # change class and name
    g = as.data.frame(g)
    h = as.data.frame(h)
    k = as.data.frame(k)
    names(g) = c("min_int", "max_int", "est", "midpt")
    names(h) = c("min_int", "max_int", "est", "midpt")
    names(k) = c("min_int", "max_int", "est", "midpt")

    # round negative estimates to 0 (this shouldn't happen)
    g$est = ifelse(sign(g$est) == -1, 0, g$est)
    h$est = ifelse(sign(h$est) == -1, 0, h$est)
    k$est = ifelse(sign(k$est) == -1, 0, k$est)

    # EIGHT: fit smooth curve
    # Loess: direct allows for extrapolation outside defined observations
    # g
    fitg = loess(est ~ midpt,
                 data = g,
                 span = g_span,
                 surface = "direct")
    predg = predict(fitg, newdata = data.frame(midpt = xg))
    # h
    if (sum(lat) != 0) {
      fith = loess(est ~ midpt,
                   data = h,
                   span = h_span,
                   surface = "direct")
      predh = predict(fith, newdata = data.frame(midpt = xh))
    } else {
      predh = rep(1, length(xh))
    }
    # k
    if (sum(marks) != 0) {
      fitk = loess((colSums(p0) - diag(p0)) ~ xk,
                   span = k_span,
                   surface = "direct")
      predk = fitk$fitted
    } else {
      predk = rep(1, length(xk))
    }

    # round negative predictions up to 0
    predg = ifelse(sign(predg) == -1, 0, predg)
    predh = ifelse(sign(predh) == -1, 0, predh)
    predk = ifelse(sign(predk) == -1, 0, predk)

    # NINE: preserve density properties of g and h
    gspline = data.frame(midpt = xg, est = predg)
    hspline = data.frame(midpt = xh, est = predh)

    # use Reimann sum to approximate integral
    # divide predictions by integral
    gspline = gspline[order(gspline$midpt), ]
    gxdif =  c(diff(gspline$midpt), 0)
    g_integral = sum(gxdif * gspline$est, na.rm = TRUE) #+ sum(gxdif*g_spline$est)/2
    predg = predg / g_integral
    g$est = g$est / g_integral
    # probably best to do this...?

    sp = ifelse(sum(xh) == 0, 0, 1)
    if (sp == 1) {
      hspline = hspline[order(hspline$midpt), ]
      hxdif =  c(diff(hspline$midpt), 0)
      h_integral = sum(hxdif * hspline$est, na.rm = TRUE) #+ sum(gxdif*g_spline$est)/2
      predh = predh / h_integral
      h$est = h$est / h_integral
    }
    predk = predk*g_integral*h_integral

    # TEN: update probability matrix
    p = update_p(p0, dist_mat, br,
                 predg, predh, predk, sp)
    predg = data.frame(midpt = xg, est = predg)
    predh = data.frame(midpt = xh, est = predh)
    predk = data.frame(midpt = xk,
                       est = predk,
                       obs = colSums(p0) - diag(p0))
    max_diff = check_p(p0, p)
    p0 = p
    n_iterations = n_iterations + 1
    if(show_progress == TRUE) {
      print(paste("n: ", n_iterations, ". max diff: ", max_diff))
    }
    md[n_iterations] = max_diff
  }

  # ELEVEN: additional calculations
  max_diag = c()
  for (i in 1:nrow(p0)) {
    if (p0[i, i] == max(p0[i, ])) {
      max_diag[i] = 1
    } else{
      max_diag[i] = 0
    }
  }
  df$mainshock = max_diag

  max_event = c()
  for (i in 1:nrow(p0)) {
    max_event[i] = which.max(p0[i, ])
  }
  df$parent = max_event

  perc_br = sum(max_diag) / length(max_diag)
  perc_diag = sum(diag(p0)) / nrow(p0)

  out = list(
    p0 = p0,
    br = br,
    br_grid = br_grid,
    perc_br = perc_br,
    perc_diag = perc_diag,
    data = df,
    #predg = predg,
    #predh = predh,
    #predk = predk,
    #g = g,
    #h = h,
    #k = k,
    # y_pix = y_pix,
    # x_pix = x_pix,
    # pix = pix,
    n_iterations = n_iterations,
    dist_mat = dist_mat,
    time_mat = time_mat,
    g_integral = g_integral,
    h_integral = h_integral,
    fitg = fitg,
    fith = fith,
    fitk = fitk,
    locs = locs,
    input = mget(names(formals()), sys.frame(sys.nframe()))
  )
  return(out)
}


# Triggering Plot ---------------------------------------------------------

#' This function plots the smooth triggering function for each component,
#' g(t), h(s), and k(m) for time, space, and marks, respectively.
#'
#' @param model the output from the ctmisd() function.
#' @param g_xlim a vector of the minimum and maximum value on the x-axis of g(t)
#' @param h_xlim a vector of the minimum and maximum value on the x-axis of h(s)
#' @param k_xlim a vector of the minimum and maximum value on the x-axis of k(m)
#' @param g_ylim a vector of the minimum and maximum value on the y-axis of g(t)
#' @param h_ylim a vector of the minimum and maximum value on the y-axis of h(s)
#' @param k_ylim a vector of the minimum and maximum value on the y-axis of k(m)
#' @param mag_label a character for the name of the marks employed
#' @param se_include TRUE if standard errors are to be included on the plots
#' @param include_pts TRUE if points, in addition to the fitted curves, are to be included on the plots
#'
#' @return \code{all_plots} all triggering component studied are plotted
#' @return \code{time_plot} only g(t) is shown
#' @return \code{space_plot} only h(s) is shown
#' @return \code{mark_plot} only k(m) is shown
#'
#' @example
#' tp = trig_plots(out_smooth_final,
#'                 g_xlim = c(0,10),
#'                 k_xlim = c(3,7),
#'                 h_xlim = c(0,25),
#'                 include_pts = F,
#'                 se_include = T)

#' @export
trig_plots = function(model,
                      g_xlim = c(min(model$fitg$x), max(model$fitg$x)),
                      h_xlim = c(min(model$h$midpt), max(model$h$midpt)),
                      k_xlim = c(min(model$k$midpt), max(model$k$midpt)),
                      g_ylim = c(0, NA),
                      h_ylim = c(0, NA),
                      k_ylim = c(0, NA),
                      mag_label = "magnitude",
                      se_include = TRUE,
                      include_pts = F) {

  # make g(t) data frame
  gfit = predict(model$fitg,
                 newdata =
                   data.frame(midpt = model$predg$midpt),
                 se = TRUE)
  g_se = gfit$se.fit / (model$g_integral ^ 2)
  g_df = data.frame(
    midpt = model$predg$midpt,
    min = model$predg$est - 2 * g_se,
    max = model$predg$est + 2 * g_se,
    est = model$predg$est
  )

  # make h(s) data frame
  hfit = predict(model$fith,
                 newdata =
                   data.frame(midpt = model$predh$midpt),
                 se = TRUE)
  h_se = hfit$se.fit / (model$h_integral ^ 2)
  h_df = data.frame(
    midpt = model$predh$midpt,
    min = model$predh$est - 2 * h_se,
    max = model$predh$est + 2 * h_se,
    est = model$predh$est
  )

  # make k(m) data frame
  kfit = predict(model$fitk,
                 newdata =
                   data.frame(midpt = model$data$marks),
                 se = TRUE)
  k_se = kfit$se.fit
  k_df = data.frame(
    midpt = model$data$marks,
    obs = model$predk$est - 2 * k_se,
    max = model$predk$est + 2 * k_se,
    est = model$predk$est
  )

  # Time Plot
  trig_g = ggplot2::ggplot(model$g, ggplot2::aes(x = midpt, y = est))

  if(se_include == TRUE){
    trig_g = trig_g + ggplot2::geom_ribbon(data = g_df,
                         ggplot2::aes(x = midpt, ymin = min, ymax = max),
                         color = "lightgray",
                         fill = "lightgray")
  }
  predg = data.frame(model$fitg$x, model$fitg$y)
  trig_g = trig_g +
    ggplot2::xlab(paste("t (time in ", model$input$time_unit, "s)", sep = "")) +
    ggplot2::ylab("g(t)") +
    ggplot2::geom_line(data = model$predg, ggplot2::aes(x = midpt, y = est)) +
    # ggplot2::geom_smooth(data = model$predg, ggplot2::aes(x = midpt, y = est),
    #                     method = "loess") +
    ggplot2::coord_cartesian(xlim = g_xlim, ylim = g_ylim) +
    ggplot2::theme(
      axis.text.x =
        ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor.x = element_blank()
    )

  if(include_pts == TRUE) {
    trig_g = trig_g + ggplot2::geom_point(alpha = 0.2)
  }

  # Magnitude Plot
  if (sum(model$data$marks) != 0) {
    trig_k = ggplot2::ggplot(model$predk, ggplot2::aes(x = midpt, y = obs))
    if(se_include == TRUE) {
      trig_k = trig_k +
        ggplot2::geom_ribbon(data = k_df,
                             ggplot2::aes(
                              x = midpt,
                              ymin = min,
                              ymax = max),
                             color = "lightgray",
                             fill = "lightgray")
    }
    trig_k = trig_k +
      ggplot2::xlab(paste("m (", mag_label, ")", sep = "")) +
      ggplot2::ylab("k(m)") +
      ggplot2::geom_line(data = model$predk, ggplot2::aes(x = midpt, y = est)) +
      ggplot2::coord_cartesian(xlim = k_xlim, ylim = k_ylim) +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid.minor.x = element_blank()
      )

    if(include_pts == TRUE) {
      trig_k = trig_k + ggplot2::geom_point(alpha = 0.2)
    }

  } else {
    trig_k = NULL
  }

  #Space Plot
  if (sum(model$dist_mat, na.rm = TRUE) != 0) {
    trig_h = ggplot2::ggplot(model$h, ggplot2::aes(x = midpt, y = est))
    if (se_include == TRUE) {
      trig_h = trig_h +
        ggplot2::geom_ribbon(data = h_df,
                             ggplot2::aes(
                               x = midpt,
                               ymin = min,
                               ymax = max
                             ),
                             color = "lightgray",
                             fill = "lightgray")

    }
    trig_h = trig_h +
      ggplot2::xlab(paste("s (distance in ", model$input$dist_unit, "s)", sep = "")) +
      ggplot2::ylab("h(s)") +
      ggplot2::geom_line(data = model$predh, ggplot2::aes(x = midpt, y = est)) +
      ggplot2::coord_cartesian(xlim = h_xlim, ylim = h_ylim) +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid.minor.x = element_blank()
      )

    if(include_pts == TRUE) {
      trig_h = trig_h + ggplot2::geom_point(alpha = 0.2)
    }

  } else{
    trig_h = NULL
  }

  # Output
  if (sum(model$dist_mat, na.rm = TRUE) == 0 &
      sum(model$data$marks) == 0) {
    all_plots = gridExtra::grid.arrange(trig_g, ncol = 1)
  } else if (sum(model$data$marks) != 0 &
             sum(model$dist_mat, na.rm = TRUE) == 0) {
    all_plots = gridExtra::grid.arrange(trig_g, trig_k, ncol = 2)
  } else if (sum(model$data$marks) == 0 &
             sum(model$dist_mat, na.rm = TRUE) != 0) {
    all_plots = gridExtra::grid.arrange(trig_g, trig_h, ncol = 2)
  } else {
    all_plots = gridExtra::grid.arrange(trig_g, trig_h, trig_k, ncol = 3)
  }
  out = list(
    all_plots = all_plots,
    time_plot = trig_g,
    space_plot = trig_h,
    mark_plot = trig_k
  )
  return(out)
}

#' Conditional Intensity
#'
#' This function calculates the conditional intensity using CTMISD.
#'
#' @param model the output from the ctmisd function.
#'
#' @return a data frame containing the time, date, latitude, longitude, mark, and
#' conditional intensity at each observed point.
#'
#' @export
cond_int_smooth = function(model) {
  times = model$data$times
  marks = model$data$marks
  lat = model$data$lat
  lon = model$data$lon

  n = nrow(model$p0)
  cond_int = c()
  # calculate br for obs 1 based on all data
  cond_int[1] = sum(diag(model$p0)) / (max(times) - min(times))

  # get overall trig component per event
  kvals = c()
  for (i in 1:n) {
    repk = rep(model$predk$est[i], i)
    kvals = c(kvals, repk)
  }
  trig = model$predg$est * model$predh$est * kvals
  br = model$br
  # get br and conditional intensity
  for (i in 2:n) {
    b = sum(1:i)#vector location of event pair (i,i)
    a = b - i + 1#vector location of event pair(1,i)
    cond_int[i] = br[i] + sum(trig[a:(b - 1)])
  }

  ci = data.frame(times, lat, lon, marks, cond_int)
  ci$dates = lubridate::ymd(model$input$ref_date) +
    lubridate::days(floor(ci$times))
  # ci$Date = as.Date(model$ref_date + ci$times)
  return(ci)
}

#'
cond_int_bin = function(model) {
  times = model$data$times
  marks = model$data$marks
  lat = model$data$lat
  lon = model$data$lon

  # binning function
  time_breaks = model$time_breaks
  space_breaks = model$space_breaks
  mark_breaks = model$mark_breaks

  bin_f <- function(u, v) {
    x <- rep(0, length(u))
    for (j in 1:length(u)) {
      for (i in 1:(length(v) - 1)) {
        if (v[i] < u[j] & u[j] <= v[i + 1]) {
          x[j] <- i
        }
      }
    }
    return(x)
  }

  cond_int = c()
  # calculate br for obs 1 based on all data
  cond_int[1] = sum(diag(model$p0)) / (max(times) - min(times))

  n = length(times)
  g_values = matrix(nrow = n,
                    ncol = n,
                    data = rep(0, n * n))
  k_values = matrix(nrow = n,
                    ncol = n,
                    data = rep(0, n * n))
  h_values = matrix(nrow = n,
                    ncol = n,
                    data = rep(0, n * n))

  # create matrix of bin values for each obs in
  # 3 components. add 1 to convert from c++
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      g_values[i, j] = model$g[model$time_bins[j, i] + 1]
    }
  }

  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      h_values[i, j] = model$h[model$dist_bins[j, i] + 1]
    }
  }

  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      k_values[i, j] = model$k[model$mark_bins[j] + 1]
    }
  }

  # get overall trig component per event
  trig = g_values * k_values * h_values
  br = model$br
  # get br and conditional intensity
  for (i in 2:n) {
    cond_int[i] = br[i] + sum(trig[i,])
  }

  ci = data.frame(times, lat, lon, marks, cond_int)
  ci$dates = lubridate::ymd(model$ref_date) +
    lubridate::days(floor(ci$times))
  return(ci)
}


# Background Rate Plot ----------------------------------------------------

#' Background Rate Plot
#'
#' Returns a heat map showing the nonstationary background rate over the spatial domain
#'
#' @param model is object exported from ctmisd function
#' @param colors is a vector of two colors representing low and high values
#' @param plot_title optional title to include
#'
#' @return a ggplot object showing the nonstationary background rate.
#'
#' @example
#' bp = br_plot(out)
#'

#' @export
br_plot = function(model, colors = c("blue", "red"), plot_title = NULL) {
  br_grid = model$br_grid
  ggplot2::ggplot(data = br_grid, aes(x = lon, y = lat, fill = br)) +
    ggplot2::scale_fill_gradient(low = colors[1], high = colors[2]) +
    ggplot2::geom_tile() +
    ggplot2::ggtitle(plot_title) +
    ggplot2::scale_x_continuous(minor_breaks = NULL) +
    ggplot2::scale_y_continuous(minor_breaks = NULL)
}


# Deviance Residual -------------------------------------------------------

#' Calculates the deviance residuals of a nonparametric Hawkes process model fit
#' using the misd or ctmisd function. May be used compare the performance of two
#' competing models.
#'
#' @param model1 the output from the misd or ctmisd function.
#' @param model2 the output from an optional second model fit with misd or ctmisd
#' @param time_pts a vector of points spanning the temporal domain.
#' The conditional intensity will be evaluated at each point.
#' @param bin1 TRUE if model 1 was estimated using nphawkes::misd, rather than ctmisd
#'
#' @return a data frame of each pixel's x and y coordinate, latitude and longitude midpoint,
#' and deviance residual
#'
#' @example
#' dr = dev_res(out)
#'

#' @export
dev_res = function(model1,
                    model2 = NULL,
                    time_pts = quantile(1:ceiling(max(model1$data$times)),
                                        probs = seq(0, 1, 0.05)),
                    bin1 = TRUE) {
  m = nrow(model1$br_grid)
  n = nrow(model1$locs)
  names(model1$locs) = c("lon_mid", "lat_mid", "x", "y", "br")
  if (bin1 == TRUE) {
    ci1 = ci_full_bin(model1, time_seq = time_pts)
    ci_obs1 = cond_int_bin(model1)
    ci_obs1 = dplyr::select(ci_obs1,-marks,-dates)
  } else {
    ci1 = ci_full(model1, time_seq = time_pts)
    ci_obs1 = cond_int_smooth(model1)
    ci_obs1 = dplyr::select(ci_obs1,-marks,-dates)
  }
  ci_obs1$x = model1$locs$x
  ci_obs1$y = model1$locs$y
  ci_obs1$br = model1$locs$br
  ci_obs1 = rename(ci_obs1, t = times)
  ci_obs1 = rename(ci_obs1, ci = cond_int)
  ci_obs1$obs = 1
  ci1$obs = 0

  # iterate over each pixel

  for (i in (min(ci1$x)):(max(ci1$x))) {
    for (j in (min(ci1$y)):(max(ci1$y))) {
      # if no points in a pixel,
      if (nrow(ci_obs1[which(ci_obs1$x == i &
                             ci_obs1$y == j), ]) == 0) {
        ci_obs1[nrow(ci_obs1) + 1, ] =
          c(
            ceiling(max(model1$data$times)),
            #t
            model1$y_pix[j + 1],
            model1$x_pix[i + 1],
            #lat, lon
            model1$br_grid[(10 * i + j + 1), 1],
            #ci assumed to be br
            # computation assumption
            i,
            j,
            #x,y
            model1$br_grid[(10 * i + j + 1), 1],
            1
          )#br, obs
      }
    }
  }
  ci1_full = rbind(ci_obs1, ci1)
  # THINK I NEED TO ORDER EVENTS BY TIME
  ci1_full$logc = log(ci1_full$ci)
  ci1_full = ci1_full[order(ci1_full$t), ]

  dr_int1 = ci1_full %>%
    dplyr::group_by(x, y) %>%
    dplyr::arrange(t) %>%
    dplyr::mutate(ci_ints = diff(c(0, t)) * ci) %>%
    dplyr::summarise(ints1 = sum(ci_ints))

  dr_sum1 = ci1_full %>%
    dplyr::filter(obs == 1) %>%
    dplyr::group_by(x, y) %>%
    dplyr::arrange(t) %>%
    dplyr::summarise(sums1 = sum(logc))

  dr1 = dplyr::left_join(dr_int1, dr_sum1, by = c("x", "y")) %>%
    dplyr::mutate(res1 = sums1 - ints1)
  # does sum log account for empty pixels?

  if (is.null(model2) == FALSE) {
    m = nrow(model2$br_grid)
    n = nrow(model2$locs)
    names(model2$locs) = c("lon_mid", "lat_mid", "x", "y", "br")
    ci2 = ci_full(model2, time_seq = time_pts)
    ci_obs2 = cond_int_smooth(model2)
    ci_obs2$x = model2$locs$x
    ci_obs2$y = model2$locs$y
    ci_obs2$br = model2$locs$br
    ci_obs2 = dplyr::rename(ci_obs2, t = times)
    ci_obs2 = dplyr::rename(ci_obs2, ci = cond_int)
    ci_obs2 = dplyr::select(ci_obs2,-marks)
    ci_obs2$obs = 1
    ci2$obs = 0

    for (i in (min(ci2$x)):(max(ci2$x))) {
      for (j in (min(ci2$y)):(max(ci2$y))) {
        if (nrow(ci_obs2[which(ci_obs2$x == i &
                               ci_obs2$y == j), ]) == 0) {
          ci_obs2[nrow(ci_obs2) + 1, ] =
            c(
              ceiling(max(model2$data$times)),
              model2$y_pix[j + 1],
              model2$x_pix[i + 1],
              model2$br_grid[(10 * i + j + 1), 1],
              i,
              j,
              model2$br_grid[(10 * i + j + 1), 1],
              1
            )
        }
      }
    }
    ci2_full = rbind(ci_obs2, ci2)
    ci2_full$logc = log(ci2_full$ci)
    ci2_full = ci2_full[order(ci2_full$t), ]
    ci2_full$prod = NA
    ci2_full$prod[1] = ci2_full$ci[1] * ci2_full$t[1]
    for (i in 2:nrow(ci2_full)) {
      ti = ci2_full$t[i]
      tdiff = ti - max(ci2_full[which(ci2_full$t < ti &
                                        ci2_full$obs == 1), 1])
      ci2_full$prod[i] = tdiff * ci2_full$ci[i]
    }
    df2 = ci2_full %>% dplyr::group_by(x, y, obs) %>%
      dplyr::summarise(sum = sum(logc),
                       int = sum(prod))
    df2 = ci2_full %>% dplyr::group_by(x, y, obs) %>%
      dplyr::summarise(sum = sum(logc),
                       int = sum(c(diff(t), 0) * ci))

    df2 = dplyr::arrange(df2, obs, x, y)
    df1 = dplyr::arrange(df1, obs, x, y)

    res = df1$sum[which(df1$obs == 1)] - df1$int[which(df1$obs == 0)] -
      (df2$sum[which(df2$obs == 1)] - df2$int[which(df2$obs == 0)])
    dev = data.frame(res = res,
                     x = model1$br_grid$lon,
                     y = model1$br_grid$lat)
    out = list(dev = dev,
               residual_sum = sum(res))
    ifelse(
      sum(res) > 0,
      print("Model 1 provides a better fit."),
      print("Model 2 provides a better fit.")
    )
  } else {
    out = dr1
  }
  return(out)
}


ci_full = function(model,
                   time_seq = seq(1, ceiling(max(model$dat$times)), 1)) {
  df = model$data
  times = model$data$times
  marks = model$data$marks
  lat = model$data$lat
  lon = model$data$lon
  xn = max(model$locs$x)
  yn = max(model$locs$y)
  x1 = min(model$locs$x)
  y1 = min(model$locs$y)
  tn = ceiling(max(times))
  #ts = seq(1, tn, 1)
  br = model$br
  x_pix = model$x_pix
  y_pix = model$y_pix
  n = length(x_pix) * length(y_pix)

  # create a data frame with a row for each pixel (midpt)
  # for each day.
  # uses the background rate for each pixel as estimated
  # in the model.
  xy = data.frame(
    lon = rep(model$br_grid$lon, length(time_seq)),
    lat = rep(model$br_grid$lat, length(time_seq)),
    x = rep(rep(x1:xn, each = length(y_pix)), length(time_seq)),
    y = rep(rep(y1:yn, length(x_pix))),
    t = rep(time_seq, each = n),
    br = rep(model$br_grid$br, length(time_seq)),
    ci = NA
  )

  # so, xy is data frame of pixels
  # and df is data frame of actual observations
  # for each daily pixel, get conditional intensity
  for (i in 1:nrow(xy)) {
    df_sub = df[which(df$times < xy$t[i]), ]
    m = nrow(df_sub)

    #for (i in 1:m) {
    sim_trig = 0
    dist_vec = rep(NA, m)
    R = 6371e3

    # get distance vector for simulated point with each previous
    # observed point
    for (j in 1:m) {
      phi1 = xy$lat[i] * pi / 180
      phi2 = df_sub$lat[j] * pi / 180
      latdiff = (df_sub$lat[j] - xy$lat[i]) * pi / 180
      londiff = (df_sub$lon[j] - xy$lon[i]) * pi / 180

      a = sin(latdiff / 2) * sin(latdiff / 2) +
        cos(phi1) * cos(phi2) * sin(londiff / 2) * sin(londiff / 2)
      b = 2 * atan2(sqrt(a), sqrt(1 - a))
      d = R * b

      dist_vec[j] = d
    }

    # convert distances to correct units
    if (model$input$dist_unit == "mile") {
      dist_vec = dist_vec * 0.000621371
    } else if (model$input$dist_unit == "kilometer") {
      dist_vec = dist_vec * 0.001
    } else {
      dist_vec = dist_vec
    }

    # get difference between simulated pixel location/time and
    # actual observations.
    xg = xy$t[i] - df_sub$times
    xh = dist_vec
    xk = df_sub$marks

    # using these new spatio-temporal differences,
    # plug into estimated g(t), h(s), k(m)
    # to calculate triggering component.
    # then add to br to get conditional intensity at that location
    # up to that point in time.
    gg = predict(model$fitg, newdata = data.frame(midpt = xg))
    hh = predict(model$fith, newdata = data.frame(midpt = xh))
    kk = predict(model$fitk, newdata = data.frame(midpt = marks))[1:m]

    gg = ifelse(gg < 0, 0, gg)
    hh = ifelse(hh < 0, 0, hh)
    kk = ifelse(kk < 0, 0, kk)

    ghk = gg * hh * kk
    sim_trig = sum(ghk)
    xy$ci[i] = xy$br[i] + sim_trig
  }
  return(xy)
}

ci_full_bin = function(model,
                       time_seq = seq(1, ceiling(max(model$dat$times)), 1)) {
  df = model$data
  times = model$data$times
  marks = model$data$marks
  lat = model$data$lat
  lon = model$data$lon
  xn = max(model$locs$x)
  yn = max(model$locs$y)
  x1 = min(model$locs$x)
  y1 = min(model$locs$y)
  tn = ceiling(max(times))
  #ts = seq(1, tn, 1)
  br = model$br
  x_pix = model$x_pix
  y_pix = model$y_pix
  n = length(x_pix) * length(y_pix)

  time_breaks = model$time_breaks
  space_breaks = model$space_breaks
  mark_breaks = model$mark_breaks

  bin_f <- function(u, v) {
    x <- rep(0, length(u))
    for (j in 1:length(u)) {
      for (i in 1:(length(v) - 1)) {
        if (v[i] < u[j] & u[j] <= v[i + 1]) {
          x[j] <- i
        }
      }
    }
    return(x)
  }

  xy = data.frame(
    lon = rep(model$br_grid$lon, length(time_seq)),
    lat = rep(model$br_grid$lat, length(time_seq)),
    x = rep(rep(x1:xn, each = length(y_pix)), length(time_seq)),
    y = rep(rep(y1:yn, length(x_pix))),
    t = rep(time_seq, each = n),
    br = rep(model$br_grid$br, length(time_seq)),
    ci = NA
  )

  for (i in 1:nrow(xy)) {
    df_sub = df[which(df$times < xy$t[i]), ]
    m = nrow(df_sub)

    #for (i in 1:m) {
    sim_trig = 0
    dist_vec = rep(NA, m)
    R = 6371e3

    # get distance vector for simulated point with each previous
    # observed point
    for (j in 1:m) {
      phi1 = xy$lat[i] * pi / 180
      phi2 = df_sub$lat[j] * pi / 180
      latdiff = (df_sub$lat[j] - xy$lat[i]) * pi / 180
      londiff = (df_sub$lon[j] - xy$lon[i]) * pi / 180

      a = sin(latdiff / 2) * sin(latdiff / 2) +
        cos(phi1) * cos(phi2) * sin(londiff / 2) * sin(londiff / 2)
      b = 2 * atan2(sqrt(a), sqrt(1 - a))
      d = R * b

      dist_vec[j] = d
    }

    # convert distances to correct units
    if (model$input$dist_unit == "mile") {
      dist_vec = dist_vec * 0.000621371
    } else if (model$input$dist_unit == "kilometer") {
      dist_vec = dist_vec * 0.001
    } else {
      dist_vec = dist_vec
    }

    xg = xy$t[i] - df_sub$times
    gb = bin_f(xg, time_breaks)
    xh = dist_vec
    hb = bin_f(xh, space_breaks)
    xk = df_sub$marks
    kb = bin_f(xk, mark_breaks)

    gg = model$g[gb]
    hh = model$h[hb]
    kk = model$k[kb]

    ghk = gg * hh * kk
    sim_trig = sum(ghk)
    xy$ci[i] = xy$br[i] + sim_trig
    i = i + 1
    # if(i%%10 == 0) {
    #   print(i/nrow(xy))
    # }
  }
  return(xy)
}


# Deviance Residual Plot --------------------------------------------------

#' Deviance residual plot
#'
#' Plot displaying the deviance residuals over the spatial domain.
#'
#' @param dev1 the deviance residuals object from the first model
#' @param dev2 the deviance residuals object from the second model
#' @param plot_title the optional character string to be used for the plot title
#' @param obs_data the observational data used in fitting the models
#' @param lat_lim a vector of the minimum latitude, maximum latitude, and height of each pixel
#' @param lon_lim a vector of the minimum longitude, maximum longitude, and width of each pixel
#' @param include_pts TRUE if observational data is plotted
#' @param include_vals TRUE if the deviance residual values are plotted in each pixel
#'
#' @example
#' dp = dev_plot(dev1 = dr1,
#'               dev2 = dr2,
#'               obs_data = hm
#'               )

#' @export
dev_plot = function(dev1,
                    dev2,
                    plot_title = NULL,
                    obs_data = NULL,
                    lon_lim = seq(min(obs_data$lon) - 0.01, max(obs_data$lon) + 0.01,
                                  (max(obs_data$lon) - min(obs_data$lon)) / 10),
                    lat_lim = seq(min(obs_data$lat) - 0.01, max(obs_data$lat) + 0.01,
                                  (max(obs_data$lat) - min(obs_data$lat)) / 10),
                    include_pts = TRUE,
                    include_vals = FALSE) {
  lon_lim = lon_lim[-length(lon_lim)]
  lat_lim = lat_lim[-length(lat_lim)]
  df = data.frame(
    dr = dev1$res1 - dev2$res1,
    x = dev1$x,
    y = dev1$y,
    lon = rep(lon_lim, each = length(unique(dev2$y))),
    lat = rep(lat_lim, length(unique(dev1$x)))
  )

  df$lat = df$lat + (max(obs_data$lat) - min(obs_data$lat)) / 20
  df$lon = df$lon + (max(obs_data$lon) - min(obs_data$lon)) / 20

  dr_plot = ggplot2::ggplot() +
    ggplot2::scale_fill_gradient2(low = "red",
                                  high = "blue",
                                  name = "Deviance \nResiduals") +
    ggplot2::geom_tile(data = df, aes(x = lon, y = lat, fill = dr)) +
    ggplot2::scale_x_continuous(ggplot2::aes(x = lon, y = lat),
                                minor_breaks = NULL) +
    ggplot2::scale_y_continuous(ggplot2::aes(x = lon, y = lat),
                                minor_breaks = NULL) +
    ggplot2::xlab("lon") +
    ggplot2::ylab("lat") +
    ggplot2::ggtitle(plot_title)


  if(include_pts == TRUE) {
    dr_plot = dr_plot + ggplot2::geom_point(data = obs_data,
                                            ggplot2::aes(x = lon, y = lat),
                                            alpha = 0.25)
  }
  if(include_vals == TRUE) {
    dr_plot = dr_plot + ggplot2::geom_text(data = df,
                                           ggplot2::aes(x = lon, y = lat,
                                                        label = round(dr)),
                                           size = 3)
  }
  return(dr_plot)
}



