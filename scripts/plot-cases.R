generate_plot_at_cti <- function(
    n_step = 30,
    AT0 = 20,
    CTI0 = 10,
    a = 0,
    b = -0.20,
    c = NULL,
    case = NULL
) {
  
  stopifnot(
    length(n_step) == 1,
    n_step >= 1,
    length(AT0) == 1,
    length(CTI0) == 1,
    length(a) == 1,
    length(b) == 1
  )
  
  constrained_c <- is.null(c)
  
  if (constrained_c) {
    c <- a - b
  }
  
  out <- tibble(
    time = 0:n_step,
    AT = NA_real_,
    CTI = NA_real_
  )
  
  out$AT[1] <- AT0
  out$CTI[1] <- CTI0
  
  for (tt in seq_len(n_step)) {
    out$AT[tt + 1] <- (1 + a) * out$AT[tt]
    out$CTI[tt + 1] <- (1 + b) * out$CTI[tt] + c * out$AT[tt]
  }
  
  if (is.null(case)) {
    case <- paste0(
      "AT0 = ", AT0,
      ", CTI0 = ", CTI0,
      ", a = ", a,
      ", b = ", b,
      ", c = ", round(c, 3)
    )
  }
  
  out |>
    mutate(
      case = case,
      D = AT - CTI,
      delta_AT = AT - lag(AT),
      delta_CTI = CTI - lag(CTI),
      delta_D = D - lag(D),
      a = a,
      b = b,
      c = c,
      c_constrained = constrained_c
    )
}

plot_at_cti_dynamics <- function(sim, facet = FALSE) {
  
  required <- c(
    "time", "AT", "CTI", "case",
    "delta_AT", "delta_CTI"
  )
  missing <- setdiff(required, names(sim))
  
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }
  
  grouping <- c("case", "plot")
  grouping <- grouping[grouping %in% names(sim)]

  sim_space <- sim |>
    group_by(across(all_of(grouping))) |>
    arrange(time, .by_group = TRUE) |>
    mutate(
      AT_next = lead(AT),
      CTI_next = lead(CTI)
    ) |>
    ungroup()

  space_limits <- range(
    c(sim_space$AT, sim_space$CTI),
    na.rm = TRUE
  )

  space_field <- tidyr::expand_grid(
    case = unique(sim_space$case),
    AT = seq(space_limits[1], space_limits[2], length.out = 6),
    CTI = seq(space_limits[1], space_limits[2], length.out = 6)
  ) |>
    left_join(
      sim_space |>
        distinct(case, a, b, c),
      by = "case"
    ) |>
    mutate(
      AT_next = AT + a * AT,
      CTI_next = CTI + b * CTI + c * AT
    )
  
  p_space <- ggplot(sim_space, aes(x = AT, y = CTI)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    geom_segment(
      data = space_field,
      aes(
        x = AT,
        y = CTI,
        xend = AT_next,
        yend = CTI_next
      ),
      color = "grey55",
      alpha = 0.65,
      arrow = arrow(
        length = grid::unit(0.08, "cm"),
        type = "closed"
      ),
      linewidth = 0.35,
      inherit.aes = FALSE
    ) +
    geom_segment(
      aes(
        xend = AT_next,
        yend = CTI_next,
        color = time
      ),
      arrow = arrow(
        length = grid::unit(0.10, "cm"),
        type = "closed"
      ),
      linewidth = 0.55,
      na.rm = TRUE
    ) +
    geom_point(
      aes(color = time),
      size = 2.2
    ) +
    scale_color_viridis_c(name = "Time") +
    coord_fixed(
      ratio = 1,
      xlim = space_limits,
      ylim = space_limits
    ) +
    labs(
      x = "Annual temperature (x)",
      y = "Community temperature index (μ)",
      title = "AT--CTI state space",
      subtitle = "Dashed line: D = AT - CTI = 0"
    ) +
    theme(
      legend.position = "right",
      strip.text = element_text(size = 10)
    )

  sim_delta_space <- sim |>
    group_by(across(all_of(grouping))) |>
    arrange(time, .by_group = TRUE) |>
    mutate(
      delta_AT_next = lead(delta_AT),
      delta_CTI_next = lead(delta_CTI)
    ) |>
    ungroup()

  delta_limits <- range(
    c(sim_delta_space$delta_AT, sim_delta_space$delta_CTI),
    na.rm = TRUE
  )

  delta_field <- tidyr::expand_grid(
    case = unique(sim_delta_space$case),
    delta_AT = seq(delta_limits[1], delta_limits[2], length.out = 6),
    delta_CTI = seq(delta_limits[1], delta_limits[2], length.out = 6)
  ) |>
    left_join(
      sim_delta_space |>
        distinct(case, a, b, c),
      by = "case"
    ) |>
    mutate(
      delta_AT_next = delta_AT + a * delta_AT,
      delta_CTI_next = delta_CTI + b * delta_CTI + c * delta_AT
    )

  p_delta_space <- ggplot(
    sim_delta_space,
    aes(x = delta_AT, y = delta_CTI)
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    geom_segment(
      data = delta_field,
      aes(
        x = delta_AT,
        y = delta_CTI,
        xend = delta_AT_next,
        yend = delta_CTI_next
      ),
      color = "grey55",
      alpha = 0.65,
      arrow = arrow(
        length = grid::unit(0.08, "cm"),
        type = "closed"
      ),
      linewidth = 0.35,
      inherit.aes = FALSE
    ) +
    geom_segment(
      aes(
        xend = delta_AT_next,
        yend = delta_CTI_next,
        color = time
      ),
      arrow = arrow(
        length = grid::unit(0.10, "cm"),
        type = "closed"
      ),
      linewidth = 0.55,
      na.rm = TRUE
    ) +
    geom_point(
      aes(color = time),
      size = 2.2,
      na.rm = TRUE
    ) +
    scale_color_viridis_c(name = "Time") +
    coord_fixed(
      ratio = 1,
      xlim = delta_limits,
      ylim = delta_limits
    ) +
    labs(
      x = expression(Delta * AT == Delta * x),
      y = expression(Delta * CTI == Delta * mu),
      title = "First-difference state space",
      subtitle = "Trajectory of first-order changes in AT and CTI"
    ) +
    theme(
      legend.position = "right",
      strip.text = element_text(size = 10)
    )
  
  sim_long <- sim |>
    select(case, time, AT, CTI) |>
    pivot_longer(
      cols = c(AT, CTI),
      names_to = "variable",
      values_to = "value"
    ) |>
    mutate(
      variable = factor(variable, levels = c("AT", "CTI"))
    )

  time_limits <- range(
    c(sim_long$time, sim_long$value),
    na.rm = TRUE
  )
  
  p_time <- ggplot(
    sim_long,
    aes(
      x = time,
      y = value,
      color = variable,
      group = interaction(case, variable, drop = TRUE)
    )
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    coord_fixed(
      ratio = 1,
      xlim = time_limits,
      ylim = time_limits
    ) +
    labs(
      x = "Time",
      y = "Value",
      color = NULL,
      title = "Temporal dynamics of AT and CTI"
    ) +
    theme(
      legend.position = "top",
      strip.text = element_text(size = 10)
    )
  
  if (facet) {
    p_space <- p_space + facet_wrap(~case)
    p_delta_space <- p_delta_space + facet_wrap(~case)
    p_time <- p_time + facet_wrap(~case)
  }
  
  list(
    space = p_space,
    delta_space = p_delta_space,
    time_value = p_time
  )
}
