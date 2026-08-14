simulate_at_cti <- function(
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
  
  required <- c("time", "AT", "CTI", "case")
  missing <- setdiff(required, names(sim))
  
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }
  
  sim_space <- sim |>
    group_by(case) |>
    arrange(time, .by_group = TRUE) |>
    mutate(
      AT_next = lead(AT),
      CTI_next = lead(CTI)
    ) |>
    ungroup()
  
  p_space <- ggplot(sim_space, aes(x = AT, y = CTI)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.7
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
    # coord_equal() +
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
  
  p_time <- ggplot(
    sim_long,
    aes(
      x = time,
      y = value,
      color = variable,
      group = variable
    )
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
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
    p_space <- p_space + facet_wrap(~case, scales = "free")
    p_time <- p_time + facet_wrap(~case, scales = "free_y")
  }
  
  list(
    space = p_space,
    time_value = p_time,
    combined = p_space | p_time
  )
}