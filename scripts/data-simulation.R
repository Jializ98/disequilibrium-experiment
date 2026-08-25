simulate_at_cti_stochastic <- function(
		n_step = 30,
		AT0 = 20,
		CTI0 = 10,
		a = 0,
		b = -0.20,
		c = NULL,
		sd_AT = 0.5,
		sd_CTI = 0.5,
		plot_AT = 0,
		plot_CTI = 0,
		case = NULL,
		seed = NULL
) {
	stopifnot(
		length(n_step) == 1,
		n_step >= 1,
		length(AT0) == 1,
		length(CTI0) == 1,
		length(a) == 1,
		length(b) == 1,
		length(sd_AT) == 1,
		sd_AT >= 0,
		length(sd_CTI) == 1,
		sd_CTI >= 0,
		length(plot_AT) == 1,
		length(plot_CTI) == 1
	)

	constrained_c <- is.null(c)

	if (constrained_c) {
		c <- a - b
	}

	if (!is.null(seed)) {
		set.seed(seed)
	}

	out <- tibble(
		time = 0:n_step,
		AT = NA_real_,
		CTI = NA_real_,
		epsilon_AT = NA_real_,
		epsilon_CTI = NA_real_
	)

	out$AT[1] <- AT0
	out$CTI[1] <- CTI0

	for (tt in seq_len(n_step)) {
		out$epsilon_AT[tt + 1] <- rnorm(1, mean = 0, sd = sd_AT)
		out$epsilon_CTI[tt + 1] <- rnorm(1, mean = 0, sd = sd_CTI)
		out$AT[tt + 1] <- (1 + a) * out$AT[tt] +
			plot_AT + out$epsilon_AT[tt + 1]
		out$CTI[tt + 1] <- (1 + b) * out$CTI[tt] +
			c * out$AT[tt] + plot_CTI + out$epsilon_CTI[tt + 1]
	}

	if (is.null(case)) {
		case <- paste0(
			"AT0 = ", AT0,
			", CTI0 = ", CTI0,
			", a = ", a,
			", b = ", b,
			", c = ", round(c, 3),
			", stochastic"
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
			sd_AT = sd_AT,
			sd_CTI = sd_CTI,
			plot_AT = plot_AT,
			plot_CTI = plot_CTI,
			c_constrained = constrained_c
		)
}

simulate_at_cti_stochastic_plots <- function(
		n_plot = 20,
		n_step = 30,
		AT0_min = 15,
		AT0_max = 25,
		CTI0_min = 5,
		CTI0_max = 35,
		a = 0,
		b = -0.20,
		c = NULL,
		sd_AT = 0.5,
		sd_CTI = 0.5,
		sigma_plot_AT = 0.5,
		sigma_plot_CTI = 0.5,
		case = "Stochastic multi-plot simulation",
		seed = NULL
) {
	stopifnot(
		length(n_plot) == 1,
		n_plot >= 1,
		length(AT0_min) == 1,
		length(AT0_max) == 1,
		AT0_min <= AT0_max,
		length(CTI0_min) == 1,
		length(CTI0_max) == 1,
		CTI0_min <= CTI0_max
		,
		length(sigma_plot_AT) == 1,
		sigma_plot_AT >= 0,
		length(sigma_plot_CTI) == 1,
		sigma_plot_CTI >= 0
	)

	if (!is.null(seed)) {
		set.seed(seed)
	}

	plot_initial_states <- tibble(
		plot = seq_len(n_plot),
		AT0 = runif(n_plot, min = AT0_min, max = AT0_max),
		CTI0 = runif(n_plot, min = CTI0_min, max = CTI0_max),
		plot_AT = rnorm(n_plot, mean = 0, sd = sigma_plot_AT),
		plot_CTI = rnorm(n_plot, mean = 0, sd = sigma_plot_CTI)
	)

	map_dfr(
		seq_len(n_plot),
		\(plot_id) {
			initial_state <- plot_initial_states[plot_id, ]

			simulate_at_cti_stochastic(
				n_step = n_step,
				AT0 = initial_state$AT0,
				CTI0 = initial_state$CTI0,
				a = a,
				b = b,
				c = c,
				sd_AT = sd_AT,
				sd_CTI = sd_CTI,
				plot_AT = initial_state$plot_AT,
				plot_CTI = initial_state$plot_CTI,
				case = case
			) |>
			mutate(
				plot = plot_id,
				AT0 = initial_state$AT0,
				CTI0 = initial_state$CTI0
			)
		}
	)
}
