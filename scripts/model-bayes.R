fit_at_cti_bayes <- function(
    sim,
    niter = 5000,
    nburnin = 2000,
    thin = 4,
    nchains = 2,
    seed = 123
) {
  required <- c("plot", "time", "AT", "CTI")
  missing <- setdiff(required, names(sim))

  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }

  if (niter <= nburnin) {
    stop("niter must be larger than nburnin.")
  }

  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop(
      "Package 'nimble' is required. Install it with install.packages('nimble')."
    )
  }

  if (!"package:nimble" %in% search()) {
    library(nimble)
  }

  chain_seeds <- if (length(seed) == 1) {
    seed + seq_len(nchains) - 1
  } else {
    seed
  }

  if (length(chain_seeds) != nchains) {
    stop("seed must have length 1 or length nchains.")
  }

  sim_fit <- sim |>
    dplyr::group_by(plot) |>
    dplyr::arrange(time, .by_group = TRUE) |>
    dplyr::mutate(
      AT_next = dplyr::lead(AT),
      CTI_next = dplyr::lead(CTI)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      is.finite(AT),
      is.finite(CTI),
      is.finite(AT_next),
      is.finite(CTI_next)
    ) |>
    dplyr::mutate(plot_id = match(plot, unique(plot)))

  if (nrow(sim_fit) == 0) {
    stop("No complete transitions are available for fitting.")
  }

  n_obs <- nrow(sim_fit)
  n_plot <- dplyr::n_distinct(sim_fit$plot_id)

  model_code <- nimble::nimbleCode({
    for (i in 1:n_obs) {
      AT_next[i] ~ dnorm(mu_AT[i], tau = tau_AT)
      mu_AT[i] <- (1 + a) * AT[i] + plot_AT[plot_id[i]]

      CTI_next[i] ~ dnorm(mu_CTI[i], tau = tau_CTI)
      mu_CTI[i] <- (1 + b) * CTI[i] + c * AT[i] + plot_CTI[plot_id[i]]
    }

    for (j in 1:n_plot) {
      plot_AT[j] ~ dnorm(0, tau = tau_plot_AT)
      plot_CTI[j] ~ dnorm(0, tau = tau_plot_CTI)
    }

    a ~ dnorm(0, tau = 1 / (0.1 ^ 2))
    b ~ dnorm(0, tau = 1 / (0.5 ^ 2))
    c ~ dnorm(0, tau = 1 / (1 ^ 2))

    sigma_AT ~ dunif(0, 10)
    sigma_CTI ~ dunif(0, 10)
    sigma_plot_AT ~ dunif(0, 10)
    sigma_plot_CTI ~ dunif(0, 10)

    tau_AT <- 1 / (sigma_AT ^ 2)
    tau_CTI <- 1 / (sigma_CTI ^ 2)
    tau_plot_AT <- 1 / (sigma_plot_AT ^ 2)
    tau_plot_CTI <- 1 / (sigma_plot_CTI ^ 2)
  })

  constants <- list(
    n_obs = n_obs,
    n_plot = n_plot,
    AT = sim_fit$AT,
    CTI = sim_fit$CTI,
    plot_id = sim_fit$plot_id
  )

  data <- list(
    AT_next = sim_fit$AT_next,
    CTI_next = sim_fit$CTI_next
  )

  make_inits <- function() {
    list(
      a = 0,
      b = 0,
      c = 0,
      sigma_AT = 0.5,
      sigma_CTI = 0.5,
      sigma_plot_AT = 0.5,
      sigma_plot_CTI = 0.5,
      plot_AT = rep(0, n_plot),
      plot_CTI = rep(0, n_plot)
    )
  }

  model <- nimble::nimbleModel(
    code = model_code,
    constants = constants,
    data = data,
    inits = make_inits()
  )
  compiled_model <- nimble::compileNimble(model)

  configure <- nimble::configureMCMC(
    model,
    monitors = c(
      "a", "b", "c",
      "sigma_AT", "sigma_CTI",
      "sigma_plot_AT", "sigma_plot_CTI"
    )
  )
  mcmc <- nimble::buildMCMC(configure)
  compiled_mcmc <- nimble::compileNimble(mcmc, project = compiled_model)

  samples <- nimble::runMCMC(
    compiled_mcmc,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    nchains = nchains,
    inits = lapply(seq_len(nchains), function(...) make_inits()),
    setSeed = chain_seeds,
    samplesAsCodaMCMC = FALSE,
    summary = FALSE,
    progressBar = FALSE,
    samples = TRUE
  )

  samples_matrix <- if (is.list(samples)) {
    do.call(rbind, samples)
  } else {
    samples
  }

  posterior_samples <- as.data.frame(samples_matrix)[, c("a", "b", "c")] |>
    mutate(q = c - (a - b))

  posterior_summary <- tibble::tibble(
    parameter = c("a", "b", "c", "q"),
    estimate = vapply(posterior_samples, median, numeric(1)),
    conf_low = vapply(
      posterior_samples,
      stats::quantile,
      numeric(1),
      probs = 0.025
    ),
    conf_high = vapply(
      posterior_samples,
      stats::quantile,
      numeric(1),
      probs = 0.975
    )
  )

  posterior_plot_data <- posterior_samples |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(
      cols = c(a, b, c, q),
      names_to = "parameter",
      values_to = "value"
    )

  posterior_plot <- ggplot2::ggplot(
    posterior_plot_data,
    ggplot2::aes(x = value)
  ) +
    ggplot2::geom_density(fill = "grey75", color = "grey25") +
    ggplot2::geom_vline(
      data = posterior_summary,
      ggplot2::aes(xintercept = estimate),
      color = "firebrick",
      linetype = "dashed",
      linewidth = 0.7
    ) +
    ggplot2::facet_wrap(~parameter, scales = "free") +
    ggplot2::labs(
      x = "Posterior value",
      y = "Density",
      title = "Posterior distributions of dynamic parameters",
      subtitle = "Dashed line: posterior median"
    ) +
    ggplot2::theme_classic()

  list(
    model = model,
    mcmc = mcmc,
    samples = samples_matrix,
    posterior_summary = posterior_summary,
    posterior_plot = posterior_plot,
    data = sim_fit
  )
}
