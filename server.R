shinyServer(function(input, output) {

    # =============================================================================================
    # BAYESIAN INFERENCE DEMO
    # =============================================================================================

    output$bayes <- renderPlot({
        alpha <- parse_number(input$bayes_alpha)
        beta <- parse_number(input$bayes_beta)
        trials <- parse_number(input$bayes_trials)
        successes <- parse_number(input$bayes_successes)
        est_prob <- successes / trials
        post_alpha <- alpha + successes
        post_beta <- beta + trials - successes
        x_vals <- seq(0, 1, length.out = 101)
        max_prior <- max(dbeta(x_vals, alpha, beta))
        max_posterior <- max(dbeta(x_vals, post_alpha, post_beta))
        max_y <- min(8, max(max_prior, max_posterior) + 0.1)
        max_binom <- dbinom(successes, trials, est_prob)
        mle <- function(p) {
            0.99 * max_y / max_binom * choose(trials, successes) * p^successes * (1 - p)^(trials - successes)
        }
        ggplot() +
            geom_function(fun = dbeta, args = list(shape1 = alpha, shape2 = beta), size = 1,
                          aes(color = "Prior")) +
            geom_function(fun = mle, xlim = c(0, 1), linetype = 2, size = 1,
                          aes(color = "Scaled Likelihood")) +
            geom_function(fun = dbeta, args = list(shape1 = post_alpha, shape2 = post_beta),
                          size = 1, aes(color = "Posterior")) +
            scale_color_manual(values = c("red1", "royalblue", "grey")) +
            labs(x = latex2exp::TeX("Parameter Value"), y = "Density", color = "Distribution") +
            theme(text = element_text(size = 17))
    })

    observeEvent(input$run_bayes, {
        output$bayes_spinner <- renderUI({
            withSpinner(plotOutput("bayes", height = "550px", width = "1000px"), type = 5, color = "#000000")
        })
    })

    # =============================================================================================
    # MONTE CARLO SAMPLING DEMO
    # =============================================================================================

    output$distribution_controls <- renderUI({
        if (input$choose_distribution == "Normal distribution") {
            splitLayout(
                textInput(inputId = "monte_carlo_normal_mean", label = "Mean:", value = 0),
                textInput(inputId = "monte_carlo_normal_sd", label = "Standard Deviation:", value = 1),
                cellArgs = list(style = "padding: 5px")
            )
        } else {
            splitLayout(
                textInput(inputId = "monte_carlo_cauchy_location", label = "Location:", value = 0),
                textInput(inputId = "monte_carlo_cauchy_scale", label = "Scale:", value = 1),
                cellArgs = list(style = "padding: 5px")
            )
        }
    })

    monte_carlo_samples <- eventReactive(input$run_monte_carlo, {
        set.seed(22)
        num_samples <- parse_number(input$num_samples)
        sample_size <- parse_number(input$monte_carlo_sample_size)
        distribution <- isolate(input$choose_distribution)
        if (distribution == "Normal distribution") {
            mean <- parse_number(input$monte_carlo_normal_mean)
            sd <- parse_number(input$monte_carlo_normal_sd)
            samples <- as.data.frame(replicate(num_samples, cummean(rnorm(sample_size, mean, sd))))
            colnames(samples) <- paste0("sample_", 1:num_samples)
            samples
        } else {
            location <- parse_number(input$monte_carlo_cauchy_location)
            scale <- parse_number(input$monte_carlo_cauchy_scale)
            samples <- as.data.frame(replicate(num_samples, cummean(rcauchy(sample_size, location, scale))))
            colnames(samples) <- paste0("sample_", 1:num_samples)
            samples
        }
    })

    output$monte_carlo <- renderPlot({
        distribution <- isolate(input$choose_distribution)
        if (distribution == "Normal distribution") {
            mean <- parse_number(input$monte_carlo_normal_mean)
            sd <- parse_number(input$monte_carlo_normal_sd)
            y_limits <- c(mean - 2 * sd, mean + 2 * sd)
        } else {
            location <- parse_number(input$monte_carlo_cauchy_location)
            scale <- parse_number(input$monte_carlo_cauchy_scale)
            y_limits <- c(location - 4 * scale, location + 4 * scale)
        }
        plot <- monte_carlo_samples() %>%
            mutate(n = row_number()) %>%
            pivot_longer(cols = paste0("sample_", 1:input$num_samples)) %>%
            ggplot(aes(x = n, y = value, color = name)) +
            geom_line() +
            ylim(y_limits) +
            scale_color_manual(values = wes_palette(n = input$num_samples, name = "Zissou1", type = "continuous")) +
            theme(legend.position = "none") +
            labs(x = "Observation Number", y = "Cumulative Sample Mean") +
            theme(text = element_text(size = 17))
        use_log_axis <- isolate(input$monte_carlo_log)
        if (use_log_axis) {
            plot + scale_x_log10()
        } else {
            plot
        }
    })

    observeEvent(input$run_monte_carlo, {
        output$monte_carlo_spinner <- renderUI({
            withSpinner(plotOutput("monte_carlo", height = "500px"), type = 5, color = "#000000")
        })
    })

    # =============================================================================================
    # MARKOV CHAINS DEMO
    # =============================================================================================

    iterate <- function(x, P, n) {
        result <- matrix(NA, n + 1, length(x))
        result[1, ] <- x
        for (i in seq_len(n)) {
            result[i + 1, ] <- x <- x %*% P
        }
        result
    }

    markov_chain_data <- eventReactive(input$run_markov_chains, {
        n <- 10
        P <- isolate(input$transition_matrix)
        state_1 <- as.data.frame(iterate(c(1, 0, 0), P, n))
        colnames(state_1) <- c("1", "2", "3")
        state_1 <- state_1 %>%
            mutate(n = row_number()) %>%
            pivot_longer(c("1", "2", "3"), names_to = "points", values_to = "state_1")
        state_2 <- as.data.frame(iterate(c(0, 1, 0), P, n))
        colnames(state_2) <- c("1", "2", "3")
        state_2 <- state_2 %>%
            mutate(n = row_number()) %>%
            pivot_longer(c("1", "2", "3"), names_to = "points", values_to = "state_2")
        state_3  <- as.data.frame(iterate(c(0, 0, 1), P, n))
        colnames(state_3) <- c("1", "2", "3")
        state_3 <- state_3 %>%
            mutate(n = row_number()) %>%
            pivot_longer(c("1", "2", "3"), names_to = "points", values_to = "state_3")
        state_1 %>%
            left_join(state_2, by = c("n", "points")) %>%
            left_join(state_3, by = c("n", "points")) %>%
            pivot_longer(cols = c("state_1", "state_2", "state_3"))
    })

    output$convergence <- renderPlot({
        markov_chain_data() %>%
            ggplot(aes(x = n, y = value, color = points, linetype = name)) +
            geom_line() +
            scale_linetype(guide = "none") +
            scale_x_continuous(breaks = 1:10) +
            scale_color_manual(values = wes_palette(n = 3, name = "Darjeeling1")) +
            labs(x = "Step", y = "Probability", color = "State") +
            theme(text = element_text(size = 17))
    })

    observeEvent(input$run_markov_chains, {
        output$convergence_spinner <- renderUI({
            withSpinner(plotOutput("convergence", height = "500px"), type = 5, color = "#000000")
        })
    })

    iterate_samples <- function(i, P, n) {
        result <- integer(n)
        for (t in seq_len(n))
            result[[t]] <- i <- sample(nrow(P), 1, pr = P[i, ])
        result
    }

    state_data <- eventReactive(input$run_markov_chains, {
        P <- isolate(input$transition_matrix)
        n <- isolate(input$steps)
        states <- iterate_samples(1, P, n)
        data.frame(state = states) %>%
            mutate(n = row_number())
    })

    output$state <- renderPlot({
         state_data() %>%
            ggplot(aes(x = n, y = state)) +
            geom_line() +
            geom_point(size = 1) +
            xlim(0, 100) +
            labs(x = "Step", y = "State") +
            scale_y_continuous(breaks = 1:3) +
            theme(text = element_text(size = 17))
    })

    observeEvent(input$run_markov_chains, {
        output$state_spinner <- renderUI({
            withSpinner(plotOutput("state", height = "500px"), type = 5, color = "#000000")
        })
    })

    output$time <- renderPlot({
        state_1 <- cummean(state_data()$state == 1)
        state_2 <- cummean(state_data()$state == 2)
        state_3 <- cummean(state_data()$state == 3)
        n <- length(state_1)
        data.frame(state_1, state_2, state_3) %>%
            mutate(n = row_number()) %>%
            pivot_longer(cols = c("state_1", "state_2", "state_3"), names_to = "state", values_to = "prob") %>%
            ggplot(aes(x = n, y = prob, color = state)) +
            geom_hline(yintercept = state_1[n], color =  wes_palette(n = 3, name = "Darjeeling1")[1], linetype = 2) +
            geom_hline(yintercept = state_2[n], color = wes_palette(n = 3, name = "Darjeeling1")[2], linetype = 2) +
            geom_hline(yintercept = state_3[n], color = wes_palette(n = 3, name = "Darjeeling1")[3], linetype = 2) +
            geom_line() +
            scale_x_log10() +
            scale_color_manual(values = wes_palette(n = 3, name = "Darjeeling1"), labels = c("1", "2", "3")) +
            labs(x = "Step", y = "Cumulative Time Spent in State", color = "State") +
            theme(text = element_text(size = 17))
    })

    observeEvent(input$run_markov_chains, {
        output$time_spinner <- renderUI({
            withSpinner(plotOutput("time", height = "500px"), type = 5, color = "#000000")
        })
    })

    # =============================================================================================
    # MARKOV CHAIN MONTE CARLO DEMO (ONE-DIMENSIONAL: GAUSSIAN MIXTURE MODEL)
    # =============================================================================================

    mcmc_sample <- eventReactive(input$run_mcmc, {
        set.seed(22)
        x_init <- parse_number(input$x_init)
        proposal_sd <- parse_number(input$proposal_sd)
        proposal_dist <- function(x) gaussian_proposal_1d(x, sd = proposal_sd)
        max_iter <- parse_number(input$max_iter)
        metropolis_hastings(x_init, gmm_posterior, proposal_dist, max_iter) %>%
            as.data.frame() %>%
            mutate(step = row_number())
    })

    output$histogram <- renderPlot({
        ggplot(mcmc_sample()) +
            geom_histogram(aes(x = theta, y = ..density..), color = "black", fill = "transparent", binwidth = 0.1) +
            geom_function(fun = gmm_posterior, color = "red") +
            xlim(c(-4, 8)) +
            labs(x = "Parameter Value", y = "Posterior Density") +
            theme(text = element_text(size = 17))
    })

    output$histogram_spinner <- renderUI({
        if (input$run_mcmc != 0) {
            withSpinner(plotOutput("histogram", height = "500px"), type = 5, color = "#000000")
        } else {
            column(
                width = 12,
                div(
                    class = "text-center",
                    h3("Click Run to draw an MCMC sample..."),
                    br(),
                    img(src = "dice.gif", width = "50%", align = "center")
                )
            )
        }
    })

    output$markov_chain <- renderPlot({
        p <- mcmc_sample() %>%
            filter(step <= 1000) %>%
            ggplot(aes(x = step, y = theta)) +
            geom_line(color = "black") +
            labs(x = "Step", y = "Sample Value") +
            theme(text = element_text(size = 17))
        x_limits <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
        theoretical_density <- ggplot() +
            geom_function(fun = gmm_posterior) +
            coord_flip() +
            lims(x = x_limits) +
            ylab("Theoretical Distribution") +
            theme(
                text = element_text(size = 17),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_blank()
            )
        grid.arrange(p, theoretical_density, ncol = 2, nrow  = 1, widths = c(4, 1))
    })

    output$markov_chain_spinner <- renderUI({
        if (input$run_mcmc != 0) {
            withSpinner(plotOutput("markov_chain", height = "500px"), type = 5, color = "#000000")
        } else {
            column(
                width = 12,
                div(
                    class = "text-center",
                    h3("Click Run to draw an MCMC sample..."),
                    br(),
                    img(src = "dice.gif", width = "50%", align = "center")
                )
            )
        }
    })

    output$decision <- renderPlotly({
        plot <- ggplot(mcmc_sample(), aes(x = step, y = proposed_theta, color = decision)) +
            geom_point(shape = 18, size = 0.5,
                       aes(text = sprintf("Step: %s<br>Proposed Value: %s<br>Acceptance Probability: %s<br>Decision: %s",
                                          step, round(proposed_theta, 3), round(acceptance_prob, 3), decision))) +
            scale_color_manual(labels = c("Accepted", "Rejected"), values = c("#049660", "#bf0000")) +
            labs(x = "Step", y = "Proposed Value", color = "Decision:") +
            theme(legend.position = "none")
        ggplotly(plot, tooltip = "text")
    })

    output$decision_spinner <- renderUI({
        if (input$run_mcmc != 0) {
            withSpinner(plotlyOutput("decision", height = "500px"), type = 5, color = "#000000")
        } else {
            column(
                width = 12,
                div(
                    class = "text-center",
                    h3("Click Run to draw an MCMC sample..."),
                    br(),
                    img(src = "dice.gif", width = "50%", align = "center")
                )
            )
        }
    })

    output$acceptance_rate <- renderPlot({
        mcmc_sample() %>%
            mutate(decision_numeric = ifelse(decision == "accepted", 1, 0)) %>%
            ggplot(aes(x = step, y = cummean(decision_numeric))) +
            geom_hline(aes(yintercept = mean(decision_numeric)), color = "red", linetype = 2) +
            geom_line(color = "black") +
            labs(x = "Step", y = "Cumulative Acceptance Rate") +
            theme(text = element_text(size = 17))
    })

    output$acceptance_rate_spinner <- renderUI({
        if (input$run_mcmc != 0) {
            withSpinner(plotOutput("acceptance_rate", height = "500px"), type = 5, color = "#000000")
        } else {
            column(
                width = 12,
                div(
                    class = "text-center",
                    h3("Click Run to draw an MCMC sample..."),
                    br(),
                    img(src = "dice.gif", width = "50%")
                )
            )
        }
    })

    # =============================================================================================
    # MARKOV CHAIN MONTE CARLO DEMO (TWO-DIMENSIONAL: SUNSPOTS AND THE GAMMA DISTRIBUTION)
    # =============================================================================================

    sample_2d <- eventReactive(input$run_2d, {
        set.seed(22)
        theta_init <- c(parse_number(input$x_init_2d_alpha), parse_number(input$x_init_2d_beta))
        proposal_sd <- c(parse_number(input$proposal_sd_2d_first), parse_number(input$proposal_sd_2d_second))
        proposal_dist <- function(x) gaussian_proposal_2d(x, sd = proposal_sd)
        max_iter <- parse_number(input$max_iter_2d)
        mean_sunspots <- subset(sunspots$mean, sunspots$mean > 0)
        results <- metropolis_hastings(theta_init, log_lik = gamma_log_lik, prior = unif_prior,
                                       proposal = proposal_dist, max_iter = max_iter, data = mean_sunspots)
        data.frame(
            x = results$theta[, 1],
            y = results$theta[, 2],
            x_prop = results$proposed_theta[, 1],
            y_prop = results$proposed_theta[, 2],
            decision = results$decision
        )
    })

    output$estimated_distribution <- renderPlot({
        sample_2d() %>%
            ggplot() +
            geom_line(aes(x = x, y = y), size = 0.8, color = "black") +
            geom_point(aes(x = x_prop, y = y_prop, color = decision), size = 2, shape = 18, alpha = 0.75) +
            scale_color_manual(labels = c("Accepted", "Rejected"), values = c("#049660", "#bf0000")) +
            labs(x = latex2exp::TeX("$\\alpha$ Parameter Value"),
                 y = latex2exp::TeX("$\\beta$ Parameter Value"),
                 color = "Decision") +
            theme(text = element_text(size = 17))
    })

    observeEvent(input$run_2d, {
        output$estimated_distribution_spinner <- renderUI({
            withSpinner(plotOutput("estimated_distribution", height = "500px"), type = 5, color = "#000000")
        })
    })
})
