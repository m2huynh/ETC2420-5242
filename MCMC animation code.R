# app.R
# install.packages(c("shiny", "ggplot2", "cowplot", "ggtext"))

library(shiny)
library(ggplot2)
library(cowplot)
library(ggtext)

# ============================================================
# Helper function: log posterior kernel
# ============================================================

log_posterior_kernel <- function(theta, a_prior, b_prior, n_trials, y_successes) {
  if (theta <= 0 || theta >= 1) return(-Inf)
  
  dbeta(theta, a_prior, b_prior, log = TRUE) +
    dbinom(y_successes, size = n_trials, prob = theta, log = TRUE)
}

# ============================================================
# Main plotting function
# ============================================================

make_mcmc_plot <- function(state, input) {
  
  theta <- state$theta
  theta_new <- state$theta_new
  r_ratio <- state$r_ratio
  alpha <- state$alpha
  u_draw <- state$u_draw
  accepted <- state$accepted
  
  t <- length(theta)
  
  old <- if (t >= 2) theta[t - 1] else theta[1]
  current_theta <- theta[t]
  
  proposal <- if (t >= 2) theta_new[t] else NA
  r <- if (t >= 2) r_ratio[t] else NA
  acc_prob <- if (t >= 2) alpha[t] else NA
  u <- if (t >= 2) u_draw[t] else NA
  is_accepted <- if (t >= 2) accepted[t] else NA
  
  trace_df <- data.frame(
    iter = seq_along(theta),
    theta = theta
  )
  
  # ----------------------------------------------------------
  # Left density panel: sideways histogram
  # ----------------------------------------------------------
  
  breaks <- seq(0, 1, by = 0.025)
  hist_obj <- hist(theta, breaks = breaks, plot = FALSE)
  
  hist_df <- data.frame(
    theta_mid = hist_obj$mids,
    count = hist_obj$counts
  )
  
  if (max(hist_df$count) == 0) {
    hist_df$width <- 0
  } else {
    hist_df$width <- hist_df$count / max(hist_df$count) * 0.85
  }
  
  density_panel <- ggplot() +
    geom_rect(
      data = hist_df,
      aes(
        xmin = 0.95 - width,
        xmax = 0.95,
        ymin = theta_mid - 0.0125,
        ymax = theta_mid + 0.0125
      ),
      fill = "grey55",
      colour = NA
    ) +
    geom_point(
      aes(x = 0.95, y = current_theta),
      size = 4,
      colour = "darkgreen"
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1),
      labels = function(x) sub("^0", "", sprintf("%.1f", x)),
      position = "right"
    ) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
      title = "Density",
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 15) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      panel.border = element_rect(
        colour = "grey35",
        fill = NA,
        linewidth = 0.7
      )
    )
  
  # ----------------------------------------------------------
  # Trace panel plus proposal density
  # ----------------------------------------------------------
  
  base_trace <- ggplot() +
    geom_line(
      data = trace_df,
      aes(x = iter, y = theta),
      colour = "grey45",
      linewidth = 0.45
    )
  
  if (t >= 2) {
    y_grid <- seq(0, 1, length.out = 300)
    proposal_density <- dnorm(y_grid, mean = old, sd = input$sigma)
    proposal_density_scaled <- proposal_density / max(proposal_density) * 12
    
    proposal_df <- data.frame(
      x = t + proposal_density_scaled,
      theta = y_grid
    )
    
    trace_panel <- base_trace +
      geom_segment(
        aes(x = t, xend = t, y = 0, yend = 1),
        colour = "black",
        linewidth = 0.8
      ) +
      geom_path(
        data = proposal_df,
        aes(x = x, y = theta),
        colour = "black",
        linewidth = 0.9
      ) +
      geom_segment(
        aes(x = t, xend = t + 12, y = old, yend = old),
        colour = "blue",
        linewidth = 1.2
      ) +
      geom_point(
        aes(x = t, y = proposal),
        colour = ifelse(is_accepted, "darkgreen", "red"),
        size = 4
      ) +
      annotate(
        "text",
        x = t + 13.5,
        y = old,
        label = "N(theta[t-1], sigma)",
        hjust = 0,
        size = 5
      )
  } else {
    trace_panel <- base_trace +
      geom_point(
        aes(x = 1, y = theta[1]),
        colour = "darkgreen",
        size = 4
      )
  }
  
  trace_panel <- trace_panel +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1),
      labels = function(x) sub("^0", "", sprintf("%.1f", x))
    ) +
    scale_x_continuous(
      limits = c(1, max(60, t + 32))
    ) +
    labs(
      title = "MCMC Iteration",
      x = NULL,
      y = expression(theta)
    ) +
    theme_classic(base_size = 15) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      panel.border = element_rect(
        colour = "grey35",
        fill = NA,
        linewidth = 0.7
      )
    )
  
  # ----------------------------------------------------------
  # Formula panel
  # ----------------------------------------------------------
  
  formula_theme <- theme_void() +
    theme(
      plot.margin = margin(8, 10, 8, 10)
    )
  
  if (t < 2) {
    
    formula_panel <- ggplot() +
      geom_richtext(
        aes(
          x = 0.04,
          y = 0.75,
          label = paste0(
            "<span style='font-size:24pt;'>",
            "Initial value:&nbsp;&nbsp;",
            "&theta;<sub>0</sub> = ",
            "<span style='color:blue;'>",
            sprintf("%.3f", theta[1]),
            "</span>",
            "</span>"
          )
        ),
        hjust = 0,
        vjust = 1,
        fill = NA,
        label.color = NA
      ) +
      geom_richtext(
        aes(
          x = 0.04,
          y = 0.45,
          label = paste0(
            "<span style='font-size:20pt;'>",
            "Click <b>Take one sample</b> or <b>Play simulations</b> ",
            "to draw from Normal(&theta;<sub>t-1</sub>, &sigma;).",
            "</span>"
          )
        ),
        hjust = 0,
        vjust = 1,
        fill = NA,
        label.color = NA
      ) +
      xlim(0, 1) +
      ylim(0, 1) +
      formula_theme
    
  } else {
    
    old_col <- "blue"
    new_col <- if (is_accepted) "darkgreen" else "red"
    
    decision_then <- paste0(
      "Then&nbsp;&nbsp;&nbsp;&nbsp;",
      "&theta;<sub>t</sub> = ",
      "&theta;<sub>new</sub> = ",
      "<span style='color:", new_col, ";'>",
      sprintf("%.3f", proposal),
      "</span>"
    )
    
    decision_otherwise <- paste0(
      "Otherwise&nbsp;&nbsp;",
      "&theta;<sub>t</sub> = ",
      "&theta;<sub>t-1</sub> = ",
      "<span style='color:", old_col, ";'>",
      sprintf("%.3f", old),
      "</span>"
    )
    
    if (is_accepted) {
      decision_then <- paste0("<b>", decision_then, "</b>")
      decision_otherwise <- paste0(
        "<span style='color:grey55;'>",
        decision_otherwise,
        "</span>"
      )
    } else {
      decision_then <- paste0(
        "<span style='color:grey55;'>",
        decision_then,
        "</span>"
      )
      decision_otherwise <- paste0("<b>", decision_otherwise, "</b>")
    }
    
    step1_left <- paste0(
      "<span style='font-size:20pt;'>",
      "<b>Step 1:</b>&nbsp;&nbsp;&nbsp;",
      "r(&theta;<sub>new</sub>, &theta;<sub>t-1</sub>)",
      "&nbsp;&nbsp;=&nbsp;&nbsp;",
      "</span>"
    )
    
    step1_frac1 <- paste0(
      "<span style='font-size:20pt;'>",
      "<span style='display:inline-block; text-align:center;'>",
      "Posterior(&theta;<sub>new</sub>)",
      "<br>",
      "<span style='border-top:2px solid black;'>",
      "Posterior(&theta;<sub>t-1</sub>)",
      "</span>",
      "</span>",
      "</span>"
    )
    
    step1_frac2 <- paste0(
      "<span style='font-size:20pt;'>",
      "<span style='display:inline-block; text-align:center;'>",
      "Beta(", input$a_prior, ",", input$b_prior, ",",
      "<span style='color:", new_col, ";'>",
      sprintf("%.3f", proposal),
      "</span>)",
      " &times; Binomial(", input$n_trials, ",", input$y_successes, ",",
      "<span style='color:", new_col, ";'>",
      sprintf("%.3f", proposal),
      "</span>)",
      "<br>",
      "<span style='border-top:2px solid black;'>",
      "Beta(", input$a_prior, ",", input$b_prior, ",",
      "<span style='color:", old_col, ";'>",
      sprintf("%.3f", old),
      "</span>)",
      " &times; Binomial(", input$n_trials, ",", input$y_successes, ",",
      "<span style='color:", old_col, ";'>",
      sprintf("%.3f", old),
      "</span>)",
      "</span>",
      "</span>",
      "</span>"
    )
    
    step2 <- paste0(
      "<span style='font-size:20pt;'>",
      "<b>Step 2:</b>&nbsp;&nbsp;&nbsp;",
      "Acceptance probability&nbsp;&nbsp;",
      "&alpha;(&theta;<sub>new</sub>, &theta;<sub>t-1</sub>)",
      "&nbsp;&nbsp;=&nbsp;&nbsp;",
      "min{r(&theta;<sub>new</sub>, &theta;<sub>t-1</sub>), 1}",
      "&nbsp;&nbsp;=&nbsp;&nbsp;",
      "min{",
      sprintf("%.3f", r),
      ", 1}",
      "&nbsp;&nbsp;=&nbsp;&nbsp;",
      sprintf("%.3f", acc_prob),
      "</span>"
    )
    
    step3 <- paste0(
      "<span style='font-size:20pt;'>",
      "<b>Step 3:</b>&nbsp;&nbsp;&nbsp;",
      "Draw u ~ Uniform(0,1)",
      "&nbsp;&nbsp;=&nbsp;&nbsp;",
      sprintf("%.3f", u),
      "</span>"
    )
    
    step4_left <- paste0(
      "<span style='font-size:20pt;'>",
      "<b>Step 4:</b>&nbsp;&nbsp;&nbsp;",
      "If&nbsp;&nbsp; u &lt; ",
      "&alpha;(&theta;<sub>new</sub>, &theta;<sub>t-1</sub>)",
      "&nbsp;&nbsp;&rarr;&nbsp;&nbsp;",
      "If&nbsp;&nbsp;",
      sprintf("%.3f", u),
      " &lt; ",
      sprintf("%.3f", acc_prob),
      "</span>"
    )
    
    formula_panel <- ggplot() +
      geom_richtext(
        aes(x = 0.02, y = 0.88, label = step1_left),
        hjust = 0,
        vjust = 0.5,
        fill = NA,
        label.color = NA
      ) +
      geom_richtext(
        aes(x = 0.31, y = 0.88, label = step1_frac1),
        hjust = 0.5,
        vjust = 0.5,
        fill = NA,
        label.color = NA,
        lineheight = 1.15
      ) +
      geom_richtext(
        aes(
          x = 0.47,
          y = 0.88,
          label = "<span style='font-size:20pt;'>&nbsp;=&nbsp;</span>"
        ),
        hjust = 0.5,
        vjust = 0.5,
        fill = NA,
        label.color = NA
      ) +
      geom_richtext(
        aes(x = 0.67, y = 0.88, label = step1_frac2),
        hjust = 0.5,
        vjust = 0.5,
        fill = NA,
        label.color = NA,
        lineheight = 1.15
      ) +
      geom_richtext(
        aes(
          x = 0.93,
          y = 0.88,
          label = paste0(
            "<span style='font-size:20pt;'>=&nbsp;&nbsp;",
            sprintf("%.3f", r),
            "</span>"
          )
        ),
        hjust = 0,
        vjust = 0.5,
        fill = NA,
        label.color = NA
      ) +
      geom_richtext(
        aes(x = 0.02, y = 0.62, label = step2),
        hjust = 0,
        vjust = 0.5,
        fill = NA,
        label.color = NA
      ) +
      geom_richtext(
        aes(x = 0.02, y = 0.40, label = step3),
        hjust = 0,
        vjust = 0.5,
        fill = NA,
        label.color = NA
      ) +
      geom_richtext(
        aes(x = 0.02, y = 0.18, label = step4_left),
        hjust = 0,
        vjust = 0.5,
        fill = NA,
        label.color = NA
      ) +
      geom_richtext(
        aes(x = 0.68, y = 0.23, label = decision_then),
        hjust = 0,
        vjust = 0.5,
        fill = NA,
        label.color = NA
      ) +
      geom_richtext(
        aes(x = 0.68, y = 0.08, label = decision_otherwise),
        hjust = 0,
        vjust = 0.5,
        fill = NA,
        label.color = NA
      ) +
      xlim(0, 1) +
      ylim(0, 1) +
      formula_theme
  }
  
  # ----------------------------------------------------------
  # Combine panels
  # ----------------------------------------------------------
  
  top <- plot_grid(
    density_panel,
    trace_panel,
    nrow = 1,
    rel_widths = c(1, 3.2),
    align = "h"
  )
  
  plot_grid(
    top,
    formula_panel,
    ncol = 1,
    rel_heights = c(1.0, 1.35)
  )
}

# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  titlePanel("Metropolis-Hastings MCMC simulator"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Model settings"),
      
      numericInput(
        "n_trials",
        "Number of trials",
        value = 10,
        min = 1,
        step = 1
      ),
      
      numericInput(
        "y_successes",
        "Number of successes",
        value = 4,
        min = 0,
        step = 1
      ),
      
      numericInput(
        "a_prior",
        "Beta prior alpha",
        value = 1,
        min = 0.1,
        step = 0.1
      ),
      
      numericInput(
        "b_prior",
        "Beta prior beta",
        value = 1,
        min = 0.1,
        step = 0.1
      ),
      
      numericInput(
        "sigma",
        "Proposal sigma",
        value = 0.08,
        min = 0.001,
        step = 0.01
      ),
      
      hr(),
      
      h4("Simulation controls"),
      
      actionButton("one_step", "Take one sample"),
      
      br(),
      br(),
      
      numericInput(
        "n_sims",
        "Number of simulations to play",
        value = 1000,
        min = 1,
        step = 1
      ),
      
      numericInput(
        "steps_per_frame",
        "MCMC steps per frame",
        value = 10,
        min = 1,
        max = 100,
        step = 1
      ),
      
      numericInput(
        "frame_delay",
        "Frame delay, milliseconds",
        value = 75,
        min = 20,
        max = 1000,
        step = 25
      ),
      
      actionButton("play", "Play simulations"),
      actionButton("stop", "Stop"),
      
      br(),
      br(),
      
      actionButton("reset", "Reset"),
      
      hr(),
      
      helpText(
        "For teaching mode, use steps per frame = 1 and delay = 250.",
        "For fast mode, use steps per frame = 10 or 50 and delay = 20 to 75."
      ),
      
      verbatimTextOutput("summary")
    ),
    
    mainPanel(
      plotOutput("mcmc_plot", height = "900px")
    )
  )
)

# ============================================================
# Server
# ============================================================

server <- function(input, output, session) {
  
  state <- reactiveValues(
    theta = c(0.5),
    theta_new = c(NA),
    r_ratio = c(NA),
    alpha = c(NA),
    u_draw = c(NA),
    accepted = c(NA),
    playing = FALSE,
    remaining = 0
  )
  
  reset_chain <- function() {
    state$theta <- c(0.5)
    state$theta_new <- c(NA)
    state$r_ratio <- c(NA)
    state$alpha <- c(NA)
    state$u_draw <- c(NA)
    state$accepted <- c(NA)
    state$playing <- FALSE
    state$remaining <- 0
  }
  
  take_one_step <- function() {
    old <- tail(state$theta, 1)
    
    proposal <- rnorm(1, mean = old, sd = input$sigma)
    
    log_r <- log_posterior_kernel(
      theta = proposal,
      a_prior = input$a_prior,
      b_prior = input$b_prior,
      n_trials = input$n_trials,
      y_successes = input$y_successes
    ) -
      log_posterior_kernel(
        theta = old,
        a_prior = input$a_prior,
        b_prior = input$b_prior,
        n_trials = input$n_trials,
        y_successes = input$y_successes
      )
    
    r <- exp(log_r)
    acc_prob <- min(r, 1)
    u <- runif(1)
    
    if (u < acc_prob) {
      next_theta <- proposal
      is_accepted <- TRUE
    } else {
      next_theta <- old
      is_accepted <- FALSE
    }
    
    state$theta <- c(state$theta, next_theta)
    state$theta_new <- c(state$theta_new, proposal)
    state$r_ratio <- c(state$r_ratio, r)
    state$alpha <- c(state$alpha, acc_prob)
    state$u_draw <- c(state$u_draw, u)
    state$accepted <- c(state$accepted, is_accepted)
  }
  
  observeEvent(input$one_step, {
    state$playing <- FALSE
    state$remaining <- 0
    take_one_step()
  })
  
  observeEvent(input$play, {
    state$playing <- TRUE
    state$remaining <- input$n_sims
  })
  
  observeEvent(input$stop, {
    state$playing <- FALSE
    state$remaining <- 0
  })
  
  observeEvent(input$reset, {
    reset_chain()
  })
  
  observeEvent(input$n_trials, {
    if (input$y_successes > input$n_trials) {
      updateNumericInput(
        session,
        "y_successes",
        value = input$n_trials
      )
    }
  })
  
  observeEvent(
    list(
      input$n_trials,
      input$y_successes,
      input$a_prior,
      input$b_prior,
      input$sigma
    ),
    {
      reset_chain()
    },
    ignoreInit = TRUE
  )
  
  # ----------------------------------------------------------
  # Faster autoplay with dynamic delay
  # ----------------------------------------------------------
  # This uses invalidateLater(), not reactiveTimer(), because
  # reactiveTimer() requires a fixed numeric interval.
  # input$frame_delay can safely change while the app is running.
  
  observe({
    invalidateLater(input$frame_delay, session)
    
    isolate({
      if (isTRUE(state$playing) && state$remaining > 0) {
        
        n_to_take <- min(input$steps_per_frame, state$remaining)
        
        for (i in seq_len(n_to_take)) {
          take_one_step()
        }
        
        state$remaining <- state$remaining - n_to_take
        
        if (state$remaining <= 0) {
          state$playing <- FALSE
        }
      }
    })
  })
  
  output$mcmc_plot <- renderPlot({
    make_mcmc_plot(state, input)
  })
  
  output$summary <- renderText({
    n_iter <- length(state$theta) - 1
    n_accept <- sum(state$accepted, na.rm = TRUE)
    
    accept_rate <- if (n_iter > 0) {
      n_accept / n_iter
    } else {
      NA
    }
    
    paste0(
      "Iterations: ", n_iter, "\n",
      "Current theta: ", sprintf("%.3f", tail(state$theta, 1)), "\n",
      "Accepted: ", n_accept, "\n",
      "Acceptance rate: ",
      ifelse(is.na(accept_rate), "NA", sprintf("%.3f", accept_rate)), "\n",
      "Remaining in play queue: ", state$remaining, "\n",
      "Steps per frame: ", input$steps_per_frame, "\n",
      "Frame delay: ", input$frame_delay, " ms"
    )
  })
}

shinyApp(ui, server)