library(dplyr)
library(ggplot2)
library(gridExtra)
library(plotly)
library(readr)
library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(shinyMatrix)
library(tidyr)
library(wesanderson)

theme_set(theme_minimal())

source("metropolis-hastings.R")

if (!file.exists("data/sunspots.csv")) {
  download.file(url = "http://www.sidc.be/silso/INFO/snmtotcsv.php", destfile = "data/sunspots.csv")
}

.variables <- c("year", "month", "year_frac", "mean", "sd", "n_obs")
sunspots <- read_delim("data/sunspots.csv", delim = ";", trim_ws = TRUE, na = c("", "-1.0", "-1"),
                       col_names = .variables, col_types = "iidddi-")

# =================================================================================================
# WELCOME PAGE
# =================================================================================================

welcome_page <- div(
  class = "welcome-panel",
  wellPanel(
    h1("Markov Chain Monte Carlo"),
    h4("Shu Amano, Maria-Cristiana Gîrjău, Leah Johnson"),
    br(),
    img(class = "center", src = "cover.png", height = "70%", width = "70%")
  )
)

# =================================================================================================
# BAYESIAN INFERENCE PAGE
# =================================================================================================

bayes_page <- tabItem(
  br(),
  tabsetPanel(
    tabPanel(
      title = "What is Bayesian Inference?",
      wellPanel(
        class = "slide",
        h1("Bayesian Inference"),
        h3("Two main interpretations of probability"),
        tags$ul(
          tags$li(HTML("<strong>Frequentist:</strong> probability is interpreted as long-run frequencies")),
          tags$li(HTML("<strong>Bayesian:</strong> probability is interpreted as a subjective degree of belief"))
        ),
        h3("What is inference?"),
        tags$ul(
          tags$li("We are often concerned with estimating a population parameter"),
          tags$li(HTML("Estimating the shape of a distribution with no <strong>closed-form expression</strong>
                       e.g., by creating a histogram")),
          tags$li("Computing Bayesian point estimates, such as the mean, median, or mode"),
          tags$li("Computing Bayesian credible intervals for a parameter of interest")
        ),
        h3("Key steps in Bayesian inference"),
        tags$ul(
          tags$li(HTML("<strong>Choose a prior distribution</strong> based on our initial beliefs about the parameters of interest")),
          tags$li(HTML("<strong>Collect data</strong> that's representative of our population of interest")),
          tags$li(HTML("<strong>Compute posterior distribution</strong> to update our prior beliefs based on the empirical observations")),
          tags$li(HTML("<strong>Repeat.</strong> Bayesian inference is potentially an iterative process"))
        )
      )
    ),
    tabPanel(
      withMathJax(),
      tags$div(HTML("<script type='text/x-mathjax-config'>
                MathJax.Hub.Config({
                tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
                });
                </script>
                ")),
      title = "Visualizing Bayesian Inference",
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h3("Bayesian Inference"),
          p(style = "text-align: justify;", "This demo illustrates the Bayesian process of updating
            your beliefs after collecting data."),
          p(style = "text-align: justify;", "The parameter of interest is the probability of success
            in a Bernoulli process."),
          tags$ul(
            tags$li(HTML("<strong>Likelihood:</strong> Binomial distribution")),
            tags$li(HTML("<strong>Prior:</strong> Beta distribution")),
            tags$li(HTML("<strong>Posterior:</strong> Beta distribution")),
          ),
          hr(),
          h4("Prior Distribution"),
          splitLayout(
            textInput(inputId = "bayes_alpha", label = "$\\alpha$ Parameter:", value = 2),
            textInput(inputId = "bayes_beta", label = "$\\beta$ Parameter:", value = 5),
            cellArgs = list(style = "padding: 5px")
          ),
          h4("Likelihood Based on Data"),
          splitLayout(
            textInput(inputId = "bayes_successes", label = "Successes:", value = 15),
            textInput(inputId = "bayes_trials", label = "Trials:", value = 20),
            cellArgs = list(style = "padding: 5px")
          ),
          div(
            class = "text-center",
            actionButton(inputId = "run_bayes", label = "Find Posterior")
          )
        ),
        mainPanel(
          width = 9,
          mainPanel(
            br(),
            uiOutput("bayes_spinner")
          )
        )
      )
    )
  )
)

# =================================================================================================
# MCMC INTRODUCTION PAGE
# =================================================================================================

introduction_page <- wellPanel(
  class = "slide",
  h1("Introduction to Markov Chain Monte Carlo"),
  h3("What is MCMC?"),
  tags$ul(
    tags$li("An ensemble of methods for sampling from high-dimensional probability distributions"),
    tags$li("Accounts for the fact that the observations may depend on one another")
  ),
  h3("Why do we need MCMC?"),
  tags$ul(
    tags$li("Useful in Bayesian inference: can sample from a posterior distribution that has no closed-form expression"),
    tags$li("Also useful in situations where we cannot easily draw independent samples from the posterior")
  ),
  h3("Key components of MCMC"),
  tags$ul(
    tags$li("Monte Carlo sampling"),
    tags$li("Markov chains")
  )
)

# =================================================================================================
# MONTE CARLO PAGE
# =================================================================================================

monte_carlo_page <- tabItem(
  br(),
  tabsetPanel(
    tabPanel(
      title = "What is Monte Carlo Sampling?",
      wellPanel(
        class = "slide",
        h1("Monte Carlo Sampling"),
        h3("How does it work?"),
        tags$ul(
          tags$li(HTML("<strong>Law of Large Numbers:</strong> the sample mean converges in probability to the theoretical mean")),
          tags$li(HTML("The Law of Large Numbers requires <strong>finite mean and variance</strong> for convergence to occur")),
          tags$li("We can use repeated sampling to estimate a probability distribution"),
          tags$li("Accuracy increases as the sample size grows")
        ),
        h3("Why is this useful?"),
        tags$ul(
          tags$li("We can estimate the shape of a distribution using repeated sampling, with no need for a closed-form expression"),
          tags$li("Using the shape of a posterior distribution, we can find Bayesian estimators")
        ),
        h3("Limitations"),
        tags$ul(
          tags$li(HTML("Does not perform well when sampling from <strong>high-dimensional</strong> distributions")),
          tags$li(HTML("Only allows us to draw <strong>independent</strong> samples from a distribution"))
        )
      )
    ),
    tabPanel(
      title = "Visualizing Monte Carlo",
      br(),
      sidebarLayout(
        sidebarPanel(
          class = "sidebar",
          width = 3,
          h3("Monte Carlo Simulation"),
          p(style = "text-align: justify;",
            HTML("This is a demonstration of the <strong>Law of Large Numbers</strong>,
                 which guarantees the convergence of Monte Carlo sampling.")),
          hr(),
          radioButtons(
            inputId = "choose_distribution",
            label = "Choose distribution:",
            choices = c("Normal distribution", "Cauchy distribution"),
            selected = "Normal distribution"
          ),
          textInput(inputId = "num_samples", label = "Number of Samples:", value = 20),
          textInput(inputId = "monte_carlo_sample_size", label = "Sample Size:", value = 10000),
          uiOutput("distribution_controls"),
          checkboxInput(inputId = "monte_carlo_log", label = "Use logarithmic x-axis", value = TRUE),
          div(
            class = "text-center",
            actionButton(inputId = "run_monte_carlo", class = "button", label = "Run")
          )
        ),
        mainPanel(
          width = 9,
          uiOutput("monte_carlo_spinner")
        )
      )
    )
  )
)

# =================================================================================================
# MARKOV CHAINS PAGE
# =================================================================================================

markov_chains_page <- tabItem(
  br(),
  tabsetPanel(
    tabPanel(
      withMathJax(),
      tags$div(HTML("<script type='text/x-mathjax-config'>
                MathJax.Hub.Config({
                tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
                });
                </script>
                ")),
      title = "What are Markov Chains?",
      wellPanel(
        class = "slide",
        h1("Markov Chains"),
        h3("High-level idea"),
        tags$ul(
          tags$li("Markov chains enable us to generate a sequence of random variables where the value of one variable depends on the value of the previous one through some fixed probability rule"),
          tags$li("To determine the next value in the sequence, we have to know the current value as well as the probability rule")
        ),
        h3("Mathematical definition"),
        tags$ul(
          tags$li("Consider a sequence of random variables $\\{X_i\\}$ indexed by a variable $i$, and defined on a state space $\\Omega$"),
          tags$li(HTML("$\\{X_i\\}$ is a Markov chain if and only if the <strong>memoryless property</strong> holds, i.e., the conditional probability of $X_{n + 1}$ depends only on $X_n$: $$P\\left(X_{n + 1} = x \\mid X_n = x_n, \\dots, X_0 = x_0\\right) = P\\left(X_{n + 1} = x \\mid X_n = x_n\\right)$$")),
          tags$li("We uniquely determine a Markov chain with a transition matrix $\\mathbf{P}$ of probabilities, defined by $$P_{ij} = P\\left(X_{i + 1} = x_j \\mid X_i = x_i\\right) \\quad \\text{where} \\quad \\sum \\limits_{j} P_{ij} = 1$$")
        ),
        h3("Example"),
        tags$ul(
          tags$li(HTML("<strong>Random Walk:</strong> the next position depends only on the current position and on the direction in which we move"))
        )
      )
    ),
    tabPanel(
      withMathJax(),
      tags$div(HTML("<script type='text/x-mathjax-config'>
                MathJax.Hub.Config({
                tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
                });
                </script>
                ")),
      title = "Visualizing Markov Chains",
      br(),
      sidebarLayout(
        sidebarPanel(
          class = "sidebar",
          width = 3,
          h3("Markov Chains"),
          p(style = "text-align: justify;", "This demo illustrates a three-state Markov process that converges to a stationary
                 distribution."),
          p(style = "text-align: justify;", "The state space is given by $\\Omega = \\{1, 2, 3\\}$."),
          hr(),
          textInput(inputId = "steps", label = "Number of Steps:", value = 5000),
          matrixInput(inputId = "transition_matrix", label = "Transition Probability Matrix:",
                      value = matrix(c(0.5, 0.25, 0.25, 0.2, 0.1, 0.7, 0.25, 0.25, 0.5),
                                     byrow = TRUE, nrow = 3),
                      class = "numeric", rows = list(names = FALSE), cols = list(names = FALSE)),
          div(
            class = "text-center",
            actionButton(inputId = "run_markov_chains", class = "button", label = "Run")
          )
        ),
        mainPanel(
          width = 9,
          tabsetPanel(
            tabPanel("Convergence", br(), uiOutput("convergence_spinner")),
            tabPanel("First 100 States", br(), uiOutput("state_spinner")),
            tabPanel("Time in Each State", br(), uiOutput("time_spinner"))
          )
        )
      )
    )
  )
)

# =================================================================================================
# METROPOLIS-HASTINGS ALGORITHM PAGE
# =================================================================================================

algorithm_page <- tabItem(
  br(),
  tabsetPanel(
    tabPanel(
      title = "Introduction to Metropolis-Hastings",
      wellPanel(
        class = "slide",
        h1("The Metropolis-Hastings Algorithm"),
        h3("What is it?"),
        tags$ul(
          tags$li("The most popular MCMC algorithm"),
          tags$li("Can generate random samples based solely on the kernel of a probability distribution")
        ),
        h3("Conditions for convergence"),
        tags$ul(
          tags$li(HTML("<strong>Existence of the stationary distribution</strong>:
                       for this condition to be satisfied, it is sufficient but not necessary that the Markov
                       chain satisfies the <strong>principle of detailed balance</strong> i.e., every transition
                       must be reversible")),
          tags$li(HTML("<strong>Uniqueness of the stationary distribution</strong>: guaranteed by the
                       ergodicity of Markov chains"))
        ),
        h3("Pseudocode"),
        tags$ul(
          tags$li("Given an initial value $x^{(0)}$, a maximum number of iterations $N$, and a proposal distribution $g\\left(x' \\mid x\\right)$"),
          tags$li("From $i = 0$ to $N$:"),
          tags$ul(
            tags$li("Choose $x^{*}$ at random from the proposal distribution $g\\left(x' \\mid x\\right)$"),
            tags$li("Calculate the acceptance probability $p_a = A\\left(x^{\\,*}, x^{(i)}\\right)$"),
            tags$li("$x^{(i + 1)} = x^{*}$ with probability $p_a$, and $x^{(i + 1)} = x^{(i)}$ with probability $1 - p_a$"),
          ),
          tags$li(HTML("Discard the <strong>burn-in period</strong> $($the first $k$ observations$)$. Our sample is given by: $\\vec{X} = x^{(k)}, x^{(k + 1)}, \\dots, x^{(N)}$"))
        )
      )
    ),
    tabPanel(
      withMathJax(),
      tags$div(HTML("<script type='text/x-mathjax-config'>
                MathJax.Hub.Config({
                tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
                });
                </script>
                ")),
      title = "Formal Derivation",
      wellPanel(
        class = "slide",
        p(
          "All states of the Markov chain are based on the principle of detailed balance:
          $$P\\left(x' \\mid x\\right) P\\left(x\\right) = P\\left(x \\mid x'\\right) P\\left(x'\\right) \\Longleftrightarrow \\frac{P\\left(x' \\mid x\\right)}{P\\left(x \\mid x'\\right)} = \\frac{P\\left(x'\\right)}{P\\left(x\\right)}$$
          Each transition from one state to the next is split up into 2 phases:"
        ),
        tags$ul(
          tags$li(HTML("<strong>The proposal:</strong> the proposal distribution $g\\left(x' \\mid x\\right)$
                       defines the probability of proposing a state $x'$ given that we are
                       currently in state $x$")),
          tags$li(HTML("<strong>The acceptance/rejection:</strong> the acceptance distribution $A\\left(x', x\\right)$
                       is the probability of accepting the proposed state $x'$ from the current state $x$.
                       If a state is accepted, we transition to it; otherwise we stay in the current state at the next step"))
        ),
        p(
          "The proposal and the acceptance/rejection steps are independent, so the transition probability
          can be written as
          $$P\\left(x' \\mid x\\right) = g\\left(x' \\mid x\\right) A\\left(x', x\\right)$$
          Substituting this expression for the transition probability into the equation for detailed balance, we obtain
          $$\\frac{A\\left(x', x\\right)}{A\\left(x, x'\\right)} = \\frac{g\\left(x \\mid x'\\right) P\\left(x'\\right)}{g\\left(x' \\mid x\\right) P\\left(x\\right)}$$
          At every step, we must choose an acceptance ratio satisfying the above equation. Under Metropolis-Hastings, we define
          $$A\\left(x', x\\right) = \\min \\left \\{1, \\frac{g\\left(x \\mid x'\\right) P\\left(x'\\right)}{g\\left(x' \\mid x\\right) P\\left(x\\right)}\\right\\}$$"
        )
      )
    )
  )
)

# =================================================================================================
# INTERACTIVE DEMO PAGE
# =================================================================================================

gmm_demo_page <- tabItem(
  br(),
  sidebarLayout(
    sidebarPanel(
      class = "sidebar",
      width = 3,
      h2("Interactive Demo"),
      p(style = "text-align: justify;", "This example uses a Gaussian mixture model as the posterior,
        and a Gaussian proposal distribution."),
      hr(),
      textInput(inputId = "max_iter", label = "Number of Steps:", value = 5000),
      textInput(inputId = "proposal_sd", label = "Std. Dev. of Proposal Distribution:", value = 4),
      textInput(inputId = "x_init", label = "Initial Value:", value = -10),
      div(
        class = "text-center",
        actionButton(inputId = "run_mcmc", class = "button", label = "Run")
      )
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Histogram of Sample", br(), uiOutput("histogram_spinner")),
        tabPanel("Sample Values Step-by-Step", br(), uiOutput("markov_chain_spinner")),
        tabPanel("Decision for Proposed Values", br(), uiOutput("decision_spinner")),
        tabPanel("Acceptance Rate", br(), uiOutput("acceptance_rate_spinner"))
      )
    )
  )
)

# =================================================================================================
# INTERACTIVE DEMO PAGE
# =================================================================================================

gamma_demo_page <- tabItem(
  br(),
  sidebarLayout(
    sidebarPanel(
      withMathJax(),
      tags$div(HTML("<script type='text/x-mathjax-config'>
                MathJax.Hub.Config({
                tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
                });
                </script>
                ")),
      class = "sidebar",
      width = 3,
      h2("Interactive Demo"),
      p(style = "text-align: justify;", "This example performs the entire Bayesian inference process
        by modeling the mean number of sunspots per month using a $\\text{Gamma}\\left(\\alpha, \\beta\\right)$
        likelihood, a uniform prior for $\\alpha$ and $\\beta$, and a Gaussian proposal."),
      hr(),
      textInput(inputId = "max_iter_2d", label = "Number of Steps:", value = 5000),
      textInput(inputId = "proposal_sd_2d_first", label = "Std. Dev. of Proposal Distribution (1):", value = 0.05),
      textInput(inputId = "proposal_sd_2d_second", label = "Std. Dev. of Proposal Distribution (2):", value = 0.5),
      splitLayout(
        textInput(inputId = "x_init_2d_alpha", label = "$\\alpha$ Initial Value:", value = 4),
        textInput(inputId = "x_init_2d_beta", label = "$\\beta$ Initial Value:", value = 10)
      ),
      div(
        class = "text-center",
        actionButton(inputId = "run_2d", class = "button", label = "Run")
      )
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Estimated Distribution", br(), uiOutput("estimated_distribution_spinner"))
     )
    )
  )
)

# =================================================================================================
# REFERENCES & FURTHER RESOURCES PAGE
# =================================================================================================

references_page <- wellPanel(
  class = "slide",
  h1("References"),
  tags$ul(
    tags$li(tags$a("https://jellis18.github.io/post/2018-01-02-mcmc-part1/")),
    tags$li(tags$a("https://www.tweag.io/blog/2019-10-25-mcmc-intro1/")),
    tags$li(tags$a("http://people.duke.edu/~ccc14/sta-663/MCMC.html")),
    tags$li(tags$a("https://github.com/Joseph94m/MCMC/blob/master/MCMC.ipynb")),
    tags$li(tags$a("https://nicercode.github.io/guides/mcmc/"))
  )
)
