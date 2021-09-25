ui <- navbarPage(
    title = "Markov Chain Monte Carlo",
    tabPanel("Welcome", welcome_page),
    tabPanel("Bayesian Inference", bayes_page),
    navbarMenu(
        "Markov Chain Monte Carlo",
        tabPanel("Introduction", introduction_page),
        tabPanel("Monte Carlo", monte_carlo_page),
        tabPanel("Markov Chains", markov_chains_page)
    ),
    tabPanel("Metropolis-Hastings Algorithm", algorithm_page),
    navbarMenu(
        "Interactive Demos",
        tabPanel("One-Dimensional Distribution", gmm_demo_page),
        tabPanel("Two-Dimensional Distribution", gamma_demo_page)
    ),
    tabPanel("References & Further Resources", references_page),
    inverse = TRUE,
    theme = "style.css",
    footer = div(id = "footer", "Shu Amano, Maria-Cristiana Gîrjău, Leah Johnson", align = "left")
)
