devtools::install_github('stablemarkets/ChiRP')
library(ChiRP)
library(BART)
library(tidyverse)
library(gtools)
library(rstanarm)
library(BayesTree)
library(rstan)
library(LaplacesDemon)
library(brms)
library(broom)
set.seed(3)

## function that takes a set N x M matrix of M posterior predictions for each of 
## the N subjects under intervention 1 (mu_a1) and intervention 0 (mu_a0)
## Performs Bayesian bootstrap and returns vector containing M posterior draws
## of Psi.

# Taken from https://github.com/stablemarkets/intro_bayesian_causal/blob/master/Nonparametrics/npbayes_ATE.R
bayes_boot = function(mu_a1, mu_a0){
    n = nrow(mu_a1)
    M = ncol(mu_a1)
    psi_post = numeric(M)
    
    for(m in 1:M){
        bb_weights = rdirichlet( 1, rep(1, n) )
        psi_post[m] = sum(bb_weights*( mu_a1[, m] - mu_a0[ , m] ))
    }
    return(psi_post)
}

# Read in data
# Make treatment 0/1 variable
nets <- read_csv("https://evalf21.classes.andrewheiss.com/data/mosquito_nets.csv") %>%
    mutate(net = ifelse(net, 1, 0)) %>%
    as.data.frame()

x_train <-
    nets %>%
    select(net, income, health, temperature) %>%
    as.data.frame()
# Create "test" data where half of data has treatment = 0 and other half has treatment = 1
x_test0 <- x_train %>% mutate(net = 0)
x_test1 <- x_train %>% mutate(net = 1)
x_test <- bind_rows(x_test0, x_test1)
y_train <- nets$malaria_risk

### First fit using standard Bayesian regression with rstanarm
stan_fit <- stan_glm(malaria_risk ~ net + income + health + temperature,
                     data = nets, seed = 123)
mu_a0 = t(posterior_predict(stan_fit, newdata = x_test0))
mu_a1 = t(posterior_predict(stan_fit, newdata = x_test1))
psi_lm = bayes_boot( mu_a1 = mu_a1, mu_a0 = mu_a0)
plot(density(psi_lm), xlab = "Average treatment effect of using a mosquito net (rstanarm)")


### Now using Dirichlet process
set.seed(2)
res=fDPMix(d_train = nets,
           formula = malaria_risk ~ net + income + health + temperature,
           d_test = x_test, 
           iter=1000, burnin=500,init_k = 10 )
N <- nrow(x_train)
psi_dp = bayes_boot( mu_a1 = res$test[(N+1):(2*N),], mu_a0 = res$test[1:N,])
plot(density(psi_dp), xlab = "Average treatment effect of using a mosquito net (rstanarm)")


### Now using BART
bart_res = wbart(x.train = x_train, y.train = y_train, 
                 x.test = x_test, ndpost = 2000, nskip = 500)
### Check convergence (looks good)
plot(bart_res$sigma, type = "l")
abline(v = 500, lwd = 2, col = "red")

bart_res_pred = t(bart_res$yhat.test)
mu_a0 = bart_res_pred[1:N, ]
mu_a1 = bart_res_pred[(N+1):(2*N), ]

psi_bart = bayes_boot( mu_a1 = mu_a1, mu_a0 = mu_a0)
plot(density(psi_bart), xlab = "Average treatment effect of using a mosquito net (BART)")


### Finally compare to propensity score approach
model_treatment <- brm(
    bf(net ~ income + temperature + health,
       decomp = "QR"),  # QR decomposition handles scaling and unscaling for us
    family = bernoulli(),  # Logistic regression
    data = nets,
    chains = 4, cores = 4, iter = 1000,
    seed = 1234, backend = "cmdstanr"
)
pred_probs_chains <- posterior_epred(model_treatment)
# Put each set of individual propensity scores into its own cell
pred_probs_nested <- pred_probs_chains %>% 
    # Convert this matrix to a data frame
    as_tibble(.name_repair = "unique") %>% 
    # Add a column for the draw number
    mutate(draw = 1:n()) %>% 
    # Make this long so that each draw gets its own row
    pivot_longer(-draw, names_to = "row", values_to = "prob") %>% 
    # Clean up the draw number 
    mutate(row = as.numeric(str_remove(row, "..."))) %>% 
    # Group by draw and nest all the scores in a cell
    group_by(draw) %>% 
    nest() %>% 
    ungroup()
outcome_models <- pred_probs_nested %>% 
    mutate(outcome_model = map(data, ~{
        # Add this version of propensity scores to the original data and calculate
        # weights. We could also do this prior to nesting everything.
        df <- bind_cols(nets, .) %>% 
            mutate(iptw = (net_num / prob) + ((1 - net_num) / (1 - prob)))
        
        # Create outcome model with this iteration of weights
        model <- lm(malaria_risk ~ net, data = df, weights = iptw)
    })) %>% 
    # Extract results
    mutate(tidied = map(outcome_model, ~tidy(.)),
           ate = map_dbl(tidied, ~filter(., term == "net") %>% pull(estimate)),
           ate_se = map_dbl(tidied, ~filter(., term == "net") %>% pull(std.error)))
plot(density(outcome_models$ate), xlab = "Average treatment effect of using a mosquito net (IPTW")
