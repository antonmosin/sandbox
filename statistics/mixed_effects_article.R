# Helper functions for the mixed effects article
library(dplyr)
library(lme4)

run_simulation <- function(n_sims,
                           n_users,
                           gamma_prob,
                           prior_alpha,
                           prior_beta,
                           test_alpha=0.05,
                           test_alternative='two-sided'){
  is_significant <- vector(mode='logical', length=n_sims)
  pval_model <- vector(mode='logical', length=n_sims)

  for (i in 1:n_sims){
    df <- generate_data(n_users = n_users,
                        gamma_prob = gamma_prob,
                        prior_alpha = prior_alpha,
                        prior_beta = prior_beta)
    group_assign <- data.frame(user_id = seq(n_users),
                               group = ifelse(runif(n_users) > 0.5, 'control', 'test'))
    df <- df %>% inner_join(group_assign, by='user_id')

    summary <- df %>%
      group_by(group) %>%
      summarise(avg_on_time = mean(is_on_time),
                n_users = length(unique(user_id)),
                n_tasks = n())

    test_res <- test_prop(p1 = summary[summary$group=='control', ] %>% pull(avg_on_time),
                          p2 = summary[summary$group=='test', ] %>% pull(avg_on_time),
                          n1 = summary[summary$group=='control', ] %>% pull(n_tasks),
                          n2 = summary[summary$group=='test', ] %>% pull(n_tasks),
                          alpha=test_alpha,
                          alternative=test_alternative)
    is_significant[i] = test_res$is_significant

    mod <- glmer(is_on_time ~ group + (1 | user_id), data=df, family=binomial)
    pval_model[i] <- summary(mod)$coefficients['grouptest', 'Pr(>|z|)']
  }
  return(list(naive_test = is_significant, mixed_model = pval_model < test_alpha))
}

generate_data <- function(n_users,
                          gamma_prob,
                          prior_alpha,
                          prior_beta){
  tasks_per_user <- rgeom(n_users, prob=gamma_prob)
  prior_prob_complet <- rbeta(n=n_users, shape1 = prior_alpha, shape2 = prior_beta)
  successes <- vector(mode="list", length=n_users)

  for (i in 1:n_users){
    if (tasks_per_user[i] == 0 ){
      next
    }
    successes[[i]] <- data.frame(is_on_time = rbinom(n=tasks_per_user[i],
                                                     p=prior_prob_complet[i],
                                                     size=1),
                                 user_id = rep(i, tasks_per_user[i]))
  }
  d <- bind_rows(successes)
  if (length(unique(d$user_id)) != sum(tasks_per_user > 0)){
    stop("Number of users with at least one task does not match!")
  }
  return(d)
}



bern_var <- function(p){
  return(p * (1 - p))
}
get_se <- function(p1, p2, n1, n2){
  return(sqrt(bern_var(p1) / n1 + bern_var(p2) / n2))
}

get_zscore <- function(p1, p2, n1, n2){
  d <- abs(p1 - p2)
  z <- d / get_se(p1, p2, n1, n2)
  return(z)
}

# normal approximation
test_prop <- function(p1, p2, n1, n2,
                      alpha=0.05,
                      alternative='two-sided'){
  z <- get_zscore(p1, p2, n1, n2)
  if (alternative == 'two-sided'){
    z_crit <- qnorm(1 - alpha / 2)
  } else {
    z_crit <- qnorm(1 - alpha)
  }
  return (list(is_significant = z > z_crit,
               z_score = z,
               z_critical = z_crit))
}
