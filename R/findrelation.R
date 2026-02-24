# Helper function for netSEMp1
## FF-1: Linear

netSEM_linear <- function(Resp, Var, x) {
  # Define significant digits
  sig_dig <- 3
  # Define formula for executing model
  formula <- glue::glue('{Resp} ~ {Var}')
  #define a simple linear model
  model <- purrr::exec("lm", formula, data = x)
  #model fitting and extracting the statistical analysis as tibble
  model_summary <- broom::glance(model) %>%
    dplyr::mutate(
      Resp = paste(enexpr(Resp)),
      Var = paste(enexpr(Var)),
      model_type = "Linear",
      .before = dplyr::everything()
    )
  # Calculating Pseudo R2
  pseudoR2 <- rcompanion::nagelkerke(model) %>%
    purrr::pluck("Pseudo.R.squared.for.model.vs.null") %>%
    tibble::as_tibble(rownames = "R2") %>%
    tidyr::pivot_longer(!'R2', values_to = 'value') %>%
    tidyr::pivot_wider(names_from = 'R2', values_from = 'value') %>%
    janitor::clean_names() %>%
    dplyr::rename(
      'pseudoR2.mcfadden' = 'mc_fadden',
      'pseudoR2.cox_snell' = 'cox_and_snell_ml',
      'pseudoR2.nagelkerke' = 'nagelkerke_cragg_and_uhler'
    ) %>%
    dplyr::transmute(pseudoR2.mcfadden, pseudoR2.cox_snell, pseudoR2.nagelkerke) %>%
    dplyr::mutate(model_type = "Linear")
  # Generating Model Equation
  equation <- model %>%
    purrr::pluck('coefficients') %>%
    tibble::as_tibble_row() %>%
    dplyr::rename(Resp = '(Intercept)') %>%
    dplyr::mutate(dplyr::across(.cols = dplyr::everything(), ~ signif(.x, sig_dig))) %>%
    dplyr::mutate(model_type = 'Linear') %>%
    dplyr::mutate(model_eq = glue::glue('y ~ {Var}x + {Resp}', Var = signif(model$coefficients[[enexpr(Var)]], sig_dig))) %>%
    dplyr::select(c(model_eq, model_type))
  # Generating output table
  suppressMessages(model_summary <- model_summary %>%
                     dplyr::inner_join(pseudoR2, by = "model_type") %>%
                     dplyr::inner_join(equation, by = "model_type"))

  return(model_summary)
}

## FF-2: Quadratic with Linear Term

netSEM_quad <- function(Resp, Var, x) {
  # Define a general quadratic model
  sig_dig <- 3
  # Define a simple linear model
  formula <-
    glue::glue('{Resp} ~ {Var} + I({Var}^2)')
  model <- purrr::exec("lm", formula, data = x)
  # Model fitting and extracting the statistical analysis as tibble
  model_summary <- broom::glance(model) %>%
    dplyr::mutate(
      Resp = paste(enexpr(Resp)),
      Var = paste(enexpr(Var)),
      model_type = "Quad",
      .before = dplyr::everything()
    )
  # Calculate Pseudo R2
  pseudoR2 <- rcompanion::nagelkerke(model) %>%
    purrr::pluck("Pseudo.R.squared.for.model.vs.null") %>%
    tibble::as_tibble(rownames = "R2") %>%
    tidyr::pivot_longer(!'R2', values_to = 'value') %>%
    tidyr::pivot_wider(names_from = 'R2', values_from = 'value') %>%
    janitor::clean_names() %>%
    dplyr::rename('pseudoR2.mcfadden' = 'mc_fadden',
                  'pseudoR2.cox_snell' = 'cox_and_snell_ml',
                  'pseudoR2.nagelkerke' = 'nagelkerke_cragg_and_uhler'
    ) %>%
    dplyr::transmute(pseudoR2.mcfadden, pseudoR2.cox_snell, pseudoR2.nagelkerke) %>%
    dplyr::mutate(model_type = "Quad")
  # Generate Equation form
  equation <- model %>%
    purrr::pluck('coefficients') %>%
    tibble::as_tibble_row() %>%
    dplyr::rename(Resp = '(Intercept)',
                  Var = enexpr(Var),
                  Var_sq = glue::glue('I({enexpr(Var)}^2)')) %>%
    dplyr::mutate(dplyr::across(.cols = dplyr::everything(), ~ signif(.x, sig_dig)))  %>%
    dplyr::mutate(model_type = 'Quad') %>%
    dplyr::mutate(model_eq = glue::glue('y ~ {Var_sq}x^2 + {Var}x + {Resp}')) %>%
    dplyr::select(c(model_eq, model_type))
  suppressMessages(model_summary <- model_summary %>%
                     dplyr::inner_join(pseudoR2, by = "model_type") %>%
                     dplyr::inner_join(equation, by = "model_type"))
  # Outpur the model summary
  return(model_summary)
}

## FF-3: Quadratic with No Linear Term

netSEM_quadNoLinear <- function(Resp, Var, x) {
  # Define a general quadratic model
  sig_dig <- 3
  # Define a
  formula <-
    glue::glue('{Resp} ~ {Var}^2')
  model <- purrr::exec("lm", formula, data = x)
  #model fitting and extracting the statistical analysis as tibble
  model_summary <- broom::glance(model) %>%
    dplyr::mutate(
      Resp = paste(enexpr(Resp)),
      Var = paste(enexpr(Var)),
      model_type = "SQuad",
      .before = dplyr::everything()
    )
  # calculate Pseudo R2
  pseudoR2 <- rcompanion::nagelkerke(model) %>%
    purrr::pluck("Pseudo.R.squared.for.model.vs.null") %>%
    tibble::as_tibble(rownames = "R2") %>%
    tidyr::pivot_longer(!'R2', values_to = 'value') %>%
    tidyr::pivot_wider(names_from = 'R2', values_from = 'value') %>%
    janitor::clean_names() %>%
    dplyr::rename('pseudoR2.mcfadden' = 'mc_fadden',
                  'pseudoR2.cox_snell' = 'cox_and_snell_ml',
                  'pseudoR2.nagelkerke' = 'nagelkerke_cragg_and_uhler'
    ) %>%
    dplyr::transmute(pseudoR2.mcfadden, pseudoR2.cox_snell, pseudoR2.nagelkerke) %>%
    dplyr::mutate(model_type = "SQuad")
  # Generate Equation form
  equation <- model %>%
    purrr::pluck('coefficients') %>%
    tibble::as_tibble_row() %>%
    dplyr::rename(Resp = '(Intercept)',
                  Var= glue::glue('{enexpr(Var)}')) %>%
    dplyr::mutate(dplyr::across(.cols = dplyr::everything(), ~ signif(.x, sig_dig)))  %>%
    dplyr::mutate(model_type = 'SQuad') %>%
    dplyr::mutate(model_eq = glue::glue('y ~ {Var}x^2 + {Resp}')) %>%
    dplyr::select(c(model_eq, model_type))
  suppressMessages(model_summary <- model_summary %>%
                     dplyr::inner_join(pseudoR2, by = "model_type") %>%
                     dplyr::inner_join(equation, by = "model_type"))
  # Outpur the model summary
  return(model_summary)
}

## FF-4: Exponential

netSEM_exponential <- function(Resp, Var, x) {
  sig_dig <- 3
  # Define a
  formula <-
    glue::glue('{Resp} ~ exp({Var})')
  #test for condition where Var data has negative or INF values which will cause lm(exp) model to fail
  Var_flag <- x[[enexpr(Var)]]
  model_flag <- all(exp(Var_flag[!is.na(Var_flag)]) != Inf)
  #if flag is TRUE, will return model statistics
  if(model_flag){
    # Define the exponential model
    model <- purrr::exec("lm", formula, data = x)
    #model fitting and extracting the statistical analysis as tibble
    model_summary <- broom::glance(model) %>%
      dplyr::mutate(
        Resp = paste(enexpr(Resp)),
        Var = paste(enexpr(Var)),
        model_type = "Exponential",
        .before = dplyr::everything()
      )
    # calculate Pseudo R2
    pseudoR2 <- rcompanion::nagelkerke(model) %>%
      purrr::pluck("Pseudo.R.squared.for.model.vs.null") %>%
      tibble::as_tibble(rownames = "R2") %>%
      tidyr::pivot_longer(!'R2', values_to = 'value') %>%
      tidyr::pivot_wider(names_from = 'R2', values_from = 'value') %>%
      janitor::clean_names() %>%
      dplyr::rename('pseudoR2.mcfadden' = 'mc_fadden',
                    'pseudoR2.cox_snell' = 'cox_and_snell_ml',
                    'pseudoR2.nagelkerke' = 'nagelkerke_cragg_and_uhler'
      ) %>%
      dplyr::transmute(pseudoR2.mcfadden, pseudoR2.cox_snell, pseudoR2.nagelkerke) %>%
      dplyr::mutate(model_type = "Exponential")
    # Generate Equation form
    equation <- model %>%
      purrr::pluck('coefficients') %>%
      tibble::as_tibble_row() %>%
      dplyr::rename(Resp = '(Intercept)',
                    Var= glue::glue('exp({enexpr(Var)})')) %>%
      dplyr::mutate(dplyr::across(.cols = dplyr::everything(), ~ signif(.x, sig_dig)))  %>%
      dplyr::mutate(model_type = 'Exponential') %>%
      dplyr::mutate(model_eq = glue::glue('y ~ {Var}exp(x) + {Resp}')) %>%
      dplyr::select(c(model_eq, model_type))
    suppressMessages(model_summary <- model_summary %>%
                       dplyr::inner_join(pseudoR2, by = "model_type") %>%
                       dplyr::inner_join(equation, by = "model_type"))
  } else {
    #sets the condition when the model is false, outputting an empty tibble
    model_summary <- tibble::tibble(NA, model_type = "Exponential")
  }
  # Outpur the model summary
  return(model_summary)
}

## FF-5: Logarithmic

netSEM_log <- function(Resp, Var, x) {
  # Define Significant Figures
  sig_dig <- 3
  # Define formula for generating model
  formula <-
    glue::glue('{Resp} ~ log({Var})')
  #subset variable column for model validity check
  Var_flag <- x[[enexpr(Var)]]
  #test for condition where Var data has negative or INF values which will cause lm(log) model to fail
  model_flag <- all(Var_flag > 0 & Var_flag < Inf, na.rm = TRUE)
  #if flag is TRUE, will return model statistics
  if(model_flag){
    # Define the log model
    model <- purrr::exec("lm", formula, data = x)
    #model fitting and extracting the statistical analysis as tibble
    model_summary <- broom::glance(model) %>%
      dplyr::mutate(
        Resp = paste(enexpr(Resp)),
        Var = paste(enexpr(Var)),
        model_type = "Log",
        .before = dplyr::everything()
      )
    # calculate Pseudo R2
    pseudoR2 <- rcompanion::nagelkerke(model) %>%
      purrr::pluck("Pseudo.R.squared.for.model.vs.null") %>%
      tibble::as_tibble(rownames = "R2") %>%
      tidyr::pivot_longer(!'R2', values_to = 'value') %>%
      tidyr::pivot_wider(names_from = 'R2', values_from = 'value') %>%
      janitor::clean_names() %>%
      dplyr::rename('pseudoR2.mcfadden' = 'mc_fadden',
                    'pseudoR2.cox_snell' = 'cox_and_snell_ml',
                    'pseudoR2.nagelkerke' = 'nagelkerke_cragg_and_uhler'
      ) %>%
      dplyr::transmute(pseudoR2.mcfadden, pseudoR2.cox_snell, pseudoR2.nagelkerke) %>%
      dplyr::mutate(model_type = "Log")
    # Generate Equation form
    equation <- model %>%
      purrr::pluck('coefficients') %>%
      tibble::as_tibble_row() %>%
      dplyr::rename(Resp = '(Intercept)',
                    Var= glue::glue('log({enexpr(Var)})')) %>%
      dplyr::mutate(dplyr::across(.cols = dplyr::everything(), ~ signif(.x, sig_dig)))  %>%
      dplyr::mutate(model_type = 'Log') %>%
      dplyr::mutate(model_eq = glue::glue('y ~ {Var}log(x) + {Resp}')) %>%
      dplyr::select(c(model_eq, model_type))
    suppressMessages(model_summary <- model_summary %>%
                       dplyr::inner_join(pseudoR2, by = "model_type") %>%
                       dplyr::inner_join(equation, by = "model_type"))
    # pseudoR2 <- rcompanion::nagelkerke(model) %>%
    #   purrr::pluck("Pseudo.R.squared.for.model.vs.null") %>%
    #   tibble::as_tibble(rownames = "R2") %>%
    #   pivot_longer(!'R2', values_to = 'value') %>%
    #   pivot_wider(names_from = 'R2', values_from = 'value') %>%
    #   janitor::clean_names() %>%
    #   dplyr::rename('pseudoR2.mcfadden' = 'mc_fadden',
    #                 'pseudoR2.cox_snell' = 'cox_and_snell_ml',
    #                 'pseudoR2.nagelkerke' = 'nagelkerke_cragg_and_uhler'
    #   ) %>%
    #   dplyr::transmute(pseudoR2.mcfadden, pseudoR2.cox_snell, pseudoR2.nagelkerke) %>%
    #   dplyr::mutate(model_type = "Log")
    # model_summary <- model_summary %>%
    #   inner_join(pseudoR2, by = "model_type")
  } else {
    #when the model is false, output an empty tibble
    model_summary <- tibble::tibble(NA, model_type = "Log")
  }
  # Outpur the model summary
  return(model_summary)
}

## FF-6: Square Root (non-linear)

netSEM_SqRt <- function(Resp, Var, x) {
  # Define Significant Figures
  sig_dig <- 3
  # Define formula for generating model
  formula <-
    glue::glue('{Resp} ~ sqrt({Var})')
  # define data to check for fail condition
  Var_flag <- x[[enexpr(Var)]]
  #test for condition where Var data has negative or INF values which will cause lm(sqrt) model to fail
  model_flag <- all(as.logical(Var_flag > 0 &
                                 Var_flag < Inf, na.rm = TRUE))
  #if flag is TRUE, will return model statistics
  if(model_flag){
    # Define the SqRoot model
    model <- purrr::exec("lm", formula, data = x)
    #model fitting and extracting the statistical analysis as tibble
    model_summary <- broom::glance(model) %>%
      dplyr::mutate(
        Resp = paste(enexpr(Resp)),
        Var = paste(enexpr(Var)),
        model_type = "SqRoot",
        .before = dplyr::everything()
      )
    # calculate Pseudo R2
    pseudoR2 <- rcompanion::nagelkerke(model) %>%
      purrr::pluck("Pseudo.R.squared.for.model.vs.null") %>%
      tibble::as_tibble(rownames = "R2") %>%
      tidyr::pivot_longer(!'R2', values_to = 'value') %>%
      tidyr::pivot_wider(names_from = 'R2', values_from = 'value') %>%
      janitor::clean_names() %>%
      dplyr::rename('pseudoR2.mcfadden' = 'mc_fadden',
                    'pseudoR2.cox_snell' = 'cox_and_snell_ml',
                    'pseudoR2.nagelkerke' = 'nagelkerke_cragg_and_uhler'
      ) %>%
      dplyr::transmute(pseudoR2.mcfadden, pseudoR2.cox_snell, pseudoR2.nagelkerke) %>%
      dplyr::mutate(model_type = "SqRoot")
    # Generate Equation form
    equation <- model %>%
      purrr::pluck('coefficients') %>%
      tibble::as_tibble_row() %>%
      dplyr::rename(Resp = '(Intercept)',
                    Var= glue::glue('sqrt({enexpr(Var)})')) %>%
      dplyr::mutate(dplyr::across(.cols = dplyr::everything(), ~ signif(.x, sig_dig)))  %>%
      dplyr::mutate(model_type = 'SqRoot') %>%
      dplyr::mutate(model_eq = glue::glue('y ~ {Var}x^0.5 + {Resp}')) %>%
      dplyr::select(c(model_eq, model_type))
    suppressMessages(model_summary <- model_summary %>%
                       dplyr::inner_join(pseudoR2, by = "model_type") %>%
                       dplyr::inner_join(equation, by = "model_type"))
  } else {
    #when the model is false, output an empty tibble
    model_summary <- tibble::tibble(NA, model_type = "SqRoot")
  }
  # Outpur the model summary
  return(model_summary)
}

## FF-7: Inverse Square Root (non-linear)

netSEM_InvSqRt <- function(Resp, Var, x) {
  # Define Significant Figures
  sig_dig <- 3
  # Define formula for generating model
  formula <-
    glue::glue('{Resp} ~ I(1/sqrt({Var}))')
  # define data to check for fail condition
  Var_flag <- x[[enexpr(Var)]]
  #test for condition where Var data has negative or INF values which will cause lm(sqrt) model to fail
  model_flag <- all(as.logical(Var_flag > 0 & Var_flag < Inf,
                               na.rm = TRUE))
  #if flag is TRUE, will return model statistics
  if(model_flag){
    # Define the Inverse SqRoot model
    model <- purrr::exec("lm", formula, data = x)
    #model fitting and extracting the statistical analysis as tibble
    model_summary <- broom::glance(model) %>%
      dplyr::mutate(
        Resp = paste(enexpr(Resp)),
        Var = paste(enexpr(Var)),
        model_type = "InvSqRoot",
        .before = dplyr::everything()
      )
    # calculate Pseudo R2
    pseudoR2 <- rcompanion::nagelkerke(model) %>%
      purrr::pluck("Pseudo.R.squared.for.model.vs.null") %>%
      tibble::as_tibble(rownames = "R2") %>%
      tidyr::pivot_longer(!'R2', values_to = 'value') %>%
      tidyr::pivot_wider(names_from = 'R2', values_from = 'value') %>%
      janitor::clean_names() %>%
      dplyr::rename('pseudoR2.mcfadden' = 'mc_fadden',
                    'pseudoR2.cox_snell' = 'cox_and_snell_ml',
                    'pseudoR2.nagelkerke' = 'nagelkerke_cragg_and_uhler'
      ) %>%
      dplyr::transmute(pseudoR2.mcfadden, pseudoR2.cox_snell, pseudoR2.nagelkerke) %>%
      dplyr::mutate(model_type = "InvSqRoot")
    # Generate Equation form
    equation <- model %>%
      purrr::pluck('coefficients') %>%
      tibble::as_tibble_row() %>%
      dplyr::rename(Resp = '(Intercept)',
                    Var = glue::glue('I(1/sqrt({enexpr(Var)}))')) %>%
      dplyr::mutate(dplyr::across(.cols = dplyr::everything(), ~ signif(.x, sig_dig)))  %>%
      dplyr::mutate(model_type = 'InvSqRoot') %>%
      dplyr::mutate(model_eq = glue::glue('y ~ {Var}1/x^0.5 + {Resp}')) %>%
      dplyr::select(c(model_eq, model_type))
    suppressMessages(model_summary <- model_summary %>%
                       dplyr::inner_join(pseudoR2, by = "model_type") %>%
                       dplyr::inner_join(equation, by = "model_type"))
  } else {
    #when the model is false, output an empty tibble
    model_summary <- tibble::tibble(NA, model_type = "InvSqRoot")
  }
  # Outpur the model summary
  return(model_summary)
}

## FF-8: NLS

netSEM_nls <- function(Resp, Var, x) {
  tryCatch({
    nls_formula <-
      glue::glue('{(Resp)} ~ a1 + a2 * exp(a3 * {(Var)})')

    resp <- x[[enexpr(Resp)]]
    var <- x[[enexpr(Var)]]

    #test for condition where Resp data has negative or INF values which will cause lm(log) model to fail
    Resp_flag <- x[[enexpr(Resp)]]
    model_flag <- all(as.logical(Resp_flag > 0 & Resp_flag < Inf,
                                 na.rm = TRUE))
    #if flag is TRUE, will return model statistics
    if(model_flag){

      a1_0 <- min(resp) / 2

      # # Take the log of nls_formula
      nls_log_formula <-
        glue::glue('log(({(Resp)} - {a1_0})) ~ {(Var)}')
      model_0 <- purrr::exec("lm", nls_log_formula, data = x)

      #NLS model
      model <- stats::nls(
        formula = nls_formula,
        data = x,
        start = (tibble::lst(a1 = a1_0, a2 = exp(coef(model_0)[1]), a3 = coef(model_0)[2])), #arbitrary starting points
        control = stats::nls.control(maxiter = 50, minFactor = 1/2^20, warnOnly = TRUE) #either converges or stops at 1000 iterations
      )

      model_coefficient <- (summary(model))$coefficients[, 1]

      a1 <- round(model_coefficient[1], 3)
      a2 <- round(model_coefficient[2], 3)
      a3 <- round(model_coefficient[3], 3)

      # Define the null model equation needed for determining the pseudo R^2 for nls-type models
      nullfunction <- function(Var, m) {m}
      # Define the null model
      null.model <- stats::nls(
        resp ~ nullfunction(var, m),
        data = x,
        start = tibble::lst(m = a1_0) #arbitrary constant
      )
      #calculates the pseudo R^2 values from three methods: McFadden, Cox and Snell and Nagelkerke and selects these models
      pseudoR2 <- rcompanion::nagelkerke(model, null = null.model) %>%
        purrr::pluck("Pseudo.R.squared.for.model.vs.null") %>%
        tibble::as_tibble(rownames = "R2") %>%
        tidyr::pivot_longer(!'R2', values_to = 'value') %>%
        tidyr::pivot_wider(names_from = 'R2', values_from = 'value') %>%
        janitor::clean_names() %>%
        dplyr::rename('pseudoR2.mcfadden' = 'mc_fadden',
                      'pseudoR2.cox_snell' = 'cox_and_snell_ml',
                      'pseudoR2.nagelkerke' = 'nagelkerke_cragg_and_uhler'
        ) %>%
        dplyr::transmute(pseudoR2.mcfadden, pseudoR2.cox_snell, pseudoR2.nagelkerke) %>%
        dplyr::mutate(model_type = "NLS")

      model_summary <- {
        #model fitting and extracting the statistical analysis as tibble
        result1 <- broom::glance(model) %>%
          dplyr::mutate(
            Resp = paste(enexpr(Resp)),
            Var = paste(enexpr(Var)),
            model_type = "NLS",
            .before = dplyr::everything()
          ) %>%
          dplyr::mutate(model_eq = glue::glue('{enexpr(Resp)} ~ {a1} + {a2} * exp({a3} * {enexpr(Var)})'))
        # coefficient of parameters in NLS model
        result2 <- broom::tidy(model) %>%
          dplyr::transmute(term, statistic, p.value) %>%
          dplyr::mutate(model_type = "NLS")
        suppressMessages(combined <- dplyr::full_join(
          result1,
          result2,
        ) %>%
          dplyr::full_join(
            pseudoR2
          ) %>%
          dplyr::select(-c('isConv', 'finTol')))
      }
    } else {
      #when the model is false, output an empty tibble
      model_summary <- tibble::tibble(NA, model_type = "NLS")
    }

    return(model_summary)},
    error = function(e){
      message('Resp = ',Resp,', Var = ',Var,' not good fit for nls.')
      model_summary <- tibble::tibble(NA, model_type = "NLS")
    }
  )
}

## FF-9: Change Point

netSEM_CP <- function(Resp, Var, x) {
  tryCatch({
    sig_dig <- 3

    resp <- x[[{{Resp}}]]
    var <- x[[{{Var}}]]
    # formula <- glue::glue('{Resp} ~ {Var}')
    # seg.formula <- glue::glue('~ {var}')
    # model.obj <- purrr::exec("lm", formula, data = x)
    # Define the change point model
    # model <- purrr::exec(segmented.lm,
    #                      obj = model.obj,
    #                      seg.Z = ~ var,
    #                      psi = NA,
    #                      control = segmented::seg.control(
    #                        K = 1,
    #                        fix.npsi = FALSE,
    #                        n.boot = 0,
    #                        it.max = 20
    #                      ))

    # Define the Change Point model
    model <- segmented::segmented(
      obj = lm(resp ~ var, data = x),
      seg.Z = ~ var,
      psi = NA,
      control = segmented::seg.control(
        K = 1,
        fix.npsi = F,
        n.boot = 0,
        it.max = 20
      )
    )

    model_coefficient <- (summary(model))$coefficients[, 1]

    #sets condition when the model is true
    model_TRUE <- broom::glance(model) %>%
      dplyr::mutate(
        Resp = paste(enexpr(Resp)),
        Var = paste(enexpr(Var)),
        model_type = "CP"
      ) %>%
      dplyr::select(Resp, Var, model_type, 1:12)

    #calculate Pseudo R2
    pseudoR2 <- rcompanion::nagelkerke(model) %>%
      purrr::pluck("Pseudo.R.squared.for.model.vs.null") %>%
      tibble::as_tibble(rownames = "R2") %>%
      tidyr::pivot_longer(!'R2', values_to = 'value') %>%
      tidyr::pivot_wider(names_from = 'R2', values_from = 'value') %>%
      janitor::clean_names() %>%
      dplyr::rename('pseudoR2.mcfadden' = 'mc_fadden',
                    'pseudoR2.cox_snell' = 'cox_and_snell_ml',
                    'pseudoR2.nagelkerke' = 'nagelkerke_cragg_and_uhler'
      ) %>%
      dplyr::transmute(pseudoR2.mcfadden, pseudoR2.cox_snell, pseudoR2.nagelkerke) %>%
      dplyr::mutate(model_type = "CP")

    #generate Equation form
    #when there is a change point
    if (length(model$coefficients) > 3) {
      # Get coefficients of model
      coefs <- coef(model)
      n_coefs <- length(coefs)

      # Count how many breakpoints exist
      # Coeffs: Intercept, var, U1.var, psi1.var, U2.var, ...
      n_psi <- (n_coefs-2)/2

      if (n_psi >= 1) {
        # Has breakpoint(s)
        intercept <- coefs["(Intercept)"]
        slope1 <- coefs["var"]

        # Start with intercept
        eq <- paste0(Resp, " = ", format(intercept, scientific = FALSE, digits = sig_dig))

        # Add first slope term
        if (slope1 >= 0) {
          eq <- paste0(eq, " + ", format(slope1, scientific = FALSE, digits = sig_dig), "*", Var)
        } else {
          eq <- paste0(eq, " ", format(slope1, scientific = FALSE, digits = sig_dig), "*", Var)
        }

        # Add breakpoint terms
        for (i in seq_len(n_psi)) {
          slope_diff <- coefs[paste0("U", i, ".var")]
          psi <- coefs[paste0("psi", i, ".var")]

          # Build psi part with _c notation
          if (psi >= 0) {
            psi_part <- paste0("(", Var, " - ", format(psi, scientific = FALSE, digits = sig_dig), ")_c")
          } else {
            psi_part <- paste0("(", Var, " + ", format(abs(psi), scientific = FALSE, digits = sig_dig), ")_c")
          }

          # Add slope_diff term
          if (slope_diff >= 0) {
            eq <- paste0(eq, " + ", format(slope_diff, scientific = FALSE, digits = sig_dig), "*", psi_part)
          } else {
            eq <- paste0(eq, " ", format(slope_diff, scientific = FALSE, digits = sig_dig), "*", psi_part)
          }
        }

        equation <- tibble::tibble(
          model_eq = eq,
          model_type = "CP",
          n_breakpoints = n_psi
        )

      } else {
        # No changepoint - simple linear
        intercept <- coefs["(Intercept)"]
        slope <- coefs["var"]

        eq <- paste0(Resp, " = ", format(intercept, scientific = FALSE, digits = sig_dig))

        if (slope >= 0) {
          eq <- paste0(eq, " + ", format(slope, scientific = FALSE, digits = sig_dig), "*", Var)
        } else {
          eq <- paste0(eq, " ", format(slope, scientific = FALSE, digits = sig_dig), "*", Var)
        }

        equation <- tibble::tibble(
          model_eq = eq,
          model_type = "CP",
          n_breakpoints = 0
        )
      }

      suppressMessages(
        model_summary <- model_TRUE %>%
          dplyr::inner_join(pseudoR2, by = "model_type") %>%
          dplyr::inner_join(equation, by = "model_type")
      )
      return(model_summary)
    }},
    error = function(e){
      print(e)
      model_summary <- tibble::tibble(NA, model_type = "CP")
    }
  )

}

#function to summarize each model
collect_netSEM_models <- function(Resp, Var, x) {

  x1 <- netSEM_linear(Resp ,Var, x)
  x2 <- netSEM_quad(Resp, Var, x)
  x3 <- netSEM_quadNoLinear(Resp, Var, x)
  x4 <- netSEM_exponential(Resp, Var, x)
  x5 <- netSEM_log(Resp, Var, x)
  x6 <- netSEM_InvSqRt(Resp, Var, x)
  x7 <- netSEM_SqRt(Resp, Var, x)
  x8 <- netSEM_nls(Resp, Var, x)
  x9 <- netSEM_CP(Resp, Var, x)

  x_final <-
    dplyr::bind_rows(x1, x2, x3, x4, x5, x6, x7, x8, x9) %>%
    dplyr::relocate("model_type")

  return(x_final)
}
