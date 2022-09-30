#function for PK model (Grimsley) based on RxODE

run_grimsley <- function(n, wgt, scr, other_covariates = NULL, verbose = F) {
  
  if (length(wgt) != n | length(scr) != n) {
    stop("covariate input does not match the number of subjects")
  }
  
  require(tidyverse)
  require(RxODE)
  
  # n = number of individual 
  # wgt = body weight (kg)
  # scr = serum creatinine (umol/l)
  
  # Dosing event
  ev <- et(
    id = seq_len(n),
    amt = 15*wgt,
    cmt = 1,
    dur = 1*60,
    addl = 2,
    ii = 24*60
  ) %>%
    et(id = seq_len(n)) %>% 
    et(seq(0, 3*24*60, 20))
  
  # Model Grimsley 
  mod <- RxODE({
    
    CL <- 3.56 * wgt / scr / 60 # convert L/h to L/min
    V1 <- 0.669 * wgt
    
    d/dt(center) <- - CL/V1*center
    
    C1 <- center/V1
    
  })
  
  # Covariate table
  covar <- data.frame(
    id = seq_len(n),
    wgt = wgt,
    scr = scr, other_covariates
  )
  
  # Run simulation
  df_sim <- rxSolve(
    mod, ev,
    iCov = covar,
    addDosing = T,
    addCov = T
  )
  
  df_sim <- df_sim %>% 
    left_join(covar %>% dplyr::select(id, all_of(names(other_covariates))), by = "id")
  
  if (verbose == T) {
    return(df_sim)
  } else {
    df_sim %>%
      dplyr::select(id, time, conc = C1, wgt, scr, all_of(names(other_covariates))) %>%
      distinct() %>% 
      return()
  }
  
}



