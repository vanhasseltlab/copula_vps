
run_grimsley <- function(n, wgt, scr, verbose = F) {
  
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
    scr = scr
  )
  
  # Run simulation
  df_sim <- rxSolve(
    mod, ev,
    iCov = covar,
    addDosing = T,
    addCov = T
  )
  
  if (verbose == T) {
    return(df_sim)
  } else {
    df_sim %>%
      select(id, time, conc = C1, wgt, scr) %>%
      distinct() %>% 
      return()
  }
  
}



