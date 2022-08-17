tosdif = function(fit1, fit2, g = NULL, alpha = NULL, rmseaa = NULL, rmseab = NULL, robust = TRUE) {
  
  # print the title
  cat("     Test of small differences in fit (McCallum et al., 2006)", "\n", "\n")
  
  # Get df from both model fits
  df1 <- fitMeasures(fit1, "df")                
  df2 <- fitMeasures(fit2, "df")  
  
  # dfa is always the model with more df
  dfa <- pmax(df1, df2)
  dfb <- pmin(df1, df2)
  
  # get n from the model fits and see if lavInspect gets the same sample size for both
  nobsA = lavInspect(fitA, what = "nobs")
  nobsB = lavInspect(fitB, what = "nobs")
  if (nobsA != nobsB) {
    warning("Warning: Unequal sample size for the two models. Preceed with caution. \n \n")
  }
  
  # if 'robust' is not a logical
  if (!is.logical(robust)) {
    stop("Argument 'robust' needs to bei either TRUE or FALSE. Calculation was stopped.")}
  
  
  # if g was not specified, set it to 1
  if(is.null(g)) {
    g = 1
  }
  
  # if alpha was not specified, then set it to 0.05
  if(is.null(alpha)) {
    alpha = 0.05
  }
  
  # if RMSEA (A or B) was not specified, set it to .06 and .05 respectively
  if(is.null(rmseaa)) {
    rmseaa = 0.06
  }
  
  if(is.null(rmseab)) {
    rmseab = 0.05
  }
  
  # if robust = TRUE or not specified, then use "chisq.scaled", else use "chisq"
  if(robust && !is.na(lavInspect(fit1, "fit")["chisq.scaled"])) {
    chidiff = abs(lavInspect(fit1, "fit")["chisq.scaled"] -
                    lavInspect(fit2, "fit")["chisq.scaled"])
  } else if (robust && is.na(lavInspect(fit1, "fit")["chisq.scaled"])) {
    chidiff = abs(lavInspect(fit1, "fit")["chisq"] -
                    lavInspect(fit2, "fit")["chisq"])
    cat("\n Warning: Robust measures were not available. Standard chi-squared was used instead. Use estimator = 'MLR' to get robust measures. \n \n")
  } else if (!robust) {
    chidiff = abs(lavInspect(fit1, "fit")["chisq"] -
                    lavInspect(fit2, "fit")["chisq"])
  }
  
  
  ddiff = dfa - dfb                                  # df difference
  fa = (dfa*rmseaa^2)/sqrt(g)                        # model A discrepancy fn value
  fb = (dfb*rmseab^2)/sqrt(g)                        # model B discrepancy fn value
  ncp = (nobsA - 1)*(fa - fb)                        # non-centrality parameter
  cval = qchisq(1 - alpha, df = ddiff, ncp = ncp)    # critical value from non-central chi^2
  sig = 1 - pchisq(chidiff, df = ddiff, ncp = ncp)   # p-value from non-central chi^2
  rm(ddiff, fa, fb, ncp, cval)                       # remove unnecessary objects
  
  # print the p-value
  cat("p =", round(sig,3), "\n \n")
  
  # Interpretation
  if(sig > 0.05) {
    cat("Interpretation: Test of small differences in fit is NOT significant. The models show at most a small difference in the model fit. Selecting the more parsimonious model (more df) is advised.")
  } else {
    cat("Interpretation: Test of small differences in fit is significant. The models show a difference in fit that is greater than the operationalised acceptable difference.")
  }
}