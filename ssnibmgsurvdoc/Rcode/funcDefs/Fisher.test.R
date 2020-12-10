# A function to get the Fisher product; to combine p-values of all studies
# gets one p-value for each study
Fisher.test <- function(p) {
  p = p[!is.na(p)]
  Xsq <- -2*sum(log(p))
  p.val <- pchisq(Xsq, df = 2*length(p), lower.tail=FALSE)
  return(c(Xsq = Xsq, p.value = p.val))
}

