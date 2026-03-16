get_coef_re <- function(object, add_main = TRUE){

  Xcoef <- as.matrix(t(object$params$Br))
  sdXcoef <- as.matrix(t(getPredictErr(object)$Br))

  Xmain <- as.matrix(object$params$B)
  Xmain <- Xmain[match(colnames(Xcoef), rownames(Xmain))]
  
  if (add_main) {
    Xcoef <- sweep(Xcoef, 2, Xmain, FUN = '+')
  }
  
  lower <- Xcoef - 1.96 * sdXcoef
  upper <- Xcoef + 1.96 * sdXcoef
  
  Xcoef <- Xcoef %>% as.data.frame() %>% mutate(sp = row.names(.)) %>%
    pivot_longer(-sp, values_to = 'mean', names_to = 'var')
  lower <- lower %>% as.data.frame() %>% mutate(sp = row.names(.)) %>%
    pivot_longer(-sp, values_to = 'lower', names_to = 'var')
  
  out.data <- upper %>% as.data.frame() %>% mutate(sp = row.names(.)) %>%
    pivot_longer(-sp, values_to = 'upper', names_to = 'var') %>% inner_join(Xcoef) %>%
    inner_join(lower)
  
  return(list(coef_data = out.data, main_eff = Xmain))
}
