# This is the MAGMA (Mark, Age, Genotype, Mixture Analysis) model created by Jim Jasper at ADF&G. 
# It was created for use with TBR and is now being adapted for use in other fisheries. 
# 
# Specifically, I am trying to make it a bit more robust when dealing with non-standard cases
# like no hatcheries, single subdistricts, etc., while retaining all previous functionality.
#
# This is the most up-to date version of the function (01/28/2019). 


magma <- function(nreps, burn, thin, nchains, wd){

#  nreps <- 4e4 ; burn <- floor(nreps / 2) ; thin <- 10 ; nchains <- 5 ;  wd <- "V:/Analysis/1_SEAK/Sockeye/Mixture/TBR/MAGMA/2017 MAGMA/New MAGMA Old Baseline" 

#  analysis_time <- magma(nreps, burn, thin, nchains, wd)

  if(!exists("update_base")){ update_base <- TRUE } # Argument for turning on/off baseline updating--useful for troubleshooting.

  if(!exists("nadapt")){ nadapt <- 100 }# A tuning parameter used to set the number of iterations to run before baseline updating begins. 

  if(!exists("sigfigs")){ sigfigs <- 16 } # Number of significant figures for writing output.  

  search_begin <- search() 

  packages <- c("abind", "rjags", "doParallel", "foreach")

  invisible(sapply(packages, function(pkg) while(!require(pkg, character.only = TRUE)) install.packages(pkg, repos = c("http://rstudio.org/_packages", "http://cran.rstudio.com"))) ) 

  setwd(wd)

  load("analysis/magma_data.RData")
 
  print(magma_data_names)

  paths <- lapply(seq(T), function(t){ 
    lapply(seq(D), function(d){ 
      lapply(seq(S[d]), function(s){ 
        lapply(seq(W), function(w){ 
          paste0("output/", years2run[t], "/D", districts[d], "/S", subdistricts[[d]][s], "/StatWeek", w) 
          }) 
        })
      }) 
    })

  build_pathes <- suppressWarnings(rapply(paths, dir.create, recursive = TRUE))

# Analysis  #####################################################################################################################################################################################################################################################################################################################################################################################################################################################################

  chains <- paste0("Chain", seq(nchains))

  cl <- makePSOCKcluster(nchains)

  registerDoParallel(cl, cores = nchains)

  beg_time <- Sys.time()

  dput(x = list("nreps" = nreps, "burn" = burn, "thin" = thin, "nchains" = nchains, "wd" = wd, "time_start" = beg_time), file = "output/MAGMA_arguments.txt")
  
  invisible(foreach(chain = chains, .packages = c("abind"))%dopar%{
  
    rdirich <- function(alpha0){ 
      if(sum(alpha0)){ vec = rgamma(length(alpha0), alpha0, 1) ; vec = vec / sum(vec) ; vec[vec == 0] = .Machine$double.xmin ; vec} 
      else{ rep(1, length(alpha0)) }}

    pops <- rownames(y)
  
    na_a <- which(is.na(metadat$a))

    na_i <- which(is.na(metadat$i))

    metadat$a <- factor(metadat$a, seq(C))

    metadat$i <- factor(metadat$i, seq(K + H))

    trait_fac <- factor(rep(names(nstates), nstates), levels = names(nstates))

    states <- paste0(paste0(trait_fac, "_"), unlist(lapply(nstates, seq.int)))

    alleles <- states[seq.int(sum(nalleles))]
  
    ages <- states[-seq.int(sum(nalleles))]

    beta <- matrix(0, nrow = nrow(y), ncol = ncol(y), dimnames = dimnames(y))

    beta[wildpops, alleles] <- matrix(rep(1 / nalleles, nalleles), 
                                      nrow = K, 
                                      ncol = sum(nalleles), 
                                      byrow = TRUE, 
                                      dimnames = list(wildpops, alleles))

    beta[pops, ages] <- matrix(rep(rep(1 / table(age_class), 
                                       table(age_class)) / max(age_class), T),
                               nrow = K + H, ncol = T * C, 
                               byrow = TRUE, 
                               dimnames = list(pops, ages))
  
    beta_prm <- y + beta

    beta_prm_prm <- beta_prm

    t_lnq <- log(apply(beta_prm, 1, function(rw){
        unlist(tapply(rw, trait_fac, function(alfa){ 
        if(sum(alfa)){ alfa / sum(alfa) } else{ 
          rep(0, length(alfa)) } }, simplify = FALSE)[names(nstates)])
      })) 

    freq <- matrix(0, 
                   nrow = nrow(x), 
                   ncol = K + H, 
                   dimnames = list(rownames(x), pops))
    
    if( length( hatcheries) > 0) {
      freq[-na_i, hatcheries] <- t(sapply(as.integer(metadat$i[-na_i]) - K, function(m){ 
        ans = rep(0L, H) ; ans[m] = 1L ; ans}))
    } # Chase Jalbert added if statement to deal with instances when there are no hatcheries. 
    
    freq[na_i, wildpops] <- exp(x[na_i, ] %*% t_lnq[, wildpops])

    p <- pPrior <- array(abind(lapply(seq(D), 
                                      function(d){abind(lapply(seq(max(S)), function(s){ 
                                        matrix((1 / table(groups[, d]) / max(groups[, d]))[groups[, d]], nrow = W, ncol = K + H, byrow = TRUE) }), along = 0)  }), along = 0), c(T, D, max(S), W, K+H))
    
    metadat$i[na_i] <- unlist(lapply(na_i, function(m){ 
      sample(K, 1, TRUE, p[metadat[m, 1], metadat[m, 2], metadat[m, 3], metadat[m, 4], seq(K)] * freq[m, seq(K)]) }))
        
    p <- aperm(apply(table(metadat[, -ncol(metadat)]) + pPrior, seq(4), rdirich), c(2:5, 1))
  
    for(rep in seq(nreps)){      

      metadat$a[na_a] <- unlist(lapply(na_a, function(m){
        sample(C, 1, TRUE, exp(t_lnq[-seq.int(sum(nalleles)), metadat$i[m]][seq.int(C * (metadat$t[m] - 1) + 1, C * metadat$t[m])])) }))

      metadat$i[na_i] <- unlist(lapply(na_i, function(m){
        sample(K, 1, TRUE, p[metadat$t[m], metadat$d[m], metadat$s[m], metadat$w[m], seq.int(K)] * freq[m, seq.int(K)])})) 

      p <- aperm(apply(table(metadat[, -ncol(metadat)]) + pPrior, seq(4), rdirich), c(2:5, 1))

      if(update_base & rep > nadapt){

        x_sum <- matrix(0L, nrow = nrow(y), ncol = ncol(y), dimnames = dimnames(y))

        x_sum[as.integer(sort(unique(metadat$i))), ] <- rowsum(x, group = metadat$i)
    
        beta_prm_prm <- beta_prm + x_sum
    
        t_lnq <- log(apply(beta_prm_prm, 1, function(rw){ 
          unlist(tapply(rw, INDEX = trait_fac, FUN = rdirich)) }) )
    
        freq[na_i, wildpops] <- exp(x[na_i, ] %*% t_lnq[, wildpops])

      }#update_base 
    
      if(rep > burn & ! rep %% thin){
      
                     invisible(lapply(seq(T), function(t){ 
                       lapply(seq(D), function(d){ 
                         lapply(seq(S[d]), function(s){ 
                           lapply(seq(W), function(w){ 

                       PiP = data.frame(diag(p[t, d, s, w, ]) %*% t(exp(t_lnq[seq.int(C * (t - 1) + 1, C * t) + sum(nalleles), ]))) ;
					   
			     file <- paste0(paths[[t]][[d]][[s]][[w]], "/PiP", t, d, s, w, "_", chain, ".csv") ;

                       if(rep == burn + thin){ 
                         write.table(as.list(ages), file = file, sep = ",", append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE) } ;               
 
                       write.table(format(PiP, trim = TRUE, digits = sigfigs, scientific = TRUE), file = file, sep = ",", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

                     }) }) }) }))

      }#if

    }#rep
  
  })#chain

  end_time <- Sys.time() - beg_time

  print(end_time)

  dput(x = list("nreps" = nreps, "burn" = burn, "thin" = thin, "nchains" = nchains, "wd" = wd, "time_start" = beg_time, "time_end" = end_time), file = "output/MAGMA_arguments.txt")
  
  stopCluster(cl)

  invisible(lapply(search()[!search() %in% search_begin], detach, character.only = TRUE))

  return(end_time)

}
