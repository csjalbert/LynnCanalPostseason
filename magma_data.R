magma_data <-
  function(wd,
           loci_name,
           wildpops_name,
           years2run,
           age_classes) {
    ## ARGS  ########################################################################################################################################################################################################################################################################################################################################################################################################################################################
    #
    #  wd <- "V:/Analysis/1_SEAK/Sockeye/Mixture/Lynn Canal Postseason/2018/MAGMA" # Working directory with a subdirectory named 'data' that contains the files named:
    #                                                 # 'baseline.RData', 'mixture.RData', 'metadata.txt', 'group_names.txt', 'groups.txt', 'harvest.txt'.
    #  loci_name <- "loci"
    #
    #  wildpops_name <- "OrderedPops"
    #
    #  years2run <- 2018
    #
    #
    #
    # 01/24/2019 - updates by Chase Jalbert to account for MAGMA used in settings other than TBR.
    # For example, TBR assumes a composite age 0X class, thus the code is hard-coded for this. In some instances, we do not
    # want a composite 0X class and may want 0's to fall in "not".
    #
    # I cleaned up a few minor formatting things and added a few if() statements as catches for:
    # no hatchery, single district, no composite 0X age class.
    #
    # Also, I removed the hard-coded age classes, so now it must be input when running "magma_data"
    # For example: age_classes = c("02", "03", "12", "13", "22", "23", "not") OR use "0X" for composite 0's
    ### Age Data ###################################################################################################################################################
    
    search_begin <- search()
    
    source("C:/R/function_gcl.r")
    
    setwd(wd)
    
    pkgs <- c("abind")
    
    require.GCL(pkgs)
    
    metadat00 <-
      read.table(
        "data/metadata.txt",
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        row.names = 1
      )
    
    euro_age <-
      setNames(metadat00$AGE_EUROPEAN, rownames(metadat00))
    
    euro_age[!is.na(euro_age)] <-
      sapply(euro_age[!is.na(euro_age)], function(chr) {
        paste(c(rep(0, 2 - nchar(chr)), chr), collapse = "")
      })
    
    fw_age_range <-
      range(as.integer(substr(
        euro_age[!is.na(euro_age)], start = 1, stop = 1
      )))
    
    sw_age_range <-
      range(as.integer(substr(
        euro_age[!is.na(euro_age)], start = 2, stop = 2
      )))
    
    fw_ages <- seq.int(fw_age_range[1], fw_age_range[2])
    
    sw_ages <- seq.int(sw_age_range[1], sw_age_range[2])
    
    euro_ages <-
      sort(apply(expand.grid(fw_ages, sw_ages), 1 , paste, collapse = ""))
    
    C <- length(euro_ages)
    
    A <- length(age_classes)
    
    age_class <- setNames(rep(A, C), euro_ages)
    
    if ('0X' %in% age_classes) {
      age_class[substring(euro_ages, 1, 1) == "0"] <- 1L
    } # Chase Jalbert added the if statement to account for instances when we do not want to hardcode age 0 composite class.
    
    age_class[euro_ages %in% age_classes] <-
      match(euro_ages[euro_ages %in% age_classes], age_classes)
    
    
    ## Mixture/Meta Data ########################################################################################################################################################################################################################################################################################################################################################################################################################################################
    
    metadat0 <- subset(metadat00, YEAR %in% years2run)
    
    attach("data/mixture.RData")
    
    mix_sillys <-
      unique(sapply(lapply(
        lapply(rownames(metadat0), strsplit, split = "_"), "[[", 1
      ), "[", 1))
    
    loci <- get(loci_name, pos = 2)
    
    nalleles <- LocusControl$nalleles[loci]
    
    invisible(PoolCollections.GCL(mix_sillys, loci, newname = "mix"))
    
    detach("file:data/mixture.RData", character.only = TRUE)
    
    invisible(RemoveIDs.GCL("mix", mix.gcl$attributes$FK_FISH_ID[!mix.gcl$attributes$SillySource %in% rownames(metadat0)]))
    
    T  <- length(years2run)
    
    year <- as.integer(factor(metadat0$YEAR, levels = years2run))
    
    age <-
      as.integer(factor(euro_age[rownames(metadat0)],  levels = euro_ages))
    
    traits <- c(loci, paste0("age", years2run))
    
    nstates <-
      c(nalleles, setNames(rep(C, T), paste0("age", years2run)))
    
    trait_fac <- factor(rep(traits, nstates), levels = traits)
    
    states_of_traits <-
      paste(trait_fac, unlist(lapply(nstates, seq)), sep = "_")
    
    x0 <- mix.gcl$counts[, loci, ]
    
    x0[is.na(x0)] <- 0L
    
    rownames(x0) <- mix.gcl$attributes$SillySource
    
    x <-
      matrix(
        0L,
        nrow = nrow(metadat0),
        ncol = sum(nstates),
        dimnames = list(rownames(metadat0), states_of_traits)
      )
    
    x[rownames(x0), seq(sum(nalleles))] <-
      Reduce(cbind, lapply(loci, function(locus) {
        x0[, locus, seq(nalleles[locus])]
      }))
    
    x[!is.na(age), seq(sum(nalleles) + 1, sum(nstates))] <-
      t(sapply(seq(sum(!is.na(
        age
      ))), function(mm) {
        ag = C * (year[!is.na(age)][mm] - 1) + age[!is.na(age)][mm]
        ags = rep(0L, T * C)
        ags[+ag] = 1L
        ags
      }))
    
    ## Baseline Data ########################################################################################################################################################################################################################################################################################################################################################################################################################################################
    
    groups0 <-
      read.table(
        "data/groups.txt",
        header = TRUE,
        row.names = 1,
        sep = "\t",
        stringsAsFactors = FALSE
      )
    
    G <- apply(groups0, 2, max)
    
    group_names0 <-
      read.table(
        "data/group_names.txt",
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      )
    
    last_group_name <-
      sapply(seq(G), function(g) {
        group_names0[G[g], g]
      })
    
    attach("data/baseline.RData")
    
    hatcheries <-
      sort(unique(metadat00$SOURCE[metadat00$SOURCE != "WILD"]))
    
    group_names <- group_names0
    
    miss_hatcheries <-
      hatcheries[!hatcheries %in% rownames(groups0)]
    
    groups <-
      rbind(groups0,
            matrix(
              rep(G, length(miss_hatcheries)),
              nrow = length(miss_hatcheries),
              ncol = length(G),
              byrow = TRUE,
              dimnames = list(miss_hatcheries, names(groups0))
            ))
    
    wildpops <- rownames(groups)[!rownames(groups) %in% hatcheries]
    
    pops <- c(wildpops, hatcheries)
    
    if (NCOL(groups) > 1) {
      groups <- groups[pops, ]
    } # Chase Jalbert added the if statement to accound for when there is only 1 district. Otherwise groups is converted to a vector, crashing downstream analyses.
    
    y0 <- FreqPop.GCL(wildpops, loci)
    
    detach("file:data/baseline.RData", character.only = TRUE)
    
    K <- length(wildpops)
    
    H <- length(hatcheries)
    
    alleles <-
      states_of_traits[-seq(sum(nalleles) + 1, sum(nstates))]
    
    y <-
      matrix(
        0L,
        nrow = K + H,
        ncol = sum(nstates),
        dimnames = list(pops, states_of_traits)
      )
    
    y[wildpops, alleles] <-
      Reduce(cbind, lapply(loci, function(locus) {
        y0[wildpops, locus, seq(nalleles[locus])]
      }))
    
    i <-
      as.integer(factor(metadat0$SOURCE, levels = hatcheries)) + K
    
    districts <-
      setNames(as.character(sort(unique(
        metadat0$DISTRICT
      ))), sort(unique(metadat0$DISTRICT)))
    
    D <- length(districts)
    
    district <-
      as.integer(factor(metadat0$DISTRICT, levels = districts))
    
    subdistricts <-
      lapply(districts, function(d) {
        setNames(sort(unique(
          subset(metadat0, DISTRICT == d)$SUBDISTRICT
        )), sort(unique(
          subset(metadat0, DISTRICT == d)$SUBDISTRICT
        )))
      })
    
    S = sapply(subdistricts, length)
    
    subdistrict <-
      as.integer(unlist(lapply(districts, function(d) {
        as.integer(factor(metadat0$SUBDISTRICT[metadat0$DISTRICT == d], levels = subdistricts[[d]]))
      })))
    
    stat_weeks <- sort(unique(metadat0$STAT_WEEK))
    
    W <- length(stat_weeks)
    
    stat_week <-
      as.integer(factor(metadat0$STAT_WEEK, levels = stat_weeks))
    
    metadat <-
      data.frame(
        t = year,
        d = district,
        s = subdistrict,
        w = stat_week,
        i = i,
        a = age,
        stringsAsFactors = FALSE,
        row.names = rownames(metadat0)
      )
    
    harvest <-
      read.table(
        "data/harvest.txt",
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      )
    
    harvest$YEAR <-
      as.integer(factor(harvest$YEAR, levels = years2run))
    
    harvest$DISTRICT <-
      as.integer(factor(harvest$DISTRICT, levels = districts))
    
    harvest$SUBDISTRICT <-
      sapply(seq(nrow(harvest)), function(rw) {
        as.integer(factor(harvest$SUBDISTRICT[rw], levels = subdistricts[[harvest$DISTRICT[rw]]]))
      })
    
    harvest$STAT_WEEK <-
      as.integer(factor(harvest$STAT_WEEK, levels = stat_weeks))
    
    harvest$HARVEST <- as.integer(harvest$HARVEST)
    
    magma_data_names <-
      c(
        "y",
        "x",
        "metadat",
        "harvest",
        "groups",
        "group_names",
        "K",
        "H",
        "A",
        "C",
        "T",
        "D",
        "S",
        "W",
        "age_class",
        "age_classes",
        "years2run",
        "districts",
        "subdistricts",
        "stat_weeks",
        "loci",
        "nalleles",
        "nstates",
        "hatcheries",
        "wildpops",
        "magma_data_names"
      )
    
    save(list = magma_data_names, file = "analysis/magma_data.RData")
    
    detach(paste0("package:", pkgs), character.only = TRUE)
    
    return(NULL)
    
  }
