#06mar2018

# makeItemParms
# Generate
# Analyze
# Summarize

# packages:
# require(lme4)
# require(CTT)
# require(dplyr)

makeItemParms <- function(startparms=NULL) {
  
  if (!is.null(startparms)) {
    
    set.seed(123)
    
    p <- read.csv(file = startparms)
    p <- p[sample(1:nrow(p), size=24, replace = FALSE),]
    p <- p[,2:3]
    p <- data.frame(item=c(1:nrow(p)),p)
    row.names(p) <- NULL
    
    # default no DIF or DIF var
    # center b values at 0 and a values at 1
    p$bdif   <- 0
    p$difvar <- 0
    p$a      <- round(p$a - mean(p$a) + 1, 3)
    p$b      <- round(p$b - mean(p$b), 3)
    
    p <- rbind(p, c(25,1,0,0,0))
    
#    write.csv(p, file="difhlm-itemparms.csv", row.names=FALSE)
    
    # list of item params for each condition
    iparms <- list()
    
    for (bval in c("low", "med", "high")) {
      for (difvar in c("no", "with")) {
        for (bdif in c("no", "with")) {
          xbval <- switch (bval,
                           low  = -0.75,
                           med  = 0.00,
                           high = 0.75
          )
          
          xdifvar <- switch (difvar,
                             no   = 0,
                             with = 0.8
          )
          
          xbdif <- switch (bdif,
                           no    = 0,
                           with  = 0.6
          )
          
          p[25,"b"]      <- xbval
          p[25,"bdif"]   <- xbdif
          p[25,"difvar"] <- xdifvar
          
          iparms[[paste(bval,bdif,difvar, sep="_")]] <- p
        }
      }
    }
    return(iparms)
  }
}

Generate <- function(condition, fixed_objects = NULL) {
  
  # icc: intraclass correlation coefficient of theta across clusters
  # J:      number of clusters
  # nj:     vector length J with n-size per cluster; if length 1 use that value 
  #         for all J clusters; if odd, still just assigns group 0/1 alternating
  #         through whole index vector when group=within;
  #         with group=between first half of all indexes are 0 then 1;
  #         this might leave out some observations at the end of total sample
  #         size odd ...
  # impact: mean difference in theta values between groups; will be set up so
  #         that mean_ref = (impact/2) and mean_focal = -(impact/2) to keep
  #         distribution centered around 0; total variance maintained at 1 by
  #         reducing the within-group-within-cluster variance
  #
  # grp:    either "within" or "between" to have the grouping variable vary 
  #         within clusters or across clusters ... how is this determined?
  # iparms: a data frame with columns "item" "a" "b" "bdif" "difvar" which are
  #         item number, discrimination, difficulty for reference group, avg.
  #         DIF value, variance of DIF value across schools
  #         
  # OPTIONAL:
  # data_list:  list of indices
  # group:      group indicator of length sum(nj)
  # rej:        cluster ability random effects
  # theta:      person within-cluster deviation ability effects
  
  # for testing
  # icc=0.2; J=2; nj=10; grp="within"; impact=0
  # iparms=data.frame(item=c(1,2,3),
  #                   a=c(.8,.9,1),
  #                   b=c(-1,0,1),
  #                   bdif=c(0,0,.2),
  #                   difvar=c(0,0,1))
  # group=NULL; rej=NULL; theta=NULL
  
  icc    <- with(condition, icc)
  J      <- with(condition, J)
  nj     <- with(condition, N)
  impact <- with(condition, impact)
  grp    <- ifelse("grp" %in% names(condition), paste(condition$grp), "within") #condition$grp, "within")
  bdif   <- condition$dif
  difvar <- condition$difvar
  bval   <- condition$bval 

  iparms <- fixed_objects[["iparms"]][[paste(bval, bdif, difvar, sep = "_")]]
  group  <- NULL
  rej    <- NULL
  theta  <- NULL
  
  if (length(nj) == 1) nj <- rep(nj, times = J)
  
  nitems  <- nrow(iparms)
  N       <- sum(nj)       # total N
  
  # cluster dif
  randif   <- numeric()
  for (i in 1:J) {
    randif <- rbind(randif, apply(iparms, 1, function(x) {
      rnorm(1, mean = 0, sd = sqrt(x["difvar"]))
    }))
  }
  
  if (is.null(group)) {
    # this creates a set of person and cluster indices for 
    # fully balanced data
    # for now, assume if data generated randomly it is all balanced group sizes
    group     <- switch(grp,
                        within  = rep(0:1, times = (N/2)),
                        between = rep(0:1, each  = (N/2)))
  }
  
  # rej: cluster ability effect; theta: student ability
  if(is.null(rej))   { rej    <- rnorm(J, mean = 0, sd = sqrt(icc)) }
  if(is.null(theta)) { 
    theta  <- rnorm(N, mean = 0,
                    sd = sqrt(1-icc-(0.5*(1-0.5)*(impact)^2))) +
      (impact*0.5) - (group*impact)
  }
  
  stopifnot(length(rej) == J & length(theta) == N & length(group) == N)
  
  # set this up as a*(theta_stud + theta_clust - b - (g * b_cluster))
  Mtheta <- matrix(rep(theta, each=nitems), nrow = N, byrow = T)
  Mrej   <- matrix(rep(rep(rej, nj), times=nitems), nrow=N, byrow=F)
  Mb     <- matrix(rep(iparms$b, each = N), nrow=N, byrow=F)
  Ma     <- matrix(rep(iparms$a, each = N), nrow=N, byrow=F)
  Mg     <- matrix(rep(group, each=nitems), nrow=N, byrow=T)
  Mbdif  <- matrix(rep(iparms$bdif,times=N), nrow=N, byrow=T)
  Mbdifj <- randif[rep(c(1:J), nj),]
  
  eta  <- Ma * (Mtheta + Mrej - Mb - (Mg * (Mbdif + Mbdifj)))
  prob <- exp(eta)/(1+exp(eta))
  y    <- prob > runif(length(prob))
  y[y==TRUE]  <- 1
  colnames(y) <- paste0("y", 1:nitems)
  
  X <- data.frame(cluster = rep(c(1:J), nj),
                  group = group,
                  y)
  
  return(list(theta2         = Mtheta[,1]+Mrej[,1], # student theta's
              cluster        = X$cluster,           # cluster ID
              cluster_effect = rej,                 # cluster effects
              cluster_dif    = randif,              # cluster DIF
              group          = X$group,             # student-level group
              resp           = X[,c(3:ncol(X))]))   # item response data
}

Analyze <- function(condition, dat, fixed_objects = NULL) {
  
  # goal of this function is to return
  # - model estimates for empty ri; ri_x_g; rc_x_g_ind; lr_x_g
  #	- total score summary stats including ICC and alpha
  #	- theta summary stats including ICC
  # - ICC of each item, using logit RI model and continuous RI model
  
  require(lme4)
  require(CTT)
  
  getModels <- function(d, g5=0, ts_type="ts", difitem) {
    # g5=0.5 to do -0.5 / 0.5 coding
    # ts_type controls matching variable; "ts"=total score
    # difitem indicates which item is being tested for DIF
    
    dat <- data.frame(y = d$resp[,paste0("y", difitem)],
                      j = d$cluster,
                      ts = scale(apply(d$resp, 1, sum)), 
                      theta = d$theta2,
                      grp = d$group-g5)
    
    fmla <- list(
      ri         = as.formula(paste0("y ~ (1|j)")),
      ri_x       = as.formula(paste0("y ~ ", ts_type, " + (1|j)")),
      ri_x_g     = as.formula(paste0("y ~ ", ts_type, "+ grp + (1|j)")),
      rc_x_g_ind = as.formula(paste0("y ~ ", ts_type,"+ grp + (1|j) + (-1+grp|j)"))
    )
    
    mods <- list()
    
    for (i in 1:length(fmla)) {
      mods[[names(fmla)[i]]] <- list(
        mod  = glmer(fmla[[i]], data = dat, family = "binomial", nAGQ = 1),
        name = names(fmla)[i]
      )
    }
    
    mods[["lr"]] <- list(
      mod  = glm(as.formula(paste0("y ~ ", ts_type, " + grp")), data = dat,
                 family = "binomial"),
      name = "lr"
    )
    
    lr_res <- data.frame()
    
    get_summary <- function(obj) {
      m    <- obj$mod
      nm   <- obj$name
      conv <- getConv(m)
      coef <- data.frame(summary(m)$coef, model=nm, conv=conv,
                         ts_type=ts_type, LL=logLik(m))
      coef$parm <- rownames(coef)
      rownames(coef) <- NULL
      if(class(m)[1] == "glmerMod"){
        vcov <- data.frame(VarCorr(m))
        vcov <- data.frame(vcov, model=nm, conv=conv, ts_type=ts_type,
                           LL=logLik(m))
        rownames(vcov) <- NULL
      } else {
        vcov <- NULL
      }
      return(list(
        coef = coef,
        vcov = vcov
      ))
    }
    
    # flag failed convergence to re-run reps
    # indicate which models didn't converge
    if(!getConv(mods[["rc_x_g_ind"]]$mod)) {
      if(!getConv(mods[["ri_x_g"]]$mod)) {
        stop(paste0("rc and ri failed to converge it", difitem))
      } else {
        stop(paste0("rc failed to converge it", difitem))
      }
    }
    if(!getConv(mods[["ri_x_g"]]$mod)) stop(paste0("ri failed to converge it", difitem))
    
    return(list(
      feff = do.call(rbind, lapply(mods, function(x) {get_summary(x)$coef})),
      reff = do.call(rbind, lapply(mods, function(x) {get_summary(x)$vcov})),
      item = difitem
    ))
    
  }
  
  getConv <- function(m) {
    if(class(m)[1] %in% c("glmerMod", "lmerMod")){
      conv <- !any(grepl("failed to converge", m@optinfo$conv$lme4$messages))
    } else {
      conv <- m$converged
    }
    return(conv)
  }
  
  getICC <- function(y, clust, type=1) {
    
    if (type==1) {
      m <- lmer(y~1+(1|clust))
      vcov <- data.frame(VarCorr(m))
      icc <- vcov[1,4]/(vcov[1,4]+vcov[2,4])
    }
    if (type==2) {
      m <- glmer(y~1+(1|clust), family="binomial", nAGQ=1)
      vcov <- data.frame(VarCorr(m))
      icc <- vcov[1,4]/(vcov[1,4]+((pi^2)/3))
    }
    return(c(
      est_icc = ifelse(getConv(m), icc, NA),
      conv = getConv(m))
    )
    
  }
  
  # get results
  
  fititem <- c(fixed_objects[["fititem"]])
  
  difmods <- list()
  
  for (it in fititem) {
    difmods[[paste0("it",it)]] <- getModels(d=dat, g5=0, ts_type="ts", difitem=it)
  }
  
  # ICC's across all items
  
  item_icc_linear <- data.frame(
    t(apply(dat$resp, 2, getICC, clust=dat$cluster, type=1)),
    item=names(dat$resp)
  )
  
  item_icc_logit <- data.frame(
    t(apply(dat$resp, 2, getICC, clust=dat$cluster, type=2)),
    item=names(dat$resp)
  )
  
  ts <- apply(dat$resp, 1, sum)
  
  ts_summary <- data.frame(
    icc   = getICC(y=ts, clust=dat$cluster, type=1)[1],
    alpha = CTT::itemAnalysis(dat$resp)$alpha,
    mean  = mean(ts),
    sd    = sd(ts),
    mean_focal = mean(ts[dat$group==1]),
    sd_focal   = sd(ts[dat$group==1]),
    mean_ref   = mean(ts[dat$group==0]),
    sd_ref     = sd(ts[dat$group==0])
  )
  
  # cleaner results
  
  get_clean_res <- function(X) {
    names(X$feff) <- c("est", "se", "z", "p", "model", "conv", "ts_type", "LL", "parm")
    fe <- with(X$feff,{
      data.frame(
        lr_int    = est[model=="lr" & parm == "(Intercept)"],
        lr_x      = est[model=="lr" & parm == "ts"],
        lr_dif    = est[model=="lr" & parm == "grp"],
        lr_dif_p  =   p[model=="lr" & parm == "grp"],
        lr_dif_se =  se[model=="lr" & parm == "grp"],
        
        ri_int    = est[model=="ri" & parm == "(Intercept)"],
        
        ri_x_int  = est[model=="ri_x" & parm == "(Intercept)"],
        ri_x_x    = est[model=="ri_x" & parm == "ts"],
        ri_x_x_se =  se[model=="ri_x" & parm == "ts"],
        
        ri_x_g_int    = est[model=="ri_x_g" & parm == "(Intercept)"],
        ri_x_g_x      = est[model=="ri_x_g" & parm == "ts"],
        ri_x_g_x_se   =  se[model=="ri_x_g" & parm == "ts"],
        ri_x_g_dif    = est[model=="ri_x_g" & parm == "grp"],
        ri_x_g_dif_se =  se[model=="ri_x_g" & parm == "grp"],
        ri_x_g_dif_p  =   p[model=="ri_x_g" & parm == "grp"],
        
        rc_x_g_int    = est[model=="rc_x_g_ind" & parm == "(Intercept)"],
        rc_x_g_x      = est[model=="rc_x_g_ind" & parm == "ts"],
        rc_x_g_x_se   =  se[model=="rc_x_g_ind" & parm == "ts"],
        rc_x_g_dif    = est[model=="rc_x_g_ind" & parm == "grp"],
        rc_x_g_dif_se =  se[model=="rc_x_g_ind" & parm == "grp"],
        rc_x_g_dif_p  =   p[model=="rc_x_g_ind" & parm == "grp"],
        
        lr_LL       =   LL[model=="lr" & parm == "(Intercept)"],
        lr_conv     = conv[model=="lr" & parm == "(Intercept)"],
        ri_LL       =   LL[model=="ri" & parm == "(Intercept)"],
        ri_conv     = conv[model=="ri" & parm == "(Intercept)"],
        ri_x_LL     =   LL[model=="ri_x" & parm == "(Intercept)"],
        ri_x_conv   = conv[model=="ri_x" & parm == "(Intercept)"],
        ri_x_g_LL   =   LL[model=="ri_x_g" & parm == "(Intercept)"],
        ri_x_g_conv = conv[model=="ri_x_g" & parm == "(Intercept)"],
        rc_x_g_LL   =   LL[model=="rc_x_g_ind" & parm == "(Intercept)"],
        rc_x_g_conv = conv[model=="rc_x_g_ind" & parm == "(Intercept)"]
        
      )
    })
    
    re <- with(X$reff, {
      data.frame(
        ri_int_var     = vcov[model=="ri"],
        ri_x_int_var   = vcov[model=="ri_x"],
        ri_x_g_int_var = vcov[model=="ri_x_g"],
        rc_x_g_int_var = vcov[model=="rc_x_g_ind" & var1=="(Intercept)"],
        rc_x_g_dif_var = vcov[model=="rc_x_g_ind" & var1=="grp"]
      )
    })
    
    return(cbind(item=X$item,re,fe))
  }
  
  clean_res <- do.call(rbind, lapply(difmods, get_clean_res))
  
  # return results
  
  ret <- list(
    item_icc_linear = item_icc_linear,
    item_icc_logit  = item_icc_logit,
    ts_summary      = ts_summary,
    res             = difmods,
    clean_res       = clean_res
  )
  
  ret
  
}

Summarize <- function(condition, results, fixed_objects = NULL) {
  # function should return one row with following outcomes
  #	- 
  #	avg DIF for RI; RC; LR
  #	significance rate for RI/RC/LR DIF
  #	avg DIF variance and DIF variance significance rate
  #	proportion of converged replications for each model (should be 1 based on above)
  #	do all calculations only on RC converged replications;
  #	do all calculations by collapsing:
  #		drop reps where RI!=conv
  #		use RI when RC!=conv
  #		call these "alt" outcomes
  
  require(dplyr)

  # var components mixture chi-square p-values
  pmix <- function(ll1, ll2, df1, df2) {
    return(
      as.numeric(0.5*pchisq(2*(ll2-ll1), df=df1, lower.tail=F) +
                   0.5*pchisq(2*(ll2-ll1), df=df2, lower.tail=F)))
  }
  
  # add additional outcomes
  
  clean_res <- do.call(rbind, lapply(results, function(x){x$clean_res}))
  
  clean_res <- within(clean_res, {
    ri_icc             = ri_int_var/(ri_int_var+((pi^2)/3))
    ri_x_g_icc         = ri_x_g_int_var/(ri_x_g_int_var+((pi^2)/3))
    ri_x_g_int_var_p   = pmix(lr_LL, ri_x_g_LL, 0, 1)
    rc_x_g_dif_var_p   = pmix(ri_x_g_LL, rc_x_g_LL, 0, 1)
    lr_dif_sig         = lr_dif_p < 0.05
    ri_x_g_dif_sig     = ri_x_g_dif_p < 0.05
    rc_x_g_dif_sig     = rc_x_g_dif_p < 0.05
    rc_x_g_dif_var_sig = rc_x_g_dif_var_p < 0.05
  })
  
  # average of all item ICCs, drop non-converged models
  
  all_icc_logits <- do.call(rbind, lapply(results, function(x) {
    x$item_icc_logit
  }))
  all_icc_linear <- do.call(rbind, lapply(results, function(x) {
    x$item_icc_linear
  }))
  
  icc_logits      <- mean(all_icc_logits[,"est_icc"], na.rm=TRUE)
  icc_logits_conv <- mean(all_icc_logits[,"conv"], na.rm=TRUE)
  icc_linear      <- mean(all_icc_linear[,"est_icc"], na.rm=TRUE)
  icc_linear_conv <- mean(all_icc_linear[,"conv"], na.rm=TRUE)
  
  # total score summary stats
  
  ts_summary <- do.call(rbind, lapply(results, function(x) {
    x$ts_summary
  }))
  ts_summary <- apply(ts_summary, 2, mean, na.rm=TRUE)
  
  # main results summaries
  # compute SD's to compare empirical to estimated SE's for DIF coef's
  
  clean_res %>% group_by(item) %>% summarise_all(funs(mean)) -> clean_res_mn
  clean_res %>% group_by(item) %>% summarise_all(funs(sd))   -> clean_res_sd
  
  # return results
  
  biglist <- list(
    mn              = clean_res_mn,
    sd              = clean_res_sd,
    icc_logits      = icc_logits,
    icc_logits_conv = icc_logits_conv,
    icc_linear      = icc_linear,
    icc_linear_conv = icc_linear_conv,
    ts              = ts_summary
  )
  
  do.call(c, biglist)
  
}



