rm(list=ls())
library(MASS)
######### global_dnds
dnds2wgs.noncoding <- function(maf, RefElement, exclsamples, negbeta, trinucMuts, outp){

  annot <- maf
  #annot <- annot[order(annot$sampleID, annot$chr, annot$pos),]
  annot <- annot[order(annot$chr, annot$pos),]
  
  ## 3. Calculate the global dN/dS
  message("[3] Estimating global rates...")
  
  Lall <- array(sapply(RefElement, function(x) x$L), dim=c(192,2,length(RefElement)))
  Nall <- array(sapply(RefElement, function(x) x$N), dim=c(192,2,length(RefElement)))
  L <- apply(Lall, c(1,2), sum)
  N <- apply(Nall, c(1,2), sum) 
  # Subfunction: fitting substitution mode
  fit_substmodel <- function(N, L, substmodel) { 
    l <- c(L); n <- c(N); r <- c(substmodel)
    n <- n[l!=0]; r <- r[l!=0]; l <- l[l!=0]
    params <- unique(base::strsplit(x=paste(r,collapse="*"), split="\\*")[[1]])
    indmat <- as.data.frame(array(0, dim=c(length(r),length(params))))
    colnames(indmat) <- params
    for (j in 1:length(r)) {
      indmat[j, base::strsplit(r[j], split="\\*")[[1]]] = 1
    }
    
    model <- glm(formula = n ~ offset(log(l)) + . -1, data=indmat, family=poisson(link=log))
    mle <- exp(coefficients(model)) # Maximum-likelihood estimates for the rate params
    ci <- exp(confint.default(model)) # Wald confidence intervals
    par <- data.frame(name=gsub("\`","",rownames(ci)), mle=mle[rownames(ci)], cilow=ci[,1], cihigh=ci[,2])
    model_qua <- glm(formula = n ~ offset(log(l)) + . -1, data=indmat, family=quasipoisson(link=log))
    mle_qua <- exp(coefficients(model_qua)) # Maximum-likelihood estimates for the rate params
    ci_qua <- exp(confint.default(model_qua)) # Wald confidence intervals
    par_qua <- data.frame(name=gsub("\`","",rownames(ci_qua)), mle_qua=mle_qua[rownames(ci_qua)], cilow_qua=ci_qua[,1], cihigh_qua=ci_qua[,2])
    
    return(list(par=par, model=model,par_qua=par_qua,model_qua=model_qua))
  }
  
  # define substitution model
  substmodel <- array(NA,dim=c(192,2))
  rownames(substmodel) <- trinucMuts
  colnames(substmodel) <- c("neutral","selection") 
  substmodel[,"neutral"] <- c(paste0("t*",trinucMuts[-length(trinucMuts)]),"t")
  substmodel[,"selection"] <- c(paste0("t*",trinucMuts[-length(trinucMuts)],"*wsel"),"t*wsel")
  poissout <- fit_substmodel(N, L, substmodel) # Original substitution model
  par <- poissout$par
  poissmodel <- poissout$model
  parmle =  setNames(par[,2], par[,1])
  globalAIC <- poissmodel$aic
  globaldeviance <- poissmodel$deviance
  mle_submodel = par
  rownames(mle_submodel) = NULL
  
  par_qua <- poissout$par_qua
  parmle_qua <- setNames(par_qua[,2],par_qua[,1])
  
  ######overdispersion calculation
  overdisp_fun <- function(model) {
    model.df <- model$df.residual
    rdf <- model.df
    rp <- residuals(model,type="pearson")
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
    c(chisq=Pearson.chisq,ratio=prat,p=pval)
  }
  overdisp_res <- as.character(overdisp_fun(poissmodel))
  
  globaldnds <- as.data.frame(matrix(c(as.character(par["wsel",c(2:4)]),globalAIC,globaldeviance,c(overdisp_res),par_qua["wsel",c(2:4)]),nrow=1))
  colnames(globaldnds) <- c("mle","ci_low","ci_high","AIC","deviance",
                            "overdis_chisq","overdis_ratio","overdis_p",
                            "mle_qua","ci_low_qua","ci_high_qua")
  
  # all the genes name
  elementname <- data.frame(element_name = sapply(RefElement, function(x) x$element_name))
  ######
  ## 4. dNdSloc: variable rate dN/dS model (gene mutation rate inferred from neutral subs in the gene only)
  genemuts = data.frame(gene_name = sapply(RefElement, function(x) x$gene_name), n_neutral=NA, n_select=NA, exp_neutral=NA, exp_select=NA)
  genemuts[,2:3] = t(sapply(RefElement, function(x) colSums(x$N)))
  mutrates = sapply(substmodel[,1], function(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]])) # Expected rate per available site
  genemuts[,4:5] = t(sapply(RefElement, function(x) colSums(x$L*mutrates)))
  numrates = length(mutrates)
  sel_loc <- data.frame()
  sel_cv <- data.frame()
  
  if (outp == 3) {
    message("[4] Running dNdSloc...")  
    selfun_loc = function(j) {
      y = as.numeric(genemuts[j,-1])
      x = RefElement[[j]]
      
      # a. Neutral model: neutral==1, select==1
      mrfold = sum(y[1:2])/sum(y[3:4]) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
      ll0 = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,negbeta),dim=c(2,numrates))), log=T)) # loglik null model
      
      # b. Missense model: wneutural==1, free select
      mrfold = max(1e-10, sum(y[c(1)])/sum(y[c(3)])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
      wfree = y[2]/y[4]/mrfold; wfree[y[2]==0] = 0# MLE of dN/dS based on the local rate (using neutral muts as neutral)
      llsel = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,wfree),dim=c(2,numrates))), log=T)) # loglik free wmis
      wfree[wfree>1e4] = 1e4
      
      p = 1-pchisq(2*(llsel-ll0),df=c(1))
      return(c(wfree,p))
    } 
    
    sel_loc = as.data.frame(t(sapply(1:nrow(genemuts), selfun_loc)))
    colnames(sel_loc) = c("wsel_loc","psel_loc")
    sel_loc$qsel_loc = p.adjust(as.numeric(sel_loc$psel_loc), method="BH")
    sel_loc = cbind(genemuts[,1:3],sel_loc)
    sel_loc = sel_loc[order(sel_loc$psel_loc,-sel_loc$wsel_loc),]
    sel_loc = na.omit(sel_loc)
  }
  
  ############
  ## 5. dNdScv: Negative binomial regression (with or without covariates) + local synonymous mutations
  
  nbreg = nbregind = NULL
  if (outp > 1) { 
    message("[5] Running dNdScv...")  
    # Covariates   
    nbrdf = genemuts[,c("n_neutral","exp_neutral")]
    #print(head(nbrdf))
    model = glm.nb(n_neutral ~ offset(log(exp_neutral)) - 1 , data = nbrdf)
    message(sprintf("Regression model for substitutions: no covariates were used (theta = %0.3g).", model$theta))   
    if (all(model$y==genemuts$n_neutral)) {
      genemuts$exp_neutral_cv = model$fitted.values
    }
    theta = model$theta
    nbreg = model
    
    # Subfunction: Analytical opt_t using only neutral subs
    mle_tcv = function(n_neutral, exp_rel_neutral, shape, scale) {
      tml = (n_neutral+shape-1)/(exp_rel_neutral+(1/scale))
      if (shape<=1) { # i.e. when theta<=1
      tml = max(shape*scale,tml) # i.e. tml is bounded to the mean of the gamma (i.e. y[9]) when theta<=1, since otherwise it takes meaningless values
      }
      return(tml)
    }
    
    # Subfunction: dNdScv per gene
    selfun_cv = function(j) {
      y = as.numeric(genemuts[j,-1])
      x = RefElement[[j]]
      exp_rel = y[3:4]/y[3]
      # Gamma
      shape = theta
      scale = y[5]/theta
      
      # a. Neutral model
      #indneut = 1:4 # vector of neutral mutation types under this model (1=synonymous, 2=missense, 3=nonsense, 4=essential_splice)
      indneut=1:2
      opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
      mrfold = max(1e-10, opt_t/y[3]) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
      ll0 = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,negbeta),dim=c(2,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik null model
      
      # b. Select model: wneutral==1, free select
      indneut = 1
      opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
      mrfold = max(1e-10, opt_t/sum(y[3])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
      wfree = y[2]/y[4]/mrfold; wfree[y[2]==0] = 0
      llsel = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,wfree),dim=c(2,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis   
      p = 1-pchisq(2*(llsel-ll0),df=c(1))
      return(c(wfree,p))
     }
    
    sel_cv = as.data.frame(t(sapply(1:nrow(genemuts), selfun_cv)))
    colnames(sel_cv) = c("wsel_cv","psel_cv")
    sel_cv$qsel_cv = p.adjust(as.numeric(sel_cv$psel_cv, method="BH"))
    sel_cv = cbind(genemuts[,1:3],sel_cv)
    sel_cv = sel_cv[order(sel_cv$psel_cv,-sel_cv$wsel_cv),]# Sorting genes in the output file
    sel_cv = na.omit(sel_cv)
  } 
  CADD_dndsWGSout <- list(globaldnds = globaldnds, sel_cv = sel_cv, sel_loc = sel_loc,
                          annotmuts = annot, genemuts = genemuts, mle_submodel = mle_submodel,
                          exclsamples = exclsamples, exclmuts = NULL, nbreg = nbreg,
                          nbregind = nbregind, poissmodel = poissmodel)
}#EOF

