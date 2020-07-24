########### Main Results of the Paper 
    
  rm(list=ls())
  
  library(tidyverse)
  library(AER)
  library(DBI)
  library(gdata)
  library(hdm)
  library(spdep)
  library(randomForest)
  library(nnet)
  library(glmnet)
  library(lfe)
  library(mgcv)
  library(mgcViz)
  library(doSNOW)
  library(parallel)
  library(doParallel)
  library(gbm)
  library(ggplot2)
  library(DHARMa)
  library(stargazer)
  library(ggmap)
  
  
  setwd("/Users/baka/academica/Burkina+Schisto/data_codes/replication")
  source("Full_DML.R")
  source("ML_algorithms.R")
  source("DML_Forest_IV.R")

  
  # cl   <- makeCluster(12, outfile="")
  # registerDoSNOW(cl)
  
  
  ############# dataset for 2009 and 2011, to replicate the main results of the paper. 
  ############# cleaned dataset with only repeated households
  
  datause <- read.csv("data_replication.csv") 
  datause$rend <- datause$rend11 + datause$rend12
  
  
  res_sc        <- matrix(0,2,6)
  colnames(res_sc) <- c("Estimate","S.E.","p-val","Mean Effect (%)","Top 5% Effect (%)","MSE")
  
  
  ######### invert prevalence to get schisto intensity
  

  lin_K      <- function(par, mu) par[1] + par[2] * mu
  invertPrev <- function(mu, prev, par, k) 1 - prev  - (1 + mu/k(par, mu))^-(k(par, mu))
  
  
  #### intensity parameter estimates from negative binomial fitting (see SI)
  
  est_par <- c(0.01014, 0.00112)
  datause$intensity <- matrix(0,dim(datause)[1],1)
  
  for (i in 1:length(datause$prev_2015map)){
    datause$intensity[i] <- uniroot(invertPrev, c(0,1000), prev =datause$prev_2015map[i] , par = est_par,k=lin_K)$root
  }
  
  
  ############ set working variables
  
  
  attach(datause) 
  
  
  datause$cult1 <- as.factor(datause$cult1)
  # fix some issues with cult2: obs have issues i.e. 107 instead of 1070
  datause$cult2 <- as.factor(ifelse(nchar(as.character(datause$cult2)) == 3, 
                                    paste(datause$cult2, 0, sep = ""), datause$cult2))
  # add 1000 as factor level for NA (no second cultivation in a crop)
  datause$cult2  <- factor(datause$cult2, levels = levels(addNA(datause$cult2)), 
                           labels = c(levels(datause$cult2), 1000), exclude = NULL)
  datause$prev   <- datause$prev_2015map
  datause$prevsq <- datause$prev^2
  datause$prev_poly1 <- poly(datause$prev,3,raw=TRUE)[,1]
  datause$prev_poly2 <- poly(datause$prev,3,raw=TRUE)[,2]
  datause$prev_poly3 <- poly(datause$prev,3,raw=TRUE)[,3]
  datause$intensitysq <- datause$intensity^2
  datause$intensity_poly1 <- poly(datause$intensity,3,raw=TRUE)[,1]
  datause$intensity_poly2 <- poly(datause$intensity,3,raw=TRUE)[,2]
  datause$intensity_poly3 <- poly(datause$intensity,3,raw=TRUE)[,3]
  
  
  # schistosomiasis variables: prevalence and intensity
  p        <- c("prev")
  p.sq     <- c("prev", "prevsq")
  p.flex   <- c("prev_poly1","prev_poly2","prev_poly3")
  
  int      <- c("intensity")
  int.sq   <- c("intensity","intensitysq")
  int.flex <- c("intensity_poly1","intensity_poly2","intensity_poly3")
  lint     <- c("log(intensity)")
  
  
  
  
  # years dummies (use only for full data analysis + lewbel IV)
  
  years <- as.data.frame(model.matrix(~factor(annee)+0))
  years[,dim(years)[2]] <- NULL
  colnames(years) <- paste("factor_years_",rep(1:dim(years)[2]),sep="")
  # factor dummies + drop last one of cult1 and first of cult2 (no 2nd cultivation)
  # necessary for the machine learning estimation, as it does not deal well with declaring factors
  

  
  cult_d1 <- as.data.frame(model.matrix(~factor(cult1)+0))
  cult_d1[,dim(cult_d1)[2]] <- NULL
  colnames(cult_d1) <- paste("factor_cult1_",rep(1:dim(cult_d1)[2]),sep="")
  cult_d2 <- as.data.frame(model.matrix(~factor(cult2)+0))
  cult_d2[,dim(cult_d2)[2]] <- NULL
  colnames(cult_d2) <- paste("factor_cult2_",rep(1:dim(cult_d2)[2]),sep="")
  cult_d <- cbind(cult_d1,cult_d2)
  # plot localization dummies (bush/camp) (drop dummy for NA+"cases")
  localisation_d <- as.data.frame(model.matrix(~factor(localisation)+0))
  localisation_d[,1] <- NULL
  colnames(localisation_d) <- paste("factor_localisation_",rep(1:dim(localisation_d)[2]),sep="")
  # antierosif (remove "aucun")
  antierosif_d <- as.data.frame(model.matrix(~factor(antierosif)+0))
  antierosif_d[,1] <- NULL
  colnames(antierosif_d) <- paste("factor_antierosif_",rep(1:dim(antierosif_d)[2]),sep="")
  # plot surface dummies (flatland/shallows/hillside - drop dummy for NA+flatland)
  relief_d <- as.data.frame(model.matrix(~factor(relief)+0))
  relief_d[,1:2] <- NULL
  colnames(relief_d) <- paste("factor_relief_",rep(1:dim(relief_d)[2]),sep="")
  # land acquisition dummies (bought/rented/gifted....)
  mode_acq_d <- as.data.frame(model.matrix(~factor(mode_acq)+0))
  mode_acq_d[,1:2] <- NULL
  colnames(mode_acq_d) <- paste("factor_mode_acq_",rep(1:dim(mode_acq_d)[2]),sep="")
  # land tenure dummies (land title/lease/APF/usage permit/acte de cession/owner/...)
  niv_acq_d <- as.data.frame(model.matrix(~factor(niv_acq)+0))
  niv_acq_d[,1:2] <- NULL
  colnames(niv_acq_d) <- paste("factor_niv_acq_",rep(1:dim(niv_acq_d)[2]),sep="")
  # contract regime dummies (owner operator/sharecropping/rental)
  # regime_contract_d <- as.data.frame(model.matrix(~factor(regime_contract)+0))
  # regime_contract_d[,1:2] <- NULL
  # colnames(regime_contract_d) <- paste("factor_regime_contract_",rep(1:dim(regime_contract_d)[2]),sep="")
  # labour type dummies (no labour/manual/plowing/motorised)
  type_labour_d <- as.data.frame(model.matrix(~factor(type_labour)+0))
  type_labour_d[,1] <- NULL
  colnames(type_labour_d) <- paste("factor_type_labour_",rep(1:dim(type_labour_d)[2]),sep="")
  # household education level dummies 
  educ_max_d <- as.data.frame(model.matrix(~factor(educ_max)+0))
  educ_max_d[,1:2] <- NULL
  colnames(educ_max_d) <- paste("factor_educ_max_",rep(1:dim(educ_max_d)[2]),sep="")
  # loss dummies
  perte_d <- as.data.frame(model.matrix(~factor(pertcult1)+0))
  perte_d[,1] <- NULL
  colnames(perte_d) <- paste("factor_perte_",rep(1:dim(perte_d)[2]),sep="")
  
  
  
  
  # datause <- cbind(datause,years,cult_d,localisation_d,relief_d,antierosif_d,mode_acq_d,niv_acq_d,regime_contract_d,type_labour_d,educ_max_d,perte_d)
  datause <- cbind(datause,years,cult_d,localisation_d,relief_d,antierosif_d,mode_acq_d,niv_acq_d,type_labour_d,educ_max_d,perte_d)
  # datause <- cbind(datause,perte_d)
  
  
  detach(datause)
  attach(datause)
  
  
  # plot/livestock/input controls
  # 
  X1   <- c("lsurf","lsurfsq","remune","entraid","tot_anim_traction","males_elev","fem_elev","totanim","morts","naissances",
            "bovins", "equins","ovins","caprins","dons","vols","autocons","ventebovins","achatbovins","venteovins","achatovins" ,
            "ventecaprins","achatcaprins","gestion","parc_recup","mode_labour","duree_jachere",
            "npk_kg","uree_kg", "phosph_kg",
            "pest_solide_g","pest_liquide_cl","herbicide_g" ,
            "herbicide_cl","fungicide_g","fungicide_cl","rodenticide_g","rodenticide_cl" )
  
  
  # household controls
  X2     <- c("nmem","agemoyen","age_chefmen","numfamilies","nchildren","nyoung","nwomen","culture_pluviale",
          "culture_maraichere" ,"arboriculture","autre_contresaison",
          "peche")
  
  
  # factor dummies - remove factors with no elements (only make estimation noisier)
  # use the second one for lewbel est (time fe)
  X3_t_pre   <- paste(grep("factor_", names(datause), value=TRUE, fixed=TRUE),sep="")
  X3_t       <- X3_t_pre[-c(93,94,96,98,100,101)]
  X3         <- X3_t[-1]
  
  # disease and climate controls
  
  X4      <- c("temp_mean","temp_day","temp_night",
             # "precip_sum",
             "precip_wet",
             "temp_dry",
             "dspell_mean","dspell_max",
             "lstD_yr_mean","lstN_yr_mean",
             "EVI_yr_mean" ,
             "NDVI_yr_mean" 
  )
  
  
  dis        <- "malaria"
  
  
  # flexible controls
  x.flex <-  c("poly(lsurf,3)","remune","poly(dons,3)",
               "poly(totanim,3)","entraid" ,"poly(tot_anim_traction,3)" ,"poly(males_elev,3)","poly(fem_elev,3)",
               "poly(morts,6)" , "naissances","bovins" ,"equins","caprins", "dons","vols",
               "autocons","poly(totanim,3)","poly(nwomen,3)","gestion", "parc_recup","mode_labour","poly(duree_jachere,3)",
               "poly(nmem,3)","poly(agemoyen,3)","poly(agemoyen,3)","age_chefmen","poly(numfamilies,3)",
               X3, paste("poly(",X4,",3)",sep=""),dis)
  
  
  
  
  
  # inverse hyperbolic sine transformation
  
  ihs            <- function(x) log(x + (x^2 + 1) ^ 0.5)
  datause$lrend  <-  ihs(datause$rend)
  datause$lsurf  <-  ihs(datause$superficie_1)
  datause$lsurfsq <- datause$lsurf^2
  datause$lsurfspl <- poly(datause$lsurf, 6)
  datause$poidscult <- datause$poidscult1 + datause$poidscult2
  datause$gestion <- ifelse(datause$gestion == 1 | is.na(datause$gestion), 0, 1)
  
  
  transf   <- c("tot_anim_traction","males_elev","fem_elev","totanim","morts","naissances",
            "bovins", "equins","ovins","caprins","dons","vols","autocons","ventebovins","achatbovins","venteovins","achatovins" ,
            "ventecaprins","achatcaprins","duree_jachere",
            "npk_kg","uree_kg", "phosph_kg",
            "pest_solide_g","pest_liquide_cl","herbicide_g" ,         
            "herbicide_cl","fungicide_g","fungicide_cl","rodenticide_g","rodenticide_cl" )
  
  
  
  dat_ihs           <- apply(datause[,c(transf,X4)],2,ihs)
  colnames(dat_ihs) <- paste("l_",colnames(dat_ihs),sep="")
  
  datause           <- cbind(datause,dat_ihs)
  
  
  # outcome: log yield
  y        <- "lrend"
  
  # log controls
  X1_log   <- c(c("lsurf","lsurfsq","remune","entraid","gestion","parc_recup","mode_labour" ),colnames(dat_ihs)[1:31])
  X4_log   <- colnames(dat_ihs)[32:42]

 
  x_t      <- c(X1, X2, X3_t, X4, dis)
  x        <- x_t[-which(x_t=="factor_years_1")]
  x_log    <- c(X1_log,X2,X3_t,X4_log,dis)
  
  
  detach(datause)
  attach(datause)
  
  
  
  # if trim households with very high intensity:
  
  # datause <- datause[-c(which(datause$intensity > quantile(datause$intensity, 0.995))),]
  # 
  # detach(datause)
  # attach(datause)
  
  
  
  
  ########### snails as IV  
  
  
  
  ##### prepare snails dataset
  
  snails_pond      <- read.csv("pond_means.csv")
  data_snails_pond <- left_join(datause,snails_pond,by=c("vid"))
  snails_river     <- read.csv("river_yearly_means.csv")
  snails2          <- snails_river %>% filter(habid=="perm_river_biomphalaria" | habid=="eph_river_bulinus")
  snails_pond_biom <- left_join(snails_pond, snails2, by=c("vid"))
  
  snails_pond_biom[is.na(snails_pond_biom)] <- 0
  data_snails_pond_biom <- left_join(datause,snails_pond_biom,by=c("vid"))
  

  ##### instruments 
  
  instr <- c("rainy","dry","mean")


  form_1stage_tfe <-  as.formula(
    paste(paste(int,"~ ",paste(paste(instr,collapse="+"),paste(x_log,collapse="+"),sep="+")),sep=""))
  
  form_2stage_tfe <-  as.formula(
    paste(paste(y,"~ ",paste("stage1_tfe",paste(x_log,collapse="+"),sep="+")),sep=""))
  

  form_2stage_tfe_gam <-  as.formula(
    paste(paste(y,"~ ",paste("s(stage1_tfe)",paste(x_log,collapse="+"),sep="+")),sep=""))
  
  data_snails_pond_biom$stage1_tfe <- fitted(lm(form_1stage_tfe,data=data_snails_pond_biom))
  # data_snails_pond2$stage1_tfe <- fitted(lm(form_1stage_tfe,data=data_snails_pond2))
  second_tfe            <- lm(form_2stage_tfe,data=data_snails_pond_biom)
    
  

  # cluster-bootstrapped errors. takes time - use the commented version for quicker results (very similar)
  set.seed(42)
  est_iv  <- coeftest(second_tfe,vcov=vcovBS,cluster = ~ vid,type="wild",R=3000)
  # est_iv <- coeftest(second_tfe,vcov=vcovCL,cluster = ~ vid)
   
  
  #### check that result holds with felm package (can't bootstrap but for validity of coefficients) 
  # formiv <- as.formula(paste(y,"~",paste(x_log,collapse="+"),"|annee |
  #                             (intensity~rainy+dry+mean)|vid"))
  # feiv <- felm(formiv,data=data_snails_pond_biom)

  
  ##### also check endogeneity with control function

  control_1 <-  as.formula(
    paste(paste(int,"~ ",paste(paste(instr,collapse="+"),paste(x_log,collapse="+"),sep="+")),sep=""))
  control_2 <-  as.formula(
    paste(paste(y,"~ ",paste("intensity+control_f",paste(x_log,collapse="+"),sep="+")),sep=""))
  data_snails_pond_biom$control_f <- lm(control_1,data=data_snails_pond_biom)$residuals
  control_fun     <- lm(control_2,data=data_snails_pond_biom)
  est_control     <- coeftest(control_fun,vcov=vcovCL,cluster = ~ vid)
  # yep
  
   

  
  res_sc[1,1:3] <- est_iv[2,c(1,2,4)]
  res_sc[1,4]   <- res_sc[1,1]*mean(intensity)
  res_sc[1,5]   <- res_sc[1,1]*quantile(intensity,0.95) 
  res_sc[1,6]   <- error(fitted(second_tfe),data_snails_pond$lrend)$err

  
  
  ##### non-parametric IV as control function
  
first_npiv_f  <- as.formula(paste(int,"~ s(rainy,dry,mean)+",paste(x_log,collapse="+")))
second_npiv_f <- as.formula(paste(y,"~ s(intensity,k=4)+npiv_resid+",paste(x_log,collapse="+")))

first_npiv  <- gam(first_npiv_f,data=data_snails_pond_biom)
npiv_resid  <- first_npiv$residuals
second_npiv <- gam(second_npiv_f,data=data_snails_pond_biom,gamma=20)

est_np_int_iv   <- getViz(second_npiv)
pl1_int_iv      <- plot( sm(est_np_int_iv,1),allTerms = FALSE) 
pdf(paste("schisto_intensity_iv_gam_",paste(unique(annee),collapse="_"),".pdf",sep=""))
pl1_int_iv + l_fitLine(colour = "blue") +l_ciLine(mul = 1.98, colour = "black", linetype = 2) +
  l_points(shape = 19, size = 1.1, alpha = 0.2) + xlim(0,100)+ scale_y_continuous(limits=c(-1,0.3))+
  labs(x="Intensity",y="Log Yield Variation"
       # title=paste("Disease Intensity on log yield/ha, ",paste(unique(annee),collapse="+"),sep="") 
  ) +
  theme_bw()+
  theme(axis.text.x= element_text( size=20),axis.text.y= element_text( size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),panel.border=element_blank(),
        plot.background = element_blank(),
        panel.grid= element_blank(),axis.line = element_line(color = 'black')
  )
dev.off()







######## double machine learning analysis (double Neyman orthogonal moment conditions)
  
  ######### see SI code for full ML analysis (comparison of methods, estimation without IV, et cetera)
  
  ######### ML + IV 
  
  
  # put tune=TRUE at your peril (takes around 6 hrs) 
  
  est_ML_base_int_mid_IV <- DML_Forest_IV(y,d=int,x_t,instr=instr,tune=F,data=data_snails_pond_biom,boot=T,id="annee",cl="vid")
  # est_ML_base_int_mid_IV <- DML_Forest_IV(y,d=int,x_t,instr=instr,tune=F,data=data_snails_pond,boot=T,id="annee",cl="vid")
  
  
  res_sc[2,1] <- est_ML_base_int_mid_IV$theta_2sls
  res_sc[2,2] <- est_ML_base_int_mid_IV$theta_2sls.se
  res_sc[2,3] <- est_ML_base_int_mid_IV$theta_2sls.p
  res_sc[2,4]   <- est_ML_base_int_mid_IV$theta_2sls*mean(intensity)
  res_sc[2,5]   <- est_ML_base_int_mid_IV$theta_2sls*quantile(intensity,0.95)
  res_sc[2,6]   <- est_ML_base_int_mid_IV$mse.y

  detach(datause)
  
  



  ##### wrap up results + plot 
  
  
  rownames(res_sc) <- c("Plot level","Plot level,\n ML-flexible")
  
  res_sc               <- round(res_sc,4)
  res_sc[,4:5]         <- 100*res_sc[,4:5]
  res_sc[,4:5]         <- round(res_sc[,4:5],2)
  
 
  
  ##### plot results
  res_inv  <- apply(res_sc,2,rev)
 
  
  cint_l   <- res_inv[,1] - 1.96*res_inv[,2]
  cint_h   <- res_inv[,1] + 1.96*res_inv[,2]
  tot_eff  <- paste(as.character(res_inv[,4])," (",as.character(res_inv[,5]),")",sep="")
  
  res_plot           <- cbind(res_inv,cint_l,cint_h)
  ticks              <- rev(c("Plot level", "Plot level,\n ML-flexible"))
  plot               <- broom::tidy(res_plot)
  

  ggplot(plot,aes(Estimate,factor(ticks,levels=ticks),stat="identity")) + 
    geom_point(size=1.9) +
    geom_errorbarh(aes(xmax = cint_l, xmin = cint_h),height=0.1,linetype=2,alpha=09,colour="black",size=0.3)+
    geom_vline(xintercept=0,linetype="solid",colour="red")+
    geom_label(aes(label=tot_eff,color=p.val),size=3)+
    xlab("Estimate") + 
    ylab("")+
    xlim(min(cint_l),max(cint_h)+0.002)  +
    theme(axis.text.y = element_text(lineheight=1,face="bold", color="black",size=, angle=0),
          legend.key.size = unit(0.3, "cm"),panel.border = element_rect(colour = "black", fill=NA, size=5))+
    scale_y_discrete(expand=c(0.1,0))+
    #uncomment for full results matrix 
    scale_color_gradient(name="p.val",low="darkblue",high="blue",space="Lab")+
    theme_classic2()
  # ggsave("est_full.pdf", width = 16, height = 9, units="cm")
  ggsave("est_full_small.pdf", height = 5, units="cm")
  # dev.off()

  
  
  
  
  ################## schistosomiasis and poverty 
  
  
  
  form_1stage_poverty <-  as.formula(
    paste(paste(int,"~ ",paste(paste(instr,collapse="+"),
                               paste(instr,":poverty",collapse="+"),
                               "poverty",paste(x_t,collapse="+"),sep="+")),sep=""))
  
  form_2stage_poverty <-  as.formula(
    paste(paste(y,"~ ",paste("stage1_poverty + stage1_poverty:poverty+poverty+",
                             paste(x_t,collapse="+"),sep="+")),sep=""))
  
  
  
  data_snails_pond_biom$stage1_poverty <- fitted(lm(form_1stage_poverty ,data=data_snails_pond_biom))
  second_poverty                      <- lm(form_2stage_poverty,data=data_snails_pond_biom)
  est_iv_poverty                      <- coeftest(second_poverty,vcov=vcovCL,cluster = ~ vid)
  est_iv_poverty 
  
  
  
  ## control function

  
  lrend_q  <- seq(0.2,0.7,0.05)
  lsurf_q  <- seq(0.2,0.7,0.05)
  
  
  res_pov      <- matrix(0,length(lrend_q),length(lsurf_q))
  res_pov_se   <- matrix(0,length(lrend_q),length(lsurf_q))
  res_pov_pv   <- matrix(0,length(lrend_q),length(lsurf_q))
  res_nopov    <- matrix(0,length(lrend_q),length(lsurf_q))
  res_nopov_pv <- matrix(0,length(lrend_q),length(lsurf_q))
  res_pov_eff  <- matrix(0,length(lrend_q),length(lsurf_q))
  res_pov_eff_se  <- matrix(0,length(lrend_q),length(lsurf_q))
  
  
  
  form_control_1_poverty       <-  as.formula(paste(paste(int,"~ ",paste(paste(instr,collapse="+"),
                                                                         "poverty",paste(x_t,collapse="+"),sep="+")),sep=""))
  form_control_2_poverty       <- as.formula(paste(paste(y,"~ ",paste("intensity + intensity:poverty+poverty+ control_f_poverty",
                                                                      paste(x_t,collapse="+"),sep="+")),sep=""))
  
  
  for (i in 1:length(lrend_q)){
    for (j in 1:length(lsurf_q)){
      data_snails_pond_biom$poverty <- ifelse(datause$lrend < quantile(datause$lrend,lrend_q[i]) &
                                                datause$lsurf < quantile(datause$lsurf,lsurf_q[j]),1,0)
      
      data_snails_pond_biom$control_f_poverty <- lm(form_control_1_poverty,data=data_snails_pond_biom)$residuals
      second_control_poverty                  <- lm(form_control_2_poverty,data=data_snails_pond_biom)
      est_iv_control_poverty                  <- coeftest(second_control_poverty,vcov=vcovCL,cluster = ~ vid)
      # very long computation time, select for bootstrapped errors
      # set.seed(42)
      # est_iv_control_poverty                  <- coeftest(second_control_poverty,vcov=vcovBS,cluster = ~ vid,type="wild",R=3000)
      res_pov[i,j]                            <- est_iv_control_poverty[dim(est_iv_control_poverty)[1],1]
      res_pov_se[i,j]                         <- est_iv_control_poverty[dim(est_iv_control_poverty)[1],2]
      res_pov_pv[i,j]                         <- est_iv_control_poverty[dim(est_iv_control_poverty)[1],4]
      res_nopov[i,j]                          <- est_iv_control_poverty[2,1]
      res_nopov_pv[i,j]                       <- est_iv_control_poverty[2,4]
      res_pov_eff[i,j]                        <- res_pov[i,j]*mean(datause$intensity[which(datause$lrend < quantile(datause$lrend,lrend_q[i]) & datause$lsurf < quantile(datause$lsurf,lsurf_q[j]))])
      res_pov_eff_se[i,j]                     <- res_pov_se[i,j]*mean(datause$intensity[which(datause$lrend < quantile(datause$lrend,lrend_q[i]) & datause$lsurf < quantile(datause$lsurf,lsurf_q[j]))])
      
    }
  }
  
  colnames(res_pov) <- lrend_q
  rownames(res_pov) <- lsurf_q
  colnames(res_pov_pv) <- lrend_q
  rownames(res_pov_pv) <- lsurf_q
  
  res_ci_u <- res_pov + 1.96*res_pov_se
  res_ci_l <- res_pov - 1.96*res_pov_se
  
  res_pov_eff    <- round(res_pov_eff*100,2)
  res_pov_eff_se <- round(res_pov_eff_se*100,2)
  res_ci_u_eff   <- res_pov_eff + 1.96*res_pov_eff_se
  res_ci_l_eff   <- res_pov_eff - 1.96*res_pov_eff_se
  
  # data_snails_pond_biom$control_f_poverty <- lm(form_control_1_poverty,data=data_snails_pond_biom)$residuals
  # second_control_poverty                  <- lm(form_control_2_poverty,data=data_snails_pond_biom)
  # est_iv_control_poverty                  <- coeftest(second_control_poverty,vcov=vcovCL,cluster = ~ vid)
  # est_iv_control_poverty
  
  
  fig <- plot_ly(showscale=F)
  fig <- fig %>% add_surface(x=lrend_q,y=lsurf_q,z=res_pov_eff) %>%
    add_surface(x=lrend_q,y=lsurf_q,z=res_ci_u_eff,opacity = 0.5)%>%
    add_surface(x=lrend_q,y=lsurf_q,z=res_ci_l_eff, opacity = 0.5) %>% 
    add_surface(x=~lrend_q,y=~lsurf_q,z=~matrix(0,length(lrend_q),length(lsurf_q)),
                colorscale = list(c(0,1),c("rgb(1,0,0)","rgb(1,0,0)")),opacity = 0.3) %>%
    layout(
      title = "",
      scene = list(
        xaxis = list(title = "Surface Quantile"),
        yaxis = list(title = "Yield Quantile"),
        zaxis = list(title = "Avg. % Extra Loss")
      ))
  fig
  #### check https://github.com/plotly/orca#installation on how to install orca
  
  orca(fig, file = "poverty_schisto.pdf")
  
  
  
  
  ########## same but total weight instead of yield
  
  
  lweight_q      <- seq(0.1,0.5,0.05)
  lsurf_q        <- seq(0.1,0.5,0.05)
  
  res_pov_w      <- matrix(0,length(lweight_q ),length(lsurf_q))
  res_pov_se_w   <- matrix(0,length(lweight_q ),length(lsurf_q))
  res_pov_pv_w   <- matrix(0,length(lweight_q ),length(lsurf_q))
  res_nopov_w    <- matrix(0,length(lweight_q ),length(lsurf_q))
  res_nopov_pv_w <- matrix(0,length(lweight_q ),length(lsurf_q))
  res_pov_eff_w  <- matrix(0,length(lweight_q ),length(lsurf_q))
  res_pov_eff_se_w <- matrix(0,length(lweight_q ),length(lsurf_q))
  
  
  
  # select 1 for cluster-bootstrapped errors (warning: very long computation time. results are very similar with non-bootstrapped clustering)

  boot_pov <- 1 
  
  for (i in 1:length(lweight_q)){
    for (j in 1:length(lsurf_q)){
      data_snails_pond_biom$poverty <- ifelse(datause$lsurf < quantile(datause$lsurf,lsurf_q[i]) &
                                                datause$poidscult < quantile(datause$poidscult,lweight_q[j]),1,0)
      data_snails_pond_biom$control_f_poverty <- lm(form_control_1_poverty,data=data_snails_pond_biom)$residuals
      second_control_poverty                  <- lm(form_control_2_poverty,data=data_snails_pond_biom)
      if (boot_pov == 1 ){ 
        set.seed(42)
        est_iv_control_poverty                  <- coeftest(second_control_poverty,vcov=vcovBS,cluster = ~ vid,type="wild",R=3000)
      }
      else {
      est_iv_control_poverty                  <- coeftest(second_control_poverty,vcov=vcovCL,cluster = ~ vid)
      }
      res_pov_w[i,j]                            <- est_iv_control_poverty[dim(est_iv_control_poverty)[1],1]
      res_pov_se_w[i,j]                         <- est_iv_control_poverty[dim(est_iv_control_poverty)[1],2]
      res_pov_pv_w[i,j]                         <- est_iv_control_poverty[dim(est_iv_control_poverty)[1],4]
      res_nopov_w[i,j]                          <- est_iv_control_poverty[2,1]
      res_nopov_pv_w[i,j]                       <- est_iv_control_poverty[2,4]
      res_pov_eff_w[i,j]                        <- res_pov[i,j]*mean(datause$intensity[which(datause$lsurf < quantile(datause$lsurf,lsurf_q[i]) & datause$poidscult < quantile(datause$poidscult,lweight_q[j]))])
      res_pov_eff_se_w[i,j]                     <- res_pov_se[i,j]*mean(datause$intensity[which(datause$lsurf < quantile(datause$lsurf,lsurf_q[i]) & datause$poidscult < quantile(datause$poidscult,lweight_q[j]))])
      
    }
  }
  
  
  
  colnames(res_pov_w) <- lweight_q
  rownames(res_pov_w) <- lsurf_q
  colnames(res_pov_pv_w) <- lweight_q
  rownames(res_pov_pv_w) <- lsurf_q
  
  res_ci_u_w <- res_pov_w + 1.96*res_pov_se_w
  res_ci_l_w <- res_pov_w - 1.96*res_pov_se_w
  
  res_pov_eff_w    <- round(res_pov_eff_w*100,2)
  res_pov_eff_se_w <- round(res_pov_eff_se_w*100,2)
  res_ci_u_eff_w   <- res_pov_eff_w + 1.96*res_pov_eff_se_w
  res_ci_l_eff_w   <- res_pov_eff_w - 1.96*res_pov_eff_se_w
  
  # data_snails_pond_biom$control_f_poverty <- lm(form_control_1_poverty,data=data_snails_pond_biom)$residuals
  # second_control_poverty                  <- lm(form_control_2_poverty,data=data_snails_pond_biom)
  # est_iv_control_poverty                  <- coeftest(second_control_poverty,vcov=vcovCL,cluster = ~ vid)
  # est_iv_control_poverty
  
  
  fig_w <- plot_ly(showscale=F)
  fig_w <- fig_w %>% add_surface(x=lweight_q,y=lsurf_q,z=res_pov_eff_w) %>%
    add_surface(x=lweight_q,y=lsurf_q,z=res_ci_u_eff_w,opacity = 0.5)%>%
    add_surface(x=lweight_q,y=lsurf_q,z=res_ci_l_eff_w, opacity = 0.5) %>% 
    add_surface(x=~lweight_q,y=~lsurf_q,z=~matrix(0,length(lweight_q),length(lsurf_q)),
                colorscale = list(c(0,1),c("rgb(1,0,0)","rgb(1,0,0)")),opacity = 0.3) %>%
    layout(
      title = "",
      scene = list(
        xaxis = list(title = "Surface Quantile"),
        yaxis = list(title = "Weight Quantile"),
        zaxis = list(title = "Avg. % Extra Loss")
      ))
  fig_w
  
  orca(fig_w, file = "poverty_schisto_weight.pdf")
  
  
  
  ########## use households who farm on average smaller plots / less productive 
  
  
  
  
  dat_groups_mid <- datause %>% group_by(mid) %>% summarise(weight_mean_mid = mean(poidscult1+poidscult2,na.rm=T) ,
                                                            rend_mean_mid = mean(rend11+rend12), 
                                                            surf_mean_mid =mean(superficie)) %>%
    ungroup()
  
  
  ### lazy way to not create a new dataframe for subsets: backup old ones and rewrite afterwards 
  ### (be careful though, and remember to rewrite again at the end of the loop)
  
  ### put 0 if no removal of cash crops necessary
  
  datause_pov <- left_join(datause,dat_groups_mid)
  
  
  cash_remove <- 1
  
  if (cash_remove == 1){
    
    datause_old <- datause
    data_snails_pond_biom_old <- data_snails_pond_biom
    
    datause <- left_join(datause,dat_groups_mid)
    datause <- datause[-which(datause$cult1==1100 | datause$cult2==1100 ),]
    data_snails_pond_biom <- data_snails_pond_biom[-which(data_snails_pond_biom$cult1==1100 | data_snails_pond_biom$cult2==1100 ),]
  }
  
  
  
  
  
  lrend_q_mid  <- c(0.05,0.1)
  lsurf_q_mid  <- c(0.05,0.1)
  
  
  res_pov_mid            <- matrix(0,length(lrend_q_mid),length(lsurf_q_mid))
  res_pov_se_mid         <- matrix(0,length(lrend_q_mid),length(lsurf_q_mid))
  res_pov_pv_mid         <- matrix(0,length(lrend_q_mid),length(lsurf_q_mid))
  res_nopov_mid          <- matrix(0,length(lrend_q_mid),length(lsurf_q_mid))
  res_nopov_se_mid       <- matrix(0,length(lrend_q_mid),length(lsurf_q_mid)) 
  res_nopov_pv_mid       <- matrix(0,length(lrend_q_mid),length(lsurf_q_mid))
  res_pov_eff_mid        <- matrix(0,length(lrend_q_mid),length(lsurf_q_mid))
  res_pov_eff_se_mid     <- matrix(0,length(lrend_q_mid),length(lsurf_q_mid))
  res_nopov_eff_mid      <- matrix(0,length(lrend_q_mid),length(lsurf_q_mid))
  
  
  
  #### same as before, put 1 for cluster-bootstrapped errors 
  boot_pov <- 1 
  
  
  for (i in 1:length(lrend_q_mid)){
    for (j in 1:length(lsurf_q_mid)){
      data_snails_pond_biom$poverty <- ifelse(datause$weight_mean_mid < quantile(datause$weight_mean_mid,lrend_q_mid[i]) &
                                                datause$surf_mean_mid < quantile(datause$surf_mean_mid,lsurf_q_mid[j]),1,0)
      
      data_snails_pond_biom$control_f_poverty <- lm(form_control_1_poverty,data=data_snails_pond_biom)$residuals
      second_control_poverty                  <- lm(form_control_2_poverty,data=data_snails_pond_biom)
      if (boot_pov == 1 ){ 
        set.seed(42)
        est_iv_control_poverty                  <- coeftest(second_control_poverty,vcov=vcovBS,cluster = ~ vid,type="wild",R=3000)
      }
      else {
        est_iv_control_poverty                  <- coeftest(second_control_poverty,vcov=vcovCL,cluster = ~ vid)
      }
      res_pov_mid[i,j]                            <- est_iv_control_poverty[dim(est_iv_control_poverty)[1],1]
      res_pov_se_mid[i,j]                         <- est_iv_control_poverty[dim(est_iv_control_poverty)[1],2]
      res_pov_pv_mid[i,j]                         <- est_iv_control_poverty[dim(est_iv_control_poverty)[1],4]
      res_nopov_mid[i,j]                          <- est_iv_control_poverty[2,1]
      res_nopov_se_mid[i,j]                       <- est_iv_control_poverty[2,2]
      res_nopov_pv_mid[i,j]                       <- est_iv_control_poverty[2,4]
      res_pov_eff_mid[i,j]                        <- res_pov_mid[i,j]*mean(datause$intensity[which(datause$weight_mean_mid < quantile(datause$weight_mean_mid,lrend_q_mid[i]) & datause$surf_mean_mid < quantile(datause$surf_mean_mid,lsurf_q_mid[j]))])
      res_pov_eff_se_mid[i,j]                     <- res_pov_se_mid[i,j]*mean(datause$intensity[which(datause$weight_mean_mid < quantile(datause$weight_mean_mid,lrend_q_mid[i]) & datause$surf_mean_mid < quantile(datause$surf_mean_mid,lsurf_q_mid[j]))])
      res_nopov_eff_mid[i,j]                      <- res_nopov_mid[i,j]*mean(datause$intensity[which(datause$weight_mean_mid > quantile(datause$weight_mean_mid,lrend_q_mid[i]) & datause$surf_mean_mid > quantile(datause$surf_mean_mid,lsurf_q_mid[j]))])
    }
  }
  
  ### rewrite over the datasets and recover the original env
  
  if (cash_remove == 1){
    datause                   <- datause_old 
    data_snails_pond_biom     <- data_snails_pond_biom_old
    remove(datause_old)
    remove(data_snails_pond_biom_old)
  }
  
  
  colnames(res_pov_mid) <- lrend_q_mid
  colnames(res_pov_mid) <- lsurf_q_mid
  
  
  res_v_pov           <- matrix(0,length(lrend_q_mid),8)
  colnames(res_v_pov) <- c("coef","se","pval","impact","coef_no","se_no","pval_no","impact_no")
  
  for (i in 1:length(lrend_q_mid)){
    res_v_pov[i,1] <- res_pov_mid[i,i] 
    res_v_pov[i,2] <- res_pov_se_mid[i,i]
    res_v_pov[i,3] <- res_pov_pv_mid[i,i]
    res_v_pov[i,4] <- res_pov_eff_mid[i,i]
    res_v_pov[i,5] <- res_nopov_mid[i,i]
    res_v_pov[i,6] <- res_nopov_se_mid[i,i]
    res_v_pov[i,7] <- res_nopov_pv_mid[i,i]
    res_v_pov[i,8] <- res_nopov_eff_mid[i,i]
  }
  
  
  res_v_pov <- round(res_v_pov,3)
  res_v_pov[,c("impact","impact_no")] <- round(res_v_pov[,c("impact","impact_no")]*100,2)
  
  
  
  ticks_pov           <- c("5%","10%")
  plot_pov            <- broom::tidy(res_v_pov)
  plot_pov[,"impact"] <-plot_pov[,"impact"]+ plot_pov[,"impact_no"]
  
  ggplot(plot_pov,aes(coef,ticks_pov))+
    geom_vline(xintercept=0,linetype="solid",colour="red")+
    geom_errorbarh(aes(xmax = coef+1.96*se, xmin =coef-1.96*se),height=0.3,
                   linetype=1,alpha=09,colour="black",size=0.4)+
    # geom_errorbarh(aes(xmax = coefdams+1.96*se, xmin =coefdams-1.96*se),height=0.3,
    #                linetype=1,alpha=09,colour="black",size=0.3)+
    geom_errorbarh(aes(xmax = coef_no+1.96*se_no, xmin =coef_no-1.96*se_no),height=0.3,
                   linetype=3,alpha=09,colour="black",size=0.4)+
    geom_label(aes(coef,label=paste("Below:\n",paste(impact,"%"),sep="")),color="blue",size=4)+
    geom_label(aes(coef_no,label=paste("Above:\n",paste(impact_no,"%"),sep="")),color="blue",size=4)+
    xlab("Estimate") + 
    ylab("Lower Quantile")+
    # scale_x_continuous(breaks=c(-0.001,0.000,0.001,0.002,0.003),limits=c(-0.003,0.0033))+
    theme(axis.text.y = element_text(lineheight=1,face="bold", color="blue", 
                                     size=9, angle=0),legend.key.size = unit(0.3, "cm"))+
    theme_classic2()
  ggsave("est_v_poverty.pdf",height=2,width=5)
  
  
  
  
#################### effect of the presence of a dam / large water reservoir


######## aggregate data at village level
  
  
  ###### to aggregate we use raw data + coordinates: can provide full details upon permission
  
  dat_spatial <- read.csv("data_replication_village.csv")
  dat_spatial$intensity <- matrix(0,dim(dat_spatial)[1],1)
  
  for (i in 1:length(dat_spatial$prev)){
    dat_spatial$intensity[i] <- uniroot(invertPrev, c(0,500), prev = dat_spatial$prev[i] , par = est_par,k=lin_K)$root
  }
  
  
  transf_v_y   <- c("rend_tot", "rend_avg")
  transf_v_x   <- c("surf_mean","dspell_mean","precip_wet", "temp_mean","dspell_max",
                    "temp_night", "temp_day",
                    "totanim","nmem","nyoung","nchildren","nwomen","agemoyen","bovins",
                    "equins","ovins","caprins","dons","vols",
                    "NDVI_yr_mean","EVI_yr_mean","lstN_yr_mean","lstD_yr_mean",
                    "tot_anim_avg","npk_avg","uree_avg","phosph_avg",
                    "pest_solide_avg","pest_liquide_avg",
                    "herbicide_cl_avg","fungicide_g_avg","fungicide_cl_avg",
                    "rodenticide_cl_avg","rodenticide_g_avg")
  
  
  
  dat_ihs_v           <- apply(dat_spatial[,c(transf_v_y,transf_v_x)],2,ihs)
  
  colnames(dat_ihs_v) <- paste("l_",colnames(dat_ihs_v),sep="")
  dat_spatial         <- cbind(dat_spatial,dat_ihs_v)
  
  
  

data_damuse   <- datause


########### provinces containing large dams: 
###### - Bagre Dam (Boulgou province, overlaps into Zoundwéogo, Kouritenga and Ganzourgou provinces): id 30, 14, 7
###### - Kompienga Dam (Kompienga province): id 35
###### - Ziga Dam (Oubritenga province): id 18
###### - Lery Dam (Nayala province): id 40

########### lakes: 
###### - Lake Bam (Bam province): id 1
###### - Lake Tengrela (Comoe province): id 6

dam_provinces   <- c(7,14,18,30,35,40) 
lake_provinces  <- c(1,6)
# only focus on dams



weights <- read.gal("data_spatial_09_11_weights.gal")
W = nb2listw(weights, style="W", zero.policy=T)


res_dam <- matrix(0,1,12)
rownames(res_dam) <- c("Village level, spatial fe")

colnames(res_dam) <- c("Interaction","s.e.inter","p.inter","Intensity","s.e.int","p.int", "Dam","s.e.dam","p.dam","MSE","mean_interaction","mean_intensity")

x            <- c(X1_log, X2, X3_t[-1],X4_log,dis)
x_t          <- c(X1_log, X2, X3_t,X4_log,dis)
dat_spatial$dam <- ifelse(dat_spatial$province %in% c(dam_provinces),1,0)
# x        <- c(X1, X2,X3[1:92],X4) 
# x          <- c(X4)

est_vil_dam_splag <- lagsarlm(l_rend_avg ~ intensity + dam + intensity:dam+l_surf_mean +
                                l_dspell_mean + l_precip_wet + 
                            l_temp_mean + l_dspell_max + l_temp_night + l_temp_day + 
                            l_totanim + l_nmem + l_nyoung + l_nchildren + l_nwomen + 
                            l_agemoyen + l_nchildren + l_bovins + l_equins + l_ovins + 
                            l_caprins + l_NDVI_yr_mean + l_EVI_yr_mean + 
                            l_lstN_yr_mean + l_lstD_yr_mean + l_tot_anim_avg+l_npk_avg+ l_uree_avg + 
                              l_phosph_avg + l_pest_solide_avg + l_pest_liquide_avg + l_herbicide_cl_avg + 
                              l_fungicide_g_avg + l_fungicide_cl_avg + l_rodenticide_cl_avg + 
                              l_rodenticide_g_avg +malaria+factor(annee), 
                            data=dat_spatial, W, type="mixed")

mse.dam_splag_v <- error(est_vil_dam_splag$fitted.values, dat_spatial$l_rend_avg)


res_dam[1,1]   <- est_vil_dam_splag$coefficients[grep("intensity:dam",names(est_vil_dam_splag$coefficients))[1]]
res_dam[1,2]   <- est_vil_dam_splag$rest.se[grep("intensity:dam",names(est_vil_dam_splag$coefficients))[1]]
res_dam[1,3]   <- 2*pt(-abs(res_dam[1,1] /res_dam[1,2]), df = nrow(dat_spatial)-1)


res_dam[1,4]   <- est_vil_dam_splag$coefficients[grep("intensity",names(est_vil_dam_splag$coefficients))[1]]
res_dam[1,5]   <- est_vil_dam_splag$rest.se[grep("intensity",names(est_vil_dam_splag$coefficients))[1]]
res_dam[1,6]   <- 2*pt(-abs(res_dam[1,4] /res_dam[1,5]), df = nrow(dat_spatial)-1)


res_dam[1,7]   <- est_vil_dam_splag$coefficients[grep("dam",names(est_vil_dam_splag$coefficients))[1]]
res_dam[1,8]   <- est_vil_dam_splag$rest.se[grep("dam",names(est_vil_dam_splag$coefficients))[1]]
res_dam[1,9]   <- 2*pt(-abs(res_dam[1,7] /res_dam[1,8]), df = nrow(dat_spatial)-1)


res_dam[1,10]  <- mse.dam_splag_v$err
res_dam[1,11]  <- res_dam[1,1]*mean(dat_spatial$intensity[which(dat_spatial$dam==1)])
res_dam[1,12]  <- res_dam[1,4]*mean(dat_spatial$intensity)


#### plot results of schisto and water resources

res_dam_r <- round(res_dam,3)
res_dam_r[,11] <- 100*res_dam_r[,11]
res_dam_r[,12] <- 100*res_dam_r[,12]
res_dam_r[,11:12] <- round(res_dam_r[,11:12],1)
ticks_dam        <- c("Plot level,\n region fe", "Plot level,\n 2SLS","Village level,\n region fe", "Village level,\n SAR fe")
plot_dam              <- broom::tidy(res_dam_r)



plot_dam_small <- plot_dam
ticks_small <- "Effect of a\n large dam"

ggplot(plot_dam_small,aes(Interaction,"Effect of a\nlarge dam\n(SAR)"))+
  geom_errorbarh(aes(xmax = Interaction+1.96*s.e.inter, xmin =Interaction-1.96*s.e.inter),height=1,
                 linetype=1,alpha=09,colour="black",size=0.3)+
  geom_vline(xintercept=0,linetype="solid",colour="red")+
  geom_errorbarh(aes(xmax = Intensity+1.9*s.e.int, xmin = Intensity-1.96*s.e.int),height=1,
                 linetype=2,alpha=09,colour="black",size=0.3)+  
  geom_label(aes(0.005,label=paste("Dam Effect:\n",paste(100*Dam,"%",sep=""),sep=""),color=p.int),size=4)+
  xlab("") + 
  ylab("")+
  geom_label(aes(Interaction,label=paste("Burden near\n a large dam:\n",paste(mean_interaction+mean_intensity,"%",
                                                                      sep=""),sep=""),color=p.inter),size=4)+
  geom_label(aes(Intensity,label=paste("Burden:\n",paste(mean_intensity,"%",sep=""),sep=""),color=p.int),size=4)+
  scale_x_continuous(breaks=c(-0.02,-0.01,0.0),limits=c(-0.018,0.007))+
  theme(axis.text.y = element_text(lineheight=1,face="bold", color="black", 
                                   size=9, angle=0),legend.key.size = unit(0.3, "cm"))+
  scale_y_discrete(expand=c(0,0.8))+
  scale_color_gradient(name="p.val",low="blue",high="blue",space="Lab",breaks=c(0.5,0.3,0.1))+
  theme_classic2()

ggsave("est_dam_small.pdf",height=2)




########### distance from dams 


res_dam_dist           <- matrix(0,1,12)
colnames(res_dam_dist) <- c("coefdistdams", "se","pval","impact","coefdams", "se","pval","impact","coefdistnodams", "se","pval","impact")
rownames(res_dam_dist) <- c("SAR")

dams_dist <- read.csv("data-dams_distance.csv")
dams_dist_use <- dams_dist[,c(5,8,9,10)]
dams_dist_use$regime_d <- ifelse(dams_dist_use$regime=="TEMPORAIRE",0,1)

dat_spatial_dist           <- left_join(dat_spatial,dams_dist_use,by=c("vid"))


est_vil_dam_splag_dist <- lagsarlm(l_rend_avg ~ intensity+dam+dam_dist+intensity:dam +
                                  intensity:dam_dist+intensity:dam:dam_dist+ dam:dam_dist+
                                l_surf_mean +
                                l_dspell_mean + l_precip_wet + 
                                l_temp_mean + l_dspell_max + l_temp_night + l_temp_day + 
                                l_totanim + l_nmem + l_nyoung + l_nchildren + l_nwomen + 
                                l_agemoyen + l_nchildren + l_bovins + l_equins + l_ovins + 
                                l_caprins + l_NDVI_yr_mean + l_EVI_yr_mean + malaria + dam:dam_dist+
                                l_lstN_yr_mean + l_lstD_yr_mean + l_tot_anim_avg+l_npk_avg+ l_uree_avg + 
                                l_phosph_avg + l_pest_solide_avg + l_pest_liquide_avg + l_herbicide_cl_avg + 
                                l_fungicide_g_avg + l_fungicide_cl_avg + l_rodenticide_cl_avg +l_rodenticide_g_avg+malaria, 
                               data=dat_spatial_dist, W, type="mixed")
 
res_dam_dist[1,1] <- est_vil_dam_splag_dist$coefficients[grep("intensity:dam:dam_dist",names(est_vil_dam_splag_dist$coefficients))[1]]
res_dam_dist[1,2] <- est_vil_dam_splag_dist$rest.se[grep("intensity:dam:dam_dist",names(est_vil_dam_splag_dist$coefficients))[1]]
res_dam_dist[1,3] <- 2*pt(-abs(res_dam_dist[1,1] /res_dam_dist[1,2]), df = dim(est_vil_dam_splag_dist$X)[1])
res_dam_dist[1,4] <- mean(dat_spatial_dist$intensity[which(dat_spatial_dist$dam == 1)])*res_dam_dist[1,1]
res_dam_dist[1,5] <- est_vil_dam_splag_dist$coefficients[grep("intensity:dam",names(est_vil_dam_splag_dist$coefficients))[1]]
res_dam_dist[1,6] <- est_vil_dam_splag_dist$rest.se[grep("intensity:dam",names(est_vil_dam_splag_dist$coefficients))[1]]
res_dam_dist[1,7] <- 2*pt(-abs(res_dam_dist[1,5] /res_dam_dist[1,6]), df = dim(est_vil_dam_splag_dist$X)[1])
res_dam_dist[1,8] <- mean(dat_spatial_dist$intensity[which(dat_spatial_dist$dam == 1)])*res_dam_dist[1,5]
res_dam_dist[1,9] <- est_vil_dam_splag_dist$coefficients[grep("intensity:dam_dist",names(est_vil_dam_splag_dist$coefficients))[1]]
res_dam_dist[1,10] <- est_vil_dam_splag_dist$rest.se[grep("intensity:dam_dist",names(est_vil_dam_splag_dist$coefficients))[1]]
res_dam_dist[1,11] <- 2*pt(-abs(res_dam_dist[1,9]/res_dam_dist[1,10]), df = dim(est_vil_dam_splag_dist$X)[1])
res_dam_dist[1,12] <- mean(dat_spatial_dist$intensity[which(dat_spatial_dist$dam == 0)])*res_dam_dist[1,9]

res_dam_dist_r           <- round(res_dam_dist,5)
res_dam_dist_r[,4*c(1:3)] <- round(100*res_dam_dist_r[,4*c(1:3)],2)

### take vil with dam=1, then take mean distance, then take mean yield among those at the mean (btw 50 and 51)

ticks_dam_dist            <- c("SAR")
plot_dam_dist              <- broom::tidy(res_dam_dist_r)
plot_dam_dist_small       <- plot_dam_dist



ggplot(plot_dam_dist_small,aes(coefdistdams,"Distance from\nwater resources\n(SAR)"))+
  geom_vline(xintercept=0,linetype="solid",colour="red")+
  geom_errorbarh(aes(xmax = coefdistdams+1.96*se, xmin =coefdistdams-1.96*se),height=0.75,
                 linetype=1,alpha=09,colour="black",size=0.2)+
  # geom_errorbarh(aes(xmax = coefdams+1.96*se, xmin =coefdams-1.96*se),height=0.3,
  #                linetype=1,alpha=09,colour="black",size=0.3)+
  geom_errorbarh(aes(xmax = coefdistnodams+1.96*se, xmin =coefdistnodams-1.96*se),height=0.75,
                 linetype=3,alpha=09,colour="black",size=0.2)+
  geom_label(aes(coefdistdams,label=paste("Distance\n(large dam)\n",paste(impact,"%"),sep="")),color="blue",size=4)+
  geom_label(aes(coefdistnodams,label=paste("Distance\n",paste(impact.2,"%"),sep="")),color="red",size=4)+
  geom_label(aes(-0.0025,label=paste("Next to\na large dam\n",paste(impact.1,"%"),sep="")),color="blue",size=4)+
  xlab("Estimate") + 
  ylab("")+
  scale_x_continuous(breaks=c(-0.001,0.000,0.001,0.002,0.003),limits=c(-0.003,0.0033))+
  theme(axis.text.y = element_text(lineheight=1,face="bold", color="blue", 
                                   size=9, angle=0),legend.key.size = unit(0.3, "cm"))+
  scale_color_gradient(name="pval",low="blue",high="red",space="Lab",breaks=c(0.5,0.3,0.1))+
  theme_classic2()
ggsave("est_dam_dist_small.pdf",height=2,width=6.8)



####### non-parametric surface of interaction disease/distance from dam


est_np_damdist <-  gam(l_rend_avg ~ s(intensity,dam_dist)+ l_surf_mean + l_dspell_mean + l_precip_wet + 
                         l_temp_mean + l_dspell_max + l_temp_night + l_temp_day + 
                         l_totanim + l_nmem + l_nyoung + l_nchildren + l_nwomen + 
                         l_agemoyen + l_nchildren + l_bovins + l_equins + l_ovins + 
                         l_caprins + l_dons + l_vols + l_NDVI_yr_mean + l_EVI_yr_mean + 
                         l_lstN_yr_mean + l_lstD_yr_mean + l_tot_anim_avg + l_uree_avg + 
                         l_phosph_avg + l_pest_solide_avg + l_pest_liquide_avg + l_herbicide_cl_avg + 
                         l_fungicide_g_avg + l_fungicide_cl_avg + l_rodenticide_cl_avg + 
                         l_rodenticide_g_avg + I(l_surf_mean^2) + factor(annee),
                       family=gaussian,data=dat_spatial_dist)

est_np_damdist   <- getViz(est_np_damdist)
pdf("intensity_distance_dams.pdf")
pl_int_dist        <-  plot( sm(est_np_damdist,1))+
                        labs(x="Intensity", y="Distance",title=" ", textsize=20) + theme_classic()+
     theme(axis.text.x= element_text( size=20),axis.text.y= element_text( size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),panel.border=element_blank(),
        plot.background = element_blank(),
        panel.grid= element_blank(),axis.line = element_line(color = 'black') ) 
pl_int_dist
dev.off()
# 

