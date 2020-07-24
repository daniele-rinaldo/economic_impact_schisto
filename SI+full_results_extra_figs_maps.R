
########### Code for the Supplementary Information Results
########### Robustness checks, subsets, crop choice analysis




############## Repeated households / villages in the full dataset around the 2010 cutoff 
############## in order to replicate Table 1 and 2 in the SI we include a partial dataset with all villages and households, but only with year and village/household id  
############## no coordinates and other covariates are included 


data    <- read.csv("data_replication_year+id.csv")

yr_map1 <- sort(c(2003,2008,2009,2010,2007,2005, 2006))
yr_map2 <- sort(c(2014 ,2015,2011, 2012, 2013, 2016, 2017))

rep_hh <- matrix(0,nrow=length(yr_map1),ncol=length(yr_map2))
rep_vil <- matrix(0,nrow=length(yr_map1),ncol=length(yr_map2))
rownames(rep_hh) <- yr_map1
colnames(rep_hh) <- yr_map2
rownames(rep_vil) <- yr_map1
colnames(rep_vil) <- yr_map2


for (i in 1:length(yr_map1)){
  for (j in 1:length(yr_map2)){
    
    dat_y <- data %>% filter(annee==yr_map1[i] | annee==yr_map2[j])
    
    n_rep_hh <- dat_y%>%
      group_by(mid) %>%
      summarise(numyears = length(unique(annee)),years=paste(unique(annee),collapse=","),
                n = n()) %>%
      arrange(desc(numyears))
    rhh <- dim(n_rep_hh[which(n_rep_hh$numyears==2),])[1]
    rep_hh[i,j] <- rhh
    
    
    n_rep_vil <- dat_y %>%
      group_by(vid) %>%
      summarise(numyears = length(unique(annee)),years=paste(unique(annee),collapse=","),
                n = n()) %>%
      arrange(desc(numyears))
    rv <- dim(n_rep_vil[which(n_rep_vil$numyears==2),])[1]
    rep_vil[i,j] <- rv
    
    remove(dat_y)
  }
}




############ Replication of the Supplementary Information results
############ Cleaned dataset for the years 2009 and 2011



res           <- matrix(0,10,6)
res_sc        <- matrix(0,2,6)
colnames(res)    <- c("Estimate","S.E.","p-val","Mean Effect (%)","Top 5% Effect (%)","MSE")
colnames(res_sc) <- c("Estimate","S.E.","p-val","Mean Effect (%)","Top 5% Effect (%)","MSE")



datause <- read.csv("data_replication.csv") 
datause$rend <- datause$rend11 + datause$rend12



################# linear regressions, fixed effects, lewbel IV/measurement error


###### fixed effects 

mid_t            <- cbind(mid,annee)
vid_t            <- cbind(vid,annee)
prov_t           <- cbind(province,annee)
com_t            <- cbind(comid,annee)
# # re-do the automatized partialling out for interacted fe
# datause$midandt          <- as.factor(paste(datause$mid,datause$annee,sep="")  )
# colnames(datause$midandt) <- midandt
# choose if you want different partialled out datasets

# fe                <- list(mid_t)
fe                <- list(prov_t, com_t, vid_t, mid_t)
# fe                <- list(midandt,prov_t)
# 
# vars to partial out
varlist           <- c(y,p.sq,int.sq, p.flex,int.flex, X1, X1_log, X2, X3, X4, X4_log, dis)

for (i in 1:length(fe)){
  
  cat('Partialling out dataset by',colnames(fe[[i]]),'\n')
  # fe.var          <- as.factor(fe[[i]][,1])
  # fe.t            <- as.factor(fe[[i]][,2])
  # feff            <-  as.data.frame(model.matrix(~fe.var+0))
  # feff[,1]        <- NULL
  # feff_t          <- as.data.frame(model.matrix(~fe.t+0))
  # feff_t[,1]        <- NULL
  
  # fixed effects
  fixed           <- colnames(fe[[i]])
  
  # partialled out dataset 
  # wdata           <- as.data.frame(cbind(mid,vid,comid,midandt,annee))
  wdata           <- as.data.frame(cbind(mid,vid,comid,annee))
  
  for(j in 1:length(varlist)){
    cat(j,'/',length(varlist),'\n')
    # use felm package for partialling out since for larger factors (mid) computation time is large
    form_fe             <- as.formula(paste(varlist[j], "~ 1 |", paste(fixed,collapse="+"),"| 0 | 0"))
    wdata[, varlist[j]] <- felm(form_fe, data=datause)$res
  }
  assign(paste("wdata",colnames(fe[[i]])[1],sep="_"),wdata)
  remove(wdata)
}



######### plot data in both years


# dat_fe_diff <- data.frame(cbind(datause$prev,wdata_mid$prev+mean(datause$prev)))
# 
# ggplot(data=dat_fe_diff, aes(x=V1)) + 
#   geom_density(alpha=0.5,adjust=2) +
#   geom_vline(xintercept=mean(dat_fe_diff$prev),color="red", linetype="dashed", size=1)
# 
# plot(density(dat_fe_diff$V1))
# lines(density(dat_fe_diff$prev))
##
# unusual result: household fe do not change the distribution almost at all
# very little heterogeneity across households
# 
dat_fe_diff <- data.frame(cbind(datause$lrend,wdata_vid$lrend+mean(datause$lrend)))

ks.test(datause$lrend,wdata_mid$lrend+mean(datause$lrend))


ggplot(dat_fe_diff)+
  geom_density(aes(x=V1),alpha=.2, fill="blue")+
  geom_density(aes(x=lrend),colour="black",fill="#FF6666",alpha=0.5)+
  xlab("Log Yield/Ha")+ylab("Density")+
  theme_bw()+
  theme(axis.text.x= element_text( size=13),axis.text.y= element_text( size=13),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),panel.border=element_blank(),
        plot.background = element_blank(),
        panel.grid= element_blank(),axis.line = element_line(color = 'black')
  )
# ggtitle(paste("Densities of original data vs. household fe, ",paste(unique(datause$annee),collapse="-"),sep=""))
ggsave(paste("density_fe_",paste(unique(datause$annee),collapse="_"),".pdf",sep=""))

# wow absolutely the same distr
# pdf(paste("density_fe_",paste(unique(datause$annee),collapse="_"),".pdf",sep=""))


##### so: fixed effects do not impact lrend at all
# they do prevalence, so the variation in the coefficient between nofe+fe is because of cross-sectional variation 
# in unobserved REAL prevalence (not the measure we have) and the individual hh characteristics


############


form_prev <- as.formula(paste(y,"~-1 +", paste(p,paste(X1_log, collapse="+"),
                                               paste(X3,collapse="+"),paste(X2,collapse="+"),
                                               paste(X4_log,collapse="+"),dis,sep="+")))

form_int <- as.formula(paste(y,"~ -1 +", paste(int,paste(X1_log, collapse="+"),
                                               paste(X3,collapse="+"),paste(X2,collapse="+"),
                                               paste(X4_log,collapse="+"),dis,sep="+")))

form_no_t_fe_prev <- as.formula(paste(y,"~ ", paste(p,paste(X1_log, collapse="+"),
                                                    paste(X3_t[1:93],collapse="+"),
                                                    paste(X2,collapse="+"),
                                                    paste(X4_log,collapse="+"),dis,
                                                    sep="+")))

form_no_t_fe_int <- as.formula(paste(y,"~ ", paste(int,paste(X1_log, collapse="+"),
                                                   paste(X3_t[1:93],collapse="+"),paste(X2,collapse="+"),
                                                   paste(X4_log,collapse="+"),dis,sep="+")))

form_large_prev <- as.formula(paste("lrend~-1+prev+",paste(x.flex,collapse="+"),sep=""))
form_large_int <- as.formula(paste("lrend~-1+intensity+",paste(x.flex,collapse="+"),sep=""))


## estimation with only time fixed effects (detrending)
# est_nofe_prev <- lm(form_no_t_fe_prev,data=datause)
est_nofe_int <- lm(form_no_t_fe_int,data=datause)

# estimation with time and household fixed effects
# est_prev_mid <- lm(form_prev, data=wdata_mid)
est_int_mid <- lm(form_int, data=wdata_mid)



# estimation with large set of controls (orthogonal polynomials on most continuous variables)
# est_mid_large_prev <- lm(form_large_prev,data=wdata_mid)
est_mid_large_int <- lm(form_large_int,data=wdata_mid)

# clustering at village level
# vcov_prev <- cluster.vcov(est_prev_mid,vid)
vcov_int_time <- cluster.vcov(est_nofe_int,vid)
vcov_int <- cluster.vcov(est_int_mid,vid)
# vcov_large_prev <- cluster.vcov(est_mid_large_prev,vid)
vcov_large_int  <- cluster.vcov(est_mid_large_int,vid)

# est_cl_prev_mid <- coeftest(est_prev_mid,vcov_prev)
est_cl_int     <- coeftest(est_nofe_int,vcov_int_time)
est_cl_int_mid <- coeftest(est_int_mid,vcov_int)
# est_cl_prev_mid_large <- coeftest(est_mid_large_prev,vcov_large_prev)
est_cl_int_mid_large <- coeftest(est_mid_large_int,vcov_large_int)


res[1,1:3] <- est_cl_int[2,c(1,2,4)]
mse.est_nofe_int <- error(est_nofe_int$fitted.values,datause$lrend)
res[1,4] <- est_cl_int[2,1]*mean(intensity)
res[1,5] <- est_cl_int[2,1]*quantile(intensity,0.95)
res[1,6] <- mse.est_nofe_int$err


res[2,1:3] <- est_cl_int_mid[1,c(1,2,4)]
mse.est_int <- error(est_int_mid$fitted.values,wdata_mid$lrend)
res[2,4] <- est_cl_int_mid[1,1]*mean(intensity)
res[2,5] <- est_cl_int_mid[1,1]*quantile(intensity,0.95)
res[2,6] <- mse.est_int$err

# 
# res[3,1:3] <- est_cl_int_mid_large[1,c(1,2,4)]
# mse.est_int_large <- error(est_mid_large_int$fitted.values,wdata_mid$lrend)
# res[3,4] <- est_cl_int_mid_large[1,1]*mean(intensity)
# res[3,5] <- est_cl_int_mid_large[1,1]*quantile(intensity,0.95)
# res[3,6] <- mse.est_int_large$err
# 




########### snails as IV  



##### prepare snails dataset

snails_pond      <- read.csv("pond_means.csv")
data_snails_pond <- left_join(datause,snails_pond,by=c("vid"))
snails_river     <- read.csv("river_yearly_means.csv")

snails2          <- snails_river %>% filter(habid=="perm_river_biomphalaria" | habid=="eph_river_bulinus")

snails_pond_biom      <- left_join(snails_pond, snails2, by=c("vid"))
snails_pond_biom[is.na(snails_pond_biom)] <- 0
data_snails_pond_biom <- left_join(datause,snails_pond_biom,by=c("vid"))



summary(felm(lrend~lsurf | annee | (intensity~dist+ mean+dry+rainy) | vid, data=data_snails_pond_biom))

##### check for villages with presence of biomphalaria snail whether mean of snails is positively corr with schisto
datacheck <- data_snails_pond_biom[which(!is.na(data_snails_pond_biom$habid)),]

snails_pond_0911 <- read.csv("snails_pond_09_11.csv")
data_snails_pond2 <- left_join(datause,snails_pond_0911,by=c("vid"))

data_snails_pond2$dry <- ifelse(data_snails_pond2$annee==2009,data_snails_pond2$dry2009,data_snails_pond2$dry2011)
data_snails_pond2$winter <- ifelse(data_snails_pond2$annee==2009,data_snails_pond2$winte2009,data_snails_pond2$winter2011)
data_snails_pond2$rainy <- ifelse(data_snails_pond2$annee==2009,data_snails_pond2$rainy2009,data_snails_pond2$rainy2011)
data_snails_pond2$yearly <- (data_snails_pond2$rainy+data_snails_pond2$winter+data_snails_pond2$dry)/3

##### instruments 

instr <- c("rainy","dry","mean")


# remove factors with almost no values (reduces noise in estimation by a lot without changing estimates)
form_1stage_tfe <-  as.formula(
  paste(paste(int,"~ ",paste(paste(instr,collapse="+"),paste(x_log,collapse="+"),sep="+")),sep=""))

form_2stage_tfe <-  as.formula(
  paste(paste(y,"~ ",paste("stage1_tfe",paste(x_log,collapse="+"),sep="+")),sep=""))


form_2stage_tfe_gam <-  as.formula(
  paste(paste(y,"~ ",paste("s(stage1_tfe)",paste(x_log,collapse="+"),sep="+")),sep=""))

data_snails_pond_biom$stage1_tfe <- fitted(lm(form_1stage_tfe,data=data_snails_pond_biom))
# data_snails_pond2$stage1_tfe <- fitted(lm(form_1stage_tfe,data=data_snails_pond2))
second_tfe            <- lm(form_2stage_tfe,data=data_snails_pond_biom)


summary(lm(intensity ~ rainy + lsurf + lsurfsq + remune + entraid + 
             gestion + parc_recup + mode_labour + l_tot_anim_traction + 
             l_males_elev + l_fem_elev + l_totanim + l_morts + l_naissances + 
             l_bovins + l_equins + l_ovins + l_caprins + l_dons + l_vols + 
             l_autocons + l_ventebovins + l_achatbovins + l_venteovins + 
             l_achatovins + l_ventecaprins + l_achatcaprins + l_duree_jachere + 
             l_npk_kg + l_uree_kg + l_phosph_kg + l_pest_solide_g + l_pest_liquide_cl + 
             l_herbicide_g + l_herbicide_cl + l_fungicide_g + l_fungicide_cl + 
             l_rodenticide_g + l_rodenticide_cl + nmem + agemoyen + age_chefmen + 
             numfamilies + nchildren + nyoung + nwomen + culture_pluviale + 
             culture_maraichere + arboriculture + autre_contresaison + 
             peche + factor_years_1 + factor_cult1_1 + factor_cult1_2 + 
             factor_cult1_3 + factor_cult1_4 + factor_cult1_5 + factor_cult1_6 + 
             factor_cult1_7 + factor_cult1_8 + factor_cult1_9 + factor_cult1_10 + 
             factor_cult1_11 + factor_cult1_12 + factor_cult1_13 + factor_cult1_14 + 
             factor_cult1_15 + factor_cult1_16 + factor_cult1_17 + factor_cult1_18 + 
             factor_cult1_19 + factor_cult1_20 + factor_cult1_21 + factor_cult1_22 + 
             factor_cult1_23 + factor_cult1_24 + factor_cult1_25 + factor_cult1_26 + 
             factor_cult1_27 + factor_cult1_28 + factor_cult1_29 + factor_cult1_30 + 
             factor_cult2_1 + factor_cult2_2 + factor_cult2_3 + factor_cult2_4 + 
             factor_cult2_5 + factor_cult2_6 + factor_cult2_7 + factor_cult2_8 + 
             factor_cult2_9 + factor_cult2_10 + factor_cult2_11 + factor_cult2_12 + 
             factor_cult2_13 + factor_cult2_14 + factor_cult2_15 + factor_cult2_16 + 
             factor_cult2_17 + factor_cult2_18 + factor_cult2_19 + factor_cult2_20 + 
             factor_cult2_21 + factor_cult2_22 + factor_cult2_23 + factor_cult2_24 + 
             factor_cult2_25 + factor_cult2_26 + factor_localisation_1 + 
             factor_localisation_2 + factor_localisation_3 + factor_relief_1 + 
             factor_relief_2 + factor_antierosif_1 + factor_antierosif_2 + 
             factor_antierosif_3 + factor_antierosif_4 + factor_antierosif_5 + 
             factor_antierosif_6 + factor_antierosif_7 + factor_antierosif_8 + 
             factor_mode_acq_1 + factor_mode_acq_2 + factor_mode_acq_3 + 
             factor_mode_acq_4 + factor_mode_acq_5 + factor_niv_acq_1 + 
             factor_niv_acq_2 + factor_niv_acq_3 + factor_niv_acq_4 + 
             factor_niv_acq_5 + factor_niv_acq_6 + factor_niv_acq_7 + 
             factor_type_labour_1 + factor_type_labour_2 + factor_type_labour_3 + 
             factor_educ_max_1 + factor_educ_max_2 + factor_educ_max_3 + 
             factor_educ_max_4 + factor_educ_max_5 + factor_educ_max_6 + 
             factor_perte_1 + factor_perte_4 + factor_perte_6 + factor_perte_8 + 
             l_temp_mean + l_temp_day + l_temp_night + l_precip_wet + 
             l_temp_dry + l_dspell_mean + l_dspell_max + l_lstD_yr_mean + 
             l_lstN_yr_mean + l_EVI_yr_mean + l_NDVI_yr_mean + malaria,data=data_snails_pond_biom))

ivnp <- gam(form_2stage_tfe_gam,data=data_snails_pond_biom)
plot(ivnp)
# data_snails_pond_biom$inter <- data_snails_pond_biom$annee*as.numeric(data_snails_pond_biom$region)
# 
# summary(felm(as.formula(paste(y,"~",paste(x_t,collapse="+"),"| annee+region| (intensity~rainy+dry+mean) | vid")),
#              data=data_snails_pond_biom))
# 
# check if this is valid for year*region fe: yes it is 


# cluster bootstrap with vid + time: extremely sig

# est_iv_tv <- coeftest(second_tfe,vcov=vcovBS,cluster = ~ annee+vid,type="wild",R=250,cores=4)
# cluster bootstrap with only vid: decent sig
set.seed(42)
est_iv  <- coeftest(second_tfe,vcov=vcovBS,cluster = ~ vid,type="wild",R=3000)
est_iv3 <- coeftest(second_tfe,vcov=vcovCL,cluster = ~ vid)

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




res[3,1:3] <- est_iv[2,c(1,2,4)]
res[3,4]   <- res[3,1]*mean(intensity)
res[3,5]   <- res[3,1]*quantile(intensity,0.95) 
res[3,6]   <- error(fitted(second_tfe),data_snails_pond$lrend)$err
res_sc[1,] <- res[3,]



######## Internal instruments as IV: Lewbel (2012)

cat("- Lewbel hetIV-GMM \n")

### ?is there heterosked?
het1 <- lmtest::bptest(est_int_mid)
het2 <- lmtest::bptest(est_nofe_int)
# omg yes, even after partialling out



y = "lrend"
endo_prev = "prev"
endo_int ="intensity" 

# select factors in non-fe case because of weird issues in inverting/decomposing the X matrix (too many 0/1s?).
# changes very little in magnitude of estimates, just to be sure
exog = paste(paste(X1,collapse="+"),paste(X2,collapse="+"),"factor_years_1",paste(colnames(perte_d)[c(1,4,6,8)],collapse="+"),paste(X4,collapse="+"),dis,sep="+")
exog_mid = paste(paste(X1,collapse="+"),paste(X2,collapse="+"),paste(X3[-1],collapse="+"),paste(X4,collapse="+"),dis,sep="+")
# try to use as instruments variables that are heteroskedastic across levels that are not partialled out:
#### i.e. use household and village-level variables if the data has household feff already removed (no plot ones)
hetv <- paste(paste(X2,collapse="+"),sep="+")

mod1 <- felm(as.formula(paste(endo_int,"~",exog_mid,sep="")) | annee+mid | 0 | 0,data=datause)
Z    <- model.matrix(as.formula(paste("~", hetv, sep=" + ")), data = datause)
Z    <- apply(Z, 2, function(x){x-mean(x)} )
colnames(Z) <- paste("Z", 1:ncol(Z), sep=".")
inst        <- matrix(apply(mod1$residuals, 2, function(x){x*Z}), nrow = dim(datause)[1])
datalew     <- cbind(datause,inst)
form_lew    <- as.formula(paste(y,"~",paste(exog_mid,collapse="+")," | annee+mid | (intensity~",paste(paste("`",seq(1:dim(Z)[2]),"`",sep=""),collapse="+")," )| parid",sep=""))
lew         <- felm(form_lew,data=datalew)


res[4,1] <- lew$coefficients[length(lew$coefficients)]
res[4,2] <- lew$se[length(lew$coefficients)]
res[4,3] <- lew$pval[length(lew$coefficients)]
res[4,4] <- res[4,1]*mean(datause$intensity)
res[4,5] <- res[4,1]*quantile(datause$intensity,0.95)
res[4,6] <- error(lew$fitted.values,datause$lrend)$err





######## double machine learning analysis (double Neyman orthogonal moment conditions)


##### full set, with K cross-validation splits and tuned methods

K            <- 3
trim         <- c(0,1)
set.seed(42)
# watch out for cv.folds option - sometimes it breaks
Boosting     <- list(bag.fraction = .5, train.fraction = 1.0, interaction.depth=2, n.trees=1000, shrinkage=.01, n.cores=1, cv.folds=1, verbose = FALSE, clas_dist= 'adaboost', reg_dist='gaussian')
Forest       <- list(clas_nodesize=1, reg_nodesize=5, ntree=500, na.action=na.omit, replace=TRUE)
RLasso       <- list(penalty = list(homoscedastic = FALSE, X.dependent.lambda =FALSE, lambda.start = NULL, c = 1.1), intercept = TRUE)
Nnet         <- list(size=2,  maxit=1000, decay=0.01, MaxNWts=5000,  trace=FALSE)
Trees        <- list(reg_method="anova", clas_method="class")
arguments    <- list(Boosting=Boosting, Forest=Forest, RLasso=RLasso, Nnet=Nnet, Trees=Trees)

# ensemble     <- list(methods=c("RLasso", "Boosting", "Forest", "Nnet"))
# methods      <- c("Boosting","Forest")
# the forest takes really forever (like 1h)
methods       <- c("Boosting","RLasso","Nnet")
methods_f     <- c("Forest")


# for neural networks, rescaling the data matrix requires removing non-numeric variables (geo ids)
nnum          <- sapply(datause,is.numeric)
datause_ML    <- datause[,nnum]
x.flex_t      <- c(x.flex,"factor_years_1")

est_full_DML        <- Full_DML(datause_ML, y, int, x_t, x.flex_t, methods=methods, nfold=K, est="plinear", arguments=arguments, ensemble=ensemble, silent=FALSE, trim)
est_full_DML_mid    <- Full_DML(wdata_mid, y, int, x, x.flex, methods=methods, nfold=K, est="plinear", arguments=arguments, ensemble=ensemble, silent=FALSE, trim)
# est_full_DML_mid_f  <- Full_DML(wdata_mid, y, int, x[1:2], x.flex, methods=methods_f, nfold=K, est="plinear", arguments=arguments, ensemble=ensemble, silent=FALSE, trim)

# pvals_DML          <- 2*pt(-abs(est_full_DML[1,c(1:(dim(est_full_DML)[2]-1))]/est_full_DML[2,c(1:(dim(est_full_DML)[2]-1))]), df = nrow(datause)-2)
pvals_DML_mid      <- 2*pt(-abs(est_full_DML_mid[1,c(1:(dim(est_full_DML_mid)[2]-1))]/est_full_DML_mid[2,c(1:(dim(est_full_DML_mid)[2]-1))]), df = nrow(datause)-2 - length(unique(datause$mid)))
# pvals_DML_mid_f      <- 2*pt(-abs(est_full_DML_mid_f[1,c(1:(dim(est_full_DML_mid_f)[2]-1))]/est_full_DML_mid_f[2,c(1:(dim(est_full_DML_mid_f)[2]-1))]), df = nrow(datause)-2 - length(unique(mid)))

# best                <- which(est_full_DML[3,1:2]==min(est_full_DML[3,1:2]))
best_mid            <- which(est_full_DML_mid[3,1:2]==min(est_full_DML_mid[3,1:2]))

# 
# res[7,1:2] <- est_full_DML_mid_f[1:2]
# res[7,3]   <- pvals_DML_mid_f
# res[7,4]   <- est_full_DML_mid_f[1]*mean(intensity)
# res[7,5]   <- est_full_DML_mid_f[1]*quantile(intensity,0.95)
# res[7,6]   <-  est_full_DML_mid_f[3]


res[5,1:2] <- est_full_DML_mid[1:2,best_mid]
res[5,3]   <- pvals_DML_mid[best_mid]
res[5,4]   <- est_full_DML_mid[1]*mean(intensity)
res[5,5]   <- est_full_DML_mid[1]*quantile(intensity,0.95)
res[5,6]   <-  est_full_DML_mid[3,best_mid]




######### ML + IV 


# put tune=TRUE at your peril (takes around 6 hrs) 

est_ML_base_int_mid_IV <- DML_Forest_IV(y,d=int,x_t,instr=instr,tune=F,data=data_snails_pond_biom,boot=T,id="annee",cl="vid")
# est_ML_base_int_mid_IV <- DML_Forest_IV(y,d=int,x_t,instr=instr,tune=F,data=data_snails_pond,boot=T,id="annee",cl="vid")


res[6,1] <- est_ML_base_int_mid_IV$theta_2sls
res[6,2] <- est_ML_base_int_mid_IV$theta_2sls.se
res[6,3] <- est_ML_base_int_mid_IV$theta_2sls.p
res[6,4]   <- est_ML_base_int_mid_IV$theta_2sls*mean(intensity)
res[6,5]   <- est_ML_base_int_mid_IV$theta_2sls*quantile(intensity,0.95)
res[6,6]   <- est_ML_base_int_mid_IV$mse.y
res_sc[2,] <- res[6,]

detach(datause)


################## aggregate at village level + spatial analysis



cat(paste("Estimation for repeated villages in",paste(unique(datause$annee),collapse=" and "), "....") )   


#### aggregate data per village


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


form_est_fe_v      <- as.formula(paste("l_rend_avg","~",
                                       int,"+",
                                       paste(colnames(dat_ihs_v)[-c(1,2)],collapse="+"),
                                       "+I(l_surf_mean^2)+malaria| annee+comid| 0 |comid",sep=""))     

est_fe_vil          <- felm(form_est_fe_v, data=dat_spatial)



res[7,1]  <- est_fe_vil$coefficients[1]
res[7,2]  <- est_fe_vil$se[1]
res[7,3]  <- est_fe_vil$pval[1]
res[7,4]  <- est_fe_vil$coefficients[1]*mean(dat_spatial$intensity)
res[7,5]  <- est_fe_vil$coefficients[1]*quantile(dat_spatial$intensity,0.95)
res[7,6]  <- error(est_fe_vil$fitted.values,dat_spatial$l_rend_avg)$err



##### only time fe does not work. probably because of spatial correlation
# since it decays quickly commune fe take care of it


##### use spatial weights to calculate Moran I for both models with and without comid fe


weights <- read.gal("data_spatial_09_11_weights.gal")
W = nb2listw(weights, style="W", zero.policy=T)



form_est_fe_t_v  <- as.formula(paste("l_rend_avg","~",
                                     int,"+",
                                     paste(colnames(dat_ihs_v)[-c(1,2)],collapse="+"),
                                     "+I(l_surf_mean^2)+malaria+ factor(annee)",sep=""))  
form_est_fe_com_v  <- as.formula(paste("l_rend_avg","~",
                                       int,"+",
                                       paste(colnames(dat_ihs_v)[-c(1,2)],collapse="+"),
                                       "+I(l_surf_mean^2)+malaria+ factor(comid)+factor(annee)",sep=""))  



est_vil       <- lm(form_est_fe_t_v , data=dat_spatial)
est_vil_comid <- lm(form_est_fe_com_v, data=dat_spatial)

m1 <- lm.morantest(est_vil, W, alternative="two.sided")
m2 <- lm.morantest(est_vil_comid, W, alternative="two.sided")

est_vil_splag <- lagsarlm(form_est_fe_t_v, data=dat_spatial,
                          W,type="mixed")

summary(est_vil_splag)

spgm(l_rend_avg~intensity,data=data_snails_pond_biom_v,listw=W,model='within',lag=TRUE)



res[8,1]  <- est_vil_splag$coefficients[2]
res[8,2]  <- est_vil_splag$rest.se[2]
res[8,3]  <-  2*pt(-abs(res[7,1] /res[7,2]), df = nrow(dat_spatial)-1)
res[8,4]  <- res[8,1]*mean(dat_spatial$intensity)
res[8,5]  <- res[8,1]*quantile(dat_spatial$intensity,0.95)
res[8,6]  <- error(est_vil_splag$fitted.values,dat_spatial$l_rend_avg)$err


######## spatial random effects 

form_est_fe_spre_v  <- as.formula(paste("l_rend_avg","~",
                                        int,"+",
                                        paste(colnames(dat_ihs_v)[-c(1,2)],collapse="+"),
                                        "+I(l_surf_mean^2)+malaria+ factor(annee)+Matern(1|lat+long)",sep=""))  


m_spatial_int_avg <- fitme(form_est_fe_spre_v,data = dat_spatial, family ="gaussian") 

res[9,1:2] <- summary(m_spatial_int_avg)[[3]][2,1:2]
res[9,3]   <- 2*pt(-abs(summary(m_spatial_int_avg)[[3]][2,1]/summary(m_spatial_int_avg)[[3]][2,2]), df = nrow(dat_spatial)-2)
res[9,4]   <- m_spatial_int_avg$fixef[2]*mean(dat_spatial$intensity)
res[9,5]   <- m_spatial_int_avg$fixef[2]*quantile(dat_spatial$intensity,0.95)
res[9,6]   <- error(predict(m_spatial_int_avg),dat_spatial$l_rend_avg)$err

sims_avg   <- simulateResiduals(m_spatial_int_avg,n=500)
pdf(paste("qqplot_spatial",paste(unique(dat_spatial$annee),collapse="_"),".pdf"))
testUniformity(sims_avg)
dev.off()


# more reasonable the intensity measure for the sum of yield
# average yield works well


MLdistMat <-  as.matrix(proxy::dist(dat_spatial[,c("lat","long")]))
hlcor     <-  HLCor(rend_avg ~ intensity+Matern(1|lat+long),data=dat_spatial,
                    distMatrix=MLdistMat,HLmethod="ML",
                    ranPars=list(nu=m_spatial_int_avg$corrPars$`1`[[1]],rho=m_spatial_int_avg$corrPars$`1`[[2]]))  



##### uncomment for plotting spatial correlation
##### one is probably enough, and it takes a while, so it's commented by default
dist     <- dist(dat_spatial[,c("lat","long")])
sp_corr  <- MaternCorr(dist, nu =m_spatial_int_avg$corrPars$`1`[[1]], rho =m_spatial_int_avg$corrPars$`1`[[2]])
spc       <- data.frame(dist=as.numeric(dist),sp=as.numeric(sp_corr))

# dev.off()
ggplot(spc,aes(dist,sp,color=sp)) +
  geom_point() +
  xlab("Village pairs distance (coordinate degrees)")+
  ylab("Estimated spatial correlation")+
  theme_bw()+
  theme(axis.text.x= element_text( size=11),axis.text.y= element_text( size=11),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13),panel.border=element_blank(),
        plot.background = element_blank(),
        panel.grid= element_blank(),axis.line = element_line(color = 'black')
  )+
  theme(legend.position="none")


# 
# 
ggsave(paste("spatial_corr_",unique(dat_spatial$annee)[1],"_",unique(dat_spatial$annee)[2],".pdf",sep=""))



#### IV at village level

data_snails_pond_v       <- left_join(dat_spatial,snails_pond,by=c("vid"))
data_snails_pond_v2      <- left_join(dat_spatial,snails_pond_0911,by=c("vid"))
data_snails_pond_biom_v  <- left_join(dat_spatial,snails_pond_biom,by=c("vid"))

data_snails_pond_v2$dry    <- ifelse(data_snails_pond_v2$annee==2009,data_snails_pond_v2$dry2009,data_snails_pond_v2$dry2011)
data_snails_pond_v2$winter <- ifelse(data_snails_pond_v2$annee==2009,data_snails_pond_v2$winte2009,data_snails_pond_v2$winter2011)
data_snails_pond_v2$rainy  <- ifelse(data_snails_pond_v2$annee==2009,data_snails_pond_v2$rainy2009,data_snails_pond_v2$rainy2011)
data_snails_pond_v2$yearly <- (data_snails_pond_v2$rainy+data_snails_pond_v2$winter+data_snails_pond_v2$dry)/3


form_firstst_v             <-  as.formula(paste(int,"~",paste(instr,collapse="+"),"+",
                                                paste(colnames(dat_ihs_v)[-c(1,2)],collapse="+"),
                                                "+I(l_surf_mean^2)+malaria+factor(annee)",sep=""))  

form_secondst_v             <-  as.formula(paste("l_rend_avg ~ vil_st1 +",
                                                 paste(colnames(dat_ihs_v)[-c(1,2)],collapse="+"),
                                                 "+I(l_surf_mean^2)+malaria+factor(annee)",sep=""))  

summary(felm(l_rend_avg~ + l_surf_mean + l_dspell_mean + 
               l_precip_wet + l_temp_mean + l_dspell_max + l_temp_night + 
               l_temp_day + l_totanim + l_nmem + l_nyoung + l_nchildren + 
               l_nwomen + l_agemoyen + l_bovins + l_equins + l_ovins + l_caprins + 
               l_dons + l_vols + l_NDVI_yr_mean + l_EVI_yr_mean + l_lstN_yr_mean + 
               l_lstD_yr_mean + l_tot_anim_avg + l_npk_avg + l_uree_avg + 
               l_phosph_avg + l_pest_solide_avg + l_pest_liquide_avg + l_herbicide_cl_avg + 
               l_fungicide_g_avg + l_fungicide_cl_avg + l_rodenticide_cl_avg + 
               l_rodenticide_g_avg + I(l_surf_mean^2) + malaria  | annee + region| 
               (intensity~ dry+yearly+winter+rainy) | 0,data_snails_pond_v2))




summary(felm(l_rend_avg~l_dspell_mean+malaria+l_precip_wet+l_temp_mean+l_dspell_max+l_temp_night+l_temp_day+l_totanim+l_nmem+l_nyoung+l_nchildren+l_nwomen+l_agemoyen+l_bovins+l_equins+l_ovins+l_caprins+l_dons+l_vols+l_NDVI_yr_mean+l_EVI_yr_mean+l_lstN_yr_mean+l_lstD_yr_mean+l_tot_anim_avg+l_npk_avg+l_uree_avg+l_phosph_avg+l_pest_solide_avg+l_pest_liquide_avg+l_herbicide_cl_avg+l_fungicide_g_avg+l_fungicide_cl_avg+l_rodenticide_cl_avg+l_rodenticide_g_avg| 
               annee+region| (intensity~dry:l_temp_mean +winter:l_temp_mean +rainy:l_temp_mean +
                                yearly:l_temp_mean +mean:l_temp_mean ) | vid, 
             data=data_snails_pond_biom_v))



data_snails_pond_biom_v$vil_st1 <- fitted(lm(form_firstst_v,data=data_snails_pond_biom_v))

second_vil                 <- lm(form_secondst_v,data=data_snails_pond_biom_v)

# est_iv_vil <- coeftest(second_vil,vcov=vcovBS,cluster = ~vid,type="wild",R=5000,cores=4)
est_iv_vil <- coeftest(second_vil,vcov=vcovCL,cluster = ~vid)
est_iv_vil[2,]

#### check with felm (no bootstrap) for consistency
# form_feiv_v <- as.formula(paste("l_rend_avg~",paste(colnames(dat_ihs_v)[-c(1,2)],collapse="+"),"+malaria | 
#                               annee | (intensity~rainy+dry+mean) | vid" ))
# feiv_v <- felm(form_feiv_v,data=data_snails_pond_biom_v)


res[10,c(1:3)] <- est_iv_vil[2,c(1,2,4)]
res[10,4]      <- res[10,1]*mean(datause$intensity)
res[10,5]      <- res[10,1]*quantile(datause$intensity,0.95)
res[10,6]      <- error(second_vil$fitted.values,dat_spatial$l_rend_avg)$err



# data_snails_pond_v$factor_years      <- ifelse(data_snails_pond_v$annee==unique(data_snails_pond_v$annee)[1],0,1)
# data_snails_pond_biom_v$factor_years <- data_snails_pond_v$factor_years 
# x_v <- c(colnames(dat_ihs_v)[-c(1,2)],"factor_years")
# 
# est_IV_DML_v <- DML_Forest_IV(y="l_rend_avg",d=int,x_v,
#                              instr=instr,tune=F,data=data_snails_pond_biom_v,boot=T,id="annee",cl="vid")
# # 
#  est_IV_DML_v
# # 

##### wrap up results


rownames(res) <- c("Plot level, time fe, full controls", "Plot level, time+hh fe, full controls", 
                   "Plot level, 2SLS, time fe","Plot level, het-IV",
                   paste("Plot level, best 4-fold Double ML, time+hh fe, best method: ",names(best_mid),sep=""),
                   "Plot level, 2SLS+Double ML, 2-fold random forest, time+hh fe",
                   "Village level, time+commune fe","Village level, SAR","Village level, spatial re",
                   "Village level, 2SLS, time fe")

rownames(res_sc) <- c("Plot level","Plot level,\n ML-flexible")


res                  <- round(res,4)
res_sc               <- round(res_sc,4)
res[,4:5]            <- 100*res[,4:5]
res[,4:5]            <- round(res[,4:5],2)
res_sc[,4:5]         <- 100*res_sc[,4:5]
res_sc[,4:5]         <- round(res_sc[,4:5],2)
results$res          <- res
last(names(results)) <- paste("res",paste(unique(datause$annee),collapse="_"),sep="_")



##### plot results

#change if you want to plot the full results (res) or the one with the main 3 (res_sc)
res_inv  <- apply(res,2,rev)
# res_inv  <- apply(res_sc,2,rev)


cint_l   <- res_inv[,1] - 1.96*res_inv[,2]
cint_h   <- res_inv[,1] + 1.96*res_inv[,2]
tot_eff  <- paste(as.character(res_inv[,4])," (",as.character(res_inv[,5]),")",sep="")

res_plot           <- cbind(res_inv,cint_l,cint_h)
ticks              <- rev(c("Plot level", "Plot level, hh fe",
                            "Plot level,\n 2SLS","Plot level,\n het-IV","Plot level, DML,\n GBM, hh fe",
                            "Plot level, 2SLS+DML, \n Random Forests, hh fe",
                            "Village level, \n commune fe","Village level, \n SAR",
                            "Village level, \nspatial re",
                            "Village level, \n2SLS"))
# ticks              <- rev(c("Plot level", "Plot level,\n ML-flexible"))
plot               <- broom::tidy(res_plot)

# require("ggrepel")

ggplot(plot,aes(Estimate,factor(ticks,levels=ticks),stat="identity")) + 
  geom_point(size=1.9) +
  geom_errorbarh(aes(xmax = cint_l, xmin = cint_h),height=0.7,linetype=2,alpha=09,colour="black",size=0.3)+
  geom_vline(xintercept=0,linetype="solid",colour="red")+
  geom_label(aes(label=tot_eff,color=p.val),size=3)+
  xlab("Estimate") + 
  ylab("")+
  xlim(min(cint_l),max(cint_h)+0.002)  +
  theme(axis.text.y = element_text(lineheight=1,face="bold", color="black",size=, angle=0),
        legend.key.size = unit(0.3, "cm"),panel.border = element_rect(colour = "black", fill=NA, size=5))+
  scale_y_discrete(expand=c(0.1,0))+
  #uncomment for full results matrix 
  scale_color_gradient(name="p.val",low="darkblue",high="red",space="Lab")+
  theme_classic2()
ggsave("est_full.pdf", width = 16, height = 9, units="cm")





######## interaction terms: channels through which schisto operates


# inter_vil_live <- paste("prev",X1,sep=":")
# inter_vil_hh <- paste("prev",X2,sep=":")
# inter_vil_fact <- paste("prev",X3,sep=":")
# inter_vil_clim <- paste("prev",X4,sep=":")
# inter_vil <- c(inter_vil_live,inter_vil_hh,inter_vil_fact,inter_vil_clim)
# inter_vil_2 <- c("prev:nwomen","prev:nchildren","prev:arboriculture",
#                "prev:autocons")
# 
# inter_vil_form <- as.formula(paste("lrend ~ ", paste(c(inter_vil,p,X1,X2,X3,X4),collapse="+"),"
#                                    | province+annee | 0 | vid",sep=""))
# 
# 
# est_inter_vil <- felm(inter_vil_form,data=datause)
# 
# summary(est_inter_vil)


# use intensity

int_vil_live_int <- paste("intensity",X1,sep=":")
int_vil_hh_int <- paste("intensity",X2,sep=":")
int_vil_fact_int <- paste("intensity",X3,sep=":")
int_vil_clim_int <- paste("intensity",X4,sep=":")
int_vil_int <- c(int_vil_live_int,int_vil_hh_int,int_vil_fact_int,int_vil_clim_int)



int_vil_form_int <- as.formula(paste("lrend ~ ", paste(c(int_vil_int,int,X1,X2,X3,X4),collapse="+")," | mid+annee | 0 | vid",sep=""))


est_int_vil_int <- felm(int_vil_form_int,data=datause)

summary(est_int_vil_int)
stargazer(est_int_vil_int,p.auto=T)
int_terms_sign <- grep("intensity:",names(est_int_vil_int$coefficients[which(est_int_vil_int$cpval<0.05),]))
stargazer(est_int_vil_int$coefficients[which(est_int_vil_int$cpval<=0.05),][int_terms_sign],title=paste(paste(unique(annee),collapse="+"),"Interaction Terms"), flip=T)


######  estimate impact of schisto on chosen inputs and check for non-significance
coefinput <- matrix(0,length(X1),2)

for (i in 1:length(X1)){
  forminput <- as.formula(paste(X1_log[i],"~", paste(paste(X1[-i],collapse="+"),paste(X2,collapse="+"),
                                                     paste(X3,collapse="+"),paste(X4_log,collapse="+"),sep="+"),"| annee | (intensity~rainy+dry+mean) | vid"))
  
  coefinput[i,] <- last(summary(felm(forminput,data=data_snails_pond_biom))$coefficients)[c(1,2)]
  rownames(coefinput) <- X1
}

coefinput <- round(coefinput,4)
coefinput[,2] <- paste("(",coefinput[,2],")",sep="")
stargazer(coefinput)

coefinput2 <- matrix(0,length(X2),2)

for (i in 1:length(X2)){
  forminput <- as.formula(paste(X2[i],"~", paste(paste(X1[-i],collapse="+"),paste(X2[-i],collapse="+"),
                                                 paste(X3,collapse="+"),paste(X4,collapse="+"),sep="+"),"| annee | (intensity ~ rainy+dry+mean ) | vid"))
  
  coefinput2[i,] <- summary(felm(forminput,data=data_snails_pond_biom))$coefficients[1,c(1,2)]
  rownames(coefinput2) <- X2
}

coefinput2 <- round(coefinput2,4)
coefinput2[,2] <- paste("(",coefinput2[,2],")",sep="")
stargazer(coefinput2)





######## crop choice: focus on cash crop (cotton)



res_c           <- matrix(0,5,5)
colnames(res_c) <- c("Estimate","S.E.","p-val","Mean Effect (%)","Top 5% Effect (%)")



### for 2009+2011 there is only 1100 as cotton (no bio/GM)
datause_cotton <- datause[which(datause$cult1==1100 | datause$cult2==1100 ),]
datause_food <- datause[-which(datause$cult1==1100 | datause$cult2==1100 ),]

wdata_mid_cotton <- wdata_mid[which(datause$cult1==1100 | datause$cult2==1100 ),]

prev_cotton    <- data.frame(lsurf=datause$lsurf,label=rep("All crops",dim(datause)[1]))
prev_nocotton  <- data.frame(lsurf=datause_cotton$lsurf,label=rep('Cotton',dim(datause_cotton)[1]))
prev_cotton_pl <- rbind(prev_cotton,prev_nocotton)

ggplot(data=prev_cotton_pl, aes(x=lsurf,fill=label)) + 
  geom_density(alpha=0.5,adjust=1) +
  scale_fill_manual(values = c("red", "lightblue"))+theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid= element_blank(),axis.line = element_line(color = 'black'))+xlab("")+ylab("")+
  xlab("Log Plot Surface")+
  ylab("Density")

ggsave("density_plot_cotton.pdf", width = 13, height = 10, units="cm")



est_fe_c <- felm(lrend ~ intensity + lsurf + lsurfsq + remune + entraid + gestion + 
                   parc_recup + mode_labour + l_tot_anim_traction + l_males_elev + 
                   l_fem_elev + l_totanim + l_morts + l_naissances + l_bovins + 
                   l_equins + l_ovins + l_caprins + l_dons + l_vols + l_autocons + 
                   l_ventebovins + l_achatbovins + l_venteovins + l_achatovins + 
                   l_ventecaprins + l_achatcaprins + l_duree_jachere + l_npk_kg + 
                   l_uree_kg + l_phosph_kg + l_pest_solide_g + l_pest_liquide_cl + 
                   l_herbicide_g + l_herbicide_cl + l_fungicide_g + l_fungicide_cl + 
                   l_rodenticide_g + l_rodenticide_cl + factor_localisation_1 + 
                   factor_localisation_2 + factor_localisation_3 + factor_relief_1 + 
                   factor_relief_2 + factor_antierosif_1 + factor_antierosif_2 + 
                   factor_antierosif_3 + factor_antierosif_4 + factor_antierosif_5 + 
                   factor_antierosif_6 + factor_antierosif_7 + factor_antierosif_8 + 
                   factor_mode_acq_1 + factor_mode_acq_2 + factor_mode_acq_3 + 
                   factor_mode_acq_4 + factor_mode_acq_5 + factor_niv_acq_1 + 
                   factor_niv_acq_2 + factor_niv_acq_3 + factor_niv_acq_4 + 
                   factor_niv_acq_5 + factor_niv_acq_6 + factor_niv_acq_7 +
                   factor_type_labour_1 + factor_type_labour_2 + factor_type_labour_3 + 
                   factor_educ_max_1 + factor_educ_max_2 + factor_educ_max_3 + 
                   factor_educ_max_4 + factor_educ_max_5 + factor_educ_max_6 + 
                   factor_perte_1 + factor_perte_4 + factor_perte_6+factor_perte_8+nmem + agemoyen + age_chefmen + 
                   numfamilies + nchildren + nyoung + nwomen + culture_pluviale + 
                   culture_maraichere + arboriculture + autre_contresaison + 
                   peche + l_temp_mean + l_temp_day + l_temp_night + l_precip_wet + 
                   l_temp_dry + l_dspell_mean + l_dspell_max +l_lstD_yr_mean + 
                   l_lstN_yr_mean + l_EVI_yr_mean + l_NDVI_yr_mean
                 | annee+mid| 0 | vid,data=datause_cotton)


res_c[1,1:3] <- summary(est_fe_c)$coefficients[1,c(1,2,4)]
res_c[1,4] <- res_c[1,1]*mean(datause_cotton$intensity)
res_c[1,5] <- res_c[1,1]*quantile(datause_cotton$intensity,0.95)

# IV
datause_cotton_biom <- data_snails_pond_biom[which(datause$cult1==1100 | datause$cult2==1100 ),]

form_1stage_tfe_cotton <-  as.formula(
  paste(paste(int,"~ ",paste(paste(instr,collapse="+"),paste(x_log,collapse="+"),sep="+")),sep=""))

form_2stage_tfe_cotton <-  as.formula(
  paste(paste(y,"~ ",paste("stage1_tfe",paste(x_log,collapse="+"),sep="+")),sep=""))


datause_cotton_biom$stage1_tfe <- fitted(lm(form_1stage_tfe_cotton,data=datause_cotton_biom))

second_tfe_cotton              <- lm(form_2stage_tfe_cotton,data=datause_cotton_biom)

set.seed(42)
est_iv_cotton  <- coeftest(second_tfe_cotton,vcov=vcovCL,cluster = ~ vid)

res_c[2,1:3] <- est_iv_cotton[2,c(1,2,4)]
res_c[2,4] <- res_c[2,1]*mean(datause_cotton$intensity)
res_c[2,5] <- res_c[2,1]*quantile(datause_cotton$intensity,0.95)


K            <- 3
trim         <- c(0,1)
set.seed(42)
# watch out for cv.folds option - sometimes it breaks
Boosting     <- list(bag.fraction = .5, train.fraction = 1.0, interaction.depth=2, n.trees=1000, shrinkage=.01, n.cores=1, cv.folds=1, verbose = FALSE, clas_dist= 'adaboost', reg_dist='gaussian')
Forest       <- list(clas_nodesize=1, reg_nodesize=5, ntree=2000, na.action=na.omit, replace=TRUE)
RLasso       <- list(penalty = list(homoscedastic = FALSE, X.dependent.lambda =FALSE, lambda.start = NULL, c = 1.1), intercept = TRUE)
Nnet         <- list(size=2,  maxit=1000, decay=0.01, MaxNWts=5000,  trace=FALSE)
Trees        <- list(reg_method="anova", clas_method="class")
arguments    <- list(Boosting=Boosting, Forest=Forest, RLasso=RLasso, Nnet=Nnet, Trees=Trees)

methods_c       <- c("Boosting","RLasso","Nnet")

nnum_c          <- sapply(datause_cotton,is.numeric)
datause_ML_cotton    <- datause_cotton[,nnum_c]
x.flex_t      <- c(x.flex,"factor_years_1")

est_full_DML_cotton               <- Full_DML(datause_ML_cotton, y, int, x_t, x.flex_t, methods=methods_c, nfold=K, est="plinear", arguments=arguments, ensemble=ensemble, silent=FALSE, trim)
est_full_DML_mid_cotton           <- Full_DML(wdata_mid_cotton, y, int, x, x.flex, methods=methods_c, nfold=K, est="plinear", arguments=arguments, ensemble=ensemble, silent=FALSE, trim)
# est_full_DML_mid_f  <- Full_DML(wdata_mid, y, int, x[1:2], x.flex, methods=methods_f, nfold=K, est="plinear", arguments=arguments, ensemble=ensemble, silent=FALSE, trim)

pvals_DML_c          <- 2*pt(-abs(est_full_DML_cotton[1,c(1:(dim(est_full_DML_cotton)[2]-1))]/est_full_DML_cotton[2,c(1:(dim(est_full_DML_cotton)[2]-1))]), df = nrow(datause_cotton)-2)
pvals_DML_mid_c      <- 2*pt(-abs(est_full_DML_mid_cotton[1,c(1:(dim(est_full_DML_mid_cotton)[2]-1))]/est_full_DML_mid_cotton[2,c(1:(dim(est_full_DML_mid_cotton)[2]-1))]), df = nrow(datause_cotton)-2 - length(unique(datause_cotton$mid)))
# pvals_DML_mid_f      <- 2*pt(-abs(est_full_DML_mid_f[1,c(1:(dim(est_full_DML_mid_f)[2]-1))]/est_full_DML_mid_f[2,c(1:(dim(est_full_DML_mid_f)[2]-1))]), df = nrow(datause)-2 - length(unique(mid)))

best_c                <- which(est_full_DML_cotton[3,1:3]==min(est_full_DML_cotton[3,1:3]))
best_mid_c            <- which(est_full_DML_mid_cotton[3,1:3]==min(est_full_DML_mid_cotton[3,1:3]))



res_c[3,1:2] <- est_full_DML_cotton[1:2,best_c]
res_c[3,3]   <- pvals_DML_c[best_c]
res_c[3,4]   <- res_c[3,1]*mean(datause_cotton$intensity)
res_c[3,5]   <- res_c[3,1]*quantile(datause_cotton$intensity,0.95)

res_c[4,1:2] <- est_full_DML_mid_cotton[1:2,best_mid_c]
res_c[4,3]   <- pvals_DML_mid_c[best_mid_c]
res_c[4,4]   <- res_c[4,1]*mean(datause_cotton$intensity)
res_c[4,5]   <- res_c[4,1]*quantile(datause_cotton$intensity,0.95)


est_ML_base_int_mid_IV_c <- DML_Forest_IV(y,d=int,x_t,instr=instr,tune=F,data=datause_cotton_biom,boot=T,id="annee",cl="vid")



res_c[5,1]   <- est_ML_base_int_mid_IV_c$theta_2sls
res_c[5,2]   <- est_ML_base_int_mid_IV_c$theta_2sls.se
res_c[5,3]   <- est_ML_base_int_mid_IV_c$theta_2sls.p
res_c[5,4]   <- res_c[5,1]*mean(datause_cotton$intensity)
res_c[5,5]   <- res_c[5,1]*quantile(datause_cotton$intensity,0.95)




res_c                  <- round(res_c,4)
res_c[,4:5]            <- 100*res_c[,4:5]
res_c[,4:5]            <- round(res_c[,4:5],2)


res_inv_c  <- apply(res_c,2,rev)
cint_l_c   <- res_inv_c[,1] - 1.96*res_inv_c[,2]
cint_h_c   <- res_inv_c[,1] + 1.96*res_inv_c[,2]
tot_eff_c   <- as.character(res_inv_c[,4])

res_plot_c           <- cbind(res_inv_c,cint_l_c,cint_h_c)
ticks_c              <- rev(c("Linear, hh fe",
                              "Linear, 2SLS","DML, Neural Nets","DML, Forest+Lasso, hh fe",
                              "DML, 2SLS+Forest"))
plot_c               <- broom::tidy(res_plot_c)



# require("ggrepel")
ggplot(plot_c,aes(Estimate,factor(ticks_c,levels=ticks_c),stat="identity")) + 
  geom_point(size=1.9) +
  geom_errorbarh(aes(xmax = cint_l_c, xmin = cint_h_c),height=0.4,linetype=2,alpha=09,colour="black",size=0.3)+
  geom_vline(xintercept=0,linetype="solid",colour="red")+
  geom_label(aes(label=tot_eff_c,color=p.val),size=3)+
  xlab("Estimate") + 
  ylab("")+
  xlim(min(cint_l_c),max(cint_h_c))  +
  theme(axis.text.y = element_text(lineheight=1,face="bold", color="black",size=, angle=0),
        legend.key.size = unit(0.3, "cm"),panel.border = element_rect(colour = "black", fill=NA, size=5))+
  scale_y_discrete(expand=c(0.1,0))+
  scale_color_gradient(name="p.val",low="red1",high="red4",space="Lab")+
  theme_classic2()
ggsave("est_full_cotton.pdf",width = 16, height = 8, units="cm")




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

dam_provinces   <- c(4,7,14,18,30,35,40) 
lake_provinces  <- c(1,6)
# only focus on dams


# data_dam <- data_damuse[data_damuse$province %in% c(lake_provinces,dam_provinces) | data_damuse$annee == yr1,]
data_dam <- data_damuse[data_damuse$province %in% c(dam_provinces),]
# data_nodam <- data_damuse[!data_damuse$province %in% c(lake_provinces,dam_provinces)| data_damuse$annee == yr1,]
data_nodam <- data_damuse[!data_damuse$province %in% c(dam_provinces),]
# data_damuse$dam <- ifelse(data_damuse$province %in% c(dam_provinces,lake_provinces),1,0)
data_damuse$dam           <- ifelse(data_damuse$province %in% c(dam_provinces),1,0)
dat_spatial$dam           <- ifelse(dat_spatial$province %in% c(dam_provinces[-1]),1,0)


#### plot disease densities per year 

prev_dens1 <- data.frame(prev=data_dam$prev,label=rep('Near a dam',dim(data_dam)[1]))
prev_dens2 <- data.frame(prev=data_damuse$prev,label=rep('Overall',dim(data_damuse)[1]))
prev_dens3 <- data.frame(prev=data_nodam$prev,label=rep('Remaining Areas',dim(data_nodam)[1]))

prev_dens = rbind(prev_dens1,prev_dens2)
prev_dens_nodam = rbind(prev_dens1,prev_dens3)
prev_dens_all = rbind(prev_dens1,prev_dens2,prev_dens3)

ggplot(data=prev_dens, aes(x=prev,fill=label)) + 
  geom_density(alpha=0.5,adjust=2) +
  geom_vline(xintercept=mean(data_dam$prev),color="red", linetype="dashed", size=1)+
  geom_vline(xintercept=mean(datause$prev),color="lightblue", linetype="dashed", size=1)+
  scale_fill_manual(values = c("red", "lightblue"))+
  xlab("Prevalence")+
  ylab("Density")
ggsave("density_dam.pdf", width = 13, height = 10, units="cm")


ggplot(data=prev_dens_nodam, aes(x=prev,fill=label)) + 
  geom_density(alpha=0.5,adjust=2) +
  geom_vline(xintercept=mean(data_dam$prev),color="red", linetype="dashed", size=1)+
  geom_vline(xintercept=mean(data_nodam$prev),color="green", linetype="dashed", size=1)+
  scale_fill_manual(values = c("red", "green"))+
  xlab("Prevalence")+
  ylab("Density")
ggsave("density_damvsnodam.pdf", width = 13, height = 10, units="cm")

ggplot(data=prev_dens_all, aes(x=prev,fill=label)) + 
  geom_density(alpha=0.5,adjust=2) +
  geom_vline(xintercept=mean(data_dam$prev),color="red", linetype="dashed", size=1)+
  geom_vline(xintercept=mean(data_nodam$prev),color="green", linetype="dashed", size=1)+
  geom_vline(xintercept=mean(datause$prev),color="lightblue", linetype="dashed", size=1)+
  scale_fill_manual(values = c("red", "green","lightblue"),name=NULL)+
  xlab("Prevalence")+
  ylab("Density")
ggsave("density_damvsnodam_all.pdf", width = 13, height = 10, units="cm")

##### show dam interaction with ols and region fe:


res_dam <- matrix(0,4,12)
rownames(res_dam) <- c("Plot level, region fe", "Plot level, 2SLS", 
                       "Village level, region fe", "Village level, spatial fe")

colnames(res_dam) <- c("Interaction","s.e.inter","p.inter","Intensity","s.e.int","p.int", "Dam","s.e.dam","p.dam","MSE","mean_interaction","mean_intensity")

x            <- c(X1_log, X2, X3_t[-1],X4_log,dis)
x_t          <- c(X1_log, X2, X3_t,X4_log,dis)
# x        <- c(X1, X2,X3[1:92],X4) 
# x          <- c(X4)

form_dam_int_reg <- as.formula(paste(y,"~",paste(paste(x,collapse="+"),paste(int,"dam+intensity:dam",sep="+"),
                                                 sep="+"), " |region+annee | 0 |region",sep=""))

est_dam_base      <- felm(lrend~intensity+dam+intensity:dam | annee | 0 |vid,data=data_damuse)
est_dam_int       <- lm(as.formula(paste(y,"~",paste(int, paste(x_t,collapse="+"),"dam+intensity:dam",sep="+"),sep="")),data=data_damuse)
# est_dam_int_hhfe  <- felm(form_dam_int_hh,data=data_damuse)
est_dam_int_regfe <- felm(form_dam_int_reg,data=data_damuse)

mse.dam_fe <- error(est_dam_int_regfe$fitted.values, data_damuse$lrend)

res_dam[1,1:3] <- summary(est_dam_int_regfe)$coefficients[dim(summary(est_dam_int_regfe)$coefficients)[1],][c(1,2,4)]
res_dam[1,4:6] <- summary(est_dam_int_regfe)$coefficients[dim(summary(est_dam_int_regfe)$coefficients)[1]-2,][c(1,2,4)]
res_dam[1,7:9] <- summary(est_dam_int_regfe)$coefficients[dim(summary(est_dam_int_regfe)$coefficients)[1]-1,][c(1,2,4)]
res_dam[1,10]  <- mse.dam_fe$err
res_dam[1,11]  <- res_dam[1,1]*mean(data_damuse$intensity[which(data_damuse$dam==1)])
res_dam[1,12]  <- res_dam[1,4]*mean(data_damuse$intensity)





instr = c("rainy", "dry"  ,"mean" )
form_1stage_tfe_dam <-  as.formula(
  paste(paste(int,"~ ",paste(paste(instr,collapse="+"),
                             paste(instr,":dam",collapse="+"),
                             "dam",paste(x_t,collapse="+"),sep="+")),sep=""))

form_2stage_tfe_dam <-  as.formula(
  paste(paste(y,"~ ",paste("stage1_tfe_dam + stage1_tfe_dam:dam + dam",paste(x_t,collapse="+"),sep="+")),sep=""))



data_snails_pond_biom$stage1_tfe_dam <- fitted(lm(form_1stage_tfe_dam,data=data_snails_pond_biom))
second_tfe_dam                       <- lm(form_2stage_tfe_dam,data=data_snails_pond_biom)



# est_iv_dam_cl  <- coeftest(est_iv_dam ,vcov=vcovBS,cluster = ~ vid,type="wild",R=3000)
est_iv_dam_cl <- coeftest(second_tfe_dam,vcov=vcovCL,cluster = ~ vid)





res_dam[2,1:3] <- est_iv_dam_cl[dim(est_iv_dam_cl)[1],][c(1,2,4)]
res_dam[2,4:6] <- est_iv_dam_cl[2,][c(1,2,4)]
res_dam[2,7:9] <- est_iv_dam_cl[3,][c(1,2,4)]
res_dam[2,10]  <- mse.dam_hetiv$err
res_dam[2,11]  <- res_dam[2,1]*mean(data_damuse$intensity[which(data_damuse$dam==1)])
res_dam[2,12]  <- res_dam[2,4]*mean(data_damuse$intensity)


############ village aggregation + dams





est_vil_dam <- felm(l_rend_avg ~ intensity + dam + intensity:dam+l_surf_mean +
                      l_dspell_mean + l_precip_wet + 
                      l_temp_mean + l_dspell_max + l_temp_night + l_temp_day + 
                      l_totanim + l_nmem + l_nyoung + l_nchildren + l_nwomen + 
                      l_agemoyen + l_nchildren + l_bovins + l_equins + l_ovins + 
                      l_caprins + l_NDVI_yr_mean + l_EVI_yr_mean + 
                      l_lstN_yr_mean + l_lstD_yr_mean + l_tot_anim_avg+l_npk_avg+ l_uree_avg + 
                      l_phosph_avg + l_pest_solide_avg + l_pest_liquide_avg + l_herbicide_cl_avg + 
                      l_fungicide_g_avg + l_fungicide_cl_avg + l_rodenticide_cl_avg + 
                      l_rodenticide_g_avg+malaria | annee +region| 0 |region, 
                    data=dat_spatial)

est_moran <- lm(l_rend_avg ~ intensity + dam + intensity:dam+l_surf_mean +
                  l_dspell_mean + l_precip_wet + 
                  l_temp_mean + l_dspell_max + l_temp_night + l_temp_day + 
                  l_totanim + l_nmem + l_nyoung + l_nchildren + l_nwomen + 
                  l_agemoyen + l_nchildren + l_bovins + l_equins + l_ovins + 
                  l_caprins + l_NDVI_yr_mean + l_EVI_yr_mean + 
                  l_lstN_yr_mean + l_lstD_yr_mean + l_tot_anim_avg+l_npk_avg+ l_uree_avg + 
                  l_phosph_avg + l_pest_solide_avg + l_pest_liquide_avg + l_herbicide_cl_avg + 
                  l_fungicide_g_avg + l_fungicide_cl_avg + l_rodenticide_cl_avg + 
                  l_rodenticide_g_avg+malaria+factor(region)+factor(annee), 
                data=dat_spatial)

mse.dam_fe_v <- error(est_vil_dam$fitted.values, dat_spatial$l_rend_avg)


m1_dam <- lm.morantest(est_moran, W, alternative="two.sided")

res_dam[3,1:3] <- summary(est_vil_dam)$coefficients[dim(summary(est_vil_dam)$coefficients)[1],][c(1,2,4)]
res_dam[3,4:6] <- summary(est_vil_dam)$coefficients[1,][c(1,2,4)]
res_dam[3,7:9] <- summary(est_vil_dam)$coefficients[2,][c(1,2,4)]
res_dam[3,10]  <- mse.dam_fe_v$err
res_dam[3,11]  <- res_dam[3,1]*mean(dat_spatial$intensity[which(dat_spatial$dam==1)])
res_dam[3,12]  <- res_dam[3,4]*mean(dat_spatial$intensity)

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
summary(est_vil_dam_splag,hausman=T)


res_dam[4,1]   <- est_vil_dam_splag$coefficients[grep("intensity:dam",names(est_vil_dam_splag$coefficients))[1]]
res_dam[4,2]   <- est_vil_dam_splag$rest.se[grep("intensity:dam",names(est_vil_dam_splag$coefficients))[1]]
res_dam[4,3]   <- 2*pt(-abs(res_dam[4,1] /res_dam[4,2]), df = nrow(dat_spatial)-1)


res_dam[4,4]   <- est_vil_dam_splag$coefficients[grep("intensity",names(est_vil_dam_splag$coefficients))[1]]
res_dam[4,5]   <- est_vil_dam_splag$rest.se[grep("intensity",names(est_vil_dam_splag$coefficients))[1]]
res_dam[4,6]   <- 2*pt(-abs(res_dam[4,4] /res_dam[4,5]), df = nrow(dat_spatial)-1)


res_dam[4,7]   <- est_vil_dam_splag$coefficients[grep("dam",names(est_vil_dam_splag$coefficients))[1]]
res_dam[4,8]   <- est_vil_dam_splag$rest.se[grep("dam",names(est_vil_dam_splag$coefficients))[1]]
res_dam[4,9]   <- 2*pt(-abs(res_dam[4,7] /res_dam[4,8]), df = nrow(dat_spatial)-1)


res_dam[4,10]  <- mse.dam_splag_v$err
res_dam[4,11]  <- res_dam[4,1]*mean(dat_spatial$intensity[which(dat_spatial$dam==1)])
res_dam[4,12]  <- res_dam[4,4]*mean(dat_spatial$intensity)


#### plot results of schisto and water resources

res_dam_r <- round(res_dam,3)
res_dam_r[,11] <- 100*res_dam_r[,11]
res_dam_r[,12] <- 100*res_dam_r[,12]
res_dam_r[,11:12] <- round(res_dam_r[,11:12],1)
ticks_dam        <- c("Plot level,\n region fe", "Plot level,\n 2SLS","Village level,\n region fe", "Village level,\n SAR fe")
plot_dam              <- broom::tidy(res_dam_r)

# require("ggrepel")



ggplot(plot_dam,aes(Interaction,ticks_dam))+
  geom_errorbarh(aes(xmax = Interaction+1.96*s.e.inter, xmin =Interaction-1.96*s.e.inter),height=0.3,
                 linetype=1,alpha=09,colour="black",size=0.3)+
  geom_vline(xintercept=0,linetype="solid",colour="red")+
  geom_errorbarh(aes(xmax = Intensity+1.96*s.e.int, xmin = Intensity-1.96*s.e.int),height=0.35,
                 linetype=2,alpha=09,colour="black",size=0.3)+  
  geom_label(aes(0.01,label=paste("Dam\n",100*Dam,sep=""),color=p.dam),size=3)+
  scale_color_gradient(position="bottom",name="p",low="blue",
                       high="red",space="Lab",breaks=c(0.01,0.1,0.2,0.3))+
  xlab("Estimate") + 
  ylab("")+
  geom_label(aes(Interaction,label=paste("I*Dam\n",mean_interaction,sep=""),color=p.inter),size=3)+
  geom_label(aes(Intensity,label=paste("I\n",mean_intensity,sep=""),color=p.int),size=3)+
  scale_x_continuous(breaks=c(-0.02,-0.01,0.0),limits=c(-0.03,0.011))+
  theme(axis.text.y = element_text(lineheight=1,face="bold", color="black", 
                                   size=9, angle=0),legend.key.size = unit(0.3, "cm"))+
  scale_y_discrete(expand=c(0,1))+
  scale_color_gradient(name="p.val",low="blue",high="red",space="Lab",breaks=c(0.5,0.3,0.1))+
  theme_classic2()

ggsave("est_dam.pdf",height=2.8)


plot_dam_small <- plot_dam[4,]
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


e# 
# 
# est_dam_spat <- fitme(l_rend_avg~ intensity:dam+intensity + dam + l_surf_mean + I(l_surf_mean^2) + l_dspell_mean + l_precip_wet +
#                         l_temp_mean + l_dspell_max + l_temp_night + l_temp_day +
#                         l_totanim + l_nmem + l_nyoung + l_nchildren + l_nwomen +
#                         l_agemoyen + l_nchildren + l_bovins + l_equins + l_ovins +
#                         l_caprins + l_dons + l_vols + l_NDVI_yr_mean + l_EVI_yr_mean +
#                         l_lstN_yr_mean + l_lstD_yr_mean + l_tot_anim_avg + l_npk_avg+l_uree_avg +
#                         l_phosph_avg + l_pest_solide_avg + l_pest_liquide_avg + l_herbicide_cl_avg +
#                         l_fungicide_g_avg + l_fungicide_cl_avg + l_rodenticide_cl_avg +
#                         l_rodenticide_g_avg +
#                         Matern(1 |lat + long), data = dat_spatial, family ="gaussian")
# 
# summary(est_dam_spat)[[3]][dim(summary(est_dam_spat)[[3]])[1],]


########### distance from dams 


res_dam_dist           <- matrix(0,2,12)
colnames(res_dam_dist) <- c("coefdistdams", "se","pval","impact","coefdams", "se","pval","impact","coefdistnodams", "se","pval","impact")
rownames(res_dam_dist) <- c("Region fe", "SAR")

dams_dist <- read.csv("data/dams_distance/data-dams_distance.csv")
dams_dist_use <- dams_dist[,c(5,8,9,10)]
dams_dist_use$regime_d <- ifelse(dams_dist_use$regime=="TEMPORAIRE",0,1)


data_snails_pond_biom_dist <- left_join(data_snails_pond_biom,dams_dist_use,by=c("vid"))
dat_spatial_dist           <- left_join(dat_spatial,dams_dist_use,by=c("vid"))


form_dam_int_reg_dist <- as.formula(paste(y,"~",paste(paste(x,collapse="+"),
                                                      paste(int,"dam+dam_dist+intensity:dam + intensity:dam_dist+intensity:dam:dam_dist",sep="+"),
                                                      sep="+"), " |region+annee | 0 |vid",sep=""))

est_dam_dist <- felm(form_dam_int_reg_dist ,data=data_snails_pond_biom_dist)
summary(est_dam_dist)

instr = c("rainy", "dry"  ,"mean" )
form_1stage_tfe_dam_dist <-  as.formula(
  paste(paste(int,"~ ",paste(paste(instr,collapse="+"),
                             paste(instr,":dam",collapse="+"),paste(instr,":dam_dist",collapse="+"),
                             paste(instr,":dam:dam_dist+dam:dam_dist",collapse="+"),
                             "dam+dam_dist",paste(x_t,collapse="+"),sep="+")),sep=""))

form_2stage_tfe_dam_dist <-  as.formula(
  paste(paste(y,"~ ",paste("stage1_tfe_dam + stage1_tfe_dam_dist:dam_dist:dam+ 
                           stage1_tfe_dam_dist:dam+stage1_tfe_dam_dist:dam_dist+dam:dam_dist+ 
                           dam_dist+dam+",paste(x_t,collapse="+"),sep="+")),sep=""))



data_snails_pond_biom_dist$stage1_tfe_dam_dist <- fitted(lm(form_1stage_tfe_dam_dist,data=data_snails_pond_biom_dist))
second_tfe_dam_dist                       <- lm(form_2stage_tfe_dam_dist,data=data_snails_pond_biom_dist)



# est_iv_dam_cl  <- coeftest(est_iv_dam ,vcov=vcovBS,cluster = ~ vid,type="wild",R=3000)
est_iv_dam_cl_dist <- coeftest(second_tfe_dam_dist,vcov=vcovCL,cluster = ~ vid)
est_iv_dam_cl_dist 



est_vil_dam_dist <- felm(l_rend_avg~intensity+dam+dam_dist+dam:dam_dist+intensity:dam + intensity:dam:dam_dist+ intensity:dam_dist+l_surf_mean + I(l_surf_mean^2) + l_dspell_mean + l_precip_wet +
                           l_temp_mean + l_dspell_max + l_dspell_mean+l_temp_night + l_temp_day +
                           l_totanim + l_nmem + l_nyoung + l_nchildren + l_nwomen +
                           l_agemoyen + l_nchildren + l_bovins + l_equins + l_ovins +
                           l_caprins + l_dons + l_vols + +malaria+l_NDVI_yr_mean + l_EVI_yr_mean +
                           l_lstN_yr_mean + l_lstD_yr_mean + l_tot_anim_avg + l_npk_avg+l_uree_avg +
                           l_phosph_avg + l_pest_solide_avg + l_pest_liquide_avg + l_herbicide_cl_avg +
                           l_fungicide_g_avg + l_fungicide_cl_avg + l_rodenticide_cl_avg +
                           l_rodenticide_g_avg +malaria| region+annee | 0 | region,data=dat_spatial_dist)
summary(est_vil_dam_dist)



res_dam_dist[1,c(1:3)]   <- last(summary(est_vil_dam_dist)$coefficients)[c(1,2,4)]
res_dam_dist[1,4]        <- mean(dat_spatial_dist$intensity[which(dat_spatial_dist$dam == 1)])*res_dam_dist[1,1]
res_dam_dist[1,c(5:7)]   <- summary(est_vil_dam_dist)$coefficients[dim(summary(est_vil_dam_dist)$coefficient)[1]-2,][c(1,2,4)]
res_dam_dist[1,8]        <- mean(dat_spatial_dist$intensity[which(dat_spatial_dist$dam == 1)])*res_dam_dist[1,5]
res_dam_dist[1,c(9:11)]   <- summary(est_vil_dam_dist)$coefficients[dim(summary(est_vil_dam_dist)$coefficient)[1]-1,][c(1,2,4)]
res_dam_dist[1,12]        <- mean(dat_spatial_dist$intensity[which(dat_spatial_dist$dam == 0)])*res_dam_dist[1,9]




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

res_dam_dist[2,1] <- est_vil_dam_splag_dist$coefficients[grep("intensity:dam:dam_dist",names(est_vil_dam_splag_dist$coefficients))[1]]
res_dam_dist[2,2] <- est_vil_dam_splag_dist$rest.se[grep("intensity:dam:dam_dist",names(est_vil_dam_splag_dist$coefficients))[1]]
res_dam_dist[2,3] <- 2*pt(-abs(res_dam_dist[2,1] /res_dam_dist[2,2]), df = dim(est_vil_dam_splag_dist$X)[1])
res_dam_dist[2,4] <- mean(dat_spatial_dist$intensity[which(dat_spatial_dist$dam == 1)])*res_dam_dist[2,1]
res_dam_dist[2,5] <- est_vil_dam_splag_dist$coefficients[grep("intensity:dam",names(est_vil_dam_splag_dist$coefficients))[1]]
res_dam_dist[2,6] <- est_vil_dam_splag_dist$rest.se[grep("intensity:dam",names(est_vil_dam_splag_dist$coefficients))[1]]
res_dam_dist[2,7] <- 2*pt(-abs(res_dam_dist[2,5] /res_dam_dist[2,6]), df = dim(est_vil_dam_splag_dist$X)[1])
res_dam_dist[2,8] <- mean(dat_spatial_dist$intensity[which(dat_spatial_dist$dam == 1)])*res_dam_dist[2,5]
res_dam_dist[2,9] <- est_vil_dam_splag_dist$coefficients[grep("intensity:dam_dist",names(est_vil_dam_splag_dist$coefficients))[1]]
res_dam_dist[2,10] <- est_vil_dam_splag_dist$rest.se[grep("intensity:dam_dist",names(est_vil_dam_splag_dist$coefficients))[1]]
res_dam_dist[2,11] <- 2*pt(-abs(res_dam_dist[2,9]/res_dam_dist[2,10]), df = dim(est_vil_dam_splag_dist$X)[1])
res_dam_dist[2,12] <- mean(dat_spatial_dist$intensity[which(dat_spatial_dist$dam == 0)])*res_dam_dist[2,9]

res_dam_dist_r           <- round(res_dam_dist,5)
res_dam_dist_r[,4*c(1:3)] <- round(100*res_dam_dist_r[,4*c(1:3)],2)

### take vil with dam=1, then take mean distance, then take mean yield among those at the mean (btw 50 and 51)

ticks_dam_dist            <- c("Region fe","SAR")
plot_dam_dist              <- broom::tidy(res_dam_dist_r)
plot_dam_dist_small       <- plot_dam_dist[2,]

ggplot(plot_dam_dist,aes(coefdistdams,ticks_dam_dist))+
  geom_vline(xintercept=0,linetype="solid",colour="red")+
  geom_errorbarh(aes(xmax = coefdistdams+1.96*se, xmin =coefdistdams-1.96*se),height=0.1,
                 linetype=1,alpha=09,colour="black",size=0.2)+
  # geom_errorbarh(aes(xmax = coefdams+1.96*se, xmin =coefdams-1.96*se),height=0.3,
  #                linetype=1,alpha=09,colour="black",size=0.3)+
  geom_errorbarh(aes(xmax = coefdistnodams+1.96*se, xmin =coefdistnodams-1.96*se),height=0.1,
                 linetype=3,alpha=09,colour="black",size=0.2)+
  geom_label(aes(coefdistdams,label=paste("I*Dam*Dist\n",paste(impact,"%"),sep="")),color="blue",size=3)+
  geom_label(aes(coefdistnodams,label=paste("I*Dist\n",paste(impact.2,"%"),sep="")),color="red",size=3)+
  geom_label(aes(-0.0025,label=paste("I*Dam\n",paste(impact.1,"%"),sep="")),color="blue",size=3)+
  xlab("Estimate") + 
  ylab("")+
  scale_x_continuous(breaks=c(-0.001,0.000,0.001,0.002,0.003),limits=c(-0.003,0.0033))+
  theme(axis.text.y = element_text(lineheight=1,face="bold", color="blue", 
                                   size=9, angle=0),legend.key.size = unit(0.3, "cm"))+
  theme_classic2()
ggsave("est_dam_dist.pdf",height=2)


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



####### non parametric


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
pdf("intensity_distance_dams_3d.pdf")
vis.gam(est_np_damdist,color="topo",view=c("intensity","dam_dist"),theta=30, 
        too.far=1,ticktype="detailed",main = " " ,xlab="Intensity",
        ylab="Distance",zlab = "Log Yield/ha",)
dev.off()
# 
# 
#   }
# }



############ summary statistics + plots + maps


prod_village_yr1 <- datause %>%  
  filter(annee==unique(annee)[1]) %>%
  group_by(vid) %>% summarise(n=n(),prev=mean(prev),rend = mean(rend),lat=mean(lat),long=mean(long),prev=unique(prev))


prod_village_yr2 <- datause %>%  
  filter(annee==unique(annee)[2]) %>%
  group_by(vid) %>% summarise(n=n(),prev=mean(prev),rend = mean(rend),lat=mean(lat),long=mean(long),prev=unique(prev))

# get the Burkina centroid and the google maps baselayer

options(repr.plot.width = 14, repr.plot.height = 10)
register_google(key="AIzaSyA6mK1VPdr0t1a1gdcgakENguRleIauiJU") # stupid google now wants an api key

bf_centroid     <- dbGetQuery(conn, "WITH cntrd as (SELECT ST_Transform(ST_Centroid(geom), 4326) as geom FROM admin.adm0) SELECT ST_X(geom), ST_Y(geom) FROM cntrd")
burkina.map_yr1    <- ggmap(get_googlemap(center = unlist(bf_centroid), scale = 4, zoom = 7), extent = "normal")
burkina.map_yr1 +  geom_point(aes(x = long, y = lat, color = rend), data = prod_village_yr1) +
  scale_color_gradient(low="yellow",high="darkgreen",name="Yield")
ggsave(paste("burkina_prod_map_",unique(annee)[1],".pdf",sep=""), width = 13, height = 10, units="cm")

burkina.map_yr2   <- ggmap(get_googlemap(center = unlist(bf_centroid), scale = 4, zoom = 7), extent = "normal")
burkina.map_yr2+  geom_point(aes(x = long, y = lat, color = rend), data = prod_village_yr2) +
  scale_color_gradient(low="yellow",high="darkgreen",name="Yield")
ggsave(paste("burkina_prod_map_",unique(annee)[2],".pdf",sep=""), width = 13, height = 10, units="cm")


burkina.map_yr1_prev    <- ggmap(get_googlemap(center = unlist(bf_centroid), scale = 4, zoom = 7), extent = "normal")
burkina.map_yr1_prev +  geom_point(aes(x = long, y = lat, color = prev), data = prod_village_yr1) +
  scale_color_gradient2(low = "darkred", high = "darkgreen", mid = "yellow", midpoint = mean(prev),name="Prevalence")
ggsave(paste("burkina_prev_map_",unique(annee)[1],".pdf",sep=""), width = 13, height = 10, units="cm")

burkina.map_yr2_prev   <- ggmap(get_googlemap(center = unlist(bf_centroid), scale = 4, zoom = 7), extent = "normal")
burkina.map_yr2_prev +  geom_point(aes(x = long, y = lat, color = prev), data = prod_village_yr2) +
  scale_color_gradient2(low = "darkred", high = "darkgreen", mid = "yellow", midpoint=mean(prev),name="Prevalence")
ggsave(paste("burkina_prev_map_",unique(annee)[2],".pdf",sep=""), width = 13, height = 10, units="cm")


## use raw dataset to use variable labels
require(foreign)
datraw <- read.spss("data/cahiers/C2_1993_2017.sav")
datraw <- data.frame(datraw)
datraw0911 <- datraw[which(datraw$ANNEE==yr1[1] | datraw$ANNEE==yr2[1]),]
cult1_r <- datraw0911$CULT1
cult2_r <- datraw0911$CULT2

cult1_r <- cult1_r[-which(cult1_r==9999)]
cult2_r <- cult2_r[-which(cult2_r==9999 | is.na(cult2_r) | cult2_r==999)]
ggplot(data.frame(cult1_r), aes(cult1_r)) + geom_bar(fill = "#FF6666") +  theme_bw()+
  theme(axis.text.x = element_text(angle = 80, hjust = 1),plot.background = element_blank(),
        panel.grid= element_blank(),axis.line = element_line(color = 'black'))+xlab("")+ylab("")
ggsave(paste("hist_cult1_",yr1[1],"-",yr2[1],".pdf",sep=""))


ggplot(data.frame(cult2_r), aes(cult2_r)) + geom_bar(fill="#FF6666") +  theme_bw()+
  theme(axis.text.x = element_text(angle = 80, hjust = 1),plot.background = element_blank(),
        panel.grid= element_blank(),axis.line = element_line(color = 'black'))+xlab("")+ylab("")
ggsave(paste("hist_cult2_",yr1[1],"-",yr2[1],".pdf",sep=""))





####### non-parametric analysis



cat("- Non-parametric estimation \n")   

# nonlinear effect of prevalence  
est_np_prev <-  gam(lrend ~ s(prev,k=4)+s(lsurf,k=5)+factor(pertcult1)+factor(cult1)+
                      factor(cult2)+factor(annee),family=gaussian,data=datause)
est_np_prev      <- getViz(est_np_prev)
pl1_prev         <- plot( sm(est_np_prev,1),allTerms = FALSE)
pdf(paste("schisto_prev_gam_",paste(unique(annee),collapse="_"),".pdf",sep=""))
pl1_prev + l_fitLine(colour = "blue") +l_ciLine(mul = 1.96, colour = "black", linetype = 2) +
  l_ciLine(mul = 1.96, colour = "black", linetype = 3) +
  l_points(shape = 19, size = 1, alpha = 0.1) + ylim(-1,0.2)+ xlim(0,0.6)+
  labs(x="Schisto Prevalence", y="Polynomial Spline, 6 knots",title=paste("Schisto on log yield/ha, ",paste(unique(annee),collapse="+"),sep=""))  +
  theme_classic()
dev.off()


# nonlinear effect of intensity

form_nl_int    <-as.formula(paste( "lrend~ s(intensity,k=3) +",
                                   paste( paste(X1, collapse="+"),paste(X2, collapse="+"),
                                          "factor(annee)+factor(cult1)+factor(cult2)+
                                          factor(localisation) + factor(relief)+ factor(antierosif)+ factor(mode_acq)+
                                          factor(niv_acq)+ factor(type_labour)+ factor(educ_max)",paste(X4, collapse="+"),sep="+")))


est_np_int <- gam(form_nl_int,family=gaussian,data=datause)
# est_np_int <-  gam(lrend ~ s(intensity)+factor(annee),family=gaussian,data=datause)
est_np_int   <- getViz(est_np_int)
pl1_int      <- plot( sm(est_np_int,1),allTerms = FALSE)
pdf(paste("schisto_intensity_full3_gam_",paste(unique(annee),collapse="_"),".pdf",sep=""))
pl1_int + l_fitLine(colour = "blue") +l_ciLine(mul = 1.96, colour = "black", linetype = 2) +
  l_points(shape = 19, size = 1, alpha = 0.1) + xlim(0,110)+ scale_y_continuous(limits=c(-1.8,0.5))+
  labs(x="Intensity",y="Log Yield"
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


# relationship with plot surface
est_np_int_surf <-  gam(lrend ~ s(intensity,k=5)+s(lsurf,k=5)+factor(cult1)+
                          factor(cult2)+factor(annee),family=gaussian,data=datause)
pdf(paste("schisto+logsurf_",paste(unique(datause$annee),collapse="_"),".pdf",sep=""))
vis.gam(est_np_int_surf,view=c("intensity","lsurf"),theta=40, too.far=1,ticktype="detailed",main = paste("Schistosomiasis and Plot Surface, ",paste(unique(annee),collapse="+"),sep=""),xlab="Intensity",
        ylab="Log Plot Surface",zlab = "Log Yield/Ha")
dev.off() 
est_np_int_surf_2 <-  gam(lrend ~  s(intensity,k=5)+s(lsurf,k=5)+ remune + entraid + 
                            tot_anim_traction + males_elev + fem_elev + totanim + morts + 
                            naissances + bovins + equins + ovins + caprins + dons + vols + 
                            autocons + ventebovins + achatbovins + venteovins + achatovins + 
                            ventecaprins + achatcaprins + gestion + parc_recup + mode_labour + 
                            duree_jachere + npk_kg + uree_kg + phosph_kg + pest_solide_g + 
                            pest_liquide_cl + herbicide_g + herbicide_cl + fungicide_g + 
                            fungicide_cl + rodenticide_g + rodenticide_cl + nmem + agemoyen + 
                            age_chefmen + numfamilies + nchildren + nyoung + nwomen + 
                            culture_pluviale + culture_maraichere + arboriculture + autre_contresaison + 
                            peche + factor(annee) + factor(cult1) + factor(cult2) + factor(localisation) + 
                            factor(relief) + factor(antierosif) + factor(mode_acq) + 
                            factor(niv_acq) + factor(type_labour) + factor(educ_max) + 
                            factor(pertcult1) + factor(pertcult2) + temp_mean + temp_day + 
                            temp_night + precip_wet + temp_dry + dspell_mean + dspell_max + 
                            lstD_yr_mean + lstN_yr_mean + EVI_yr_mean + NDVI_yr_mean,
                          family=gaussian,data=datause)


est_np_int_surf_22   <- getViz(est_np_int_surf_2)
pdf(paste("schisto+logsurf_int_2_",paste(unique(datause$annee),collapse="_"),".pdf",sep=""))
pl_int_surf        <-  plot( sm(est_np_int_surf_22,1))  + labs(x="Schistosomiasis Intensity", 
                                                               y="Log Plot Area", title=" ") + 
  guides(fill=guide_legend(title="Resp.",reverse=T,nrow=4))+ theme_classic()
pl_int_surf
dev.off()

pdf(paste("schisto+logsurf_heat_",paste(unique(datause$annee),collapse="_"),".pdf",sep=""))
vis.gam(est_np_int_surf_2,view=c("intensity","lsurf"),plot.type="contour",color="heat",xlab="Intensity",
        ylab="Log Plot Surface",main = " ")
dev.off()  

# tensor+malaria (yr1)
# dat_mal <- datause[which(datause$annee==unique(annee)[1]),]
est_np_tens_mal  <-  gam(lrend ~ s(prev,malaria,k=6)+s(lsurf,k=5)+
                           + remune + entraid + 
                           tot_anim_traction + males_elev + fem_elev + totanim + morts + 
                           naissances + bovins + equins + ovins + caprins + dons + vols + 
                           autocons + ventebovins + achatbovins + venteovins + achatovins + 
                           ventecaprins + achatcaprins + gestion + parc_recup + mode_labour + 
                           duree_jachere + npk_kg + uree_kg + phosph_kg + pest_solide_g + 
                           pest_liquide_cl + herbicide_g + herbicide_cl + fungicide_g + 
                           fungicide_cl + rodenticide_g + rodenticide_cl + nmem + agemoyen + 
                           age_chefmen + numfamilies + nchildren + nyoung + nwomen + 
                           culture_pluviale + culture_maraichere + arboriculture + autre_contresaison + 
                           peche +factor(cult1) + factor(cult2) + factor(localisation) + 
                           factor(relief) + factor(antierosif) + factor(mode_acq) + 
                           factor(niv_acq) + factor(type_labour) + factor(educ_max) + 
                           factor(pertcult1) + factor(pertcult2) + temp_mean + temp_day + 
                           temp_night + precip_wet + temp_dry + dspell_mean + dspell_max + 
                           lstD_yr_mean + lstN_yr_mean + EVI_yr_mean + NDVI_yr_mean,
                         family=gaussian,data=datause)
est_np_tens_mal   <- getViz(est_np_tens_mal)
pdf(paste("schisto_malaria_gam_",paste(unique(annee),collapse="_"),".pdf",sep=""))
pl2_s         <-  plot( sm( est_np_tens_mal,1)) + labs(x="Schisto Prevalence", 
                                                       y="Median Malaria Prevalence ", 
                                                       title=paste("Interaction Schisto/Malaria, ",paste(unique(annee),collapse="+"),sep=""))   + theme_classic()
pl2_s
dev.off()

pdf(paste("schisto_malaria_3d_2_gam_",paste(unique(annee),collapse="_"),".pdf",sep=""))
vis.gam(est_np_tens_mal,color="heat",view=c("prev","malaria"),theta=150,too.far=2,ticktype="detailed",
        # main = "Interaction schistosomiasis/malaria" ,
        xlab="\n Schistosomiasis",
        ylab="\n Malaria",zlab = "\n Log Yield/Ha",cex.lab=1.6,cex.axis=1.6)
dev.off()







# pdf(pl2_s,paste("schisto+malaria_",unique(dat_mal$annee),".pdf",sep=""))

# using all data

form_all    <- as.formula(paste(y, "~ s(prev_2015map, k=5) +",
                                paste( paste(X1, collapse="+"),paste(X2, collapse="+"),
                                       "factor(annee)+factor(cult1)+factor(cult2)+
                                       factor(localisation) + factor(relief)+ factor(antierosif)+ factor(mode_acq)+
                                       factor(niv_acq)+ factor(type_labour)+ factor(educ_max)",paste(X4, collapse="+"),sep="+")))

# est_np_all   <-  gam(lrend ~ s(prev_2015map,k=6)+log(superficie)+
#                        factor(annee)+factor(pertcult1),family=gaussian,data=data)
est_np_all   <-  gam(form_all,family=gaussian,data=data)
est_np_all   <- getViz(est_np_all)
pl1          <- plot( sm(est_np_all,1),allTerms = FALSE)
pdf("schisto_prev_gam_allyears.pdf")
pl1 + l_fitLine(colour = "blue") +l_ciLine(mul = 1.96, colour = "black", linetype = 2) +
  l_ciLine(mul = 1.96, colour = "black", linetype = 3) +
  l_points(shape = 19, size = 1, alpha = 0.1) + ylim(-1,0.3)+
  labs(x="Prevalence", y="Polynomial Spline, 6 knots",title="2003-2017")  +
  theme_classic()
dev.off()





############ study of prevalence instead of intensity

res_prev           <- matrix(0,5,5)
colnames(res_prev) <- c("Estimate","S.E.","p-val","Mean Effect (%)","Top 5% Effect (%)")

est_fe_prev <- felm(as.formula(paste(y,"~",paste(p,paste(x_log,collapse="+"),sep="+"),"|annee+mid|0|vid",sep="")),data=datause)


res_prev[1,1:3] <- summary(est_fe_prev)$coefficients[1,c(1,2,4)]
res_prev[1,4] <- res_prev[1,1]*mean(datause$prev)
res_prev[1,5] <- res_prev[1,1]*quantile(datause$prev,0.95)

# IV

form_1stage_tfe_prev <-  as.formula(
  paste(paste(p,"~ ",paste(paste(instr,collapse="+"),paste(x_log,collapse="+"),sep="+")),sep=""))

form_2stage_tfe_prev <-  as.formula(
  paste(paste(y,"~ ",paste("stage1_tfe_prev",paste(x_log,collapse="+"),sep="+")),sep=""))


data_snails_pond_biom$stage1_tfe_prev <- fitted(lm(form_1stage_tfe_prev,data=data_snails_pond_biom))

second_tfe_prev             <- lm(form_2stage_tfe_prev,data=data_snails_pond_biom)

set.seed(42)
est_iv_prev  <- coeftest(second_tfe_prev,vcov=vcovCL,cluster = ~ vid)

res_prev[2,1:3] <- est_iv_prev[2,c(1,2,4)]
res_prev[2,4] <- res_prev[2,1]*mean(datause$prev)
res_prev[2,5] <- res_prev[2,1]*quantile(datause$prev,0.95)


K            <- 3
trim         <- c(0,1)
set.seed(42)
# watch out for cv.folds option - sometimes it breaks
Boosting     <- list(bag.fraction = .5, train.fraction = 1.0, interaction.depth=2, n.trees=1000, shrinkage=.01, n.cores=1, cv.folds=1, verbose = FALSE, clas_dist= 'adaboost', reg_dist='gaussian')
Forest       <- list(clas_nodesize=1, reg_nodesize=5, ntree=2000, na.action=na.omit, replace=TRUE)
RLasso       <- list(penalty = list(homoscedastic = FALSE, X.dependent.lambda =FALSE, lambda.start = NULL, c = 1.1), intercept = TRUE)
Nnet         <- list(size=2,  maxit=1000, decay=0.01, MaxNWts=5000,  trace=FALSE)
Trees        <- list(reg_method="anova", clas_method="class")
arguments    <- list(Boosting=Boosting, Forest=Forest, RLasso=RLasso, Nnet=Nnet, Trees=Trees)

methods       <- c("Boosting","RLasso","Nnet")


est_full_DML_prev             <- Full_DML(datause_ML, y, p, x_t, x.flex_t, methods=methods, nfold=K, est="plinear", arguments=arguments, ensemble=ensemble, silent=FALSE, trim)
est_full_DML_mid_prev         <- Full_DML(wdata_mid, y,p, x, x.flex, methods=methods, nfold=K, est="plinear", arguments=arguments, ensemble=ensemble, silent=FALSE, trim)
# est_full_DML_mid_f  <- Full_DML(wdata_mid, y, int, x[1:2], x.flex, methods=methods_f, nfold=K, est="plinear", arguments=arguments, ensemble=ensemble, silent=FALSE, trim)

pvals_DML_prev           <- 2*pt(-abs(est_full_DML_prev [1,c(1:(dim(est_full_DML_prev)[2]-1))]/est_full_DML_prev[2,c(1:(dim(est_full_DML_prev)[2]-1))]), df = nrow(datause)-2)
pvals_DML_mid_prev       <- 2*pt(-abs(est_full_DML_mid_prev[1,c(1:(dim(est_full_DML_mid_prev)[2]-1))]/est_full_DML_mid_prev[2,c(1:(dim(est_full_DML_mid_prev)[2]-1))]), df = nrow(datause)-2 - length(unique(datause$mid)))
# pvals_DML_mid_f      <- 2*pt(-abs(est_full_DML_mid_f[1,c(1:(dim(est_full_DML_mid_f)[2]-1))]/est_full_DML_mid_f[2,c(1:(dim(est_full_DML_mid_f)[2]-1))]), df = nrow(datause)-2 - length(unique(mid)))

best_prev                <- which(est_full_DML_prev[3,1:3]==min(est_full_DML_prev[3,1:3]))
best_mid_prev            <- which(est_full_DML_mid_prev[3,1:3]==min(est_full_DML_mid_prev[3,1:3]))



est_ML_base_int_mid_IV_prev <- DML_Forest_IV(y,d=p,x_t,instr=instr,tune=F,data=data_snails_pond_biom,boot=T,id="annee",cl="vid")



res_prev[3,1]   <- est_ML_base_int_mid_IV_prev$theta_2sls
res_prev[3,2]   <- est_ML_base_int_mid_IV_prev$theta_2sls.se
res_prev[3,3]   <- est_ML_base_int_mid_IV_prev$theta_2sls.p
res_prev[3,4]   <- res_prev[3,1]*mean(datause$prev)
res_prev[3,5]   <- res_prev[3,1]*quantile(datause$prev,0.95)




est_prev_v <- felm(l_rend_avg~ + l_surf_mean + l_dspell_mean + 
                     l_precip_wet + l_temp_mean + l_dspell_max + l_temp_night + 
                     l_temp_day + l_totanim + l_nmem + l_nyoung + l_nchildren + 
                     l_nwomen + l_agemoyen + l_bovins + l_equins + l_ovins + l_caprins + 
                     l_dons + l_vols + l_NDVI_yr_mean + l_EVI_yr_mean + l_lstN_yr_mean + 
                     l_lstD_yr_mean + l_tot_anim_avg + l_npk_avg + l_uree_avg + 
                     l_phosph_avg + l_pest_solide_avg + l_pest_liquide_avg + l_herbicide_cl_avg + 
                     l_fungicide_g_avg + l_fungicide_cl_avg + l_rodenticide_cl_avg + 
                     l_rodenticide_g_avg + I(l_surf_mean^2) + malaria  | annee  | 
                     (prev ~ dry+rainy+mean) | vid,data_snails_pond_biom_v)




res_prev                  <- round(res_prev,4)
res_prev[,4:5]            <- 100*res_prev[,4:5]
res_prev[,4:5]            <- round(res_prev[,4:5],2)


res_inv_prev   <- apply(res_prev[1:3,],2,rev)
cint_l_prev   <- res_inv_prev[,1] - 1.96*res_inv_prev[,2]
cint_h_prev   <- res_inv_prev[,1] + 1.96*res_inv_prev[,2]
tot_eff_prev   <- as.character(res_inv_prev[,4])

res_plot_prev           <- cbind(res_inv_prev,cint_l_prev,cint_h_prev)
ticks_prev              <- rev(c("Linear",
                                 "Linear, 2SLS",
                                 "Flexible,\nML-2SLS+Forest"))
plot_prev               <- broom::tidy(res_plot_prev)



# require("ggrepel")
ggplot(plot_prev,aes(Estimate,factor(ticks_prev,levels=ticks_prev),stat="identity")) + 
  geom_point(size=1.9) +
  geom_errorbarh(aes(xmax = cint_l_prev, xmin = cint_h_prev),height=0.2,linetype=2,alpha=09,colour="black",size=0.3)+
  geom_vline(xintercept=0,linetype="solid",colour="red")+
  geom_label(aes(label=tot_eff_prev,color=p.val),size=3)+
  xlab("Estimate (Prevalence)") + 
  ylab("")+
  xlim(min(cint_l_prev),max(cint_h_prev))  +
  theme(axis.text.y = element_text(lineheight=1,face="bold", color="black",size=, angle=0),
        legend.key.size = unit(0.3, "cm"),panel.border = element_rect(colour = "black", fill=NA, size=5))+
  scale_y_discrete(expand=c(0.1,0))+
  scale_color_gradient(name="p.val",low="blue",high="red",space="Lab",breaks=c(0.15,0.1,0.05))+
  theme_classic2()
ggsave("est_full_prev.pdf",height = 6, units="cm")






