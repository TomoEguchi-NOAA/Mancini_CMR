rm(list=ls())
library(jagsUI)
library(tidyverse)
library(lubridate)
library(reshape)
library(bayesplot)
library(ggridges)
library(loo)

# library(RMark)
# library(R2ucare)

source("Mancini_functions.R")
save.fig <- F

MCMC.params <- list(n.samples = 50000,
                    n.burnin = 30000,
                    n.thin = 5,
                    n.chains = parallel::detectCores())

community.names <- c("BKS", "BMA", "GNO", 
                     "IES", "LSI", "MUL",
                     "PAO", "PLM")

# the best models for communities, determined through AICc in MARK:
models.MARK <- data.frame(community = community.names,
                          ID = c(2, 10, 11,
                                 11, 2, 1,
                                 2, 11),
                          model = c("Phidot_pt", "PhiTSM_pdot", "PhiTSM_pt",
                                    "PhiTSM_pt", "Phidot_pt", "Phidot_pdot",
                                    "Phidot_pt", "PhiTSM_pt"))



#dat.1 <- get.data("data/GTC_20190725_Tomo_v2.csv")

dat.1.Cm <- get.data.Cm("data/GTC_Cm Data_updated_2020-04-28_TE_v2.csv")

dat.1.Cm %>% 
  select(ID, CCL) %>% 
  na.omit() %>%
  group_by(ID) %>%
  summarise(min_CCL = min(CCL, na.rm = T)) %>% #-> tmp
  filter(!is.infinite(min_CCL)) -> data.CCL

# some communities don't have enough recaptures to do CMR modeling
c <- 0
k <- 2
k1 <- k2 <- 1


for (k in 1:length(community.names)){
  
  dat.1.Cm.community <- filter(dat.1.Cm, community == community.names[k])
  CJS.data <- dat2CJS(dat.1.Cm.community, save.file = FALSE)
  n.GT2.cap <- length(which(rowSums(CJS.data$data) > 1))
  CJS.data$data %>% rownames_to_column(var = "ID") -> CH.1 #data.CJS
  
  cap.dates <- paste0(colnames(CJS.data$data), "-01")
  delta.dates <- signif(as.numeric(as.Date(cap.dates[2:length(cap.dates)]) -
                                     as.Date(cap.dates[1:(length(cap.dates)-1)]))/365, 1)
  # 
  
  
  # capture history
  CH <- as.matrix(CJS.data$data)
  
  nInd <- nrow(CJS.data$data)
  
  ns <- colSums(CJS.data$data)
  
  # find the first capture date
  get.first <- function(x) min(which(x != 0))
  first.cap <- apply(CH, 1, get.first)
  
  m <- matrix(data = 2, nrow = nrow(CH), ncol = ncol(CH))
  for (k2 in 1:nrow(m)){
    m[k2, first.cap[k2]] <- 1
  }
  
  jags.data <- list(y = CH, 
                    f = first.cap, 
                    nind = nInd, 
                    n.occasions = dim(CH)[2],
                    n = ns, 
                    T = ncol(CJS.data$data),
                    dt = c(0, delta.dates),
                    mean.K = 2000,
                    m = m)
  
  ## parameters to monitor - when this is changed, make sure to change
  ## summary statistics index at the end of this script. 
  parameters <- c("beta", "gamma", "N", "prop.trans", 
                  "r", "K",  "beta.p",
                  "deviance")   # removed "loglik" because not doing LOOIC. 2021-03-26
                # too many parameters are monitored so removing p for now. 2021-03-29 
                # also removing "mu_N", "sigma_N","sigma_logitP",
                # "mu", "mu1", "mu2", "mu.p",

  
  jags.input <- list(raw.data = dat.1.Cm.community,
                     CJS.data = CJS.data,
                     CH.1 = CH.1,
                     jags.data = jags.data,
                     parameters.to.save = parameters,
                     run.date = Sys.Date())
  
  if (!file.exists(paste0("RData/CJS_Cm_jags_input_", community.names[k], ".rds")))
    saveRDS(jags.input, 
            file = paste0("RData/CJS_Cm_jags_input_", community.names[k], ".rds"))        
  
  # models need to change according to the output from Mark, or do we do 
  # a similar model selection process using Pareto K?
  MCMC.params$model.file <- paste0("models/Model_CJS_", filter(models.MARK, 
                                                        community == community.names[k])["model"],
                                   ".txt")
  
  if (!file.exists(paste0("RData/CJS_Cm_jags_M",
                          filter(models.MARK, 
                                 community == community.names[k])["ID"], "_", 
                          community.names[k], ".rds"))){
    
    jm <- jags(data = jags.data,
               parameters.to.save= parameters,
               model.file = MCMC.params$model.file,
               n.chains = MCMC.params$n.chains,
               n.burnin = MCMC.params$n.burnin,
               n.thin = MCMC.params$n.thin,
               n.iter = MCMC.params$n.samples,
               DIC = T, 
               parallel=T)
    
    saveRDS(jm, 
            file = paste0("RData/CJS_Cm_jags_M",
                          filter(models.MARK, 
                                 community == community.names[k])["ID"], "_", 
                          community.names[k], ".rds"))
  }  
  
  # loo.out[[k1]] <- compute.LOOIC(loglik = jm$sims.list$loglik, 
  #                                MCMC.params = MCMC.params, 
  #                                data.vector = as.vector(jags.data$y))
  rm(list = c("jm"))

}

#saveRDS(loo.out, file = paste0("RData/CJS_Cm_jags_", community.names[k], "_loo.rds"))
        

