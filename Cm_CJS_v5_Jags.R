rm(list=ls())
library(jagsUI)
library(tidyverse)
library(lubridate)
library(reshape)
library(bayesplot)
library(RMark)
#library(ggridges)
library(loo)

# library(RMark)
# library(R2ucare)

source("Mancini_functions.R")
save.fig <- F

MCMC.params <- list(n.samples = 50000,
                    n.burnin = 30000,
                    n.thin = 5,
                    n.chains = 5)

# Changed from community to site in Sept 2024 as per Agnes's suggestion
site.names <- c("BKS", "CI_IES", "CL_BMA", "GNO", "MUL", "PAO")
# community.names <- c("BKS", "BMA", "GNO", 
#                      "IES", "LSI", "MUL",
#                      "PAO", "PLM")

#Need to figure out which models were considered best and second best in Mark
best.model.ID <- best.2.model.ID <- vector(mode = "numeric", length = length(site.names))
best.model.name <- best.2.model.name <- vector(mode = "character", length = length(site.names))

k <- 1                                                
for (k in 1:length(site.names)){            
  Cm.results.Mark <- readRDS(file = paste0("RData/CJS_Cm_2023_RMark_", site.names[k], ".rds"))
  
  models.df <- model.table(Cm.results.Mark)
  
  best.model.ID[k] <- as.numeric(row.names(models.df)[1])
  best.model.name[k] <- models.df[1, "model"]
  
  best.2.model.ID[k] <- as.numeric(row.names(models.df)[2])
  best.2.model.name[k] <- models.df[2, "model"]
}

tmp1 <- str_replace_all(best.model.name, "p", "_p")
tmp2 <- str_replace_all(tmp1, "\\(~1\\)", "dot")
tmp3 <- str_replace_all(tmp2, "\\(~time\\)", "t")
tmp4 <- str_replace(tmp3, "\\(~tsm\\)", "TSM")
best.model.names <- str_replace(tmp4, "\\(~effort\\)", "effort")

tmp1.2 <- str_replace_all(best.2.model.name, "p", "_p")
tmp2.2 <- str_replace_all(tmp1.2, "\\(~1\\)", "dot")
tmp3.2 <- str_replace_all(tmp2.2, "\\(~time\\)", "t")
tmp4.2 <- str_replace_all(tmp3.2, "\\(~tsm\\)", "TSM")
tmp5.2 <- str_replace(tmp4.2, "\\(~effort\\)", "effort")
best.2.model.names <- str_replace_all(tmp5.2, "\\(~tsm \\+ effort\\)", "TSM-effort")

# These are the ones used in Cm_CJS_Mark_jags_community Aug 2021 v2.Rmd
# For GNO, Phi[TSM]p[t] resulted in near 1 survival rate estimate, so used the second best
# model, which was Phi[~1]p[t]
models.MARK <- data.frame(site = site.names,
                          ID.1 = best.model.ID,
                          ID.2 = best.2.model.ID,
                          model.1 = best.model.names,
                          model.2 = best.2.model.names)

#dat.1 <- get.data("data/GTC_20190725_Tomo_v2.csv")
k <- 2
# Use the input files from Mark analysis.
for (k in 1:length(site.names)){
  
  Cm.inputs <- readRDS(file = paste0("RData/CJS_Cm_2023_RMark_input_", 
                                     site.names[k], ".rds"))

  Cm.output <- readRDS(file = paste0("RData/CJS_Cm_2023_RMark_", 
                                     site.names[k], ".rds"))
  CH <- Cm.inputs$CH.1 %>% select(-ID) %>% as.matrix()
  
  nInd <- sum(Cm.inputs$CH.R2ucare$effY)
  
  ns <- colSums(Cm.inputs$CH.1 %>% select(-ID))
  
  cap.dates <- paste0(colnames(CH), "-01")
  delta.dates <- signif(as.numeric(as.Date(cap.dates[2:length(cap.dates)]) -
                                     as.Date(cap.dates[1:(length(cap.dates)-1)]))/365, 1)

  # find the first capture date
  get.first <- function(x) min(which(x != 0))
  first.cap <- apply(CH, 1, get.first)
  
  m <- matrix(data = 2, nrow = nrow(CH), ncol = ncol(CH))
  for (k2 in 1:nrow(m)){
    m[k2, first.cap[k2]] <- 1
  }

  ## parameters to monitor 
  parameters <- c("beta", "gamma", "N", "prop.trans", 
                  "r", "K",   
                  "mean.phi",  "phi",
                  "beta.p", "mean.p", "p",    
                  "deviance")
  
  jags.data <- list(y = CH, 
                    f = first.cap, 
                    nind = nInd, 
                    n.occasions = dim(CH)[2],
                    n = ns, 
                    T = ncol(CH),
                    dt = c(0, delta.dates),
                    mean.K = 2000,
                    n.caught = as.vector(colSums(CH)),
                    m = m)
  k3 <- 2
  for (k3 in 1:2){
    if (length(grep("effort", models.MARK[k, paste0("model.", k3)])) > 0){
      dat.1.Cm <- get.data.Cm.2023("data/240925_Base datos_Cm_2001-2023 TE.csv")
      dat.1.Cm.site <- filter(dat.1.Cm, 
                              Macro_site == site.names[k]) 
      
      first.sample. <- first.sample(site.names[k])    
      
      dat.2.Cm.site <- filter(dat.1.Cm.site, 
                              DATE >= first.sample.$first.sample.date)
      
      dat.2.Cm.site %>% select(season, "DATE") %>% 
        group_by(season) %>% #-> tmp3
        summarise(effort = n_distinct(DATE)) -> effort.season
      
      jags.data$x <- as.vector(effort.season$effort)
      
    }
    
    
    jags.input <- list(CH = CH,
                       jags.data = jags.data,
                       parameters.to.save = parameters,
                       run.date = Sys.Date())
    
    out.input.file.name <- paste0("RData/CJS_Cm_2023_jags_input_v5_M", 
                                  filter(models.MARK, 
                                         site == site.names[k])[paste0("ID.", k3)], "_", 
                                  site.names[k], ".rds")
    
    if (!file.exists(out.input.file.name))
      saveRDS(jags.input, 
              file = out.input.file.name)        
  
  # models need to change according to the output from Mark, or do we do 
  # a similar model selection process using Pareto K?

    MCMC.params$model.file <- paste0("models/Model_CJS_", 
                                     filter(models.MARK, 
                                            site == site.names[k])[paste0("model.", k3)],
                                     ".txt")

    out.file.name <- paste0("RData/CJS_Cm_2023_jags_v5_M",
                            filter(models.MARK, 
                                   site == site.names[k])[paste0("ID.", k3)], "_", 
                            site.names[k], ".rds")
    
    if (!file.exists(out.file.name)){
      
      begin.time <- Sys.time()
      jm <- jags(data = jags.data,
                 parameters.to.save= parameters,
                 model.file = MCMC.params$model.file,
                 n.chains = MCMC.params$n.chains,
                 n.burnin = MCMC.params$n.burnin,
                 n.thin = MCMC.params$n.thin,
                 n.iter = MCMC.params$n.samples,
                 DIC = T, 
                 parallel=T)
      
      run.time <- Sys.time() - begin.time
      jags.out <- list(jm = jm,
                       MCMC.params = MCMC.params,
                       run.date = Sys.Date(),
                       run.time = run.time,
                       System = Sys.getenv(),
                       Model = read_delim(MCMC.params$model.file, "\n"))
      
      saveRDS(jm, 
              file = out.file.name)
      
      # loo.out[[k1]] <- compute.LOOIC(loglik = jm$sims.list$loglik, 
      #                                MCMC.params = MCMC.params, 
      #                                data.vector = as.vector(jags.data$y))
      rm(list = c("jm"))
    }
  
  }  
  
}

#saveRDS(loo.out, file = paste0("RData/CJS_Cm_jags_", site.names[k], "_loo.rds"))
        

