
first.sample <- function(site.name){
  first.sample.season <- switch(site.name,
                                "BKS" = "2015-01",
                                "CL_BMA" = "2001-08",
                                "GNO" = "2001-08",
                                "CI_IES" = "2006-01",
                                #"LSI" = "2002-08-01",
                                "MUL" = "2008-08",
                                "PAO" = "2001-08")
  
  year <- unlist(str_split(first.sample.season, "-") )[1] %>% as.numeric()
  season.begin <- unlist(str_split(first.sample.season, "-") )[2]
  first.sample.date <- as.Date(ifelse(season.begin == "01",
                                      as.Date(paste0(year-1, "-10-01")),
                                      as.Date(paste0(year, "-05-01"))),
                               format = "%Y-%m-%d") 
  
  return(list(first.sample.season = first.sample.season,
              first.sample.date = first.sample.date))
  
} 

Nhats_comparison <- function(loc, Model.ID.Mark, Model.ID.Jags,
                             N.Phi.p.hats.Mark, N.Phi.hats.Jags, save.fig,
                             fig.height, fig.width){
  
  #Abundance estimates between the two models were similar?
  Nhats <- data.frame(N_mean_Jags = N.Phi.hats.Jags$Nhats$mean,
                      N_median_Jags = N.Phi.hats.Jags$Nhats$median,
                      N_lcl_Jags = N.Phi.hats.Jags$Nhats$lcl,
                      N_ucl_Jags = N.Phi.hats.Jags$Nhats$ucl,
                      N_mean_Mark = c(N.Phi.p.hats.Mark$N.hats$Nhat, NA),
                      N_lcl_Mark = c(N.Phi.p.hats.Mark$N.hats$lcl_Nhat, NA),
                      N_ucl_Mark = c(N.Phi.p.hats.Mark$N.hats$ucl_Nhat, NA))
  
  p.Nhats.comparison <- ggplot(data = Nhats) +
    geom_point(aes(x = N_mean_Jags, y = N_mean_Mark)) +
    geom_errorbar(aes(x = N_mean_Jags, ymin = N_lcl_Mark, ymax = N_ucl_Mark)) +
    geom_errorbarh(aes(y = N_mean_Mark, xmin = N_lcl_Jags, xmax = N_ucl_Jags)) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    ylab(paste0("MLE (M", Model.ID.Mark, ")")) + 
    xlab(paste0("Bayesian (M", Model.ID.Jags, ")")) +
    labs(title = loc)
  
  if (save.fig)
    ggsave(plot = p.Nhats.comparison, 
           filename = paste0("figures/Cm_2023_Nhats_comp_", loc, 
                             "_M", Model.ID.Mark, "_M", Model.ID.Jags, ".png"),
           height = fig.height, width = fig.width,
           device = "png", dpi = 600)
}

plot.Nhats.Jags <- function(loc, 
                            Model.ID,
                            N.Phi.p.hats.Jags, 
                            save.fig, 
                            fig.height, fig.width, ext = "jags_v5.png"){
  p.Nhats <- ggplot(data = N.Phi.p.hats.Jags$Nhats) +
    geom_point(aes(x = season, y = (mean))) +
    geom_errorbar(aes(x = season, ymin = (lcl), ymax = (ucl))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ylab("Abundance (95% CI)") +
    labs(title = paste0(loc, " (Bayesian, Model ", Model.ID, ")"))
  
  if (save.fig)
    ggsave(plot = p.Nhats, 
           filename = paste0("figures/Cm_2023_Nhats_", loc, "_M", Model.ID, "_", ext),
           height = fig.height, width = fig.width,
           device = "png", dpi = 600)
  
  #return(p.Nhats)
}

plot.Nhats.Mark <- function(loc, 
                            Model.ID, 
                            N.Phi.p.hats.Mark, 
                            save.fig, 
                            fig.height, fig.width){
  
  p.Nhats <- ggplot(data = N.Phi.p.hats.Mark$N.hats) +
    geom_point(aes(x = season, y = (Nhat))) +
    geom_errorbar(aes(x = season, ymin = (lcl_Nhat), ymax = (ucl_Nhat))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ylab("Abundance (95% CI)") +
    labs(title = paste0(loc, " (MLE, Model ", Model.ID, ")"))
  
  if (save.fig)
    ggsave(plot = p.Nhats, 
           filename = paste0("figures/Cm_2023_Nhats_", loc, "_M", Model.ID, ".png"),
           height = fig.height, width = fig.width,
           device = "png", dpi = 600)
  
  #return(p.Nhats)
}

capture.recapture.stats <- function(loc, 
                                    dat.1.Cm.community, 
                                    #Cm.inputs, 
                                    save.fig = FALSE, 
                                    fig.height = 4, 
                                    fig.width = 6){
  
  CJS.data <- dat2CJS(dat.1.Cm.community, save.file = FALSE)
  n.occ <- ncol(CJS.data$data)
  n.caught <- colSums(CJS.data$data)
  
  p.captures <- ggplot(dat.1.Cm.community) + 
    geom_point(aes(x = DATE,  y = ID, color = CCL)) +
    geom_path(aes(x = DATE, y = ID)) +
    theme(axis.text.y = element_blank()) +
    labs(title = loc)
  
  if (save.fig)
    ggsave(filename = paste0("figures/Cm_capture_history_", loc, ".png"),
           plot = p.captures, device = "png", dpi = 600,
           height = fig.height, width = fig.width)
  
  n.caps <- rowSums(CJS.data$data) %>% 
    data.frame() %>% 
    rownames_to_column(var = "ID") 
  
  colnames(n.caps) <- c("ID", "n")   # there has to be a better way but... 
  
  p.n.captures <- ggplot(data = n.caps) + 
    geom_histogram(aes(x = n), binwidth = 1) +
    labs(title = loc)
  
  if (save.fig)
    ggsave(filename = paste0("figures/Cm_capture_histogram_", loc, ".png"),
           plot = p.n.captures, device = "png", dpi = 600,
           height = fig.height, width = fig.width)
  
  tmp.counts <- t(apply(CJS.data$data, MARGIN = 1, FUN = n.captures))
  tmp.counts.df <- data.frame(count.captures(tmp.counts)) 
  rownames(tmp.counts.df) <- 1:max(tmp.counts)
  colnames(tmp.counts.df) <- colnames(CJS.data$data)
  tmp.counts.df <- rownames_to_column(tmp.counts.df, var = "n.recaptures")
  tmp.counts.df %>% 
    pivot_longer(!n.recaptures, 
                 names_to = "Season", 
                 values_to = "Freq") -> counts.freq

  p.n.recap <- ggplot(counts.freq) + 
    geom_col(aes(x = Season, 
                 y = Freq, 
                 fill = n.recaptures)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "top") +
    ylab("# individuals") + labs(fill = "# recaptures", title = loc)
  
  
  if (save.fig)
    ggsave(filename = paste0("figures/Cm_recaptures_", loc, ".png"),
           plot = p.n.recap, device = "png", dpi = 600,
           height = fig.height, width = fig.width)
  
  loc.stats <- list(CJS.data = CJS.data,
                    n.occ = n.occ,
                    n.caught = n.caught,
                    n.caps = n.caps,
                    plots = list(captures = p.captures,
                                 n.captures = p.n.captures,
                                 n.recap = p.n.recap))
    
  return(loc.stats)
}

# Counts the number of individuals that are caught x times at each capture occasion. Use this with output from n.captures.
count.captures <- function(x){
  max.n <- max(x)
  out.table <- matrix(nrow = max.n, ncol = ncol(x))
  for (k in 1:ncol(x)){
    for (k1 in 1:max.n){
      out.table[k1, k] <- sum(x[,k] == k1)
    }
  }
  return(out.table)
}


# this function enumerate the capture sequence for each individual. Rather than having 0s and 1s, they are turned into 0s and 1s, 2s, 3s, etc. Cumulative sum but keeping zeros as zeros. Use it with apply to a capture history dataframe. When used with apply, the output is individuals in columns and capture occasions in rows. So, transform the output to make it comparable to the capture history matrix.
#t(apply(X, MARGIN = 1, FUN = n.captures)), where X is the capture history matrix/dataframe
n.captures <- function(x){
  #x.1 <- vector(mode = "numeric", length = ncol(x))
  x.1 <- x[1]
  for (k in 2:length(x)){
    if (x[k] > 0) { 
      x[k] <- x.1 + x[k]
      x.1 <- x[k]
    }
  }
  
  return(x)
}

stats.Bayesian_v5 <- function(model.ID, Phi.spec, p.spec, loc){
  #models.MARK %>% filter(community == loc) %>% select(ID) -> model.ID
  
  Cm.results <- readRDS(file = paste0("RData/CJS_Cm_2023_jags_v5_M", 
                                      model.ID, "_", loc, ".rds"))
  
  Cm.inputs <- readRDS(file = paste0("RData/CJS_Cm_2023_jags_input_v5_M", 
                                     model.ID,"_", loc, ".rds"))
  
  Phi.par <- switch(Phi.spec,
                    "~1" = "mean.phi",
                    "~time" = "phi[",
                    "~tsm" = "phi[")
  
  p.par <- switch(p.spec,
                  "~1" = "mean.p",
                  "~time" = "p[",
                  "~effort" = "p[")

  real.estimates <- Cm.results$summary %>% 
    as.data.frame() %>% 
    rownames_to_column("parameter")
  
  params <- list(phi = Phi.par,
                 p = p.par,
                 p.trans = "prop.trans")
  
  N.Phi.hats <- extract.Nhats.jags_v5(Cm.inputs, real.estimates, params)
  
  N.Phi.hats$Model.ID <- model.ID
  return(N.Phi.hats)  
}

stats.Bayesian <- function(models.MARK, loc){
  models.MARK %>% filter(community == loc) %>% select(ID) -> model.ID
  
  Cm.results <- readRDS(file = paste0("RData/CJS_Cm_jags_M", model.ID, "_", loc, ".rds"))
  
  Cm.inputs <- readRDS(file = paste0("RData/CJS_Cm_jags_input_", loc, ".rds"))
  
  real.estimates <- Cm.results$summary %>% as.data.frame() %>% rownames_to_column("parameter")
  
  N.Phi.hats <- extract.Nhats.jags(Cm.inputs, real.estimates)
  return(N.Phi.hats)  
}

extract.Nhats <- function(Cm.inputs, Cm.results, real.estimates){
  data.0 <- Cm.inputs$CJS.data$data
  
  p.hats <- real.estimates[grep("p", real.estimates$parameter),]
  phi.hats <- real.estimates[grep("Phi", real.estimates$parameter),]
  
  p.hats[which(p.hats[,"estimate"] < 0.001), c("estimate", "se", "lcl", "ucl")] <- NA
  if (nrow(p.hats) == 1){
    p.hats %>% 
      slice(rep(1, each = (length(colnames(data.0))-1))) %>%
      mutate(season = colnames(data.0[1:(ncol(data.0)-1)]),
             n.caught = colSums(data.0[,1:(ncol(data.0)-1)])) -> p.hats
    
    
    
  } else if (nrow(p.hats) != (ncol(data.0)-1)){  # When effort is affecting p
    effort.p.hats <- data.frame(effort = unique(Cm.inputs$effort$effort),
                               p.hat = real.estimates[grep("p", real.estimates$parameter),]) %>%
      select(-c(p.hat.fixed, p.hat.note))
    
    Cm.inputs$effort %>% 
      left_join(effort.p.hats, by = "effort") %>%
      left_join(data.frame(season = colnames(Cm.inputs$CJS.data$data),
                           n.caught = colSums(Cm.inputs$CJS.data$data)),
                by = "season") %>%
      dplyr::rename(parameter = p.hat.parameter,
                    estimate = p.hat.estimate,
                    se = p.hat.se,
                    lcl = p.hat.lcl,
                    ucl = p.hat.ucl) -> p.hats
    
  } else if (nrow(p.hats) == (ncol(data.0)-1)){
    p.hats <- select(p.hats, -c(fixed, note)) %>%
      mutate(season = colnames(data.0[1:(ncol(data.0)-1)]),
             n.caught = colSums(data.0[,1:(ncol(data.0)-1)]))
  }
  
  #phats$season <- colnames(data.0)[1:(ncol(data.0)-1)]
  
  #n.caught <- colSums(data.0)
  
  model.averaged.Phi <- model.average(Cm.results, parameter = "Phi")
  
  p.hats %>% 
    mutate(Nhat = n.caught/p.hats$estimate,
           SE_Nhat = (n.caught/estimate) * se/estimate,
           lcl_Nhat = n.caught/ucl,
           ucl_Nhat = n.caught/lcl) -> N.hats.df
  
  # Nhats.df <- data.frame(season = colnames(data.0)[1:(ncol(data.0)-1)],
  #                        Nhat = (n.caught[1:(length(n.caught) - 1)]/phats$estimate) ) %>%
  #   mutate(SE_Nhat = (n.caught[1:(length(n.caught) - 1)]/phats$estimate) * phats$se/phats$estimate,
  #          #lcl  = (n.caught[2:length(n.caught)]/phats$lcl) * p.residents,
  #          #ucl = (n.caught[2:length(n.caught)]/phats$ucl) * p.residents,
  #          lcl = (n.caught[1:(length(n.caught)-1)]/phats$estimate)  - 1.96 * SE_Nhat,
  #          ucl = (n.caught[1:(length(n.caught)-1)]/phats$estimate)  + 1.96 * SE_Nhat,
  #          lcl2 = ifelse(lcl < 0, 0, lcl))
  
  out.list <- list(p.hats = p.hats,
                   N.hats = N.hats.df,
                   Phi.hats = phi.hats,
                   Phi.hat_avg = model.averaged.Phi)
  return(out.list)
}

extract.Nhats.jags_v5 <- function(Cm.inputs, real.estimates, params){
  #data.0 <- Cm.inputs$CJS.data$data
  data.0 <- Cm.inputs$CH
  Nhats <- real.estimates[grep("N[", real.estimates$parameter, fixed = T),] %>%
    transmute(parameter = parameter,
              mean = mean,
              sd = sd,
              lcl = `2.5%`,
              median = `50%`,
              ucl = `97.5%`,
              Rhat = Rhat,
              season = colnames(data.0))
  
  Phihats <- real.estimates[grep(params$phi, real.estimates$parameter, fixed = T),] %>%
    transmute(parameter = parameter,
              mean = mean,
              sd = sd,
              lcl = `2.5%`,
              median = `50%`,
              ucl = `97.5%`,
              Rhat = Rhat)
  
  phats <- real.estimates[grep(params$p, real.estimates$parameter, fixed = T),] %>%
    transmute(parameter = parameter,
              mean = mean,
              sd = sd,
              lcl = `2.5%`,
              median = `50%`,
              ucl = `97.5%`,
              Rhat = Rhat)

  Gammahats <- real.estimates[grep("gamma", real.estimates$parameter, fixed = T),] %>%
    transmute(parameter = parameter,
              mean = mean,
              sd = sd,
              lcl = `2.5%`,
              median = `50%`,
              ucl = `97.5%`,
              Rhat = Rhat)
  
  p.trans <- real.estimates[grep("prop.trans", real.estimates$parameter, fixed = T),] %>%
    transmute(parameter = parameter,
              mean = mean,
              sd = sd,
              lcl = `2.5%`,
              median = `50%`,
              ucl = `97.5%`,
              Rhat = Rhat)
  
  out.list <- list(Nhats = Nhats,
                   Phihats = Phihats,
                   phats = phats,
                   Gammahats = Gammahats,
                   p.trans = p.trans)
  return(out.list)
}

extract.Nhats.jags <- function(Cm.inputs, real.estimates){
  data.0 <- Cm.inputs$CJS.data$data
  #data.0 <- Cm.inputs$CH
  Nhats <- real.estimates[grep("N[", real.estimates$parameter, fixed = T),] %>%
    transmute(parameter = parameter,
              mean = mean,
              sd = sd,
              lcl = `2.5%`,
              median = `50%`,
              ucl = `97.5%`,
              Rhat = Rhat,
              season = colnames(data.0))
  
  Phihats <- real.estimates[grep("phi", real.estimates$parameter, fixed = T),] %>%
    transmute(parameter = parameter,
              mean = mean,
              sd = sd,
              lcl = `2.5%`,
              median = `50%`,
              ucl = `97.5%`,
              Rhat = Rhat)
  
  Gammahats <- real.estimates[grep("gamma", real.estimates$parameter, fixed = T),] %>%
    transmute(parameter = parameter,
              mean = mean,
              sd = sd,
              lcl = `2.5%`,
              median = `50%`,
              ucl = `97.5%`,
              Rhat = Rhat)
  
  out.list <- list(Nhats = Nhats,
                   Phihats = Phihats,
                   Gammahats = Gammahats)
  return(out.list)
}

do_analysis <- function(dp, ddl)
{
  # create formulas for Phi
  # tsm is time-since-marking; check for transient effects
  Phi.dot <-  list(formula = ~ 1)  
  #Phi.weight <- list(formula= ~ min_weight)   # many missing data 
  #Phi.t <- list(formula = ~ time)             # we never have this model worked for turtles... 
  #Phi.season <- list(formula = ~ sum_win)      # this also is unlikely... 
  #Phi.transience <- list(formula = ~ Transient)
  Phi.tsm <- list(formula = ~ tsm)
  
  #create formulas for p
  p.dot <- list(formula = ~ 1)
  p.t <- list(formula = ~ time)
  p.tsm <- list(formula = ~ tsm)
  #p.transience <- list(formula = ~ Transient)
  #p.tsm.transience <- list(formula = ~ tsm + Transient)
  #p.t.transience <- list(formula = ~ time + Transient)
  p.effort <- list(formula = ~ effort)
  p.season <- list(formula = ~ sum_win)
  p.tsm.season <- list(formula = ~ tsm + sum_win)
  p.tsm.effort <- list(formula = ~ tsm + effort)
  
  # create all combinations 
  cml <- create.model.list("CJS")
  
  # run all all models and return as a list with class marklist
  results <- mark.wrapper(model.list = cml,
                          data=dp,
                          ddl=ddl,
                          output=FALSE,
                          silent=TRUE)
  return(results)
}


Cm.vonBert.jags.data <- function(dat.1){
  # remove rows with is.na(CCL) is true:
  dat.1 %>% filter(!is.na(CCL)) -> dat.2
  
  n.cap.ID <- table(dat.2$ID)
  recap.ID <- data.frame(n.cap.ID[n.cap.ID > 2])
  colnames(recap.ID) <- c("ID", "Freq")
  
  recap.ID %>% left_join(dat.2, by = "ID") -> recap.data
  
  # Make length and capture date matrices
  unique.ID <- recap.ID$ID
  size.mat <- date.mat <- matrix(nrow = length(unique.ID),
                                 ncol = max(recap.data$Freq))
  
  date.1 <- structure(numeric(length(unique.ID)), class = "Date")
  n.vec <- vector(mode = "numeric", length = length(unique.ID))
  
  k <- 1
  for (k in 1:length(unique.ID)){
    tmp.ID <- filter(recap.data, ID == as.character(unique.ID[k]))
    size.mat[k, 1:nrow(tmp.ID)] <- tmp.ID$CCL
    date.mat[k, 1:nrow(tmp.ID)] <- tmp.ID$DATE - min(tmp.ID$DATE)
    date.1[k] <- min(tmp.ID$DATE)
    n.vec[k] <- nrow(tmp.ID)
  }
  
  date.mat <- date.mat[, 2:ncol(date.mat)]/365
  
  jags.data <- list(nIndiv = length(unique.ID),
                    n = n.vec,
                    L = size.mat,
                    t = date.mat)
  
  return(list(jags.data = jags.data,
              ID = unique.ID))
}



vonBert.jags.data <- function(dat.1, sp.code){
  # remove rows with is.na(CCL) is true:
  dat.1 %>% filter(species == sp.code) %>%
    filter(!is.na(CCL)) -> dat.2
  
  n.cap.ID <- table(dat.2$ID)
  recap.ID <- data.frame(n.cap.ID[n.cap.ID > 2])
  colnames(recap.ID) <- c("ID", "Freq")
  
  recap.ID %>% left_join(dat.2, by = "ID") -> recap.data
  
  # Make length and capture date matrices
  unique.ID <- recap.ID$ID
  size.mat <- date.mat <- matrix(nrow = length(unique.ID),
                                 ncol = max(recap.data$Freq))
  
  date.1 <- structure(numeric(length(unique.ID)), class = "Date")
  n.vec <- vector(mode = "numeric", length = length(unique.ID))
  
  k <- 1
  for (k in 1:length(unique.ID)){
    tmp.ID <- filter(recap.data, ID == as.character(unique.ID[k]))
    size.mat[k, 1:nrow(tmp.ID)] <- tmp.ID$CCL
    date.mat[k, 1:nrow(tmp.ID)] <- tmp.ID$DATE - min(tmp.ID$DATE)
    date.1[k] <- min(tmp.ID$DATE)
    n.vec[k] <- nrow(tmp.ID)
  }
  
  date.mat <- date.mat[, 2:ncol(date.mat)]/365
  
  jags.data <- list(nIndiv = length(unique.ID),
                    n = n.vec,
                    L = size.mat,
                    t = date.mat)
  
  return(list(jags.data = jags.data,
              ID = unique.ID))
}

get.data.Cm <- function(filename){
  col.def <- cols(Event_ID = col_character(),
                  Turtle_no = col_integer(),
                  Count = col_integer(),
                  Season = col_character(),
                  Year = col_integer(),
                  Month = col_integer(),
                  Day = col_integer(),
                  Sp_code = col_character(),
                  Turtle_ID = col_character(),
                  Recapture = col_character(),
                  Community_code = col_character(),
                  Start_date = col_date(format = "%m/%d/%Y"),
                  Start_time = col_time(format = "%H:%M"),
                  End_time = col_time(format = "%H:%M"),
                  Tot_hours = col_double(),
                  Tot_hours_estimated = col_double(),
                  Monitoring_type = col_character(),
                  Monitoring_technique = col_character(),
                  Region = col_character(),
                  Site_type_general = col_character(),
                  Type_site_specific = col_character(),
                  Lat = col_double(),
                  Long = col_double(),
                  Capture_date = col_date(format = "%m/%d/%Y"),
                  Species = col_character(),
                  SCL = col_double(),
                  SCW = col_double(),
                  CCL = col_double(),
                  CCW = col_double(),
                  BD = col_double(),
                  PL = col_double(),
                  TTL = col_double(),
                  Weight = col_double(),
                  Sex = col_character(),
                  Right_tag_new = col_character(),
                  Left_tag_new = col_character(),
                  Right_tag_old = col_character(),
                  Left_tag_old = col_character())
  
  dat.1 <- read_csv(file = filename, col_types = col.def)
  
  dat.1 %>% mutate(ID = as.factor(Turtle_ID),
                   CDATE = Capture_date) %>% 
    transmute(ID = ID,
              detect = 1,
              DATE = CDATE,
              season = as.factor(Season),
              SCL = SCL,
              CCL = CCL,
              weight_kg = Weight,
              sex = Sex,
              species = as.factor(Sp_code),
              community = Community_code)-> dat.1
  
  return(dat.1)
}

get.data.Cm.2023 <- function(filename){
  col.def <- cols(Season = col_character(),
                  Year = col_integer(),
                  Month = col_integer(),
                  Day = col_integer(),
                  Region = col_character(),
                  Macro_site = col_character(),
                  Site_code = col_character(),
                  Date = col_date(format = "%m/%d/%Y"),
                  TurtleID = col_character(),
                  Site_name = col_character(),
                  Recapture = col_character(),
                  Rtag_new = col_character(),
                  Ltag_new = col_character(),
                  Rtag_old = col_character(),
                  Ltag_old = col_character(),
                  SCL = col_double(),
                  CCL = col_double(),
                  Weight = col_double(),
                  Data_ownership = col_character())
  
  dat.1 <- read_csv(file = filename, col_types = col.def)
  
  dat.1 %>% mutate(ID = as.factor(TurtleID)) %>% 
    transmute(ID = ID,
              detect = 1,
              DATE = Date,
              SCL = SCL,
              CCL = CCL,
              weight_kg = Weight,
              season = Season,
              community = as.factor(Site_code),
              Macro_site = as.factor(Macro_site))-> dat.1
  
  return(dat.1)
}

get.data <- function(filename){
  col.def <- cols(monitoring_event = col_character(),
                  value = col_integer(),
                  season = col_character(),
                  year = col_integer(),
                  month = col_integer(),
                  day = col_integer(),
                  species = col_character(),
                  turtle_code = col_character(),
                  recapture = col_character(),
                  community = col_character(),
                  start_date = col_date(format = "%d/%m/%Y"),
                  tot_hours = col_double(),
                  type_monitoring = col_character(),
                  methodology = col_character(),
                  longitude_net = col_number(),
                  type_site_gen = col_character(),
                  type_site_spec = col_character(),
                  site_name = col_character(),
                  latitude = col_double(),
                  longitude = col_double(),
                  SCL_CCL = col_double(),
                  SCL = col_double(),
                  SCW = col_double(),
                  CCL = col_double(),
                  CCW = col_double(),
                  BD = col_double(),
                  PL = col_double(),
                  TTL = col_double(),
                  weight_kg = col_double(),
                  sex = col_character(),
                  tipe_monitoring = col_character(),
                  BCI = col_double(),
                  Condition = col_character(),
                  Recapture_new = col_character())
  
  dat.1 <- read_csv(file = filename, col_types = col.def)
  
  dat.1 %>% mutate(ID = as.factor(turtle_code),
                   CDATE = as.Date(paste0(year, "-", month, "-", day))) %>% 
    transmute(ID = ID,
              detect = 1,
              species = as.factor(species),
              DATE = CDATE,
              season = as.factor(season),
              SCL = SCL,
              CCL = CCL,
              weight_kg = weight_kg,
              sex = sex)-> dat.1
  
  # "indefindio" species are greens per Agnese 
  # email on 2019-09-04
  dat.1[dat.1$species == "In", "species"] <- "Cm"
  
  return(dat.1)
}

get.data.Ei <- function(filename){
  col.def <- cols(Event_ID = col_character(),
                  Season = col_character(),
                  Species = col_character(),
                  Turtle_ID = col_character(),
                  Recapture = col_character(),
                  Location = col_character(),
                  Start_time = col_time(format = "%H:%m"),
                  End_time = col_time(format = "%H:%m"),
                  Total_hrs = col_double(),
                  Type = col_character(),
                  Metodologia = col_character(),
                  Longitud_lanceos = col_character(),
                  Site_type_general = col_character(),
                  Site_type_specific = col_character(),
                  site_name = col_character(),
                  Lat_Corr = col_double(),
                  Long_Corr = col_double(),
                  Capture_date = col_date(format = "%m/%d/%Y"),
                  Capture_time = col_time(format = "%H:%m"),
                  Turlte_name = col_character(),
                  SCL = col_double(),
                  SCW = col_double(),
                  CCL = col_double(),
                  CCW = col_double(),
                  BD = col_double(),
                  PL = col_double(),
                  TTL = col_double(),
                  Weight = col_double(),
                  Sex = col_character(),
                  tag_r_new = col_character(),
                  tag_l_new = col_character(),
                  tag_r_old = col_character(),
                  tag_l_old = col_character(),
                  Pit_tag = col_character(),
                  pit_new = col_character(),
                  Pit_old = col_character())
  
  dat.1 <- read_csv(file = filename, col_types = col.def)
  
  dat.1 %>% mutate(ID = as.factor(Turtle_ID),
                   CDATE = Capture_date) %>% 
    transmute(ID = ID,
              season = Season,
              detect = 1,
              DATE = CDATE,
              SCL = SCL,
              CCL = CCL,
              weight_kg = Weight,
              sex = Sex)-> dat.1
  
  return(dat.1)
}


dat2CJS <- function(dat.1, save.file = FALSE, filename = "not.saved"){
  
  # Create ID by Date and assign 1s
  tmp <-melt(dat.1, 
             id.vars = c("ID", "season"), 
             measure.vars = "detect")
  
  # make a table with ID by season
  dat.01 <- cast(tmp, 
                 formula = ID ~ season,
                 fun.aggregate = length)
  
  # replace > 1 with ones
  dat.01 <- as.data.frame(dat.01) %>%
    remove_rownames() %>%
    column_to_rownames(var = "ID")
  
  dat.01[(dat.01 > 1)] <- 1
  
  # save file for later
  if (save.file){
    out.name <- filename
    write.csv(dat.01, 
              file = out.name, 
              row.names = T,
              quote = F)
    
  } else {
    out.name <- "not.saved"
  }
  
  out <- list(filename = out.name,
              data = dat.01)
  
  return(out)
  
}


known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

dat2dat01.year <- function(dat.1, year, save.file = FALSE){
  dat.year <- filter(dat.1, YEAR == year)
  
  # Create ID by Date and assign 1s
  tmp <-melt(dat.year, 
             id.var = c("ID", "DATE"), 
             measure.var = "detect")
  
  # make a table with ID by Date
  dat.01.year <- cast(tmp, ID ~ DATE)
  
  # replace NAs with zeros
  dat.01.year[is.na(dat.01.year)] <- 0
  dat.01.year <- as.data.frame(dat.01.year)
  
  # save file for later
  if (save.file){
    out.name = paste0("data/Cm_01_", year, ".csv")
    write.csv(dat.01.year, 
              file = out.name, 
              row.names = F,
              quote = F)
    
  } else {
    out.name <- "not.saved"
  }
  
  out <- list(filename = out.name,
              data = dat.01.year)
  
  return(out)
  
}


dat2dat01 <- function(dat.1, save.file = FALSE){
  
  # Create ID by Date and assign 1s
  tmp <-melt(dat.1, 
             id.var = c("ID", "season"), 
             measure.var = "detect")
  
  # make a table with ID by year
  dat.01 <- cast(tmp, ID ~ season)
  
  # replace > 1 with ones
  dat.01 <- as.data.frame(dat.01) %>%
    remove_rownames() %>%
    column_to_rownames(var = "ID")
  
  dat.01[(dat.01 > 1)] <- 1
  
  # save file for later
  if (save.file){
    out.name <- "data/Cm_01_all.csv"
    write.csv(dat.01, 
              file = out.name, 
              row.names = T,
              quote = F)
    
  } else {
    out.name <- "not.saved"
  }
  
  out <- list(filename = out.name,
              data = dat.01)
  
  return(out)
  
}


dat2CJS_covCCL <- function(dat.1){
  
  # Create ID by Date and assign 1s
  tmp <-melt(dat.1, 
             id.vars = c("ID", "season"), 
             measure.vars = "CCL")
  
  # make a table with ID by season
  dat.CCL <- reshape2::dcast(tmp, 
                             formula = ID ~ season,
                             value.var = "value",
                             fun.aggregate = mean)
  
  out <- dat.CCL
  return(out)
  
}


compute.LOOIC <- function(loglik, data.vector, MCMC.params){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
  loglik.vec <- as.vector(loglik)
  loglik.mat <- matrix(loglik.vec[!is.na(data.vector)], 
                        nrow = MCMC.params$n.chains * n.per.chain)
  
  Reff <- relative_eff(exp(loglik.mat),
                       chain_id = rep(1:MCMC.params$n.chains,
                                      each = n.per.chain),
                       cores = 4)
  
  loo.out <- rstanarm::loo(loglik.mat, 
                           r_eff = Reff, 
                           cores = 4, k_threshold = 0.7)
  
  out.list <- list(Reff = Reff,
                   loo.out = loo.out)
  
  return(out.list)  
}

# extract looic and pareto k statistics. the first input (loo.out) should come from compute.LOOIC above. 
pareto.looic.fcn <- function(loo.out, models){
  pareto.k <- lapply(loo.out, 
                     FUN = function(x) x$loo.out)
  
  # find maximum pareto k values
  max.pareto.k <- unlist(lapply(pareto.k,
                                FUN = function(x) max(x$diagnostics$pareto_k)))
  
  # find the models that have max(pareto k) < 0.7
  good.models <- models[which(max.pareto.k < 0.7)]
  good.models.pareto.k <- pareto.k[which(max.pareto.k < 0.7)]
  
  looic.esimates <- lapply(lapply(loo.out[which(max.pareto.k < 0.7)], 
                                  FUN = function(x) x$loo.out),
                           FUN = function(x) x$estimates)
  
  looic <- unlist(lapply(looic.esimates, 
                         FUN = function(x) x["looic", "Estimate"]))
  
  loo.out.list <- lapply(loo.out[which(max.pareto.k < 0.7)], 
                         FUN = function(x) x$loo.out)
  
  # calculate model weights
  model.weights <- loo_model_weights(loo.out.list)
  out.list <- list(model.weights = model.weights,
                   looic = looic,
                   good.models = good.models,
                   good.models.pareto.k = good.models.pareto.k)
}



# Extracting posterior samples of deviance or any other variable from jags output:
extract.samples <- function(varname, zm){
  dev <- unlist(lapply(zm, FUN = function(x) x[, varname]))
  return(dev)
}

