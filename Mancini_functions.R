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
              season = season,
              SCL = SCL,
              CCL = CCL,
              weight_kg = weight_kg,
              sex = sex)-> dat.1
  
  dat.1[dat.1$species == "In", "species"] <- "Cm"
  
  return(dat.1)
}


dat2CJS <- function(dat.1, save.file = FALSE){
  
  # Create ID by Date and assign 1s
  tmp <-melt(dat.1, 
             id.vars = c("ID", "season"), 
             measure.vars = "detect")
  
  # make a table with ID by season
  dat.01 <- cast(tmp, ID ~ season)
  
  # replace > 1 with ones
  dat.01 <- as.data.frame(dat.01) %>%
    remove_rownames() %>%
    column_to_rownames(var = "ID")
  
  dat.01[(dat.01 > 1)] <- 1
  
  # save file for later
  if (save.file){
    out.name <- "data/Cm_01_CJS_pt.csv"
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



