## ----knitr, echo = FALSE-------------------------------------------------
library(knitr)
opts_chunk$set(message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8)

## ----netrc_req, echo = FALSE---------------------------------------------
# This chunk is only useful for BioConductor checks and shouldn't affect any other setup
ISR_login <- Sys.getenv("ISR_login")
ISR_pwd <- Sys.getenv("ISR_pwd")
if(ISR_login != ""  &  ISR_pwd != ""){
  netrc_file <- tempfile("ImmuneSpaceR_tmp_netrc")
  netrc_string <- paste("machine www.immunespace.org login", ISR_login, "password", ISR_pwd)
  write(x = netrc_string, file = netrc_file)
  labkey.netrc.file <- netrc_file
}

## ------------------------------------------------------------------------
library(ImmuneSpaceR)
library(data.table)
library(ggplot2)

## ------------------------------------------------------------------------
con <- CreateConnection("SDY144")
flow <- con$getDataset("fcs_analyzed_result")
hai  <- con$getDataset("hai")
vn   <- con$getDataset("neut_ab_titer")

## ----subset--------------------------------------------------------------
pb <- flow[population_name_reported %in% c("Plasma cells,Freq. of,B lym CD27+",
                                           "Plasmablast,Freq. of,Q3: CD19+, CD20-")]
pb <- pb[, population_cell_number := as.numeric(population_cell_number)]
pb <- pb[study_time_collected == 7 & study_time_collected_unit == "Days"] #13 subjects
pb <- pb[, list(participant_id, population_cell_number, population_name_reported)]

## ----FC------------------------------------------------------------------
# HAI
hai <- hai[,response:=value_reported/value_reported[study_time_collected==0],
                 by="virus,cohort,participant_id"][study_time_collected==30]
hai <- hai[, list(participant_id, virus, response)]
dat_hai <- merge(hai, pb, by = "participant_id", allow.cartesian = TRUE)
# VN
vn <- vn[, response:=value_reported/value_reported[study_time_collected==0],
                 by="virus,cohort,participant_id"][study_time_collected==30]
vn <- vn[, list(participant_id, virus, response)]
dat_vn <- merge(vn, pb, by = "participant_id", allow.cartesian = TRUE)

## ----HAI, dev='png'------------------------------------------------------
ggplot(dat_hai, aes(x = population_cell_number, y = response)) +
  geom_point() + geom_smooth(method = "lm") +
  facet_grid(virus~population_name_reported, scale = "free") +
  xlab("Number of cells") + ylab("HI fold-increase Day 30 vs. baseline") + theme_IS()

## ----VN, dev='png'-------------------------------------------------------
ggplot(dat_vn, aes(x = population_cell_number, y = response)) +
  geom_point() + geom_smooth(method = "lm") +
  facet_grid(virus~population_name_reported, scale = "free") +
  xlab("Number of cells") + ylab("VN fold-increase Day 30 vs. baseline") + theme_IS()

