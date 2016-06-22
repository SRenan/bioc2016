## ----knitr-opts, echo = FALSE--------------------------------------------
library(knitr)
opts_chunk$set(message = FALSE, warning = FALSE)
opts_chunk$set(fig.align = "center", fig.width = 10)

## ----usecase-------------------------------------------------------------
library(GEOquery)
gse <- getGEO("GSE29619")

## ------------------------------------------------------------------------
gse

## ------------------------------------------------------------------------
names(gse)

## ------------------------------------------------------------------------
library(Biobase)
es_LAIV <- gse[[1]]
head(pData(es_LAIV), 3)

## ----combine-es, error = TRUE--------------------------------------------
combine(gse[[1]], gse[[2]])

## ----CreateConnection269-------------------------------------------------
library(ImmuneSpaceR)
sdy269 <- CreateConnection("SDY269")
sdy269

## ----getGEM_TIV2008------------------------------------------------------
TIV2008 <- sdy269$getGEMatrix("TIV_2008")
TIV2008

## ----pdata_TIV2008-------------------------------------------------------
head(pData(TIV2008))

## ----getGEM_multicohorts-------------------------------------------------
es269 <- sdy269$getGEMatrix(c("TIV_2008", "LAIV_2008"), summary = TRUE)
es269

## ----getGEM_cohortnames--------------------------------------------------
es269 <- sdy269$getGEMatrix(cohort = c("TIV Group 2008", "LAIV group 2008"), summary = TRUE)

## ------------------------------------------------------------------------
all <- CreateConnection("") #All studies
es <- all$getGEMatrix(c("TIV_2007", "TIV_2008", "LAIV_2008"), summary = TRUE)
head(pData(es))

## ----listds--------------------------------------------------------------
library(data.table)
sdy269$listDatasets()
hai <- sdy269$getDataset("hai")
hai

## ----getDataset_cross----------------------------------------------------
ahai <- all$getDataset("hai")

## ----makeFilter----------------------------------------------------------
library(Rlabkey)
virus_filter <- makeFilter(c("virus", "CONTAINS", "H1N1"))
hai_f <- sdy269$getDataset("hai", colFilter = virus_filter)
virus_filter2 <- makeFilter(c("virus", "EQUAL", "A/Brisbane/59/2007 (H1N1)"))
hai_f <- sdy269$getDataset("hai", colFilter = virus_filter2)

#multiple filters can be specified
analyte_filter <- makeFilter(c("Analyte", "EQUAL", "IFNg"), c("Study time collected", "IN", "0;7"))
elisa <- sdy269$getDataset("elisa", colFilter = analyte_filter)

## ----cross-assay---------------------------------------------------------
# Elispot
analyte_filter2 <- makeFilter(c("Analyte", "EQUAL", "IgG"), c("Study time collected", "EQUAL", "7"))
elispot <- sdy269$getDataset("elispot", colFilter = analyte_filter2, reload = TRUE)
elispot <- elispot[, elispot_response := spot_number_reported + 1]
elispot <- elispot[, list(participant_id, elispot_response)]

# Flow
fcs <- sdy269$getDataset("fcs_analyzed_result")
fcs <- fcs[, fcs_response := (as.double(population_cell_number) + 1) / as.double(base_parent_population)][study_time_collected == 7]
res <- merge(elispot, fcs, by = "participant_id")

library(ggplot2)
ggplot(res, aes(x = as.double(fcs_response), y = elispot_response, color = cohort)) +
  geom_point() + scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") + 
  xlab("Total plasmablasts (%)") + ylab("Influenza specific cells\n (per 10^6 PBMCs)") +
  theme_IS()

## ----quick-plot----------------------------------------------------------
sdy269$quick_plot("hai", normalize = FALSE)
sdy269$quick_plot("hai", filter = virus_filter2, normalize = FALSE, color = "Age", shape = "Gender")

## ----qp_cross------------------------------------------------------------
virus_filter3 <- makeFilter(c("cohort", "contains", "TIV"), c("study_time_collected", "IN", "0;21;28;30;180"))
all$quick_plot("hai", filter = virus_filter3, normalize = TRUE, color = "Age")

