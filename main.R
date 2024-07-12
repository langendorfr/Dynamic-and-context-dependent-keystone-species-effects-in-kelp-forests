#### Author: Ryan E. Langendorf
#### License: GPLv3
#### Reproduces main text results in 'Dynamic and context-dependent keystone species effects in kelp forests'

## Primitives
# rm(list = ls()) ## Use to clear all variables so the code runs in a clean environment
# set.seed(42) ## Use to replicate results in the paper

## Packages
library("doParallel")
library("dplyr")
library("foreach")
library("GenSA")
library("ggplot2"); theme_set(theme_bw())
library("ggraph")
library("lme4")
library("magrittr")
library("Matrix")
library("MuMIn")
library("nlme")
library("optimx")
library("piecewiseSEM")
library("readr")
library("reshape2")
library("reticulate") ## Requires a python installation
library("stringr")
library("tibble")
library("tidygraph")

## Primary Directory
# pd <- "/home/ryan/SanNic" ## NOTE: Set to directory with all files and scripts
setwd(pd) ## NOTE: This should throw an error until the previous line is uncommented and set to the directory with all the files and scripts

## Custom functions
source("net_mat.R")

## UUID (1 second resolution) for each model run so we do not accidentally save over previous model fit output
tu_id <- gsub(x = Sys.time(), pattern = "-|:|[.]| ", replacement = "")

## Runtime Parameters
verbose <- TRUE
transform <- "Untransformed" #"Kenner"
normalization <- "unity" #"percent"
weighting_scheme <- "AllCausers"
model_kind <- "piecewiseSEM"
model_type <- "FixedEffects" #"DredgeMixed" #"MixedEffects"
intra_inter_priority <- "inter" #"intra"
model_selection <- "AICc" #"BIC"
location <- "SN" #"BC" ## NOTE: In the paper SN = SNI and BC = WVI
model_penalty_upscale <- 1e3 ## Care more about adhering to the bounds
model_penalty_enforcement <- "ForcedBounds" #"UnforcedBounds"
optim_method <- "optimSAoptim" #"gensaoptim" #"gensa" #"optimx" #"gensa" #"dfd"
optim_iters <- 3
optim_bound <- 10
noise <- 0 #1e-10
theta <- "pcapsem"
date <- "2024-5-10"
removedabsentinteractors <- "RemovedAbsentInteractors"
state <- "pca" #umap #smap
cutoff <- 999 ## Set the maximum number of steps in interaction chains considered when calculating net effects

## SN (SNI) data
data_sn_untransformed <- readr::read_csv("INPUT_SNI_Untransformed.csv") ## Used in the paper
data_sn_kenner <- readr::read_csv("INPUT_SNI_KennerTransformed.csv")

## BC (WVI) data
data_w_untransformed <- readr::read_csv("INPUT_WVI_Untransformed.csv") ## Used in the paper
data_w_kenner <- readr::read_csv("INPUT_WVI_KennerTransformed.csv")

## Combined network layout node coordinates
layout_coords_combined <- readr::read_csv("Network_NodeCoordinates.csv")
layout_coords_sn <- layout_coords_combined
layout_coords_bc <- layout_coords_combined

## Data processing
switch (location,
    "SN" = {
        switch (transform,
            "Untransformed" = {
                data <- data_sn_untransformed
                colnames(data)[which(colnames(data) == "otters_AOS")] = "otter_AOS"

                ## Added on 2023-5-23 so `data = dplyr::select(data, all_of(colnames(data_sn_kenner)))` works, because the kenner transformed version of the SN data uses the older bivalve group with more than just crassadoma_gigantea
                colnames(data)[which(colnames(data) == "crassadoma_gigantea")] = "bivalve"

                ## Modified on 2023-5-23 to include abalone which are not in data_sn_kenner (Ellen only updated the untransformed version)
                data = dplyr::select(data, all_of(colnames(data_sn_kenner)), abalone, pterygophera)
            },
            "Kenner" = {
                data <- data_sn_kenner
            }
        )
    },
    "BC" = {
        switch (transform,
            "Untransformed" = {
                data <- data_w_untransformed
                data = dplyr::select(data, all_of(colnames(data_w_kenner)))
            },
            "Kenner" = {
                data <- data_w_kenner
            }
        )
    },
    stop("ERROR: Unknown location. Must be `SN` or `BC`.")
)

## Missing baysizes, set to zero to prevent any apply crashes
data[is.na(data)] = 0

## Sort again just to be safe
switch (location,
    "SN" = {
        data_arranged <- dplyr::arrange(data, station, year)
    },
    "BC" = {
        data_arranged <- dplyr::arrange(data, SITE, YEAR)
    },
    stop("ERROR: Unknown location. Must be `SN` or `BC`.")
)

## Remove alternative kinds of otter data
switch (location,
    "SN" = {
        data <- dplyr::select(data_arranged, -c("otter_AOS", "otters2km_smooth", "ottersAOS_smooth", "otters_2km"))
        colnames(data)[which(colnames(data) == "otters_2km_density")] = "otters"
        data_otter_kind <- "Otters2kmDensity"
    },
    "BC" = {
        data <- dplyr::select(data_arranged, -c("OTTER_2km", "OTTER2km_smooth", "OTTER_baylevel", "OTTERbaydensity_smooth"))
        colnames(data)[which(colnames(data) == "OTTER_density")] = "OTTER"
        data_otter_kind <- "OTTERdensity"
    },
    stop("ERROR: Unknown location. Must be `SN` or `BC`.")
)

## If using Kenner transform, need to log otter data too
if (transform == "Kenner") {
    data_otter_kind = paste0(data_otter_kind, "-log")

    switch (location,
        "SN" = {
            data$otters = log(1 + data$otters)
        },
        "BC" = {
            data$OTTER = log(1 + data$OTTER)
        },
        stop("ERROR: Unknown location. Must be `SN` or `BC`.")
    )
}

## Apply unity normalization, and salt if desired
switch (location,
    "SN" = {
        data <- dplyr::select(data, -c("kelletia"))
        data <- dplyr::select(data, -c("megathura"))
        data <- dplyr::select(data, -c("year", "season", "period", "station", "otterlevel", "UsableArea"))
    },
    "BC" = {
        data <- dplyr::select(data, -c("YEAR", "SITE", "OTTERLEVEL", "BAYSIZE"))
    },
    stop("ERROR: Unknown location. Must be `SN` or `BC`.")
)

data <- apply(data, 2, function(x) {
    if (var(x) == 0) {
        return(runif(length(x)) * noise)
    } else {
        noise_data <- runif(length(x)) * noise
        salted_data <- x + noise_data
        scaled_data <- salted_data/max(abs(salted_data)) ## Add abs(), really just for ENSO, which already is higher in positive values so this is just to make a point and does not affect the running of this particular dataset

        if (normalization == "percent") {
            scaled_data = scaled_data * 100
        }

        return(scaled_data)
    }
}) %>% as_tibble()

switch (location,
    "SN" = {
        data <- dplyr::bind_cols(dplyr::select(data_arranged, c("year", "season", "period", "station", "otterlevel", "UsableArea")), data)
    },
    "BC" = {
        data <- dplyr::bind_cols(dplyr::select(data_arranged, c("YEAR", "SITE", "OTTERLEVEL", "BAYSIZE")), data)
    },
    stop("ERROR: Unknown location. Must be `SN` or `BC`.")
)

## Create lags, and remove observations with no lag (first observation at each location)
switch (location,
    "SN" = {
        stations <- data$station %>% unique() %>% sort()
    },
    "BC" = {
        stations <- data$SITE %>% unique() %>% sort()
    },
    stop("ERROR: Unknown location. Must be `SN` or `BC`.")
)

for (s in seq_along(stations)) {

    switch (location,
        "SN" = {
            data_temp <- dplyr::filter(data, station == stations[s])
        },
        "BC" = {
            data_temp <- dplyr::filter(data, SITE == stations[s])
        },
        stop("ERROR: Unknown location. Must be `SN` or `BC`.")
    )

    rows <- nrow(data_temp)
    data_t0 <- data_temp[2:rows,]
    data_t1 <- data_temp[1:(rows-1),]
    colnames(data_t1) = paste0(colnames(data_t1), "_t1")

    data_temp_lagged = dplyr::bind_cols(data_t0, data_t1)

    if (s == 1) {
        data_lagged <- data_temp_lagged
    } else {
        data_lagged = dplyr::bind_rows(data_lagged, data_temp_lagged)
    }

}

## Remove rows without an immediate lag (1 year since previous observation)
keep_id <- rep(NA, nrow(data_lagged))
keep_id[1] = TRUE

for (id in 2:length(keep_id)) {
    switch (location,
        "SN" = {
            same_station <- data_lagged$station[id] == data_lagged$station_t1[id]
            correct_lag <- data_lagged$year[id] == 1 + data_lagged$year_t1[id]
        },
        "BC" = {
            same_station <- data_lagged$SITE[id] == data_lagged$SITE_t1[id]
            correct_lag <- data_lagged$YEAR[id] == 1 + data_lagged$YEAR_t1[id]
        },
        stop("ERROR: Unknown location. Must be `SN` or `BC`.")
    )

    if ( all(same_station, correct_lag) ) {
        keep_id[id] = TRUE
    } else {
        keep_id[id] = FALSE
    }
}

if (!all(keep_id[1:5])) {
    print("WARNING: Lag cleanup might have failed.")
}

data_lagCleanup <- dplyr::filter(data_lagged, keep_id)
data = data_lagCleanup

## Rescale ENSO to have same expected variance as the rest of the variables in the [0,1] range
data$ENSO = data$ENSO * 0.5
data$ENSO_t1 = data$ENSO_t1 * 0.5

## Make names of urchin species identical between matrices and data
if (location == "SN") {
    col_id_sf <- which(colnames(data) == "strongylo_franciscanus")
    col_id_sp <- which(colnames(data) == "strongylo_purpuratus")
    colnames(data)[col_id_sf] = "s_franciscanus"
    colnames(data)[col_id_sp] = "s_purpuratus"

    col_id_sf <- which(colnames(data) == "strongylo_franciscanus_t1")
    col_id_sp <- which(colnames(data) == "strongylo_purpuratus_t1")
    colnames(data)[col_id_sf] = "s_franciscanus_t1"
    colnames(data)[col_id_sp] = "s_purpuratus_t1"

    col_id_patiria <- which(colnames(data) == "patiria_miniata")
    colnames(data)[col_id_patiria] = "patiria"

    col_id_patiria <- which(colnames(data) == "patiria_miniata_t1")
    colnames(data)[col_id_patiria] = "patiria_t1"

    col_id_pisaster <- which(colnames(data) == "pisaster_giganteus")
    colnames(data)[col_id_pisaster] = "pisaster"

    col_id_pisaster <- which(colnames(data) == "pisaster_giganteus_t1")
    colnames(data)[col_id_pisaster] = "pisaster_t1"

    col_id_cucum <- which(colnames(data) == "other_cucum")
    colnames(data)[col_id_cucum] = "other_cucs"

    col_id_cucum <- which(colnames(data) == "other_cucum_t1")
    colnames(data)[col_id_cucum] = "other_cucs_t1"

    col_id_cucum <- which(colnames(data) == "othercucum")
    colnames(data)[col_id_cucum] = "other_cucs"

    col_id_cucum <- which(colnames(data) == "othercucum_t1")
    colnames(data)[col_id_cucum] = "other_cucs_t1"

    col_id_enso <- which(colnames(data) == "ENSO")
    colnames(data)[col_id_enso] = "enso"

    col_id_enso <- which(colnames(data) == "ENSO_t1")
    colnames(data)[col_id_enso] = "enso_t1"

    col_id_stars <- which(colnames(data) == "otherstars")
    colnames(data)[col_id_stars] = "other_stars"

    col_id_stars <- which(colnames(data) == "otherstars_t1")
    colnames(data)[col_id_stars] = "other_stars_t1"
}

## Model from expert adjacency matrix
switch (location,
    "SN" = {
        matrix_model_raw <- read.table("matrix_model_SNI.csv", sep = ",", header = TRUE, row.names = 1)
    },
    "BC" = {
        matrix_model_raw <- read.table("matrix_model_WVI.csv", sep = ",", header = TRUE, row.names = 1)
    },
    stop("ERROR: Unknown location. Must be `SN` or `BC`.")
)

matrix_model <- matrix_model_raw
matrix_model = as.matrix(matrix_model)

## Store the full set of interactors
interactors_causes_master <- paste0(colnames(matrix_model), "_t1")
interactors_effects_master <- rownames(matrix_model)
interactors_all_master <- c(interactors_effects_master, substr(interactors_causes_master, 1, nchar(interactors_causes_master)-3)) %>% unique() %>% sort()

## Add state, here measured with PCA, to data
switch (location,
    "SN" = {
        pca_data <- data %>% dplyr::select(ends_with("_t1")) %>% dplyr::select(-c("year_t1", "season_t1", "period_t1", "station_t1", "otterlevel_t1", "UsableArea_t1"))
    },
    "BC" = {
        pca_data <- data %>% dplyr::select(ends_with("_t1")) %>% dplyr::select(-c("YEAR_t1", "SITE_t1", "OTTERLEVEL_t1", "BAYSIZE_t1"))
    },
    stop("ERROR: Unknown location. Must be `SN` or `BC`.")
)

pca_full <- stats::prcomp(pca_data)
data_weighting_pca <- pca_full$x %>% as.matrix()
pca_importance <- as.numeric(summary(pca_full)$importance[2,])
pca_importance = pca_importance/sum(pca_importance)
pca1 <- pca_full$x[, 1]
pca2 <- pca_full$x[, 2]

## Use UMAP instead of PCA
if (state == "umap") {
    umap_2d <- umap::umap(pca_data, n_components = 1)$layout
    pca1 <- umap_2d[,1]
    pca2 <- umap_2d[,2]
}

pca1_squared <- pca1^2
pca2_squared <- pca2^2

data = data %>% transmute(zca1 = pca1, zca1_squared = pca1_squared, zca2 = pca2, zca2_squared = pca2_squared, across(everything()))
model_data <- data


## Create a Primary Directory for each run of the model
dir.create(paste0(pd, "/Output_", date), showWarnings = FALSE) ## Hide warnings because every run of the model after the first will warn that the directory already exists
pd_foreach <- paste0(pd, "/Output_", date, "/", location, "_", tu_id)

## Create and move into pd_foreach
dir.create(pd_foreach)
dir.create(paste0(pd_foreach, "/Models"))
dir.create(paste0(pd_foreach, "/Matrices"))
dir.create(paste0(pd_foreach, "/Networks"))
dir.create(paste0(pd_foreach, "/Tables"))
dir.create(paste0(pd_foreach, "/Figures"))
setwd(pd_foreach)

matrix <- matrix_model

## Dimensions of the system
num_rows <- nrow(matrix)
num_cols <- ncol(matrix)

## Names of the variables in the system
names_rows <- rownames(matrix)
names_cols <- colnames(matrix)

## Separate matrices for the interactions and their bounds
matrix_model <- matrix(NA, nrow = num_rows, ncol = num_cols)
matrix_lowerBound <- matrix(NA, nrow = num_rows, ncol = num_cols)
matrix_upperBound <- matrix(NA, nrow = num_rows, ncol = num_cols)

## Add names of the input matrix to the model and bounds matrices
rownames(matrix_model) = names_rows
colnames(matrix_model) = names_cols

rownames(matrix_lowerBound) = names_rows
colnames(matrix_lowerBound) = names_cols

rownames(matrix_upperBound) = names_rows
colnames(matrix_upperBound) = names_cols

## Spread the information about how each interaction is bounded across the three kinds of matrices
matrix_model[matrix == "U"] = 1
matrix_model[matrix == "P"] = 1
matrix_model[matrix == "N"] = 1

## If having convergence issues can try using values close to zero instead of zero
matrix_lowerBound[matrix == "P"] = 0 #1e-10
matrix_upperBound[matrix == "N"] = 0 #-1e-10

## Initialize parallel backend
cluster <- makeCluster(3, outfile = "") ## Careful not to overwhelm your system! Three cores is probably safe but may take hours or days. Usually there are only a few models that take time to converge so even getting up to 10 cores should dramatically improve runtime.
registerDoParallel(cluster)

## Build the model
## Add equations one effect variable (row) at a time
model_fit_foreach <- foreach (row = 1:nrow(matrix_model), .packages = c("tibble", "dplyr", "piecewiseSEM", "lme4", "nlme", "Matrix", "optimx", "MuMIn")) %dopar% {

    ## Initialize vectors of interactors
    interactors <- {}
    interactors_causes <- {}
    interactors_effects <- {}

    tbl_interactions <- tibble(
        Effect = character(),
        Cause = character(),
        zca1 = integer(),
        zca1_squared = integer(),
        zca2 = integer(),
        zca2_squared = integer(),
        `zca1*zca2` = integer(),
        `zca1_squared*zca2_squared` = integer()
    )

    effect_row <- rownames(matrix_model)[row]

    if(verbose) { print(paste0("Row ", row, " of ", nrow(matrix_model), " | Effect = ", effect_row) ) }

    ## Check that the effect variable has any causes (non-empty columns)
    if (sum(matrix_model[row,], na.rm = TRUE) != 0) {

        ## All models are of the generic form `effect ~ cause`
        model_row = paste0("(", effect_row, " ~ ")

        ## Add the effect (row) variable's name to the growing list of interactors
        interactors = c(interactors, effect_row)
        interactors_effects = c(interactors_effects, effect_row)

        ## Add the cause (column) variable's name to the growing list of interactors
        causes <- colnames(matrix_model)[which(matrix_model[row,] == 1)]
        ## Save column names for creating bounds later
        causes_row <- causes
        causes = paste0(causes, "_t1")
        causes_row_model <- causes
        interactors = c(interactors, causes)
        interactors_causes = c(interactors_causes, causes)

        ## Combine the coefficient with the cause variable with the generic additive formula
        causes = paste0(causes, collapse = " + ")
        model_row = paste0(model_row, causes)

        model_row = paste0(model_row, ", data = model_data")

        #Fit a basic linear model and set up a deviance function (slightly trickier)
        lmmod <- paste0("lmmod <- lm", model_row, ")")
        eval(parse(text = lmmod))

        ## define fixed parameter bounds
        ## length = number of cause variables + 1 for the intercept, which comes first
        lowerBound <- c(-Inf, matrix_lowerBound[row, causes_row])
        upperBound <- c(Inf, matrix_upperBound[row, causes_row])

        ## Convert NA, meaning no bound, to infinity
        lowerBound[is.na(lowerBound)] = -Inf
        upperBound[is.na(upperBound)] = Inf

        ## Build Residual Sum of Squares model
        RSSfnmod <- paste0("RSSfnmod <- function(par) {  model_data$", effect_row, " - (par[1]")

        RSSfnmod_testdata <- RSSfnmod
        RSSfnmod_var <- list()
        par_counter <- 1

        ## Intraspecific causes most go first so their relation to state can be held constant in intraspecific estimations
        intraspecific_id <- which(causes_row_model == paste0(effect_row, "_t1"))
        causes_row_model = c(causes_row_model[intraspecific_id], causes_row_model) %>% unique() ## This works because unique() keeps the first instance

        ## Need to reorder lower and upper bounds with the changed order of the causes to keep the intraspecific causes first
        lowerBound_row <- c(lowerBound[1], lowerBound[1+intraspecific_id], lowerBound[-c(1, 1+intraspecific_id)])
        upperBound_row <- c(upperBound[1], upperBound[1+intraspecific_id], upperBound[-c(1, 1+intraspecific_id)])

        for (var in seq_along(causes_row_model)) {
            par_counter = par_counter + 1
            RSSfnmod_var_addition <- paste0("(par[", par_counter, "]")
            RSSfnmod_var_testdata_addition <- RSSfnmod_var_addition

            ## Intraspecific causes have no second term in the model (would repeat their first term)
            if (var == 1) {
                gmod_formula <- paste0(effect_row, " ~ ", causes_row_model[var], " * (zca1*zca2+zca1_squared*zca2_squared) - (zca1*zca2+zca1_squared*zca2_squared)")
            } else {
                if (intraspecific_interaction == "NO_EFFECT_OF_STATE") {
                    gmod_formula <- paste0(effect_row, " ~ ", causes_row_model[var], " * (zca1*zca2+zca1_squared*zca2_squared) + ", effect_row, "_t1 - (zca1*zca2+zca1_squared*zca2_squared)")
                } else {
                    gmod_formula <- paste0(effect_row, " ~ ", causes_row_model[var], " * (zca1*zca2+zca1_squared*zca2_squared) + ", effect_row, "_t1 * (", intraspecific_interaction, ") - (zca1*zca2+zca1_squared*zca2_squared)")
                }
            }

            gmod <- lm(as.formula(gmod_formula), data = model_data, na.action = "na.fail")

            if (var == 1) {
                dredgemod <- MuMIn::dredge(gmod, fixed = paste0(effect_row, "_t1"), beta = "none", rank = model_selection)
            } else {
                dredgemod <- MuMIn::dredge(gmod, beta = "none", rank = model_selection)
            }

            ## Make sure to include the intraspecific interactions, except for when defining this (e.g. spp ~ spp_t1*fun(state))
            if (var == 1) {
                topmod <- get.models(dredgemod, subset = 1)[[1]]$coefficients
            } else {
                necessary_cols <- startsWith(names(dredgemod), paste0(effect_row, "_t1:"))
                topmod_id <- which(rowSums(as.data.frame(is.na(dredgemod[, necessary_cols]))) == 0)[1] %>% as.numeric()
                topmod <- get.models(dredgemod, subset = topmod_id)[[1]]$coefficients

                ## Fixes the issue with names of interactions being reordered
                topmod_dredgemod <- as.matrix(dredgemod)[topmod_id,]
                topmod_dredgemod = topmod_dredgemod[!is.na(topmod_dredgemod)]
                topmod_dredgemod = topmod_dredgemod[1:length(topmod)]
                topmod = topmod_dredgemod
            }

            fixed_pars <- which(names(topmod) %in% c("(Intercept)", paste0(effect_row, "_t1"), causes_row_model[var]))

            ## Now remove intrascpecific interactions, except for when defining this (e.g. spp ~ spp_t1*fun(state))
            if (var != 1) {
                names_topmod <- names(topmod)
                for (name in seq_along(names_topmod)) {
                    if (startsWith(names_topmod[name], paste0(effect_row, "_t1"))) {
                        fixed_pars = c(fixed_pars, name)
                    }
                }
            }
            fixed_pars = unique(fixed_pars) ## Direct (not a function of state) intrascpecific interactions will be repeated

            topmod = topmod[-fixed_pars]

            tbl_topmod_var_addition <- {}

            if (length(topmod) > 0) {
                for (topmod_var in seq_along(topmod)) {
                    topmod_var_addition <- strsplit(names(topmod)[topmod_var], ":")[[1]][-1]

                    ## E.g. three-way interactions
                    if (length(topmod_var_addition) > 1) {
                        topmod_var_addition_modeldataadded <- paste(topmod_var_addition, collapse = "*model_data$")
                        topmod_var_addition = paste(topmod_var_addition, collapse = "*")
                    } else {
                        topmod_var_addition_modeldataadded <- topmod_var_addition
                    }

                    par_counter = par_counter + 1

                    RSSfnmod_var_testdata_addition = paste0(RSSfnmod_var_testdata_addition, " + par[", par_counter, "]*", topmod_var_addition)
                    RSSfnmod_var_addition = paste0(RSSfnmod_var_addition, " + par[", par_counter, "]*model_data$", topmod_var_addition_modeldataadded)

                    tbl_topmod_var_addition = paste0(tbl_topmod_var_addition, topmod_var_addition, " + ")

                }
            } else {
                ## Add the plus sign at the end so it will be trimmed off like the normal pca terms that accrue in the loop above
                tbl_topmod_var_addition = "NO_EFFECT_OF_STATE + "
            }

            RSSfnmod_var_addition = paste0(RSSfnmod_var_addition, ")")
            RSSfnmod_var_testdata_addition = paste0(RSSfnmod_var_testdata_addition, ")")

            RSSfnmod_var[[var]] = RSSfnmod_var_addition

            ## Add upper and lower bounds for each var
            RSSfnmod = paste0(RSSfnmod, " + model_data$", causes_row_model[var], "*", RSSfnmod_var_addition)

            ## Same as above but for the localparams estimation below, so without `model_data$`
            RSSfnmod_testdata = paste0(RSSfnmod_testdata, " + ", causes_row_model[var], "*", RSSfnmod_var_testdata_addition)

            tbl_topmod_var_addition = substr(tbl_topmod_var_addition, 1, nchar(tbl_topmod_var_addition)-3)

            ## Save the intraspecific interaction so can be constant in interspecific estimations
            if (var == 1) {
                intraspecific_interaction <- tbl_topmod_var_addition
            }

            if (intra_inter_priority == "inter") {
                intraspecific_interaction = "NO_EFFECT_OF_STATE"
            }

            tbl_topmod_var_addition_split <- strsplit(x = tbl_topmod_var_addition, split = " [+] ")[[1]]
            tbl_interactions_addition_estimates <- (colnames(tbl_interactions)[-c(1,2)] %in% tbl_topmod_var_addition_split) %>% as.matrix() %>% t() %>% as.data.frame()
            names(tbl_interactions_addition_estimates) = colnames(tbl_interactions)[-c(1,2)]

            tbl_interactions_addition_estimates = as.matrix(tbl_interactions_addition_estimates)[1,]

            ## This works because the FALSE entries becomes zeros
            tbl_interactions_addition_estimates_pre <- tbl_interactions_addition_estimates
            tbl_interactions_addition_estimates[tbl_interactions_addition_estimates] = topmod

            tbl_interactions_addition = cbind(
                tibble(
                    Effect = effect_row,
                    Cause = causes_row_model[var]
                ),
                as.list(tbl_interactions_addition_estimates)
            )

            tbl_interactions = dplyr::bind_rows(tbl_interactions, tbl_interactions_addition)
        }

        RSSfnmod = paste0(RSSfnmod, ")")
        RSSfnmod_testdata = paste0(RSSfnmod_testdata, ")")

        RSSfnmod_suffix <- " }"

        RSSfnmod = paste0(RSSfnmod, RSSfnmod_suffix)
        RSSfnmod_testdata = paste0(RSSfnmod_testdata, RSSfnmod_suffix)

        RSSfnmod_string <- RSSfnmod
        RSSfnmod_string_testdata <- RSSfnmod_testdata

        RSSfnmod = eval(parse(text = RSSfnmod))

        ## Use quadratic penalty on bounds violations to enforce expert knowledge about sign of interactions
        RSSfnmod_penalizedBounds <- function(par, verbose = FALSE, penalty = FALSE) {
            model_fit_FixedEffects <- RSSfnmod(par)

            model_fit =  sum(model_fit_FixedEffects^2)

            model_penalty <- 0
            model_penalty_raw <- model_penalty ## For keeping track of bounds violations
            for (v in seq_along(RSSfnmod_var)) {
                v_values <- eval(parse(text = RSSfnmod_var[[v]]))

                ## Need to add something in case no bounds are violated, so model_penalty + ___ doesn't throw an error
                model_penalty_addition <- 0 + sum(abs(v_values[v_values < lowerBound_row[-1][v]]))
                model_penalty_raw = model_penalty_raw + model_penalty_addition
                model_penalty = model_penalty + (model_penalty_addition)^2

                model_penalty_addition <- 0 + sum(abs(v_values[v_values > upperBound_row[-1][v]]))
                model_penalty_raw = model_penalty_raw + model_penalty_addition
                model_penalty = model_penalty + (model_penalty_addition)^2
            }

            ## Allows penalty to matter more so can ensure bounds are being observed
            model_penalty = model_penalty * model_penalty_upscale

            model_fit_penalty <- model_fit + model_penalty

            if (verbose) { print(paste0("Row = ", row, " | Fit = ", model_fit, " | Penalty = ", model_penalty, " | Raw Penalty = ", model_penalty_raw)) }

            if (penalty == TRUE) {
                return(model_penalty_raw)
            } else {
                return(model_fit_penalty)
            }
        }

        ## Estimate model parameters using penalized RSS via optimization
        modresult_iters <- list()
        for (iter in 1:optim_iters) {

            switch(optim_method,
                dfd = {
                    modresult <- optim(par = rnorm(par_counter), fn = RSSfnmod_penalizedBounds, method=c("BFGS"))
                },
                optimx = {
                    modresult_optimx <- optimx::optimx(par = runif(par_counter, max = 1e-5) + rep(-0.001, par_counter), fn = RSSfnmod_penalizedBounds, lower = rep(-optim_bound, par_counter), upper = rep(optim_bound, par_counter), method = c("L-BFGS-B"))
                    modresult <- list(par = modresult_optimx[1, 1:par_counter] %>% unlist())
                },
                gensa = {
                    modresult <- GenSA::GenSA(fn = RSSfnmod_penalizedBounds, lower = rep(-optim_bound, par_counter), upper = rep(optim_bound, par_counter))
                    modresult$convergence = NA
                },
                gensaoptim = {
                    modresult_GenSA <- GenSA::GenSA(fn = RSSfnmod_penalizedBounds, lower = rep(-optim_bound, par_counter), upper = rep(optim_bound, par_counter))
                    modresult <- optim(par = modresult_GenSA$par, fn = RSSfnmod_penalizedBounds, method=c("BFGS"))
                },
                optimSAoptim = {
                    modresult_optim <- optim(par = rnorm(par_counter), fn = RSSfnmod_penalizedBounds, method=c("BFGS"))
                    modresult_GenSA <- GenSA::GenSA(par = modresult_optim$par, fn = RSSfnmod_penalizedBounds, lower = rep(-optim_bound, par_counter), upper = rep(optim_bound, par_counter))
                    modresult <- optim(par = modresult_GenSA$par, fn = RSSfnmod_penalizedBounds, method=c("BFGS"))
                },
                {
                    stop("ERROR: Unkonwn `optim_method` variable.")
                }
            )

            modresult$penalty_raw = RSSfnmod_penalizedBounds(modresult$par, penalty = TRUE)

            modresult_iters[[iter]] = modresult

        }

        optim_values <- sapply(modresult_iters, function(x){x$value})
        optim_convergence <- sapply(modresult_iters, function(x){x$convergence})

        ## If none of the optimizations converged normally then use the best fit, but otherwise use the best fit from those that did converge
        if(sum(optim_convergence == 0, na.rm = TRUE) == 0) {
            optim_id <- which(optim_values == min(optim_values))[1]
        } else {
            optim_id_ignore <- which(optim_convergence != 0)
            optim_values[optim_id_ignore] = NA
            optim_id <- which(optim_values == min(optim_values, na.rm = TRUE))[1]
        }

        modresult <- modresult_iters[[optim_id]]

        if (verbose) { print(paste0("Convergence flag of best fit = ", modresult$convergence, " | All convergence flags = ", paste0(sapply(modresult_iters, function(x){x$convergence}), collapse = " "))) }
        if (verbose) { print(paste0("Best fit = ", modresult$value, " | All model fits = ", paste0(sapply(modresult_iters, function(x){x$value}), collapse = " "))) }

        ## Just to get the print statement to check if the bounds are being observed
        RSSfnmod_penalizedBounds(modresult$par, verbose = TRUE)

        model_fit_FixedEffects <- RSSfnmod(modresult$par)

        model_fit =  sum(model_fit_FixedEffects^2)

        RsquaredValue <- 1 - ( model_fit / sum((model_data[[effect_row]] - mean(model_data[[effect_row]]))^2) )
        df_n <- nrow(model_data)
        df_p <- length(modresult$par)-1

        RSSfnmod_string_testdata_row = RSSfnmod_string_testdata
        RSSfnmod_suffix_row = RSSfnmod_suffix
        par_counter_row = par_counter
        RSquaredPerSpecies_row = RsquaredValue
        AdjustedRSquaredPerSpecies_row = 1 - ( (1-RsquaredValue) * ( (df_n-1) / (df_n-df_p-1) ) )

        output <- list(
            modresult_list = modresult,
            modresult_convergence = modresult$convergence,
            modresult_penalty = modresult$penalty_raw,
            RSSfnmod_string_testdata_list = RSSfnmod_string_testdata_row,
            RSSfnmod_suffix_list = RSSfnmod_suffix_row,
            lmmod_list = lmmod,
            par_counter_list = par_counter_row,
            RSquaredPerSpecies = RSquaredPerSpecies_row,
            AdjustedRSquaredPerSpecies = AdjustedRSquaredPerSpecies_row,
            tbl_interactions = tbl_interactions
        )

    ## Row variable has no causes
    } else {
        output <- list(
            modresult_list = NA,
            modresult_convergence = NA,
            modresult_penalty = NA,
            RSSfnmod_string_testdata_list = NA,
            RSSfnmod_suffix_list = NA,
            lmmod_list = NA,
            par_counter_list = NA,
            RSquaredPerSpecies = NA,
            AdjustedRSquaredPerSpecies = NA,
            tbl_interactions = tbl_interactions
        )
    }

    return(output)
}

## Assign elements of output list to their own variables
modresult_list <- lapply(model_fit_foreach, function(x){x$modresult_list})
modresult_convergence <- lapply(model_fit_foreach, function(x){x$modresult_convergence})
modresult_penalty <- lapply(model_fit_foreach, function(x){x$modresult_penalty})
RSSfnmod_string_testdata_list <- lapply(model_fit_foreach, function(x){x$RSSfnmod_string_testdata_list})
RSSfnmod_suffix_list <- lapply(model_fit_foreach, function(x){x$RSSfnmod_suffix_list})
lmmod_list <- lapply(model_fit_foreach, function(x){x$lmmod_list})
par_counter_list <- lapply(model_fit_foreach, function(x){x$par_counter_list})
RSquaredPerSpecies <- sapply(model_fit_foreach, function(x){x$RSquaredPerSpecies})
AdjustedRSquaredPerSpecies <- sapply(model_fit_foreach, function(x){x$AdjustedRSquaredPerSpecies})
tbl_interactions <- do.call(rbind, lapply(model_fit_foreach, function(x){x$tbl_interactions}))

RSquaredPerSpecies = tibble::enframe(RSquaredPerSpecies)
colnames(RSquaredPerSpecies) = c("EffectVariable_RowInMatrixModel", "R_squared")

AdjustedRSquaredPerSpecies = tibble::enframe(AdjustedRSquaredPerSpecies)
colnames(AdjustedRSquaredPerSpecies) = c("EffectVariable_RowInMatrixModel", "Adjusted_R_Squared")

RSquaredPerSpecies = dplyr::full_join(RSquaredPerSpecies, AdjustedRSquaredPerSpecies, by = "EffectVariable_RowInMatrixModel")
RSquaredPerSpecies = dplyr::transmute(RSquaredPerSpecies, EffectSpecies = rownames(matrix_model), across(everything()))

setwd(paste0(pd_foreach, "/Tables"))
    write.table(tbl_interactions, file = paste0("tbl_interactions_", location, ".csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    write.table(RSquaredPerSpecies, file = paste0("RsquaredPerEffect_", location, ".csv"), sep = ",", col.names = TRUE, row.names = FALSE)
setwd(pd_foreach)

## Save all parameters to help compare models
par_list <- lapply(modresult_list, \(x){x[[1]]})
saveRDS(object = par_list, file = "par_list.rds")

## Save equations (in string form)
saveRDS(object = RSSfnmod_string_testdata_list, file = "RSSfnmod_string_testdata_list.rds")

## Save the models so can infer interactions at any point in the state space
setwd(paste0(pd_foreach, "/Tables"))
    saveRDS(object = lmmod_list, file = "lmmod_list.rds")
setwd(pd_foreach)














## Now use state-based model to estimate interactions for each observation

## For making movies of networks changing
figure_names <- rep(NA, nrow(data))

## NOTE: This parallel loop will use the same backend cluster created for the parallel loop above
model_estimate_foreach <- foreach (observation = 1:nrow(data), .combine = "rbind", .packages = c("tibble", "dplyr", "piecewiseSEM", "lme4", "optimx", "ggraph", "tidygraph", "magrittr")) %dopar% {
    print(paste0("psem observation ", observation))
    setwd(pd_foreach)

    ## Initialize vectors of interactors
    model <- "model <- piecewiseSEM::psem("

    for (row in 1:nrow(matrix_model)) {

        ## Check that the effect variable has any causes (non-empty columns)
        if (sum(matrix_model[row,], na.rm = TRUE) != 0) {

            causes <- colnames(matrix_model)[which(matrix_model[row,] == 1)]
            causes = paste0(causes, "_t1")
            causes_row_model <- causes

            effect_row <- rownames(matrix_model)[row]

            ## Intraspecific causes most go first so their relation to state can be held constant in intraspecific estimations
            intraspecific_id <- which(causes_row_model == paste0(effect_row, "_t1"))
            causes_row_model = c(causes_row_model[intraspecific_id], causes_row_model) %>% unique() ## This works because unique() keeps the first instance

            zca1val <- model_data$zca1[observation]
            zca1squaredval <- model_data$zca1_squared[observation]
            zca2val <- model_data$zca2[observation]
            zca2squaredval <- model_data$zca2_squared[observation]

            testdata <- matrix(rep(c(rep(0,length(causes_row_model)), zca1val, zca1squaredval, zca2val, zca2squaredval), 1+length(causes_row_model)), ncol = 4+length(causes_row_model), byrow = TRUE)
            diag(testdata[-1, -c(ncol(testdata)-1, ncol(testdata))]) = 1
            colnames(testdata) = c(causes_row_model, "zca1", "zca1_squared", "zca2", "zca2_squared")
            testdata = as.data.frame(testdata)

            modresult = modresult_list[[row]]

            par <- as.numeric(modresult$par[1:par_counter_list[[row]]])

            prefix_length <- 43 + nchar(effect_row)
            estis <- with(testdata, eval(parse(text = substr(RSSfnmod_string_testdata_list[[row]], prefix_length, nchar(RSSfnmod_string_testdata_list[[row]])-2 )    )))

            ## These are the estimated local interaction, after accounting for the intercept
            localparams <- c(estis[1], estis[2:length(estis)]-estis[1])

            # finally, sub into the global model:
            lmmod = lmmod_list[[row]]
            names(localparams) = c("(Intercept)", causes_row_model)

            lmmod$coefficients = localparams[names(lmmod$coefficients)]

            ## Save this effect variable's model to feed into the piecewise SEM
            model_row_name <- paste0("model_", row)
            assign(model_row_name, lmmod)
            model = paste0(model, model_row_name, ", ")

        }
    }

    ## Remove final semicolon
    model = substr(model, 1, (nchar(model)-2))

    ## Add final parenthesis
    model = paste0(model, ")")

    eval(parse(text = model))

    ## Save model
    setwd(paste0(pd_foreach, "/Models"))

    saveRDS(object = model,
            file = paste0("StateBasedSEM_", location, "_Observation-", observation ,".rds"))
    setwd(pd_foreach)

    coefficients <- piecewiseSEM::coefs(model)[, c(1:3)] %>% as_tibble()
    coefficients$Response = coefficients$Response %>% as.character()
    coefficients$Predictor = coefficients$Predictor %>% as.character()
    coefficients$Estimate = coefficients$Estimate %>% as.numeric()

    sem_mat <- matrix(0, nrow = nrow(matrix_model_raw), ncol = ncol(matrix_model_raw))
    colnames(sem_mat) = paste0(colnames(matrix_model_raw), "_t1")
    rownames(sem_mat) = rownames(matrix_model_raw)

    for (row in 1:nrow(coefficients)) {
        sem_mat[coefficients$Response[row], coefficients$Predictor[row]] = coefficients$Estimate[row]
    }

    ## State-dependency allows nonzero interaction strengths when the causer was unobserved, so manually set these interactions to zero
    if (removedabsentinteractors == "RemovedAbsentInteractors") {
        for (cause in colnames(sem_mat)) {
            if (unlist(data %>% dplyr::select(all_of(cause)))[observation] <= 1e-5) {
                sem_mat[,cause] = 0
            }
        }

        for (effect in rownames(sem_mat)) {
            if (unlist(data %>% dplyr::select(all_of(effect)))[observation] <= 1e-5) {
                sem_mat[effect,] = 0
            }
        }
    }

    ## Force expert knowledge bounds on interactions (positive or negative) to be observed
    if (model_penalty_enforcement == "ForcedBounds") {

        ## The NAs will NOT be replaced with zero here, but just ignored (R logic)
        sem_mat[sem_mat > matrix_upperBound] = 0
        sem_mat[sem_mat < matrix_lowerBound] = 0
    }

    setwd(paste0(pd_foreach, "/Matrices"))

    utils::write.table(sem_mat,
        file = paste0("Matrix_", location, "_Observation-", observation, ".csv"),
        sep = ",",
        col.names = TRUE,
        row.names = TRUE
    )
    setwd(pd_foreach)

    ## Create linked list to plot network
    ll <- sem_mat %>% as.matrix() %>% reshape2::melt() %>% as_tibble()
    graph_table <- tibble(From = sapply(ll$Var2 %>% as.character(), function(x){substr(x, 1, nchar(x)-3)}),
                        To = ll$Var1 %>% as.character(),
                        Strength = ll$value,
                        Sign = NA,
                        Color = NA,
                        Group = NA)

    for (i in 1:nrow(graph_table)) {
        graph_table$Group[i] = 1
        if (graph_table$Strength[i] > 0) {
            graph_table$Sign[i] = 1
            graph_table$Color[i] = "firebrick1"
        } else if (graph_table$Strength[i] == 0) {
            graph_table$Sign[i] = 0
            graph_table$Color[i] = "white"
        } else {
            graph_table$Sign[i] = -1
            graph_table$Color[i] = "dodgerblue"
        }
    }

    graph_table = dplyr::arrange(graph_table, From, To)
    node_names_all <- c(graph_table$From, graph_table$To) %>% unique() %>% sort()
    graph_table = dplyr::filter(graph_table, Strength != 0)
    node_names <- c(graph_table$From, graph_table$To) %>% unique() %>% sort()
    node_names_filtered <- node_names_all[!(node_names_all %in% node_names)]

    for (name in 1:length(node_names_all)) {
        new_row <- tibble(From = node_names_all[name], To = node_names_all[name], Strength = 2, Sign = 0, Color = "white", Group = NA)
        graph_table = dplyr::bind_rows(graph_table, new_row)
    }

    ## Remove tiny interaction strengths, which could be tolerance of the bounds on known (expert opinions) absent interactions
    graph_table = dplyr::filter(graph_table, abs(Strength) > 1e-5)
    graph_table = dplyr::arrange(graph_table, From, To)

    graph_table = dplyr::mutate(graph_table, Interaction_Strength = abs(Strength), Interaction_Sign = as.factor(as.character(Sign)))

    ## Swap in compact names
    if (location == "SN") {
        graph_table$From = sapply(graph_table$From, function(x){layout_coords_sn$Name[which(layout_coords_sn$Name_Data == x)]})
        graph_table$To = sapply(graph_table$To, function(x){layout_coords_sn$Name[which(layout_coords_sn$Name_Data == x)]})

        layout_coords <- layout_coords_sn
    } else if (location == "BC") {
        graph_table$From = sapply(graph_table$From, function(x){layout_coords_bc$Name[which(layout_coords_bc$Name_Data == x)]})
        graph_table$To = sapply(graph_table$To, function(x){layout_coords_bc$Name[which(layout_coords_bc$Name_Data == x)]})

        layout_coords <- layout_coords_bc
    } else {
        stop("ERROR: Unknown location. Must be `SN` or `BC`.")
    }

    graph_table %<>% dplyr::filter(From != "ENSO", To != "ENSO")
    graph_table %<>% dplyr::filter(From != "enso", To != "enso")

    layout_coords %<>% dplyr::arrange(match(Name, names((graph_table %>% dplyr::arrange(From, To) %>% tidygraph::as_tbl_graph())[1])))

    layout_coords %<>%
        dplyr::filter(Name != "ENSO") %>%
        dplyr::filter(Name != "enso")

    layout_coords %<>% dplyr::filter(Location == location) ## Added for this combined version (same layout of SN and BC) because ggraph layouts cannot have duplicate node rows

    graph_table_tidy <- graph_table %>% tidygraph::as_tbl_graph()

    lc_align <- graph_table_tidy %>% as_tibble()
    colnames(lc_align) = "Name"
    layout_coords = dplyr::left_join(lc_align, layout_coords, by = "Name")

    lc <- layout_coords %>%
        dplyr::select(c("X_coord", "Y_coord")) %>%
        as.matrix()

    if (all(graph_table$Sign == 1)) {
        color_high <- "dodgerblue"
        color_low <- "dodgerblue"
    } else if (all(graph_table$Sign == -1)) {
        color_high <- "firebrick1"
        color_low <- "firebrick1"
    } else {
        color_high <- "dodgerblue"
        color_low <- "firebrick1"
    }

    fig_network_title <- switch (location,
        "SN" = {
            paste0(data$station[observation], " | ", data$year[observation], paste0(" | Otters = ", round(data$otters_t1[observation], 2)*100, "% -> ", round(data$otters[observation], 2)*100, "% | PCA PSEM" ))
        },
        "BC" = {
            paste0(data$SITE[observation], " | ", data$YEAR[observation], paste0(" | Otters = ", round(data$OTTER_t1[observation], 2)*100, "% -> ", round(data$OTTER[observation], 2)*100, "% | PCA PSEM" ))
        },
        stop("ERROR: Unknown location. Must be `SN` or `BC`.")
    )

    png_name <- paste0("Network_", location, "_Observation-", observation, ".png")
    figure_names[observation] = png_name
    setwd(paste0(pd_foreach, "/Networks"))
        grDevices::png(filename = png_name, width = 1000, height = 1000) #, units = "cm", res = 500)
            fig <- ggraph(graph_table_tidy, layout = lc) +
                    theme_void() +
                    geom_edge_loop(aes(width = Interaction_Strength,
                                    color = Interaction_Sign,
                                    span = 150,
                                    direction = 180,
                                    strength = 0.1,
                                    ),
                                arrow = arrow(type = "closed",
                                            length = unit(3, "mm")
                                            )) +
                    geom_node_point(shape = 21, size = 20, fill = "white", alpha = 1) +
                    geom_edge_arc(aes(width = Interaction_Strength,
                                    color = Interaction_Sign,
                                    start_cap = label_rect(node1.name),
                                    end_cap = label_rect(node2.name)
                                    ),
                                start_cap = circle(8, "mm"),
                                end_cap = circle(8, "mm"),
                                strength = 0.01,
                                arrow = arrow(type = "closed",
                                            length = unit(2, "mm")
                                            )) +
                    geom_node_point(shape = 21, size = 20, fill = layout_coords$Color, alpha = 0.5) +
                    geom_node_text(aes(label = name),
                                size = 5) +
                    scale_edge_colour_manual(values = c(color_low, "white", color_high)) +
                    theme(plot.margin = grid::unit(c(10,10,10,10), "mm"),
                        legend.position = "bottom",
                        legend.justification = "center",
                        legend.margin = margin(50,50,50,50)) +
                    scale_edge_width(range = c(0.01, 5)) +
                    coord_cartesian(clip = "off") +
                    labs(title = fig_network_title)
        print(fig)
    dev.off()
    setwd(pd_foreach)

    melted_mat <- reshape2::melt(as.matrix(sem_mat)) %>% as_tibble()
    melted_mat = dplyr::mutate(melted_mat, interaction = paste0(melted_mat[,1] %>% unlist(), ".", melted_mat[,2] %>% unlist()))
    melted_mat = dplyr::select(melted_mat, interaction, value)

    output_wide <- melted_mat$value
    names(output_wide) = melted_mat$interaction

    return(output_wide)
}

output_wide <- model_estimate_foreach %>% as_tibble()

switch (location,
    "SN" = {
        metadata <- dplyr::select(data, c("year", "season", "period", "station", "otterlevel", "UsableArea"))
        metadata_addition <- dplyr::select(data, -c("year", "season", "period", "station", "otterlevel", "UsableArea"))
    },
    "BC" = {
        metadata <- dplyr::select(data, c("YEAR", "SITE", "OTTERLEVEL", "BAYSIZE"))
        metadata_addition <- dplyr::select(data, -c("YEAR", "SITE", "OTTERLEVEL", "BAYSIZE"))
    },
    stop("ERROR: Unknown location. Must be `SN` or `BC`.")
)

ordination_input <- dplyr::bind_cols(metadata, output_wide)
regular_input <- dplyr::bind_cols(metadata, metadata_addition)
StateBasedSEM <- dplyr::bind_cols(data, output_wide)

setwd(paste0(pd_foreach, "/Tables"))
    write.table(StateBasedSEM, file = paste0("StateBasedSEM_", location, ".csv"), sep = ",", col.names = TRUE, row.names = FALSE)
setwd(pd)

## Finally, save metadata
metadata_output <- tibble(
    tu_id = tu_id,
    date = date,
    location = location,
    data_otter_kind = data_otter_kind,
    transform = transform,
    normalization = normalization,
    model_selection = model_selection,
    intra_inter_priority = intra_inter_priority,
    model_kind = model_kind,
    model_type = model_type,
    model_penalty_upscale = model_penalty_upscale,
    optim_method = optim_method,
    optim_iters = optim_iters,
    optim_bound = optim_bound,
    modresult_convergence = paste0(modresult_convergence, collapse = " "),
    modresult_penalty = paste0(modresult_penalty, collapse = " "),
    intercept = TRUE,
    nevermissing = TRUE,
    state = state,
    weighting_scheme = weighting_scheme,
    removedabsentinteractors = removedabsentinteractors,
    model_penalty_enforcement = model_penalty_enforcement
)

setwd(pd_foreach)
write.table(metadata_output, file = paste0("metadata_", tu_id, ".csv"), sep = ",", col.names = TRUE, row.names = FALSE)














## Now calculate net effects by tracing chains of interactions

data <- StateBasedSEM

setwd(paste0(pd_foreach, "/Matrices"))
mat_files <- list.files(pattern = paste0("Matrix_", location, "_Observation"))

matrices <- list()
ids <- 1:nrow(data)

dir.create(paste0(pd_foreach, "/Matrices/NetEffects-", cutoff), showWarnings = FALSE)

## NOTE: This parallel loop will use the same backend cluster created for the parallel loops above
net_effects <- foreach (m = ids, .packages = c("tibble", "dplyr", "reshape2", "igraph", "reticulate")) %dopar% {
    if (verbose) { print(m) }

    file <- paste0("Matrix_", location, "_Observation-", m, ".csv")

    setwd(paste0(pd_foreach, "/Matrices"))

    sem_mat = read.table(file, sep = ",")

    ## Sort rows and columns alphabetically
    sem_mat = sem_mat[sem_mat %>% rownames() %>% sort(), sem_mat %>% colnames() %>% sort()]

    ## Create vector of interactors, without the column "_t1" suffix
    interactors <- c(rownames(sem_mat), stringr::str_sub(sem_mat %>% colnames, end=-4)) %>% unique() %>% sort()

    mat_complexity = sum(sem_mat != 0, na.rm = TRUE) / (nrow(sem_mat) * ncol(sem_mat))

    net_effects_row <- net_mat(direct_effects = sem_mat, cutoff = cutoff, rowVScol = "col", node_specific = FALSE, verbose = FALSE)

    melted_mat <- reshape2::melt(as.matrix(net_effects_row)) %>% as_tibble()

    melted_mat = dplyr::mutate(melted_mat, interaction = paste0(melted_mat[,1] %>% unlist(), ".", melted_mat[,2] %>% unlist()))
    melted_mat = dplyr::select(melted_mat, interaction, value)
    colnames(melted_mat)[2] = m %>% as.character()

    return(
        list(
            interaction_names = melted_mat[[1]],
            interaction_values = melted_mat[[2]]
        )
    )

}

## Shut down parallel backend and free up resources
parallel::stopCluster(cluster)
# closeAllConnections() ## WARNING: Only use if having trouble freeing up cores. This could close other jobs running in parallel.

net_ts <- do.call(rbind, lapply(net_effects, function(x){x[[2]]})) %>% as_tibble()
colnames(net_ts) = net_effects[[1]][1] %>% unlist()
colnames(net_ts) = paste0("net_", colnames(net_ts))

NetEffects_OUTPUT <- dplyr::bind_cols(data, net_ts)

setwd(paste0(pd_foreach, "/Tables"))
write.table(NetEffects_OUTPUT,
            file = paste0("NetEffects_", location, "_cutoff-", cutoff, ".csv"),
            sep = ",",
            col.names = TRUE,
            row.names = FALSE)
setwd(pd_foreach)
