library(caret)

estimate_regressor_performance <- function(y_actual, y_predicted){
  rmse      <- sqrt(sum((y_actual - y_predicted)^2) / length(y_actual))
  pearsonr  <- cor(y_actual, y_predicted, method = "pearson")
  c(pearsonr, rmse)
}

build_xgb_regressor <- function(x, y){
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    allowParallel = FALSE
  )
  
  grid <- expand.grid(
    nrounds = c(500, 1000, 1500, 2000),
    max_depth = c(2, 4, 6, 8),
    eta = c(0.001, 0.01, 0.1, 0.2, 0.3),
    gamma = 0.05,
    subsample = 0.5,
    colsample_bytree = 1,
    min_child_weight = 1
  )

  set.seed(47567)
  suppressWarnings(
    xgb_fit <- train(
      x = x, 
      y = y,
      method = "xgbTree",
      eval_metric = "rmse",
      verbose = 0, 
      trControl = ctrl,
      tuneGrid = grid,
      objective = "reg:squarederror"
    )
  )

  xgb_fit
}

read_year_data <- function(year){
  x_amino_acid <- read.csv(
    file = paste0("../input/", year, "/CASF_", year, "_aa.csv"),
    header = TRUE, 
    row.names = 1)
  x_ligand_properties <- read.csv(
    file = paste0("../input/", year, "/CASF_", year, "_ligprop.csv"),
    header = TRUE, 
    row.names = 1)
  y_train <- read.csv(
    file = paste0("../input/", year, "/CASF_", year, "_activities_train.csv"),
    header = FALSE, 
    row.names = 1)
  y_test <- read.csv(
    file = paste0("../input/", year, "/CASF_", year, "_activities_test.csv"),
    header = FALSE, 
    row.names = 1)
  
  train_samples <- rownames(y_train)
  y_train <- y_train[, 1]
  names(y_train) <- train_samples

  test_samples <- rownames(y_test)
  y_test <- y_test[, 1]
  names(y_test) <- test_samples

  list(amino_acid = x_amino_acid, ligand_properties = x_ligand_properties,
       y_train = y_train, y_test = y_test)
}

build_and_predict <- function(x_train, x_test, y_train, y_test){
  model <- build_xgb_regressor(x_train, y_train)
  y_hat <- predict(model, x_test)
  estimate_regressor_performance(y_test, y_hat)
}

model_data <- function(year_data){
  y_train <- year_data$y_train
  y_test  <- year_data$y_test

  train_samples <- names(y_train)
  test_samples  <- names(y_test)

  # Amino
  x_train_amino <- year_data$amino_acid[train_samples, ]
  x_test_amino  <- year_data$amino_acid[test_samples, ]
  amino_perf    <- build_and_predict(x_train_amino, x_test_amino, y_train, y_test)
  
  #------------------------------------------------------------------
  # RDKit
  x_train_ligand <- year_data$ligand_properties[train_samples, ]
  x_test_ligand  <- year_data$ligand_properties[test_samples, ]
  ligand_perf    <- build_and_predict(x_train_ligand, x_test_ligand, y_train, y_test)

  #------------------------------------------------------------------
  # RDKit + Amino
  x_train_ligand_amino <- cbind(x_train_ligand, x_train_amino)
  x_test_ligand_amino  <- cbind(x_test_ligand, x_test_amino)
  ligand_amino_perf    <- build_and_predict(
    x_train_ligand_amino, 
    x_test_ligand_amino, 
    y_train, y_test)

  #------------------------------------------------------------------

  performance <- rbind(amino_perf, ligand_perf, ligand_amino_perf)
  rownames(performance) <- c("Amino", "Ligand", "Amino+Ligand")
  colnames(performance) <- c("PearsonR", "RMSE")

  performance
}

main <- function(){
  years <- c(2007, 2013, 2016)
  for (year in years){
    year_data <- read_year_data(year)
    message(year)
    results <-  model_data(year_data)
    write.csv(results, file = paste0("../output/", year, ".csv"))
  }
}

main()