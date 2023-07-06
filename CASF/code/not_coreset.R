library(caret)

generate_cv_fold_indices <- function(n, folds){
  set.seed(13579)
  split(sample(1:n), rep(1:folds, length = n))
}

estimate_regressor_performance <- function(y_actual, y_predicted){
  rmse      <- sqrt(sum((y_actual - y_predicted)^2) / length(y_actual))
  pearsonr  <- cor(y_actual, y_predicted, method = "pearson")
  rsquared  <- 1 - (sum((y_actual - y_predicted)^2) / 
                    sum((y_actual - mean(y_actual))^2))
  c(pearsonr, rmse, rsquared)
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

read_data <- function(){
  amino <- read.csv(file = "../input/2019/refined_2019_aa.csv", header = TRUE, row.names = 1)
  ligand <- read.csv(file = "../input/2019/CASF_2019_ligprop.csv", header = TRUE, row.names = 1)
  activities <- read.csv(file = "../input/2019/PDBind_2019_activities.csv", header = FALSE, row.names = 1)

  samples <- intersect(rownames(activities), rownames(amino))

  y <- activities[samples, ]
  names(y) <- samples

  amino <- amino[samples, ]
  ligand <- ligand[samples, ]

  list(amino_acid = amino, ligand_properties = ligand, y = y)
}

build_and_predict <- function(x_train, x_test, y_train, y_test){
  model <- build_xgb_regressor(x_train, y_train)
  y_hat <- predict(model, x_test)
  estimate_regressor_performance(y_test, y_hat)
}

model_data <- function(data, fold_idx){
  y_train <- data$y[-fold_idx]
  y_test  <- data$y[fold_idx]

  train_samples <- names(y_train)
  test_samples  <- names(y_test)

  # Amino
  x_train_amino <- data$amino_acid[train_samples, ]
  x_test_amino  <- data$amino_acid[test_samples, ]
  amino_perf    <- build_and_predict(x_train_amino, x_test_amino, y_train, y_test)
  
  #------------------------------------------------------------------
  # RDKit
  x_train_ligand <- data$ligand_properties[train_samples, ]
  x_test_ligand  <- data$ligand_properties[test_samples, ]
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


main <- function(nfolds = 5){
  data <- read_data()
  folds <- generate_cv_fold_indices(length(data$y), nfolds)

  for (i in 1:nfolds){
    message(i)
    fold_idx <- folds[[i]]
    results <- model_data(data, fold_idx)
    write.csv(results, file = paste0("../output/2019_fold_", i, ".csv"))
  }
}

main()