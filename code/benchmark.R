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

perform_experiments <- function(x, y, nfolds = 5){
  results <- matrix(nrow = 5, ncol = 3)
  colnames(results) <- c("PearsonR", "RMSE", "RSquared")

  folds <- generate_cv_fold_indices(length(y), nfolds)
  for (i in 1:nfolds){
    message(paste0("  fold: ", i))
    fold <- folds[[i]]
    x_train <- x[-fold, ]
    y_train <- y[-fold]
    x_test  <- x[fold, ]
    y_test  <- y[fold]

    model <- build_xgb_regressor(x_train, y_train)
    y_hat <- predict(model, x_test)
    results[i, ] <- estimate_regressor_performance(y_test, y_hat)
  }

  results
}

write_results <- function(results, exp_index, dataset){
  write.csv(
    results, 
    file = paste0("../output/experiment", exp_index, "_", dataset, ".csv"),
    row.names = FALSE)
}

main <- function(){
  for (exp_index in 1:6){
    message(paste0("Experiment dataset ", exp_index))
    expdf <- read.csv(file = paste0("../input/experiment", exp_index, ".csv"),
                      header = TRUE, row.names = 1)

    y <- expdf[, 1]

    # Ligand
    message("Ligand")
    x_ligand <- expdf[, 2:19]
    ligand_results <- perform_experiments(x_ligand, y)
    write_results(ligand_results, exp_index, "ligand")

    # Protein
    message("Protein")
    x_protein <- expdf[, 20:39]
    protein_results <- perform_experiments(x_protein, y)
    write_results(protein_results, exp_index, "protein")
    
    # Ligand and protein
    message("Ligand and Protein")
    x_ligand_protein <- expdf[, 2:39]
    ligand_protein_results <- perform_experiments(x_ligand_protein, y)
    write_results(ligand_protein_results, exp_index, "ligand_and_protein")
  }
}

main()