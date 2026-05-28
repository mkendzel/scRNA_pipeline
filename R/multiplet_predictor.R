#' Predict Multiplet Rate from Cells Recovered
#'
#' Predicts multiplet rates based on cell recovery counts using
#' a pre-defined linear model calibrated with internal data. Modeled data provided by
#' 10X Genomics Chromium Next GEM single Cell 5' reagent kits v2 (dual index).
#' Reference manual: CG000330 Rev G
#'
#' @param new_data A data frame containing new samples to predict. Must contain
#'   a column named 'sample' with sample identifiers and a column named 
#'   'CellsRecovered' with cell recovery counts.
#'
#' @return A list with two components:
#' \itemize{
#'   \item \code{detailed_predictions}: A data frame containing the original
#'         sample data with predicted multiplet rates and confidence intervals
#'   \item \code{quick_lookup_matrix}: A matrix with predicted values, lower,
#'         and upper confidence bounds, with sample names as row names
#' }
#'
#' @examples
#' # Create new data for prediction
#' new_samples <- data.frame(
#'   sample = c("mock_1", "mock_2", "TX_1"),
#'   CellsRecovered = c(19649, 30123, 25354)
#' )
#'
#' # Predict multiplet rates
#' results <- predict_multiplet_rate(new_samples)
#' 
#' # View detailed predictions
#' print(results$detailed_predictions)
#' 
#' # View quick lookup matrix
#' print(results$quick_lookup_matrix)
#'
#' @export

predict_multiplet_rate <- function(new_data) {
  # Calibration data (internal to the function)
  calibration_data <- data.frame(
    MultipletRate = c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9, 4.6, 5.4, 6.1, 6.9, 7.6),
    CellsLoaded = c(870, 1700, 3500, 5300, 7000, 8700, 10500, 12200, 14000, 15700, 17400),
    CellsRecovered = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
  )
  calibration_data$MultipletRate <- calibration_data$MultipletRate / 100
  
  # Fit the model
  model_simple <- lm(MultipletRate ~ CellsRecovered, data = calibration_data)
  
  # Make predictions
  predictions <- cbind(new_data, as.data.frame(
    predict(model_simple, newdata = new_data["CellsRecovered"], interval = "confidence")
  ))
  
  # Create named matrix
  pred_mat <- as.matrix(predictions[, c("fit", "lwr", "upr")])
  rownames(pred_mat) <- predictions$sample
  
  return(list(
    detailed_predictions = predictions,
    quick_lookup_matrix = pred_mat
  ))
}