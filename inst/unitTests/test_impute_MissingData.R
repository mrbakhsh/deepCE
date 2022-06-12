test_impute_MissingData = function() {
  data("HelaCE")
  Imputed_data <-
    impute_MissingData(HelaCE,bound_left = 2,bound_right = 2)
  checkEqualsNumeric(sum(as.matrix(Imputed_data)), 335366.9, tolerance = 1e-2)
}
