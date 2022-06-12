test_scaling_m = function() {
  data("HelaCE")
  s_data <-
    scaling_m(HelaCE)
  checkEqualsNumeric(sum(as.matrix(s_data)), 7326.892, tolerance = 1e-2)
}
