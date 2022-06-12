test_filter_ConsecutivePep = function() {
  data("HelaCE")
  filt_dat <- filter_ConsecutivePep(HelaCE,10)
  checkEqualsNumeric(sum(as.matrix(filt_dat)), 236920, tolerance = 1e-2)
}
