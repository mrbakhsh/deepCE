test_build_trainingData = function() {
  data("refcpx")
  data("HelaCE")
  tdata <- build_trainingData(HelaCE, refcpx)
  checkEquals(ncol(tdata$train_d), 624L)
  checkEquals(length(tdata$train_l), 2758L)

}
