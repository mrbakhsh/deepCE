test_similarity_score= function() {
  data('HelaCE')
  scored_Data <- similarity_score(HelaCE, metric = "pearson")
  checkEqualsNumeric(sum(as.matrix(scored_Data)), 15676.69, tolerance = 1e-2)
}
