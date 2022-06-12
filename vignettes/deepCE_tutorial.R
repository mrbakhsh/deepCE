## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----warning = FALSE, message=FALSE-------------------------------------------
# Loading packages required for data handling & visualization
library(ggplot2)
library(tidyr)
library(dplyr)
library(igraph)
library(fgsea)

# Loading deepCE package
library(deepCE)
# Loading the demo data
data(HelaCE)
dim(HelaCE)

## -----------------------------------------------------------------------------
dat_p1 <- impute_MissingData(HelaCE,bound_left = 2,bound_right = 2)

## -----------------------------------------------------------------------------
dat_p2 <-
  RemoveSinglePeak(dat_p1)

## -----------------------------------------------------------------------------
dat_p3 <-
  filter_ConsecutivePep(dat_p2,min_stretch_length=2)

## -----------------------------------------------------------------------------
dat_p4 <-
  scaling_m(dat_p3)

## ----message = FALSE, warning = FALSE-----------------------------------------
conc_dat <- getPPI(dat_p4)

## -----------------------------------------------------------------------------
# Load reference set 
data("refcpx")
t_data <- build_trainingData(dat_p4, refcpx)

## -----------------------------------------------------------------------------
dat_p2 <-
  RemoveSinglePeak(dat_p1)

## -----------------------------------------------------------------------------
dat_p3 <-
  filter_ConsecutivePep(dat_p2,min_stretch_length=2)

## -----------------------------------------------------------------------------
dat_p4 <-
  scaling_m(dat_p3)

## ----message = FALSE, warning = FALSE-----------------------------------------
conc_dat <- getPPI(dat_p4)

## -----------------------------------------------------------------------------
# Load reference set 
data("refcpx")
t_data <- build_trainingData(dat_p4, refcpx)

## ----message = FALSE, warning = FALSE-----------------------------------------
# set the seed to ensure reproducible output
set.seed(100)
MLP_interactions <- 
  MLP_m(conc_dat, #concatenated co-elution profiles
      t_data$train_d, #training data matrix
      t_data$train_l, #training data label
      cv_fold = 5) 
head(MLP_interactions)

## ----message = FALSE, warning = FALSE, eval = FALSE---------------------------
#  # set the seed to ensure reproducible output
#  set.seed(101)
#  oneDCNN_interactions <-
#    oneD_CNN(conc_dat, #concatenated co-elution profiles
#        t_data$train_d, #training data matrix
#        t_data$train_l, #training data label
#        cv_fold = 5)

## ----message = FALSE, warning = FALSE-----------------------------------------
# set the seed to ensure reproducible output
set.seed(100)
MLP_interactions <- 
  MLP_m(conc_dat, #concatenated co-elution profiles
      t_data$train_d, #training data matrix
      t_data$train_l, #training data label
      cv_fold = 5) 
head(MLP_interactions)

## ----message = FALSE, warning = FALSE, eval = FALSE---------------------------
#  # set the seed to ensure reproducible output
#  set.seed(101)
#  oneDCNN_interactions <-
#    oneD_CNN(conc_dat, #concatenated co-elution profiles
#        t_data$train_d, #training data matrix
#        t_data$train_l, #training data label
#        cv_fold = 5)

## ----results="hide", warning = FALSE, fig.show="hide"-------------------------
set.seed(102)
tuning_result <-
    MLP_tunning(t_data$train_d, 
                t_data$train_l,
                nlayers = 2,
                powerto1 = c(4,6,7),
                powerto2 = c(4,6,7), 
                b_size = 128, 
                metrics = "accuracy",
                epochs = 50, k = 0.3)

## -----------------------------------------------------------------------------
f<- 
  tuning_result %>%
  filter(metric == "loss") 

min_val_error <- 
  f %>%
  filter(data == "validation") %>%
  group_by(run) %>% summarise(min_value.error = 
                                round(min(value, na.rm = TRUE), 3))
head(min_val_error)

## ----results="hide", warning = FALSE------------------------------------------
ggplot(f, aes(epoch, value, color = data)) +
  geom_point(size = 0.1) +
  geom_line() +
  theme_bw() +
  facet_wrap(.~ run)+
  geom_label(data = min_val_error, aes(label=min_value.error), 
             x = Inf, y = -Inf, hjust=1, vjust=0,
             inherit.aes = FALSE)


## -----------------------------------------------------------------------------
pred_cpx <- get_clusters(csize = 3, 
                         d = 0.3, p = 2,mx_overlap = 0.8,
                         tpath =file.path(system.file("extdata", 
                                                      package = "deepCE")))

## -----------------------------------------------------------------------------
data("refcpx")
set.seed(103)
Clust_tuning_result <-
  Clust_tuning(refcpx, csize = 3, 
                d = c(0.3,0.4),
                p = c(2, 2.5),
                mx_overlap = c(0.6,0.7),
                tpath =
                  file.path(system.file("extdata", package = "deepCE")))

## -----------------------------------------------------------------------------
ggplot(Clust_tuning_result, aes(tune_names, compScore)) +
  geom_point(size = 3) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size=15)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(
    color = "black", # label color and size
    size = 12))  +
  theme(axis.ticks = element_line(
    colour = "black",
    size = 0.75, linetype = "solid"
  )) +
  theme(axis.ticks.length=unit(0.3, "cm"))


## ----results="hide", message = FALSE, warning = FALSE-------------------------
data("HelaCE")
data("refcpx")
MLP_prediction_outputs <- 
  predPPI_MLP(HelaCE,
              refcpx,
              cv_fold = 5,
              tpath = tempdir())

## ----message = FALSE, warning = FALSE, eval=FALSE-----------------------------
#  data("HelaCE")
#  data("refcpx")
#  oneDCNN_prediction_outputs <-
#    predPPI_1D.CNN(HelaCE,
#                refcpx,
#                cv_fold = 5)

## ---- eval=FALSE--------------------------------------------------------------
#  ig <-
#     graph_from_data_frame(MLP_prediction_outputs$ppi_input_ClusterOne)
#  
#  createNetworkFromIgraph(ig,"myIgraph")

## ---- warning=FALSE, message=FALSE, results='hide', eval = FALSE--------------
#  # extract the predicted complexes
#  predcpx <- MLP_prediction_outputs$predicted_cpx
#  
#  enrich_result <-
#    enrichfindCPX(predcpx,
#                  threshold = 0.05,
#                  sources = "GO:BP",
#                  p.corrction.method = "bonferroni",
#                  org = "hsapiens")

## ----message=FALSE, results='hide', warning=FALSE, fig.show="hide"------------
# profile 1
data("m1")
# profile 2
data("m2")
# biological term (here is known complexes)
data("refcpx")
diff_output <- diffPPI(m1, 
                       m2, 
                       refcpx,
                       minSize = 2,
                       maxSize = 10000,
                       tpath = tempdir())

## -----------------------------------------------------------------------------
plotEnrichment(diff_output$term_list[["F1F0-ATP synthase, mitochondrial"]],
               diff_output$vec_rank) +
  labs(title="F1F0-ATP synthase, mitochondrial") +
   theme_bw() +
  theme(text = element_text(size=15)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(
    color = "black", # label color and size
    size = 12))  +
  theme(axis.ticks = element_line(
    colour = "black",
    size = 0.75, linetype = "solid"
  )) +
  theme(axis.ticks.length=unit(0.3, "cm"))


## -----------------------------------------------------------------------------
scored_Data <- similarity_score(HelaCE, metric = "pearson")

## -----------------------------------------------------------------------------
sessionInfo()

