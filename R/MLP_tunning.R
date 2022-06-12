


  #' MLP_tunning
  #' @title Tune Hyperparamters for Multi-layer Perception (MLP)
  #' @param train_d A matrix of training data containing numerical features,
  #' generated from \code{\link{build_trainingData}}.
  #' @param train_l A vector of binary categorical label (0-1),
  #' generated from  \code{\link{build_trainingData}}.
  #' @param nlayers a vector on integers, describing the number
  #' of hidden layers.
  #' @param powerto1 a vector of integers, describing the number of neurons
  #' in the first hidden layer as defined by two to the power of this value.
  #' @param powerto2 a vector of integers, describing the number of neurons
  #' in the subsequent hidden layer as defined by two to the power of
  #' this value.
  #' @param drate  A vector of numbers, describing the dropout rates range.
  #' @param optimizer A vector of strings, describing the name of the
  #' optimizer.
  #' @param b_size A vector of integers, describing
  #' number of samples per gradient update.
  #' @param metrics List of metrics to be used by the model during k-fold
  #' cross-validation. Defaults to 'accuracy'.
  #' See \code{\link[keras]{compile.keras.engine.training.Model}} for
  #' more details.
  #' @param epochs Number of epochs to train the model.
  #' @param k Float between 0 and 1. Fraction of the training data to be used
  #' as validation data.
  #' @return A data.frame containing model performance across different
  #' combination of parameters.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom dplyr left_join
  #' @importFrom keras keras_model_sequential
  #' @importFrom keras layer_dense
  #' @importFrom keras fit
  #' @importFrom keras compile
  #' @importFrom keras layer_batch_normalization
  #' @importFrom keras layer_dropout
  #' @importFrom keras callback_early_stopping
  #' @importFrom keras regularizer_l2
  #' @importFrom tidyr unite
  #' @importFrom tibble rowid_to_column
  #' @importFrom purrr map
  #' @importFrom purrr map_df
  #' @importFrom purrr pmap
  #' @description This function tune different combinations of hyperparameters
  #' for Multi-layer Perception (MLP) model.
  #' @export
  #' @examples
  #' t_data <- build_trainingData(HelaCE, refcpx)
  #' MLP_tunning <-
  #' MLP_tunning(t_data$train_d,
  #' t_data$train_l,
  #' nlayers = 2,
  #' powerto1 = c(4,6),
  #' powerto2 = c(4,6),
  #' drate = 0.3,
  #' optimizer = "rmsprop",
  #' b_size = 50,
  #' metrics = "AUC",
  #' epochs = 5,
  #' k = 0.3)


  MLP_tunning <-
    function(train_d,
             train_l,
             nlayers = c(2,3),
             powerto1 = c(4,6),
             powerto2 = c(4,6), drate = 0.3,
             optimizer = "rmsprop",
             b_size = 50, metrics = "accuracy",epochs = 5, k = 0.3){

      Var1 <- NA
      Var2 <- NA
      Var3 <- NA
      Var4 <- NA
      Var5 <- NA
      Var6 <- NA
      param <- NA

      if(nrow(train_d) != length(train_l)) {
        stop("Size of the class lables must be equal to the number of
             interactions in the training set")

      }

      if(!is.matrix(train_d)) train_d <- as.matrix(train_d)

      if(is.character(train_d) == TRUE){
        stop("Training matrix must include numerical variables")
      }

      if(is.character(train_l) == TRUE){
        stop("Class labels must be a numerical vector")
      }

      if(!is.vector(train_l)) train_l <- as.vector(train_l)



      tdata <- train_d
      ldata <- train_l

      tune_f <- function(nlayers = nlayers,
                         powerto1 = powerto1,
                         powerto2 = powerto2,
                         drate = drate,
                         optimizer = optimizer,
                         b_size = b_size) {

        # define the exmaple nework
        network <- keras_model_sequential() %>%
          layer_dense(units = 64, activation = "relu",
                      input_shape = c(ncol(tdata)),
                      kernel_regularizer = regularizer_l2(0.001)) %>%
          layer_batch_normalization() %>%
          layer_dropout(rate = 0.3) %>%
          layer_dense(units = 1, activation = "sigmoid")


        # perform tuning
        netwrok <- keras_model_sequential() %>%
          layer_dense(units = 2^powerto1, activation = "relu",
                      input_shape = c(ncol(tdata)),
                      kernel_regularizer = regularizer_l2(0.001)) %>%
          layer_batch_normalization() %>%
          layer_dropout(rate = drate)


        if(nlayers > 1) {
          map(2:nlayers, ~ network %>%
                layer_dense(units = 2^powerto2, activation = "relu",
                            kernel_regularizer = regularizer_l2(0.001)) %>%
                layer_batch_normalization() %>%
                layer_dropout(rate = drate)
          )

        }

        network %>%
          layer_dense(units = 1, activation = "sigmoid")

        network %>% compile(
          optimizer = optimizer,
          loss = "binary_crossentropy",
          metrics = c(metrics))

        history <- network %>% fit(
          tdata,
          ldata,
          epochs = epochs,
          batch_size = b_size,
          validation_split = k,
          callbacks = list(callback_early_stopping(monitor = "loss",
                                                   patience = 5,
                                                   restore_best_weights = TRUE))

        )
      }

      n_layer <-  nlayers
      powerto_1 <- powerto1
      powerto_2 <- powerto2
      r1 <- drate
      opt <- optimizer
      bc <- b_size

      param_comb <-
        expand.grid(n_layer,powerto_1,powerto_2,r1,opt,bc)

      pmap(list(param_comb$Var1,param_comb$Var2, param_comb$Var3,
                param_comb$Var4, param_comb$Var5,
                param_comb$Var6),tune_f)  %>%
        map_df(data.frame, .id = "run") -> df_runs


      p <- param_comb
      p[,1] <- paste0("nlayer:", p[,1])
      p[,2] <- paste0("n_layer1:", p[,2])
      p[,3] <- paste0("n_layers:", p[,3])
      p[,4] <- paste0("rate:", p[,4])
      p[,5] <- paste0("opt:", p[,5])
      p[,6] <- paste0("b_size:", p[,6])




      # collapse the paramters
      p <- unite(p, param,
                 c(Var1,Var2,Var3,Var4,Var5,Var6), sep = "~")

      p <- rowid_to_column(p, "run")
      p$run <- as.character(as.numeric(p$run))

      output <- left_join(p, df_runs, by = "run")
      return(output)


    }







