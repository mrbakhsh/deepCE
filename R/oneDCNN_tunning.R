


  #' oneDCNN_tunning
  #' @title Tune Hyperparamters for 1D-CNN
  #' @param train_d A matrix of training data containing numerical features,
  #' generated from \code{\link{build_trainingData}}.
  #' @param train_l A vector of binary categorical label (0-1),
  #' generated from  \code{\link{build_trainingData}}.
  #' @param nlayers A vector of integers, describing
  #' number of hidden layers; defaults to 3.
  #' @param filters_1 A vector of integers, describing
  #' the dimensionality of the output space
  #' (i.e. the number of output filters in the first convolution).
  #' @param filters_2 A vector of integers, describing the
  #' dimensionality of the output space
  #' (i.e. the number of output filters in the subsequent convolution).
  #' @param kernel_size A vector of integers, specifying the
  #' height and width of the 2D convolution window. Can be a single integer to
  #' specify the same value for all spatial dimensions.
  #' @param strides A vector of integers, specifying the
  #' strides of the convolution along the height and width. Can be a single
  #' integer to specify the same value for all spatial dimensions.
  #' @param pool_size A vector of integers, describing the pool_size.
  #' @param powerto A vector of integers, describing
  #' the number of neurons in the last layer as defined by two to
  #' the power of this value.
  #' @param drate A vector of numbers, describing the dropout rates.
  #' @param optimizer A vector of strings, specifying the name of optimizer.
  #' For most models, this defaults to "rmsprop".
  #' @param b_size A vector of integers, describing,
  #' number of samples per gradient update.
  #' @param epochs Number of epochs to train the model.
  #' @param k Float between 0 and 1. Fraction of the training data to be used
  #' as validation data.
  #' @param metrics List of metrics to be used by the model during k-fold
  #' cross-validation. Defaults to 'accuracy'.
  #' @return A data.frame containing model performance across different
  #' combination of parameters.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @description This function tune different combinations of hyperparameters
  #' for one-dimensional convolutional neural network (1D-CNN) 1D-CNN model.
  #' @export








  oneDCNN_tunning <-
    function(train_d, train_l,
             nlayers = 2,
             filters_1 = c(32, 64),
             filters_2 = c(32, 64),
             kernel_size = c(1,2),
             strides = 2,
             pool_size = 2,
             powerto = 6,
             drate = 0.1,
             optimizer = "rmsprop",
             b_size = 64,
             metrics = "accuracy",epochs = 10, k = 0.1){

      Var1 <- NA
      Var2 <- NA
      Var3 <- NA
      Var4 <- NA
      Var5 <- NA
      Var6 <- NA
      Var7 <- NA
      Var8 <- NA
      Var9 <- NA
      Var10 <- NA
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
                         filters_1 = filters_1,
                         kernel_size = kernel_size,
                         strides = strides,
                         pool_size = pool_size,
                         drate = drate,
                         filters_2 = filters_2,
                         powerto = powerto,
                         optimizer = optimizer,
                         b_size = b_size) {

        # define the exmaple nework
        network <- keras_model_sequential() %>%
          layer_conv_1d(filters = 32,
                        kernel_size = 2,
                        strides = 2,
                        activation = "relu",
                        input_shape = c(ncol(tdata),1)) %>%
          layer_max_pooling_1d(pool_size = 2) %>%
          layer_batch_normalization() %>%
          layer_dropout(rate = drate)




        # perform tuning
        netwrok <-
          keras_model_sequential() %>%
          layer_conv_1d(filters = filters_1,
                        kernel_size = kernel_size,
                        strides = strides,
                        activation = "relu",
                        input_shape = c(ncol(tdata),1)) %>%
          layer_max_pooling_1d(pool_size = pool_size) %>%
          layer_batch_normalization() %>%
          layer_dropout(rate = drate)


        if(nlayers > 1) {
          map(2:nlayers, ~ network %>%
                layer_conv_1d(filters = filters_2,
                              kernel_size = kernel_size,
                              strides = strides,
                              activation = "relu",
                              input_shape = c(ncol(tdata),1)) %>%
                layer_max_pooling_1d(pool_size = pool_size) %>%
                layer_batch_normalization() %>%
                layer_dropout(rate = drate)
          )

        }


        network %>%
          layer_flatten() %>%
          layer_dense(units = 2^powerto, activation = "relu") %>%
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
      filters_t1 <- filters_1
      kernel_size_t = kernel_size
      strides_t = strides
      pool_size_t = pool_size
      r1_t <- drate
      filters_t2 <- filters_2
      powerto_t <- powerto
      opt_t <- optimizer
      bc_t <- b_size


      param_comb <-
        expand.grid(n_layer,filters_t1,kernel_size_t,strides_t,pool_size_t,
                    r1_t,filters_t2,powerto_t,opt_t,bc_t)

      pmap(list(param_comb$Var1,param_comb$Var2,
                param_comb$Var3, param_comb$Var4,
                param_comb$Var5,param_comb$Var6,
                param_comb$Var7, param_comb$Var8,
                param_comb$Var9, param_comb$Var10),tune_f)  %>%
        map_df(data.frame, .id = "run") -> df_runs


      p <- param_comb
      p[,1] <- paste0("nlayer:", p[,1])
      p[,2] <- paste0("filters_1L:", p[,2])
      p[,3] <- paste0("kernel_size:", p[,3])
      p[,4] <- paste0("strides:", p[,4])
      p[,5] <- paste0("pool_size:", p[,5])
      p[,6] <- paste0("rate:", p[,6])
      p[,7] <- paste0("filters_2L:", p[,7])
      p[,8] <- paste0("powerto:", p[,8])
      p[,9] <- paste0("opt:", p[,9])
      p[,10] <- paste0("b_size:", p[,10])




      # collapse the paramters
      p <- unite(p, param,
                 c(Var1,Var2,Var3,Var4,Var5,Var6,Var7,Var8,Var9,Var10), sep = "~")

      p <- rowid_to_column(p, "run")
      p$run <- as.character(as.numeric(p$run))

      output <- left_join(p, df_runs, by = "run")

      return(output)

    }







