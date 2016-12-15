#' Train gene expression model from methylation profiles
#'
#' \code{train_classification_model} trains a regression model for predicting gene
#' expression levels by taking as input the higher order methylation features
#' extracted from specific genomic regions.
#'
#' @param formula An object of class \code{\link[stats]{formula}}, e.g. see
#'   \code{\link[stats]{lm}} function. If NULL, the simple linear regression
#'   model is used.
#' @param model_name A string denoting the regression model. Currently,
#'   available models are: \code{"svm"}, \code{"randomForest"}, \code{"rlm"},
#'   \code{"mars"} and \code{"lm"}.
#' @param train The training data.
#' @param is_summary Logical, print the summary statistics.
#'
#' @return A list containing the following elements: \itemize{ \item
#'   \code{formula}: The formula that was used. \item \code{gex_model}: The
#'   fitted model. \item \code{train_pred} The predicted values for the training
#'   data. \item \code{train_errors}: The training error metrics. }
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{predict_model_gex}}
#'
#' @importFrom stats formula predict
#'
#' @export
train_classification_model <- function(formula = NULL, model_name = "svm", train,
                                      is_summary = TRUE){
    if (is.null(formula)){
        formula <- y ~ .
    }

    wts <- 1 - (table(train$y) / (sum(table(train$y))))
    if (model_name == "randomForest"){
        model <- randomForest::randomForest(formula = formula,
                                            data = train,
                                            ntree = 1001,
                                            classwt = wts,
                                            nodesize = 3)

        # Make predictions
        train_pred <- predict(object = model,
                              newdata = train[, 1:(NCOL(train) - 1), drop = FALSE],
                              type = "response")
    }else if (model_name == "glm"){
        model <- stats::glm(formula = formula,
                            data = train,
                            family = stats::binomial(link = "logit"))
        # Make predictions
        train_pred <- predict(object = model,
                              newdata = train[, 1:(NCOL(train) - 1), drop = FALSE],
                              type = "response")
    }else{
        model <- e1071::svm(formula = formula,
                            data = train,
                            kernel = "radial",
                            cross = 10,
                            probability = TRUE,
                            cost = 1,
                            class.weights = wts)

        # Make predictions
        train_pred <- predict(object = model,
                              newdata = train[, 1:(NCOL(train) - 1), drop = FALSE],
                              type = "response",
                              probability = TRUE)
    }

    if (length(train_pred) != length(train$y)){
        warning("The classification model returned NAs")
        train_errors <- NULL
    }else{
        # Calculate model errors
        if(is.factor(train_pred)){
            fitted_results <- ifelse(as.numeric(levels(train_pred))[train_pred] > 0.5, 1, 0)
        }else{
            fitted_results <- ifelse(train_pred > 0.5, 1, 0)
        }
        train_errors <- mean(fitted_results != train$y)
        confusion_matrix <- caret::confusionMatrix(as.factor(fitted_results), train$y)

        if (is_summary){
            message("-- Train Errors --")
            print(paste("Accuracy", 1 - train_errors))
        }
    }

    out <- list(formula      = formula,
                gex_model    = model,
                train_pred   = train_pred,
                train_errors = train_errors,
                confusion_matrix = confusion_matrix)
    return(out)
}



#' Predict gene expression model from methylation profiles
#'
#' \code{predict_classification_model} makes predictions of gene expression
#' levels using a model trained on higher order methylation features extracted
#' from specific genomic regions.
#'
#' @param model The fitted regression model, i.e. the output of
#'   \code{\link{train_model_gex}}.
#' @param test The testing data.
#' @param is_summary Logical, print the summary statistics.
#'
#' @return A list containing the following elements: \itemize{ \item
#'   \code{test_pred}: The predicted values for the test data. \item
#'   \code{test_errors}: The test error metrics. }
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{train_model_gex}}
#'
#' @importFrom stats lm predict
#'
#' @export
predict_classification_model <- function(model, test, is_summary = TRUE){
    # Convert to a data.frame
    test <- as.data.frame(test)

    # Make predictions
    test_pred <- predict(object  = model,
                         newdata = test[, 1:(NCOL(test) - 1), drop = FALSE],
                         type = "response",
                         probability = TRUE)

    # Calculate model errors
    if(is.factor(test_pred)){
        fitted_results <- ifelse(as.numeric(levels(test_pred))[test_pred] > 0.5, 1, 0)
    }else{
        fitted_results <- ifelse(test_pred > 0.5, 1, 0)
    }
    test_errors <- mean(fitted_results != test$y)
    confusion_matrix <- caret::confusionMatrix(as.factor(fitted_results), test$y)

    if (is_summary){
        message("-- Test Errors --")
        print(paste("Accuracy", 1 - test_errors))
    }

    out <- list(test_pred   = test_pred,
                test_errors = test_errors,
                confusion_matrix = confusion_matrix)
    return(out)
}
