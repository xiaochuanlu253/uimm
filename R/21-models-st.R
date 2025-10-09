methods::setOldClass("rpart")
#' @exportClass STModel
setClass("STModel",
         slots = list(IMM = "MortalityModel", model = "rpart", tnodes = "character")
)

# helper
.get_tree <- function(formula, data) LTRCtrees::LTRCART(formula, data = data)

#' @export
STModel <- function(formula, data = NULL, data_train = NULL, data_test = NULL, IMM = NULL) {
  if (is.null(IMM)) IMM <- MortalityModel(data = data, data_train = data_train, data_test = data_test)
  tree <- .get_tree(formula, IMM@data_train)
  train_copy <- IMM@data_train
  train_copy$node <- as.vector(tree$where)
  tnodes <- levels(factor(train_copy$node))
  test_copy <- IMM@data_test
  if (!is.null(test_copy)) {
    test_copy$group <- predict(partykit::as.party(tree), newdata = test_copy, type = "node")
  }
  IMM@data <- list(train = train_copy, test = test_copy)
  IMM@data_train <- train_copy; IMM@data_test <- test_copy
  methods::new("STModel", IMM = IMM, model = tree, tnodes = tnodes)
}

#' @export
setMethod("show", "STModel", function(object) {
  cat("Training data size:", nrow(object@IMM@data_train), "\n")
  cat("Test data size:", if (is.null(object@IMM@data_test)) 0L else nrow(object@IMM@data_test), "\n")
  cat("Terminal nodes:", paste(object@tnodes, collapse = ", "), "\n")
})

#' @export
setMethod("predict", "STModel",
          function(object, newdata = NULL,
                   type = c("survival", "node"),
                   age_range = 80:115) {
            type <- match.arg(type)
            if (is.null(newdata)) newdata <- object@IMM@data_test
            if (type == "survival") {

            }
            if (type == "node")
              return(predict(object@model %>% partykit::as.party(), newdata = newdata, type = "node"))
          })

#' @export
setMethod("evaluate", "STModel",
          function(object, newdata = NULL, ...) {
            if (is.null(newdata)) newdata <- object@IMM@data_test
            preds <- predict(object, newdata, type = "node")
            list(nodes = object@tnodes,
                 n_by_node = as.integer(table(factor(preds, levels = object@tnodes))))
          })
