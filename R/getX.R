getX <- function(model){UseMethod("getX")}

getX.lmerModgANOVA <- function(model){getME(model,"X")}

getX.list <- function(model){model$X}
