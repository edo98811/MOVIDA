FEANames_wrapper <- function(...) {
  tryCatch(
    DeeDeeExperiment::FEANames(...),
    error = function(e) {
      warning("getFEA_list_wrapper error: ", conditionMessage(e))
      NULL
    }
  )
}

DEANames_wrapper <- function(...) {
  tryCatch(
    DeeDeeExperiment::DEANames(...),
    error = function(e) {
      warning("getDEA_list_wrapper error: ", conditionMessage(e))
      NULL
    }
  )
}

FEA_wrapper <- function(...) {
  tryCatch(
    DeeDeeExperiment::FEA(...),
    error = function(e) {
      warning("FEA_wrapper error: ", conditionMessage(e))
      NULL
    }
  )
}

DEA_wrapper <- function(...) {
  tryCatch(
    DeeDeeExperiment::DEA(...),
    error = function(e) {
      warning("DEA_wrapper error: ", conditionMessage(e))
      NULL
    }
  )
}
