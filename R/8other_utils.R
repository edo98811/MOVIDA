fea_names_wrapper <- function(...) {
  tryCatch(
    DeeDeeExperiment::fea_names(...),
    error = function(e) {
      warning("get_fea_list_wrapper error: ", conditionMessage(e))
      NULL
    }
  )
}

dea_names_wrapper <- function(...) {
  tryCatch(
    DeeDeeExperiment::dea_names(...),
    error = function(e) {
      warning("get_dea_list_wrapper error: ", conditionMessage(e))
      NULL
    }
  )
}

fea_wrapper <- function(...) {
  tryCatch(
    DeeDeeExperiment::fea(...),
    error = function(e) {
      warning("fea_wrapper error: ", conditionMessage(e))
      NULL
    }
  )
}

dea_wrapper <- function(...) {
  tryCatch(
    DeeDeeExperiment::dea(...),
    error = function(e) {
      warning("dea_wrapper error: ", conditionMessage(e))
      NULL
    }
  )
}
