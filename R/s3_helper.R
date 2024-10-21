get_flatnames3d <- function(name, dim = NULL, first_index = 1){

  if(is.null(dim)){
    return(name)
  }

  indices <- lapply(dim, \(ind){first_index:(ind + first_index - 1)})

  indices <- tidyr::expand_grid(!!!indices) |>
    rev() # for column-major order

  indices <- indices |>
    tidyr::unite("ind", sep = ", ")

  names_vec <- paste0(name, "[",  indices$ind, "]")

  return(names_vec)

}
