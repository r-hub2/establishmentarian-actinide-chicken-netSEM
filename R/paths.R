#' paths data frame
#'
#' @description The function finds the best paths through netSEM modeling.
#'
#' @details paths_analysis uses a recursive algorithm to rank all different paths through a netSEM model by the lowest R2adj in the pathway (the weakest link). Function returns the best paths in a data frame.
#'
#' @param a a netSEM object containing every variable
#' @param top_n an integer. How many paths from each variable should the function look at. Default value is 3.
#' @param p_val a numeric. The threshold below which a p-value is considered significant, and thus the path is considered valid. The default is 0.05, meaning that paths with p-values greater than 0.05 will not be considered significant. This parameter can be adjusted based on the specifics of your analysis.
#'
#' @return a data frame with two columns: Paths, containing the best paths through the netSEM as a string and min_R2 containing the smallest R2adj value between two variables in the path.
#' @export
#'
#' @examples
#' \dontrun{
#' ## Load the sample acrylic data set
#' data(acrylic)
#'
#' ## Run netSEM
#' ans <- netSEMp1(acrylic, "IrradTot", "YI")
#'
#' path_data <- paths(ans, 3, 0.05)
#' }

paths <- function(a, top_n = 3, p_val = 0.05) {

  find_path <- function(current_var, output_var, cor_mat, path, min_R2adj) {
    path <- c(path, current_var)

    next_vars <- find_next_var(current_var, cor_mat, path)
    next_vars_vec <- as.vector(next_vars$v1)
    x <- 1
    if (current_var == output_var | length(next_vars_vec) == 0) {
      return(list(list(path = path, 
                       min_R2 = min_R2adj)))  # return a list of paths
    } else {
      paths <- list()
      for (var in next_vars_vec) {
        new_min_R2adj <- min(next_vars$R2adj[x], min_R2adj)
        paths <- c(paths, 
                   find_path(var, 
                             output_var, 
                             cor_mat, 
                             path, 
                             new_min_R2adj))  # append lists of paths
        x <- x + 1
      }
      return(paths)
    }
  }

  find_next_var <- function(current_var, cor_mat, path) {
    # Exclude variables already in the path
    diff_vars <- unique(setdiff(cor_mat$var1, path))

    if(length(diff_vars) == 0) {
      return(data.frame(v1 = character(), R2adj = numeric())) # return an empty dataframe if no variables left
    } else {
      # Filter cor_mat for rows where var2 is the current_var
      cor_mat_current <- cor_mat[cor_mat$var2 == current_var, ]

      # Join this with diff_vars
      available_vars <- merge(x = data.frame(v1 = diff_vars), 
                              y = cor_mat_current, 
                              by.x = "v1", 
                              by.y = "var1", all.x = TRUE)
      available_vars <- available_vars %>% 
        dplyr::filter(p_value < p_val)

      # Find the top_n variables most correlated with the current one
      available_vars <- available_vars %>%
        dplyr::arrange(dplyr::desc(R2adj)) %>%
        head(n = top_n)
      return(available_vars)
    }
  }

  output <- a[["bestModels"]]$Resp[1]
  imput <- a[["bestModels"]]$Var[1]

  #create data frame with all the needed netSEM values
  group_data <- data.frame(id = 1:length(a$adj.r.squared), 
                           var1 = NA, 
                           var2 = NA, 
                           R2adj = NA)

  #fill data frame from netSEM
  group_data <- a[["bestModels"]] %>%
    dplyr::rename(
      var2 = Var,
      var1 = Resp,
      R2adj = adj.r.squared,
      model = best_model,
      p_value = p.value
    )

  group_data <- dplyr::select(group_data, var1, var2, R2adj, model, p_value)

  group_data <- group_data %>%
    dplyr::group_by(var1, var2, model) %>%
    dplyr::slice_max(order_by = R2adj, n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(var1, var2) %>%
    dplyr::slice_max(order_by = R2adj, n = 1) %>%
    dplyr::ungroup()

  # Create the lookup table
  lookup_table <-
    setNames(group_data$R2adj,
             paste(group_data$var1, 
                   group_data$var2, 
                   group_data$R2adj))

  # Start with a very high initial minimum R^2
  #paths <- find_path(input_var, output_var, cor_mat, c(), 1)
  paths <- find_path(a[["bestModels"]]$Var[1],
                     output,
                     group_data,
                     c(a[["bestModels"]]$Var[1]), 1)

  # Convert to data frame as before
  paths_df <- do.call(rbind, lapply(paths, function(x) {
    data.frame(path = paste(x$path, collapse = " -> "),
               length = length(x$path) - 2,
               min_R2 = x$min_R2, stringsAsFactors = FALSE)
  }))

  #clean the recursive leftover
  paths_df <- paths_df %>% 
    dplyr::mutate(path = trimws(sub("^[^>]*> ", "", path))) %>%
    dplyr::arrange(length, desc(min_R2))

  return(paths_df)
}