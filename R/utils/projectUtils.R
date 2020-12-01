#######################################################################
# Script used for saving, building data used all over the project
######################################################################
#' Internal Project use:
#' Function to build h5 hierarchal data using list dataframes
#'
#' @param data list of named dataframes.
#' @param path_to_fname path to desired saved h5 data.
#' @param h5_group grouping inside h5
#' @examples
#' build_list_df_to_h5(list(age = tapA, gender = tapG, education = tapE), "cct_rest.h5", "rest/gender")
build_list_df_to_h5 <- function(data,
                                path_to_fname,
                                h5_group){
  if(is.null(data)){
    return("ERROR: PLEASE PARSE NON-NULL DATA")
  }
  if(!file.exists(path_to_fname)){
    h5createFile(path_to_fname)
  }
  # delete groups when rewriting for consistency and repeatability when testing
  if((h5_group %in% h5ls(path_to_fname)$name) &
     nrow(h5ls(path_to_fname)) != 0){
    h5delete(path_to_fname, h5_group)
  }
  h5createGroup(path_to_fname, h5_group)
  for(cols in names(data)){
    h5write(data[[cols]],
            path_to_fname,
            file.path(h5_group, cols))
  }
}

#' Internal Project use:
#' Function to get the current running file
get_this_file <- function(){
    tryCatch({
      return(basename(rstudioapi::getSourceEditorContext()$path))},
      error = function(e){
        this.file <- commandArgs() %>%
          tibble::enframe(name=NULL) %>%
          tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
          dplyr::filter(key == "--file") %>%
          dplyr::pull(value)
        return(basename(this.file))
      })
}

#' utility function to set features as row names
#' this is used to shape .h5 file inability to store rownames
#' @param data: dataframe containing feature-rows
#' @return index feature as rows
index_features <-function(data){
  data <- data %>% 
    magrittr::set_rownames(.$feature) %>% 
    dplyr::select(-feature)
  return(data)
}

#' utility function to index features on each
#' hierarchial structure of the n of 1 matrix data
#' @param tod_data: list of list of each user and their metrics on the n of 1 analysis
#' @return n of 1 data with indexed features
index_user_metrics <- function(tod_data){
  for(metrics in names(tod_data)){
    tod_data[[metrics]] <- 
      map(tod_data[[metrics]], 
          function(x){
            x <- x %>% index_features(.)
          })
  }
  return(tod_data)
}

