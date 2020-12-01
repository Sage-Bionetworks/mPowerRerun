#' function get tapping worse hand association
#' @param datTap tap data with healthcode and numTaps
#' @return a summarized version of user based on their worse performing
#' hands in terms of number of tapping
worse.tap.hands.association <- function(datTap){
    tap.feat.list <- list()
    tap.feat.list$left.tap <- datTap %>%
        filter((tappingHands == "left")) %>% 
        dplyr::group_by(healthCode) %>% 
        dplyr::summarise(left.median.num.taps = median(numberTaps, na.rm = TRUE))
    tap.feat.list$right.tap <- datTap %>% 
        filter((tappingHands == "right"))%>% 
        dplyr::group_by(healthCode) %>% 
        dplyr::summarise(right.median.num.taps = median(numberTaps, na.rm = TRUE))
    tap.feat.list$null.tap <- datTap %>% 
        filter(is.na(tappingHands)) %>% 
        dplyr::group_by(healthCode) %>% 
        dplyr::summarise()
    user.worse.hand.association <- tap.feat.list$right.tap %>% 
        full_join(tap.feat.list$left.tap, by = "healthCode") %>% 
        mutate(right.median.num.taps = replace_na(right.median.num.taps, 0),
               left.median.num.taps = replace_na(left.median.num.taps, 0)) %>%
        mutate(inferred.worse.hand = case_when(
            right.median.num.taps == left.median.num.taps ~ "right", ## initial assumptions
            right.median.num.taps > left.median.num.taps ~ "left",
            right.median.num.taps < left.median.num.taps ~ "right")) %>% 
        full_join(tap.feat.list$null.tap %>% dplyr::select(healthCode), by = c("healthCode")) %>% 
        dplyr::select(-c(right.median.num.taps, left.median.num.taps))
    return(user.worse.hand.association)
}