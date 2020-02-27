compareRhythms_rain <- function(y, expDesign, Tau=24){

 expDesign <- cbind(expDesign, colNumber = seq(nrow(expDesign)))
 group_IDs <- base::unique(expDesign$group)
 expDesign_A <- base::subset(expDesign, group == group_IDs[1])

 expDesign_A <- expDesign_A[base::order(expDesign_A$time),]

 expDesign_B <- base::subset(expDesign, group == group_IDs[2])

 expDesign_B <- expDesign_B[base::order(expDesign_B$time),]

 deltat_A <- min(diff(unique(expDesign_A$time)))


 time_A <- base::seq(min(expDesign_A$time),
             max(expDesign_A$time),
             by = deltat_A)

 measure.sequence_A <- base::table(expDesign_A$time)
 measure.sequence_A <- base::sapply(time_A, function(t) ifelse(any(names(measure.sequence_A) == t),
                                                         measure.sequence_A[names(measure.sequence_A) == t],
                                                         0))
 deltat_B <- min(diff(unique(expDesign_B$time)))


 time_B <- base::seq(min(expDesign_B$time),
               max(expDesign_B$time),
               by = deltat_B)

 measure.sequence_B <- table(expDesign_B$time)
 measure.sequence_B <- sapply(time_B, function(t) ifelse(any(names(measure.sequence_B) == t),
                                                         measure.sequence_B[names(measure.sequence_B) == t],
                                                         0))
 return(expDesign_B)
}
