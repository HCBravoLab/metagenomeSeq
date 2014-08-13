test_cumNormStat <- function() {
    data(lungData)
    checkEquals(cumNormStat(lungData),0.7014946)
 }