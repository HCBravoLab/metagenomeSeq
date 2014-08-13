test_cumNorm <- function() {
    data(mouseData)
    checkEquals(cumNorm(mouseData,p=.5), mouseData)
}