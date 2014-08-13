test_cumNormStatFast <- function() {
    data(lungData)
    checkEquals(cumNormStatFast(lungData),0.7014946)
}