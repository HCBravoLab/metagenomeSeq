test_calcNormFactors <- function() {
    data(lungData)
    point25 = c(29,2475,2198,836,722,1820,79,1171,1985,710,145,742,848,89,1981)
    point = c(43,2475,2198,836,722,1820,119,1171,1985,710,145,742,848,89,1981)
    point100=as.numeric(unlist(libSize(lungData[,1:15])))
    checkEquals(as.numeric(unlist(calcNormFactors(lungData[,1:15]))),point)
    checkEquals(as.numeric(unlist(calcNormFactors(lungData[,1:15],p=.25))),point25)
    checkEquals(as.numeric(unlist(calcNormFactors(lungData[,1:15],p=1))),point100)
 }