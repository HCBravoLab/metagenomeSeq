rm ./metagenomeSeq.1/.DS_Store
rm ./metagenomeSeq.1/data/.DS_Store
rm ./metagenomeSeq.1/inst/.DS_Store
rm ./metagenomeSeq.1/vignettes/.DS_Store
rm ./metagenomeSeq.1/R/.DS_Store
rm ./metagenomeSeq.1/man/.DS_Store
rm ./metagenomeSeq.1/vignettes/.Rhistory
rm ./metagenomeSeq.1/vignettes/Rplots.pdf
rm ./metagenomeSeq.1/man/.Rhistory
rm ./metagenomeSeq.1/man/.Rapp.history 

 
R CMD build metagenomeSeq.1
#R CMD check metaR_0.99.tar.gz

#R CMD build --resave-data metaR.1/
#R CMD check metagenomeSeq_0.99.0.tar.gz

