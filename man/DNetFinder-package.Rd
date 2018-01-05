\name{DNetFinder-package}
\alias{DNetFinder-package}
\alias{DNetFinder}
\docType{package}
\title{
\packageTitle{DNetFinder}
}
\description{
\packageDescription{DNetFinder}
}
\details{
\packageIndices{DNetFinder}
}
\author{
\packageAuthor{DNetFinder}

Maintainer: \packageMaintainer{DNetFinder}
}
\references{
Li, X., Zhao, T., Yuan, X., Liu, H. (2015). The flare Package for High Dimensional Linear Regression and Precision Matrix Estimation in R. Journal of Machine Learning Research, 16:553-557

Liu, H., Lafferty, J., Wasserman, L. (2009). The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. Journal of Machine Learning Research, 10:2295-2328

Liu, W. (2017). Structural Similarity and Difference Testing on Multiple Sparse Gaussian Graphical Models. Annals of Statistics, 45(6):2680-2707

Tibshirani, R. (1996). Regression Shrinkage and Selection via the Lasso. Journal of the Royal Statistical Society Series B, 58(1):267-288

Zhang, Q. (2017). Structural Difference Testing on Multiple Nonparanormal Graphical Models with False Discovery Rate Control. Preprint.
}


\keyword{Differential network}
\keyword{Nonparanormal graphical model}
\keyword{False discovery rate control}
\seealso{
lassoGGM(), lassoNPN(), DNetGGM(), DNetNPN()
}
\examples{
library(flare)
library(DNetFinder)
Data1=read.table(system.file("extdata","Data1.txt",package="DNetFinder"),header=FALSE)
Data2=read.table(system.file("extdata","Data2.txt",package="DNetFinder"),header=FALSE)
BetaGGM1=read.table(system.file("extdata","BetaGGM1.txt",package="DNetFinder"),header=FALSE)
BetaGGM2=read.table(system.file("extdata","BetaGGM2.txt",package="DNetFinder"),header=FALSE)
BetaNPN1=read.table(system.file("extdata","BetaNPN1.txt",package="DNetFinder"),header=FALSE)
BetaNPN2=read.table(system.file("extdata","BetaNPN2.txt",package="DNetFinder"),header=FALSE)
est_coefGGM=lassoGGM(Data1)
est_coefNPN=lassoNPN(Data1)
est_DNGGM=DNetGGM(Data1,Data2,BetaGGM1,BetaGGM2,alpha=0.1)
est_DNNPN=DNetNPN(Data1,Data2,BetaNPN1,BetaNPN2,alpha=0.1)
}
