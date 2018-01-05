\name{DNetNPN}
\alias{DNetNPN}
\title{
Testing for the structural difference between two NPNGMs
}
\description{
The function "DNetNPN" tests for the structural difference between two nonparanormal graphical models with false discovery rate control.
}
\usage{
DNetNPN(Data_mat1,Data_mat2,Beta_mat1,Beta_mat2,alpha)
}
\arguments{
  \item{Data_mat1}{
An n1 by p data matrix for the first NPNGM, where each row represents one observation 
}
  \item{Data_mat2}{
An n2 by p data matrix for the second NPNGM, where each row represents one observation 
}
  \item{Beta_mat1}{
A p-1 by p coefficient matrix for the first NPNGM, where each column contains the regression coefficients of one variable on the other p-1 variables.
}
  \item{Beta_mat2}{
A p-1 by p coefficient matrix for the second NPNGM. See \code{Beta_mat1} for details.
}
 \item{alpha}{
User-specified FDR level
}
}
\details{
The multiple testing procedure asymptotically controls the false discovery rate. See Zhang (2017) for details.
}
\value{
Estimated differential network, where "1" represents a differential edge and "0" represents a common edge (or no edge) between two NPNGMs.
}
\references{
Li, X., Zhao, T., Yuan, X., Liu, H. (2015). The flare Package for High Dimensional Linear Regression and Precision Matrix Estimation in R. Journal of Machine Learning Research, 16:553-557

Liu, H., Lafferty, J., Wasserman, L. (2009). The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. Journal of Machine Learning Research, 10:2295-2328

Liu, W. (2017). Structural Similarity and Difference Testing on Multiple Sparse Gaussian Graphical Models. Annals of Statistics, 45(6):2680-2707

Tibshirani, R. (1996). Regression Shrinkage and Selection via the Lasso. Journal of the Royal Statistical Society Series B, 58(1):267-288

Zhang, Q. (2017). Structural Difference Testing on Multiple Nonparanormal Graphical Models with False Discovery Rate Control. Preprint.
}
\author{
Qingyang Zhang
}
\note{
Besides lasso, other estimators such as Dantzig selector or square-root lasso can also be used. See detailed discussion in Liu (2017) and Zhang (2017).
}

\seealso{
DNetGGM()
}
\examples{
Data1=read.table(system.file("extdata","Data1.txt",package="DNetFinder"),header=FALSE)
Data2=read.table(system.file("extdata","Data2.txt",package="DNetFinder"),header=FALSE)
BetaNPN1=read.table(system.file("extdata","BetaNPN1.txt",package="DNetFinder"),header=FALSE)
BetaNPN2=read.table(system.file("extdata","BetaNPN2.txt",package="DNetFinder"),header=FALSE)
est_DNNPN=DNetNPN(Data1,Data2,BetaNPN1,BetaNPN2,alpha=0.1)
}
\keyword{Nonparanormal graphical model}
\keyword{False discovery rate}
\keyword{Differential network}
