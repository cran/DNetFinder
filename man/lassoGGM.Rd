\name{lassoGGM}
\alias{lassoGGM}
\title{
Estimating the regression coefficients in GGMs with lasso 
}
\description{
The function "lassoGGM" computes the lasso estimates of the regression coefficents in GGMs for constructing the test statistic. 
}
\usage{
lassoGGM(Data_mat)
}

\arguments{
  \item{Data_mat}{
A n by p data matrix, where each row represents one observation 
}
}
\details{
The tuning parameter in the lasso regression is chosen as in Liu (2017).
}
\value{
The estimated coefficient matrix by lasso
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
Other estimators such as Dantzig selector or square-root lasso can also be used. See detailed discussion in Liu (2017) and Zhang (2017).
}

\seealso{
lassoNPN()
}
\examples{
Data1=read.table(system.file("extdata","Data1.txt",package="DNetFinder"),header=FALSE)
est_coefGGM=lassoGGM(Data1)
}
\keyword{Gaussian graphical model}
\keyword{Lasso estimate}
