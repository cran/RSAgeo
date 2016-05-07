\name{gdata}
\docType{data}
\alias{gdata}
\title{Dataset for RSA}
\description{
Simulated date matrix with response variable distributed with 'powered exponential' covariance matrix, 
with \eqn{\phi=25, \tau=1, \kappa=1, \sigma=1}.
}
\usage{data(gdata)}
\format{
  A sample dataset with 1000 rows and 4 columns. Each row represents one observation. The first two column of the matrix consists of x1 and x2, where x1 and x2 provide the 2D location for each obsevation. The third column of the matrix gives the y value, which is the response value. Covariate is included in Column 4. This is the default form for the variable Data in function RSA.
}
\details{

####code used to generate gdata####\cr
require("geoR")\cr
require("RandomFields")\cr
DataNum=1000;\cr
gData=grf(DataNum,grid="irreg",DataNum,DataNum,xlims=c(0,100),ylims=c(0,100),nsim=1,mean=0,
          cov.mode="powered.exponential",cov.par=c(1,25),nugget=1,kappa=1)\cr
l=gData$coords\cr
y=gData$data\cr
x=rnorm(DataNum)\cr
y2=y+.5+x\cr
gdata=cbind(l,y2,x)\cr
}
\examples{
data(gdata)
## This is an illustration of the function. 
## In practice, Total_Iteration is recommended to be set to 2500 or higher.
## N_subset is recommemded to be set to 300 or higher.
RSA(gdata,N_subset=50,Stepscale=40,Total_Iteration=100,Warm=20)

}
\keyword{datasets}
