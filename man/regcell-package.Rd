\name{regcell-package}
\alias{regcell-package}
\alias{regcell}
\docType{package}
\title{
\packageTitle{regcell}
}
\description{
\packageDescription{regcell}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{regcell}
\packageIndices{regcell}
robust sparse regression under cellwise contamination
}
\author{
\packageAuthor{regcell}

Maintainer: \packageMaintainer{regcell}
}
\references{
will be on arxiv soon
}
cellwise contamination; sparse regression; robust variable selection
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
data = genevar()
y = data$y
x = data$x
fit = sregcell_std(y,x)

}
