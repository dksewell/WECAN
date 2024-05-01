\name{eClust_GEM_wn}
\alias{eClust_GEM_wn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Edge clustering Algorithm for weighted network with a noisy component.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eClust_GEM_wn(A, K = 20,
              p = 3,
              maxIter=1e3,
              maxIterVB=100,
              a_0 = 1,
              A_0 = 3,
              B_0 = 0.5,
              nu_0 = 2,
              eta_0 = 5,
              QNSteps=25,
              CGSteps=25,
              distribution = "Normal",
              if_log = FALSE,
              eps=1e-5,
              upsilon_SR1 = 2,
              upsilon_SR2 = 2,
              upsilon_UV = 2,
              Phi_SR1 = toeplitz((2:1)/2),
              Phi_SR2 = toeplitz((2:1)/2),
              Phi_UV = toeplitz((2:1)/2),
              noise_mean = 0.1,
              alpha_init = NULL,
              a_a = 1,
              b_a = 200,
              t_0 = 0.05,
              alpha_0 = 1.1,
              beta_0 = 25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
   The adjcent matrix of the network.
 }
  \item{p}{
  Denmision of the latent features.
}
  \item{noise_mean}{
  The mean weight of the "noisy edges".
}
  \item{distribution}{
  Exponenital distribution family for the weights of the edges.
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{estimates$Z_est}{
 Estimated edge cluster.
}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Haomin Li
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
