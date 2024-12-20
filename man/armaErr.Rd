% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/error.R
\docType{class}
\name{armaErr}
\alias{armaErr}
\title{\code{armaErr} class generator using \code{\link{arima.sim}}}
\description{
\code{armaErr} class generator using \code{\link{arima.sim}}

\code{armaErr} class generator using \code{\link{arima.sim}}
}
\examples{

# set up a stationary arma model
testErr <- armaErr$new(model = list(ar=c(.5), ma=c(.2)))
testErr$checkModel()

# show testing for an unstationary ARMA model
testFail <- armaErr$new(model = list(ar=c(-.8, .5), ma=c(.2)))
testFail$checkModel()

# show that r* functions from gamlss.family distribution functions can be
# passed to arima.sim
errBeta <- arima.sim(list(ar=c(.5)), 1000, innov = rBE(1000, mu=.5, sigma=.2))

# note that the rand.gen option only allows the defaults for the function
# passed to rand.gen
#errBeta <- arima.sim(list(ar=c(.5)), 1000, rand.gen = rBE)
# that said, rand.gen's parameters can be passed via ...
#errBeta <- arima.sim(list(ar=c(.5)), 1000, rand.gen = rBE, mu=.9, sigma=.2)

# note that even though beta innovations are used, the resulting data is not
# constrained to be in (0,1)
hist(errBeta)
plot(errBeta)
auto.arima(ts(errBeta))
}
\author{
Stephen Tueller \email{stueller@rti.org}
}
\section{Super class}{
\code{\link[PersonAlyticsPower:err]{PersonAlyticsPower::err}} -> \code{errARMA}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{model}}{The parameters for \code{FUN}.}

\item{\code{errVar}}{The error variance.}

\item{\code{FUN}}{The function for simulating error terms.}

\item{\code{fam}}{A \code{\link{gamlss.family}} family for the error term distribution.
After data are simulated using \code{FUN} and its \code{model}, it can be
transformed to have a distribution different from those available to \code{FUN}.}

\item{\code{famParms}}{Parameters to be passed to \code{fam}, see \code{\link{gamlss.family}}.}
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{model}}{The parameters for \code{FUN}.}

\item{\code{errVar}}{The error variance.}

\item{\code{FUN}}{The function for simulating error terms.}

\item{\code{fam}}{A \code{\link{gamlss.family}} family for the error term distribution.
After data are simulated using \code{FUN} and its \code{model}, it can be
transformed to have a distribution different from those available to \code{FUN}.}

\item{\code{famParms}}{Parameters to be passed to \code{fam}, see \code{\link{gamlss.family}}.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-errARMA-new}{\code{armaErr$new()}}
\item \href{#method-errARMA-checkModel}{\code{armaErr$checkModel()}}
\item \href{#method-errARMA-makeErrors}{\code{armaErr$makeErrors()}}
\item \href{#method-errARMA-clone}{\code{armaErr$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="PersonAlyticsPower" data-topic="err" data-id="print"><a href='../../PersonAlyticsPower/html/err.html#method-err-print'><code>PersonAlyticsPower::err$print()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-errARMA-new"></a>}}
\if{latex}{\out{\hypertarget{method-errARMA-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{armaErr$new(
  model = list(ar = c(0.5), ma = c(0)),
  fam = "NO",
  famParms = list(mu = 0, sigma = 1)
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-errARMA-checkModel"></a>}}
\if{latex}{\out{\hypertarget{method-errARMA-checkModel}{}}}
\subsection{Method \code{checkModel()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{armaErr$checkModel(
  seed = 1234,
  n = 1000,
  doPlot = TRUE,
  doStats = TRUE,
  debug = FALSE
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-errARMA-makeErrors"></a>}}
\if{latex}{\out{\hypertarget{method-errARMA-makeErrors}{}}}
\subsection{Method \code{makeErrors()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{armaErr$makeErrors(n, nObservations, seed = 123, debug = FALSE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-errARMA-clone"></a>}}
\if{latex}{\out{\hypertarget{method-errARMA-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{armaErr$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
