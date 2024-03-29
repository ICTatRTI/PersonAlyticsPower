% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/error.R
\docType{class}
\name{Err}
\alias{Err}
\title{\code{Err} class generator}
\description{
\code{Err} class generator

\code{Err} class generator
}
\author{
Stephen Tueller \email{stueller@rti.org}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{model}}{The parameters for \code{FUN}.}

\item{\code{errVar}}{The error variance.}

\item{\code{FUN}}{The function for simulating error terms.}

\item{\code{fam}}{A random variable generator from the
\code{\link{gamlss.family}} for the error term distribution.
After data are simulated using \code{FUN} and its \code{model}, it can be
transformed to have a distribution implied by \code{fam}.}

\item{\code{famParms}}{Parameters to be passed to \code{fam}, see \code{\link{gamlss.family}}.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Err$new()}}
\item \href{#method-print}{\code{Err$print()}}
\item \href{#method-clone}{\code{Err$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Err$new(model = NULL, errVar = 1, FUN = NULL, fam = NULL, famParms = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Err$print(...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Err$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
