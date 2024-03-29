% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dotplot.R
\name{seqdotplot}
\alias{seqdotplot}
\title{Draw sequence dotplot}
\usage{
seqdotplot(
  s1,
  s2,
  width = 10,
  step = 1,
  n_mismatches = 0,
  geom = c("auto", "point", "tile"),
  threads = getOption("mc.cores", 1L)
)
}
\arguments{
\item{s1}{the first sequence.}

\item{s2}{the second sequence.}

\item{width}{the width of the window for comparing matches.}

\item{step}{the step size to move the window forward.}

\item{n_mismatches}{the allowed number of mismatches within window.}

\item{geom}{the geom used to draw the dots.}

\item{threads}{the number of threads to use.}
}
\value{
ggplot object containing the sequence dotplot
}
\description{
Draws a sequence dotplot with forward and reverse complement matching. The
longer sequence will be represented on the X-axis.
}
