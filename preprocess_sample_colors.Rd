\name{preprocess_sample_colors}
\alias{preprocess_sample_colors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pre-processing of PSI input data for plotting with colors
}
\description{
\code{preprocess_sample_colors} re-orders PSI sample columns to a specified order
and generates a corresponding color code sequence. The sample order and
corresponding colors are taken from a pre-defined tab-delimited file 
(see details).
}
\usage{
preprocess_sample_colors(psi, database)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psi}{
  a \emph{n} x \emph{m} data frame of PSI values where \emph{n} is the number of AS events
and \emph{m} is the number of samples.
  }
  \item{database}{
  filename of the samples database in tab-delimited format.
  }
}
\details{
\code{preprocess_sample_colors} depends on a pre-defined "samples inventory"
database file in tab-delimited format. This file is species-specific and 
consists of four columns: \emph{Order, SampleName, GroupName, RColorCode}. This
header is required in the file.

An example entry looks like:
\preformatted{Order SampleName  GroupName RColorCode
1 Ooctye  EarlyDev  36}

The \emph{SampleName} must match the column names in \code{psi}.
\cr

The RColorCode value corresponds to the index of the vector produced by
\code{colors()}. For example, RColorCode = 36 corresponds to:
\preformatted{> cols <- colors()
> mycolour <- cols[36]
> mycolour
[1] "brown4"}

}
\value{
A list consisting of:
  \item{data}{data frame of PSI values with columns re-ordered}
  \item{col}{vector of colors corresponding to the re-ordered columns}
  \item{group.index}{list of column indices corresponding to each \emph{GroupName}}
  \item{group.col}{vector of colors corresponding to each \emph{GroupName}}
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
Kevin Ha <k.ha@mail.utoronto.ca>, Blencowe Lab
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

See Plot_PSIs_NEW_KH.R for usage example.

}

