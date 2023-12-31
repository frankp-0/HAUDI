% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_tracts.R
\name{make_Tracts}
\alias{make_Tracts}
\title{Convert FLARE output to RFMix output}
\usage{
make_Tracts(anc_vcf_file, chunk_size, out_file)
}
\arguments{
\item{anc_vcf_file}{a file path to a VCF file with phased genotypes and estimated local ancestry, i.e. the output of flare (https://doi.org/10.1101/2022.08.02.502540)}

\item{chunk_size}{an integer indicating the max number of VCF records to read at a time}

\item{out_file}{a file path for storing output}
}
\description{
Convert FLARE output to RFMix output
}
\author{
Frank Ockerman
}
