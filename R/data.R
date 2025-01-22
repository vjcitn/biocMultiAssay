#' @name miniACC
#'
#' @title Adrenocortical Carcinoma (ACC) MultiAssayExperiment
#'
#' @docType data
#'
#' @references Zheng S *et al.*: Comprehensive Pan-Genomic Characterization of
#' Adrenocortical Carcinoma. Cancer Cell 2016, 29:723-736.
#'
#' @description
#' A [`MultiAssayExperiment`] object providing a reduced version of
#' the TCGA ACC dataset for all 92 patients. RNA-seq, copy number, and somatic
#' mutations are included only for genes whose proteins are included in the
#' reverse-phase protein array. The MicroRNA-seq dataset is also included,
#' with infrequently expressed microRNA removed. Clinical, pathological, and
#' subtype information are provided by \code{colData(miniACC)}, and some
#' additional details are provided by metadata(miniACC).
#'
#' @keywords data
#'
#' @format A `MultiAssayExperiment` with 5 experiments, providing:
#' \describe{
#'     \item{RNASeq2GeneNorm}{RNA-seq count data: an `ExpressionSet`
#'     with 198 rows and 79 columns}
#'     \item{gistict}{Reccurent copy number lesions identified by GISTIC2:
#'     a `SummarizedExperiment` with 198 rows and 90 columns}
#'     \item{RPPAArray}{Reverse Phase Protein Array: an `ExpressionSet`
#'     with 33 rows and 46 columns. Rows are indexed by genes,
#'     but protein annotations are available from
#'     \code{featureData(miniACC[["RPPAArray"]])}. The source of these
#'     annotations is noted in \code{abstract(miniACC[["RPPAArray"]])}}
#'     \item{Mutations}{Somatic mutations: a `matrix` with 223 rows and
#'     90 columns. 1 for any kind of non-silent mutation, zero for silent
#'     (synonymous) or no mutation.}
#'     \item{miRNASeqGene}{microRNA sequencing: an `ExpressionSet` with
#'     471 rows and 80 columns. Rows not having at least 5 counts in at least
#'     5 samples were removed.}
#' }
#'
#' @author Levi Waldron \email{lwaldron.research@gmail.com}
#'
#' @source https://github.com/waldronlab/multiassayexperiment-tcga
#'
#' @usage data("miniACC")
#'
#' @examples
#'
#' data("miniACC")
#' metadata(miniACC)
#' colnames(colData(miniACC))
#' table(miniACC$vital_status)
#' longForm(
#'     miniACC["MAPK3", , ],
#'     colDataCols = c("vital_status", "days_to_death")
#' )
#'
#' wideFormat(
#'     miniACC["MAPK3", , ],
#'     colDataCols = c("vital_status", "days_to_death")
#' )
#'
"miniACC"
