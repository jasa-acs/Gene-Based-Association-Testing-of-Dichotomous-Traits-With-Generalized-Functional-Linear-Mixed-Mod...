# Copyright 2020, Georgetown University and University of Pittsburgh.  All Rights Reserved.
#
# This code is dervied from the 'mega2pedgene.R' file from the
#   Mega2R package.
#   Contributors to Mega2R: Robert V. Baron and Daniel E. Weeks.
#
#' load Mega2 SQLite database and perform initialization for PedGFLMM usage
#'
#' @description
#'  This populates the \bold{R} data frames from the specified \bold{Mega2} SQLite database.
#'
#' @param db specifies the path of a \bold{Mega2} SQLite database containing study data.
#'
#' @param verbose TRUE indicates that diagnostic printouts should be enabled.
#'  This value is saved in the returned environment.
#' @param traitname Name of the affection status trait to use to set the case/control values; by default, "default"
#'
#' @return "environment" containing data frames from an SQLite database and some computed values.
#'
#' @import Mega2R
#' @importFrom Mega2R applyFnToRanges
#' @importFrom utils read.table write.table
#' @export
#'
#' @note
#'  \emph{init_PedGFLMM} sets up the schaidPed and pedPer data frames that are used later in the \emph{DOPedGFLMM} calculation.
#'  In addition, it initializes a matrix to aid
#'   in translating a genotype allele matrix to a genotype count matrix.
#'
#'  It also initializes the results data frame \emph{envir$PedGFLMM_results} to zero rows.
#'
#' @seealso \code{\link{DOPedGFLMM}}, \code{\link{Mega2PedGFLMM}}
#'
#' @examples
#' db = system.file("exdata", "seqsimmGFLMM.db", package="PedGFLMM")
#' ENV = init_PedGFLMM(db, traitname = "default")
#' ls(ENV)
#'
init_PedGFLMM = function (db = NULL, verbose = FALSE, traitname = "default") {

    if (is.null(db))
        stop("You must specify a database argument!\n", call. = FALSE)

    envir = dbmega2_import(db, verbose = verbose)

    fam = mkfam(envir = envir, traitname = traitname)
#   fam = fam[fam$trait != 0, ]  # b,c vs a
#   vs
    fam$trantrait = NA
    fam$trantrait[fam$trait == 1] = 0
    fam$trantrait[fam$trait == 2] = 1
    fam$trait = NULL

    setfam(fam, envir = envir)  # also updates unified_genotype_table
    envir$traitname <- traitname
    envir$schaidPed = envir$fam[ , c(-1, -2)]
    colnames(envir$schaidPed) = c("ped", "person", "father", "mother", "sex", "trait")
    envir$pedPer = envir$schaidPed[ , 1:2]
#   envir$mt = matrix(c(11, 12, 21, 22, 0,    0, 1, 1, 2, 0), nrow = 5, ncol = 2)
    envir$mt1 = c(0x10001, 0x10002, 0x20001, 0x20002, 0)
    envir$mt2 = c(      0,       1,       1,       2, 0)
    envir$PedGFLMM_results <- data.frame(chr = character(0), gene = character(0),
                                        nvariants = numeric(0),
                                        start = numeric(0), end = numeric(0),
                                        LRT.bsm.Bsp = numeric(0),
                                        LRT.bsm.fsp = numeric(0),
                                        stringsAsFactors = FALSE)
    return (envir)
}


#' Execute the PedGFLMM_beta_smooth_only function on a transcript ranges
#'
#' @description
#'  This example function illustates how to use functions from the \code{Mega2R} R package to
#'  iterate over defined gene ranges, computing the \code{PedGFLMM_beta_smooth_only} statistics
#'  for each gene that contains more than two polymorphic markers.
#'
#' Execute the PedGFLMM_beta_smooth_only function on the first \emph{gs} default gene transcript ranges (gs = 1:100).
#'  Update the \emph{envir$PedGFLMM_results} data frame with the results.
#"
#' @param gs a subrange of the default transcript ranges over which to calculate the \emph{DOPedGFLMM} function.
#'
#' @param genes a list of genes over which to calculate the \emph{DOPedGFLMM} function.
#'  The value, "*", means use all the transcripts in the selected Bioconductor database.
#'  If genes is NULL, the gs range of the internal \emph{refRanges} will be used.
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return None
#'  the data frame with the PedGFLMM_beta_smooth_only results is stored in the environment and named \emph{PedGFLMM_results},
#'  viz. envir$PedGFLMM_results
#'
#' @import Mega2R
#' @importFrom Mega2R applyFnToRanges
#' @export
#'
#' @seealso \code{\link{init_PedGFLMM}}
#'
#' @examples
#' db = system.file("exdata", "seqsimmGFLMM.db", package="PedGFLMM")
#' ENV = init_PedGFLMM(db)
#' ENV$verbose = TRUE
#' Mega2PedGFLMM(gs = 50:60)
#'
Mega2PedGFLMM = function (gs = 1:100, genes = NULL, envir = ENV) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

    if (is.null(genes))
        Mega2R::applyFnToRanges(DOPedGFLMM, envir$refRanges[gs, ], envir$refIndices, envir = envir)
    else
        Mega2R::applyFnToGenes(DOPedGFLMM, genes, envir = envir)
}

#' PedGFLMM_beta_smooth_only call back function
#'
#' @description
#'  First, ignore call backs that have less than two polymorphic markers.  Second, convert the genotypesraw()
#'  patterns of 0x10001, 0x10002 (or 0x20001), 0x20002, 0 from the genotype matrix
#'  to the numbers 0, 1, 2, 0 for each marker. (Reverse, the order if allele "1" has the
#'  minor allele frequency.)  Next, prepend the pedigree and person columns of the family data
#'  to this modified genotype matrix.  Finally, invoke \code{PedGFLMM} with the family data and
#'  genotype matrix to compute the \code{PedGFLMM_beta_smooth_only} statistics.  Save the p-values for each
#'  statistic in the \emph{envir$PedGFLMM_results} data frame.
#'
#' @param markers_arg a data.frame with the following 5 observations:
#' \describe{
#' \item{locus_link}{is the ordinal ranking of this marker among all loci}
#' \item{locus_link_fill}{is the position of corresponding genotype data in the
#' \emph{unified_genotype_table}}
#' \item{MarkerName}{is the text name of the marker}
#' \item{chromosome}{is the integer chromosome number}
#' \item{position}{is the integer base pair position of marker}
#'  }
#'
#' @param range_arg one row of a ranges_arg.  The latter is a data frame of at least three
#'  integer columns.  The columns indicate a range:
#'  a chromosome number, a start base pair value, and an end base pair value.
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return None
#' @import Mega2R
#' @importFrom Mega2R applyFnToRanges
#' @export
#'
#' @note
#'  This function computes the PedGFLMM_beta_smooth_only statistics and appends the output to the data frame, \emph{envir$PedGFLMM_results}.  It will
#'  print out the lines as they are generated if \emph{envir$verbose} is TRUE. The data frame
#'  \emph{envir$PedGFLMM_results} is initialized by \emph{init_PedGFLMM}, and is appended to
#'  each time \emph{DOPedGFLMM} is run.
#'
#' @seealso \code{\link{init_PedGFLMM}}
#'
#' @examples
#' db = system.file("exdata", "seqsimmGFLMM.db", package="PedGFLMM")
#' ENV = init_PedGFLMM(db)
#' ENV$verbose = TRUE
#' Mega2R::applyFnToRanges(DOPedGFLMM, ENV$refRanges[50:60,], ENV$refIndices)
#'
#' \donttest{
#' # Not run
#' Mega2R::applyFnToGenes(DOPedGFLMM, genes_arg = c("CEP104"))
#' }
#'
DOPedGFLMM = function(markers_arg, range_arg, envir = ENV) {

    if (is.null(range_arg))
        stop("DOPedGFLMM: range is not defined.", calls. = FALSE)

    geno_arg = getgenotypesraw(markers_arg, envir);

    markerNames = markers_arg$MarkerName
    gene  = as.character(range_arg[,envir$refCol[4]])

    di = dim(geno_arg)
    geno = matrix(0, nrow = (di[1]), ncol = di[2])
    for (k in 1:(di[2])) {
        vec = envir$mt2[match(as.integer(geno_arg[ , k]), envir$mt1)]
        g0 = sum(vec == 0)
        g1 = sum(vec == 1)
        g2 = sum(vec == 2)
        if (envir$verbose)
            cat(gene, markerNames[k], g0, g1, g2, "\n")
        if (g0 < g2) {
           geno[ , k] = 2 - vec
        } else {
           geno[ , k] =     vec
        }
    }

    geno = matrix(geno, nrow = di[1])
    maf = colMeans(geno)
    pos = markerNames[maf > 0]
    pheno <- mkphenotype()
    cov <- pheno[,-1*which(names(pheno)==envir$traitname)]
    names(cov)[1] <- "ped"
    names(cov)[2] <- "person"

    if (length(pos) >= 2) {       # at least 2 polymorphic variants #
        geno <- geno[ , maf > 0]     # remove non-polymorphic variants #
        nsnp    <- ncol(geno)
        # weight <- rep(1, ncol(geno))

        pedgeno <- cbind(envir$pedPer, geno)
        ID=1:nrow(envir$schaidPed)
        IDped <- cbind(ID,envir$schaidPed)

        betabasis_Bsp = 3
        genobasis_Bsp = 3
        order  =   3

        bsmooth.Bsp <- PedGFLMM_beta_smooth_only(ped = IDped, geno = pedgeno,
                 pos = markers_arg$position[maf>0], order = order, beta_basis=betabasis_Bsp,
                 covariate = cov,
                 base = "bspline")
        LRT.bsm.Bsp <- bsmooth.Bsp$LRT
#       Wald.bsmooth.Bsp <- bsmooth.Bsp$Wald

        betabasis_Fsp = 3
        genobasis_Fsp = 3

        bsmooth.fsp=PedGFLMM_beta_smooth_only(ped = IDped, geno = pedgeno,
                pos = markers_arg$position[maf>0], order = order, beta_basis=betabasis_Fsp,
                covariate = cov,
                base = "fspline")
        LRT.bsm.fsp <- bsmooth.fsp$LRT
#       Wald.bsm.fsp <- bsmooth.fsp$Wald


        ## read out the results ##
        chr   <- as.character(range_arg[,envir$refCol[1]])
        start <- range_arg[,envir$refCol[2]]
        end   <- range_arg[,envir$refCol[3]]

        result = list(chr, gene, nsnp, start, end,
                      LRT.bsm.Bsp, LRT.bsm.fsp)
        lastp1 = nrow(envir$PedGFLMM_results) + 1
        envir$PedGFLMM_results[lastp1,] = result
        if (envir$verbose) {
            print(envir$PedGFLMM_results[lastp1, ])
        }

    } else {
        if (envir$verbose)
            message("Only one marker in range.  Ignored!\n")
    }
}
