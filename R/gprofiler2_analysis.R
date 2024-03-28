#' Imports a standard resulting data.frame from DESeq2
#'
#' @param nombre an \code{integer} that can take the value 1 or 2 specifying if you want to
#' import only one \code{data.frame} or a \code{list} of two data.frames
#' Default: 1
#'
#' @return if no argument given the object returned is a .csv file that contains a \code{data.frame} named df1 with:
#' df1 specifications:
#' 25 rows with a padj < 0.05 (significant) & a abs(log2FoldChange) > log2(1.5)
#' 25 rows with a padj < 0.05 (significant) & a abs(log2FoldChange) < log2(1.5)
#' 25 rows with a padj > 0.05 (not significant) & a abs(log2FolChcange) > log2(1.5)
#' 25 rows with a padj > 0.05 (not significant) & a abs(log2FoldChange) < log2(1.5)
#' if argument given (number=2) the object returned is a \code{list} of two .csv files that contain a \code{data.frame} named df1 and df2.
#' df2 specifications:
#' 50 rows with a padj < 0.05 (significant) & a abs(log2FoldChange) > log2(2)
#' 50 rows with a padj < 0.05 (significant) & a abs(log2FoldChange) < log2(2)
#' 50 rows with a padj > 0.05 (not significant) & a abs(log2FolChcange) > log2(2)
#' 50 rows with a padj > 0.05 (not significant) & a abs(log2FoldChange) < log2(2)
#'
#' @examples
#' de_res <- get_demo_de_res()
#' de_res <- get_demo_de_res(2)
#'
#' @importFrom readr read_csv
#'
#' @export

get_demo_de_res <- function(number=1){
  stopifnot(number %in% c(1, 2))
  if(number  == 2){
    return(list(df1=read_csv(system.file("extdata/demo_de_res.csv", package = "rnaseq")), df2=read_csv(system.file("extdata/demo_de_res2.csv", package = "rnaseq"))))
  } else {
    return(df1=read_csv(system.file("extdata/demo_de_res.csv", package = "rnaseq")))
  }
}




#' Produces a gene profile from DESeq2 results, converts the result for clusterProfiler visualization
#'
#' @param de_res either a \code{data.frame} or a list of \code{data.frame} returned by the
#' \code{DESeq2::results} function.
#' @param p_threshold the numeric threshold of the adjusted p-value to be considered as significant.
#' the value has to be between 0 and 1.
#' Default: 0.05
#' @param fc_threshold the numeric threshold of the fold change to be considered as significant.
#' the value has to be positive.
#' Default: 1.5
#' @param organism the string that defines corresponding source organism for the gene list.
#' The organism names are usually constructed by concatenating the first letter of the name
#' and the family name, e.g human - hsapiens.
#' Default : "hsapiens"
#'
#' @return a list of two lists. list$gp is the \code{list} returned by \code{gprofiler::gost} that contains
#' either one or many gene profiles(gp) objects depending on the parameter de_res.
#' list$cp is the \code{list} returned by the \code{gprofiler2_convert} function.
#' the list$cp can contain one or many clusterProfiler gene profiles (cp) depending on the parameter de_res.
#' Each cp has one enrichResult object returned by \code{multienrichjam::enrichDF2enrichResult} and one \code{data.frame}
#' df generated from the result \code{list} in the enrichResult object
#'
#' @examples
#' de_res <- get_demo_de_res()
#' enr <- gprofiler2_analysis(de_res)
#'
#' @importFrom gprofiler2 gost
#' @importFrom clusterProfiler merge_result
#' @importFrom clusterProfiler dotplot
#' @importFrom purrr map
#' @importFrom purrr walk
#' @importFrom tidyr tibble
#'
#' @export

gprofiler2_analysis <- function(de_res, p_threshold = 0.05, fc_threshold = 1.5, organism = "hsapiens"){

  # Testing the input

  ##de_res
  validate_de_res <- function(x){
    stopifnot(nrow(x) > 0)
    expected_values <- c("ensembl_gene", "padj", "log2FoldChange")
    stopifnot(expected_values %in% colnames(x))
    stopifnot(is.character(x$ensembl_gene))
    stopifnot(!any(is.na(x$ensembl_gene)))
  }

  if(is.data.frame(de_res)){
    validate_de_res(de_res)
    de_res <- list(de_res = de_res)
  } else if(is.list(de_res)){
    stopifnot(length(names(de_res)) == length(unique(names(de_res))))
    stopifnot(!("" %in% allNames(de_res)))
    walk(de_res, ~ validate_de_res(.x))
  } else {
    stop("de_res should be data.frame or list")
  }

  ##p_threshold
  stopifnot(is.numeric(p_threshold))
  stopifnot(0 < p_threshold & p_threshold < 1)

  ##fc_threshold
  stopifnot(is.numeric(fc_threshold))
  stopifnot(fc_threshold > 0)

  ##organism
  stopifnot(is.character(organism))

  functional_enrichment <- function(current_de_res){

    #(1.0) Bloc that filters data according to set parameters
    signif_de_res <- subset(current_de_res, padj < p_threshold & abs(log2FoldChange) >= log2(fc_threshold))


    #(2.0) Bloc that does the Hypergeometric analysis using gprofiler2

    if(nrow(signif_de_res) == 0){
      empty_df <- setNames(data.frame(matrix(ncol = ncol(signif_de_res), nrow = 0)), colnames(signif_de_res))
      return(empty_df)
    }

    gp <- gprofiler2::gost(signif_de_res$ensembl_gene,
                           organism = organism,
                           ordered_query = FALSE,
                           multi_query = FALSE,
                           significant = TRUE,
                           evcodes = TRUE,
                           user_threshold = p_threshold,
                           sources = c("GO:BP", "REAC", "KEGG", "WP"))
    if(is.null(gp)){
      gp <- list(result=as.data.frame(tibble(query=character(),
                                             significant=logical(),
                                             p_value=double(),
                                             term_size=integer(),
                                             query_size=integer(),
                                             intersection_size=integer(),
                                             precision=double(),
                                             recall=double(),
                                             term_id=character(),
                                             source=character(),
                                             term_name=character(),
                                             effective_domain_size=integer(),
                                             source_order=integer(),
                                             parents=list(),
                                             evidence_codes=character(),
                                             intersection=character(),
                                             GeneRatio=character())), meta=list())

      cp <- list(enrichResult = NULL, df=as.data.frame(tibble(ID=character(),
                                                              query=character(),
                                                              significant=logical(),
                                                              pvalue=double(),
                                                              term_size=integer(),
                                                              query_size=integer(),
                                                              intersection_size=integer(),
                                                              precision=double(),
                                                              recall=double(),
                                                              source=character(),
                                                              Description=character(),
                                                              effective_domain_size=integer(),
                                                              source_order=integer(),
                                                              parents=list(),
                                                              evidence_codes=character(),
                                                              geneID=character(),
                                                              GeneRatio=character(),
                                                              p.adjust=double(),
                                                              Count=integer(),
                                                              setSize=double())))

    } else{
      gp$result$GeneRatio = paste0(gp$result$intersection_size, "/", gp$result$query_size)

      cp <- list(enrichResult=gprofiler2_convert(gp))

      cp$df = cp$enrichResult@result

    }
    return(list(gp=gp, cp=cp))

  }
  res <- map(de_res, ~ functional_enrichment(.x))
  final_res <- list(gp = map(res, ~ .x$gp), cp = map(res, ~ .x$cp))
  return(final_res)
}




#' Converts the result of gprofiler2::gost in a compatible object for clusterProfiler
#'
#' @param gp the \code{list} returned by the \code{gprofiler::gost}
#' function wich is a list of two lists.
#' result contains a data.frame with the results from the hypergeometric test.
#' meta contains lists of the inputs and the queries for reproductibility.
#'
#' @return cp is a converted \code{list} from the \code{list} returned by the \code{gprofiler::gost}
#' function and is compatible with the vizualisation function from the code{clusterProfiler} package.
#'
#' @importFrom multienrichjam enrichDF2enrichResult
#'
#' @export

gprofiler2_convert <- function(gp){

  # Testing the input

  ##gp
  stopifnot(is.list(gp))
  stopifnot(length(gp) == 2)
  stopifnot((!any(is.na(gp))))

  ##gp$result
  stopifnot(is.list(gp$result))
  stopifnot(length(gp$result) == 17)
  expected_values <- c("query", "significant", "p_value",
                       "term_size", "query_size", "intersection_size",
                       "precision", "recall", "term_id", "source",
                       "term_name", "effective_domain_size", "source_order",
                       "parents", "evidence_codes", "intersection")
  stopifnot(expected_values %in% colnames(gp$result))

  ##gp$result$parents
  stopifnot(is.list(gp$result$parents))

  ##gp$meta
  stopifnot(is.list(gp$meta))

  #(1.0) Bloc that converts the result from gprofiler2::gost in a compatible object for clusterProfiler visualization
  cp <- multienrichjam::enrichDF2enrichResult(enrichDF = gp$result, keyColname = "term_id",
                              geneColname = "intersection", pvalueColname = "p_value",
                              descriptionColname = "term_name", pvalueCutoff = 0.05, geneRatioColname = "GeneRatio")
  return(cp)
}
