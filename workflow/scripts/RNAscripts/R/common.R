#' Splits string at -vs- to extract the two contrast groups.
#'
#' @param cont string containing the contrast groups
#' @return A character vector of length two with the contrast split apart
#' @export
split_contrast_groups <- function(cont) {
  stringr::str_split(cont, "-vs-", simplify = T) %>% unlist()
}

#' Joins the DEseq table and the fpkm read table together
#'
#' @param DEseq_tb Deseq result table from DEseq. Saved as TSV table by deseq rule
#' @param fpkm FPKM table from snakemake pipeline. You can use any fpkm table.
#' @return sorted combined tibble between the input files.
#' @importFrom magrittr %>%
#' @export
join_tables <- function(DEseq_tb, fpkm) {

  diff_exp <- DEseq_tb

  filer <- fpkm %>% dplyr::filter(gene %in% diff_exp$gene_id) %>% dplyr::filter(!duplicated(gene))
  joined_df <- dplyr::inner_join(filer, DEseq_tb, by = c('gene' = 'gene_id')) %>% dplyr::arrange(padj)

  return(joined_df)
}
#' Calculates mean values from the Four groups, which are still hardcoded. and returns
#' the mean values plus gene names and their adjusted p values
#'
#' @param mat Object that can be safely coerced to tibble with col names from Vascaging
#' @param extra_col Extra column that should be extracted and added
#' @param contrast_groups contrast_groups
#' @param s_map Sample map of samples
#' @export
#' @return tibble with mean values, gene names and padj values
#' @importFrom magrittr %>%
mean_tibble_from_mat <- function(mat, extra_col = NULL, contrast_groups = contrast_groups, s_map) {
  mat_tb <- tibble::as_tibble(mat)

  mean_tb <- tibble::tibble(.rows = nrow(mat_tb))
  for( gr in contrast_groups) {
    col_n <- s_map %>% dplyr::filter(condition == gr) %>% dplyr::pull(sample)
    mean_tb[[gr]] <- rowMeans(dplyr::select(mat_tb, col_n ))
  }
  if("gname" %in% colnames(mat_tb)) {
    mean_tb["gname"] <- mat_tb %>% dplyr::select('gname')
  }
  if("padj" %in% colnames(mat_tb)) {

    mean_tb["padj"] <- mat_tb %>% dplyr::select("padj")
  }
  if (!is.null(extra_col)) {
    mean_tb[extra_col] <- mat_tb %>% dplyr::select(extra_col)
  }

  mean_tb
}


#' test
#'
#' @param info Info field to check
#' @param index index that is to be checked
#'
#' @return
#' @export
#'
#' @examples
#' NULL
check_against_index <- function(info, index) {
  return(as.numeric(index %in% info))
}


#' Create dynamic table from significant gene names
#'
#' @param sig_gene_names import list of genes to put in matrix
#'
#' @return DT::datatable of genes across comparison groups
#' @export
#'
#' @examples
#' create_overlap_matrix(list("g1" = c("gene1", "gene2, "gene3"), "g2" = c("gene1")))
#' 
create_overlap_matrix <- function(sig_gene_names) {
  
  gene_name_index <- unlist(sig_gene_names) %>% unique() %>% sort()
  


  gene_matrix <- sapply(sig_gene_names, check_against_index, index = gene_name_index)
  rownames(gene_matrix) <- gene_name_index  
  
  row_sum <- rowSums(gene_matrix)
  
  gene_matrix <- cbind(gene_matrix, total=rowSums(gene_matrix))
  
  ordered_matrix <- gene_matrix[order(gene_matrix[,'total'], decreasing = T),]
  
  DT::datatable(ordered_matrix)
  
}
#' filters table so only contrasted samples remain
#'
#' @param rld DESeq2 object that has been normalized with either (rlog2 or vst)
#' @param contrast_groups character vecctor of contrast groups being analyzed.
#' @param j_df joined tibble of deseq + fpkm values to match genes and gnames
#' @param reorder_genes whether the order of genes should be orded by logfolgchange
#' @return matrix of epxression values
#' @export
filter_split_table <- function(rld, contrast_groups, j_df, reorder_genes = TRUE) {
  expression_values <- as.data.frame(SummarizedExperiment::assay(rld)[j_df$gene,rld$condition %in% contrast_groups, drop = F])
  rownames(expression_values) <- j_df$gname
  if (reorder_genes) {
    expression_values <- expression_values[order(abs(j_df$logFoldChange)),]
  }
  expression_values
}

#' Save pheatmapplot to svg  file
#' 
#' @param x plot object to save
#' @param filename name where to save file
#' @param width plotwidth
#' @param height plotheight
#' @return nothing
#' @importFrom grDevices dev.off svg
#' @export
save_pheatmap_svg <- function(x, filename,width = 7 , height = 7 ) {
  svglite::svglite(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
} 

#' Save any plot to svg  file
#' 
#' @param x plot object to save
#' @param filename name where to save file
#' @param width plotwidth
#' @param height plotheight
#' @return nothing
#' @importFrom grDevices dev.off svg
#' @export
save_plot_svg <- function(x, filename, width = 7 , height = 7 ) {
  svglite::svglite(filename, width = width, height = height)
  x 
  dev.off()
}

#' Save pheatmapplot to svg  file
#' 
#' @param x plot object to save
#' @param filename name where to save file
#' @param width plotwidth
#' @param height plotheight
#' @return nothing
#' @importFrom grDevices dev.off svg
#' @import svglite
#' @export
save_cheatmap_svg <- function(x, filename, width = 7 , height = 7 ) {
    svglite::svglite(filename, width = width, height = height)
  ComplexHeatmap::draw(x)
  dev.off()
}
