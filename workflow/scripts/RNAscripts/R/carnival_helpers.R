#' Convert Human gene list
#'
#' @param x character vektor of gene symbols to translate (DEPrecated)
#' @param hmart hunan biomart
#' @param mmart mouse biomart
#'
#' @return
#' @export
#'
#' @examples
#' NULL
convertHumanGeneList <- function(x, hmart, mmart) {
  genesV2 <- biomaRt::getLDS(
    attributes = c("hgnc_symbol"),
    filters = "hgnc_symbol",
    values = x,
    mart = hmart,
    attributesL = c("mgi_symbol"),
    martL = mmart,
    uniqueRows = T
  )
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  # print(head(humanx))
  return(humanx)
}

#' convert between human and mouse homologs using homolog table
#'
#' @param x List of genes
#' @param jax_database Database from outside to reduce http calls
#' @param input_organism Organism being converted from (Homo sapiens or Mus musculus)
#'
#' @return
#' @export
#'
#' @examples
#' convertHumanGeneHomologs(c("APLNR", "APLNR"))
convertHumanGeneHomologs <-
  function(x,
           jax_database = NULL,
           input_organism = "Homo sapiens") {
    if (is.null(jax_database)) {
      mouse_human_homologs <-
        readr::read_tsv(
          "http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt",
          col_names = c("hgene", "hID", "mgene", "mID", "lcol")
        )
    } else {
      mouse_human_homologs <- jax_database
    }
    if (input_organism == "Homo sapiens") {
      gene_list <- mouse_human_homologs %>%
        dplyr::filter(hgene %in% x) %>%
        dplyr::pull(mgene)
    } else if (input_organism == "Mus musculus") {
      gene_list <- mouse_human_homologs %>%
        dplyr::filter(mgene %in% x) %>%
        dplyr::arrange(match(mgene, x)) %>%
        dplyr::pull(hgene)
    } else {
      stop("input_organism not supported choose Mus musculus or Homo sapiens")
    }

    gene_list
  }


#' Generate TF list from Saezlab Transcript tutorial
#'
#' @param df Dataframe from Dorothea
#' @param top Top n genes to
#' @param access_idx Col to access
#'
#' @return
#' @export
#'
#' @examples
#' NULL
generateTFList <- function(df = df,
                           top = 50,
                           access_idx = 1) {
  if (top == "all") {
    top <- nrow(df)
  }
  if (top > nrow(df)) {
    warning(
      "Number of to TF's inserted exceeds the number of actual TF's in the\n            data frame. All the TF's will be considered."
    )
    top <- nrow(df)
  }
  ctrl <- intersect(x = access_idx, y = 1:ncol(df))
  if (length(ctrl) == 0) {
    stop("The indeces you inserted do not correspond to \n              the number of columns/samples")
  }
  returnList <- list()
  for (ii in seq_along(ctrl)) {
    tfThresh <- sort(x = abs(df[, ctrl[ii]]), decreasing = TRUE)[top]
    temp <- which(abs(df[, ctrl[ii]]) >= tfThresh)
    currDF <- matrix(
      data = ,
      nrow = 1,
      ncol = top
    )
    colnames(currDF) <- rownames(df)[temp[1:top]]
    currDF[1, ] <- df[temp[1:top], ctrl[ii]]
    currDF <- as.data.frame(currDF)
    returnList[[length(returnList) + 1]] <- currDF
  }
  names(returnList) <- colnames(df)[ctrl]
  return(returnList)
}
#' assigny progeny scores to genes via OmnipathR
#'
#' @param progeny table of Progeny output
#' @param progenyMembers List of Members from Progeny
#' @param id ID of gnees
#' @param access_idx Score Access ID
#'
#'
#' @return
#' @export
#'
#' @examples
#' NULL
assignPROGENyScores <-
  function(progeny = progeny,
           progenyMembers = progenyMembers,
           id = "gene",
           access_idx = 1) {
    if (id == "uniprot") {
      idx <- which(names(progenyMembers) == "uniprot")
      progenyMembers <- progenyMembers[[idx]]
    } else {
      idx <- which(names(progenyMembers) == "gene")
      progenyMembers <- progenyMembers[[idx]]
    }
    members <- matrix(data = , nrow = 1, ncol = 2)
    pathways <- colnames(progeny)
    ctrl <- intersect(x = access_idx, y = 1:nrow(progeny))
    if (length(ctrl) == 0) {
      stop("The indeces you inserted do not correspond to \n              the number of rows/samples")
    }
    for (ii in seq_along(pathways)) {
      mm <- progenyMembers[[which(names(progenyMembers) ==
        pathways[ii])]]
      for (jj in seq_along(mm)) {
        members <- rbind(members, c(pathways[ii], mm[jj]))
      }
    }
    members <- members[-1, ]
    scores <-
      matrix(
        data = ,
        nrow = nrow(progeny),
        ncol = nrow(members)
      )
    colnames(scores) <- members[, 2]
    rownames(scores) <- rownames(progeny)
    members <- unique(members)
    for (i in 1:ncol(scores)) {
      for (j in 1:nrow(scores)) {
        scores[j, i] <- as.numeric(progeny[j, members[which(members[
          ,
          2
        ] == colnames(scores)[i]), 1]])
      }
    }
    pxList <- list()
    for (ii in seq_along(access_idx)) {
      pxList[[length(pxList) + 1]] <-
        as.data.frame(t(as.matrix(scores[access_idx[ii], ])))
    }
    names(pxList) <- rownames(progeny)[ctrl]
    return(pxList)
  }

#' \code{df_to_viper_regulon}
#'
#' This function is designed to generate a ready to use regulon object for viper
#' from a 3 column dataframe representation of a target set collection.
#'
#' @param df a dataframe of n*3 dimension. The first column corresponds the targets,
#' and the second column indicates which regulon does each target belongs to.
#' The third column corresponds to the weight and sign of the interaction between
#' a regulon and its targets.
#'
#' @return a list where each element is a regulon in the viper format. This list
#' is ready to be used as a regulon set in viper.
df_to_viper_regulon <- function(df) {
  names(df) <- c("feature", "pathway", "sign")
  df <- df[complete.cases(df), ]

  pathway_regulon <- list(0)
  i <- 1
  for (pathway in unique(df$pathway))
  {
    pathway_feature_list <- list(0)
    features <- df[df$pathway == pathway, 3]
    names(features) <- df[df$pathway == pathway, 1]
    pathway_feature_list[[1]] <- features
    pathway_feature_list[[2]] <- rep(1, length(features))
    names(pathway_feature_list) <- c("tfmode", "likelihood")

    pathway_regulon[[i]] <- pathway_feature_list
    i <- i + 1
  }
  names(pathway_regulon) <- unique(df$pathway)
  return(pathway_regulon)
}

#' rename ENSG IDS from OTP using fpkm table
#'
#' @param count_table Table with ENSG IDs with version as rownames
#' @param fpkm_table FPKM table from rna-star-deseq pipeline `fpkm/all.tsv` as tibble
#' @param fpkm_col column name of genes in fpkm table
#' @param tibble_col If count_table is tibble give tibble col
#' @return
#' @export
#'
#' @examples
#' NULL
rename_count_rownames <- function(count_table, fpkm_table, fpkm_col = "gene", tibble_col = NULL) {

  # Check if tibble or data.frame
  if(tibble::is_tibble(count_table)) {
    stopifnot("tibble_col must be defined when count_table is tibble" = 
                !is.null(tibble_col))
    filer <- fpkm_table %>% dplyr::filter(!duplicated(gname)) %>% dplyr::select(gene, gname)                               

    count_table %>% dplyr::filter(!!sym(tibble_col) %in% filer$gene) -> count_table
    by_vec <- fpkm_col
    names(by_vec) <- tibble_col
    dplyr::left_join(count_table, filer, by = by_vec) -> count_table
  } else {
    filer <- fpkm_table %>%
    dplyr::filter(gene %in% rownames(count_table)) %>%
    dplyr::filter(!duplicated(gname)) %>% dplyr::select(gene, gname)
    count_table <- count_table[filer$gene, ]
    rownames(count_table) <- filer$gname
  }
  count_table
}


#' Parse and change variable types for carnival
#'
#' @param carnival_res Carnival Results object
#'
#' @return
#' @export
#'
#' @examples
process_carnival <- function(carnival_res) {
  carnival_res$weightedSIF <-
    data.frame(carnival_res$weightedSIF, stringsAsFactors = F)
  carnival_res$weightedSIF$Sign <-
    as.numeric(carnival_res$weightedSIF$Sign)
  carnival_res$weightedSIF$Weight <-
    as.numeric(carnival_res$weightedSIF$Weight)

  carnival_res$nodesAttributes <-
    data.frame(carnival_res$nodesAttributes, stringsAsFactors = F)
  carnival_res$nodesAttributes$ZeroAct <-
    as.numeric(carnival_res$nodesAttributes$ZeroAct)
  carnival_res$nodesAttributes$UpAct <-
    as.numeric(carnival_res$nodesAttributes$UpAct)
  carnival_res$nodesAttributes$DownAct <-
    as.numeric(carnival_res$nodesAttributes$DownAct)
  carnival_res$nodesAttributes$AvgAct <-
    as.numeric(carnival_res$nodesAttributes$AvgAct)
  carnival_res
}


#' get organism ID used in omnipath for Mus musculus or Homo sapiens
#'
#' @param organism_name Name of organism (Mus musculus, Homo sapiens)
#'
#' @return
#' @export
#'
#' @examples NULL
get_organism_omnipath_id <- function(organism_name) {
  if (organism_name == "Homo sapiens") {
    org_num <- 9606
  } else if (organism_name == "Mus musculus") {
    org_num <- 10090
  } else {
    stop(
      "Organism not supported. Please select a supported organism \n Mus musculus or Homo sapiens"
    )
  }
  org_num
}

#' get organism name used in omnipath for Mus musculus or Homo sapiens
#'
#' @param organism_name Name of organism (Mus musculus, Homo sapiens)
#'
#' @return
#' @export
#'
#' @examples NULL
get_organism_omnipath_name <- function(organism_name) {
  if (organism_name == "Homo sapiens") {
    org_num <- "Human"
  } else if (organism_name == "Mus musculus") {
    org_num <- "Mouse"
  } else {
    stop(
      "Organism not supported. Please select a supported organism \n Mus musculus or Homo sapiens"
    )
  }
  org_num
}

#' get organism name used in ensembl for Mus musculus or Homo sapiens
#'
#' @param organism_name Name of organism (Mus musculus, Homo sapiens)
#'
#' @return
#' @export
#'
#' @examples NULL
get_organism_ensembl_name <- function(organism_name) {
  if (organism_name == "Homo sapiens") {
    org_name <- "hsapiens"
  } else if (organism_name == "Mus musculus") {
    org_name <- "mmusculus"
  } else {
    stop(
      "Organism not supported. Please select a supported organism \n Mus musculus or Homo sapiens"
    )
  }
  org_name
}
