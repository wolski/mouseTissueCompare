
#'
#'
#'Selects those columns with same strain information
#'
#'@param all matrix starting wit columns named either Mito_ or Tissue_
#'
#'@export
#'
selectMatchStrains <- function(all){
  rownames(all) = sub("^1\\/","",rownames(all))
  tiss <- all[,grep("Tissue_",colnames(all))]
  colnames(tiss)<-gsub("Tissue_","", colnames(tiss))

  mito <- all[,grep("Mito_",colnames(all))]
  colnames(mito)<-gsub("Mito_","", colnames(mito))

  common<-intersect(colnames(tiss),colnames(mito))
  tiss <- tiss[,colnames(tiss) %in% common]
  mito <- mito[,colnames(mito) %in% common]
  tiss <- tiss[,order(colnames(tiss))]
  mito <- mito[,order(colnames(mito))]

  stopifnot(colnames(tiss) == colnames(mito))

  colnames(tiss)<- paste("Tissue_",colnames(tiss),sep="")
  colnames(mito)<- paste("Mito_",colnames(mito),sep="")
  all <- cbind("Tissue_"=tiss, "Mito_"=mito)
}

#' filter ms exeperiment and aggreage peptides proteins
#' @param exp - msexperiment
#' @param scorethresh - threshold for mscore
#' @param n - number of samples with mscore below scorethreshold
#' @param TOPPeptidesPerProtein - number of peptides per protein to use for quantification
#' @export
filterAndAggregate <- function(exp,scorethresh =0.005, n = 2,TOPPeptidesPerProtein=5){
  exp = filterByScore(exp, scorethresh =scorethresh, n = n)
  exp = aggregatePeptides(exp)
  exp = orderByRT(exp)
  # remove NA's
  tmp = apply(exp$Intensity, 1, function(x){sum(is.na(x)) == 0 })
  exp = subset( exp, tmp )
  # select top 5 peptides
  exp = selectTopPeptidesPerProtein(exp,peptop = TOPPeptidesPerProtein)
  exp = aggregateProteins(exp)
  return(exp$Intensity)
}

#' keep samples containing string
#' @param tiss substring the sample must contain, it is colled tiss because it is the tissue name usually
#' @param exp ms experiment
#' @export
#'
keepSample <- function(exp,tiss="Brain"  ){
  tmp <- colnames(exp$score)[!grepl(tiss, colnames(exp$score))]
  exp <- removeSamples(exp,tmp)
  return(exp)
}

#' load eperiment, aggregate peptides from top transitions etc.
#'
#' removes decoys
#' keeps only proteotypic peptides
#' @param SpecLib path to openswath output
#' @param maxNRTransitions max transitions per peptide
#' @param maxNRPeptides nr peptides per protein to keeps
#' @export
#'
prefilterNewDataSet <- function(SpecLib, maxNRTransitions = 3, maxNRPeptides = 10){
  obj = fread(SpecLib)
  # number of transitions / peptide
  nrt = maxNRTransitions
  # nr of peptides per protein
  peptop = maxNRPeptides

  obj=prepareDF(obj)
  colnames(obj)
  msexpori = loadTransitonsMSExperiment(obj,nrt=nrt,peptop=peptop)
  length(unique(msexpori$pepinfo$transition_group_id))
  stopifnot(nrt==length(msexpori$pepinfo$transition_group_id) / length(unique(msexpori$pepinfo$transition_group_id)))
  msexpori = removeDecoys(msexpori)
  tmp = grep("^1\\/",msexpori$pepinfo$ProteinName)
  msexpori = subset(msexpori, tmp)
  return(msexpori)
}

#' annotate sample
#'
#' @param msexpori msexperiment
#' @param FileSampleMapping mapping of files to meaningfull label
#' @param Samples2Remove sample to remove
#' @export
#'
annotateSamples <- function(msexpori,FileSampleMapping,Samples2Remove=""){
  msexpori = removeSamples(msexpori, sample=Samples2Remove)
  ncn = read.table(FileSampleMapping,header=T,stringsAsFactors = FALSE)
  ncn = ncn[ ncn$Filename != Samples2Remove, ]
  ncnOrdering <- match(colnames(msexpori$Intensity),ncn$Filename )
  stopifnot(ncn$Filename[ncnOrdering] == colnames(msexpori$Intensity))
  msexpori = mycolnames(msexpori,ncn[ncnOrdering,1])
  return(msexpori)
}


