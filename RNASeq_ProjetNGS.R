#Projet NGS (Next Generation of Sequencing) 2022
#auteurs: P.Lucie, R.Vanessa, T.V.C.Valérie, S.Roussaya
#version 10/02/2022
#########################################################################################
#    Pipeline d'analyses d'expressions différentielle de gènes issus de données RNASeq
#
#
# Chargement des différentes librairies. NB: utilisé packages_istall.R pour les installer

library(limma) #Data analysis, linear models and differential expression for microarray and RNASeq data

library(edgeR) #Differential expression analysis of RNA-seq expression profiles with biological replication. 
#Implements a range of statistical methodology based on the negative binomial distributions, 
#including empirical Bayes estimation, exact tests, generalized linear models and quasi-likelihood tests.

library(Rsubread) #Includes functionality for read mapping, read counting, SNP calling, structural variant detection 
#and gene fusion discovery

library(SRAdb) #The Sequence Read Archive (SRA) is the largest public repository of sequencing data from the NGS

########### Etape 1: Importer les données ########################
url <- "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-9285/E-MTAB-9285.sdrf.txt"
data <- read.delim(url)

########### Etape 2: Contrôles qualités et filtration ######################## https://bioinformatics-core-shared-training.github.io/RNAseq-R/align-and-count.nb.html
# Extract quality scores for 100 reads
qs <- qualityScores(filename="",nreads=100)
# Check dimension of qs
dim(qs)
# Check first few elements of qs with head
head(qs)


########### Etape 3: Gene annotation (if necessary) #############

########### Etape 4: Alignement des read ########################
FASTQFiles <- paste0(Target$indivudial,"_",Target$center, ".fastq.gz")
BAMFiles <- paste0(Target$individual, "_", Target$center, "subread.bam")
buldindex(basename="hg19_subread",reference="GRCh37.fa")
align ("hg19_subread", readfile1 = FastQFiles, input_format="gzFASTQ", output_file=BAMFiles)
gene <- feateureCounts(BAMFiles, useMetaFeature=TRUE, annot.inbuilt="hg19", allowMultiOverlaps=TRUE)
propmap <- propmapped(bamFiles)

########### Etape 5: Matrice des read counts ####################
# filter out excluded samples (marked as "X")]...
#TODO: Why?? adapt to our project....
sampleSelection <- sampleSelection[which(sml != "X")]
# gs is a  Factor (= labelized vector, for quantitative data manipulation  )
gs <- factor(sampleSelection)
levels(gs) <- groups

# build the design matrix having size : 'selectedSampleIdList' rows  X 'groups' columns
# first we use the transpose dataFrameCounts, because we need GSM as variables ... using matrix representation for transpose t()
dataFrameCounts <- data.frame(
  matrix(vector(), 0, length(sampleIdList),  dimnames=list(c(), c( sampleIdList))),
  stringsAsFactors=F)
print( summary(dataFrameCounts) )


mtdataFrameCounts <- t(data.matrix(dataFrameCounts))
mtdataFrameCounts <- as.data.frame(mtdataFrameCounts)
# then add column 'group' to mtdataFrameCounts
mtdataFrameCounts$group <- gs

# finally build design
design <- model.matrix(~group + 0, mtdataFrameCounts)
colnames(design) <- levels(gs)
# samples 2 groups
sample2group <- apply(design, 1, function(v)  colnames(design)[which(v==1)] )

#objet DGEList
dge <- DGEList(counts=dataFrameCounts, genes=rownames(dataFrameCounts))

########### Etape 6: Normaliser #######################
dge <- calcNormFactors(dge)

########### Etape 7: Définir les contraste ######################

########### Etape 8: Conversion logarithmique (log2) ############
l2 <- log2(dataframe +2)

########### Etape 9: Appliquer modèles statistiques #############

########### Etape 10: Génerer les visualisations ################





