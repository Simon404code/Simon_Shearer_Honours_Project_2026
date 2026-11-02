
library(GENESPACE)
library(BiocGenerics)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(GenomeInfoDbData)
library(rtracklayer)
library(data.table)
library(XVector)
library(Biostrings)
library(ggplot2)
library(devtools)
library(parallel)
library(fpc)
library(dbscan) 
library(dplyr)

genomeRepo <- "~/h_project/50_genespace_run/rss/"
wd <- "/scratch/genespace2"
path2mcscanx <- "~/h_project/50_genespace_run/MCScanX/"
genomes2run <- c("RK05HAP1", "RK05HAP2", "RK17HAP1", "RK17HAP2", "RK205HAP1", "RK205HAP2", "RK25HAP1", "RK25HAP2", "RK313HAP1", "RK313HAP2", "RK320HAP1", "RK320HAP2", "RK339HAP1", "RK339HAP2", "RK388HAP1", "RK388HAP2", "RK436HAP1", "RK436HAP2", "RK441HAP1", "RK441HAP2", "RK465HAP1", "RK465HAP2", "RK498HAP1", "RK498HAP2", "RK50HAP1", "RK50HAP2", "RK561HAP1", "RK561HAP2", "RK614HAP1", "RK614HAP2", "RK687HAP1", "RK687HAP2", "RK715HAP1", "RK715HAP2", "RK74HAP1", "RK74HAP2", "RK792HAP1", "RK792HAP2", "RK84HAP1", "RK84HAP2", "RK86HAP1", "RK86HAP2", "RM01HAP1", "RM01HAP2", "RM02HAP1", "RM02HAP2", "RM03HAP1", "RM03HAP2", "TAG0014HAP1", "TAG0014HAP2")

parsedPaths <- parse_annotations(rawGenomeRepo=genomeRepo, genomeDirs=genomes2run, genomeIDs= genomes2run, presets = "phytozome", genespaceWd=wd)

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx)

out <- run_genespace(gpar, overwrite = T)

ripd <- plot_riparian(
  gsParam = out,
  refGenome = "TAGHAP1", 
  useRegions = FALSE)

#############Anneri's code below 
genomeRepo <- "50_genespace_run/other_rss/files_fa/"
faFiles <- list.files(path = genomeRepo, recursive = T, pattern = ".fa", full.names = TRUE)

names(faFiles) <- gsub("Eucalyptus_grandis_var_|UP_FAM[0-9]+_|\\.mainGenome.fasta", "", basename(faFiles))


source("clean_windows.R") ### Load updated library because source functions do not respect namespace
source("riparian_windows.R")
wnd <- init_windows(wd = wd, nCores = 6, faFiles = faFiles, path2mcscanx = path2mcscanx)

out <- clean_windows(wnd)

bcf <- list.files(path = file.path(wd, "output"), pattern = "synBlkCoords", full.names = T)
rip <- riparian_windows(faFiles = faFiles, blockCoordFiles = bcf)