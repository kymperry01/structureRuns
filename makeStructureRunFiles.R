# ##############################################################################
# Run structure analysis on Px and Pa populations
# K.Perry, 19/12/2016

# This script does:
# make sed commands for running in unix terminal (tizard1),
# which convert template subfiles and mainparams files to 
# create independent structure runs file for all values of 'k' and 'run'.
# 
# Instructions:
# template 'mainparams', 'extraparams' and structure subfiles must reside in the same linux directory (specified by tizDir below)
# note, the template sub must have the correct linux working directory in the 'cd' command ... was too much trouble to code this in.
# run the two R functions to write sed commands to a .txt file on local directory
# navigate to the linux directory in the unix terminal
# copy and paste sed commands into terminal to generate the new files and queue them.
# ##############################################################################

library(dplyr)
library(stringr)


# #################################################
# define functions
# #################################################

# 
# --- A function to configure the extraparams file ---

# ... the extraparams file contains the model assumptions for the structure analysis
# ... call this function once intially create a named extraparams file for your analysis
# ... thereafter, only run if you want to adjust the parameter set
# ... note: here, I've only added parameters that I may wish to change from the defaults
# ... to alter other parameters, code them into the function

# This function prints a sed command to screen.
# ... to create the new extraparams file, paste the sed command into your linux terminal
# ... in the directory with your extraparams template file.

setExtraParams <- function(noadmix = 0, # 0 assumes admixture model 
                           locprior = 1, # 1 assumes locprior mode
                           freqscorr = 1, # 1 assumes correlated allele frequencie
                           lambda = 1.0, # set lambda = 1. Works ok for most data)
                           tizDir = "/scratch/kperry/structure/",
                           epTemplate = "extraparams_template",
                           outPrefix) { # desired identifier of parameter file/param set. e.g. <prefix>_extraparams
  
  out <- paste0(outPrefix, ".extraparams")
  paste0("sed ",
         "-e ", "'s/NOADMIX 0", "/NOADMIX ", noadmix, "/' ",
         "-e ", "'s/LOCPRIOR 1", "/LOCPRIOR ", locprior, "/' ",
         "-e ", "'s/FREQSCORR 1", "/FREQSCORR ", freqscorr, "/' ",
         "-e ", "'s/LAMBDA 1.0", "/LAMBDA ", format(lambda, digits = 2, nsmall = 1), "/' ",
         paste0(tizDir, epTemplate), " > ", paste0(tizDir, out), ";")
  
} # End function

# --- A function to create sed commands for pasting directly into linux terminal 
# ... it take template mainparams and subfiles and write new ones for each k and r
sedStructureFiles <- function(minK, maxK, startRunNum, numRuns,
                             numInds, numLoci, burnins, mcmc,
                             inFile, # structure input file name
                             outPrefix, # prefix string for mainparams and subfiles
                             mpTemplate, # name of mainparams template file in tizard dir
                             extraparams, # name of extraparams file in tizard dir. Should be the output file of setExtraParams()
                             subTemplate, # name of subfile template in tizard dir
                             locOutDir, # local dir to write sed cmds and subfiles
                             # whether to write the `queue subfile` commands by K value (the default) or run number (TRUE).
                             # ordering by run is useful when you want to queue the files in run batches over time, to ensure the random seed numbers are different for each run
                             # (the structure randomize==1 takes the random seed from the system clock)
                             # ordering by run simply makes cut/paste of command minot terminal more convenient 
                             orderByRun = FALSE) { # wh order to write the `queue subfile` commands

  
  # 1. set up the first mainparams file for k = minK, run = startRunNum
  # ... then make sed commands to create a file for each 'k' and 'run'
  kSeq <- seq(minK, maxK)
  rSeq <- seq(startRunNum, startRunNum + numRuns - 1)
  kr <- paste0("_k", minK, "_", rSeq[1]) # unique combination of 'k' and 'r' values for each run
  
  mp1 <- paste0(outPrefix, kr, ".mainparams") # name of the first mainparams file
  setMainParams <- paste0("sed ", 
                          "-e ", "'s/MAXPOPS 00", "/MAXPOPS ", minK, "/' ",
                          "-e ", "'s/BURNIN 00", "/BURNIN ", format(burnins, scientific = FALSE), "/' ",
                          "-e ", "'s/NUMREPS 00", "/NUMREPS ", format(mcmc, scientific = FALSE), "/' ",
                          "-e ", "'s/INFILE infile", "/INFILE ", inFile, "/' ",
                          "-e ", "'s/OUTFILE outfile", "/OUTFILE ", mp1, "/' ",
                          "-e ", "'s/NUMINDS 00", "/NUMINDS ", numInds, "/' ",
                          "-e ", "'s/NUMLOCI 00", "/NUMLOCI ", numLoci, "/' ",
                          mpTemplate, " > ",
                          mp1, ";")
  
  # 2. Set first subfile, create sed commands to make the rest
  sub1 <- paste0(outPrefix, kr, ".sub") # name of the first sub file 
  setSubfile <- paste0("sed ",
                       "-e ", "'s/jobname/", paste0(outPrefix, kr), "/' ",
                       "-e ", "'s/mainparams/", mp1, "/' ",
                       "-e ", "'s/extraparams/", extraparams, "/' ",
                       "-e ", "'s/infile/", inFile, "/' ",
                       "-e ", "'s/outfile/", paste0(outPrefix, kr), "/' ",
                       subTemplate, " > ", 
                       sub1, ";")
  
  sed <- lapply(kSeq, function(k) {
    
    # for each run in rSeq
    lapply(rSeq, function(r){
      krNew <- paste0("_k", k, "_", r)
      mpOut <- paste0(outPrefix, krNew, ".mainparams")
      subOut <- paste0(outPrefix, krNew, ".sub")
      mainparams <- paste0("sed ", 
                           "-e ", "'s/MAXPOPS ", minK, "/MAXPOPS ", k, "/' ",
                           "-e ", "'s/", kr, "/", krNew, "/g' ",
                           mp1, " > ", mpOut, ";")
      subfiles <- paste0("sed -e ", "'s/", kr, "/", krNew, "/g' ",
                         sub1, " > ", subOut, ";")
      qsubs <- paste0("sbatch ", subOut, ";")
      data_frame(mainparams, subfiles, qsubs)
      
      }) %>%
      bind_rows()
    }) %>%
    bind_rows()
  
  orderQsubByRun <- function(x) {
    run <- str_extract(x, "\\d{1,3}\\.sub") %>% 
      gsub("\\.sub", "", .) %>% 
      as.integer()
    cbind(data.frame(run),
          data.frame(x)) %>%
      arrange(run) %>%
      dplyr::select(x) %>%
      unlist() %>%
      as.character()
  }
  
  if (orderByRun) {
    queue <- orderQsubByRun(sed$qsubs)
  } else {
    queue <- sed$qsubs
  }
  
  outpath <- paste0(locOutDir, "makeStructureFilesCommands.", outPrefix, ".txt")
  write.table(c(setMainParams, sed$mainparams[-1]), outpath,
              row.names = FALSE, quote = FALSE,
              col.names = "# Run these commands in terminal to create the mainparams and subfiles files for STRUCTURE (one file for each 'k' value and 'run')")
  write.table(c(setSubfile, sed$subfiles[-1]), outpath,
              row.names = FALSE, quote = FALSE, col.names = FALSE,
              append = TRUE) 
  write.table(queue, outpath,
              row.names = FALSE, quote = FALSE, col.names = FALSE,
              append = TRUE) 
  
  }  # End function

# ####################################
# Pa structure analysis
# 27/12/2016

# Here, we'll run full structure analysis using param set from RAD1 (locprior, correlate allele frequencies, lambda = 1) 

# In unix (tizard) terminal, navigate to /scratch/kperry/structure
# Generate extraparams file containing the model assumptions (paste the command into unix terminal)
setExtraParams(tizDir = "/scratch/kperry/structure/",
               epTemplate = "extraparams_template",
               outPrefix = "ps1") # param set 1

# **********************************************
# Pa, 943 bi-allelic SNPs only, max missing 0.7 
# **********************************************
# Create sed commands to make mainparams and subfiles for each 'k' and 'run'  
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 10,
                 numInds = 45, numLoci = 943,
                 burnins = 250000, mcmc = 500000,
                 inFile = "Pa.943.SNPs.unix.structure", 
                 outPrefix = "pa.943",
                 mpTemplate = "mainparams_template", 
                 extraparams = "ps1.extraparams",
                 subTemplate = "structure_template.sub",
                 locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/")


# **************************
# papx4pops, bi-allelic SNPs
# **************************

sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 10,
                  numInds = 89, numLoci = 843,
                  burnins = 250000, mcmc = 500000,
                  inFile = "pxpa4pops.843.SNPs.unix.structure", 
                  outPrefix = "pxpa.843",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub",
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/")




# ****************************
# make files for px first run;
# do short clustering runs for k1-15, 5 reps, 100000 burnins/reps
# 26/12/2016
# *****************************

# Do some short clustering runs using 995 biallelic SNP sites to determine optimal k
# (afterwards, do longer runs at optimal k +- 1)

sedStructureFiles(minK = 1, maxK = 20, startRunNum = 1, numRuns = 5,
                  numInds = 857, numLoci = 995,
                  burnins = 10000, mcmc = 10000, 
                  inFile = "Px.variants.hardFilt.maxMiss0.7.SNPs.unix.structure", 
                  outPrefix = "px.995",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub",
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/")

#> optimal k = 11 (!!), but there was a fair bit of noise in the data.
#> download the plots, view plots in clumpp/distruct out of interest

# Next 
# a. ... do long runs of K1-15 for all Px
# b. ... do short clustering runs using variants file that include 2 Pa

# ******************************************************************
# Long clustering runs, Px 'all', 1008 biallelic SNPs, with duplicate samples removed.
# 14/1/2017
sedStructureFiles(minK = 13, maxK = 15, startRunNum = 1, numRuns = 10,
                  numInds = 842, numLoci = 1008,
                  burnins = 250000, mcmc = 500000, 
                  inFile = "Px.1008.SNPs.unix.structure", 
                  outPrefix = "px1008",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub",
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/")


# ***************************************
# Now, run structure separately for 2014 and 2015 samples

# 2014, 440 ind
# (15/2/2017: final run, 500k burn, 10^6k mc)
# Edit: 8/10/2017: repeat, this time 20 runs using a different random seed each time (using OrderByRun = TRUE, and pasting commands in for each run wth a time delay different times)
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 11, numRuns = 10,
                  numInds = 440, numLoci = 1008,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.1008SNPs.440ind.2014.unix.structure", 
                  outPrefix = "px1008.2014",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub",
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)


# 2015, 402ind
# 20/2/2017: run with more burnins and mcmc reps
# Edit: 8/10/2017: repeat, this time 20 runs using a different random seed each time (using OrderByRun = TRUE, and pasting commands in for each run wth a time delay different times)
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 11, numRuns = 10,
                  numInds = 402, numLoci = 1008,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.1008SNPs.402ind.2015.unix.structure", 
                  outPrefix = "px1008.2015",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub",
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)
# edt:31/1/2017. Just need to run the commands for 2015, then queue the files.

# *******************************************
# Run structure for Paus.5pops.53ind, K = 1-7
# 12/2/2017
# *******************************************

sedStructureFiles(minK = 1, maxK = 7, startRunNum = 1, numRuns = 10,
                  numInds = 53, numLoci = 978,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Paus.5pops.53ind.978SNPs.unix.structure", 
                  outPrefix = "pa978",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub",
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/")
# edt:31/1/2017. Just need to run the commands for 2015, then queue the

#*********************************************************************************
# Run structure for pxpa.5pops.101ind, K = 1-10, 500000 burnin, 10^6 reps, 10 runs
# 13/2/2017
# ********************************************************************************
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 10,
                  numInds = 101, numLoci = 708,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "pxpa.5pops.101ind.708SNPs.unix.structure", 
                  outPrefix = "pxpa708",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub",
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/")

# #######################################################################################
# Run structure analysis for Plutella australiana populations, excluding sample `calca-7`
# 17/6/2017
# 19/6/2017: at 10 runs, two vlues of k had zero standard errors, precluding evanno.
# ... therefore, perform another 10 runs for each dataset.
# 20/6/2017: I accidentally did 20 instead of 10 extra runs for each K, making 30 in total
# #######################################################################################
# pxpa, 100ind, 707 SNPs
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 10,
                  numInds = 100, numLoci = 707,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "pxpa.5pops.100ind.707SNPs.unix.structure", 
                  outPrefix = "pxpa707",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub",
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)

# Paus, 52ind, 976 SNPs
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 10,
                  numInds = 52, numLoci = 976,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Paus.5pops.52ind.976SNPs.unix.structure", 
                  outPrefix = "pa976",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub",
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)



# ###########################################################################################
# 25/9/2017
# Run structure for the Px 2015 dataset, excluding SouthEnd samples
# Px, 386 ind, 27 pops
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 20,
                  numInds = 386, numLoci = 1008,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.27pops.386ind.2015.unix.structure", 
                  outPrefix = "px386",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub",
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)



# 26/9/2017
# Run structure for the Px 2014 dataset, excluding 3 divergent/unusual Esperance samples
# Px, 437 ind, 31 pops
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 20,
                  numInds = 437, numLoci = 1008,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.31pops.437ind.2014.unix.structure", 
                  outPrefix = "px437",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub",
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)


# 10/10/2017
# Run STRUCTURE for all P. xylostella samples, 59 pops, 842 ind.
sedStructureFiles(minK = 1, maxK = 15, startRunNum = 1, numRuns = 20,
                  numInds = 842, numLoci = 1008,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.59pops.842ind.1008SNPs.unix.structure", 
                  outPrefix = "px1008.842ind",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)

# Now run STRUCTURE for all P. xylostella samples EXCLUDING Southend: 58 pops, 826 ind.
sedStructureFiles(minK = 1, maxK = 15, startRunNum = 1, numRuns = 20,
                  numInds = 826, numLoci = 1008,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.58pops.826ind.1008SNPs.rmSouthend.unix.structure", 
                  outPrefix = "px1008.826ind",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)

# ################
# Run STRUCTURE for P. xylostella samples from resampled populations
# 15/10/2017

# first, the 6 resampled locations from vegetables (12 pops) 
sedStructureFiles(minK = 1, maxK = 7, startRunNum = 1, numRuns = 20,
                  numInds = 159, numLoci = 1008,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.12pops.159ind.1008SNPs.unix.structure", 
                  outPrefix = "px159",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)

# now, the 9 resampled locations < 10km distance (18 pops) 
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 20,
                  numInds = 248, numLoci = 1008,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.18pops.248ind.1008SNPs.unix.structure", 
                  outPrefix = "px248",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)

#########################################################
# Run structure for the final datasets for the Px study #
# 22/12/2017                                            #
######################################################### 

# The final datasets are:
# `Px.31pops.434ind.2014.maxMiss0.7.1009SNPs.structure`
# `Px.28pops.399ind.2015.maxMiss0.7.1009SNPs.structure`
# `Px.31pops.434ind.2014.maxMiss0.8.859SNPs.structure`
# `Px.28pops.399ind.2015.maxMiss0.8.859SNPs.structure`

# ------------
# max miss 0.7
# ------------
# `Px.31pops.434ind.2014.maxMiss0.7.1009SNPs.structure`
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 15,
                  numInds = 434, numLoci = 1009,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.31pops.434ind.2014.maxMiss0.7.1009SNPs.unix.structure", 
                  outPrefix = "px434.1009SNPs",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)

# `Px.28pops.399ind.2015.maxMiss0.7.1009SNPs.structure`
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 15,
                  numInds = 399, numLoci = 1009,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.28pops.399ind.2015.maxMiss0.7.1009SNPs.unix.structure", 
                  outPrefix = "px399.1009SNPs",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)

# -------------------
# max missing 0.8:
# -------------------
# `Px.31pops.434ind.2014.maxMiss0.8.859SNPs.structure`
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 15,
                  numInds = 434, numLoci = 859,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.31pops.434ind.2014.maxMiss0.8.859SNPs.unix.structure", 
                  outPrefix = "px434.859SNPs",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)


# `Px.28pops.399ind.2015.maxMiss0.8.859SNPs.structure`
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 15,
                  numInds = 399, numLoci = 859,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.28pops.399ind.2015.maxMiss0.8.859SNPs.unix.structure", 
                  outPrefix = "px399.859SNPs",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)


# --------------------------
# Px 3 pops, 42 individuals
# --------------------------

# `Px.3pops.42ind.1554SNPs.structure` # esperance, southend and dalby QLD.
sedStructureFiles(minK = 7, maxK = 8, startRunNum = 1, numRuns = 10,
                  numInds = 42, numLoci = 1554,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.3pops.42ind.1554SNPs.unix.structure", 
                  outPrefix = "px42ind.1554SNPs",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)

# --------------------------------------
# now max miss 0.7, maf 0.01, 1225 SNPs
# --------------------------------------

# 2014: `Px.31pops.434ind.2014.maxMiss0.7.maf0.01.1225SNPs.unix.structure`
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 15,
                  numInds = 434, numLoci = 1225,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.31pops.434ind.2014.maxMiss0.7.maf0.01.1225SNPs.unix.structure", 
                  outPrefix = "px434ind.1225SNPs",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)

# 2015: `Px.28pops.399ind.2015.maxMiss0.7.maf0.01.1225SNPs.unix.structure`
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 15,
                  numInds = 399, numLoci = 1225,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.28pops.399ind.2015.maxMiss0.7.maf0.01.1225SNPs.unix.structure", 
                  outPrefix = "px399ind.1225SNPs",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)

# 2015 excluding southend: `Px.27pops.383ind.2015.maxMiss0.7.maf0.01.1225SNPs.unix.structure`
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 15,
                  numInds = 383, numLoci = 1225,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.27pops.383ind.2015.maxMiss0.7.maf0.01.1225SNPs.unix.structure", 
                  outPrefix = "px383.1225SNPs",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)


# --------------------------------------------------
# now maxmiss 0.8, maf 0.01, 1032 SNPs (31/12/2017)
# --------------------------------------------------

# px 31 pops, 434ind, 2014, 1032 SNPs
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 15,
                  numInds = 434, numLoci = 1032,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.31pops.434ind.2014.maxMiss0.8.maf0.01.1032SNPs.unix.structure", 
                  outPrefix = "px434.1032SNPs",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)

# px 28 pops, 399ind, 2015, 1032 SNPs
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 15,
                  numInds = 399, numLoci = 1032,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.28pops.399ind.2015.maxMiss0.8.maf0.01.1032SNPs.unix.structure", 
                  outPrefix = "px399.1032SNPs",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)

# px 27 pops, 383ind, 2015, 1032 SNPs, exluding southend population
sedStructureFiles(minK = 1, maxK = 10, startRunNum = 1, numRuns = 15,
                  numInds = 383, numLoci = 1032,
                  burnins = 500000, mcmc = 1000000, 
                  inFile = "Px.27pops.383ind.2015.maxMiss0.8.maf0.01.1032SNPs.unix.structure", 
                  outPrefix = "px383.1032SNPs",
                  mpTemplate = "mainparams_template", 
                  extraparams = "ps1.extraparams",
                  subTemplate = "structure_template.sub", 
                  locOutDir = "C:/UserData/Kym/PhD/RAD2/structure/",
                  orderByRun = TRUE)



# End script
# ####################################






