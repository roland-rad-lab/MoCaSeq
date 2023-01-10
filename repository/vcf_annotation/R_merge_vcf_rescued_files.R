# GOAL  : Combine VCFrescued files
# USAGE : Rscript scripts/R_merge_vcf_rescued_files.R path/to/rescued/directory path/to/output/directory
# SOURCE: Modified Niklas's original script just to add command line arguments

# Load libraries
suppressPackageStartupMessages(library('rapportools'))
suppressPackageStartupMessages(library('data.table'))
suppressPackageStartupMessages(library('stringr'))


# Command Arguments
warnings         <- warnings();
args             <- commandArgs(TRUE);
inputDir         <- args[1] # Input folder to search all the vcf files 
outputDir        <- args[2]
callerBase       <- args[3]
callerRescue     <- args[4]

if(length(callerBase) < 1){
  callerBase    <- callerBase
}else{
  callerBase    <- 'Mutect2'
}

if(length(callerRescue) < 1){
  callerRescue    <- callerRescue
}else{
  callerRescue    <- 'Mutect2'
}

main <- function() {


  # Create output dir if not exists
  system(paste("mkdir -p ", outputDir, sep=''))

  # 2. provide a list of all rescued folders (with full path)
  # 2a) automatically find the folders on the server
  # allfolders <- findFolders(inputDir)

  # 2b) provide them manually
  #allfolders <- c("/home/rad/Downloads/rescued1/","/home/rad/Downloads/rescued2/")
  allfolders <- c(paste0(inputDir, "/results/rescued"))
  
  # 3. start  
  startRun(allfolders)
}

findFolders <- function(serverpath){
  samplenames <- list.dirs(serverpath, recursive = F, full.names = F)
  allfolders <- c()
  for(sample in samplenames){
    somefolder <- paste(serverpath,sample,"results/rescued/",sep='/')
    # only keep of the subfolder results/rescued/ is found
    if(dir.exists(somefolder)){
      allfolders <- c(allfolders, somefolder)
    }
  }
  return(allfolders)
}

startRun <- function(allfolders){
  # loop over every input folder
  for(filepath in allfolders){

    # get all files in that folder
    allfiles <- list.files(filepath, full.names = T, pattern = ".txt")
    # print(allfiles)
    
    if(length(allfiles) == 0){
      stop("Rescue folder does not exist or is empty!")
    }
    
    # loop over all files in the given folder
    outDT <- data.table() # init data
    for(afile in allfiles){
      
      #afile <- allfiles[2]
      #afile <- "/run/user/1000/gvfs/smb-share:server=imostorage.med.tum.de,share=fastq/Studies/AGRad_AAVLiver/AAVL1.1iLL/results/rescued//AAVL1.1iLL.Mutect2.AAVL1.1iLL.Mutect2.txt"
      
      # get the filename
      fileTag <- basename(afile)
      
      # check if there are . in the filename, if yes replace it
      splitregex <- "^(.*)\\.(Mutect2|Strelka)\\.(.*)\\.(Mutect2|Strelka)(.*)\\.txt"
      base <- gsub(splitregex, "\\1", fileTag)
      callerBase <- gsub(splitregex, "\\2", fileTag)
      rescue <- gsub(splitregex, "\\3", fileTag)
      callerRescue <- gsub(splitregex, "\\4", fileTag)
      
      fileEnding <- gsub(splitregex, "\\5", fileTag)
      # with this regex there will be a . before the "unique" in this case we remove it 
      if(fileEnding == ".unique"){
        fileEnding <- "unique"
      }

      # split the input filename to get the variables determining e.g. rescued by what with which caller
      # OLD
      #fileTag <- strsplit(fileTag, "\\.")[[1]]
      #base <- fileTag[1]
      #callerBase <- fileTag[2]
      #rescue <- fileTag[3]
      #callerRescue <- fileTag[4]
      #fileEnding <- fileTag[5]
      
      # read the file
      DT <- fread(afile)
      
      # the classic file (base) will always be base=rescue with the same caller
      if(base == rescue & callerBase == callerRescue){
        
        # ignore the empty .unique.txt file
        if(fileEnding != "unique"){
          DT[, RESCUED := "no"]
        } else {
          next()
        }
        
      } else {
        # these are the rescue files (first part of filename is not equal to last part)
        # e.g. DS21_8661_LivMet-1.Mutect2.DS21_8661_LungMet-1.Mutect2.unique.txt
        # here we only want the "unique" files
        
        if(fileEnding == "unique"){
          rescuedBy <- paste0(rescue,"-", callerRescue)
          DT[, RESCUED := rescuedBy]
        } else {
          next()
        }
      }
      
      # bind to the final output data
      outDT <- rbind(outDT, DT)
    }

    # print(outDT)
    # optional filtering
    # outDT <- outDT[`GEN[Normal].AD` == 0]
    
    # formating
    outDT[, ID := paste(CHROM, POS, REF, ALT, sep = "_")]
    outDT[, SearchID := paste(CHROM, POS, sep = ":")]
    outDT[, fileID := base]
    outDT[, totalCov := `GEN[Tumor].RD` + `GEN[Tumor].AD`]
    outDT[, VAF := `GEN[Tumor].AD` / `totalCov`]
    setcolorder(outDT, c("fileID", "ID","SearchID", "CHROM", "POS", "REF", "ALT", "GEN[Tumor].VF","VAF", "totalCov"))

    outDT[, `ANN[*].IMPACT` := NULL]
    
    #c(colnames(outDT[, -"ANN[*].FEATUREID"]))
    collapcols <- c("ANN[*].GENE", "ANN[*].EFFECT", "ANN[*].FEATUREID", "ANN[*].FEATUREID", "ANN[*].HGVS_C", "ANN[*].HGVS_P")

    allothercols <- colnames(outDT[, -collapcols, with=F])
    outDT <- outDT[, .(paste0(`ANN[*].GENE`, collapse = ","),
                      paste0(`ANN[*].EFFECT`, collapse = ","), 
                      paste0(`ANN[*].FEATUREID`, collapse = ","),
                      paste0(`ANN[*].HGVS_C`, collapse = ","),
                      paste0(`ANN[*].HGVS_P`, collapse = ",")), by=allothercols]
    setnames(outDT, "V1", "GENE")
    setnames(outDT, "V2", "EFFECT")
    setnames(outDT, "V3", "FEATUREID")
    setnames(outDT, "V4", "HGVS_C")
    setnames(outDT, "V5", "HGVS_P")
    
    # collapse duplicated rescue rows (the same mutation can be rescued by different samples)
    allothercols <- colnames(outDT[, -c("RESCUED")])
    outDT <- outDT[, paste0(RESCUED, collapse = ","), by=c(allothercols)]
    setnames(outDT, "V1", "RESCUED") 
    
    # Select and rename the relevant columns
    outDT <- outDT[, .(FileID=fileID, VariantID=ID, SearchID, Chrom=CHROM, GenomicPos=POS, Ref=REF, Alt=ALT, 
              TumorVAFMutect2=`GEN[Tumor].VF`, TumorVAFCustom=VAF, TumorTotalCov=totalCov, 
              TumorRefCov=`GEN[Tumor].RD`, TumorAltCov=`GEN[Tumor].AD`, 
              NormalRefCov=`GEN[Normal].RD`, NormalAltCov=`GEN[Normal].AD`,
              GeneName=GENE, VariantLocationInGene=EFFECT, GeneID=FEATUREID, HGVSc=HGVS_C, HGVSp=HGVS_P, Rescued=RESCUED)]
    
    # Rename the header
    #setnames(outDT, c("FileID", "VariantID","SearchID", "Chrom", "GenomicPos", "Ref", "Alt", "TumorVAFMutect2", "TumorVAFCustom", "TumorTotalCov", "TumorRefCov", "TumorAltCov", "NormalRefCov", "NormalAltCov", 
    #                  "GeneName", "VariantLocationInGene", "GeneID", "HGVSc", "HGVSp", "Rescued"))

    # combine the given output folder with a new results filename
    outFile <- paste0(outputDir, "/", base, "_", callerBase, ".txt")
    
    # print results
    fwrite(outDT, outFile, sep = "\t")

    cat(paste0("\t- ",outFile,"\n"))
  }
}


# Print session info as log file formatted in tabular format
# Source: https://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
print_session_info <- function(session_logfile){
	suppressPackageStartupMessages(library("devtools"))
	suppressPackageStartupMessages(library("knitr"))

	# Get all the session info to the variable
	my_session_info <- devtools::session_info()

	# Print it in the tabular format using knitr
	writeLines(text = {
	    paste(sep = "\n", collapse = "",
	          paste0(rep("-", 80), collapse = ""),
	          paste(paste0(rep("-", 32), collapse = ""),
	                "R environment",
	                paste0(rep("-", 33), collapse = "")),
	          paste0(rep("-", 80), collapse = ""),
	          paste(knitr::kable(data.frame(setting = names(my_session_info$platform),
	                                  value = as.character(my_session_info$platform))), collapse = "\n"),
	          paste0(rep("-", 80), collapse = ""),
	          paste(paste0(rep("-", 35), collapse = ""),
	                "packages",
	                paste0(rep("-", 35), collapse = "")),
	          paste0(rep("-", 80), collapse = ""),
	          paste(knitr::kable(my_session_info$packages), collapse = "\n")
	          )
	}, con = session_logfile)
}

# Run the main function
main()





