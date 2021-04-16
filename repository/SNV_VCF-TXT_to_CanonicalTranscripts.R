library(data.table)
args <- commandArgs(TRUE)

mutfile = args[1]
species = args[2]

#species <- "Human"
#mutfile <- "/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD1/hMANEC_combined/1N7T85/results/Mutect2/1N7T85.Tumor.Mutect2.txt"

if(species == "Human"){
  transcriptTerm <- "ENST"
} else if(species == "Mouse"){
  transcriptTerm <- "ENSMUST"
} else {
  stop(paste0("Error: invalid species name used! Valid: Human or Mouse. Found: ", species))
}

mutDT <- fread(mutfile)
featureDT <- mutDT[, "ANN[*].FEATUREID"]
names(featureDT) <- "transcript"

# clean up the list and only keep valid transcript IDs
featureDT <- featureDT[transcript != ""]
featureDT <- featureDT[grep(transcriptTerm, transcript)]

# remove the version number
featureDT[, transcript := gsub("\\..*$","", transcript)]

outfile <- paste0(dirname(mutfile), "/", gsub(".txt", ".canonicalTranscripts.txt", basename(mutfile)))

fwrite(featureDT, outfile, col.names = F, quote = F)