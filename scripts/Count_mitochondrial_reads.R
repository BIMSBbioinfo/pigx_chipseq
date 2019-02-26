# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Count_chrM')

# ---------------------------------------------------------------------------- #


Count_chrM <- function(
    infiles, 
    genome_chrlen, 
    outfile, 
    chrM_givenName=NULL
) {
    
    library(Rsamtools)
    library(GenomeInfoDb)
    library(stringr)
    
    # ------------------------------------------------------------------------ #
    
    ## define some default names for chrM
    chrM_defaultName <-  c("chrM","MT","MtDNA","chMT")
    chrM_defaultName <- unique(c(chrM_defaultName,chrM_givenName))
    
    infiles <- gsub(pattern = '.bai',replacement = '',infiles)
    
    ## read in genome chrom lengths
    genome <- read.delim(genome_chrlen,header = FALSE,stringsAsFactors = FALSE)
    names(genome) <- c("seqnames","seqlengths")
    
    ## convert to seqinfo 
    genome <- as(genome,"Seqinfo")
    chr_names <- seqnames(genome)
    
    ## check if chrM is included
    if (!any(seqnames(genome) %in% chrM_defaultName)) {
      warning("No chrM detected, returning 0 counts.")
      gr <- GRanges(seqnames = "chrM",ranges = IRanges(start = 1,end = 1e6))
    } else {
      message(chr_names[chr_names %in% chrM_defaultName]," detected.")
      gr <- as(genome[chr_names[chr_names %in% chrM_defaultName]],"GRanges")
    }
    
    res <- list()
    for(i in 1:length(infiles)){
      
      infile = infiles[i]
      name = str_replace(basename(infile),'.sorted.bam','')
      message(name)
      count = tryCatch(Rsamtools::countBam(file = infile,
                                  param = Rsamtools::ScanBamParam(which = gr))$records,
                       error=function(e) 0)
      stat = data.frame("value"=count, "stat"="reads.chrM", sample_name = name,row.names = NULL)
      res[[name]] = stat
      
    }
    resdf <- do.call("rbind",res)
    saveRDS(resdf, file = outfile)
  
}

# ---------------------------------------------------------------------------- #
# function call
Count_chrM(
    infiles     = argv$input[['bamfiles']],
    genome_chrlen = argv$params[['chrlen']],
    chrM_givenName = argv$params[['chrM_givenName']],
    outfile     = argv$output[['outfile']]
)


