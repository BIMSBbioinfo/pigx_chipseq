# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Count_chrM')

# ---------------------------------------------------------------------------- #


Count_chrM <- function(
    infiles, 
    bowtie2log, 
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
    
    ## load bowtie2 log
    cnts.stat = readRDS(bowtie2log)
    
    res <- list()
    for(i in 1:length(infiles)){
      
      infile = infiles[i]
      name = basename(dirname(infile))
      message(name)
      # subset for sample
      cnts.stat.sub = cnts.stat[cnts.stat$sample_name == name,]
      print(cnts.stat.sub)
      mapped.total = cnts.stat.sub$value[cnts.stat.sub$stat == 'mapped.total']
      count = tryCatch(Rsamtools::countBam(file = infile,
                                  param = Rsamtools::ScanBamParam(which = gr))$records,
                       error=function(e) 0)
      stat = data.frame("value"=c(count,count/mapped.total*100), 
                        "stat"=c("reads.chrM","perc.chrM"), 
                        sample_name = name,
                        row.names = NULL)
      res[[name]] = stat
      
    }
    resdf <- do.call("rbind",res)
    saveRDS(resdf, file = outfile)
  
}

# ---------------------------------------------------------------------------- #
# function call
Count_chrM(
    infiles        = argv$input[['bamfiles']],
    bowtie2log     = argv$input[['bowtie2log']],
    genome_chrlen  = argv$params[['chrlen']],
    chrM_givenName = argv$params[['chrM_givenName']],
    outfile        = argv$output[['outfile']]
)


