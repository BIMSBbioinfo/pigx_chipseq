# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extend_Regions', parseFlags = TRUE)

# ---------------------------------------------------------------------------- #
Extend_Regions = function(
    input_bamfile,
    input_stats,
    outpath,
    extend          = NULL,
    bam_name        = NULL,
    scale_index     = FALSE,
    scale_factor    = 1e6,
    chunk_size      = 1e7,
    library_type    = 'single',
    spikein         = 'No',
    deduplicate     = FALSE,
    minmapq         = 0
){

    # ------------------------------------------------------------------------ #
    suppressPackageStartupMessages({
        library(GenomicRanges)
        library(GenomicAlignments)
        library(rtracklayer)
        library(Rsamtools)
    })

    # ------------------------------------------------------------------------ #
    if(!is.numeric(extend))
        stop('Extend_Regions: extend parameter needs to be a number')

    if(!file.exists(input_bamfile))
        stop('Extend_Regions: input file does not exist')

    # ------------------------------------------------------------------------ #
    param   = ScanBamParam(flag=scanBamFlag(isPaired = (library_type == 'paired'),
                                            isProperPair = library_type == 'paired',
                                            isDuplicate = deduplicate),
                           mapqFilter = minmapq)

    # ------------------------------------------------------------------------ #
    bamfile = BamFile(input_bamfile, yieldSize=chunk_size)
    open(bamfile)
    lcov = NULL
    total = 0
    while (length(chunk <- readGAlignments(bamfile, param=param))) {
        gchunk = trim(resize(granges(chunk), width=extend))
        cov = coverage(gchunk)
        if(is.null(lcov)){
            lcov = cov
        }else{
            lcov = lcov + cov
        }
        total = total+length(gchunk)
    }
    close(bamfile)

    # ------------------------------------------------------------------------ #
    # checks whether to scale the bigWig
    if(scale_index == TRUE || toupper(scale_index) == 'YES'){

      # checks for spike-in normalization
      if(toupper(spikein) == 'YES'){
        stats        = readRDS(input_stats)
        stats        = subset(stats, sample_name == bam_name)
        mapped_main  = subset(stats, genome_type=='Main' & stat == 'mapped.total')$value
        mapped_spike = subset(stats, genome_type=='Spike-in' & stat == 'mapped.total')$value
        scale_factor = mapped_main / mapped_spike
        lcov = round(lcov * scale_factor, 3)

      # normalization to total number of counts
      }else{
        lcov = round(lcov * (scale_factor/total),2)
      }
    }

    # ------------------------------------------------------------------------ #
    export.bw(lcov, outpath)
}

# ---------------------------------------------------------------------------- #
Extend_Regions(
  input_bamfile   = argv$input[['bamfile']],
  input_stats     = argv$input[['stats']],
  outpath         = argv$output[['outfile']],
  extend          = argv$params[['extension']],
  bam_name        = argv$params[['sample_name']],
  scale_index     = argv$params[['scale']],
  library_type    = argv$params[['library_type']],
  spikein         = argv$params[['spikein']],
  minmapq         = argv$params[['minmapq']],
  deduplicate     = argv$params[['deduplicate']]
)
