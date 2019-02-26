# ----------------------------------------------------------------------------- #
rule count_mitochondrial_reads:
    input:
        bamfiles = expand(os.path.join(PATH_MAPPED, "{name}", "{name}" + BAM_SUFFIX + ".bai"), name=NAMES)
    output:
        outfile  = os.path.join(PATH_RDS, "Count_chrM.rds")
    params:
        chrlen = rules.index_to_chrlen.output.outfile,
        chrM_givenName = PARAMS['count_mitochondrial_reads']['chrM_givenName'],
        threads   = 1,
        mem       = '1G',
        Rscript   = SOFTWARE['Rscript']['executable']
    log:
       logfile = os.path.join(PATH_LOG, 'count_chrM.log')
    message:
        """
            Running: count_mitochondrial_reads:
                output:  {output.outfile}
        """
    run:
        RunRscript(input, output, params, log.logfile, 'Count_mitochondrial_reads.R')
