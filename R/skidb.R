#############################################################################r
## Marcin Imielinski
## The Broad Institute of MIT and Harvard / Cancer program.
## marcin@broadinstitute.org
##
## Weill-Cornell Medical College
## mai9037@med.cornell.edu
##
## New York Genome Center
## mimielinski@nygenome.org
##
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

#' @import tools
#' @import Biostrings
#' @import parallel
#' @import GenomicRanges
#' @import gUtils
#' @importFrom data.table data.table setkey setkeyv := setnames is.data.table as.data.table


#' @name skidb_env
#' @title skidb_env
#'
#' @description
#'
#' Get environment variables of all flat file paths associated with skidb, or specific file path.  T
#'
#' This can be used to
#' populate the file pointers with the appropriate files or rewire the paths using environment variables
#' having the same name as the skidb_env environment variable.  skidb_env is used by other skidb functions
#' to grab files and paths to software.
#'
#'
#' gets all the individual / pair and sample info for an individual / pair set
#' @author Marcin Imielinski
#' @export
skidb_env = function(query = NULL)
    {
        DB_ROOT = Sys.getenv('SKI_DB_ROOT')
        if (nchar(DB_ROOT)==0)
            {
                warning('SKI_DB_ROOT not set, so using ~/DB as a default root path')
                DB_ROOT = normalizePath('~/DB/')
            }
        GAF.DIR = paste(DB_ROOT, '/GAF/', sep = "/")
        SOFTWARE.DIR = Sys.getenv('SKI_SOFTWARE_ROOT')
        if (nchar(SOFTWARE.DIR)==0)
            SOFTWARE.DIR = normalizePath('~/Software/')
        UCSC.DIR = paste(DB_ROOT, "UCSC", sep = "/")
        env = list(
            DB_ROOT = DB_ROOT,
            MSIGDB_PATH = paste(DB_ROOT, '/GSEA/msigdb.v5.0.symbols.gmt', sep = ""),
            UNIPROT.DAT = paste(DB_ROOT, '/Uniprot/uniprot_sprot_human.dat', sep = ""),
            UNIPROT.ROOT = paste(DB_ROOT, '/Uniprot/', sep = ""),
            UNIPROT.FT = paste(DB_ROOT, '/Uniprot/uniprot_sprot_human.ft.txt', sep = ""),
            UNIPROT.RS=  paste(DB_ROOT, '/Uniprot/uniprot_sprot_human.rs.txt', sep = ""),
            UNIPROT.FA = paste(DB_ROOT, '/Uniprot/HUMAN.fasta', sep = ""),
            UNIPROT.ACMAP = paste(DB_ROOT, '/Uniprot/uniprot_sprot_human.acmap.txt', sep = ""),
            UNIPROT.FEATUREKEYS = paste(DB_ROOT, '/Uniprot/uniprot_feature_keys.txt', sep = ""),
            PFAM.HUMAN = paste(DB_ROOT, '/Pfam/Pfam-A.full.human.lean', sep = ""),
            PFAM.HUMAN.RDS = paste(DB_ROOT, '/Pfam/Pfam-A.full.human.lean.rds', sep = ""),
            PFAM.HUMAN.META.RDS = paste(DB_ROOT, '/Pfam/Pfam-A.full.human.meta.rds', sep = ""),
            PFAM.FA = paste(DB_ROOT, '/Pfam/Pfam-A.fasta', sep = ""),
            GAF.DIR = paste(DB_ROOT, '/GAF/', sep = "/"),
            GAF.SOURCE = paste(GAF.DIR, 'TCGA.hg19.June2011.gaf', sep = "/"),
            CCDS.MAP.FILE = paste(DB_ROOT, '/NCBI/CCDS2Sequence.current.txt', sep = ''),
            CCDS.FA.FILE = paste(DB_ROOT, '/NCBI/CCDS_nucleotide.current.fna', sep = ''),
            HG19.FFT = paste0(DB_ROOT, '/ffTracks/hg19.rds'),
            HG19.CONTEXT.FFT = paste0(DB_ROOT, '/ffTracks/hg19.context64.rds'),
            REFGENE.FILE.HG19 = paste(DB_ROOT, '/UCSC/refGene.hg19.txt', sep = ""),
            REFGENE.FILE.HG19.GR =  paste(DB_ROOT, '/UCSC/refGene.hg19.gr.rds', sep = ""),
            REFGENE.FILE.HG19.GRL =  paste(DB_ROOT, '/UCSC/refGene.hg19.grl.rds', sep = ""),
            GENCODE.FILE.HG19.GR = paste(DB_ROOT, '/GENCODE/gencode.v19.annotation.gtf.gr.rds', sep = ''),
            REFGENE.FILE.HG19.CDS.GRL =  paste(DB_ROOT, '/UCSC/refGene.hg19.cds.grl.rds', sep = ""),
            REFGENE.FILE.HG18 = paste(DB_ROOT, "/UCSC/refGene.txt", sep = ""),
            SNP.MARKER.MAP = paste(DB_ROOT, "/SNP.ANNOT/Affymetrix.hg19.NA31.markers.data.table.rds", sep = ""),
            ENSEMBL.GENE = paste0(DB_ROOT, '/Ensembl/Ensembl64.gene.txt'),
            ENSEMBL.TX = paste0(DB_ROOT, '/Ensembl/Ensembl64.transcriptome.map'),
            HG.SIZE = 2897310462,
            KG.FILE.HG19 = paste(DB_ROOT, '/UCSC/knownGene.hg19.txt', sep = ""),
            KG.FILE.HG18 = paste(DB_ROOT, '/UCSC/knownGene.hg18.txt', sep = ""),
            KG2RG.FILE.HG19 = paste(DB_ROOT, '/UCSC/knownToRefSeq.hg19.txt', sep = ""),
            KG2RG.FILE.HG18 = paste(DB_ROOT, '/UCSC/knownToRefSeq.txt', sep = ""),
            GENOMIC.FEATURES.ROOT = paste(DB_ROOT, '/GenomicFeatures/', sep = "/"),
            UCSC.CHROM.SIZES.PATH = paste(DB_ROOT, "/UCSC/hg19.ucsc.chrom.sizes.txt", sep = ""),
            UCSC.DIR = UCSC.DIR,
            REPEATMASKER.HG19 = paste(UCSC.DIR, 'hg19.repeatmask.txt', sep = '/'),
            REPEATMASKER.HG19.GR = paste(UCSC.DIR, 'hg19.repeatmask.rds', sep = '/'),
            ENTREZ.GENE.FN = paste(DB_ROOT, "/Entrez/Homo_sapiens.gene_info", sep = ""),
            SOFTWARE.DIR = SOFTWARE.DIR,
            BLAST.DIR = paste(SOFTWARE.DIR, '/ncbi-blast-2.2.26+/bin/', sep = ""),
            MUSCLE.PATH = paste(SOFTWARE.DIR, '/MUSCLE/MUSCLE', sep = ''),
            TMP.DIR = paste0(DB_ROOT, '/tmp/'),
            FH.ENTITY.URL = "http://firehose:8080/cga/ws/entity/getAnnotations/%s?entityNames=%s&filterSetType=Individual_Set&workspaceName=%s",
            FH.LIST.URL = "http://firehose:8080/cga/ws/list/%s?workspaceName=%s",
            FH.GET.ENTITIES.URLS = "http://firehose:8080/cga/ws/entity/getAnnotations/%s?workspaceName=%s&%s",
            FH.WKFLOW.URL = "http://firehose:8080/cga/ws/workflow/execute/Individual_Set/%s/%s?workspaceName=%s",
            FH.START.URL = "http://firehose:8080/cga/ws/pipeline/start/Individual_Set/%s/%s?workspaceName=%s",
            FH.IMPORT.URL = "http://firehose:8080/cga/ws/list/Individual_Set?importFromFile=&workspaceName=%s",
            FH.UNPW = '~/.fh.config',
            FH.GETCONTAINERS.URL = "http://firehose:8080/cga/ws/list/%s?getContainers=&%s&containerType=%s&workspaceName=%s",
            FISS = '/xchip/tcga/Tools/gdac/bin/fiss',
            FISS_TMP = paste0(DB_ROOT, '/tmp.fiss')
        )


        if (is.null(query))
            return(env)

        if (nchar(Sys.getenv(query))>0)
            return(Sys.getenv(query))
        else
            {
                if (!file.exists(env[[query]]))
                    warning(cat('File', env[[query]], 'does not exist!\n'))
                return(env[[query]])
            }
    }

##  ________  __    __        __    __    __      __  __
## |        \|  \  |  \      |  \  |  \  |  \    |  \|  \
## | $$$$$$$$| $$  | $$      | $$  | $$ _| $$_    \$$| $$
## | $$__    | $$__| $$      | $$  | $$|   $$ \  |  \| $$
## | $$  \   | $$    $$      | $$  | $$ \$$$$$$  | $$| $$
## | $$$$$   | $$$$$$$$      | $$  | $$  | $$ __ | $$| $$
## | $$      | $$  | $$      | $$__/ $$  | $$|  \| $$| $$
## | $$      | $$  | $$       \$$    $$   \$$  $$| $$| $$
##  \$$       \$$   \$$        \$$$$$$     \$$$$  \$$ \$$


##############
#' @name fiss_get
#' @title fiss_get
#'
#' @description
#' replacement for get fh
#'
#' gets all the individual / pair and sample info for an individual / pair set
#' @author Marcin Imielinski
#' @export
##############
fiss_get = function(set = NULL,
    space = Sys.getenv('FIREHOSE_WORKSPACE'),
    wkspace = space, type = 'Pair', # can also be "individual", "Sample"
    sample_type_col = "sample_type", ## can be changed eg to "type"
    fuse = TRUE, ## if TRUE then will return FUSE path instead
    verbose = F, get.containers = T, mc.cores = 1, domain = 'cga')
{
  if (!(domain %in% c('cga', 'tcga-gdac', 'gtex')))
    stop('domain not recognized')

  path_type = ifelse(fuse, 'FusePath', 'SystemPath')

  FISS = paste(skidb_env('FISS'), '-d ', domain)
  if (grepl(' ', wkspace))
      {
          warning('Workspace contains spaces .. trimming')
          wkspace = gsub(' .*', '', wkspace)
      }

  if (grepl('(pair)|(pset)', type, ignore.case = T))
    {
      if (verbose)
        cat('Getting pairs\n')

      if (domain == 'tcga-gdac')
          str = paste(FISS, ' annot_get ', wkspace, ' pair pset=', set, sep = '')
      else
          str = paste(FISS, ' annot_get ', " -p ", path_type, " ", wkspace, ' pair pset=', set, sep = '')

      if (verbose)
        cat('Running:', str, '\n')

      p = pipe(str);

      pairs= read.delim(p, strings = F)
      rownames(pairs) = pairs$pair_id;

      if (any(grepl('There\\.are\\.validation\\.errors', colnames(pairs))))
        stop('Pair set / workspace combination not found')

      if (verbose)
        cat('Getting samples\n')

      if (domain == 'tcga-gdac')
          str = paste(FISS,  ' annot_get ', wkspace, ' sample pset=', set, sep = '')
      else
          str = paste(FISS,  ' annot_get ', " -p ", path_type, " ", wkspace, ' sample pset=', set, sep = '')

      if (verbose)
        cat('Running:', str, '\n')

      p = pipe(str);
      samples = read.delim(p, strings = F)

      if (verbose)
        cat('Matching up pairs and samples\n')

      sample_fields = setdiff(names(samples), c('individual_id'))

      pairs[, paste('Tumor', sample_fields, sep = "_")] = samples[match(pairs$case_sample, samples$sample_id), sample_fields]
      pairs[, paste('Normal', sample_fields, sep = "_")] = samples[match(pairs$control_sample, samples$sample_id), sample_fields]

      if (get.containers)
        {
          if (verbose)
            cat('Getting containers\n')
          containers = fiss_getcontainers(samples$sample_id, 'sample', wkspace, verbose = verbose, mc.cores = mc.cores, domain = domain)
          s2p = structure(containers$Individual, names = containers$Sample)
          s2p = s2p[intersect(names(s2p), c(pairs$case_sample, pairs$control_sample))]
          pairs$individual_id[((match(names(s2p), c(pairs$case_sample, pairs$control_sample))-1) %% nrow(pairs))+1] = s2p

          ## now fetch individual level annotations for all individuals in pair set
          if (domain == 'tcga-gdac')
              str = paste(FISS,  ' annot_get ', wkspace, ' indiv=', paste(pairs$individual_id, collapse = ','), sep = '')
          else
              str = paste(FISS,  ' annot_get ', " -p ", path_type, " ", wkspace, ' indiv=', paste(pairs$individual_id, collapse = ','), sep = '')

          if (verbose)

          p = pipe(str)

          CHUNKSIZE = 50
          ghetto.chunks = c(seq(1, nrow(pairs), CHUNKSIZE), nrow(pairs)+1); ## FH can't handle too large get container requests

          indiv_annots = do.call('rrbind', mclapply(ghetto.chunks[1:(length(ghetto.chunks)-1)], function(i)
            {
              if (verbose)
                cat('Getting indiv chunk ', i ,'\n')

              if (verbose)
                cat('Running:', str, '\n')

              if (domain == 'tcga-gdac')
                  str = paste(FISS, ' annot_get ',  wkspace, '  indiv=',
                      paste(pairs$individual_id[i:min(nrow(pairs), i+CHUNKSIZE-1)], collapse = ','), sep = '')
              else
                  str = paste(FISS, ' annot_get ', " -p ", path_type, " ", wkspace, '  indiv=',
                      paste(pairs$individual_id[i:min(nrow(pairs), i+CHUNKSIZE-1)], collapse = ','), sep = '')

              p = pipe(str)
              out = read.delim(p, strings = F, header = T)
              return(out);
            }, mc.cores = mc.cores))

          ix = match(indiv_annots$individual_id, pairs$individual_id)
          cnames = setdiff(names(indiv_annots), names(pairs))
          if (length(cnames)>0)
            {
              pairs[, cnames] = NA
              pairs[ix, cnames] = indiv_annots[, cnames]
            }

        }
      return(pairs)
    }
  else if (grepl('(indiv)|(iset)', type, ignore.case = T))  # individual
    {
      if (verbose)
        cat('Getting individuals\n')

      if (domain == 'tcga-gdac')
          str = paste(FISS,  ' annot_get ',wkspace, ' indiv iset=', set, sep = '')
      else
          str = paste(FISS,  ' annot_get ', " -p ", path_type, " ", wkspace, ' indiv iset=', set, sep = '')

      p = pipe(str)

      individuals= read.delim(p, strings = F)
      rownames(individuals) = individuals$individual_id;

      if (any(grepl('There\\.are\\.validation\\.errors', colnames(individuals))))
        stop('Individual set / workspace combination not found')

      if (verbose)
        cat('Getting samples\n')

      if (domain == 'tcga-gdac')
          p = pipe(paste(FISS,  ' annot_get ', wkspace, ' sample iset=', set, sep = ''))
      else
          p = pipe(paste(FISS,  ' annot_get ', " -p ", path_type, " ", wkspace, ' sample iset=', set, sep = ''))

      samples = read.delim(p, strings = F)

      if (get.containers)
        {
          if (verbose)
            cat('Getting containers\n')
          containers = fiss_getcontainers(samples$sample_id, 'sample', wkspace = wkspace, domain = domain, mc.cores = mc.cores)
          s2p = structure(containers$Individual, names = containers$Sample)
        }
      else # hack
        {
          if (!is.null(samples$Individual_Id))
            s2p = structure(samples$Individual_Id, names = samples$sample_id)
          else if (!is.null(samples$individual_id))
            s2p = structure(samples$individual_id, names = samples$sample_id)
          else
            s2p = structure(gsub('\\-Tumor|\\-Normal', '', samples$sample_id), names = samples$sample_id)
        }

      if (verbose)
        cat('Matching up individuals and samples\n')

      sample_fields = setdiff(names(samples), c('individual_id'))

      if ((sample_type_col %in% names(samples))) ## "type" is flex
          samples$sample_type = samples[, sample_type_col]
      else
          samples$sample_type = "Sample"

      ix = is.na(samples$sample_type) | nchar(samples$sample_type)==0
      if (any(ix))
          samples$sample_type[ix] = "Sample"

      ## now we want "sample_type" to aggressively dedup sample --> individual
      ## mappings, since we will create an extra set of columns for
      ## each unique value of sample_type, we want to avoid overwriting
      ## in case there are several samples of the same type for the same individuals
      samples$sample_type = sapply(strsplit(dedup(paste(s2p[samples$sample_id], samples$sample_type, sep = '\t'), suffix = ''), '\t'), function(x) x[2])

      sample.split = split(samples, samples$sample_type)

      for (stype in names(sample.split))
        {
          ix = match(s2p[sample.split[[stype]]$sample_id], individuals$individual_id);
          fields = paste(stype, sample_fields, sep = "_")
          individuals[, fields] = NA;
          individuals[ix[!is.na(ix)], fields] = sample.split[[stype]][!is.na(ix), sample_fields]
        }

      return(individuals)
  }
  else if (grepl('(samp)|(sset)', type, ignore.case = T))  # sample
      {

          if (domain == 'tcga-gdac')
              str = paste(FISS,  ' annot_get ', wkspace, ' sample sset=', set, sep = '')
          else
              str = paste(FISS,  ' annot_get ', " -p ", path_type, " ", wkspace, ' sample pset=', set, sep = '')

          if (verbose)
              cat('Running:', str, '\n')
          p = pipe(str)

          samples = read.delim(p, strings = F)

          if (any(grepl('There\\.are\\.validation\\.errors', colnames(samples))))
              if (domain == 'tcga-gdac')
                  stop('Sample set / workspace combination not found')
              else
                  {
                      str = paste(FISS,  ' annot_get ', " -p ", path_type, " ", wkspace, ' sample iset=', set, sep = '')
                      p = pipe(str)
                      samples = read.delim(p, strings = F)

                      if (any(grepl('There\\.are\\.validation\\.errors', colnames(samples))))
                          stop('Sample set / workspace combination not found')
                  }

          if (domain != 'tcga-gdac')
              if (nrow(samples)>0)
                  {
                      containers = tryCatch(fiss_getcontainers(samples$sample_id, 'sample', wkspace, verbose = verbose, mc.cores = mc.cores, domain = domain), error = function(e) data.frame())

                      if (nrow(containers)>0)
                          samples$individual_id =  structure(containers$Individual, names = containers$Sample)[samples$sample_id]

                      rownames(samples) = samples$sample_id
                  }

          return(samples)
      }
  else
    stop('entity type not recognized')
}


#' @name fiss_getcontainers
#' @title fiss_getcontainers
#'
#' @description
#'
#' low level function to get containers
#'
#'
#' @author Marcin Imielinski
#' @export
fiss_getcontainers = function(id, type = 'sample', wkspace = Sys.getenv('FIREHOSE_WORKSPACE'), verbose = F, mc.cores = 1, domain = 'cga')
    {

    if (!(domain %in% c('cga', 'tcga-gdac', 'gtex')))
      stop('domain not recognized')
    FISS = paste(skidb_env('FISS'), '-d ', domain)

    CHUNKSIZE = 50
    ghetto.chunks = c(seq(1, length(id), CHUNKSIZE), length(id)+1); ## FH can't handle too large get container requests

    containers = do.call('rbind', mclapply(ghetto.chunks[1:(length(ghetto.chunks)-1)], function(i)
      {
        if (verbose)
          cat('Getting container chunk ', i, '\n')
        p = pipe(paste(FISS, ' entity_get_containers ', wkspace, ' ', type, '=',
          paste(id[i:min(length(id), i+CHUNKSIZE-1)], collapse = ' '), sep = ''))
        out = read.delim(p, strings = F, header = F)
        return(out)
      }, mc.cores = mc.cores))

    # ugly loop to parse out the different container data which are outputted as a two column table

    if (nrow(containers)==0)
        return(containers)

    require(Hmisc)
#    header.ix = c(which(capitalize(containers[,1],un = T) %in% capitalize(type, un = T)), nrow(containers)+1)
    header.ix = c(which(skitools::capitalize(containers[,1]) %in% capitalize(type)), nrow(containers)+1)
    sub.tabs = cbind(header.ix[-length(header.ix)], header.ix[-1]-1, NA)[1:(length(header.ix)-1), ]

    if (is.null(dim(sub.tabs)))
      return(data.frame())

    sub.tabs[,3] = sub.tabs[,2]-sub.tabs[,1]+1

    sub.tab.id = unlist(lapply(1:nrow(sub.tabs), function(x) rep(x, sub.tabs[x, 3])))
    tmp.split = lapply(split(containers, sub.tab.id), function(x) {nm = x[1,]; y = x[-1, ]; colnames(y) = nm; y})

    tmp.split = tmp.split[sapply(tmp.split, nrow)>0]

    if (length(tmp.split)==0)
      return(data.frame())

    ## first combine tmp.split subtables with same column 2 name
    tmp.split.colnames = sapply(tmp.split, function(x) colnames(x)[2])
    unames = unique(tmp.split.colnames)

    tmp.split = lapply(unames, function(x) do.call('rbind', tmp.split[which(tmp.split.colnames %in% x)]))

    containers = tmp.split[[1]]
    if (length(tmp.split)>1)
      for (i in 2:length(tmp.split))
        containers = merge(containers, tmp.split[[i]], by = names(containers)[1], all = T)

    return(containers)
  }


#' @name fiss_delete
#' @title fiss_delete
#'
#' @description
#'
#' Delete entities via fiss
#'
#' @author Marcin Imielinski
#' @export
fiss_delete = function(
  entities = NULL,
  type = 'Individual', # can also be pair
  space = Sys.getenv('FIREHOSE_WORKSPACE'),
  wkspace = space,
  char = F,
  domain = 'cga')
  {
    if (!(domain %in% c('cga', 'tcga-gdac', 'gtex')))
      stop('domain not recognized')
    FISS = paste(skidb_env('FISS'), '-d ', domain)

    if (!char)
        {
            if (grepl('indiv', type, ignore.case = T))
                system(paste(FISS, ' indiv_delete ', wkspace, paste(entities, collapse = ' ')))
            else if (grepl('pair', type, ignore.case = T))
                system(paste(FISS, ' pair_delete ', wkspace, paste(entities, collapse = ' ')))
            else if (grepl('sample', type, ignore.case = T))
                system(paste(FISS, ' sample_delete ', wkspace, paste(entities, collapse = ' ')))
            else if (grepl('sset', type, ignore.case = T))
                system(paste(FISS, ' sset_delete ', wkspace, paste(entities, collapse = ' ')))
            else if (grepl('iset', type, ignore.case = T))
                system(paste(FISS, ' iset_delete ', wkspace, paste(entities, collapse = ' ')))
            else if (grepl('pset', type, ignore.case = T))
                system(paste(FISS, ' pset_delete ', wkspace, paste(entities, collapse = ' ')))
            else if (grepl('annot', type, ignore.case = T))
                system(paste(FISS, ' annot_delete ', wkspace, paste(entities, collapse = ' ')))
        }
    else
        {
            if (grepl('indiv', type, ignore.case = T))
                return(paste(FISS, ' indiv_delete ', wkspace, paste(entities, collapse = ' ')))
            else if (grepl('pair', type, ignore.case = T))
                return(paste(FISS, ' pair_delete ', wkspace, paste(entities, collapse = ' ')))
            else if (grepl('sample', type, ignore.case = T))
                return(paste(FISS, ' sample_delete ', wkspace, paste(entities, collapse = ' ')))
            else if (grepl('sset', type, ignore.case = T))
                return(paste(FISS, ' sset_delete ', wkspace, paste(entities, collapse = ' ')))
            else if (grepl('iset', type, ignore.case = T))
                return(paste(FISS, ' iset_delete ', wkspace, paste(entities, collapse = ' ')))
            else if (grepl('pset', type, ignore.case = T))
                return(paste(FISS, ' pset_delete ', wkspace, paste(entities, collapse = ' ')))
            else if (grepl('annot', type, ignore.case = T))
                system(paste(FISS, ' annot_delete ', wkspace, paste(entities, collapse = ' ')))
        }
  }


#' @name read.xml
#' @title read.xml
#' @description
#'
#' reads a bunch of xml in files into data.table, flattens, and
#' and concatenates the files via rbindlist fill
#'
#' @export
read.xml = function(files, pattern = NULL, as.data.table = TRUE, na.sort = TRUE, mc.cores = 1)
{
    out = rbindlist(mclapply(files, function(x)
        {
            row = unlist(XML::xmlToList(x))
            names(row) = gsub('\\.+', '.', names(row))

            if (!is.null(pattern))
                row = row[grep(pattern, names(row))]
            
            return(as.data.table(as.list(row)))
        }, mc.cores = mc.cores), fill = TRUE)

    if (na.sort)
        {
            ix = rev(as.numeric(names(sort(table(colSums(!is.na(row)))))))
            out = out[, ix, with = FALSE]
        }
    
    if (!as.data.table)
        out = as.data.frame(out)

    return(out)                       
}

#' @name fiss_task
#' @title fiss_task
#'
#' @description
#'
#' Run a task within firehose via fiss
#'
#' @author Marcin Imielinski
#' @export
fiss_task = function(
  name = NULL, ## (scalar) task name, if null will output a vector of task names in workspace
  entity = NULL, ## entity name
  type = 'pset', ## (scalar) entity type on which task operates (pset, iset, sset, indiv, sample, pair)
  start = F, ## if true will "start"
  stop = F, ## if true will "stop" task (overrides stop)
  export = FALSE, ## if true will export xml of task to a character vector
  verbose = F,
  domain = 'cga',
  space = Sys.getenv('FIREHOSE_WORKSPACE'),
  wkspace = space)
  {
    if (!(domain %in% c('cga', 'tcga-gdac', 'gtex')))
      stop('domain not recognized')
    FISS = paste(skidb_env('FISS'), '-d ', domain)

    if (is.null(name))
      {
        p = pipe(paste(FISS, ' task_list ', wkspace))
        out = readLines(p);
        close(p)
      }
    else if (export)
        {
            rfile = paste(skidb_env('FISS_TMP'), paste('tmp', runif(1), '.xml', sep = ''), sep = '/')
            str = sprintf('fiss task_export %s "%s" %s', space, gsub(" ", "%20", name), rfile)
            system(str)
            out = readLines(rfile)
            system(paste('rm', rfile))
            return(out)
        }
    else
      {
        if (is.null(entity))
          stop('entity must be specified')

        if (!(type %in% c('pset', 'iset', 'sset', 'indiv', 'sample', 'pair')))
          stop('type must be one of the following: ', paste(c('pset', 'iset', 'sset', 'indiv', 'sample', 'pair'), collapse = ', '))

        if (stop)
          cmd = 'task_stop'
        else if (start)
          cmd = 'task_start'
        else
          cmd = 'task_status'

        out = do.call('rbind', lapply(entity, function(x)
          {
            str = paste(FISS, ' ', cmd, ' ', wkspace, ' ', type, '=', x, ' ', '\\"', name, '\\"', sep = '')
            if (verbose)
              cat('Running: ', str, '\n')
            p = pipe(str)
            out = read.delim(p, strings = F, header = F)
            names(out) = c('entity', 'status', 'task', 'version', 'timestamp')
            return(out)
          }))
      }

    return(out)
  }


#' @name fiss_list
#' @title fiss_list
#'
#' @description
#'
#' List all entities of a given type in a given workspace.
#'
#' @author Marcin Imielinski
#' @export
fiss_list = function(
  type = 'pset', ## (scalar) entity type on which task operates (pset, iset, sset, indiv, sample, pair)
  set = NULL,  ## set name, eg sample set, pair set, individual set (only relevant for type = indiv, sample, or pair)
  space = Sys.getenv('FIREHOSE_WORKSPACE'),
  wkspace = space,
  verbose = FALSE,
  domain = 'cga'
  )
  {
    if (!(domain %in% c('cga', 'tcga-gdac', 'gtex')))
      stop('domain not recognized')
    FISS = paste(skidb_env('FISS'), '-d ', domain)

    ALLOWABLE.TYPES = c('pset', 'iset', 'sset', 'indiv', 'sample', 'pair', 'task', 'space')

    if (!(type %in% ALLOWABLE.TYPES))
      stop('type must be one of the following: ', paste(ALLOWABLE.TYPES, collapse = ', '))

    cmd = paste(type, 'list', sep = '_')

    if (type %in% c('pair', 'indiv', 'sample'))
        {
            if (!is.null(set) | type != 'pair')
                str = paste(FISS, ' ', cmd, ' ', wkspace, set)
            else
                str = paste(FISS, ' ', cmd, wkspace, 'all')
        }

    else
      str = paste(FISS, ' ', cmd, ' ', wkspace)

    if (verbose)
        print(str)

    p = pipe(str)
    out = readLines(p)
    close(p)

    return(out)
  }


#####

#####

#' @name fiss_import
#' @title fiss_import
#'
#' @description
#'
#' import individual, sample, or pair annotation to firehose workspace via fiss, writing to temp tables
#' (NOTE: importing a sample annotation = annotating a sample, or setting an annotation for a sample)
#'
#' Value of dat depends on value of type
#'
#' if type = 'sample' or 'pair' then dat is data.frame whose rows are named by the individual or pair_id, depending the value of type, i.e.
#' "pair_id" if type = 'Pair', 'sample_id' if type = 'Sample', and the remaining columns are annotations to add
#' NOTE: annotations already have to be declared from within firehose, not currently possibly to do from FISS
#'
#' if type = 'iset', 'pset', or 'sset' then dat is a two column data.frame with first column having set name and second column having the
#' individuals / pairs / samples belonging to it.
#'
#' NOTE: if type = 'wkspace' then dat is just a named vector specifying workspace annotations
#'
#' @author Marcin Imielinski
#' @export
fiss_import = function(dat,
  type = 'pair', ## can be
  space = Sys.getenv('FIREHOSE_WORKSPACE'),
  wkspace = space,
  domain = 'cga')
  {
    if (is.data.table(dat))
        dat = as.data.frame(dat)

    system(paste('mkdir -p', skidb_env('FISS_TMP')))
    if (!(domain %in% c('cga', 'tcga-gdac', 'gtex')))
      stop('domain not recognized')

    dat[is.na(dat)] = ''
    FISS = paste(skidb_env('FISS'), '-d ', domain)

    if (grepl('pair', type, ignore.case = T))
      {
        ix <- max(grep('^pair_id$', colnames(dat), ignore.case = T), 1)
        pair_id = dat[, ix]
        dat = dat[, -ix, drop = F]

        if (ncol(dat)>0)
          tabs = lapply(colnames(dat), function(x) cbind(data.frame(Pair_Id = pair_id), dat[, x, drop = FALSE]))
        else
          tabs = list(data.frame(Pair_Id = pair_id))

        rnames = paste(skidb_env('FISS_TMP'), '/tmp', runif(length(tabs)), sep = '')
        mapply(write.tab, tabs, rnames)
        system(paste(FISS, 'pair_import ', wkspace, paste(rnames, collapse = ' ')))
        system(paste('rm', paste(rnames, collapse = ' ')))
      }
    else if (grepl('indiv', type, ignore.case = T))
      {
        ix <- max(grep('^individual_id$', colnames(dat), ignore.case = T), 1)
        individual_id = dat[, ix]
        dat = dat[, -ix, drop = F]

        if (ncol(dat)>0)
          tabs = lapply(colnames(dat), function(x) cbind(data.frame(Individual_Id = individual_id), dat[, x, drop = FALSE]))
        else
          tabs = list(data.frame(Individual_Id = individual_id))

        rnames = paste(skidb_env('FISS_TMP'), '/tmp', runif(length(tabs)), sep = '')
        mapply(write.tab, tabs, rnames)
        system(paste(FISS, 'indiv_import ', wkspace, paste(rnames, collapse = ' ')))
        system(paste('rm', paste(rnames, collapse = ' ')))
      }
    else if (grepl('samp', type, ignore.case = T))
      {
        ix <- max(grep('^sample_id$', colnames(dat), ignore.case = T), 1)
        sample_id = dat[, ix]
        dat = dat[, -ix, drop = F]

        if (ncol(dat)>0)
          tabs = lapply(colnames(dat), function(x) cbind(data.frame(Sample_Id = sample_id), dat[, x, drop = FALSE]))
        else
          tabs = list(data.frame(Sample_Id = sample_id))

        rnames = paste(skidb_env('FISS_TMP'), '/tmp', runif(length(tabs)), sep = '')
        mapply(write.tab, tabs, rnames)
        system(paste(FISS, 'sample_import ', wkspace, paste(rnames, collapse = ' ')))
        system(paste('rm', paste(rnames, collapse = ' ')))
      }
    else if (grepl('iset', type, ignore.case = T))
      {
        if (ncol(dat)!=2)
          stop('fiss_import expects two column table, first column set and second column entity when type = "iset", "pset", or "sset"')
        rnames = paste(skidb_env('FISS_TMP'), '/tmp', runif(1), sep = '')
        colnames(dat) = c('Individual_Set_Id', 'Individual_Id')
        write.tab(dat, rnames)
        system(paste(FISS, 'iset_import ', wkspace, rnames))
        system(paste('rm', rnames))
      }
    else if (grepl('pset', type, ignore.case = T))
      {
        if (ncol(dat)!=2)
          stop('fiss_import expects two column table, first column set and second column entity when type = "iset", "pset", or "sset"')
        rnames = paste(skidb_env('FISS_TMP'), '/tmp', runif(1), sep = '')
        colnames(dat) = c('Pair_Set_id', 'Pair_Id')
        write.tab(dat, rnames)
        system(paste(FISS, 'pset_import ', wkspace, rnames))
        system(paste('rm', rnames))
      }
    else if (grepl('sset', type, ignore.case = T))
      {
        if (ncol(dat)!=2)
          stop('fiss_import expects two column table, first column set and second column entity when type = "iset", "pset", or "sset"')
        rnames = paste(skidb_env('FISS_TMP'), '/tmp', runif(1), sep = '')
        colnames(dat) = c('Sample_Set_Id', 'Sample_Id')
        write.tab(dat, rnames)
        system(paste(FISS, 'sset_import ', wkspace, rnames))
        system(paste('rm', rnames))
      }
    else if (grepl('(wkspace)|(workspace)', type, ignore.case = T))
      {
        if (!is.vector(dat))
          stop('input for wkspace annotation import must be a named vector')
        sapply(names(dat), function(x)
               system(paste(FISS, 'annot_set ', wkspace, x, dat[x])))
      }
    else
      stop(paste('type argument value', type, 'not implemented'))
  }

#' @name fiss_new_set
#' @title fiss_new_set
#'
#' creates individual, pair, sample, pair set, iset entities for an input sample set in given wkspace
#'
#' Input table samples must have fields
#' $sample_id
#' $individual_id
#'
#' and can have optional fields
#' $pair_id, $case (logical), each pair_id must occur with one $case = TRUE, and one $case = FALSE
#' $pset_id
#' $iset_id
#'
#' as well as other sample level annotations associated with the sample (eg bam_file_wgs)
#' @description
#'
#' List all entities of a given type in a given workspace.
#'
#' @author Marcin Imielinski
#' @export
fiss_new_set = function(samples, space, verbose = T)
  {
    if (is.null(samples$sample_id) | is.null(samples$individual_id))
      stop('Mandatory fields are sample_id and individual_id')


    samples = as.data.frame(samples)
    individuals = data.frame(individual_id = unique(samples$individual_id))

    if (!is.null(samples$iset_id))
      isets = samples[!duplicated(samples$iset_id, samples$individual_id), c('iset_id', 'individual_id'), ]

    pairs = isets = psets = NULL

    if (!is.null(samples$pair_id))
      if(!is.null(samples$case))
        {
          case.samples = samples[samples$case, ]
          control.samples = samples[!samples$case, ]

          paired = intersect(case.samples$pair_id, control.samples$pair_id)
          if (length(paired>0))
            pairs = data.frame(pair_id = paired, individual_id = case.samples$individual_id[match(paired, case.samples$pair_id)],
              case_sample = case.samples$sample_id[match(paired, case.samples$pair_id)],
              control_sample = control.samples$sample_id[match(paired, control.samples$pair_id)])

          if (verbose)
            cat('Making', length(paired), 'pairs\n')
        }
      else
        warning('must include logical $case field if including pair_id, otherwise will be ignored')

    individuals = data.frame(individual_id = unique(samples$individual_id))

    if (!is.null(samples$iset_id))
      isets = samples[!duplicated(samples$iset_id, samples$pair_id), c('iset_id', 'individual_id')]

    if (!is.null(samples$pset_id) & !is.null(pairs))
      psets = samples[!duplicated(samples$pset_id, samples$pair_id) & samples$pair_id %in% pairs$pair_id, c('pset_id', 'pair_id')]

    scols = c('sample_id', 'individual_id', colnames(samples)[!(colnames(samples) %in% c('sample_id', 'individual_id', 'pair_id', 'iset_id', 'case', 'pset_id'))])
    samples = samples[, scols]
    samples = samples[!duplicated(samples), ]

    ## have to import in this order ..
    ## first individuals --> samples --> pairs --> pair_sets (+/- isets)

    if (verbose)
      cat('Importing individuals \n')

    fiss_import(individuals, space = space, type = 'individual')

    if (!is.null(isets))
      {
        if (verbose)
          cat('Importing individual sets\n')

        fiss_import(isets, space = space, type = 'iset')
      }

    if (verbose)
      cat('Importing samples\n')
    fiss_import(samples, space = space, type = 'sample')

    if (!is.null(pairs))
      {
        if (verbose)
          cat('Importing pairs\n')
        fiss_import(pairs, space = space, type = 'pair')
      }

    if (!is.null(psets))
      {
        if (verbose)
          cat('Importing pair sets\n')
        fiss_import(psets, space = space, type = 'pset')
      }
  }



#' @name fiss_status
#' @title fiss_status
#'
#' @description
#'
#' Gets workflow status for iset in workspace.
#'
#' NOTE: requires that a file ~/.fh.config exists in home directory with a single line: "-u username:password"
#' (where username = fh username, password = fh password) --> this file should be "chmod 700" to maintain (some) pw security
#' List all entities of a given type in a given workspace.
#'
#' @author Marcin Imielinski
#' @export
fh_status = function(workflow = "Capture_Workflow",  # eg Capture_QC_Workflow, WGS_Workflow, WGS_QC_Workflow
  isetname = DEFAULT.ISET,
  wkspace = Sys.getenv('FIREHOSE_WORKSPACE'))
  {
    if (!file.exists(skidb_env('FH.UNPW')))
      stop(sprintf('Need to create file %s, containing single line "-u username:password"', skidb_env('FH.UNPW')))
    else
      fh.unpw = gsub('\\-u ', '', readLines(skidb_env('FH.UNPW')));

    tmp = getURL(sprintf(skidb_env('FH.WKFLOW.URL'), isetname, workflow, wkspace), userpwd = fh.unpw);
    con = textConnection(tmp)
    out= read.delim(con, strings = F, skip = 2, header = TRUE);
    close(con);
    rownames(out) = out$Pipeline;
    out$Pipeline = NULL;

    return(out);
  }


#' @name fh_start
#' @title fh_start
#'
#' @description
#'
#' Pushes button of firehose pipeline "pipeline" in individual set "isetnmae" in workspace "wkspace"
#'
#' NOTE: requires that a file ~/.fh.config exists in home directory with a single line: "-u username:password"
#' (where username = fh username, password = fh password) --> this file should be "chmod 700" to maintain (some) pw security
#'
#' @author Marcin Imielinski
#' @export
fh_start = function(pipeline, #
  isetname = DEFAULT.ISET,
  wkspace = Sys.getenv('FIREHOSE_WORKSPACE'))
  {
    if (!file.exists(skidb_env('FH.UNPW')))
      stop(sprintf('Need to create file %s, containing single line "-u username:password"', skidb_env('FH.UNPW')))
    else
      fh.unpw = gsub('\\-u ', '', readLines(skidb_env('FH.UNPW')));

    out = getURL(sprintf(skidb_env('FH.START.URL'), isetname, gsub(' ', '%20', pipeline), wkspace), userpwd = fh.unpw);
    return(out);
  }


#' @name fh_import
#' @title fh_import
#'
#' @description
#'
#' Import individual set with name iset name "name" and individuals "inds" into firehose
#'
#' NOTE: requires that a file ~/.fh.config exists in home directory with a single line: "-u username:password"
#' (where username = fh username, password = fh password) --> this file should be "chmod 700" to maintain (some) pw security
#'
#' @author Marcin Imielinski
#' @export
fh_import = function(name, inds,  #
  wkspace = Sys.getenv('FIREHOSE_WORKSPACE'))
  {
    if (!file.exists(skidb_env('FH.UNPW')))
      stop(sprintf('Need to create file %s, containing single line "-u username:password"', skidb_env('FH.UNPW')))
    else
      fh.unpw = gsub('\\-u ', '', readLines(skidb_env('FH.UNPW')));

    tmpname = paste('~/.tmp.', round(1e6*runif(1)), '.txt', sep = "");
    write.table(data.frame(individual_set_id = name, individual_id = inds, stringsAsFactors = F), tmpname, row.names = F, sep = "\t", quote = F);
    out = postForm(uri = sprintf(skidb_env('FH.IMPORT.URL'), wkspace), importFile=paste('', file_path_as_absolute(tmpname), sep = ""), .opts = curlOptions(userpwd = fh.unpw));
    system(paste('rm', tmpname));
    return(out);
  }



#' @name get_fh
#' @title get_fh
#'
#' @description
#'
#' Deprecated (pre-fiss) function of accessing firehose
#'
#' Retrieves sample or individual table from firehose for individual set isetname and workspace wkspace.
#' NOTE: isetname can be "All", in which case all individuals in that workspace are retrieved.
#'
#' If isetname, wkspace don't exist will yield a two row table that contains errors.
#'
#' NOTE: requires that a file ~/.fh.config exists in home directory with a single line: "-u username:password"
#' (where username = fh username, password = fh password) --> this file should be "chmod 700" to maintain (some) pw security
#' @author Marcin Imielinski
#' @export
get_fh = function(isetname,
  wkspace = "An_LUAD",
  type = 'Individual', # can also be "Sample"
  merged = TRUE,
  hack.sample.mapping = T, # hack for individual--> sample mapping where we strip -Tumor or -Normal from sample ID, to get individual ID
                           # if F will perform mapping using correct (but currently very slow way) ie FH getcontainers query
  verbose = FALSE)
  {

  #  if (is.null(type))
  #    type = 'Individual'
  #  else
  #    merged = FALSE;

    if (grepl('indiv', type, ignore.case = TRUE))
      type = "Individual"
    else if (grepl('pair', type, ignore.case = TRUE))
      type = 'Pair'
    else if (grepl('samp', type, ignore.case = TRUE))
      type = "Sample"
    else
      stop("Variable 'type' not recognized, should be either sample or individual");

    if (!merged)
      {
        if (!file.exists(skidb_env('FH.UNPW')))
          stop(sprintf('Need to create file %s, containing single line "-u username:password"', skidb_env('FH.UNPW')))
        else
          fh.unpw = gsub('\\-u ', '', readLines(skidb_env('FH.UNPW')));

        if (isetname == 'All')
          {
            ## list all entities in workspace
            tmp = getURL(sprintf(skidb_env('FH.LIST.URL'), type, wkspace), userpwd = fh.unpw);
            con = textConnection(tmp)
            entities = readLines(con);
            close(con);

            ## given entity list, grab annotations all entities in workspace
            CHUNK.SIZE = 100; # bolus size for firehose URL
            entities = entities[nchar(entities)>0];
            out = NULL;
            for (i in seq(1, length(entities), CHUNK.SIZE))
              {
                this.chunk = i:min(length(entities), i+CHUNK.SIZE-1);
                url = sprintf(skidb_env('FH.GET.ENTITIES.URL'), type, wkspace, paste(paste('entityNames=', entities[this.chunk], sep = ""), collapse = '&'));
                tmp = getURL(url, userpwd = fh.unpw);
                con = textConnection(tmp)
                out = rrbind(out, read.delim(con, strings = F, row.names = NULL), union = TRUE);
                close(con);
              }
          }
        else
          {
            tmp = getURL(sprintf(skidb_env('FH.ENTITY.URL'), type, isetname, wkspace), userpwd = fh.unpw);
            con = textConnection(tmp)
            out= read.delim(con, strings = F, row.names = NULL);
            close(con);
          }

        if (is.null(out$sample_id) & is.null(out$individual_id))
          stop('Empty table returned.  Check individual set');

        if (type == "Sample")
          {
            if (hack.sample.mapping)  # the quick and dirty way
              out$individual_id = gsub('(\\-Tumor)|(\\-Normal)$', '', out$sample_id)
            else # the correct and slow way
              out$individual_id = fh_getcontainers(out$sample_id, id.type = 'Sample', container.type = 'Individual', wkspace = wkspace, verbose = verbose);
            rownames(out) = out$sample_id;
          }
        else
          rownames(out) = out$individual_id;

        names(out) = gsub('^sample_id', 'Sample_Id', gsub('^individual_id', 'Individual_Id', names(out)));
        return(out);
      }
    else
      {
        individuals = get_fh(isetname, wkspace = wkspace, type = type, merged = F, verbose = verbose);
        samples = get_fh(isetname, wkspace = wkspace, type = 'Sample', merged = F, verbose = verbose);

        sample_fields = setdiff(names(samples), c('Individual_Id'))

        samples.tumor = samples[samples$sample_type == "Tumor", ]
        samples.normal = samples[samples$sample_type == "Normal", ];

        individuals[, paste('Tumor_', sample_fields, sep = "")] = samples.tumor[match(rownames(individuals), samples.tumor$Individual_Id), sample_fields]
        individuals[, paste('Normal_', sample_fields, sep = "")] =  samples.normal[match(rownames(individuals), samples.normal$Individual_Id), sample_fields]
        return(individuals)
      }
f  }


#' @name fh_containers
#' @title fh_containers
#'
#' @description
#'
#' (deprecated pre-fiss function)
#'
#' call to getcontainers firehose API function for a set of entity ids (eg Sample Id)
#'
#' id entity type should match the entity ids (ie if sample ids should be "Sample", if individual ids should be "Individual")
#'
#' @author Marcin Imielinski
#' @export
fh_getcontainers = function(ids, id.type = 'Sample', container.type = 'Individual', wkspace = Sys.getenv('FIREHOSE_WORKSPACE'), verbose = FALSE)
  {
    if (!file.exists(skidb_env('FH.UNPW')))
      stop(sprintf('Need to create file %s, containing single line "-u username:password"', skidb_env('FH.UNPW')))
    else
      fh.unpw = gsub('\\-u ', '', readLines(skidb_env('FH.UNPW')));

    FH.STOMACH.SIZE = 200;

    out = c()

    for (i in seq(1, length(ids), FH.STOMACH.SIZE))
         {
           ix = i:min(i+FH.STOMACH.SIZE-1, length(ids))

           if (verbose)
             print(sprintf('%s: Getting containers for ids %s:%s', as.character(Sys.time()), min(ix), max(ix)))

           url = sprintf(skidb_env('FH.GETCONTAINERS.URL'), id.type, paste(paste('entityNames', ids[ix], sep = '='), collapse = '&'), container.type, wkspace)
           tmp = getURL(url, userpwd = fh.unpw);
           con = textConnection(tmp)
           tmp.out = readLines(con);
           names.ix = grep('^#', tmp.out)
           out[ix] = tmp.out[names.ix+1]
           names(out)[ix] = gsub('^#Sample\\: ', '', tmp.out[names.ix])
           close(con);
         }

    return(out[ids])
  }



##  _______                                       __        _______   _______   __
## |       \                                     |  \      |       \ |       \ |  \
## | $$$$$$$\  ______    ______    ______    ____| $$      | $$$$$$$\| $$$$$$$\| $$_______
## | $$__/ $$ /      \  /      \  |      \  /      $$      | $$  | $$| $$__/ $$ \$/       \
## | $$    $$|  $$$$$$\|  $$$$$$\  \$$$$$$\|  $$$$$$$      | $$  | $$| $$    $$  |  $$$$$$$
## | $$$$$$$\| $$   \$$| $$  | $$ /      $$| $$  | $$      | $$  | $$| $$$$$$$\   \$$    \
## | $$__/ $$| $$      | $$__/ $$|  $$$$$$$| $$__| $$      | $$__/ $$| $$__/ $$   _\$$$$$$\
## | $$    $$| $$       \$$    $$ \$$    $$ \$$    $$      | $$    $$| $$    $$  |       $$
##  \$$$$$$$  \$$        \$$$$$$   \$$$$$$$  \$$$$$$$       \$$$$$$$  \$$$$$$$    \$$$$$$$



############################

############################


#' @name get_alignment_summaries
#' @title get_alignment_summaries
#'
#' @description
#'
#' get_alignment_summaries
#'
#' Using firehose table grabs the picard alignment summaries and pulls out samplewise statistics for tumor and normal
#' including pf_reads, pf_aligned_bases etc.
#'
#' Also computes ballpark X coverage for wgs and wes using back of the envelope sizes for each (34MB, 3GB)
#'
#' Returns table with Tumor_ and Normal_ columns added at the end
#' call to getcontainers firehose API function for a set of entity ids (eg Sample Id)
#'
#' @author Marcin Imielinski
#' @export
get_alignment_summaries = function(iset, wgs = T, short = FALSE)
{
  if (wgs)
    {
      samples = data.frame(sample_id = c(iset$tumor_sample_id, iset$normal_sample_id),
        bam = c(iset$tumor_clean_bam_file_wgs, iset$normal_clean_bam_file_wgs), stringsAsFactors = F)
      denom = skidb_env('HG.SIZE');
    }
  else
    {
      samples = data.frame(sample_id = c(iset$tumor_sample_id, iset$normal_sample_id),
        bam = c(iset$tumor_clean_bam_file_capture, iset$normal_clean_bam_file_capture), stringsAsFactors = F)
      denom =  36557327;
    }

  samples$asm = gsub('bam$', 'alignment_summary_metrics', samples$bam)
  asm.sum = NULL;
  for (i in 1:nrow(samples))
    {
      if (file.exists(samples$asm[i]))
        {
          asm.data = read.delim(samples$asm[i], strings = F, skip = 6)
          asm.data = asm.data[asm.data$CATEGORY == "PAIR", ][1, , drop = FALSE]
          asm.sum = rbind(asm.sum, data.frame(
            total_reads = sum(as.numeric((asm.data$TOTAL_READS))),
            pf_reads = sum(as.numeric((asm.data$PF_READS))),
            pf_aligned_reads = sum(as.numeric((asm.data$PF_READS_ALIGNED))),
            pf_hq_reads = sum(as.numeric((asm.data$PF_HQ_ALIGNED_READS))),
            pf_aligned_bases = sum(as.numeric((asm.data$PF_ALIGNED_BASES))),
            pf_hq_aligned_bases = sum(as.numeric((asm.data$PF_HQ_ALIGNED_BASES))),
            pf_hq_aligned_q20_bases = sum(as.numeric((asm.data$PF_HQ_ALIGNED_Q20_BASES))),
            mean_read_length = mean(asm.data$MEAN_READ_LENGTH),
            estimated_X_aligned = sum(as.numeric((asm.data$PF_ALIGNED_BASES)))/denom,
            estimated_X_hq_aligned = sum(as.numeric((asm.data$PF_HQ_ALIGNED_BASES)))/denom
            ))
          rownames(asm.sum)[nrow(asm.sum)] = samples$sample_id[i]
        }
    }

  tum.colnames = paste('tumor', colnames(asm.sum), c('wes', 'wgs')[wgs+1], sep = "_")
  norm.colnames = paste('normal',colnames(asm.sum), c('wes', 'wgs')[wgs+1], sep = "_")
  iset[, tum.colnames] = asm.sum[iset$tumor_sample_id, colnames(asm.sum)]
  iset[, norm.colnames] = asm.sum[iset$normal_sample_id, colnames(asm.sum)]

  if (short == TRUE)
    return(iset[, c('Individual_id', tum.colnames, norm.colnames)])
  else
    return(iset)
}


##  _______              ______                                                                    _______   _______   __
## |       \            /      \                                                                  |       \ |       \ |  \
## | $$$$$$$\  ______  |  $$$$$$\ ______    ______    ______   _______    _______   ______        | $$$$$$$\| $$$$$$$\| $$_______
## | $$__| $$ /      \ | $$_  \$$/      \  /      \  /      \ |       \  /       \ /      \       | $$  | $$| $$__/ $$ \$/       \
## | $$    $$|  $$$$$$\| $$ \   |  $$$$$$\|  $$$$$$\|  $$$$$$\| $$$$$$$\|  $$$$$$$|  $$$$$$\      | $$  | $$| $$    $$  |  $$$$$$$
## | $$$$$$$\| $$    $$| $$$$   | $$    $$| $$   \$$| $$    $$| $$  | $$| $$      | $$    $$      | $$  | $$| $$$$$$$\   \$$    \
## | $$  | $$| $$$$$$$$| $$     | $$$$$$$$| $$      | $$$$$$$$| $$  | $$| $$_____ | $$$$$$$$      | $$__/ $$| $$__/ $$   _\$$$$$$\
## | $$  | $$ \$$     \| $$      \$$     \| $$       \$$     \| $$  | $$ \$$     \ \$$     \      | $$    $$| $$    $$  |       $$
##  \$$   \$$  \$$$$$$$ \$$       \$$$$$$$ \$$        \$$$$$$$ \$$   \$$  \$$$$$$$  \$$$$$$$       \$$$$$$$  \$$$$$$$    \$$$$$$$


#' @name get_alignment_summaries
#' @title get_alignment_summaries
#'
#' @description
#'
#' process_uniprot_ft
#'
#' Takes uniprot data file and outputs a "seg" file style data frame / txt of protein features file with columns $ID, $begin, $end,  $feature.name, $feature.description
#'
#' @author = functioun Marcin Imielinski
#' @export
process_uniprot_ft = function(out.file = skidb_env('UNIPROT.FT'), uniprot.dat = skidb_env('UNIPROT.DAT'))
  {
    p = pipe(paste('grep -P "^((FT)|(ID))"', uniprot.dat, ' | wc -l '))
    numfeatures = as.numeric(readLines(p)); # read only FT and ID lines of uniprot file
    close(p)

    feature.name.and.begin = feature.description = ID = rep("", numfeatures)
    begin = end = rep(-1, numfeatures)

    #get field widths
    p = pipe(paste('grep -P "^((ID))"', uniprot.dat));
    open(p);
    fw.prot = get.field.widths(readLines(p,1));
    close(p)

    #get exact field widths for (bizarre doubly horizontally justified) fixed width protein feature row
    p = pipe(paste('grep -P "^((FT))"', uniprot.dat));
    open(p);
    line = paste(readLines(p,1), ' ', sep ="");
    match = sort(c(as.numeric(gregexpr('\\w+[ ]+', line, perl = T)[[1]]), as.numeric(gregexpr('\\s+\\w+', line, perl = T)[[1]]), nchar(line)))
    fw.feature = c(match[3]-match[1], match[6]-match[3], match[8]-match[6], nchar(line)-match[8])
    close(p)
    # iterate through file
    p = pipe(paste('grep -P "^((FT)|(ID))"', uniprot.dat)); # read only FT and ID lines of uniprot file
    open(p)
    k = 0;
    current.protein = NA;
    line = readLines(p,1)
    while (length(line)==1)
      {
        k = k+1;
        if (k%%1000 ==0 ) {print(paste(k, ID[k-1], feature.name.and.begin[k-1], feature.description[k-1], end[k-1], sep = "  %%%% "))}
        if (substr(line, 1, 2) == "ID")
          {
            current.protein = substr(line, fw.prot[1]+1, fw.prot[1]+fw.prot[2])
            ID[k] = current.protein;
            fields = strsplit.fwf(line, fw.prot)
            feature.name.and.begin[k] = gsub('.*(\\d+) AA\\.', 'ProteinBoundary 1', line);
            feature.description[k] = '';
            end[k] = gsub('.* (\\d+) AA\\.', '\\1', line)
          }
        else
          {
            ID[k] = current.protein;
            fields = strsplit.fwf(line, fw.feature)
            feature.name.and.begin[k] = fields[2];  # pt 2 of hack to deal with wierd fixed width left vs right justified uniprot text file format
            feature.description[k] = fields[4];
            end[k] = fields[3];
          }
        line = readLines(p,1)
      }

    tmp = strsplit(gsub('\\s+', '\t', feature.name.and.begin, perl = T), '\t') # pt 3 of hack to deal with wierd fixed width left vs right justified uniprot text file format
    out = data.frame(ID = trim(ID), feature.name = trim(sapply(tmp, function(x) x[[1]])), feature.description = trim(feature.description),
      begin = as.numeric(sapply(tmp, function(x) if (length(x)>1) x[[2]] else NA )), end = as.numeric(end), stringsAsFactors = F)

    # clean up multiline protein feature descriptions
    na = which(is.na(out$end));
    tmp = out$feature.description;
    for (i in rev(na))
      tmp[i-1] = paste(tmp[i-1], tmp[i])
    out$feature.description = tmp;
    out = out[-na, ];

    # match up feature_names with feature_tpes
    up.keys = read.delim(skidb_env('UNIPROT.FEATUREKEYS'), strings = F);
    out$feature.type = up.keys$feature_type[match(out$feature.name, up.keys$feature)]

    if (!is.null(out.file))
      write.table(out, out.file, sep = "\t", row.names = F, quote = F)

    return(out)
  }

#' @name process_uniprot_acmap
#' @title process_uniprot_acmap
#'
#' @description
#'
#' Takes uniprot data file and outputs a two column table mapping id's to accession numbers
#'
#' @author Marcin Imielinski
#' @export
process_uniprot_acmap = function(out.file = skidb_env('UNIPROT.ACMAP'), uniprot.dat = skidb_env('UNIPROT.DAT'))
  {
    p = pipe(paste('grep -P "^AC"', uniprot.dat, ' | grep -Po ";" | wc -l'))
    numids = as.numeric(readLines(p)); # count up number of unique accession id's (to preallocate)
    close(p)

    ## get field widths for efficient parsing later on
    p = pipe(paste('grep -P "^((ID))"', uniprot.dat));
    open(p);
    fw.id = get.field.widths(readLines(p,1));
    close(p)

    p = pipe(paste('grep -P "^((AC))"', uniprot.dat));
    open(p);
    fw.ac = get.field.widths(readLines(p,1));
    close(p)

    ID = AC = rep("", numids); #preallocate output table columns (for some reason R is much faster at updating big vectors than big data frames)

    # iterate through file
    p = pipe(paste('grep -P "^((AC)|(ID))"', uniprot.dat)); # read only AC and ID lines of uniprot file
    open(p)
    k = 1;
    current.protein = NA;
    line = readLines(p,1)
    while (length(line)==1)
      {
        if (substr(line, 1, 2) == "ID")
            current.protein = substr(line, fw.id[1]+1, fw.id[1]+fw.id[2])
        else
          {
            ac = strsplit(substr(line, fw.ac[1],nchar(line)), ';')[[1]]
            these.ind = k:(k+length(ac)-1);
            ID[these.ind] = current.protein;
            AC[these.ind] = ac;
            k = k+length(ac);
          }
        line = readLines(p,1)
      }

    out = data.frame(ID = trim(ID), AC = trim(AC), stringsAsFactors = F)

    if (!is.null(out.file))
      write.table(out, out.file, sep = "\t", row.names = F, quote = F)

    return(out)
  }

#' @name process_pfam
#' @title process_pfam
#'
#' @description
#'
#' Given pfam msa file (stockholm format)
#' output data as list of character matrices with rownames corresponding to sequence data
#' and list names corresponding to pfam domain IDs.
#'
#' @author Marcin Imielinski
#' @export
process_pfam= function(pfam.file = skidb_env('PFAM.HUMAN'), out.file = skidb_env('PFAM.HUMAN.RDS'), lim = Inf)
  {
    pali = read_stockholm(pfam.file);

    ix = sapply(pali, function(x) length(rownames(x)))>0
    if (any(!ix))
      warning(paste('Some pfam entries without rownames:', paste(names(pali)[!ix], collapse = ', ')))

    pali = pali[which(ix)]
    l.ix = unlist(lapply(1:length(pali), function(x) rep(x, nrow(pali[[x]]))))
    pali.name = names(pali)[l.ix]
    tmp.txt = unlist(lapply(pali, function(x) trim(rownames(x))))
    tmp.mat = t(matrix(unlist(strsplit(tmp.txt, '[-\\/]')), nrow = 3))
    keep.ix = tmp.mat[,1] %in% seqlevels(up.ft)
    tmp.mat = tmp.mat[keep.ix, ]
    keep.ix.l = split(keep.ix, l.ix)
    pname= names(pali);
    pali = lapply(1:length(pali), function(x) pali[[x]][keep.ix.l[[x]], , drop = F])
    names(pali) = pname;
    pali = pali[sapply(1:length(pali), function(x) any(keep.ix.l[[x]]))]
    l.ix = l.ix[keep.ix]
    pali.name = pali.name[keep.ix]
    gr = suppressWarnings(gr.fix(gr.fix(GRanges(tmp.mat[,1], IRanges(as.numeric(tmp.mat[,2]), as.numeric(tmp.mat[,3])), strand = '+'), uniprot_seqinfo())))
    grl = split(gr, pali.name)[names(pali)]
    up2pf = read_pfam_meta()
    up2pf$Uniprot = up2pf$UniprotID
    up2pf$Pfam = up2pf$fullacc
    out = list(grl = grl, pali = pali, up2pf = up2pf[, c('Uniprot', 'Pfam', 'desc')])

    saveRDS(out, out.file)
  }

#' @name process_pfam_metaa
#' @title process_pfam_meta
#'
#' @description
#'
#' processes pfam fasta from pfam fasta and main pfam file
#'
#' @author Marcin Imielinski
#' @export
process_pfam_meta = function(pfam.fa.file = skidb_env('PFAM.FA'), outfile = skidb_env('PFAM.HUMAN.META.RDS'), verbose = T)
{
  if (verbose)
    cat('Reading pfam fasta\n')

  if (verbose)
    cat('parsing pfam main\n')
  pfam = read_pfam()
  bla = lapply(pfam, rownames);
  pfam.names = sort(trim(structure(unlist(bla), names = unlist(lapply(1:length(pfam), function(x) rep(names(pfam)[x], length(bla[[x]])))))))
  pfam.names = pfam.names[!duplicated(pfam.names)]
  tmp = strsplit(pfam.names, '[ \\/\\]' )

  pfam.meta = as.data.frame(t(matrix(unlist(tmp), ncol = length(tmp))))
  pfam.meta[,3] = sapply(strsplit(as.character(pfam.meta[,2]), '-'), function(x) x[[2]])
  pfam.meta[,2] = sapply(strsplit(as.character(pfam.meta[,2]), '-'), function(x) x[[1]])
  pfam.meta[,4] = names(pfam.names)
  pfam.meta[,5] = sapply(strsplit(as.character(pfam.meta[,4]), '[;\\.]'), function(x) x[[1]])
  names(pfam.meta) = c('UniprotID', 'begin', 'end', 'fullacc', 'hmmacc')

  if (verbose)
    cat('parsing pfam fasta\n')
  pfam.fa = readDNAStringSet('~/links/tcga_marcin/DB/Pfam/Pfam-A.fasta')
  pfam.fa = pfam.fa[grep('HUMAN', names(pfam.fa))]
  tmp = strsplit(names(pfam.fa), '[ \\/\\]' )
  pfam.meta2 = as.data.frame(t(matrix(unlist(tmp), ncol = length(tmp))))
  pfam.meta2[,3] = sapply(strsplit(as.character(pfam.meta2[,2]), '-'), function(x) x[[2]])
  pfam.meta2[,2] = sapply(strsplit(as.character(pfam.meta2[,2]), '-'), function(x) x[[1]])
  pfam.meta2[,5] = sapply(strsplit(as.character(pfam.meta2[,4]), '[;\\.]'), function(x) x[[1]])
  names(pfam.meta) = c('UniprotID', 'begin', 'end', 'desc', 'hmmacc')

  pfam.meta$desc = pfam.meta2$desc[match(pfam.meta$hmmacc, pfam.meta2$hmmacc)]
  saveRDS(pfam.meta, outfile)
}


#' @name read_stockholm
#' @title read_stockholm
#'
#' @description
#'
#' Reads Pfam stockholm format into a list of matrices named by each record accession
#' Each row of each matrix comprises a sequence alignment.
#'
#' @author Marcin Imielinski
#' @export
read_stockholm = function(pfam.file, lim = Inf)
  {
    p = pipe(paste('grep -P "^//"', pfam.file, ' | wc -l'))
    numdom = as.numeric(readLines(p)); # count up number of unique accession id's (to preallocate)
    close(p)

    numdom = pmin(lim, numdom);

    ## scan once to pre-allocate matrices
    f = file(pfam.file, 'r')
    line = readLines(f,1)

    dom.dim = matrix(0, ncol = 2, nrow = numdom);
    text.pos = matrix(0, ncol = 2, nrow = numdom);
    dom.names = rep('', numdom);
    i = 0;
    waiting = FALSE;
    while (length(line)==1 & i <= lim)
      {
        if (substr(line, 1, 11) == '# STOCKHOLM')
          {
            if ((i %% 1000)==0)
              cat(sprintf('Read %s Pfam records ..\n', i))
            i = i+1;
            waiting = TRUE
          }
        else if (substr(line, 1, 7) == '#=GF AC')
          dom.names[i] = substr(line, 11, nchar(line))
        else if (substr(line, 1, 1) != '#' & substr(line, 1, 1) != '/')
        {
          if (waiting)
            {
              sp.ix = which(strsplit(line, '')[[1]]==' ')
              if (length(sp.ix)>0)
                {
                  text.pos[i,1] = min(sp.ix)
                  text.pos[i,2] = max(sp.ix)+1
                  dom.dim[i, 2] = nchar(line)-text.pos[i,2]+1
                }
              waiting = FALSE;
            }
          dom.dim[i, 1] = dom.dim[i, 1] + 1;
        }
        last = line;
        line = readLines(f, 1);
      }
    close(f)

    # pre allocate matrices
    out = lapply(1:numdom, function(x) matrix(NA, nrow = dom.dim[x, 1], ncol = dom.dim[x, 2]))
    names(out) = dom.names;

    # start over
    f = file(pfam.file, 'r')
    line = readLines(f,1)

    i = 0;
    waiting = FALSE;
    while (length(line)==1 & i <= lim)
      {
        if (substr(line, 1, 11) == '# STOCKHOLM')
          {
            if ((i %% 1000)==0)
              cat(sprintf('Populated %s Pfam records ..\n', i))

            i = i+1;
            j = 1;

            if (i>1)
              rownames(out[[i-1]]) = seq.names

            if (i<=nrow(dom.dim))
              seq.names = rep('', dom.dim[i, 1]);
          }
        else if (substr(line, 1, 1) != '#' & substr(line, 1, 1) != '/')
        {
          seq.names[j] = substr(line, 1, text.pos[i,2]-1)
          out[[i]][j, ] = strsplit(substr(line, text.pos[i,2], nchar(line)), '')[[1]]
          j = j+1;
        }
        line = readLines(f, 1);
      }
    close(f)

    rownames(out[[length(out)]]) = seq.names

    out = out[which(sapply(out, ncol)!=0)]
    return(out)
  }


#' @name process_uniprot_rs
#' @title process_uniprot_rs
#'
#' @description
#'
#' Takes uniprot data file and outputs a two column table mapping id's to refseq numbers
#'
#' @author Marcin Imielinski
#' @export
process_uniprot_rs = function(out.file = skidb_env('UNIPROT.RS'), uniprot.dat = skidb_env('UNIPROT.DAT'))
  {
    p = pipe(paste('grep -P "^DR\\s+RefSeq"', uniprot.dat, ' | grep -Po ";" | wc -l'))
    numids = as.numeric(readLines(p)); # count up number of unique accession id's (to preallocate)
    close(p)

    ## get field widths for efficient parsing later on
    p = pipe(paste('grep -P "^((ID))"', uniprot.dat));
    open(p);
    fw.id = get.field.widths(readLines(p,1));
    close(p)

    p = pipe(paste('grep -P "^((DR\\s+RefSeq;))"', uniprot.dat));
    open(p);
    fw.rs = get.field.widths(readLines(p,1));
    close(p)

    ID = RS = rep("", numids); #preallocate output table columns (for some reason R is much faster at updating big vectors than big data frames)

    # iterate through file
    p = pipe(paste('grep -P "^((DR\\s+RefSeq)|(ID))"', uniprot.dat)); # read only DR RefSeq  and ID lines of uniprot file
    open(p)
    k = 1;
    current.protein = NA;
    line = readLines(p,1)
    while (length(line)==1)
      {
        if (substr(line, 1, 2) == "ID")
            current.protein = substr(line, fw.id[1]+1, fw.id[1]+fw.id[2])
        else
          {
            rs = strsplit(substr(line, fw.rs[1],nchar(line)-1), ';')[[1]]
            rs = rs[2:length(rs)];
            these.ind = k:(k+length(rs)-1);
            ID[these.ind] = current.protein;
            RS[these.ind] = rs;
            k = k+length(rs);
          }
        line = readLines(p,1)
      }

    out = data.frame(ID = trim(ID), RS = trim(RS), stringsAsFactors = F)

    if (!is.null(out.file))
      write.table(out, out.file, sep = "\t", row.names = F, quote = F)

    return(out)
  }

#' @name read_uniprot_ft
#' @title read_uniprot_ft
#'
#' @description
#'
#' reads uniprot ft from raw flat file dump
#'
#' @author Marcin Imielinski
#' @export
read_uniprot_ft = function(granges = F, grl = T)
{
  up.ft = read.delim(skidb_env('UNIPROT.FT'), strings = F)
  if (granges | grl)
    {
      up.ft.tmp = up.ft[!is.na(up.ft$begin), ]
      si.tmp = aggregate(formula = end ~ ID, data = up.ft.tmp, FUN = max);
      si = Seqinfo(si.tmp[,1], si.tmp[,2]);
      up.ft = GRanges(up.ft.tmp$ID, IRanges(up.ft.tmp$begin, up.ft.tmp$end), strand = '*', seqlengths = seqlengths(si))
      values(up.ft) = up.ft.tmp[, setdiff(names(up.ft.tmp), c("seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "genome", "start", "end", "width", "element"))];
      up.ft$label = up.ft.tmp$feature.description
      names(up.ft) = NULL

      if (grl)
        up.ft = split(up.ft, seqnames(up.ft))

    }

  return(up.ft)
}


#' @name read_pfam
#' @title read_pfam
#'
#' @description
#'
#' reads pfam dump from saved .rds path location
#'
#' @author Marcin Imielinski
#' @export
read_pfam = function()
{
  return(readRDS(skidb_env('PFAM.HUMAN.RDS')))
}


#' @name process_repeatmasker
#' @title process_repeatmasker
#'
#' @description
#'
#' process repeat masker file into GRanges
#'
#' @author Marcin Imielinski
#' @export
process_repeatmasker = function(infile = skidb_env('REPEATMASKER.HG19'), outfile = skidb_env('REPEATMASKER.HG19.GR'), chrsub = T)
{
  tmp = read.delim(infile, strings = F)
  if (chrsub)
    tmp$genoName = gsub('chr', '', tmp$genoName)
  out = GRanges(tmp$genoName, IRanges(tmp$genoStart + 1, tmp$genoEnd + 1), strand = tmp$strand)
  values(out) = tmp[, c('repName', 'repClass', 'repFamily', 'swScore', 'id')]
  saveRDS(out, outfile)
}


#' @name read_pfam_meta
#' @title read_pfam_meta
#'
#' @description
#'
#' reads meta dat for pfam from cached rds file
#'
#' @author Marcin Imielinski
read_pfam_meta = function(gr = F)
{
  pfam.meta = readRDS(skidb_env('PFAM.HUMAN.META.RDS'))

  if (gr)
    {
      out = GRanges(pfam.meta$UniprotID, IRanges(as.numeric(pfam.meta$begin), as.numeric(pfam.meta$end)), desc = pfam.meta$desc, hmmacc = pfam.meta$hmmacc, fullacc = pfam.meta$fullacc);
      out = gr.fix(out, uniprot_seqinfo(), drop = T)
    }
  else
    return(pfam.meta)
}


#' @name read_uniprot_map
#' @title read_uniprot_map
#'
#' @description
#'
#' reads uniprot map from file
#'
#' @author Marcin Imielinski
read_uniprot_map = function(mapto = NULL)
{
  if (is.null(mapto))
    stop(sprintf('mapto argument options are:\n %s', paste(gsub('map_Uniprot2(.*)\\.txt$', '\\1', dir(skidb_env('UNIPROT.ROOT'), 'map_Uniprot2')), collapse = ", ")))

  return(read.delim(sprintf('%s/map_Uniprot2%s.txt', skidb_env('UNIPROT.ROOT'), mapto), strings = F))
}

#' @name read_hg_bionano
#' @title read_hg_bionano
#'
#' @description
#'
#' reads bionano human genome into GRanges format
#'
#' @author Marcin Imielinski
#' @export
read_hg_bionano = function(fn = CMAP.HG19)
  {
    cmap = .read_bionano(fn)

    cmap$CMapId = names(hg_seqlengths())[cmap$CMapId]
    cmap = cmap[!is.na(cmap$CMapId), ]
    out = GRanges(cmap$CMapId, IRanges(cmap$Position, width = 1),  strand = '*', seqlengths = hg_seqlengths())
    values(out) = cmap
    return(out)
  }

#' @name read_lumpy
#' @title read_lumpy
#'
#' @description
#'
#' reads LUMPY sv bedpe format, returning rearrangement junctions as GRangesList
#' with strands pointing AWAY from junction
#'
#' @author Marcin Imielinski
#' @export
read_lumpy = function(file, sl = hg_seqlengths(), samples = c('tumor', 'normal'))
  {
    df = read.delim(file, strings = F, header = F)
    names(df) = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'id', 'score', 'str1', 'str2', 'type', 'info', 'strcfg', 'bppos', 'bpint')

    if (any(grepl('IDS', df$info)))
      {
        tmp = strsplit(gsub('IDS:', '', df$info), ';')
        tmp.i = rep(1:length(tmp), sapply(tmp, length))
        tmp2 = strsplit(unlist(tmp), ',')
        tmp.j = as.numeric(sapply(tmp2, function(x) x[[1]]))
        tmp.v = as.numeric(sapply(tmp2, function(x) x[[2]]))
        dat = sparseMatrix(tmp.i, tmp.j, x = tmp.v)
        for (i in 1:length(samples))
          {
            this.col = paste(samples[i], 'reads', sep = '_')
            df[, this.col] = 0
            if (i<=ncol(dat))
                df[, this.col] = dat[,i]
          }
      }

    new.sn = setdiff(c(df$chr1, df$chr2), names(sl))
    if (length(new.sn)!=0)
      sl[new.sn] = NA

    ra1 = GRanges(df$chr1, IRanges(df$start1, df$end1), strand = ifelse(df$str1=='+', '-', '+'), seqlengths = sl)
    ra2 = GRanges(df$chr2, IRanges(df$start2, df$end2), strand = ifelse(df$str2=='+', '-', '+'), seqlengths = sl)

    out = grl.pivot(GRangesList(ra1, ra2))
    values(out) = df

  }

#' @name read_uniprot_seqinfo
#' @title read_uniprot_seqinfo
#'
#' @description
#'
#' extracts seqinfo from uniprot fa file
#'
#' @author Marcin Imielinski
#' @export
uniprot_seqinfo = function()
  {
    up.fa = read_uniprot_fa()
    tmp = structure(names = names(up.fa), width(up.fa));
    tmp = tmp[!duplicated(names(tmp))]
    return(Seqinfo(names(tmp), tmp))
  }

#' @name read_uniprot_fa
#' @title read_uniprot_fa
#'
#' @description
#'
#' reads uniprot fasta file into BioStrings BStringSet format
#'
#' @author Marcin Imielinski
#' @export
read_uniprot_fa = function(trim.names = T #if false will output full identifiers as stored in uniprot fasta, otherwise only uniprot id's (e.g. U2AF1_HUMAN)
  )
{
  out = as(readBStringSet(skidb_env('UNIPROT.FA'), format = "fasta"), 'AAStringSet')

  if (trim.names)
    names(out) = gsub('^[^\\|]+\\|[^\\|]+\\|(\\w+)+ .*$', '\\1', names(out))

  return(out)
}

#' @name read_rsem
#' @title read_rsem
#'
#' @description
#'
#' process rsem files from UNC
#'
#' @author Marcin Imielinski
#' @export
read_rsem = function(file)
  {
    out = read.delim(file, strings = F, skip = 1)

    if (!is.null(out$gene_id))
      {
        out$gene_sym = sapply(strsplit(out$gene_id, '\\|'), function(x) x[[1]])
        out$gene_id = sapply(strsplit(out$gene_id, '\\|'), function(x) x[[2]])
      }

    rownames(out) = NULL

    return(out)
  }

#' @name read_snpdata
#' @title read_snpdata
#'
#' @description
#'
#' reads snp data given array definition file in either rds data table format
#' or txt format
#'
#' returns GRanges or list with $gr and $td for GRanges and length(2) trackData of total and allelic copy number
#'
#' process rsem files from UNC
#'
#' @author Marcin Imielinski
#' @export
read_snpdata = function(file, markermap = skidb_env('SNP.MARKER.MAP'), td = FALSE)
{
    if (!is.data.table(markermap))
        snpmarkers = tryCatch(readRDS(markermap), error = function(e) NULL)
    else
        snpmarkers = markermap

    if (is.null(snpmarkers))
        {
            warning("Reloading markermap file from scratch")
            snpmarkers = fread(markermap, header = F)
            setnames(snpmarkers, c("name", "chr", "pos"))
            setkey(snpmarkers, name)
        }

    snpdata = fread(file, skip = 1);
    setnames(snpdata, c("name", "A", "B"))
    snpdata[, ":="(chr = snpmarkers[name, chr], pos = snpmarkers[name, pos])]
    gr = seg2gr(snpdata)
    gr$high = pmax(gr$A, gr$B); gr$low = pmin(gr$A, gr$B); gr$total = gr$A + gr$B;

    if (td) ## return both trackData and GRanges object
        {
            gr.low = gr.high = gr[, c()];
            gr.low$signal = gr$low;
            gr.high$signal = gr$high;
            gr.low$type = "low";
            gr.high$type = "high";
            gr.allele = c(gr.low, gr.high);
            td.allele = trackData(gr.allele, y.field = 'signal', colormaps = list(type = c(high = "red", low = "blue")));
            td.total = trackData(gr, y.field = 'total')
            td = c(td.allele, td.total)
            return(list(gr = gr, td = td))
        }
    else
        return(gr)
}


#' @name read_refgene
#' @title read_refgene
#'
#' @description
#'
#' if exons = T, returns an expanded data frame with 1 row per exon
#' if uniprot = T, maps to uniprot annotations, and limits only to those refseq genes that map to at least
#' one uniprot record
#'
#' @author Marcin Imielinski
#' @export
read_refGene = function(hg19 = TRUE, exons = F, cds = F, granges = F, grl = F, chr.numeric = F, chr = F, uniprot = T, ccds = T, force.make = F)
{
  if (!force.make & hg19) ## cached versions
      {
          if (granges & !exons & file.exists(skidb_env('REFGENE.FILE.HG19.GR')))
              return(readRDS(skidb_env('REFGENE.FILE.HG19.GR')))
          if (grl & !cds &  file.exists(skidb_env('REFGENE.FILE.HG19.GRL')))
              return(readRDS(skidb_env('REFGENE.FILE.HG19.GRL')))
          if (grl & cds &  file.exists(skidb_env('GENCODE.FILE.HG19.CDS.GRL')))
              return(readRDS(skidb_env('GENCODE.FILE.HG19.CDS.GRL')))
      }

  if (hg19)
      rg = read.delim(REFGENE.FILE.HG19, row.names=NULL, header=0, strings = F)
  else
      rg = read.delim(skidb_env('REFGENE.FILE.HG18'), row.names=NULL, header=0, strings = F)

  colnames(rg) = c("num", "NM_id", "chr", "strand", "s1", "e1", "s2", "e2", "n_exons", "exon_starts", "exon_ends", "num2", "gene_sym", "cmpl1", "cmpl2", "exon_frame")
  if (chr.numeric)
    {
      chrs = rg[,"chr"]
      nchrs = gsub('MT', '25', gsub('Y', 24, gsub('X', '23', gsub("chr", "", chrs))))
      keep_chrs = as.character( c(1:25) )
      keep_ix = (nchrs %in% keep_chrs)
      rg = rg[ keep_ix, ]
      rg[,"chr"] = as.numeric( nchrs[keep_ix] )
    }

  if (!chr)
    rg$chr = gsub('chr', '', rg$chr)

  rg$start = rg$s1+1; ## we don't like 0 based coordinates for ranges
  rg$end = rg$e1;
  rg$transcript.length = rg.length(rg)

  if (uniprot)
    {
      up.rs = read_uniprot_map('RefSeq_NT')
      up.rs$NM_id = sapply(strsplit(up.rs$CrossRefID, '\\.'), function(x) x[[1]])
      up.kb = read_uniprot_map('UniProtKB-ID');
      up.kb$UP_id  = up.kb[,2]
      up.final = merge(up.kb, up.rs, by = 'Uniprot.ID', all = T)[, c('UP_id', 'NM_id')]
      up.final = up.final[rev(order(up.final$UP_id %in% names(up.kb))), ]
      up.final = up.final[!is.na(up.final$NM_id), ]
      up.final = up.final[order(grepl('^Q', up.final$UP_id)), ]
      rg$Uniprot = up.final$UP_id[match(rg$NM_id, up.final$NM_id)]
      rg = rg[!is.na(as.vector(rg$Uniprot)), ]
    }

  if (ccds)
    {
      ccds.map = read.delim(skidb_env('CCDS.MAP.FILE'), strings = F)
      ccds.map$nid = gsub('\\..*', '', ccds.map$nucleotide_ID)
      rg$CCDS = ccds.map[match(rg$NM_id, ccds.map$nid), ]$X.ccds
    }

  if (grl | exons | cds)
    {
      rg$uid = 1:nrow(rg)
      rg = rg2exons(rg, cds = cds);
    }

  if (granges | grl)
    {
      tmp = GRanges(rg$chr, IRanges(rg$start, rg$end), strand = rg$strand, seqlengths = hg_seqlengths(include.junk = T, chr = chr, hg19 = hg19))
      values(tmp) = rg[, setdiff(names(rg), c('chr', 'start', 'end', 'strand'))]
      rg = tmp;

      if (!exons & !cds)
        saveRDS(rg, skidb_env('REFGENE.FILE.HG19.GR'))

      if (grl)
        {
          # separate tx.level and exon level attributes
          tx.cols = intersect(c('transcript.length', 'NM_id', 'uid',  'n_exons', 's1', 's2', 'e1', 'e2', 'gene_sym', 'Uniprot', 'CCDS', 'exon_starts', 'exon_ends', 'num', 'chr', 'pos1', 'pos2'), names(values(rg)))
          tx.vals = values(rg)[, tx.cols]
          uix = !duplicated(rg$uid)
          tx.str = structure(as.character(strand(rg)[uix]), names = unique(rg$uid))
          tx.chr = structure(as.character(seqnames(rg)[uix]), names = unique(rg$uid))
          tx.vals2 = tx.vals
          tx.vals = tx.vals[!duplicated(tx.vals$uid), ]
          values(rg) = values(rg)[, c('uid', setdiff(names(values(rg)), tx.cols))]
          ix = order(c(-1, 1)[1+as.numeric(strand(rg)=='+')] * start(rg))
          rg = split(rg[ix], rg$uid[ix])
          values(rg) = tx.vals[match(names(rg), tx.vals$uid),
                  setdiff(tx.cols, c('exon_starts', 'exon_ends', 'exon_frame', 'num', 'chr', 'pos1', 'pos2'))]
          values(rg)$chr = tx.chr[values(rg)$uid]
          values(rg)$str = tx.str[values(rg)$uid]
          values(rg)$s1 = values(rg)$s1+1 ## fix 0 based coordinates for output, otherwise these will be confusing downstream
          values(rg)$s2 = values(rg)$s2+1
          values(rg)$uid = NULL;

          if (!cds)
            saveRDS(rg, skidb_env('REFGENE.FILE.HG19.GRL'))
          else
            saveRDS(rg, skidb_env('REFGENE.FILE.HG19.CDS.GRL'))
        }
    }

  return(rg)
}

#' @name read_ensembl
#' @title read_ensembl
#'
#' @description
#'
#' reads ensembl transcripts into genomic ranges
#'
#' @author Marcin Imielinski
#' @export
read_ensembl = function(
  transcript = F, # if FALSE will only output genes (not transcript ranges)
  flat = F, # if TRUE will flatten exons across genes and output a gr with exon labels
  grl = F, # if TRUE will output grl of exons for each transcript
  types = '')
  {
    e.gene = read.delim(skidb_env('ENSEMBL.GENE'), strings = F)
    other.fields = setdiff(colnames(e.gene),
      c('Chromosome.Name', 'Gene.Start..bp.','Gene.End..bp.', 'Transcript.Start..bp.',
        'Transcript.End..bp', 'Ensembl.Transcript.ID', 'Ensembl.Gene.ID'))

    if (any(nchar(types)>0))
      e.gene = e.gene[e.gene$Gene.Biotype %in% types | e.gene$Transcript.Biotype %in% types, ]

    if (transcript | grl | flat)
      {
        if (grl | flat)
          {
            e.tx = read.delim(skidb_env('ENSEMBL.TX'), strings = F, header = F)
            e.tx = e.tx[e.tx[,1] %in% e.gene$Ensembl.Transcript.ID, ]
            tmp = lapply(strsplit(e.tx[,2], ','), function(x) strsplit(x, ':|\\-'))
            nexon = sapply(tmp, length)
            tx.id = Rle(e.tx[,1], nexon)
            ix = match(e.tx[,1], e.gene$Ensembl.Transcript.ID)
            str = Rle(c('-', '+')[1+(e.gene$Strand[ix]>0)], nexon)
            exons = matrix(unlist(tmp), ncol = 3, byrow = T)
            out = split(GRanges(exons[,1], IRanges(as.numeric(exons[,2]), as.numeric(exons[,3])), strand = str), tx.id);
            values(out)$tx.id = e.gene$Ensembl.Transcript.ID[ix]
            values(out)$gene.id = e.gene$Ensembl.Gene.ID[ix]
            values(out)$gene = e.gene$Associated.Gene.Name[ix]
            names(out) = e.gene$Ensembl.Transcript.ID[ix]

            if (flat)
              {
                tmp.out = unlist(out)
                ix = match(names(tmp.out), names(out))
                values(tmp.out) = values(out)[ix, ]
                out = unique(tmp.out)
                names(out) = NULL
              }
          }
        else
          {
            out = GRanges(e.gene$Chromosome.Name, IRanges(as.numeric(e.gene$Transcript.Start..bp.), as.numeric(e.gene$Transcript.End..bp.)),
              strand = c('-', '+')[1+(e.gene$Strand>0)])
            values(out) = e.gene[, other.fields]
            out$tx.id = e.gene$Ensembl.Transcript.ID
            out$gene.id = e.gene$Ensembl.Gene.ID
            out$gene = e.gene$Associated.Gene.Name;
            names(out) = e.gene$Ensembl.Transcript.ID
          }
      }
    else
      {
        e.gene = e.gene[!duplicated(e.gene$Ensembl.Gene.ID), ]
        out = GRanges(e.gene$Chromosome.Name, IRanges(e.gene$Gene.Start..bp., e.gene$Gene.End..bp.), strand = c('-', '+')[1+(e.gene$Strand>0)])
        values(out) = e.gene[, other.fields]
        out$gene.id = e.gene$Ensembl.Gene.ID
        out$gene = e.gene$Associated.Gene.Name;
        names(out) = e.gene$Ensembl.Gene.ID
      }

    out = gr.fix(out, gr.sub(seqinfo2gr(read_hg()), 'chr', ''))

    return(out)
  }


#' @name read_kg
#' @title read_kg
#'
#' @description
#'
#' reads knownGene and refGene and crosses the two (left join) to yield an uber table
#'
#' @author Marcin Imielinski
#' @export
read_kg = function(hg19 = T, chr = F, gr = F,
  grl = F) # grl = T will override gr = T
{
  rg = read_refGene(hg19)[, c('NM_id', 'gene_sym')]

  if (hg19)
    {
      kg = read.delim(skidb_env('KG.FILE.HG19'), strings = F)
      kg2rg = read.delim(skidb_env('KG2RG.FILE.HG19'), strings = F, header = F)
    }
  else
    {
      kg = read.delim(skidb_env('KG.FILE.HG18'), strings = F)
      kg2rg = read.delim(skidb_env('KG2RG.FILE.HG18'), strings = F, header = F)
    }

  names(kg)[1] = 'name';
  names(kg2rg) = c('name', 'NM_id')

  kg = merge(kg, kg2rg, by = 'name', all.x = T);
  kg = merge(kg, rg, by = 'NM_id', all.x = T);

  if (!chr)
    kg$chrom = gsub('chr', '', kg$chrom);

  if (grl)
    {
      val = kg[, setdiff(names(kg), c('strand', 'exonStarts', 'exonEnds'))];
      val$str = kg$strand;
      st = strsplit(kg$exonStarts, ',')
      en = strsplit(kg$exonEnds, ',')
      chr = lapply(1:length(st), function(x) rep(kg$chrom[x], length(st[[x]])))
      str = lapply(1:length(st), function(x) rep(kg$strand[x], length(st[[x]])))
      id = lapply(1:length(st), function(x) rep(x, length(st[[x]])));
      tmp = GRanges(unlist(chr), IRanges(as.numeric(unlist(st)), as.numeric(unlist(en))), strand = unlist(str))
      kg = split(tmp, id);
      values(kg) = val;
      names(kg) = val$name;
    }
  else if (gr)
    {
      val = kg[, setdiff(names(kg), c('chrom', 'strand', 'txStart', 'txEnd', 'exonStarts', 'exonEnds'))];
      kg = GRanges(kg$chrom, IRanges(kg$txStart, kg$txEnd), strand = kg$strand);
      values(kg) = val;
      names(kg) = val$name;
    }

  return(kg)
}

#' @name process_gaf
#' @title process_gaf
#'
#' @description
#'
#' processes gaf from source (ie obtained from TCGA) into skidb_env('GAF.DIR')
#' @author Marcin Imielinski
#' @export
process_gaf= function(gaf.file = skidb_env('GAF.SOURCE'), verbose = T)
  {
    if (verbose)
      print(sprintf('Reading gaf from %s and splitting across FeatureType', gaf.file))

    gaf = read.delim(gaf.file, strings = F)
    gafs = split(gaf, gaf$FeatureType);

    if (verbose)
      print(sprintf('Writing split gaf subtables (%s) as rds files to central directory (%s)', paste(names(gafs), collapse = ", "), skidb_env('GAF.DIR')))

#    lapply(names(gafs), function(x) saveRDS(gafs[[x]], paste(skidb_env('GAF.DIR'), "/", x, ".rds", sep = "")))

    if (verbose)
      print(sprintf('Making genomic ranges version of gaf genes and dumping as rds file to %s/gene.gr.rds', skidb_env('GAF.DIR')))

    rownames(gafs$gene) = gafs$gene$FeatureID
    grl = gaf2gr(gafs$gene)

    saveRDS(grl, paste(skidb_env('GAF.DIR'), "/gene.gr.rds", sep = ""))

    if (verbose)
      print(sprintf('Making genomic ranges version of gaf transcript and dumping as rds file to %s/transcript.gr.rds', skidb_env('GAF.DIR')))

    grl = gaf2gr(gafs$transcript[gafs$transcript$CompositeType == 'genome', ])

    names(grl) = gafs$transcript[gafs$transcript$CompositeType == 'genome', ]$FeatureID
    values(grl)$gene = gsub('\\|.*', '', gafs$transcript[gafs$transcript$CompositeType == 'genome', ]$Gene)

    saveRDS(grl, paste(skidb_env('GAF.DIR'), "/transcript.gr.rds", sep = ""))

    if (verbose)
      print('done')
  }


#' @name read_gencode
#' @title read_gencode
#'
#' @description
#'
#' reads GENCODE file dump from .rds of GRanges.
#'
#' If that file doesn't exist
#'
#' @author Marcin Imielinski
#' @export
read_gencode = function(type = NULL, by = NULL, fn.rds = skidb_env('GENCODE.FILE.HG19.GR'), fn = GENCODE.FILE.HG19)
    {
        TYPES = c('exon', 'gene', 'transcript', 'CDS')
        BY = c('transcript_id', 'gene_id')

        if (!is.null(by))
            by = toupper(by)

        if (!is.null(type))
            type = toupper(type)

        if (file.exists(fn.rds))
            ge = readRDS(fn.rds)
        else
            {
                require(gUtils)
                cat('Importing read_gencode file from', GENCODE.FILE.HG19, 'and writing to', skidb_env('GENCODE.FILE.HG19.GR'))
                ge = import.ucsc(GENCODE.FILE.HG19)
                saveRDS(ge, skidb_env('GENCODE.FILE.HG19.GR'))
            }

        if (!is.null(type))
            {
                type = grep(type, TYPES, value = TRUE, ignore.case = TRUE)[1]
                if (!all(type %in% TYPES))
                    stop(sprintf('Type should be in %s', paste(TYPES, collapse = ',')))
                tx = ge[ge$type %in% 'transcript']
                ge = ge[ge$type %in% type]

                if (type == 'CDS' & is.null(by))
                    return(.gencode_transcript_split(ge, tx))
            }

        if (!is.null(by))
            {
                by = grep(by, BY, value = TRUE, ignore.case = TRUE)[1]

                if (!(by %in% BY))
                    stop(sprintf('Type should be in %s', paste(TYPES, collapse = ',')))

                if (by == 'transcript_id')
                    return(.gencode_transcript_split(ge, tx))
                else
                    return(split(ge, values(gene)[, by]))
            }
        return(ge)
    }


#' @name .gencode_transcript_split
#' @rdname gencode_transcript_split
#' @title .gencode_transcript_split
#'
#' @description
#'
#' splits gencode gr into transcript, taking care of some junky issues in the meantime
#'
#' @author Marcin Imielinski
.gencode_transcript_split = function(gr, tx)
    {
        require(Hmisc)
#        gr.span = seg2gr(gr2dt(gr)[, list(seqnames = seqnames[1], start = min(start), end = max(end), strand = strand[1], transcript_name = transcript_name[1], gene_name = gene_name[1], gene_status = gene_status[1], gene_type = gene_type[1]), by = transcript_id], seqlengths = seqlengths(gr))

        ## keep.ix = !is.na(values(gr.span)[, by.og])
        ## gr.span = gr.span[keep.ix,]
        ## tmp.id = paste(values(gr.span)[, by.og], seqnames(gr.span), sep = '_')
        ## if (any(ix <- duplicated(values(gr.span)[, by.og])))
        ##     values(gr.span)[, by.og][ix] = tmp.id[ix]

        ## gr = gr[!is.na(values(gr)[, by.og])]
        ## tmp.id.gr = paste(values(gr)[, by.og], seqnames(gr), sep = '_')
        ## gr$transcript_id = values(gr.span)[, by.og][match(tmp.id.gr, tmp.id)]

        gr$start.local = gr2dt(gr)[, id := 1:length(gr)][, tmp.st := 1+c(0, cumsum(width)[-length(width)]), by = transcript_id][, keyby = id][, tmp.st]

        gr$end.local = gr$start.local + width(gr) -1
        grl = split(gr, gr$transcript_id)
        tmp.val = as.data.frame(tx)[match(names(grl), tx$transcript_id), ]
        rownames(tmp.val) = tmp.val$transcript_id
        ## rownames(tmp.val) = tmp.val$splitkey
        ## names(tmp.val)[match('splitkey', names(tmp.val))] = by
        names(tmp.val) = capitalize(names(tmp.val))
        values(grl) = tmp.val
        return(grl)
    }


#' @name read_gaf
#' @title read_gaf
#'
#' @description
#'
#' reads any one of the gaf rds files stored in skidb_env('GAF.DIR') by process_gaf
#'
#' @author Marcin Imielinski
#' @export
read_gaf = function(type = 'gene', gr = F  # choices include  "transcript", "componentExon", "compositeExon", "junction", "miRNA", "pre-miRNA"
  )
  {
    if (gr)
      {
        LEGAL.TYPES = c('transcript', 'gene', 'miRNA')
        if (type %in% LEGAL.TYPES)
          return(readRDS(paste(skidb_env('GAF.DIR'), '/', type, '.gr.rds', sep = "")))
        else
          stop(sprintf('GRanges only supported for: %s', paste(LEGAL.TYPES, collapse = ', ')))
      }

    out = readRDS(paste(skidb_env('GAF.DIR'), '/', type, '.rds', sep = ""))
    out$Gene = sapply(strsplit(out$Gene, '\\|'), function(x) if (length(x)>0) x[[1]] else NA)

    if (type == 'gene')
      {
        rownames(out) = out$FeatureID;
        txt = paste(gsub('(\\w+)\\:(\\d+)\\-(\\d+)\\:([\\+\\-]$)', '\\1\t\\2\t\\3\t\\4', out$GeneLocus, perl = T), collapse = "\n")
        con = textConnection(txt);
        gene.info = read.delim(con, strings = F, header = F, col.names = c('chr', 'pos1', 'pos2', 'strand'));
        gene.info$chr = gsub('chr', '', gene.info$chr)
        out = cbind(out, gene.info);
        close(con)
      }

    return(out)
  }

#' @name read_gaf
#' @title read_gaf
#'
#' @description
#'
#' returns GRanges of genome coordinates of gaf features or gChain mapping feature coordinates to composite coordinates
#'
#' @author Marcin Imielinski
#' @export
read_gaf = function(file = 'gene', ## can be a path or also a keyword corresponding to any gaf feature eg gene, componentExon, junction, in which case the directory is queried for the appropriate file
  gchain = F)
  {
    GAF.COLNAMES = c('EntryNumber', 'FeatureID', 'FeatureType', 'FeatureDBSource', 'FeatureDBVersion', 'FeatureDBDate', 'FeatureSeqFileName', 'Composite', 'CompositeType', 'CompositeTypeDBSource', 'CompositeDBVersion', 'CompositeDBDate', 'AlignmentType', 'FeatureCoordinates', 'CompositeCoordinates', 'Gene', 'GeneLocus', 'FeatureAliases', 'FeatureInfo'); mirna.genome = read.delim('~/DB/GAF/miRNA.genome.v3_0.gaf', strings = F, header = F, col.names = GAF.COLNAMES)
    mirna.genome = read.delim(filep, strings = F, header = F, col.names = GAF.COLNAMES)
    return(mirna.genome)
  }

#' @name read_ccds_fa
#' @title read_ccds_fa
#
#' @description
#'
#' reads ccds sequences
#'
#' @author Marcin Imielinski
#' @export
read_ccds_fa = function()
  {
    ccds.fa = readDNAStringSet(skidb_env('CCDS.FA.FILE'))
    ccds.nm = sapply(strsplit(names(ccds.fa), '\\|'), function(x) x[[1]]);
    ix = !duplicated(ccds.nm)
    ccds.fa = ccds.fa[ix]
    names(ccds.fa) = ccds.nm[ix]
    return(ccds.fa)
  }


#' @name read_hg
#' @title read_hg
#'
#' @description
#' Wrapper around BSgenome call
#'
#' Retreives either the BSgenome hg18 or hg19 genome by default.  Requires packages
#' BSgenome.Hsapiens.UCSC.hg19 for hg19 and BSgenome.Hsapiens.UCSC.hg19 for hg18.
#'
#' If fft = TRUE, can also also return the hg19 ffTrack (requires that the file exists)
#' Requires the existence of environment variable HG.FFT pointing to ffTrack .rds file..
#'
#' @param hg19 Logical whether to return hg18 or hg19 BSgenome. Default TRUE
#' @param fft Logical whether to return an ffTrack. Default FALSE
#' @return BSGenome or ffTrack of the genome
#' @export
read_hg = function(hg19 = T, fft = F)
  {
      if (fft)
          {
              if (file.exists(Sys.getenv('HG.FFT')))
                  REFGENE.FILE.HG19.FFT = Sys.getenv('HG.FFT')
              else if (file.exists('~/DB/ffTracks/hg19.rds'))
                  REFGENE.FILE.HG19.FFT = '~/DB/ffTracks/hg19.rds'
              else
                  stop("Need to supply environment variable to FFtracked genome or load BSGenome. Env Var: HG.FFT")

              return(readRDS(REFGENE.FILE.HG19.FFT))
          }
    else
        {
            require(BSgenome)
            if (hg19)
                library(BSgenome.Hsapiens.UCSC.hg19)
            else
                library(BSgenome.Hsapiens.UCSC.hg18)
            return(Hsapiens)
        }
  }


#' @name read_context
#' @title read_hg_context
#'
#' @description
#'
#' reads hg context file
#'
#' @author Marcin Imielinski
#' @export
read_hg_context = function()
    {
      return(readRDS(skidb_env('HG19.CONTEXT.FFT')))
  }


#' @name gaf2gr
#' @title gaf2gr
#'
#' @description
#'
#' converts one or more gaf data frame rows to GRanges list (works only for gaf rows with composite type genome) using "CompositeCoordinates"
#' named by FeatureID
#'
#' @author Marcin Imielinski
#' @export
gaf2gr = function(gaf, strip.chr = T)
{
  ix = gaf$CompositeType == 'genome'
  if (any(!ix))
    {
      warning(sprintf('Removing %s rows of input gaf that are not of CompositeType genome', length(which(ix))))
      gaf = gaf[ix, ]
    }

  if (nrow(gaf)==0)
    stop('No rows supplied')

  gaf$CompositeCoordinates = gsub('UNKNOWN;', '', gaf$CompositeCoordinates)
  tmp = lapply(strsplit(gsub('[^\\:]\\-', ':', gaf$CompositeCoordinates), ';'), function(x) strsplit(x, ':'))
  gene.info = as.data.frame(matrix(unlist(tmp), ncol = 4, byrow = T, dimnames = list(NULL, c('chr', 'pos1', 'pos2', 'strand'))), stringsAsFactors = F)

  gaf$GeneLocus = gsub('UNKNOWN;', '', gaf$GeneLocus)
  tmp = lapply(strsplit(gsub('[^\\:]\\-', ':', gaf$CompositeCoordinates), ';'), function(x) strsplit(x, ':'))
  gene.info = as.data.frame(matrix(unlist(tmp), ncol = 4, byrow = T, dimnames = list(NULL, c('chr', 'pos1', 'pos2', 'strand'))), stringsAsFactors = F)

  if (strip.chr)
    gene.info$chr = gsub('chr', '', gene.info$chr)

  exon.pos = lapply(strsplit(gene.info$pos, ','), function(x) t(sapply(strsplit(x, '-'), as.numeric)))

  grl = do.call('GRangesList',
    lapply(1:nrow(gene.info),
           function(x)  GRanges(seqnames = Rle(gene.info$chr[x], nrow(exon.pos[[x]])),
                                strand = gene.info$strand[x], seqlengths = hg_seqlengths(),
                                ranges = IRanges(start = exon.pos[[x]][,1], end = exon.pos[[x]][,2]))))

  if (!is.null(rownames(gaf)))
    names(grl) = rownames(gaf)
  else
    names(grl) = gaf$FeatureID

  return(grl)
}

#' @name rg2exons
#' @title rg2exons
#'
#' extracts a df list of exons from one or more transcripts of the rg data frame
#' reorient = T, numbers exons relative to transcript start.
#' @description
#'
#' reads knownGene and refGene and crosses the two (left join) to yield an uber table
#'
#' @author Marcin Imielinski
rg2exons = function(rg, reorient = FALSE,
  cds = FALSE, # if TRUE will map start and end to starts / ends of cds (otherwise exon boundary)
  gr = FALSE # if TRUE will return granges list instead of df
  )
  {
    exon_starts = strsplit(rg$exon_starts, ',');
    exon_ends = strsplit(rg$exon_ends, ',');
    exon_frames = strsplit(rg$exon_frame, ',');

    exon_starts = lapply(exon_starts, as.numeric)
    exon_ends = lapply(exon_ends, as.numeric)
    exon_ix = lapply(1:length(exon_starts), function(x) rep(x, length(exon_starts[[x]])))
    exon_id = lapply(1:length(exon_starts), function(x) if(rg$strand[x][1] == "-") rev(1:length(exon_starts[[x]])) else 1:length(exon_starts[[x]]))
    cds_starts = lapply(1:length(exon_starts), function(x) pmax(exon_starts[[x]]+1, rg$s2[x]+1));  ## pesky 0 coordinates again ..
    cds_ends = lapply(1:length(exon_ends), function(x) pmin(exon_ends[[x]], rg$e2[x]));

    if (reorient)
      {
        cds_starts = lapply(1:length(cds_starts), function(x) cds_starts[[x]]-rg$s1[x]+1);
        cds_ends = lapply(1:length(cds_ends), function(x) cds_ends[[x]]-rg$s1[x]+1);
        exon_starts = lapply(1:length(exon_starts), function(x) exon_starts[[x]]-rg$s1[x]+1);
        exon_ends = lapply(1:length(exon_ends), function(x) exon_ends[[x]]-rg$s1[x]+1);
      }

    # modifying exon_start coordinate here to 1 based (prob should do it above but don't want to screw things up)
    out = cbind(rg[unlist(exon_ix), setdiff(names(rg), c('exon_starts', 'exon_ends', 'exon_frame'))],
      data.frame(exon_id = as.numeric(unlist(exon_id)), exon_start = as.numeric(unlist(exon_starts))+1, exon_end = as.numeric(unlist(exon_ends)),
                 cds_start = as.numeric(unlist(cds_starts)), cds_end = as.numeric(unlist(cds_ends)), exon_frame = as.numeric(unlist(exon_frames)), stringsAsFactors = F))

    if (cds)
      {
        out$start = out$cds_start
        out$end = out$cds_end
        out = out[(out$end-out$start)>=0, ];  ## cds may induce zero ranges in UTR's
      }
    else
      {
        out$start = out$exon_start
        out$end = out$exon_end
      }

    if (gr)
      {
        out = split(seg2gr(out), out$NM_id);
        values(out)$hugo = sapply(out, function(x) values(x)$gene_sym[1])
        if (!is.null(values(out[[1]])$Uniprot))
          values(out)$uniprot = sapply(out, function(x) values(x)$Uniprot[1])
      }

    return(out)
}


#' @name read_entrez
#' @title read_entrez
#'
#' @description
#'
#' reads knownGene and refGene and crosses the two (left join) to yield an uber table
#'
#' @author Marcin Imielinski
#' @export
read_entrez = function()
{
   eg = read.delim(skidb_env('ENTREZ.GENE.FN'), row.names=NULL, header=0, strings = F )
   colnames(eg) = c("tax_id", "GeneID", "Symbol", "LocusTag", "Synonyms", "dbXrefs", "chr", "MapLocation", "Description", "GeneType", "OfficialSymbol", "OfficialFullName", "NomenclatureStatus", "OtherDesignation", "ModificationDate")
   eg = eg[!duplicated(eg$Symbol), ];
   eg = eg[!is.na(eg$Symbol), ]
   rownames(eg) = eg$Symbol;
   return(eg)
}

#' @name read_msigdb
#' @title read_msigdb
#'
#' @description
#'
#' reads knownGene and refGene and crosses the two (left join) to yield an uber table
#'
#' @author Marcin Imielinski
#' @export
read_msigdb = function(path = skidb_env('MSIGDB_PATH'))
{
  tmp = sapply(readLines(path), function(x) gsub(' ', '', strsplit(gsub('///', '\t', x), '\t')[[1]]))
  col = vector(mode = 'character', length = length(tmp));
  out = list(flat = data.frame(name = col, url = col, genes = col, stringsAsFactors = F), lists = list());
  for (i in 1:length(tmp))
    {
      out$flat[i, 'name'] = tmp[[i]][1];
      out$flat[i, 'url'] = tmp[[i]][2];
      out$flat[i, 'genes'] = paste(tmp[[i]][3:length(tmp[[i]])], collapse = ", ");
      out$lists[i] = list(tmp[[i]][3:length(tmp[[i]])]);
    }
  names(out$lists) = out$flat$name;
  return(out)
}


#' @name read_gmt
#' @title read_gmt
#'
#' @description
#'
#' reads "GMT" format eg reactome
#'
#' @author Marcin Imielinski
#' @export
read_gmt = function(fn)
{
    lines = strsplit(readLines(fn), '\t')
    out = lapply(lines, function(x) x[3:length(x)])
    names(out) = sapply(lines, function(x) x[1])
    return(out)
}


#' @name read_txdb
#' @title read_txdb
#'
#' @description
#'
#' reads UCSC based txdb in GenomicFeatures directory
#'
#' @author Marcin Imielinski
#' @export
read_txdb = function(txdbname = 'hg19.knownGene')
  {
    fn = paste(skidb_env('GENOMIC.FEATURES.ROOT'), txdbname, ".sqlite", sep = "")

    if (!file.exists(fn))
      stop(sprintf('DB not found, available DBs include:\n %s. \nUse make_txdb if desired DB not available',
                   paste(gsub('\\.sqllite', '', dir(skidb_env('GENOMIC.FEATURES.ROOT'))), collapse = ", ")))

    return(loadFeatures(fn))
  }


#' @name make_txdb
#' @title make_txdb
#'
#' @description
#'
#' reads UCSC based txdb in GenomicFeatures directory
#'
#' @author Marcin Imielinski
#' @export
make_txdb = function(db = 'hg19', tablename = 'knownGene')
  {
    txdb = makeTranscriptDbFromUCSC(genome = 'hg19', tablename = 'knownGene')
    saveFeatures(txdb, paste(skidb_env('GENOMIC.FEATURES.ROOT'), "/", db, ".", tablename, ".sqlite", sep = ""))
  }

#' @name blastp
#' @title blastp
#'
#' @description
#'
#' Run blastp (via command line) on (single) sequence against uniprot DB (default) or arbitrary DB
#'
#' @author Marcin Imielinski
#' @export
blastp = function(seq, # can be either a character string of amino acids or a Biostrings Xstring
  db = skidb_env('UNIPROT.FA'))
{
  require(Biostrings)

  if (is.character(seq))
    seq = AAString(seq)
  else if (inherits(tmp, 'XString'))
    seq = as(seq, 'XStringSet');

  if (length(seq)>1)
    seq = seq[1];

  tmp.fasta.path = paste(skidb_env('TMP.DIR'), 'blastp.marcin.', round(runif(1)*1e15), '.fasta', sep = "")
  write.XStringSet(seq, tmp.fasta.path)

  blast.cmd = paste(skidb_env('BLAST.DIR'), '/blastp', ' -query ', tmp.fasta.path,  ' -outfmt 3 ', ' -db ', skidb_env('UNIPROT.FA'), sep = "")

  p = pipe(blast.cmd);
  txt = readLines(p)
  close(p)
  blank.ix = which(txt=="")
  nonblank.ix = which(txt!="")

  # assuming that the query results stats are 2 lines from res.ix
  res.start.ix = grep('Sequences producing significant alignments', txt)+2;  blank.ix = which(txt=="");
  res.end.ix = min(blank.ix[blank.ix>res.start.ix]);

# assuming that the query results stats are 2 lines from res.ix
  res.start.ix = grep('Sequences producing significant alignments', txt)+2;  blank.ix = which(txt=="");
  res.end.ix = min(blank.ix[blank.ix>res.start.ix])-1;

  #ugly fixed width processing
  w = get.field.widths(txt[res.start.ix])
  w = c(sum(w[1:(length(w)-2)]), w[(length(w)-1):length(w)])

  out = as.data.frame(matrix(unlist(lapply(txt[res.start.ix:res.end.ix], function(x) trim(strsplit.fwf(x, w)))), ncol = 3, byrow = T), stringsAsFactors = F)
  names(out) = c('match', 'score', 'E')
  out$E = as.numeric(out$E)
  out$score = as.numeric(out$score)
  out$protein =  gsub('^[^\\|]+\\|[^\\|]+\\|(\\w+)+ .*$', '\\1', out$match)

  #more ugly fixed width processing
  align.start.ix = min(nonblank.ix[nonblank.ix>res.end.ix])
  align.end.ix = min(blank.ix[blank.ix > align.start])-1
  w = get.field.widths(txt[align.start.ix])
  w = c(sum(w[1:(length(w)-2)]), w[(length(w)-1):length(w)])

  out = cbind(out, as.data.frame(matrix(unlist(lapply(txt[align.start.ix:align.end.ix], function(x) trim(strsplit.fwf(x, w)))), ncol = 3, byrow = T), stringsAsFactors = F))
  names(out) = c('match', 'score', 'E')
  out$E = as.numeric(out$E)
  out$score = as.numeric(out$score)
  out$protein =  gsub('^[^\\|]+\\|[^\\|]+\\|(\\w+)+ .*$', '\\1', out$match)

  text.con = textConnection(gsub(txt[res.start.ix:res.end.ix]))
  con = read.delim(text.con, strings = F, header = F)
  close(text.con)
}

#' @name tcga_bcr2four
#' @title tcga_bcr2four
#'
#' @description
#'
#' Extracts "four digit" code from tcga "bcr" or "dcc" style barcode
#'
#' @author Marcin Imielinski
#' @export
tcga_bcr2four = function(x)
  {
    ix = grep('TCGA', x)
    x[ix] = gsub('.*(TCGA.*)', '\\1', x[ix])

    return(mapply(function(x,y) if (length(x)>2) if (nchar(x[[1]]) == 4 & nchar(x[[2]]) == 2 & nchar(x[[3]])==4) x[[3]] else y else y, strsplit(gsub('[\\_\\]', '-', x), '\\-'), x))
  }

#' @name get_ucsc_chroms
#' @title get_ucsc_chroms
#'
#' @description
#'
#' returns segments corresponding to coordinates of ucsc hg18 / hg19 chromosomesfour
#'
#' "clean" will return without all of the unmapped segments (e.g Un_gl000221)
#' reads knownGene and refGene and crosses the two (left join) to yield an uber table
#'
#' @author Marcin Imielinski
#' @export
get_ucsc_chroms = function(hg19 = T, clean = T, numeric = F, add.chr = F)
{
  if (hg19)
    fn = paste(skidb_env('UCSC.DIR'), 'hg19.chrom.sizes.txt', sep = "/")
  else
    fn = paste(skidb_env('UCSC.DIR'), 'hg18.chrom.sizes.txt', sep = "/")

  chroms = read.delim(fn, strings = F, header = F);
  names(chroms) = c('chr', 'pos2')
  chroms$pos1 = 1;

  if (clean)
    chroms = chroms[!grepl('(random)|(Un\\_)|(hap)', chroms$chr), ]

  if (add.chr)
    chroms$chr = paste('chr', chroms$chr, sep = "")

  if (numeric)
    chroms$chr = chr2num(chroms$chr)

  return(chroms[, c('chr', 'pos1', 'pos2')])
}

#' @name get_ucsc_bands
#' @title get_ucsc_bands
#'
#' @description
#' returns segments corresponding to coordinates of ucsc hg18 / hg19 chrom bands
#'
#'
#'
#' @author Marcin Imielinski
#' @export
get_ucsc_bands = function(hg19 = T, arms = F, add.chr = F, gr = T)
  {
    if (hg19)
      fn = paste(skidb_env('UCSC.DIR'), 'hg19.cytoband.txt', sep = "/")
    else
      fn = paste(skidb_env('UCSC.DIR'), 'hg18.cytoband.txt', sep = "/")

    bands = read.delim(fn, strings = F)
    names(bands) = c('chr', 'pos1', 'pos2', 'name', 'stain')

    if (add.chr == F)
      bands$chr = gsub('chr', '', bands$chr)

    bands$name = paste(bands$chr, bands$name, sep = "")

    if (arms)
      {
        bands$name = gsub('([pq]).*', '\\1', bands$name)
        tmp.pos1 = aggregate(bands$pos1, by = list(bands$name), FUN = min)
        tmp.pos2 = aggregate(bands$pos2, by = list(bands$name), FUN = max)
        tmp.chr = aggregate(bands$chr, by = list(bands$name), FUN = unique)
        bands = data.frame(name = tmp.pos1[,1], chr = tmp.chr[,2],
          pos1 = tmp.pos1[,2], pos2 = tmp.pos2[,2], stringsAsFactors = F)
      }

    bands = bands[order(chr2num(bands$chr), bands$pos2), ]
    rownames(bands) = bands$name;

    if (gr)
      {
        sl = aggregate(formula = pos2 ~ chr, data = bands, FUN = max);
        bands = GRanges(seqnames = bands$chr, IRanges(pmax(1, bands$pos1), bands$pos2), band = bands$name, stain = bands$stain, seqlengths = structure(sl[,2]+1, names = sl[,1]))
      }

    return(bands)
  }


#' @name muscle
#' @title muscle
#'
#' @description
#'
#' runs muscle program on XStringSet (e.g. AAStringSet, DNAStringSet) object "seq" and
#' return output as XStringSet
#'
#' reads knownGene and refGene and crosses the two (left join) to yield an uber table
#'
#' @author Marcin Imielinski
#' @export
muscle = function(seq, sub.mat = NULL, gap.open = -10, gap.extend = -4, html.out = NULL)
{
  tmp.name = paste(skidb_env('TMP.DIR'), '/muscle_tmp', gsub('0\\.', '', as.character(runif(1))), sep = '')
  tmp.name.out = paste(tmp.name, '.out', sep = '')
  writeXStringSet(seq, tmp.name)
  if (!is.null(html.out))
  {
    if (nchar(html.out)>0)
    {
      cat('Writing html results to ', html.out, '.\n')
      html.out = paste('-htmlout', html.out)
    }
  } else
    html.out = '';

  cmd = paste(skidb_env('MUSCLE.PATH'), '-in', tmp.name, '-fastaout', tmp.name.out, html.out)
  if(!is.null(sub.mat))
    {
      tmp.name.BLOSUM = paste(tmp.name, ".BLOSUM", sep = "")
    sub.mat = rbind(dimnames(sub.mat)[[2]], sub.mat)
    write.fwf(as.data.frame(sub.mat), tmp.name.BLOSUM, rownames = T, colnames = F)
    cmd = paste(cmd, "-matrix", tmp.name.BLOSUM, "-gapopen", gap.open, "-gapextend", gap.extend)
  }
  cmd = paste(cmd, '2>/dev/null')
  system(cmd)

  if (inherits(seq, 'AAStringSet'))
    return(readAAStringSet(tmp.name.out, format = 'fasta'))
  else
    return(readDNAStringSet(tmp.name.out, format = 'fasta'))
}

#' @name run_prego
#' @title run_prego
#'
#' @description
#'
#' writes RA and CN output to directory, so it can be run from command line
#'
#' intervals = GRanges tiling of genome with $tum.cov and $norm.cov fields (obtain using bam.cov.gr)
#' ra = output of ra_breaks, GRangesList of signed breakpoint pairs pointing in the direction of the joined segment
#'
#' @author Marcin Imielinski
#' @export
read_prego = function(fn)
  {
    res.fn = paste(fn, '.results', sep = '')
    res.tmp = readLines(res.fn)
    res = structure(lapply(split(res.tmp, cumsum(grepl('edges', res.tmp))), function(x) read.delim(textConnection(x), strings = F, skip = 1, header = F,
      col.names = c('node1', 'chr1', 'pos1', 'node2', 'chr2', 'pos2', 'cn'))), names = gsub(':', '', grep('edges', res.tmp, value = T)))
    res[[1]]$tag = paste(res[[1]]$node1, ':', res[[1]]$node2, sep = '')

    out = list()
    segstats = GRanges(res[[1]]$chr1, IRanges(res[[1]]$pos1, res[[1]]$pos2), strand = '+',
      cn = res[[1]]$cn, left.tag = res[[1]]$node1, right.tag = res[[1]]$node2)
    segstats = gr.fix(c(segstats, gr.flip(segstats)))

    neg.ix = which(strand(segstats) == '-')
    tag1 = segstats$right.tag;
    tag1[neg.ix] = segstats$left.tag[neg.ix]
    tag2 = segstats$left.tag;
    tag2[neg.ix] = segstats$right.tag[neg.ix]

    out$adj.cn = matrix(0, nrow = length(segstats), ncol = length(segstats), dimnames = list(tag1, tag2))
    out$adj.cn[cbind(res[[2]]$node1, res[[2]]$node2)] = res[[2]]$cn
    out$adj.cn[cbind(res[[2]]$node2, res[[2]]$node1)] = res[[2]]$cn
    out$adj.cn[cbind(res[[3]]$node1, res[[3]]$node2)] = res[[3]]$cn
    out$adj.cn[cbind(res[[3]]$node2, res[[3]]$node1)] = res[[3]]$cn

    out$adj.type = matrix('', nrow = length(segstats), ncol = length(segstats), dimnames = list(tag1, tag2))
    out$adj.type[cbind(res[[2]]$node1, res[[2]]$node2)] = 'ref'
    out$adj.type[cbind(res[[2]]$node2, res[[2]]$node1)] = 'ref'
    out$adj.type[cbind(res[[3]]$node1, res[[3]]$node2)] = 'ab'
    out$adj.type[cbind(res[[3]]$node2, res[[3]]$node1)] = 'ab'
    out$segstats = segstats;

    return(out)
  }

#' @name read_ucsc_files_table
#' @title read_ucsc_files_table
#'
#' @description
#' parses contents of url pointing to files.txt file on UCSC server
#' corresponding to track metadata and the corresponding files (which are
#' assumed to lie in the same directory
#'
#' e.g. "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/files.txt"#'
#' @author Marcin Imielinski
#' @export
read_ucsc_files_table = function(url)
  {
    files = read.delim(url, strings = F, header = F)
    names(files) = c('path', 'metadata')
    tmp = lapply(strsplit(files$metadata, ';'), function(x) strsplit(trim(x), '='))
    tmp2 = lapply(tmp, function(x) structure(sapply(x, function(x) x[[2]]), names = sapply(x, function(x) x[[1]])))
    all.names = unique(unlist(sapply(tmp2, function(x) names(x))))
    base.path = gsub('\\/[^\\/]+$', '', url)
    out = data.frame(url = paste(base.path, files$path, sep = '/'), file = files$path, stringsAsFactors = F)

    out[, all.names] = NA
    for (i in 1:nrow(out))
      out[i, names(tmp2[[i]])] = tmp2[[i]]

    return(out)
  }

#' @name read_hicdomain
#' @title read_hicdomain
#'
#' @description
#'
#' reads hic domain output from Rao .. Lieberman Cell 2014
#'
#' @author Marcin Imielinski
#' @export
read_hicdomain = function(fn)
{
    dat = read.delim(fn, strings = FALSE)

    gr1 = dat[, 1:3]
    gr2 = dat[, 3 + 1:3]
    names(gr1) = names(gr2) = c('chr', 'start', 'end')
    val = dat[, -c(1:6)]
    out = grl.pivot(GRangesList(seg2gr(gr1), seg2gr(gr2)))
    values(out) = val
    return(out)
}

#' @name wig2bw
#' @title wig2bw
#'
#' @description
#'
#' Runs UCSC command tool wig2big on a set of wig files (files.wig) and outputs to the corresponding bigwig files (files.bw)
#' wig files can be zipped as well
#'
#' If files.bw is not provided then will just replace wig with big wig in files.wig name)
#'
#' @author Marcin Imielinski
#' @export
wig2bw = function(files.wig, files.bw = NULL,
  chrom.path = skidb_env('UCSC.CHROM.SIZES.PATH'),
  fork = 0 #if greater than >0 then will mclapply across "fork" parallele processes
  )
{
  require(parallel)
  if (is.null(files.bw))
    files.bw = paste(files.wig, '.bw', sep = "")

  if (fork>0)
    {
      require(parallel)
      fork = min(fork, 20) ; # can't get too crazy

      mclapply(1:length(files.wig),
             function(x) system(sprintf('wigToBigWig %s %s %s', files.wig[x],chrom.path, files.bw[x])), mc.cores = fork)
    }
  else
    sapply(1:length(files.wig),
           function(x) system(sprintf('wigToBigWig %s %s %s', files.wig[x],chrom.path, files.bw[x])))
}

#' @name bw2rle
#' @title bw2rle
#'
#' @description
#'
#' convert each big wig (bw) file in files.bw into an RDS file storing a SimpleRleList
#' across list of regions (GRanges or seg data frame)
#'
#' if regions is NULL will process all data stored in the bw file (warning can be large)
#'
#' These files can be loaded with readRDS
#'
#' @author Marcin Imielinski
#' @export
bw2rle = function(files.bw, regions = NULL, files.rle = NULL,
  precision = 2, # decimal precision which to keep (current limitation of import.bw), maximum is ~8
  fork = 0, ## if true will MCL apply over (max 20) processors
  return.objects = TRUE, # if false, will return only filenames
  dump.files = TRUE, # if false, will not dump rle Rdata files
  verbose = FALSE
  )
  {
    require(parallel)
    if (!is.null(regions))
      if (inherits(regions, 'data.frame')) ## we try to convert to genomic ranges if we see data frame
        regions = seg2gr(regions)

    if (is.null(files.rle))
      files.rle = paste(files.bw, '.rle.rds', sep = "")

    get.fun = function(x)
      {
        if (file.exists(files.bw[x]))
          {
            if (verbose == T)
              print(paste("Processing", files.bw[x]))

            tmp = import.bw2(files.bw[x], regions, fill = NA, output.rle = T);

            if (dump.files)
              saveRDS(tmp, file = files.rle[x])
            if (!return.objects) # only return the file name
              tmp = files.rle[x]

            gc()
          }
        else
          tmp = NA
        return(tmp);
      }

    if (fork==0)
      out = sapply(1:length(files.bw), get.fun)
    else
      {
        require(parallel)
        fork = min(fork, 20) ; # can't get too crazy
        out = mclapply(1:length(files.bw), get.fun, mc.cores = fork)
      }

    names(out) = files.bw;

    return(out)
  }

#' @name wig2rlea
#' @title wig2rle
#'
#' @description
#'
#' convert each big wig (bw) file in files.bw into an RDS file storing a SimpleRleList
#' across list of regions (GRanges or seg data frame)
#'
#' These files can be loaded with readRDS
#'
#' @author Marcin Imielinski
#' @export
wig2rle = function(files.wig, files.rle,
  fork = 0, ## if true will MCL apply over (max 20) processors
  return.objects = TRUE, # if false, will return only filenames
  dump.files = TRUE, # if false, will not dump rle Rdata files
  verbose = FALSE
  )
  {
    if (inherits(regions, 'data.frame')) ## we try to convert to genomic ranges if we see data frame
      regions = seg2gr(regions)

    if (is.null(files.rle))
      files.rle = paste(files.wig, '.rle.rds', sep = "")

    BLOCK.SIZE = 1e5;

    .wig2rle = function(x)
      {
        fh = file(files.wig[x])
        open(fh)

        chunk = readLines(fh, n = BLOCK.SIZE)
        chunk = chunk[-grep('^track type', chunk)];

        while (length(chunk)>0)
          {
            step.ix = grep('(fixedStep)|(variableStep)', chunk)

#            step.names = sapply(strsplit(gsub('^fixedStep ', '', bla[ix]), ' '), function(x) sapply(strsplit(x, '='), function(y) y[1]))
            step.vals = sapply(strsplit(gsub('^fixedStep ', '', bla[ix]), ' '), function(x) sapply(strsplit(x, '='), function(y) y[2]))


          }
      }

  }


##########################
#' @name bam2bw
#' @title bam2bw
#'
#' @description
#'
#'
#' Runs BEDTools genomeCoverageBed and then UCSC tool BedGraphToWigWig
#' on a set of bam files to output per base coverage in bigwig format.
#'
#' If files.bw is not provided then will just replace .bam with .bedtools.coverage.bw in files.wig name)
#'
#' NOTE: assumes that BedGraphToBigWig and genomeCoverageBed already exist on file system
#' and skidb_env('UCSC.CHROM.SIZES.PATH') should be valid.   Also assumes bam is sorted (ie samtools sort bam)
#'
#' reads knownGene and refGene and crosses the two (left join) to yield an uber table
#'
#' @author Marcin Imielinski
#' @export
bam2bw = function(files.bam, files.bw = NULL,
  fork = 0 #if greater than >0 then will mclapply across "fork" parallel processes
  )
{
  if (is.null(files.bw))
    files.bw = paste(files.bam, '.bedtools.coverage.bw', sep = "")


  .bam2bw = function(x) { cmd = paste(
                                    sprintf('genomeCoverageBed -ibam %s -g %s -bg > %s.tmp.bedgraph',
                                            files.bam[x], skidb_env('UCSC.CHROM.SIZES.PATH'), files.bw[x]),
                                    sprintf('bedGraphToBigWig %s.tmp.bedgraph %s %s', files.bw[x], skidb_env('UCSC.CHROM.SIZES.PATH'), files.bw[x]),
                          #          sprintf('rm %s.tmp.bedgraph', files.bw[x]),
                                    sep = "\n");
                         writeLines(paste("Running:", cmd));
                         system(cmd);
                        }

  if (fork>0)
    {
      require(parallel)
      fork = min(fork, 20) ; # can't get too crazy

      mclapply(1:length(files.bam), .bam2bw, mc.cores = fork)
    }
  else
    sapply(1:length(files.bam), .bam2bw)

}


#' @name wig2covered
#' @title wig2covered
#'
#' @description
#' takes wigfile of coverage with path fn (which is either a text, gz, or bz2 file) and
#' outputs all "covered" positions (chr pos) to out.fn as space separated file.
#'
#' NOTE: assumes all positions without 0 value are covered
# NOTE: assumes wig file is in fixedStep format
#' reads knownGene and refGene and crosses the two (left join) to yield an uber table
#'
#' @author Marcin Imielinski
#' @export
wig2covered = function(fn, out.fn = paste(fn, '.covered.txt', sep = ''), verbose = T,
  gz = T ## whether to gzip after dumping
  )
  {
    type = gsub('.*\\.([^\\.]+)', '\\1', fn)

    if (!(type %in% c('txt', 'wig', 'bz2', 'gz')))
      stop('Input file must be wig file in text format (extension .wig, .txt), bz2 format (extension .bz2), or gz format (extension .gz')


    if (type == 'gz')
      {
        p = pipe(paste('gunzip -c', fn, " | grep -n fixedStep | sed 's/\\(\\:fixedStep\\)\\|\\(chrom\\=\\)\\|\\(start\\=\\)\\|\\(step\\=\\)//g'"), 'r')
        fs = read.delim(p, sep = ' ', header = F, strings = F, col.names = c('line', 'chr', 'start', 'step'))
        close(p)
      }
    else if (type == 'bz2')
      {
        p = pipe(paste('bunzip2 -c', fn, " | grep -n fixedStep | sed 's/\\(\\:fixedStep\\)\\|\\(chrom\\=\\)\\|\\(start\\=\\)\\|\\(step\\=\\)//g'"), 'r')
        fs = read.delim(p, sep = ' ', header = F, strings = F, col.names = c('line', 'chr', 'start', 'step'))
        close(p)
      }
    else
      {
        p = pipe(paste("grep -n fixedStep ", fn, " | sed 's/\\(\\:fixedStep\\)\\|\\(chrom\\=\\)\\|\\(start\\=\\)\\|\\(step\\=\\)//g'"), 'r')
        fs = read.delim(p, sep = ' ', header = F, strings = F, col.names = c('line', 'chr', 'start', 'step'))
        close(p)
      }

    con.out = file(out.fn, 'w')
    l = fs$line
    chr = fs$chr
    st = fs$start
    if (type == 'gz')
      wigf = gzfile(fn, 'r')
    else if (type == 'bz2')
      wigf = bz2file(fn, 'r')
    else
      wigf = file(fn, 'r')

    line = readLines(wigf, l[1]-1)
    last = 0;
    for (i in 2:length(l))
      {
        lines = suppressWarnings(readLines(wigf, l[i]-l[i-1]))
        ix = lines!='0'
        ix[1] = FALSE
        if (any(ix))
          {
            ix = which(ix)
                                        #        out[(last + 1):(last + length(ix))] =  paste(chr[i-1], st[i-1]+ix)
            writeLines(paste(chr[i-1], st[i-1]+ix), con.out)
            last = last + length(ix)
          }
        if ((i %% 100)==0 & verbose)
          cat(i, lines[1], '\n')
      }
    close(con.out)
    close(wigf)

    if (gz)
      system(paste('gzip', out.fn))
  }

#' @name read_maf
#' @title read_maf
#'
#' @description
#'
#' read maf taking care of headers
#'
#' @author Marcin Imielinski
#' @export
read_maf = function(fn, nskip = NULL, cols = NULL, dt = FALSE, add.path = FALSE, verbose = FALSE)
    {
        if (verbose)
            cat('Reading', fn, '\n')
        if (is.null(nskip))
            {
                header = readLines(fn, 20)
                nskip = sum(grepl('^#', header))

                if (length(header)==0)
                    return(data.frame())

                if (nskip == length(header))
                    {
                        header = readLines(fn, 50)
                        nskip = sum(grepl('^#', header))
                    }

                if (nskip == length(header))
                    {
                        header = readLines(fn, 500)
                        nskip = sum(grepl('^#', header))
                    }

                if (nskip == length(header))
                    stop('Header is extra long check file')
            }

        if (dt)
            {
                out = tryCatch(fread(fn, skip = nskip), error = function(e) 'Error')
                if (is.data.table(out))
                    if (!is.null(cols))
                        out = out[, intersect(cols, names(out)), with = FALSE]
                if (add.path)
                    out$source.path = fn
            }
        else
            {
                out <- tryCatch(read.delim(fn, skip = nskip, strings = FALSE), error = function(e) 'Error')
                if (is.character(out))
                    {
                        warning('Error getting maf file')
                        return(NULL)
                    }

                if (is.data.frame(out))
                    if (!is.null(cols))
                        out = out[, cols]

                if (add.path)
                    out$source.path = fn
            }
        return(as.data.table(out))
    }


#' @name read_vcf
#' @title read_vcf
#'
#' @description
#'
#' wrapper around variatnAnnotation reads VCF into granges or data.table format
#'
#' @author Marcin Imielinski
#' @export
read_vcf = function(fn, gr = NULL, hg = 'hg19', geno = NULL, swap.header = NULL, verbose = FALSE, add.path = FALSE, tmp.dir = '~/temp/.tmpvcf', ...)
    {
        require(VariantAnnotation)
        in.fn = fn

        if (verbose)
            cat('Loading', fn, '\n')

        if (!is.null(gr))
            {
                tmp.slice.fn = paste(tmp.dir, '/vcf_tmp', gsub('0\\.', '', as.character(runif(1))), '.vcf', sep = '')
                cmd = sprintf('bcftools view %s %s > %s', fn,  paste(gr.string(gr.stripstrand(gr)), collapse = ' '), tmp.slice.fn)
                if (verbose)
                    cat('Running', cmd, '\n')
                system(cmd)
                fn = tmp.slice.fn
            }

        if (!is.null(swap.header))
            {
                if (!file.exists(swap.header))
                    stop(sprintf('Swap header file %s does not exist\n', swap.header))

                system(paste('mkdir -p', tmp.dir))
                tmp.name = paste(tmp.dir, '/vcf_tmp', gsub('0\\.', '', as.character(runif(1))), '.vcf', sep = '')
                if (grepl('gz$', fn))
                    system(sprintf("zcat %s | grep '^[^#]' > %s.body", fn, tmp.name))
                else
                    system(sprintf("grep '^[^#]' %s > %s.body", fn, tmp.name))

                if (grepl('gz$', swap.header))
                    system(sprintf("zcat %s | grep '^[#]' > %s.header", swap.header, tmp.name))
                 else
                    system(sprintf("grep '^[#]' %s > %s.header", swap.header, tmp.name))

                system(sprintf("cat %s.header %s.body > %s", tmp.name, tmp.name, tmp.name))
                vcf = readVcf(tmp.name, hg, ...)
                system(sprintf("rm %s %s.body %s.header", tmp.name, tmp.name, tmp.name))
            }
        else
            vcf = readVcf(fn, hg, ...)
        
        out = granges(vcf)

        if (!is.null(values(out)))
            values(out) = cbind(values(out), info(vcf))
        else
            values(out) = info(vcf)
        
                                    
        if (!is.null(geno))
        {
            if (geno)
                for (g in  names(geno(vcf)))
                {
                    geno = names(geno(vcf))
                    warning(sprintf('Loading all geno field:\n\t%s', paste(geno, collapse = ',')))
                }
            
            gt = NULL
            for (g in geno)
            {
                m = as.data.frame(geno(vcf)[[g]])
                names(m) = paste(g, names(m), sep = '_')
                if (is.null(gt))
                    gt = m
                else
                    gt = cbind(gt, m)
            }
            values(out) = cbind(values(out), as(gt, 'DataFrame'))
        }

        if (!is.null(gr))
            system(paste('rm', tmp.slice.fn))

        if (add.path)
            values(out)$path = in.fn

        return(out)
    }

#' @name read_10x
#' @title tally 10x fragments across windows
#' @description
#'
#' pull reads from 10x .bam across windows and return matrix barcodes x windows or window x window co-coverage matrix
#'
#' @param bam bam file
#' @param wins GRanges of windows to query and tally
#' @param cocov logical flag whether to return cocov length(win) x length(win) matrix or matrix of fragments x windows
#'
#'
read_10x = function(bam = NULL, wins, reads = NULL, cocov = TRUE, readcount = FALSE, verbose = FALSE)
{
    require(Rsamtools)
    require(data.table)

    win = streduce(wins)
    if (verbose)
        cat('Extracting reads\n')

    if (!is.null(bam))
        reads = read.bam(bam, win, tag = c('RX', 'BX'), verbose = verbose, as.data.table = TRUE, pairs.grl = FALSE)
    else if (is.null(reads))
        stop('provide either reads or bam as input')

    if (all(is.na(reads$RX)))
        reads[, RX := BX]

    if (verbose)
        cat('Finished extracting reads\n')

    # reads2 = seg2gr(reads[!is.na(MD), ][is.dup(RX), ][, span := diff(range(start)), by = RX])
    reads2 = seg2gr(reads[is.dup(RX), ][ , span := diff(range(start)), by = RX][!is.na(RX), ])
    #    reads2 = seg2gr(reads[!is.na(MD), ][, span := diff(range(start)), by = RX])

    if (length(reads2)==0)
        {
            warning('No reads found with RX flag')
            return(sparseMatrix(1, 1, x = 0, dims = c(length(wins), length(wins))))
        }

    r2 = gr2dt(reads2[, c('RX')] %*% wins)

    if (nrow(r2)==0)
        return(sparseMatrix(1, 1, x = 0, dims = c(length(wins), length(wins))))

     r2 = r2[!is.na(RX), ]

    setkey(r2, RX)

    if (readcount)
        r2[, weight := 1, by = RX]
    else
        r2[, weight := width/sum(width), by = RX]

    if (verbose)
        cat('Finished making data.table\n')

    require(Matrix)
    urx = unique(r2$RX)
    M = sparseMatrix(as.integer(factor(r2$RX), urx), r2$subject.id, x = r2$weight)
    rownames(M) = as.character(urx)

    if (verbose)
        cat('Finished tallying\n')

    if (cocov)
        {
            out = (t(M) %*% M)
            if (nrow(out) != length(wins) | ncol(out) != length(wins))
                {
                    tmp = out
                    out = sparseMatrix(1,1, x = 0, dims = c(length(wins), length(wins)))
                    out[1:nrow(tmp), 1:ncol(tmp)] = tmp
                }
            return(out)
        }
    else
        return(M)
}


#############
# run meme
#
#
#############
meme = function(sequences, basedir = './',
        tilim = 18000,
        maxsize = 60000,
        nmotifs = 8,
        minw = 3,
        maxw = 50,
        revcomp = TRUE,
        MEMEPATH = '/broad/software/free/Linux/redhat_6_x86_64/pkgs/meme_4.7.0/bin/meme')
    {

        outdir = paste(basedir, 'memeres', sep = '/')
        system(paste('mkdir -p', outdir))

        if (is.character(sequences))
            sequences = DNAStringSet(sequences)

        if (is.null(names(sequences)))
            names(sequences) = as.character(1:length(sequences))

        writeXStringSet(sequences, paste0(outdir, '/sequences.fa'))

        cmd = sprintf('cd %s; %s sequences.fa -dna -oc . -nostatus -time %s -maxsize %s -mod zoops -nmotifs %s -minw %s -maxw 50 -V', normalizePath(outdir), MEMEPATH, tilim, maxsize, nmotifs, minw, maxw)
        if (revcomp)
            cmd = paste(cmd, '-revcomp')

        cat('Running:', cmd, '\n')
        system(cmd)
        cat('MEME processing done check results at', outdir, '\n')
    }


#############
# run dreme
#
#
#############
dreme = function(case, control = NULL,
    basedir = './',
    tilim = 1000,
    mink = 3,
    maxk = 8,
    e.thresh = 0.05,
    m.thresh = NULL,
    ngen = 100,
    MEMEPATH = '/broad/software/free/Linux/redhat_6_x86_64/pkgs/meme_4.7.0/bin/dreme')
    {

        outdir = paste(basedir, 'dreme_out', sep = '/')
        system(paste('mkdir -p', outdir))

        if (is.factor(case))
            case = as.character(case)

        if (is.factor(control))
            control = as.character(control)

        if (is.character(case))
            case = DNAStringSet(case)

        if (is.null(names(case)))
            names(case) = as.character(1:length(case))

        writeXStringSet(case, paste0(outdir, '/case.fa'))

        if (!is.null(control))
            {if (is.character(control))
                 control = DNAStringSet(control)

                if (is.null(names(control)))
                    names(control) = as.character(1:length(control))
                writeXStringSet(control, paste0(outdir, '/control.fa'))
            }

        if (is.null(m.thresh))
            cmd = sprintf('cd %s; %s -p dreme_out/case.fa  -n dreme_out/control.fa -t %s -mink %s -maxk %s -g %s -v 5 -e %s', normalizePath(basedir), MEMEPATH, tilim, mink, maxk, ngen, e.thresh)
        else
            cmd = sprintf('cd %s; %s -p dreme_out/case.fa  -n dreme_out/control.fa -t %s -mink %s -maxk %s -g %s -v 5 -m %s', normalizePath(basedir), MEMEPATH, tilim, mink, maxk, ngen, m.thresh)

        cat('Running:', cmd, '\n')
        system(cmd)
        cat('DREME processing done check results at', outdir, '\n')
    }


#' @name dump_dranger_for_IGV_dump
#' @title _dranger_for_IGV_snapshots
#' @description
#'
#' dumps "snapshot_queue" file from dranger table into igv.dir for IGV screen shots of translocations
#' Takes as input:
#' dranger = data frame obtained by reading in dranger output
#' iset= firehose iset name OR data.frame of individuals table from firehose obtained by doing get_fh(...,  merged = TRUE) which contains columns $Individual_Id,
#' $Tumor_clean_bam_file_wgs, $Normal_clean_Bam_file_wgs;
#' wkspace = firehose workspace (needs to be provided if iset is a name of a firehose iset)
#'
#' (should be followed by command line call to "make snapshots" or "make snapshots_batch" via Makefile in snapshots directory to take snapshots on the farm
#' of queued mutations on the farm)
#'
#' @export
#' @author Marcin Imielinski
dump_dranger_for_IGV_snapshots = function(dranger, iset, wkspace = "Lung", igv.dir = "./", mutation.log.name = 'mutation_log', snapshots.queue.filename = 'snapshot_queue', dump.mutation.log = TRUE, chr.include = T)
{
  SNAPSHOT.SUBDIR = "images";
  DRANGER.MAKEFILE.PATH = '/home/unix/marcin/Scripts/Makefiles/Makefile.dranger';

  system(sprintf('mkdir -p %s/%s', igv.dir, SNAPSHOT.SUBDIR));

  if (is.character(iset))
    iset = get_fh(iset, wkspace = wkspace, merged = TRUE);

  dranger[, c('Tumor_bam', 'Normal_bam')] = iset[match(dranger$individual, iset$Individual_Id), c('Tumor_clean_bam_file_wgs', 'Normal_clean_bam_file_wgs')]
  dranger$label <- paste(dranger$gene1, dranger$gene2, sep = "-");
  dranger$label[dranger$gene1 == dranger$gene2] = dranger$gene1[dranger$gene1 == dranger$gene2];
  dranger$label = paste(dranger$individual, dranger$label, sep = "_")

  if (chr.include == T & !any(grepl('chr', dranger$chr1)))
    {
      dranger$chr1 = paste('chr', dranger$chr1, sep = "");
      dranger$chr2 = paste('chr', dranger$chr2, sep = "");
    }

  if (chr.include == F & any(grepl('chr', dranger$chr1)))
    {
      dranger$chr1 = gsub('chr', '', dranger$chr1)
      dranger$chr2 = gsub('chr', '', dranger$chr2)
OA    }

  dranger$snapshot_path = paste(file_path_as_absolute(paste(igv.dir, SNAPSHOT.SUBDIR, sep = "/")),
    paste(paste(dranger$label, dranger$chr1, dranger$min1, dranger$max1, dranger$chr2, dranger$min2, dranger$max2, sep = "_"), '.png', sep = ""), sep = "/");
  dranger$Variant_Classification = gsub('([^\\(]+)[\\(]?.*$', '\\1', dranger$fusion)
  write.tab(dranger[,c('individual', 'chr1', 'min1', 'max1', 'chr2', 'min2', 'max2', 'snapshot_path', 'Tumor_bam', 'Normal_bam')], paste(igv.dir, snapshots.queue.filename, sep = "/"), col.names = F)

  dranger$reviewed = "";

  # dump mutation log so that mutation_station can eventually annotate
  mutation_log = dranger[, c('label', 'chr1', 'min1', 'max1', 'chr2', 'min2', 'max2', 'Variant_Classification', 'reviewed', 'snapshot_path')]
  if (dump.mutation.log == T)
    write.table(mutation_log,
                paste(igv.dir, mutation.log.name, sep = "/"),
                quote = F, sep = "\t", row.names = F, col.names = TRUE)

  # transfer makefile to new igv snapshot directory to make images
  system(paste('cp', DRANGER.MAKEFILE.PATH, paste(igv.dir, '/Makefile', sep = "")))

  writeLines(sprintf('Instructins: Now go into shell, cd into directory "%s", and type either "make snapshots" or "make snapshots_batch" to generate screen shots either locally or on LSF, respectively.', igv.dir))

  return(dranger)
}

##############
#' @name dump_mutations_for_IGV_snapshots
#' @title dump_mutations_for_IGV_snapshots
#'
#' @description
#' dumps "snapshot_queue" file from maf into igv.dir for IGV screen shots
#' (should be followed by command line call to "make snapshots" or "make snapshots_batch" via Makefile in snapshots directory to take snapshots on the farm
#' of queued mutations on the farm)
#'
#' Optionally dumps "mutation_log" which is a text file input to Mutation Station (Eran Hodis 2011) .. this tdf provides
#' bsolute paths to pngs (resulting from the predicted outputs of the above) and columns for mutation review.
#'
#' Should set dump.mutation.log to TRUE for first run and then for subsequent re-runs (ie to re-do LSF failures) should run with
#' dump.mutation.log = FALSE, so as to not rewrite the mutation log.
#' @export
#' @author Marcin Imielinski
##############
dump_mutations_for_IGV_snapshots = function(maf, igv.dir,
  iset = "LUAD_all_exomes", wkspace = "Lung", # wkspace only needs to be provided if iset is a string
  mutation.log.name = 'mutation_log', snapshots.queue.filename = 'snapshot_queue', dump.mutation.log = TRUE,
  capture = TRUE) # capture = F  if you want to take images from whole genome bams)
{
  SNAPSHOT.SUBDIR = "images";
  MUTATION.MAKEFILE.PATH = '/home/unix/marcin/Scripts/Makefiles/Makefile.mutations';
  system(sprintf('mkdir -p %s/%s', igv.dir, SNAPSHOT.SUBDIR));

  if (is.character(iset))
    {
      print(sprintf('Pulling firehose individuals table %s from workspace %s ...', iset, wkspace))
      iset = get_fh(iset, wkspace = wkspace);
    }

  if (is.null(maf$patient_name))
    maf$patient_name = gsub('\\-Tumor', '', maf$Tumor_Sample_Barcode);

  # match up maf's with tumor and normal bam file paths
  upatient = data.frame(id = unique(maf$patient_name), tumorbam = '', normalbam = '', stringsAsFactors = F);
  maf$patient_ix = match(maf$patient_name, upatient$id);

  if (capture)
    upatient[, c('Tumor_bam', 'Normal_bam')] = iset[match(upatient$id, iset$Individual_Id), c('Tumor_clean_bam_file_capture', 'Normal_clean_bam_file_capture')]
  else
    upatient[, c('Tumor_bam', 'Normal_bam')] = iset[match(upatient$id, iset$Individual_Id), c('Tumor_clean_bam_file_wgs', 'Normal_clean_bam_file_wgs')]
  maf[, c('Tumor_bam', 'Normal_bam')] = upatient[maf$patient_ix, c('Tumor_bam', 'Normal_bam')]

  if (is.null(maf$Start_position))
    {
      maf$Start_position = as.numeric(gsub('g\\.chr\\w+\\:(\\d+).*', '\\1', maf$Genome_Change));
      maf$End_position = as.numeric(gsub('g\\.chr\\w+\\:(\\d+).*', '\\1', maf$Genome_Change));
    }

  if (is.null(maf$Chromosome))
    maf$Chromosome = as.numeric(gsub('Y', '24', gsub('X', '23', gsub('g\\.chr(\\w+)\\:\\d+.*', '\\1', maf$Genome_Change))));

  maf = maf[order(maf$patient_name), ]; # will optimize time for igv snapshots
  maf$rowid = 1:nrow(maf);
  maf$rowid = maf$Protein_Change;
  maf$label = maf_label(maf, short = TRUE, use.hugo.symbol=T);

  maf$snapshot_path = paste(file_path_as_absolute(paste(igv.dir, SNAPSHOT.SUBDIR, sep = "/")),
    paste(maf$label, '.png', sep = ""), sep = "/");

  maf$chr = paste('chr', maf$Chromosome, sep = "");

  out = maf[, c('chr', 'Start_position', 'End_position', 'snapshot_path', 'Tumor_bam', 'Normal_bam')]

  bad = which(is.na(out$Tumor_bam) | is.na(out$Normal_bam));

  if (length(bad)>0)
    {
      warning(sprintf('There are %d events without a tumor or normal bam link .. removing', length(bad)))
      out = out[-bad, ];
    }

  print(sprintf('Dumping %d total events for snapshotting', nrow(out)));

  write.table(out,
              paste(igv.dir, snapshots.queue.filename, sep = "/"),
              quote = F, sep = "\t", row.names = F, col.names = F)
  maf$reviewed = 'Unreviewed';

  mutation_log = maf[, c('rowid', 'Hugo_Symbol', 'patient_name', 'chr', 'Start_position', 'End_position', 'Variant_Classification', 'reviewed', 'label', 'snapshot_path')]
  names(mutation_log) = c('mutsigrank','gene','individual','chrom','start_pos','end_pos','type','ManualReviewJudgement','label', 'snapshot_path')

  if (dump.mutation.log == T)
      write.table(mutation_log,
                  paste(igv.dir, mutation.log.name, sep = "/"),
                  quote = F, sep = "\t", row.names = F, col.names = TRUE)

  mutation_log$label = maf$label;

  # transfer makefile to new igv snapshot directory to make images
  system(paste('cp', MUTATION.MAKEFILE.PATH, paste(igv.dir, '/Makefile', sep = "")))

  writeLines(sprintf('Instructions: Now go into shell, cd into directory "%s", and type either "make snapshots" or "make snapshots_batch" to generate screen shots either locally or on LSF, respectively.', igv.dir))
  
  return(mutation_log)
}


#' @name read_peaks
#' @title read_peaks
#'
#' @description
#'
#' reads UCSC / ENCODE / Roadmap broad and narrow peak data into granges object 
#'
#' @author Marcin Imielinski
#' @export
read_peaks = function(fn, chr.sub = TRUE)
    {
        ln = grep('^((track)|(browser)|\\#)', readLines(fn), value = TRUE, invert = TRUE)
        tab = fread(paste(ln, collapse = '\n'), header = FALSE)
        nm = c('chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak')
        setnames(tab, nm[1:ncol(tab)])
        if (chr.sub)
            tab[, chrom := gsub('chr', '', chrom)]
        return(seg2gr(tab))
    }
