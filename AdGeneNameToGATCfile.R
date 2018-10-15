# GATC in gene or not
#' @author Dominic Ritler
#' @import biomart, ggplot2
 
library(biomaRt)
library(ggplot2)


# ###########################  function part ###################################

#' Import GATC file and add additional column 
#' 
#' @param file.p the file path to the GATC file.
#' @return a data.frame 
#' @examples ReadGatcFile("myDamIDruns/results/123-gatc.txt")
ReadGatcFile <- function(file.p) {
  # read the GACT file origin from RDamID
  #
  # Args:
  #   file.p: (string) the path to the GATC file 
  #   
  # Returns:
  #   a data frame: [chromosome, start position of GATC, number of 
  #                 methylated by each sample] it also adds a additional column
  #                 to the data.frame filled with 0 for the genome information
  #                 to be added. 
  
  gatc.hits <- read.csv(file.p, row.names=1, sep="")
  gatc.hits$inGene <- rep(0, length(gatc.hits[,1]))
  
  return(gatc.hits)
}

#' download gene information from Biomart 
#' 
#' @param genome.name the ensembl biomart species name.
#' @return a data.frame containing chromosome, start stop and the gene name 
#' @examples getGeneInformation("celegans_gene_ensembl")
getGeneInformation <- function(genome.name) {
  
  database <- "ensembl"
  genome.name <- genome.name
  
  mart.querry <- useMart(database, dataset=genome.name)
  
  atr        <- c("ensembl_gene_id", "chromosome_name", "start_position", 
                  "end_position")
  filter     <- c("chromosome_name")
  
  
  querry.res <- data.frame(getBM(attributes=atr, #filters = filter, 
                                 mart=mart.querry))  #values = list(chr1),
  
  # ad chr to chromosome name as it is only I II III ... to match gatc file
  querry.res$chromosome_name <- paste("chr", querry.res$chromosome_name, sep="")
  
  
  return(querry.res)
}

#' add gene name to each gatc site if gene exist at this position 
#' 
#' @param gatc.hits the gatc data.frame
#' @param querry.res the biomart gene data.frame
#' @return the gact data.drame containing the genome information  
#' @examples AddGeneInfoToGatc(gatc.df, bimart.df)
AddGeneInfoToGatc <- function(gatc.hits, querry.res) {
  
  # make progress bar
  pb = txtProgressBar(min = 0, max = length(querry.res[,1]), initial = 0) 
  
  for (chro in unique(querry.res$chromosome_name)) {
    cat("\n", "Run chromosome: ", chro, "\n") # cat
    
    chr.subset <- which(querry.res$chromosome_name %in% chro)
    gatc.subset <- which(gatc.hits$chr %in% chro)
    
    for (i in chr.subset) {
      #chr.n <- querry.res$chromosome_name[i]
      g.sta <- querry.res$start_position[i]
      g.sto <- querry.res$end_position[i]
      
      
      #ch.match <- gatc.hits$chr %in% chr.n
      pos.match <- gatc.hits$start.position[gatc.subset] >= g.sta & gatc.hits$start.position[gatc.subset] <= g.sto
      #sum(pos.match)
      
      if ( sum(pos.match) != 0 ) {
        # get position and override 
        gatc.hits$inGene[gatc.subset[pos.match]] <- querry.res$ensembl_gene_id[i]
      } 
      setTxtProgressBar(pb, i)
    }
  }
  cat("\n")
  return(gatc.hits)
}

#' plot and print results
#' 
#' @param merged.df the gatc + genome information data.frame
#' @param result.name the file path to save the results
#' @param group.name the names of the group to compare
#' @param group.clust numeric grouping index
#' @return plots and saved files   
#' @examples MakeResult(gatc.df, bimart.df)
MakeResult <- function(merged.df, result.name, group.name, group.clust,
                       plot.res) {
  
  
  # save list
  write.csv2(merged.df, file = paste(result.name, "gene_list.csv"))
  
  # calculate presentage of total GACT in gene region  
  in.g <- sum(merged.df$inGene != 0)
  
  # calculate frequency
  freq.in.g <- in.g / length(merged.df$inGene) 
  
  pdf(paste(result.name, "GATC_in_gene_frequency.pdf"))
  barplot(c(freq.in.g, 1 - freq.in.g), 
          ylab = c("GATC frequency"),
          main = "Total GATC in gene",  
          names.arg = c(paste("in gene: ", round(freq.in.g, digits = 5)," %" ), 
                        paste("outside gene: ", round(1 - freq.in.g, digits = 5), " %")), 
          xlab = '')
  
  dev.off()
  
  
  
  # DAU test
  if (plot.res) {
    if (is.null(group.name)) {
      print("DamID samples from GATC file must be in group.name")
      return(NULL)
    }
    if (is.null(group.clust)) {
      print("groups for DamID samples from GATC files must be set")
      return(NULL)
    }
    
  # calculate numbers of group clusters 
  g.clusts <- unique(group.clust)
  
  # calculate methylated gatc in gene regions vs outside gene regions
  if (length(g.clusts) > 0) {
    for (gr in g.clusts) {
      
      # merge group
      gru.n <- group.name[which(group.clust == gr)]
      gatc.gr <- rowSums(merged.df[, names(merged.df) %in% gru.n])
      
      # find in gene
      temp <- data.frame(merged.df$inGene, gatc.gr)
      temp <- temp[which(temp[,2] != 0),] 
      
      in.g <- sum(temp[,1] != 0) / length(temp[,1]) # make frequency
      
      pdf(paste(result.name, "_", gr, "_group_GATC_in_gene_ratio.pdf"))
      barplot(c(in.g, 1 - in.g),  ylab = c("Met GATC frequency"),
              main = paste("Total methilated GATC in: ", toString(gru.n)), 
              names.arg = c(paste("in gene: ", round(in.g, digits = 5)," %" ), 
                            paste("outside gene: ", round(1 - in.g, digits = 5), " %")), 
              xlab = '')
      
      
      dev.off()
      
    }
  }
  
  
  met.in.g <- NULL
  met.in.intr <- NULL
  
}
}

#' add genome information to each GATC site in a GATC file originating from
#' a RDamID run.
#' 
#' @param file.p the file path to the GATC file.
#' @param genome.name ghe ensembl biomart species name.
#' @param group.name the names of the group to compare
#' @param group.clust numeric grouping index
#' @param result.name the file name to save the result in
#' @param plot.res if TRUE plot result 
#' @return the gatc file with the additional genome information column
#' @examples AddGenomeToGatc("myDamIDruns/results/123-gatc.txt")
AddGenomeToGatc <- function(file.p, genome.name, group.name=NULL, 
                            group.clust=NULL, result.name="res.txt", 
                            plot.res=FALSE) {
  
  # import gatc file
  print("Import GATC file:")
  gatc.df <- ReadGatcFile(file.p)
  
  # import genome information
  print("Get genome:")
  genome.df <- getGeneInformation(genome.name)
  
  # add genome information to gatc file
  print("make gatc gene list:")
  merged.df <- AddGeneInfoToGatc(gatc.df, genome.df)
  
  # plot data
  print("save files:")
  MakeResult(merged.df, result.name, group.name, group.clust, plot.res)
  
}


############################## run part ########################################


# add GATC file path
file.p <- "reads-t-gatc-1536588248.txt"
# add species tag
genome.name <- "celegans_gene_ensembl"
# do you want to have gatc ration in vs outside gene of your samples
plot.res <- TRUE
# add names of samples you want to plot the gatc in gene ratio
group.name <- c("gfpDamS1", "gfpDamS2", "TS1", "TS2")
# add groups (all group.name with the same number are in one group)
group.clust <- c(1, 1, 2, 2)
# add result file path and name
result.name <- "gatc_gene_test"



# run script
AddGenomeToGatc(file.p = file.p, genome.name = genome.name, 
                group.name = group.name, group.clust = group.clust,
                result.name = result.name, plot.res = plot.res)

