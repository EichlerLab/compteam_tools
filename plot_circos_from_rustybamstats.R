library("circlize")
library("stringr")
library("colorRamps")
library("tidyverse")
library("data.table")
library("argparser")

#################################
# plot_circos_from_rustybamstats.R
#
# written 12/18/2021 by kmiyamot (Katy Munson)
# last edited 1/11/2022 by kmiyamot
#
# Generates circos plots using the circlize package in R. Designed for visualizing whole-
# genome alignments against a reference or another high-quality assembly.
#
# Input: .bed file of alignments as generated from rustybam stats
# Output: .pdf files of 1) all query contigs vs. all ref contigs
#    2) Each ref contig against all querys that match to it
#    3) Each query contig against all refs that match to it
#
# Output is filtered based on lengths of contigs and alignments.
#
# Usage: Rscript plot_circos_from_rustybamstats.R -i <file name> -o <output file prefix>
#
#################################

#INPUT REQUIREMENTS
#Before starting, align query to ref with minimap2, e.g.
#
#module load minimap2/2.24
#minimap2 -x asm20 -c --eqx --secondary=no -t 48 {ref} {query} > alignments.paf
#
#Then prepare data for plot using rustybam
#
#module load rustybam/0.1.23
#rustybam stats --paf < alignments.paf > alignments_circos.bed

#################################
#
# Parameters to edit:
#
#################################

## Filtering Params
# Whole genome plot
# show links at least this long (bp)
GEN_MATCH_FILTER <- 200000
# show contigs at least this long (bp)
GEN_LENGTH_FILTER <- 10000000

# Individual ref/query contigs plots
# show links at least this long (bp)
TIG_MATCH_FILTER <- 50000
# show contigs at least this long (bp)
TIG_LENGTH_FILTER <- 1000000

## Output Formats
plot_whole_genome = TRUE
plot_chr_deck = TRUE
plot_tig_deck = TRUE

#################################
#
# Code follows:
#
#################################


## Input and output file names
p <- arg_parser("IO")
p <- add_argument(p, "--input", help="Input rustybam stats bed file", default="alignments_circos.bed")
p <- add_argument(p, "--output", help="Prefix for output files", default="circos_out")

my_file_names <- parse_args(p)

file_name = my_file_names$input
output_name = my_file_names$output


#Read rustybam stats file
mydata = read.table(
  file_name,
  header = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
)

#Calculate additional stats
mydata$ref_match_len <-
  mydata$reference_end - mydata$reference_start
mydata$query_match_len <- mydata$query_end - mydata$query_start

#Filter out information we won't possibly use
rough_match_filter <- min(GEN_MATCH_FILTER, TIG_MATCH_FILTER)
rough_length_filter <- min(GEN_LENGTH_FILTER, TIG_LENGTH_FILTER)

mydata <- mydata[which (
  mydata$ref_match_len > rough_match_filter &
    mydata$query_match_len > rough_match_filter &
    mydata$query_length > rough_length_filter &
    mydata$reference_length > rough_length_filter
) , ]


#select chromosome info for plot

ref_chrs = mydata %>% dplyr::select(c(X.reference_name, reference_length)) %>% unique()
colnames(ref_chrs) = c("chr", "end")
query_chrs = mydata %>% dplyr::select(c(query_name, query_length)) %>% unique()
colnames(query_chrs) = c("chr", "end")

#Sort ref chromosomes in size order.
ref_chrs <- ref_chrs %>% dplyr::arrange(desc(end))

#Calculate best ref match for query chrs
temp_data <- as.data.table(mydata)

pt <- temp_data[, .(match_len=sum(query_match_len)), by = c("query_name", "X.reference_name")]
shared_lengths <- pt %>% tidyr::pivot_wider(names_from = "X.reference_name", values_from = "match_len") %>% as.data.frame()
rownames(shared_lengths) = shared_lengths$query_name
shared_lengths <- shared_lengths[ , -which(names(shared_lengths) %in% c("query_name"))]

#order query chrs by best match ref size order
maxes <-
  cbind(rownames(shared_lengths), colnames(shared_lengths)[max.col(replace(shared_lengths, is.na(shared_lengths), -Inf), ties.method =
                                                                     "first")])
colnames(maxes) = c("query_chr", "ref_chr")
maxes <- merge(maxes, ref_chrs, by.x = "ref_chr", by.y = "chr")
colnames(maxes) <- gsub('^end', 'ref_end', names(maxes))
maxes <- merge(maxes, query_chrs, by.x = "query_chr", by.y = "chr")
colnames(maxes) <- gsub('^end', 'query_end', names(maxes))

query_chrs <-
  maxes %>% dplyr::arrange(ref_end) %>% dplyr::select(query_chr, query_end)
colnames(query_chrs) <- gsub('query_', '', names(query_chrs))
query_chrs$chr <- as.character(query_chrs$chr)



#Filter files
if (plot_whole_genome) {
  gendata <- mydata[which (
    mydata$ref_match_len > GEN_MATCH_FILTER &
      mydata$query_match_len > GEN_MATCH_FILTER &
      mydata$query_length > GEN_LENGTH_FILTER &
      mydata$reference_length > GEN_LENGTH_FILTER
  ) , ]
  
}
if (plot_tig_deck | plot_chr_deck) {
  tigdata <- mydata[which (
    mydata$ref_match_len > TIG_MATCH_FILTER &
      mydata$query_match_len > TIG_MATCH_FILTER &
      mydata$query_length > TIG_LENGTH_FILTER &
      mydata$reference_length > TIG_LENGTH_FILTER
  ) , ]
  
}

#clear original data table
rm(mydata)




##Plot asm vs genome circos
if (plot_whole_genome) {
  #Create chr pseudo-bed

  gen_ref_chrs <- ref_chrs[ref_chrs$chr %in% unique(gendata$X.reference_name),]
  gen_query_chrs <- query_chrs[query_chrs$chr %in% unique(gendata$query_name),]  
  all_chrs <- rbind(gen_ref_chrs, gen_query_chrs)
  all_chrs$start = 1
  all_chrs <- all_chrs[, c("chr", "start", "end")]

  
  #choose colors for links
  
  gen_ref_chrs$color = paste(primary.colors(length(gen_ref_chrs$chr), steps = , no.white = TRUE), "80", sep = "")
  gendata <-
    merge(gendata, gen_ref_chrs, by.x = "X.reference_name", by.y = "chr")
    
  #Create data pseudo-bed files
  
  bed1 <-
    gendata[, c("X.reference_name", "reference_start", "reference_end")]
  bed2 <- gendata[, c("query_name", "query_start", "query_end")]
  

  ### graphing code
  
  pdf(paste(output_name, "genome.pdf", sep="_"))
  

  # count ref and query contigs and set gaps between sets
  
  howmany_refchrs = length(gen_ref_chrs$chr)
  howmany_querychrs = length(gen_query_chrs$chr)

  gaps = c(rep(1,howmany_refchrs-1), 5, rep(1, howmany_querychrs-1), 5)

  circos.par(gap.after = gaps)

  # Increase size of canvas to avoid clipping of contig labels

  circos.par("canvas.xlim" = c(-1.5, 1.5), "canvas.ylim" = c(-1.5, 1.5))


  # initialize plot with no labels
  c = circos.genomicInitialize(all_chrs, plotType = NULL)

  # add contig labels, rotated so they don't overlap each other 

  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
        CELL_META$sector.index, cex = 0.6, facing = "clockwise", niceFacing = TRUE, adj = c(-0.1, 0))
  }, 
  track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)


# Add genomic axis (length ticks)

  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    circos.genomicAxis(h = "top", direction = "outside", labels.cex = 0.25) }, cell.padding = c(0, 0, 0, 0), bg.border = NA)


  # Plot chord diagram
  circos.genomicLink(bed1, bed2, col = gendata$color, border = "#00000000")
  circos.clear()

  # close plot
  dev.off()

  
}


if (plot_chr_deck | plot_tig_deck) {
  all_chrs <- rbind(ref_chrs, query_chrs)
  all_chrs$start = 1
  all_chrs <- all_chrs[, c("chr", "start", "end")]
}



##graph deck by chr

if (plot_chr_deck) {
  my_chrs <- ref_chrs$chr

  pdf(paste(output_name, "bychr.pdf", sep="_"))

  for (this_chr in my_chrs) {
    #subset links to match this chr
    match_data <-
      subset(tigdata [which (tigdata$X.reference_name == this_chr), ])
    
    #choose chrs that have valid links
    these_query_chrs <- unique(as.character(match_data$query_name))
    these_chrs <- append(this_chr, these_query_chrs)
    
    #generate new chr file
    match_chrs <- all_chrs[which (all_chrs$chr %in% these_chrs), ]
    
    #choose colors
    new_cols <-
      cbind(these_query_chrs, paste(primary.colors(
        length(these_query_chrs),
        steps = ,
        no.white = TRUE
      ), "80", sep = ""))
    colnames(new_cols) = c("query_name", "color_q")
    match_data <- merge(match_data, new_cols)
    match_data$color_q <- as.character(match_data$color_q)

    # count query contigs and set gaps between sets
    howmany_querychrs = length(these_query_chrs)
    gaps = c(5, rep(1, howmany_querychrs-1), 5)
    circos.par(gap.after = gaps)
    
    #plot this circos

    c = circos.genomicInitialize(match_chrs)
    
    bed1 <-
      match_data[, c("X.reference_name", "reference_start", "reference_end")]
    bed2 <-
      match_data[, c("query_name", "query_start", "query_end")]
    
    circos.genomicLink(bed1, bed2, col = match_data$color_q)
    circos.clear()
  }

  
  #close plot
  dev.off()

 
 
}



##graph deck by contig
if (plot_tig_deck) {
  my_chrs <- query_chrs$chr
  
  pdf(paste(output_name, "bytig.pdf", sep="_"))

  for (this_chr in my_chrs) {
    #subset links to match this chr
    match_data <-
      subset(tigdata [which (tigdata$query_name == this_chr), ])
    
    #choose chrs that have valid links
    these_ref_chrs <-
      unique(as.character(match_data$X.reference_name))
    these_chrs <- append(this_chr, these_ref_chrs)
    
    #generate new chr file
    match_chrs <- all_chrs[which (all_chrs$chr %in% these_chrs), ]
    
    #choose colors
    new_cols <-
      cbind(these_ref_chrs, paste(primary.colors(
        length(these_ref_chrs),
        steps = ,
        no.white = TRUE
      ), "80", sep = ""))
    colnames(new_cols) = c("X.reference_name", "color_r")
    match_data <- merge(match_data, new_cols)
    match_data$color_r <- as.character(match_data$color_r)
    
    # count ref contigs and set gaps between sets
    howmany_refchrs = length(these_ref_chrs)
    gaps = c(rep(1, howmany_refchrs-1), 5, 5)
    circos.par(gap.after = gaps)


    #plot this circos
    
    c = circos.genomicInitialize(match_chrs)
    
    bed1 <-
      match_data[, c("X.reference_name", "reference_start", "reference_end")]
    bed2 <-
      match_data[, c("query_name", "query_start", "query_end")]
    
    circos.genomicLink(bed1, bed2, col = match_data$color_r)
  }
  
  #close plot
  dev.off()
  circos.clear()
  
}
