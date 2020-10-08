# Author :Shanzida jahan Siddique 

## Set a working directory 
setwd ("/Users/shanzida/Documents/BIFX 552 Bioinformatics application /class 10:10")

## visualize working directory
getwd()

## install bioconductor package ,bioLite.R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
## Storing Generic Ranges with IRanges(page 270)
## installation of IRanges package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IRanges")
library(IRanges)
##Create IRanges objects.
-------------------------------------
## Each IRanges object has the two basic components of any range: a start and end position. We can create ranges with the IRanges() function. For example, a range starting at position 4 and ending at position 13 would be created with:
rng <- IRanges(start=4, end=13)
rng
# create ranges by specifying their width, and either start or end position:
IRanges(start=4, width=3)
IRanges(end=5, width=5)
# the IRanges() constructor (a function that creates a new object) can take vector arguments, creates an IRanges object containing many ranges:
x <- IRanges(start=c(4, 7, 2, 20), end=c(13, 7, 5, 23))
x
#Like many R objects, each range can be given a name. This can be accomplished by setting the names argument in IRanges, or using the function names():
names(x) <- letters[1:4]
x
# access the start positions, end positions, and widths of each range in this object with the methods start(), end(), and width():
start(x)
end(x)
width(x)
# These functions also work with <- to set start, end, and width position. For example, we could increment a range’s end position by 4 positions with:
end(x) <- end(x) + 4
x
# The range() method returns the span of the ranges kept in an IRanges object:
range(x)
#  subset IRanges just as we would any other R objects (vectors, dataframes, matrices), using either numeric, logical, or character (name) index
x[2:3]
start(x) < 5
x[start(x) < 5]
x[width(x) > 8]
x['a']
# Ranges can also be easily merged using the function c()
a <- IRanges(start=7, width=4)
b <- IRanges(start=2, end=5)
c(a, b)
## Basic Range Operations: Arithmetic, Transformations, and Set Operations
-------------------------------------------------------------------------------
  #IRanges objects can be grown or shrunk using arithmetic operations like +, -, and * (the division operator, /, doesn’t make sense on ranges, so it’s not supported). Growing ranges is useful for adding a buffer region. For example, we might want to include a few kilobases of sequence up and downstream of a coding region rather than just the coding region itself. With IRanges objects, addition (subtraction) will grow (shrink) a range symmetrically by the value added (subtracted) to it:
x <- IRanges(start=c(40, 80), end=c(67, 114))
x
x+4L
x-10L
# The IRanges package method restrict() cuts a set of ranges such that they fall inside of a certain bound
y <- IRanges(start=c(4, 6, 10, 12), width=13)
y


restrict(y, 5, 10)
# Another important transformation is flank(), which returns the regions that flank (are on the side of) each range in an IRanges object. flank() is useful in creating ranges upstream and downstream of protein coding genes that could contain pro‐ moter sequences.
x
flank(x, width=7)
# By default, flank() creates ranges width positions upstream of the ranges passed to it. Flanking ranges downstream can be created by setting start=FALSE:
flank(x, width=7, start=FALSE)
# Another common operation is reduce(). the reduce() operation takes a set of possi‐ bly overlapping ranges and reduces them to a set of nonoverlapping ranges that cover the same positions. Any overlapping ranges are merged into a single range in the result.
set.seed(0) # set the random number generator seed
alns <- IRanges(start=sample(seq_len(50), 20), width=5)
head(alns, 4)
reduce(alns)
# A similar operation to reduce() is gaps(), which returns the gaps (uncovered por‐ tions) between ranges.
gaps(alns)
# Another class of useful range operations are analogous to set operations. Each range can be thought of as a set of consecutive integers, so an IRange object like IRange(start=4, end=7) is simply the integers 4, 5, 6, and 7. This opens up the abil‐ ity to think about range operations as set operations like difference (setdiff()), intersection (intersect()), union (union()), and complement (which is simply the function gaps() we saw earlier)
a <- IRanges(start=4, end=13)
b <- IRanges(start=12, end=17) 
intersect(a, b)
setdiff(a, b)
union(b, a)
union(a, b)
##  Finding Overlapping Ranges
----------------------------------
# the basic task of finding overlaps between two sets of IRanges objects using the findOverlaps() function. findOverlaps() takes query and subject IRanges objects as its first two arguments
qry <- IRanges(start=c(1, 26, 19, 11, 21, 7), end=c(16, 30, 19, 15, 24, 8), names=letters[1:6])
sbj <- IRanges(start=c(1, 19, 10), end=c(5, 29, 16), names=letters[24:26]) 
qry
sbj
# Calling findOver laps(qry, sbj) returns an object with class Hits, which stores these overlaps:
hts <- findOverlaps(qry, sbj)
hts
# Access indices by using the accessor functions queryHits() and subjectHits().Find the names of each query and subject range with an overlap:
names(qry)[queryHits(hts)]
names(sbj)[subjectHits(hts)]
# limit the overlap results to only include query ranges that fall entirely within subject ranges with type=within 
hts_within <- findOverlaps(qry, sbj, type="within")
hts_within
# the options "first", "last", and "arbitrary" all lead findOverlaps() to return only one overlapping subject range per query (or NA if no overlap is found), results are returned in an integer vector where each element corresponds to a query range in qry:
findOverlaps(qry, sbj, select="first")
findOverlaps(qry, sbj, select="last")
findOverlaps(qry, sbj, select="arbitrary")
# Creating an IntervalTree object from an IRanges object:
sbj_it <- IntervalTree(sbj)
sbj_it
class(sbj_it)
#Using this sbj_it object illustrates we can use findOverlaps() with IntervalTree objects just as we would a regular IRanges object—the interfaces are identical:
findOverlaps(qry, sbj_it)
# Use hits objects to extract infor‐ mation from the overlapping ranges
as.matrix(hts)
countQueryHits(hts)
setNames(countQueryHits(hts), names(qry))
countSubjectHits(hts)
setNames(countSubjectHits(hts), names(sbj))
ranges(hts, qry, sbj)
# The functions subsetByOverlaps() and countOverlaps() simplify some of the most common operations performed on ranges once overlaps are found: keeping only the subset of queries that overlap subjects, and counting overlaps.
countOverlaps(qry, sbj)
subsetByOverlaps(qry, sbj)
## Finding Nearest Ranges and Calculating Distance
  ------------------------------------------------------
# Another common set of operations on ranges focuses on finding ranges that neighbor query ranges. In the IRanges package, there are three functions for this type of opera‐ tion: nearest(), precede(), and follow(). The nearest() function returns the near‐ est range, regardless of whether it’s upstream or downstream of the query. precede() and follow() return the nearest range that the query is upstream of or downstream of, respectively
qry <- IRanges(start=6, end=13, name='query')
sbj <- IRanges(start=c(2, 4, 18, 19), end=c(4, 5, 21, 24), names=1:4)
qry
sbj
nearest(qry, sbj)
precede(qry, sbj)
follow(qry, sbj)
# This family of functions for finding nearest ranges also includes distan ceToNearest() and distance(), which return the distance to the nearest range and the pairwise distances between ranges. We’ll create some random ranges to use in this example:
qry <- IRanges(sample(seq_len(1000), 5), width=10) 
sbj <- IRanges(sample(seq_len(1000), 5), width=10)
qry
sbj
# use distanceToNearest() to find neighboring ranges. It works a lot like findOverlaps()—for each query range, it finds the closest subject range, and returns everything in a Hits object with an additional column indicating the distance:
distanceToNearest(qry, sbj)
distance(qry, sbj)
  ## Run Length Encoding and Views
  -------------------------------------
# Run-length encoding and coverage()
#compress these runs using a scheme called run-length encoding. Run-length encoding compresses this sequence, storing it as: 3 fours, 2 threes, 1 two, 5 ones, 7 zeros, 3 ones, 7 fours.
x<-as.integer(c(4,4,4,3,3,2,1,1,1,1,1,0,0,0,
                  0, 0, 0, 0, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4))
xrle <- Rle(x)
xrle
#  revert back to vector form with as.vector():
as.vector(xrle)
# Run-length encoded objects support most of the basic operations that regular R vec‐ tors do, including subsetting, arithemetic and comparison operations, summary functions, and math functions:
xrle+4L
xrle/2
xrle>3
xrle[xrle > 3]
sum(xrle)
summary(xrle)
round(cos(xrle), 2)
# access an Rle object’s lengths and values using the functions run Lengths() and runValues():
runLength(xrle)
runValue(xrle)
#One place where we encounter run-length encoded values is in working with cover age(). The coverage() function takes a set of ranges and returns their coverage as an Rle object (to the end of the rightmost range). Simulating 70 random ranges over a sequence of 100 positions:
set.seed(0)
rngs <- IRanges(start=sample(seq_len(60), 10), width=7)
names(rngs)[9] <- "A" 
rngs_cov <- coverage(rngs)
rngs_cov
rngs_cov > 3 # where is coverage greater than 3?
rngs_cov[as.vector(rngs_cov) > 3] # extract the depths that are greater than 3
# subset Rle objects directly with IRanges objects:
rngs_cov[rngs['A']]
# mean coverage within this range
mean(rngs_cov[rngs['A']])
# Going from run-length encoded sequences to ranges with slice()
#take our coverage Rle object rngs_cov and slice it to create ranges corresponding to regions with more than 2x coverage:
min_cov2 <- slice(rngs_cov, lower=2)
min_cov2
# extract out the underlying ranges:
ranges(min_cov2)
#Advanced IRanges: Views
#summarize the views we created earlier using slice() using functions like viewMeans(), viewMaxs(), and even viewApply(), which applies an arbitrary function to views:
viewMeans(min_cov2)
viewMaxs(min_cov2)
viewApply(min_cov2, median)
# calculate the average coverage for windows 5-positions wide:
length(rngs_cov)
bwidth <- 5L
end <- bwidth * floor(length(rngs_cov) / bwidth)
windows <- IRanges(start=seq(1, end, bwidth), width=bwidth)
cov_by_wnd <- Views(rngs_cov, windows)
head(cov_by_wnd)
viewMeans(cov_by_wnd)
head(windows)
  ## Storing Genomic Ranges with GenomicRanges
  ------------------------------------------------
# create GRanges objects
library(GenomicRanges)
gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"),
              ranges=IRanges(start=5:8, width=10),
              strand=c("+", "-", "-", "+"))
gr
# specify the sequence lengths in the GRanges constructor, or set it after the object has been created using the seqlengths() function:
seqlens <- c(chr1=152, chr2=432, chr3=903)
gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"),
              ranges=IRanges(start=5:8, width=10),
              strand=c("+", "-", "-", "+"),
              gc=round(runif(4), 3),
              seqlengths=seqlens)
seqlengths(gr) <- seqlens # another way to do the same as above
gr
# access data in GRanges objects
start(gr)
end(gr)
width(gr)
# For the GRanges-specific data like sequence name and strand, there are new accessor functions—seqnames and strand:
seqnames(gr)
strand(gr)
# The returned objects are all run-length encoded. If we wish to extract all IRanges ranges from a GRanges object, we can use the ranges accessor function:
ranges(gr)
# GRanges has a length that can be accessed with length(), and supports names:
length(gr)
names(gr) <- letters[1:length(gr)]
gr
# All ranges with a start position greater than 7:
start(gr) > 7
gr[start(gr) > 7]
# Using the seqname() accessor, we can count how many ranges there are per chromosome and then subset to include only ranges for a particular chromosome:
table(seqnames(gr))
gr[seqnames(gr) == "chr1"]
# The mcols() accessor is used access metadata columns:
mcols(gr)
# The usual syntactic shortcut for accessing a column
mcols(gr)$gc
gr$gc
# compute the average GC content of all ranges on chr1:
mcols(gr[seqnames(gr) == "chr1"])$gc
mean(mcols(gr[seqnames(gr) == "chr1"])$gc)

  ## Grouping Data with GRangesList
  -------------------------------------
# GRanges Lists can be created manually:
gr1 <- GRanges(c("chr1", "chr2"), IRanges(start=c(32, 95), width=c(24, 123)))
gr2 <- GRanges(c("chr8", "chr2"), IRanges(start=c(27, 12), width=c(42, 34)))
grl <- GRangesList(gr1, gr2)
grl
unlist(grl)
doubled_grl <- c(grl, grl)
length(doubled_grl)
doubled_grl[2]
doubled_grl[[2]]
seqnames(grl)
start(grl)
# reate some random GRanges data, and demonstrate splitting by sequence name:
chrs <- c("chr3", "chr1", "chr2", "chr2", "chr3", "chr1")
gr <- GRanges(chrs, IRanges(sample(1:100, 6, replace=TRUE),
                            width=sample(3:30, 6, replace=TRUE)))

head(gr)
gr_split <- split(gr, seqnames(gr))
gr_split[[1]]
names(gr_split)
# use lapply() and sapply() to iterate through all elements and apply a function.
lapply(gr_split, function(x) order(width(x)))
sapply(gr_split, function(x) min(start(x)))
sapply(gr_split, length)
elementLengths(gr_split)
# reduce() called on a GRangesList object automatically works at the list-element level:
reduce(gr_split)

  ##Working with Annotation Data: GenomicFeatures and rtracklayer
  -------------------------------------------------------------------
  # install GenomicFeatures
library(BiocInstaller)
biocLite("GenomicFeatures")
biocLite("TxDb.Mmusculus.UCSC.mm10.ensGene")
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
txdb
class(txdb)
# access all gene regions in Mus musculus (in this version of Ensembl annotation). There’s a simple accessor function for this, unsurprisingly named genes():
mm_genes <- genes(txdb)
head(mm_genes)
head(mm_genes)
length(mm_genes)
# retrieve all exons grouped by transcript or gene:
mm_exons_by_tx <- exonsBy(txdb, by="tx") 
mm_exons_by_gn <- exonsBy(txdb, by="gene")
length(mm_exons_by_tx)
#  use a subset of chro‐ mosomes by setting which sequences our transcriptDb should query using the following approach:
seqlevels(txdb, force=TRUE) <- "chr1"
chr1_exons <- exonsBy(txdb, "tx")
all(unlist(seqnames(chr1_exons)) == "chr1")
# get all genes within this expanded region with:
qtl_region <- GRanges("chr8", IRanges(123260562, 123557264))
qtl_region_expanded <- qtl_region + 10e3
transcriptsByOverlaps(txdb, qtl_region_expanded)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")
library(rtracklayer)
mm_gtf <- import("/Users/shanzida/Documents/BIFX 552 Bioinformatics application /class 10:10/Mus_musculus.GRCm38.75_chr1.gtf.gz")
colnames(mcols(mm_gtf)) # metadata columns read in
getwd
set.seed(0)
pseudogene_i <- which(mm_gtf$gene_biotype == "pseudogene" &
                        mm_gtf$type == "gene")
pseudogene_sample <- sample(pseudogene_i, 5)
export(mm_gtf[pseudogene_sample], con="five_random_pseudogene.gtf",
       format="GTF")
bed_data <- mm_gtf[pseudogene_sample]
mcols(bed_data) <- NULL # clear out metadata columns
export(bed_data, con="five_random_pseudogene.bed", format="BED")
## Retrieving Promoter Regions: Flank and Promoters
------------------------------------------------------
table(mm_gtf$gene_biotype)
chr1_pcg <- mm_gtf[mm_gtf$type == "gene" &
                     mm_gtf$gene_biotype == "protein_coding"]
summary(width(chr1_pcg))
length(chr1_pcg)
chr1_pcg_3kb_up <- flank(chr1_pcg, width=3000)
chr1_pcg_3kb_up2 <- promoters(chr1_pcg, upstream=3000, downstream=0)
identical(chr1_pcg_3kb_up, chr1_pcg_3kb_up2)
# load the BSgenome.Mmusculus.UCSC.mm10 package and poke around
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome")
library(BSgenome)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome)
mm_gm <- BSgenome.Mmusculus.UCSC.mm10
organism(mm_gm)
providerVersion(mm_gm)
provider(mm_gm)
# BSge‐ nome packages contain sequences for each chromosome, stored in a list-like struc‐ ture we can access using indexing:
seqinfo(mm_gm)
mm_gm$chrM
mm_gm[[22]]
# BSgenome objects can be searched using the string-matching and alignment functions in the Bio strings packages.
library(Biostrings)
matchPattern("GGCGCGCC", mm_gm$chr1)
all(seqlevels(chr1_pcg_3kb_up) %in% seqlevels(mm_gm))

# With our chromosome names consistent between our GRanges promoter regions and the mouse BSgenome package, it’s easy to grab the sequences for particular regions kept in a GRanges object:
chr1_3kb_seqs <- getSeq(mm_gm, chr1_pcg_3kb_up)
chr1_3kb_seqs

  ## Getting Intergenic and Intronic Regions: Gaps, Reduce, and Setdiffs in Practice
------------------------------------------------------------------------------------------
  # With genomic ranges, gaps are calculated on every combination of strand and sequence. Here’s an example of this on a toy GRanges object so it’s clear what’s happening:
gr2 <- GRanges(c("chr1", "chr2"), IRanges(start=c(4, 12), width=6), strand=c("+", "-"), seqlengths=c(chr1=21, chr2=41))
gr2 # so we can see what these ranges look like
gaps(gr2)  
gr3<-gr2
strand(gr3) <- "*"
gaps(gr3)[strand(gaps(gr3)) == "*"]

 ## Finding and Working with Overlapping Ranges
-------------------------------------------------
# Using rtracklayer load bed file:
library(rtracklayer)
dbsnp137 <- import("mm10_snp137_chr1_trunc.bed.gz")
#subset to look at chromosome 1 exons (because our variants are only from chro‐ mosome 1):
collapsed_exons <- reduce(exons(txdb), ignore.strand=TRUE)
chr1_collapsed_exons <- collapsed_exons[seqnames(collapsed_exons) == "chr1"]
  # look at the length distribution of  variants:
summary(width(dbsnp137))
dbsnp137$name[which.max(width(dbsnp137))]
# count  zero-width features too, we’ll resize using the resize() function:
dbsnp137_resized <- dbsnp137
zw_i <- width(dbsnp137_resized) == 0
dbsnp137_resized[zw_i] <- resize(dbsnp137_resized[zw_i], width=1)
# tell findOverlaps() to ignore strand:
hits <- findOverlaps(dbsnp137_resized, chr1_collapsed_exons, ignore.strand=TRUE)
hits
length(unique(queryHits(hits)))
length(unique(queryHits(hits)))/length(dbsnp137_resized)
# use the method subsetByOverlaps():
subsetByOverlaps(dbsnp137_resized, chr1_collapsed_exons, ignore.strand=TRUE)
  
##Calculating Coverage of GRanges Objects
  --------------------------------------------
#using R’s sample() function:
set.seed(0)
chr19_len <- seqlengths(txdb)['chr19']
chr19_len
start_pos <- sample(1:(chr19_len-150), 2047719, replace=TRUE)
reads <- GRanges("chr19", IRanges(start=start_pos, width=150))
#use the coverage() method from GenomicRanges to calculate the coverage of these random reads:

cov_reads <- coverage(reads)
#calculate mean coverage per chromosome
mean(cov_reads)
#It’s also easy to calculate how much of this chromosome has no reads covering it (this will happen with shotgun sequencing, due to the random nature of read placement). We can do this two ways. First, we could use == and table():
table(cov_reads == 0)
# use some run-length encoding tricks (these are faster, and scale to larger data better):
sum(runLength(cov_reads)[runValue(cov_reads) == 0])
406487/chr19_len


  
