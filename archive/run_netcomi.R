library(NetCoMi)
library(phyloseq)

# Load data sets
data("amgut1.filt") # ASV count matrix
data("amgut2.filt.phy") # phyloseq objext

# Agglomerate to genus level
amgut_genus <- tax_glom(amgut2.filt.phy, taxrank = "Rank6")

# Rename taxonomic table and make Rank6 (genus) unique
amgut_genus_renamed <- renameTaxa(amgut_genus,
                                  pat = "<name>",
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Rank6")

head(tax_table(amgut_genus))

# Using Spearman correlation instead of spring to avoid segfault
net_res <- netConstruct(amgut_genus_renamed,
                           taxRank = "Rank6",
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 50),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           measure = "sparcc",
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 2,
                           seed = 123456)

head(net_res$edgelist1)
