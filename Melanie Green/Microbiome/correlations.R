library(propr)

alltaxa <- read.delim("H:\\Endocrinology\\Green\\Microbiome\\PCO vs control\\Data\\alltaxa_cts_29May2018.txt")
family <- read.delim("H:\\Endocrinology\\Green\\Microbiome\\PCO vs control\\Data\\family_cts_29May2018.txt")

# metadata looks like clinical data without any microbiome data
metadata <- read.delim("H:\\Endocrinology\\Green\\Microbiome\\PCO vs control\\Data\\metadata_16July2018.txt")
order <- read.delim("H:\\Endocrinology\\Green\\Microbiome\\PCO vs control\\Data\\order_cts_29May2018.txt")
phylum <- read.delim("H:\\Endocrinology\\Green\\Microbiome\\PCO vs control\\Data\\phylum_cts_29May2018.txt")
quant <- read.delim("H:\\Endocrinology\\Green\\Microbiome\\PCO vs control\\Data\\quantdata_23Aug2018b.txt")

clinical <- read.csv("H:\\Endocrinology\\Green\\Microbiome\\PCO vs control\\Data\\clinical data-9-2-2018.csv")
