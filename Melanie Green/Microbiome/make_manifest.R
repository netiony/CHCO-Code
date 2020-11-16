# Read in files, get forward and reverse
files <- list.files("/Volumes/som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Laura Tim projects/Melanie Green/Microbiome/Data_Cleaned/fastq",full.names = T)
# Need unzipped fastq files. Remove compressed from manifest
files <- files[-c(grep("\\.bz2",files))]
forward <- files[grep("-R1-",files)]
reverse <- files[grep("-R2-",files)]
# Make manifest DF
manifest <- as.data.frame(matrix(nrow = length(forward),ncol = 3))
colnames(manifest) <- c("sample-id","forward-absolute-filepath","reverse-absolute-filepath")
manifest$`sample-id` <- sub("-R.*","",basename(forward))
manifest$`forward-absolute-filepath` <- forward
manifest$`reverse-absolute-filepath` <- reverse
# Write
write.table(manifest,file = "/Volumes/som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Laura Tim projects/Melanie Green/Microbiome/Data_Cleaned/manifest.txt",
            quote=FALSE, sep='\t',row.names = F)