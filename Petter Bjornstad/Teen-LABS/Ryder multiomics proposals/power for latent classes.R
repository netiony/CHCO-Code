t <- read_xlsx("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Grants/Ryder Teen-LABS multiomics/10yr Pct BMI change trajectory groupings 08_30_2023.xlsx")

table(t$GROUP)

# 1  2  3  4 
# 46 87 99 28 
# If we combine 1+2 and 3+4, this gives group sizes of 133 and 127