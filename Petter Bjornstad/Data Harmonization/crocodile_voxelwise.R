library(R.matlab)
# Import
controls = readMat("/Users/timvigers/Dropbox/Work/Petter Bjornstad/CROCODILE/Data_Raw/IndivValuesVoxelWise/Controls.mat")
t1d = readMat("/Users/timvigers/Dropbox/Work/Petter Bjornstad/CROCODILE/Data_Raw/IndivValuesVoxelWise/T1D.mat")
# Format
controls_indiv = unlist(controls[[3]],recursive = F)
names(controls_indiv) = rep(rep(1:10,each = 4),2)
