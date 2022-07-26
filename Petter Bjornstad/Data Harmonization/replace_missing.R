df = read.csv("/Users/timvigers/Downloads/PENGUIN_DATA_2022-07-26_1225.csv",na.strings = "")

# Check for missing
df[df == -97] = -9997
df[df == -98] = -9998
df[df == -99] = -9999

df[df == -997] = -9997
df[df == -998] = -9998
df[df == -999] = -9999

write.csv(df,"~/penguin.csv",na="",row.names = F)
