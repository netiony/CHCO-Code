library(readxl)
library(dplyr)
library(rtf)

f <- "/Users/pylell/Dropbox/TODAY DKD manuscript [shared]/CJASN response 3/TODAY DKD supplemental tables 1-12.xlsx"

# read in existing tables
t1 = read_excel(f,sheet = "Sup T1. Albuminuria")
t2 = read_excel(f, sheet = "Sup T2. Hyperfiltration")
t3 = read_excel(f, sheet = "Sup T3. Rapid eGFR decline")
t4 = read_excel(f, sheet = "Sup T4. Albuminuria H")
t5 = read_excel(f, sheet = "Sup T5. Hyperfiltration H")
t6 = read_excel(f, sheet = "Sup T6. Rapid eGFR decline H")
t7 = read_excel(f, sheet = "Sup T7. Albuminuria NHB")
t8 = read_excel(f, sheet = "Sup T8. Hyperfiltration NHB")
t9 = read_excel(f, sheet = "Sup T9. Rapid eGFR decline NHB")
t10 = read_excel(f, sheet = "Sup T10. Albuminuria NHW")
t11 = read_excel(f, sheet = "Sup T11. Hyperfiltration NHW")
t12 = read_excel(f, sheet = "Sup T12. Rapid eGFR decline NHW")

t1 <- t1 %>% select(AptName, Target, TargetFullName, UniProt, EntrezGeneID, EntrezGeneSymbol, 
                    estimate, std.error, p.value, adj.p.value, conf.low, conf.high)
t1$estimate <- round(t1$estimate, 2)
t1$std.error <- round(t1$std.error, 2)
t1$p.value <- ifelse(t1$p.value < 0.001, "<0.001", round(t1$p.value, 2))
t1$adj.p.value <- ifelse(t1$adj.p.value < 0.001, "<0.001", round(t1$adj.p.value, 2))
t1$conf.high <- round(t1$conf.high, 2)
t1$conf.low <- round(t1$conf.low, 2)

rtffile <- RTF("/Users/pylell/Dropbox/TODAY DKD manuscript [shared]/CJASN response 3/supplemental tables 1-12.doc")  # this can be an .rtf or a .doc
addTable(rtffile, t1)
done(rtffile)







