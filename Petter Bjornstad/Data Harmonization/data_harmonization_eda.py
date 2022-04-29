import sweetviz as sv 
import pandas as pd
df = pd.read_csv('/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/merged_dataset_2022-04-29.csv')
my_report = sv.analyze(df,pairwise_analysis='off')
my_report.show_html(filepath='/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/eda.html',open_browser=False)