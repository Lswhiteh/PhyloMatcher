import pandas as pd
import numpy as np

manual_df = pd.read_csv("../data/Logan/Mammalia_Equivalence_List.csv")
manual_df.drop(manual_df.columns[len(manual_df.columns) - 1], axis=1, inplace=True)
manual_df.reset_index(inplace=True)
taxo_df = pd.read_csv("../data/mammalia_species.csv", sep="\t")

man_list = list(manual_df.reset_index().to_numpy().flatten())
taxo_list = list(taxo_df.to_numpy().flatten())

man_list = list(set([i for i in man_list if str(i) != "nan"]))
taxo_list = list(set([i for i in taxo_list if str(i) != "nan"]))

unique = []
matches = []

for i in taxo_list:
    if i in man_list:
        matches.append(i)
    else:
        unique.append(i)

print("Unique:", len(unique))
print("Matches:", len(matches))
