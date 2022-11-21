"""
Helper functions for data harmonization.
"""


def combine_checkboxes(df, base_name="", levels=[], sep=" & ",
                       drop=True):
    # Libraries
    import pandas as pd
    import numpy as np
    # A list of integer strings (one for each level)
    num = [str(s) for s in range(1, len(levels) + 1)]
    # Get column names
    cols = [base_name + "___" + s for s in num]
    # Which columns are checked
    values = df[cols].apply(lambda x: " & ".join(
        [str(s + 1) for s in np.where((x == 1) | (x == "1") | (x == "Checked"))[0]]), axis=1)
    # Dictionary for number to level translation
    column_translation = dict(zip(num, levels))
    # Replace numbers with levels
    for word, replacement in column_translation.items():
        values = [s.replace(word, replacement) for s in values]
    # Add to dataframe
    df[base_name] = values
    # Drop old columns
    if drop:
        df.drop(cols, axis=1, inplace=True)
    # Return df
    return df


def find_duplicate_columns(df, score_thresh=80):
    # Libraries
    import pandas as pd
    from fuzzywuzzy import fuzz
    from fuzzywuzzy import process
    from collections import Counter
    # Exact matches
    exact = [k for k, v in Counter(df.columns).items() if v > 1]
    # Remove underscores from column names for better similarity scoring.
    cols = [c.replace("_", " ") for c in df.columns]
    # Find each column name's closest match
    matches = []
    for var in cols:
        closest = process.extractBests(var, set(cols) - {var}, limit=1,
                                       score_cutoff=score_thresh)
        if len(closest) > 0:
            closest = (var,) + closest[0]
            matches.append(closest)
    # To dataframe
    matches = pd.DataFrame(matches,
                           columns=["original_variable", "closest_match", "score"])
    # Return
    return matches


def find_different_column_types(df1, df2):
    overlap = list(set(df1.columns).intersection(set(df2.columns)))
    diff = []
    for c in overlap:
        t1 = df1[c].dtypes
        t2 = df2[c].dtypes
        if t1 != t2:
            diff.append(c)
