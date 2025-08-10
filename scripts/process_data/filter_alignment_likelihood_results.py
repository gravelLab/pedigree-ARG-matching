import pandas as pd

from scripts.utility.basic_utility import get_filepath

# Load CSV
input_filepath = get_filepath("Specify the filepath: ")
df = pd.read_csv(input_filepath)

# Filter where last column (number_of_alignments) is >= 100
filtered = df[df["number_of_alignments"] >= 100]

# Calculate mean and max of the second column (identity_percentile)
mean_val = filtered["identity_percentile"].mean()
median_val = filtered["identity_percentile"].median()
max_val = filtered["identity_percentile"].max()

print(f"Mean of second column (filtered): {mean_val}")
print(f"Max of second column (filtered): {max_val}")
print(f"Median: {median_val}")
