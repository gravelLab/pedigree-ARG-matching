import pandas as pd
import matplotlib.pyplot as plt

from scripts.utility import get_filepath

# Read the CSV file
filepath = get_filepath("Specify the path to the data file:")
df = pd.read_csv(filepath)

# Define x-axis and values to plot
x = df['proband_number']
y_columns = ['no_solutions', 'individual_and_spouses', 'individual_and_non_spouse',
             'no_individual_spouse', 'neither_individual_nor_spouse']
cumulative_sum = df[y_columns].cumsum(axis=1)

# Plot the stacked area chart
plt.figure(figsize=(10, 6))

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']  # Example colors for each area

# Plot areas between polylines
for i, col in enumerate(y_columns):
    if i == 0:
        plt.fill_between(x, 0, df[col], color=colors[i], label=col, alpha=0.7)
    else:
        plt.fill_between(x, cumulative_sum.iloc[:, i-1], cumulative_sum.iloc[:, i], color=colors[i], label=col, alpha=0.7)

# Add labels and legend
plt.xlabel("Proband Number")
plt.ylabel("Values")
plt.title("Stacked Area Chart")
plt.legend(
    title="Categories",
    loc="upper center",
    bbox_to_anchor=(0.5, -0.1),  # Centered below the plot
    ncol=3,  # Number of columns in the legend
)

# Set all x-values as ticks
plt.xticks(ticks=x, labels=x)  # Rotate labels for better visibility if needed

# Adjust layout to prevent overlap
plt.tight_layout()
plt.subplots_adjust(bottom=0.2)  # Add extra space at the bottom

# Display the plot
plt.show()
