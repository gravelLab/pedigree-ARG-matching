import pandas as pd
import matplotlib.pyplot as plt
from scripts.utility import get_filepath, get_non_existing_path


def main():
    filepath = get_filepath("Specify the path to the data file:")
    result_name = get_non_existing_path("Specify the resulting figure's filepath (without extension):")
    title = input("Specify the title of the figure:")

    df = pd.read_csv(filepath)

    # Define x-axis and values to plot
    x = df['proband_number']
    y_columns = ['no_solutions', 'individual_and_spouses', 'individual_and_non_spouse',
                 'no_individual_spouse', 'neither_individual_nor_spouse']

    # Normalize the values to percentages
    df['total'] = df[y_columns].sum(axis=1)
    for col in y_columns:
        df[col] = (df[col] / df['total']) * 100

    # Calculate cumulative sum for stacked plot
    cumulative_sum = df[y_columns].cumsum(axis=1)

    plt.figure(figsize=(27.16, 14))

    colors = ['#898989', '#00452c', '#00965f', '#f1c338', '#be5103']
    handles = []

    for i, col in enumerate(y_columns):
        if i == 0:
            fill = plt.fill_between(x, 0, df[col], color=colors[i], label=col, alpha=0.7)
        else:
            fill = plt.fill_between(x, cumulative_sum.iloc[:, i - 1], cumulative_sum.iloc[:, i], color=colors[i],
                                    label=col, alpha=0.7)
        handles.append(fill)

    plt.xlabel("Proband Number", fontsize=20)
    plt.ylabel("Percentage (%)", fontsize=20)
    plt.title(title, fontsize=20)

    # Set y-axis limits explicitly
    plt.ylim(0, 100)

    # Define custom legend
    legend_labels = [
        "No solutions",
        "Individual present, the rest are the spouses",
        "Individual present, and a non-spouse assignment present",
        "Individual not present, spouse present",
        "Neither individual nor spouse present",
    ]
    grouped_handles = handles

    plt.legend(
        grouped_handles,
        legend_labels,
        title="Categories",
        title_fontsize=28,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=3,
        columnspacing=1.5,
        frameon=False,
        fontsize=28
    )

    plt.xticks(ticks=x, labels=x, fontsize=18)
    # plt.gca().axes.get_xaxis().set_visible(False)

    plt.yticks(fontsize=18)
    # plt.xlim(min(x), max(x))
    plt.xlim(4, max(x))
    # Use tight_layout for better adjustment
    plt.tight_layout()

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)

    # Save and display the plot
    plt.savefig(fname=f"{result_name}.svg")
    plt.show()


if __name__ == "__main__":
    main()
