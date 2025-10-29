import csv
from pathlib import Path

from matplotlib import pyplot as plt

from scripts.utility.basic_utility import get_filepath, get_number_input_with_lower_bound, get_number_input_in_bounds, \
    get_basename_without_extension

axis_options = ("Types of axis:\n"
                "1) Linear axis\n"
                "2) Log-scale axis\n")
log_axis_options = 2


def main():
    input_csv_filepath = get_filepath("Specify the input csv file path:")
    input_filename_no_ext = get_basename_without_extension(input_csv_filepath)
    output_filename = f"{input_filename_no_ext}.png"
    output_filepath = Path(input_csv_filepath).parent / output_filename

    x_axis_column_index = get_number_input_with_lower_bound(
        "Specify the x-axis column index in the csv file:", 0
    )
    y_axis_column_index = get_number_input_with_lower_bound(
        "Specify the y-axis column index in the csv file:", 0
    )

    #figure_title = input("Specify the figure title: ")
    #x_axis_title = input("Specify the x-axis title: ")
    #y_axis_title = input("Specify the y-axis title: ")

    print(axis_options)
    x_axis_option = get_number_input_in_bounds("Choose the x-axis type:", 1, 2)
    y_axis_option = get_number_input_in_bounds("Choose the y-axis type:", 1, 2)

    x_vals, y_vals = [], []
    with open(input_csv_filepath, newline="") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            try:
                x_vals.append(float(row[x_axis_column_index]))
                y_vals.append(float(row[y_axis_column_index]))
            except (ValueError, IndexError):
                continue  # skip bad rows

    # Make the figure larger — increase figsize as needed
    plt.figure(figsize=(14, 10))  # larger canvas

    # Create scatter plot
    plt.scatter(x_vals, y_vals, s=70)  # you can increase `s` for larger points

    # Log scale options
    if x_axis_option == 2:
        plt.xscale("log")
    if y_axis_option == 2:
        plt.yscale("log")

    # Axis labels and title with larger font sizes
    #plt.xlabel(x_axis_title, fontsize=20, labelpad=10)
    #plt.ylabel(y_axis_title, fontsize=20, labelpad=10)
    #plt.title(figure_title, fontsize=24, pad=20)

    # Make tick labels larger
    plt.tick_params(axis='both', which='major', labelsize=30, pad=20)
    plt.tick_params(axis='both', which='minor', labelsize=30, pad=20)

    # Grid for better readability
    plt.grid(True, which='major', linestyle='--', linewidth=1)

    # Save with higher DPI for clarity
    plt.tight_layout()
    plt.savefig(output_filepath, dpi=400)
    print(f"Scatter plot saved to {output_filepath}")


if __name__ == "__main__":
    main()
