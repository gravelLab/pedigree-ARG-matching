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
    x_axis_column_index = get_number_input_with_lower_bound("Specify the x-axis column index in the csv file:",
                                                            0)
    y_axis_column_index = get_number_input_with_lower_bound("Specify the y-axis column index in the csv file:",
                                                            0)
    y_axis_column_index_2 = get_number_input_with_lower_bound("Specify the y-axis column index in the csv file:",
                                                            0)
    figure_title = input("Specify the figure title: ")
    x_axis_title = input("Specify the x-axis title: ")
    y_axis_title = input("Specify the y-axis title: ")
    print(axis_options)
    x_axis_option = get_number_input_in_bounds("Choose the x-axis type:", 1, 2)
    y_axis_option = get_number_input_in_bounds("Choose the y-axis type:", 1, 2)
    x_vals, y_vals = [], []
    with open(input_csv_filepath, newline="") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            try:
                x_vals.append(float(row[x_axis_column_index]))
                y_vals.append(float(row[y_axis_column_index]) / float(row[y_axis_column_index_2]))
            except (ValueError, IndexError):
                continue  # skip bad rows
    plt.figure(figsize=(12, 9))
    plt.scatter(x_vals, y_vals)
    if x_axis_option == 2:
        plt.xscale("log")
    if y_axis_option == 2:
        plt.yscale("log")

    plt.xlabel(x_axis_title)
    plt.ylabel(y_axis_title)
    plt.title(figure_title)
    # Plot the y=0.95x line
    # x_vals = np.array(x_vals)
    # x_line = np.linspace(x_vals.min(), x_vals.max(), 1000)
    # y_line = x_line
    # plt.plot(x_line, y_line, color="blue", linestyle="-", label="y = x")
    plt.grid(True, which='major', linestyle='--', linewidth=1)
    plt.savefig(output_filepath, dpi=300)
    print(f"Scatter plot saved to {output_filepath}")


if __name__ == "__main__":
    main()
