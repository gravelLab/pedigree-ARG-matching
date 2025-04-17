from scripts.utility import *


def main():
    # filename = get_file_path("Specify the path to the statistics file:\n")
    # filename = "pedigrees/likelihood_comparison_simulation/full"
    filename2 = "pedigrees/likelihood_comparison_simulation/likelihood"
    # result = input("Specify the name of the resulting image:\n")
    # result = f"{result}.png"
    result = "likelihood.png"

    # Assuming your data is stored in a file named 'data.txt'
    # with open(filename, 'r') as file:
    #     data = [line.split(": ") for line in file]
    #
    # x_values = [int(pair[0].strip()) for pair in data]
    # y_values = [int(pair[1].strip()) for pair in data]

    # Read the second dataset
    with open(filename2, 'r') as file2:
        data2 = [line.split(": ") for line in file2]

    x_values2 = [int(pair[0].strip()) for pair in data2]
    y_values2 = [int(pair[1].strip()) for pair in data2]

    plt.figure(figsize=(8, 8))

    # plt.scatter(x_values, y_values, color='blue')
    plt.scatter(x_values2, y_values2, color='red')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Proband number')
    plt.ylabel('Number of alignments')
    plt.title('Subclade Alignments')
    plt.savefig(result)
    plt.show()


if __name__ == "__main__":
    main()
