import matplotlib.pyplot as plt
from utility import *

# filename = get_file_path("Specify the path to the statistics file:\n")
filename = "pedigrees/simulation/results.txt"
# result = input("Specify the name of the resulting image:\n")
# result = f"{result}.png"
result = "results.png"

# Assuming your data is stored in a file named 'data.txt'
with open(filename, 'r') as file:
    data = [line.split(": ") for line in file]

x_values = [int(pair[0].strip()) for pair in data]
y_values = [int(pair[1].strip()) for pair in data]

plt.figure(figsize=(8, 8))

plt.scatter(x_values, y_values)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Proband number')
plt.ylabel('Number of alignments')
plt.title('Subclade Alignments')

# plt.xticks(x_values)
#y_ticks = list(range(0, 101, 25)) + list(range(100, max(y_values) + 1, 100))
#plt.yticks(y_ticks)
#x_ticks = list(range(0, 101, 25)) + list(range(100, max(x_values) + 1, 1000))
#plt.xticks(x_ticks)

plt.savefig(result)
plt.show()
