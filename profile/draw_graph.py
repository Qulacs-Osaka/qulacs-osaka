import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Change extension of `file_name` to png.
def default_output_file_name(file_name):
    splitted = file_name.split(".")
    splitted[-1] = "png"
    return ".".join(splitted)

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="CSV File to draw graph")
parser.add_argument("-o", "--output", help="Output file name")
args = parser.parse_args()

file_name = args.file
data = pd.read_csv(file_name)
if args.output is None:
    output_file_name = default_output_file_name(file_name)

x = np.array(data["qubit_count"])
y = np.array(data["memory_usage"])
plt.xlabel("Qubit Count")
plt.ylabel("Memory Usage(MB)")
plt.title("Memory Usage over Qubit Count")
plt.plot(x, y)
plt.savefig(output_file_name)