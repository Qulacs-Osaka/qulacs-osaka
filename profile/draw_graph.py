import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

def distinguish_method(file_name):
    if file_name.find("arnoldi") > -1:
        return "Arnoldi"
    elif file_name.find("lanczos") > -1:
        return "Lanczos"
    elif file_name.find("power") > -1:
        return "Power"
    else:
        return ""

# Change extension of `file_name` to png.
def default_output_file_name(file_name):
    splitted = file_name.split(".")
    splitted[-1] = "png"
    return ".".join(splitted)

def draw_graph(file_name):
    data = pd.read_csv(file_name)
    x = np.array(data["qubit_count"])
    y = np.array(data["memory_usage"])
    plt.xlabel("Qubit Count")
    plt.gca().get_xaxis().set_major_locator(ticker.MultipleLocator(1))
    plt.ylabel("Memory Usage(MB)")
    plt.yscale("log")

    method = distinguish_method(file_name)
    plt.title("Memory Usage over Qubit Count in {} method".format(method))
    plt.plot(x, y)
    output_file_name = default_output_file_name(file_name)
    plt.savefig(output_file_name)
    plt.clf()
    plt.close()

parser = argparse.ArgumentParser()
parser.add_argument("files", nargs="+", help="CSV Files to draw graph.")
args = parser.parse_args()
files = args.files

for file in files:
    draw_graph(file)
