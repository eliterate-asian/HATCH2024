import argparse
import numpy as np
import os
import time
import matplotlib.pyplot as plt


def slideplot(seqA, seqB, window, thresh): #seqA should ALWAYS be ref genome
    if len(seqB) >= len(seqA):
        print("[ERROR] Fasta Sequence Equivalent or Larger than the Reference!")
        return 1
    
    slice_size = len(seqB)
    

    return dp

def dotplot(seqA, seqB, w, s):

    # determine size of dotplot matrix


    # Initialize the dotplot matrix with zeros
    dp = np.zeros((len(seqA), len(seqB)), dtype=int)

    # Iterate over all positions in seqA and seqB
    for i in range(len(seqA)):
        for j in range(len(seqB)):
            # Compute the window around position i in seqA and position j in seqB
            windowA = seqA[max(i-w, 0):min(i+w+1, len(seqA))]
            windowB = seqB[max(j-w, 0):min(j+w+1, len(seqB))]

            # Count the number of matching symbols in the window
            matches = sum([1 for x, y in zip(windowA, windowB) if x == y])

            # Set the matrix element to 1 if the number of matches is at least s
            if matches >= s:
                dp[i,j] = 1

    return dp

def dotplot2Ascii(dp, seqA, seqB, heading, filename):
    # Open the output file
    with open(filename, 'w') as f:
        # Write the heading
        f.write(heading + "\n")
        # Write the x-axis labels
        f.write(" " * 5 + seqB + "\n")
        # Write the dotplot matrix
        for i in range(len(seqA)):
            f.write(seqA[i] + " | ")
            for j in range(len(seqB)):
                if dp[i][j]:
                    f.write("*")
                else:
                    f.write(".")
            f.write("\n")

def dotplot2Graphics(dp, labelA, labelB, heading, filename):
    # create a new figure
    fig, ax = plt.subplots()

    # plot the dots using a scatter plot
    rows, cols = np.where(dp)
    ax.scatter(cols, rows, marker='.', color='black')

    # set the labels and title
    ax.set_xlabel(labelB)
    ax.set_ylabel(labelA)
    ax.set_title(heading)

    # set the tick positions and labels
    xticks = np.arange(0.5, dp.shape[1], 1)
    yticks = np.arange(0.5, dp.shape[0], 1)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(np.arange(1, dp.shape[1]+1))
    ax.set_yticklabels(np.arange(1, dp.shape[0]+1)[::-1])

    # save the figure to a file and display it on screen
    plt.savefig(filename)
    plt.show()

      

seq1 = './Tools/HomoSapien_Chr1.fasta'
seq2 = './Tools/testb.fasta' #Astro2_MTHFR.fasta' 
#'./Tools/HomoSapien_Chr1.fasta'

with open(seq1) as fileA, open(seq2) as fileB:
        seqRef = "".join([line.strip() for line in fileA if not line.startswith(">")])
        seqFas = "".join([line.strip() for line in fileB if not line.startswith(">")])
seqRef = seqRef.replace("N", "")

print(len(seqRef))
# starttime = time.time()
# print("Executing Matrix. Start.")
# dp = dotplot(seqRef, seqFas, 250, 205)
# plottitle = str("TestPlot-" + str(starttime)[-6:] +".png")
# dotplot2Graphics(dp, 'Sequence A', 'Sequence B', "SandboxTable", plottitle)
# executiontime = time.time() - starttime
# executiontime = executiontime/60
# print("Execution complete. End in Minutes: ", executiontime)
