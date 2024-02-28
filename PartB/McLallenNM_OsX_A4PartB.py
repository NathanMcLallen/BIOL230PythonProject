# Makes a graph comparing coding region length with exon number for both augustus and genscan gene predictions of the D.melanogaster genome
# Also outputs the gene number, average exons, SD of exons, and highest exon number for both predictions in a csv file

# Dependencies: matplotlib, pandas, scikitlearn, statsmodels, numpy

import matplotlib.pyplot as plt
import pandas as pd
from sklearn import datasets
import statsmodels.api as sm
from statsmodels.formula.api import ols
import numpy

def main():

    # Extracts the needed columns from the augustus gene data CSV file and calls the function to calculate transcript lengths
    augData = pd.read_csv("./genes-augustus.csv", usecols=["txStart", "txEnd", "exonCount"])
    #print (augData.head())
    fullAugData = addTxLength(augData)
    #print (fullAugData.head())

    # Does the same for the genscan gene data
    GSData = pd.read_csv("./genes-genscan.csv", usecols=["txStart", "txEnd", "exonCount"])
    fullGSData = addTxLength(GSData)
    #print (fullGSData.head())
    

    # Puts required columns into their own dataframes to be graphed
    augLength = fullAugData.txLength
    augExons = fullAugData.exonCount

    GSLength = fullGSData.txLength
    GSExons = fullGSData.exonCount

    # Gets stats for exons and prints to csv file
    augStats = getBasicStats(augExons)
    GSStats = getBasicStats(GSExons)
    augStats.insert(0, "Augustus")
    GSStats.insert(0, "GenScan")

    statsDF = pd.DataFrame()
    statsDF = statsDF.assign (StatValue = ["Gene predicter", "Number of genes", "Average exons per gene", "Standard deviation of exons per gene", "Highest exon number"])
    statsDF = statsDF.assign (Augustus = augStats)
    statsDF = statsDF.assign (GenScan = GSStats)

    outFile = open("ExonData.csv", "w")
    #outFile.write(statsDF)
    statsDF.to_csv("ExonData.csv", sep="\t", index=False)
    outFile.close()

    # Makes a scatterplot with all important features
    plt.scatter(augLength, augExons)
    plt.scatter(GSLength, GSExons, alpha=0.4)

    plt.xlabel("Predicted gene transcript length (bp)")
    plt.ylabel("Predicted number of exons")
    plt.legend(["Augustus predictions", "Genscan predictions"])
    plt.title("Exons vs gene length for Augustus and GenScan predicitions of D.melanogaster genes")


    # Saves scatterplot to pdf
    plt.savefig("ExonsVsLength.pdf", format="pdf", bbox_inches="tight")
    #plt.show()
    


# adds a new column to a dataframe for transcript lengths by using the transcript start and end columns
def addTxLength(mainDataFrame):
    txEnds = mainDataFrame.txEnd
    txStarts = mainDataFrame.txStart
    codingLengths = []
    for i in range(len(mainDataFrame.index)):
        length = txEnds[i] - txStarts[i]
        codingLengths.append(length)
    toReturn = mainDataFrame.assign (txLength=codingLengths)
    return toReturn


# Given a single dataframe column, returns the number of datapoints, average, SD, and maximum in a list in that order
def getBasicStats(dataFrame):
    toReturn = []
    total = 0
    count = 0
    maximum = 0
    for i in range(len(dataFrame.index)):
        count += 1
        total = total + dataFrame[i]
        if (maximum < dataFrame[i]):
            maximum = dataFrame[i]
    average = total / count
    sdTop = 0
    for i in range(len(dataFrame.index)):
        toAdd = dataFrame[i] - average
        toAdd = toAdd * toAdd
        sdTop = sdTop + toAdd
    standardD = sdTop / count
    standardD = numpy.sqrt(standardD)
    toReturn.append(count)
    toReturn.append(average)
    toReturn.append(standardD)
    toReturn.append(maximum)
    #print (toReturn)
    return toReturn

    



if __name__ == "__main__":
    main()
    #End program
    print ("\nGraphs created")
    exit