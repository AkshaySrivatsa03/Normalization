# A collection of different normalization techniques for Data Cleaning and preparation

import csv
import math
import copy

originalInfo = {}
min_values = []
max_values = []
mean_values = []
std_vals = []
no_of_rows = 0

# values of means for each attribute when sample size is 10
sample_mean_values = [4.4, 127.3, 59.8, 19.9, 89.3, 27.709999999999997, 0.5078000000000001, 35.9]

# values of standard deviations for each attribute when sample size is 10
sample_std_values = [10.554619841567009, 120.68761328321975, 143.1401411205117, 152.96404806358908, 532.265065545354,
                     533.3357188488317, 533.3393289338412, 534.5275856189277]

# This method calculates the mean values for each attribute

def calcMeans(normalized_DataSet):

    global mean_values
    global no_of_rows

    mean_values = []

    done = False

    for i in range(len(normalized_DataSet["1"])):

        total = 0.0
        count = 0

        for key in normalized_DataSet:
            total += (normalized_DataSet[key])[i]
            count += 1

            if not done:
                no_of_rows = count
                done = True

        total = float(total)/count

        mean_values.append(total)
    return None

# This method calculates Standard Deviation for each attribute

def stdDev(normalized_DataSet):

    global mean_values
    global std_vals
    global no_of_rows

    std_vals = []

    calcMeans(normalized_DataSet)

    variance = 1.0
    std_val = 1.0

    for i in range(len(normalized_DataSet["1"])):

        for key in normalized_DataSet:
            variance += ((normalized_DataSet[key])[i] - mean_values[i] ) **2

        variance = float(variance) / no_of_rows
        std_val = math.sqrt(variance)

        std_vals.append(std_val)

    return None

# This method prints out the normalized dataset and labels it as per the method that called it

def printNewDataSet(newDataSet, Name):

    Name += ".csv"

    with open(Name, 'w', newline = '') as csvfile:
        filewriter = csv.writer(csvfile, delimiter = ',',
                            quotechar = '|', quoting = csv.QUOTE_MINIMAL)
        filewriter.writerow(['Patient', 'No.times pregnant','Plasma glucose concentration', 'Diastolic blood pressure'
                             , 'Tricepts skin fold thickness', '2-hr serum insulin','Body mass index',
                             'Diabetes pedigree function', 'Age'])

        for key in newDataSet:

            lis = newDataSet[key]
            listToWrite = []
            listToWrite.append(str(key))

            for i in range(len(lis)):
                listToWrite.append(lis[i])

            filewriter.writerow(listToWrite)

    return None

# This method finds the minimum of the attribute and returns that value

def minAttribute():

    global originalInfo
    global min_values

    for i in range(len(originalInfo['1'])):

        min_of_the_attribute = float("inf")

        for key in originalInfo:

            if min_of_the_attribute > (originalInfo[key])[i]:
                min_of_the_attribute = (originalInfo[key])[i]

        min_values.append(min_of_the_attribute)

    return None

# This method finds the maximum of the attribute and returns that value

def maxAttribute():

    global originalInfo
    global max_values

    for i in range(len(originalInfo["1"])):

        max_of_the_attribute = 0.0

        for key in originalInfo:

            if max_of_the_attribute < (originalInfo[key])[i]:
                max_of_the_attribute = (originalInfo[key])[i]

        max_values.append(max_of_the_attribute)

    return None

# This method gives us the standard normal distribution dataset

def stdNormalization():

    global originalInfo
    global mean_values
    global std_vals

    normalized_DataSet = copy.deepcopy(originalInfo)

    stdDev(normalized_DataSet)

    new_value = 1.0

    for i in range(len(normalized_DataSet["1"])):

        for key in normalized_DataSet:
            new_value = (normalized_DataSet[key])[i] - mean_values[i]
            new_value = float(new_value) / std_vals[i]
            (normalized_DataSet[key])[i] = new_value

    printNewDataSet(normalized_DataSet, "Standard Normalized Dataset")

    return None

# This method implements min max normalization

def minMaxNormalization():

    global originalInfo
    global max_values
    global min_values

    normalized_DataSet = copy.deepcopy(originalInfo)

    for i in range(len(normalized_DataSet["1"])):

        for key in normalized_DataSet:

            numerator = ( (normalized_DataSet[key])[i] - min_values[i] )
            denominator = ( max_values[i] - min_values[i] )
            if denominator > 0:
                new_value = float(numerator)/denominator
            else:
                new_value = float(numerator)

            (normalized_DataSet[key])[i] = new_value

    printNewDataSet(normalized_DataSet, "Min Max Normalized Dataset")
    return None

# This method implements the method to perform a variant of the min max normalization technique

def variantNormalization():

    global originalInfo
    global max_values
    global min_values

    new_max = 10.0
    new_min = 0.0

    normalized_DataSet = copy.deepcopy(originalInfo)

    for i in range(len(normalized_DataSet["1"])):

        for key in normalized_DataSet:

            numerator = ( (normalized_DataSet[key])[i] - min_values[i] )
            denominator = ( max_values[i] - min_values[i] )
            new_value = 1.0
            if denominator > 0:
                new_value = float(numerator)/denominator
            else:
                new_value = float(numerator)
            factor = float (new_max - new_min)
            new_value = (new_value * factor) + new_min

            (normalized_DataSet[key])[i] = new_value

    printNewDataSet(normalized_DataSet, "Variant Min Max Normalized Dataset")
    return None

# This method executes Decimal Scaling Normalization for the dataset

def decimalScaling():

    global originalInfo
    global max_values

    new_values = []

    normalized_DataSet = copy.deepcopy(originalInfo)

    for i in range(len(normalized_DataSet["1"])):

        new_value = 2.0
        k = 0

        while new_value > 1.0:
            new_value = max_values[i] / (10**k)
            k += 1
        new_values.append(k)

    for i in range(len(normalized_DataSet["1"])):
        for key in normalized_DataSet:
            (normalized_DataSet[key])[i] = float( (normalized_DataSet[key])[i] ) /(10** (new_values[i]) )

    printNewDataSet(normalized_DataSet, "Decimal Scaling Normalized Dataset")
    return None

# This method implements the sigmoidal normalization process

def sigmoidal():

    global mean_values
    global std_vals

    normalized_DataSet = copy.deepcopy(originalInfo)

    stdDev(normalized_DataSet)

    new_value = 1.0

    for i in range(len(normalized_DataSet["1"])):

        for key in normalized_DataSet:
            new_value = (normalized_DataSet[key])[i] - mean_values[i]
            new_value = float(new_value) / std_vals[i]
            numerator = (1 - (1.0 / math.exp(new_value)) )
            denominator = (1 + (1.0 / math.exp(new_value)) )
            new_value = float(numerator)/ denominator
            (normalized_DataSet[key])[i] = new_value

    printNewDataSet(normalized_DataSet, "Sigmoidal Normalized Dataset")
    return None

# This method calculates the Softmax normalization of the dataset

def Softmax():

    global mean_values
    global std_vals

    normalized_DataSet = copy.deepcopy(originalInfo)

    stdDev(normalized_DataSet)

    new_value = 1.0

    for i in range(len(normalized_DataSet["1"])):

        for key in normalized_DataSet:
            new_value = (normalized_DataSet[key])[i] - mean_values[i]
            new_value = float(new_value) / std_vals[i]
            denominator = (1 + (1.0 / math.exp(new_value)) )
            new_value = 1.0 / denominator
            (normalized_DataSet[key])[i] = new_value

    printNewDataSet(normalized_DataSet, "Softmax Normalized Dataset")

    return None

# This method implements the Student's t-distribution normaliztaion technique

def StudentDistrib():

    global originalInfo
    global sample_mean_values
    global sample_std_values
    global mean_values
    global no_of_rows

    normalized_DataSet = copy.deepcopy(originalInfo)

    calcMeans(normalized_DataSet)

    for i in range(len(normalized_DataSet["1"])):

        new_value = 1.0

        for key in normalized_DataSet:
            numerator = (sample_mean_values[i] - mean_values[i])
            denominator = ( float(std_vals[i]) / math.sqrt(no_of_rows) )
            new_value = float(numerator) / denominator
            (normalized_DataSet[key])[i] = new_value

    printNewDataSet(normalized_DataSet, "Student T Distribution Dataset")
    return None

def readData():

    global originalInfo

    data = open('Diabetes.csv', 'rt')
    # data = open('Test.csv', 'rt')
    read = csv.reader(data)

    for eachRow in read:

        if eachRow[0] != "Patient":
            listToAdd = []
            for i in range(1, len(eachRow)-1):
                listToAdd.append(float(eachRow[i]))
            originalInfo[eachRow[0]] = listToAdd

    return None

def Main():

    readData()

    minAttribute()
    maxAttribute()

    stdNormalization()
    minMaxNormalization()
    variantNormalization()
    decimalScaling()
    sigmoidal()
    Softmax()
    StudentDistrib()

    return None

Main()