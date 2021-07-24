import sys
import csv
import json
from tabulate import tabulate
from beautifultable import BeautifulTable
import pandas as pd
import sys, getopt
import pymongo
from pymongo import MongoClient
# pprint library is used to make the output look more pretty
from pprint import pprint

# Define constants
h2o = 18.010565
glycerol = 92.047344

# database setting
password = "chemistry"
dbname = "Cluster0"
client = pymongo.MongoClient("mongodb+srv://dly1994:" + password + "@cluster0.nzxnf.mongodb.net/" + dbname + "?retryWrites=true&w=majority")
db = client.test


# read the user input
# ---------------------------------------------
def read_user_input(confirm, glycerolipids, totalMass, ion, tolerance):
    while not confirm or int(confirm) != 1:
        # reset Glycerolipids
        glycerolipidsNum = "-1"
        while not glycerolipidsNum or int(glycerolipidsNum) < 1 or int(glycerolipidsNum) > 3:
            glycerolipidsNum = input("\nPlease choose a Glycerolipids from the following list: \n"
                    + "1. Monoacylglycerol (MAG) \n"
                    + "2. Diacylglycerol (DAG) \n"
                    + "3. Triacylglycerol (TAG) \n\n")
        # get the glycerolipid choice
        glycerolipidsNum = int(glycerolipidsNum)
        if glycerolipidsNum == 1:
            glycerolipids = "MAG"
        elif glycerolipidsNum == 2:
            glycerolipids = "DAG"
        elif glycerolipidsNum == 3:
            glycerolipids = "TAG"
    
        # reset totalMass
        totalMass = ""
        # get the totalMass from user
        while not totalMass or float(totalMass) < 0:
            totalMass = float(input("\nPlease enter total mass of your acids: "))

        # reset ionsNum value
        ionsNum = "-1"
        # check for valid selection
        while not ionsNum or int(ionsNum) < 0 or int(ionsNum) > 6:
            ionsNum = input("\nPlease choose an ion from the following list: \n" 
                    + "0. Radical (-) \n"
                    + "1. Minus Hydrogen (-H+) \n" 
                    + "2. Hydrogen (H+) \n"
                    + "3. Ammonia (NH4+) \n"
                    + "4. Sodium (Na+) \n"
                    + "5. Lithium (Li+) \n"
                    + "6. Potassium (K+) \n\n")
        # get the ions value
        ionsNum = int(ionsNum)
        if ionsNum == 0:
            ion = "Radical (-)"
        elif ionsNum == 1:
            ion = "Minus_Hydrogen (-H+)"
        elif ionsNum == 2:
            ion = "Hydrogen (H+)"
        elif ionsNum == 3:
            ion = "Ammonia (NH4+)"
        elif ionsNum == 4:
            ion = "Sodium (Na+)"
        elif ionsNum == 5:
            ion = "Lithium (Li+)"
        elif ionsNum == 6:
            ion = "Potassium (K+)"

        # reset tolerance 
        tolerance_choice = "-1"
        # 0.5  0.2  0.1  0.05  0.01  0.005  0.001  0.0005
        while not tolerance_choice or int(tolerance_choice) < 1 or int(tolerance_choice) > 8:
            tolerance_choice = input("\nPlease choose the mass tolerance for your result from the following list: \n"
                        + "1. +/- 0.5\n"
                        + "2. +/- 0.2\n"
                        + "3. +/- 0.1\n"
                        + "4. +/- 0.05\n"
                        + "5. +/- 0.01\n"
                        + "6. +/- 0.005\n"
                        + "7. +/- 0.001\n"
                        + "8. +/- 0.0005\n\n")
        tolerance_choice = int(tolerance_choice)
        if tolerance_choice == 1:
            tolerance = 0.5
        if tolerance_choice == 2:
            tolerance = 0.2
        if tolerance_choice == 3:
            tolerance = 0.1
        if tolerance_choice == 4:
            tolerance = 0.05
        if tolerance_choice == 5:
            tolerance = 0.01
        if tolerance_choice == 6:
            tolerance = 0.005
        if tolerance_choice == 7:
            tolerance = 0.001
        if tolerance_choice == 8:
            tolerance = 0.0005

        # reset confirm value
        confirm = "0"
        # ask user confirmation
        while not confirm or int(confirm) < 1 or int(confirm) > 2:
            confirm = input("\nDo you want to get all fatty acid combinations for " + glycerolipids + " for a Total Mass of: " + str(totalMass) 
                        + " with ion: " + ion 
                        + " and a mass tolerance of +/- " + str(tolerance) + "? \n"
                        + "1. YES \n"
                        + "2. NO \n\n")

    return [confirm, glycerolipids, totalMass, ion, tolerance]


# Given a list & a targetMass, 
# get all triplet in the list that adds up to that targetMass
# ****todo need to handle the case if multiple acids have same mass****
# -------------------------------------------------
def get_all_combination(acidList, targetMass, tolerance, ionMass, numOfAcid):
    if numOfAcid == 3:
        result = find_three_acid_combination(acidList, targetMass, tolerance, ionMass)
    elif numOfAcid == 2:
        result = find_two_acid_combination(acidList, targetMass, tolerance, ionMass)
    elif numOfAcid == 1:
        result = find_one_acid(acidList, targetMass, tolerance, ionMass)
    
    if len(result) == 0:
        print("No combination found")
    else:
        print("Found combination!")
        count = 1
        for item in result:
            table = BeautifulTable(maxwidth=140, numeric_precision=4)
            table.column_headers = ["Compound", "Formula", "Neutral Mass", "Common Name", "IUPAC Name", "Description"]
            
            firstAcid, secondAcid, thirdAcid = {}, {}, {}
            stardardMass = -1
            if numOfAcid > 0:
                firstAcid = item[0]
                table.append_row([firstAcid["Compound"], firstAcid["Formula"], firstAcid["Neutral Mass"], firstAcid["Common Name"], firstAcid["IUPAC Name"], firstAcid["Description"]])
            if numOfAcid > 1:
                secondAcid = item[1]
                table.append_row([secondAcid["Compound"], secondAcid["Formula"], secondAcid["Neutral Mass"], secondAcid["Common Name"], secondAcid["IUPAC Name"], secondAcid["Description"]])
            if numOfAcid > 2:
                thirdAcid = item[2]
                table.append_row([thirdAcid["Compound"], thirdAcid["Formula"], thirdAcid["Neutral Mass"], thirdAcid["Common Name"], thirdAcid["IUPAC Name"], thirdAcid["Description"]])

            standardMass = item[numOfAcid]
            percentageError = item[numOfAcid+1]

            print(str(count) + '.')
            # debug
            print('Standard Mass - ' + str(round(standardMass, 4)))
            print('Percentge Error - ' + str(round(percentageError, 4))) 
            count += 1
            print(table)
    # return result #todo


def find_three_acid_combination(acidList, targetMass, tolerance, ionMass):
    # testing value: TAG, 835.7749, Hydrogen
    result = []
    # sort the dictionary by value
    acidList.sort(key = lambda x: float(x['Neutral Mass']))

    N = len(acidList)
    for i in range(N):
        firstAcid = acidList[i]
        firstAcidMass = float(firstAcid['Neutral Mass'])

        startIndex, endIndex = i, N-1
        while startIndex <= endIndex:
            secondAcid = acidList[startIndex]
            secondAcidMass = float(secondAcid['Neutral Mass'])

            thirdAcid = acidList[endIndex]
            thirdAcidMass = float(thirdAcid['Neutral Mass'])

            tempTargetMass = float(firstAcidMass + secondAcidMass + thirdAcidMass)
            
            if abs(tempTargetMass - targetMass) <= tolerance:
                standardMass = tempTargetMass + float(ionMass) + glycerol - h2o*3
                totalMass = targetMass + float(ionMass) + glycerol - h2o*3
                percentageError = (standardMass - totalMass) / standardMass * 1000000
                result.append([firstAcid, secondAcid, thirdAcid, standardMass, percentageError])
                if startIndex == endIndex-1 and secondAcidMass == thirdAcidMass:
                    result.append([firstAcid, secondAcid, secondAcid, standardMass, percentageError])
                    result.append([firstAcid, thirdAcid, thirdAcid, standardMass, percentageError])
                    break
                else: # check if the acid mass at next index is the same as the current one
                    if startIndex < endIndex and secondAcidMass == float(acidList[startIndex+1]['Neutral Mass']):
                        startIndex += 1
                    else:
                        endIndex -= 1
            elif tempTargetMass < targetMass:
                startIndex += 1
            else:
                endIndex -= 1

    return result


def find_two_acid_combination(acidList, targetMass, tolerance, ionMass):
    # testing value: DAG, 597.5452, Hydrogen
    result = []
    # sort the dictionary by value
    acidList.sort(key = lambda x: float(x['Neutral Mass']))

    startIndex = 0
    endIndex = len(acidList)-1
    while startIndex <= endIndex:
        firstAcid = acidList[startIndex]
        firstAcidMass = float(firstAcid['Neutral Mass'])

        secondAcid = acidList[endIndex]
        secondAcidMass = float(secondAcid['Neutral Mass'])

        tempTargetMass = float(firstAcidMass + secondAcidMass)
            
        if abs(tempTargetMass - targetMass) <= tolerance:
            standardMass = tempTargetMass + float(ionMass) + glycerol - h2o*2
            totalMass = targetMass + float(ionMass) + glycerol - h2o*2
            percentageError = (standardMass - totalMass) / standardMass * 1000000
            result.append([firstAcid, secondAcid, standardMass, percentageError])
            if startIndex == endIndex-1 and firstAcidMass == secondAcidMass:
                result.append([firstAcid, firstAcid, standardMass, percentageError])
                result.append([secondAcid, secondAcid, standardMass, percentageError])
                break
            else:
                if startIndex < endIndex and firstAcidMass == float(acidList[startIndex+1]['Neutral Mass']):
                    startIndex += 1
                else:
                    endIndex -= 1
        elif tempTargetMass < targetMass:
            startIndex += 1
        else:
            endIndex -= 1

    return result


def find_one_acid(acidList, targetMass, tolerance, ionMass):
    # testing value: MAG, 331.2843, Hydrogen
    result = []
    # sort the dictionary by value
    acidList.sort(key = lambda x: float(x['Neutral Mass']))

    N = len(acidList)
    for i in range(N):
        firstAcid = acidList[i]
        firstAcidMass = float(firstAcid['Neutral Mass'])
            
        if abs(firstAcidMass - targetMass) <= tolerance:
            standardMass = firstAcidMass + float(ionMass) + glycerol - h2o
            totalMass = targetMass + float(ionMass) + glycerol - h2o
            percentageError = (standardMass - totalMass) / standardMass * 1000000
            result.append([firstAcid, standardMass, percentageError])

    return result


def import_csv_to_mongodb(file, header, type):
    # reset database
    if type == "acid":
        db.acid.drop()
    elif type == "ion":
        db.ion.drop()

    reader = csv.DictReader(file)
    for each in reader:
        print(each)
        row={}
        for field in header:
            row[field]=each[field]
        if type == "acid":
            db.acid.insert(row)
        elif type == "ion":
            db.ion.insert(row)


def get_ion_mass(ionFormula):
    if ionFormula == "0": 
        return 0
    document = db.ion.find_one({'Formula': ionFormula})
    return document['Mass']
