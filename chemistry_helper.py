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

# database setting
password = "chemistry"
dbname = "Cluster0"
client = pymongo.MongoClient("mongodb+srv://dly1994:" + password + "@cluster0.nzxnf.mongodb.net/" + dbname + "?retryWrites=true&w=majority")
db = client.test


# read the user input
# ---------------------------------------------
def read_user_input(confirm, totalMass, ion):
    while not confirm or int(confirm) != 1:
        # reset totalMass
        totalMass = ""
        # get the totalMass from user
        while not totalMass or float(totalMass) < 0:
            totalMass = float(input("\nPlease enter total mass of your acids: "))
        print(totalMass)
        # reset ionsNum value
        ionsNum = "-1"
        # check for valid selection
        while not ionsNum or int(ionsNum) < 0 or int(ionsNum) > 6:
            ionsNum = input("\nPlease choose an ion from the following list: \n" 
                    + "0. Netrual \n"
                    + "1. Radical (-) \n" 
                    + "2. Hydrogen (H+) \n"
                    + "3. Ammonia (NH4+) \n"
                    + "4. Sodium (Na+) \n"
                    + "5. Lithium (Li+) \n"
                    + "6. Potassium (K+) \n\n")
        # get the ions value
        ionsNum = int(ionsNum)
        if ionsNum == 0:
            ion = "Netrual (0)"
        elif ionsNum == 1:
            ion = "Radical (-)"
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
        
        # reset confirm value
        confirm = "0"
        # ask user confirmation
        while not confirm or int(confirm) < 1 or int(confirm) > 2:
            confirm = input("\nDo you want to get all fatty acid combination for a Total Mass of: " + str(totalMass) 
                        + " with ion: " + ion + "? \n"
                        + "1. YES \n"
                        + "2. NO \n\n")

    return [confirm, totalMass, ion]


# Given a list & a totalMass, 
# get all triplet in the list that adds up to that totalMass
# ****todo need to handle the case if multiple acids have same mass****
# -------------------------------------------------
def get_all_combination(acidList, totalMass):
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

            tempTotalMass = float(firstAcidMass + secondAcidMass + thirdAcidMass)
            
            if round(tempTotalMass, 4) == round(totalMass, 4):
                result.append([firstAcid, secondAcid, thirdAcid])
                if startIndex == endIndex-1 and secondAcidMass == thirdAcidMass:
                    result.append([firstAcid, secondAcid, secondAcid])
                    result.append([firstAcid, thirdAcid, thirdAcid])
                    break
                else: # check if the acid mass at next index is the same as the current one
                    if secondAcidMass == float(acidList[startIndex+1]['Neutral Mass']):
                        startIndex += 1
                    else:
                        endIndex -= 1
            elif tempTotalMass < totalMass:
                startIndex += 1
            else:
                endIndex -= 1

    if len(result) == 0:
        print("No combination found")
    else:
        print("Found combination!")
        count = 1
        for item in result:
            table = BeautifulTable()
            table.column_headers = ["Compound", "Formula", "Neutral Mass", "Common Name", "IUPAC Name", "Description"]
            firstAcid = item[0]
            secondAcid = item[1]
            thirdAcid = item[2]
            table.append_row([firstAcid["Compound"], firstAcid["Formula"], firstAcid["Neutral Mass"], firstAcid["Common Name"], firstAcid["IUPAC Name"], firstAcid["Description"]])
            table.append_row([secondAcid["Compound"], secondAcid["Formula"], secondAcid["Neutral Mass"], secondAcid["Common Name"], secondAcid["IUPAC Name"], secondAcid["Description"]])
            table.append_row([thirdAcid["Compound"], thirdAcid["Formula"], thirdAcid["Neutral Mass"], thirdAcid["Common Name"], thirdAcid["IUPAC Name"], thirdAcid["Description"]])
            print(str(count) + '.')
            count += 1
            print(table)
    # return result #todo


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
