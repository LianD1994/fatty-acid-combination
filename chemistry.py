import sys
from chemistry_helper import read_user_input
from chemistry_helper import get_all_combination
from chemistry_helper import get_ion_mass
import pymongo
from pymongo import MongoClient
# pprint library is used to make the output look more pretty
from pprint import pprint

# database setting
password = "chemistry"
dbname = "Cluster0"
client = pymongo.MongoClient("mongodb+srv://dly1994:" + password + "@cluster0.nzxnf.mongodb.net/" + dbname + "?retryWrites=true&w=majority")
db = client.test

# Define constants
h2o = 18.010565
glycerol = 92.047344

while True:
    # read the user input
    userInput = read_user_input("0", "", "", "-1")
    totalMass = userInput[1]
    ion = userInput[2].split()
    tolerance = userInput[3]
    ionName = ion[0]
    ionMass = get_ion_mass(ion[1][1:-1])

    print("Calculating combinations for: \n" 
                + "Total Mass: " + str(totalMass) + "\n"
                + "ion: " + ionName + "(" + str(ionMass) + ")\n"
                + "tolerance: " + str(tolerance) + "\n")

    # Get list of fatty acids from database 
    cursor = db.acid.find({})
    acidList = []
    for document in cursor:
        acidList.append(document)

    # Get list of fatty acids from database 
    cursor = db.ethyl_ester.find({})
    ethyl_ester_list = []
    for document in cursor:
        ethyl_ester_list.append(document)

    # Calculate all possible combinations 
    get_all_combination(acidList, ethyl_ester_list, totalMass, tolerance, ionMass)
