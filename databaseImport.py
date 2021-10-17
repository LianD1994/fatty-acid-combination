import csv
import json
import pandas as pd
import sys, getopt
import pymongo
from pymongo import MongoClient
# pprint library is used to make the output look more pretty
from pprint import pprint
from chemistry_helper import import_csv_to_mongodb

# read csv file
fatty_acid_csvfile = open('/Users/lian/code/fatty-acid-combination/csv/fatty_acid.csv', 'r', encoding='utf-8-sig')
# solution_ion_csvfile = open('/Users/lian/Documents/Chemistry/csv/solution_ion.csv', 'r')
ethyl_ester_csvfile = open('/Users/lian/code/fatty-acid-combination/csv/ethyl_ester.csv', 'r', encoding='utf-8-sig')

fatty_acid_header = ["Compound", "Formula", "Neutral Mass", "Common Name", "IUPAC Name", "Description"]
# solution_ion_header = ["Solution ", "Solution Name", "Formula", "Mass"]
ethyl_ester_header = ["Compound", "Formula", "Neutral Mass", "Common Name", "IUPAC Name", "Description"]

# import_csv_to_mongodb(fatty_acid_csvfile, fatty_acid_header, "acid")
import_csv_to_mongodb(ethyl_ester_csvfile, ethyl_ester_header, "ethyl_ester")
# import_csv_to_mongodb(solution_ion_csvfile, solution_ion_header, "ion")

# Issue the serverStatus command and print the results
# serverStatusResult=db.command("serverStatus")
# pprint(serverStatusResult)