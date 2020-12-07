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
fatty_acid_csvfile = open('/Users/lian/Documents/Chemistry/csv/fatty_acid.csv', 'r', encoding='utf-8-sig')
solution_ion_csvfile = open('/Users/lian/Documents/Chemistry/csv/solution_ion.csv', 'r')
fatty_acid_header = ["Compound", "Formula", "Neutral Mass", "Common Name", "IUPAC Name", "Description"]
solution_ion_header = ["Solution ", "Solution Name", "Formula", "Mass"]

import_csv_to_mongodb(fatty_acid_csvfile, fatty_acid_header, "acid")
import_csv_to_mongodb(solution_ion_csvfile, solution_ion_header, "ion")

# Issue the serverStatus command and print the results
# serverStatusResult=db.command("serverStatus")
# pprint(serverStatusResult)