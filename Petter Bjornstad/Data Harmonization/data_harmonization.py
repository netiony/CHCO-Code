"""
This code is designed to pull data from multiple REDCap projects and harmonize
the data in a single dataset. Some studies are cross-sectional but include
measures at multiple visits, and some studies are longitudinal. So, this code
outputs data in a "semi-long" format with one row per study procedure, and a
visit column for longitudinal clustering.
"""
__author__ = "Tim Vigers"
__credits__ = ["Tim Vigers"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"

# Libraries
import os
os.chdir(
    "/Users/timvigers/Documents/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
from casper import clean_casper
from coffee import clean_coffee
from crocodile import clean_crocodile
from improve import clean_improve
from penguin import clean_penguin
from renal_heir import clean_renal_heir
# Use individual data functions to import cleaned DFs
casper = clean_casper()
coffee = clean_coffee()
crocodile = clean_crocodile()
improve = clean_improve()
penguin = clean_penguin()
renal_heir = clean_renal_heir()
