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
import redcap
import pandas as pd

# REDCap import variables
tokens = pd.read_csv(
    "~/Documents/Work/CHCO/Petter Bjornstad/Data Harmonization/api_tokens.csv")
uri = "https://redcap.ucdenver.edu/api/"

# Import CROCODILE
proj = redcap.Project(url=uri, your_rc_api_tkn)
