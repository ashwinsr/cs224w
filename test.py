# Usage: python test.py
# This should print the FEVD matrices for three variables

from fevd import FEVD

DATA_FILE = 'test_data.csv'
DATA_COLUMNS = ['BND.Open', 'GLD.Open', 'IYF.Open']

fevd = FEVD(DATA_FILE, DATA_COLUMNS)
print fevd.run()
