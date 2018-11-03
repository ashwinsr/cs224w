# Usage: python test.py
# This should print the FEVD matrices for three variables

from fevd import FEVD

DATA_FILE = 'fevd_test_data.csv'
#DATA_COLUMNS = ['BND.Open', 'GLD.Open', 'IYF.Open']
DATA_COLUMNS = ['GLD.Open', 'BND.Open', 'IYF.Open']
#DATA_COLUMNS = ['IYF.Open', 'BND.Open', 'GLD.Open']

fevd = FEVD(DATA_FILE, DATA_COLUMNS)
print fevd.get_matrix(10)
