################################################################################
# Forward Error Variance Decomposition
# ------------------------------------
# Given K time series variables, each of which have observations at T time
# periods, this computes the forward error variance decomposition between the
# variables.
#
# Usage:
#     from fevd import FEVD
#     fevd = FEVD('data.csv', ['Column_A', 'Column_B'])
#
# This will open the 'data.csv' file and extract the two columns specified and
# use them as the variables for FEVD. The class also accepts one optional
# parameter: lag_order.
#     lag_order: This is the number of values prior to the current timestep
#                that should be taken into account when computing the VAR. 
#                Default value is X because.
#
#     matrix = fevd.get_matrix(look_ahead)
#
#     look_ahead: This is the number of time periods ahead that the forecast
#                 willl estimate 


#FIXME. Default value is 15 because that's what Diebold and Yilmaz use in 2011. This sholud be pertrubed.
# This will return a 3-tensor of dimensions (look_ahead, K, K). For each of the
# look_ahead values, it will compute the K x K FEVD matrix. Note, that for each
# value in the K x K matrix, the value at position (i, j) represent the 

# at that look_ahead value.
################################################################################

import pandas as pd
import numpy as np
from statsmodels.tsa.api import VAR

class FEVD:
    def __init__(self, data_file_name, data_columns, lag_order=10):
        # Read raw csv and save columns we want
        data = pd.read_csv(data_file_name)
        data = data[data_columns]


        model = VAR(data.values)
        self.results = model.fit(lag_order)

    def get_matrix(self, look_ahead=10):
        fevd = self.results.fevd(look_ahead)
        
        return np.swapaxes(fevd.decomp, 0, 1)[look_ahead-1]
