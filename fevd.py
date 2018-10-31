################################################################################
# Forward Error Variance Decomposition
# ------------------------------------
# Given K time series variables, each of which have observations at T time
# periods, this computes the forward error variance decomposition between the
# variables.
#
# Usage:
# TODO: explain lag order and default value
# TODO: explain look ahead and default value
# TODO: explain return (what is directionality of edges in matrix)
# TODO: explain usage
################################################################################

import pandas as pd
import numpy as np
from statsmodels.tsa.api import VAR

class FEVD:
    def __init__(self, data_file_name, data_columns, lag_order=10, look_ahead=10):
        # Read raw csv and save columns we want
        data = pd.read_csv(data_file_name)
        self.data = data[data_columns]

        # Save hyperparameters
        self.lag_order = lag_order
        self.look_ahead = look_ahead

    def run(self):
        model = VAR(self.data.values)
        results = model.fit(self.lag_order)
        fevd = results.fevd(self.look_ahead)
        
        return np.swapaxes(fevd.decomp, 0, 1)
