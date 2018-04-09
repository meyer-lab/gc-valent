import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

output = pd.read_csv('IL2_IL15_chain-0.csv') # import as pandas array
headers = list(output)
log_posts = headers[0:17]
posts = headers[17:34]
IL2_residuals = headers[34:50]
IL15_residuals = headers[50:66]
traf_residuals = headers[66:82]


mat = output.as_matrix() # convert output to matrix
#print(mat)
print(posts)
for ii in range(17,34):
    plt.hist(mat[:,ii])
    plt.title(headers[ii])
    plt.show()