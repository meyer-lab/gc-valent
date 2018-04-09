import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

output = pd.read_csv('IL2_IL15_chain-0.csv') # import as pandas array
headers = list(output)
log_posts = headers[0:17]
posts = headers[17:34]
IL2_residuals = headers[34:50]
IL15_residuals = headers[50:66]
traf_residuals = headers[66:94]


mat = output.as_matrix() # convert output to matrix
#print(mat)
print(posts)
for ii in range(17,34):
    plt.hist(mat[:,ii])
    plt.title(headers[ii])
    plt.show()
    
print(mat.shape)
sq_err = np.zeros((500, 50))
for ii in range(34,94):
    sq_err[:, ii-34] = (mat[:,ii])**2

IL2_sq_err = sq_err[:,0:15]
IL15_sq_err = sq_err[:,15:30]
traf_sq_err = sq_err[:,30:50]

IL2_sum_sq_err = np.sum(IL2_sq_err, axis=1)
IL15_sum_sq_err = np.sum(IL15_sq_err, axis=1)
traf_sum_sq_err = np.sum(traf_sq_err, axis=1)
print(headers)
print(traf_residuals)