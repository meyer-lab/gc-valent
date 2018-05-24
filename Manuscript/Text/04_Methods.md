# Methods

All analysis was implemented in Python, and can be found at <https://github.com/meyer-lab/type-I-ckine-model>, release 1.0 (doi: [00.0000/arc0000000](https://doi.org/doi-url)).



## Model

### Base model



### Parameters and assumptions



### Model fitting

### Tensor Generation and Factorization
### Tensor generation
After fitting the model and obtaining the posterior distributions of the unknown parameters, we generated a multi-dimensional dataset for the receptor-complex concentrations at different timepoints for various combinations of receptor expression rates and initial ligand concentrations. The receptor expression rates are those of IL2Ra, IL2Rb, gc, IL15Ra, IL7R, and IL9R, and they range from 10^-3 to 10^2 nM/min. The initial ligand concentrations to stimulate the cells with include IL2, IL15, IL7, and IL9, and they cover from 10^-3 to 10^3 nM. Each combination of these initial conditions, when inputted into the model, results in a matrix of size 1000 timepoints by 56 concentration values. In order to interpret the concentration values at different timepoints, we reduced the 56 values to 16, four of which include the amount of surface and endosomal active species per cytokine, six include the surface receptor concentrations, and another six encompass the total receptor amounts (both surface and endosomal). As a result, from running X combinations, the final dataset is a tensor of dimensions X combinations by 1000 timepoints by 16 values. From there, further tensor decomposition methods such as tensor factorization were employed in order to achieve intuitive representation and handling of higher dimensional data without losing the structural characteristics of the data itself. 

### Tensor factorization
Before deconstructing the dataset and applying any decomposition methods, the tensor was scaled along the dimension of the 16 values in order to mitigate the different scales available in the varied results. This z-scoring ensures the data is put on a similar scale during analysis. 
Tensor factorization analysis was performed using a python package called TensorLy. Given our tensor of order 3, we applied Canonical Polyadic Decomposition (also known as CP or PARAFAC) to obtain one matrix (or factor) per mode of the tensor. This decomposition follows the scheme that just as a matrix can be expressed as the sum of the outer product of two vectors, a tensor of order 3 can be expressed as the sum of the ‘R’ outer product of three vectors. Each vector then belongs to a decomposition matrix, called a factor, along the corresponding dimension. From there, each of the three factor matrices will have R columns. 
 
In order to determine the number of components (R) necessary to explain 90 percent of the variance in the original tensor, an R2X metric was computed and compared for decomposing the tensor using different number of components. This method first required to decompose the tensor into the sum of R vector products utilizing tensorly.parafac, which results in 3 factor matrices whose second dimension is R. Following, the tensor was reconstructed with tensorly.kruskal_to_tensor by multiplying the factors out. Then, the R2X was computed using
R^2X = 1 - (variance(X_r-X)/variance(X)).
Here, Xr is the reconstructed tensor and X is the original tensor. By iterating over different possible R component numbers, the percent explained variance is obtained. 
 
From the R2X values, we determined that 8 components were needed to explain 90 percent of the variance. The original tensor was then decomposed to 8 components and the columns within each factor matrix were plotted against each other to help visualize correlation and relationships between variables. 

### Correlation Coefficients
As mentioned above, we generated the tensor using 10 input variables. To calculate the correlation coefficients between these input variables and the 8 components to which the tensor was decomposed to, we use the Pearson product-moment correlation coefficient. This is done through numpy.corrcoef between all possible values of each input, i.e. IL2Ra, IL2Rb, gc, IL15Ra, IL7R, and IL9R expression rates and IL2, IL15, IL7, and IL9 ligand concentrations, and each column of the combination decomposition factor matrix. 
