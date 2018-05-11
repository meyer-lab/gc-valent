# Methods

All analysis was implemented in Python, and can be found at <https://github.com/meyer-lab/type-I-ckine-model>, release 1.0 (doi: [00.0000/arc0000000](https://doi.org/doi-url)).



## Model

### Base model

Cytokine binding to receptors was modeled using a kinetic model.  The ligands of interest were four cytokines: IL2, IL15, IL7, and IL9. Each ligand can bind to the common gamma chain receptor (gc) and a private family of receptors specific to each ligand. IL2's private receptors were IL2Ra and IL2Rb; IL15's private receptors were IL15Ra and IL2Rb; IL7's private receptor was IL7Ra; IL9's private receptor was IL9R. After a ligand binds to a receptor, it forms a complex that can bind to other free receptors to make larger complexes. Every complex has exactly one ligand and has at least one receptor. Each of these dimerization steps had forward and reverse reaction rate constants based on the dissociation constants of the complexes involved.

Our model also accounted for endosomal trafficking. All species were allowed to be endocytosed, recycled, and degraded. All species  Binding events were allowed to occur inside the endosome. The surface and endosomal ligand concentrations were assumed to be equal for each cytokine. The number densities of each complex and receptor were different between the cell surface and endosome. 

The rate of change for each free receptor and complex was modeled through ordinary differential equations (ODEs). There were 26 ODEs that represented the rate of change for free receptors and complexes on the cell surface, 26 ODEs that represented the rate of change for free receptors and complexes in the endosome, and 4 ODEs that represented the rate of change of ligand concentration. For a complex (C) that forms through the binding of ligand (L) to free receptor (R) with a forward reaction rate (kf) and a reverse reaction rate (kr), the ODE would be:

dC/dt=(kf*L*R)-(kr*C) 

In this example, dC/dt is in units of number / cell * min, L is in units of molarity, R is in number / cell, C is in number / cell, kf is in 1 / nM * min, and kr is in 1 / min. For our model, all of the free receptors and complexes were measured in units of number per cell and all ligands were measured in units of concentration (nM). Due to these unit choices for our species, the rate constants for ligand binding to a free receptors had units of 1 / nM * min, rate constants for the forward dimerization of free receptor to complex had units of 1 / min * number per cell, and the dissociation of a complex into another complex and free receptor had units of 1 / min.

Type I cytokine signaling follows the JAK-STAT signaling pathway which is initiated when two JAK subunits come in contact with one another. JAK proteins are found on the intracellular regions of the gc, IL2Rb, IL7Ra and IL9R receptors; therefore all complexes which contained at least two of those receptors were deemed to be active species.

unit tests?

ODE solver? ... absolute tolerance of 1E-9 and relative tolerance of 1E-12.


### Parameters and assumptions

All ligand-receptor binding processes had a rate constant of ‘kfbnd’ which was assumed to be on rate of 10^7 M-1 sec-1; the only exception to this assumption was the binding of IL2 to gc which was assumed to be on rate of 10^6 M-1 sec-1 (Voss, et al (1993). PNAS. 90, 2428–2432). All forward reaction rates for dimerization were represented by ‘kfwd’. All reverse reaction rates were unique. Each forward and reverse reaction rate was related through a dissociation constant. Many of these dissociation constants were found in published surface plasmon resonance and isothermal calorimetry experiments.

[insert a table of dissociation constants?]
IL2:
k1rev / kfbnd = 10 nM (doi:10.1016/j.jmb.2004.04.038)
k2rev / kfbnd = 144 nM (doi:10.1016/j.jmb.2004.04.038)
k3rev / k3fwd = 50000 nM (doi:10.1016/j.jmb.2004.04.038 -- says that there is no binding of IL2 to gc which means Kd is large)
k10rev / kfwd = 12 nM (doi:10.1016/j.jmb.2004.04.038)
k11rev / kfwd = 63 nM (doi:10.1016/j.jmb.2004.04.038)
k5rev / kfwd = 1.5 nM (doi:10.1016/j.jmb.2004.04.038) -- not yet included in code

IL15:
k13rev / kfbnd = 0.065 nM (multiple papers suggesting 30-100 pM)
k14rev / kfbnd = 438 nM (doi:10.1038/ni.2449)
k15rev / kfbnd = 50000 nM (reference unknown) // gc-IL15 binding very weak (> 50 uM)

IL7:
k25rev / kfbnd = 59 nM (DOI:10.1111/j.1600-065X.2012.01160.x)
k26rev / kfbnd = 50000 nM (reference unknown) // General assumption that cytokine doesn't bind to free gc

IL9:
k30rev / kfbnd = 50000 nM (reference unknown) // cytokine doesn't bind to free gc

Not dissociation constants but measured in literature:
k10rev / kfwd = 12 nM (doi:10.1016/j.jmb.2004.04.038)
k11rev / kfwd = 63 nM (doi:10.1016/j.jmb.2004.04.038)



We used detailed balance to eliminate 12 unknown rate constants. Detailed balance loops were based on formation of full complexes and initial assembly.

The rate of endocytosis is quantified by a rate constant of ‘activeEndo’ for active complexes and ‘endo’ for all other species. The fraction of all endosomal species sent to lysosomes is ‘sortF’. All endosomal species not sent to lysosomes are recycled back to the cell surface. The rate constants to quantify degradation and recycling are ‘kDeg’ and ‘kRec’, respectively. There is no autocrine ligand produced by the cells. Receptors can be synthesized by the cells and placed on the cell surface; these receptor synthesis rates are specific to each receptor. The volume of the entire endosome was 623 (units unknown; got info from TAM paper). The surface area of the  endosome is half the size of the cell surface (got info from same TAM paper).


### Model fitting

We used Markov chain Monte Carlo to fit the unknown parameters in our model to experimental data. The data we used was from experiments performed by Ring, et al. (https://doi.org/10.1038/ni.2449). The three data sets we fit our model to were p-STAT5 activity levels when exposed to varying concentrations of IL-2, p-STAT5 activity levels when exposed to varying concentrations of IL-15, and percent of total IL2Rb that is located on the cell surface over the course of 100 minutes. YT-1 cells were used for all three data-sets and each experiment was repeated with YT-1 cells absent of IL2Ra. The surface IL2Rb experiments were performed for IL2 and IL15 concentrations of 1 nM and 500 nM (4 total per cell type). In order to fit our model to the data, we had to create functions that would calculate percent of maximal p-STAT5 activity and percent of surface IL2Rb under various cellular conditions. These functions allowed us to compare our base model to the experimental data directly.

When fitting our model to the data, we used PyMC model that incorporates Bayesian Statistics in order to store what the likelihood of the model is for a given point. Within our PyMC model, we assumed a lognormal distribution and a standard deviation of 1 for all reaction rate parameters and trafficking parameters; the only exception to these distributions was ‘sortF’ in which we assumed a beta distribution with shape parameters of alpha=2 and beta=7. Executing this fitting process yielded distributions for the likelihood of each unknown parameter and the residuals produced at each experimental data-point that was fit. 
