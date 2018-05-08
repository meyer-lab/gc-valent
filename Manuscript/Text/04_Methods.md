# Methods

All analysis was implemented in Python, and can be found at <https://github.com/meyer-lab/type-I-ckine-model>, release 1.0 (doi: [00.0000/arc0000000](https://doi.org/doi-url)).



## Model

### Base model

Cytokine binding to receptors was modeled using an equilibrium model of ligand-receptor interaction.  The ligands in this model were the four cytokines: IL2, IL15, IL7, and IL9. The model allows each ligand molecule to bind to the common gamma chain receptor (gc) and a private family of receptors specific to each ligand. The private receptors for each ligand were as follows – IL2 has IL2Ra and IL2Rb; IL15 has IL15Ra and IL2Rb; IL7 has IL7Ra; IL9 has IL9R. All of these ligand-receptor binding processes had a rate constant of ‘kfbnd’ which was assumed to be on rate of 10^7 M-1 sec-1; the only exception to this assumption was the binding of IL2 to gc which was assumed to be on rate of 10^6 M-1 sec-1 (Voss, et al (1993). PNAS. 90, 2428–2432). Once ligands were bound to receptors, they can bind to other free receptors to make larger complexes. Each of these dimerization steps had forward and reverse reaction rate constants based on the dissociation constants of the complexes involved. While the forward reaction rates for dimerization were represented by ‘kfwd’, the reverse reaction rates, known as dissociation rates, were unique. Many of these dissociation constants were found in published surface plasmon resonance and isothermal calorimetry experiments.

The rate of change for each free receptor and complex was modeled through ordinary differential equations (ODEs). There were 26 ordinary differential equations that represented the rate of change for each surface species. For a complex (C) that forms through the binding of ligand (L) to free receptor (R) with a forward reaction rate (kf) and a reverse reaction rate (kr), the ODE would be:

dC/dt=(kf*L*R)-(kr*C) 

In this example, dC/dt is in units of number / cell * min, L is in units of molarity, R is in number / cell, C is in number / cell, kf is in 1 / molarity * min, and kr is in 1 / min. For our model, all of the free receptors and complexes were measured in units of number per cell and all ligands were measured in units of concentration (nM). Due to these unit choices for our species, the rate constants for ligand binding to a free receptors had units of 1 / nM * min (?), rate constants for the forward dimerization of free receptor to complex had units of 1 / min * number per cell, and the dissociation of a complex into another complex and free receptor had units of 1 / min.

We used detailed balance to eliminate unknown rate constants. Detailed balance eliminates unknown rate constants because the network of ligand-receptor and dimerization reactions form “loops”. Within each “loop”, we can express the conversion of one complex to another using two paths and each path involves a unique combination of rate constants. Thus we were able to solve for one unknown rate constant in terms of other rate constants with each loop in our reaction network.

Type I cytokine signaling follows the JAK-STAT signaling pathway which is initiated when two JAK subunits come in contact with one another. Since JAK proteins are found on the intracellular regions of the gc, IL2Rb, IL7Ra and IL9R receptors, all complexes which contained at least two of those receptors were deemed to be active. 

Another component our model accounted for was endosomal trafficking. Each free ligand, free receptor, and complex can move from the cell surface into the endosome. Binding events are allowed to occur inside the endosome. Since the concentration of species in the endosome differs from that at the cell surface, an additional 26 ODEs must be constructed to keep track of all the endosomal ligands, receptors, and complexes. This engulfing process is quantified by a rate constant of ‘activeEndo’ for active complexes and ‘endo’ for all other species. Once in endosomal vesicles, the contents are either sent to lysosomes for degradation or back to the cell surface for recycling. ‘SortF’ represents the fraction of all endosomal species that are sent to lysosomes for degradation; the remaining endosomal species (1 – ‘sortF’) are recycled back to the cell surface. The rate constants to quantify degradation and recycling are ‘kDeg’ and ‘kRec’, respectively. While there is no autocrine ligand produced by the cells, each receptor can be produced by the cells and placed back on the cell surface. Each receptor has its own rate constant associated with the synthesis process.


### Parameters and assumptions

In this model we assumed that there was no autocrine ligand; however, all receptors were synthesized by the cells. All units


### Model fitting

We used Markov chain Monte Carlo to fit the unknown parameters in our model to experimental data. The data we used was from experiments performed by Ring, et al. (https://doi.org/10.1038/ni.2449). The three data sets we fit our model to were p-STAT5 activity levels when exposed to varying concentrations of IL-2, p-STAT5 activity levels when exposed to varying concentrations of IL-15, and percent of total IL2Rb that is located on the cell surface over the course of 100 minutes. YT-1 cells were used for all three data-sets and each experiment was repeated with YT-1 cells absent of IL2Ra. The surface IL2Rb experiments were performed for IL2 and IL15 concentrations of 1 nM and 500 nM (4 total per cell type). In order to fit our model to the data, we had to create functions that would calculate percent of maximal p-STAT5 activity and percent of surface IL2Rb under various cellular conditions. These functions allowed us to compare our base model to the experimental data directly.

When fitting our model to the data, we used PyMC model that incorporates Bayesian Statistics in order to store what the likelihood of the model is for a given point. Within our PyMC model, we assumed a lognormal distribution and a standard deviation of 1 for all reaction rate parameters and trafficking parameters; the only exception to these distributions was ‘sortF’ in which we assumed a beta distribution with shape parameters of alpha=2 and beta=7. Executing this fitting process yielded distributions for the likelihood of each unknown parameter and the residuals produced at each experimental data-point that was fit. 
