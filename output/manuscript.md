---
author-meta:
- Brian Orcutt-Jahns
- Zoe S. Kim
- Scott M. Carlson
- Aaron S. Meyer
date-meta: '2020-04-16'
header-includes: '<!--

  Manubot generated metadata rendered from header-includes-template.html.

  Suggest improvements at https://github.com/manubot/manubot/blob/master/manubot/process/header-includes-template.html

  -->

  <meta name="dc.format" content="text/html" />

  <meta name="dc.title" content="Fitting a mechanistic model of IL-2/IL-15 response to single cell measurements using moment propagation" />

  <meta name="citation_title" content="Fitting a mechanistic model of IL-2/IL-15 response to single cell measurements using moment propagation" />

  <meta property="og:title" content="Fitting a mechanistic model of IL-2/IL-15 response to single cell measurements using moment propagation" />

  <meta property="twitter:title" content="Fitting a mechanistic model of IL-2/IL-15 response to single cell measurements using moment propagation" />

  <meta name="dc.date" content="2020-04-16" />

  <meta name="citation_publication_date" content="2020-04-16" />

  <meta name="dc.language" content="en-US" />

  <meta name="citation_language" content="en-US" />

  <meta name="dc.relation.ispartof" content="Manubot" />

  <meta name="dc.publisher" content="Manubot" />

  <meta name="citation_journal_title" content="Manubot" />

  <meta name="citation_technical_report_institution" content="Manubot" />

  <meta name="citation_author" content="Brian Orcutt-Jahns" />

  <meta name="citation_author_institution" content="Department of Bioengineering, University of California, Los Angeles" />

  <meta name="citation_author" content="Zoe S. Kim" />

  <meta name="citation_author_institution" content="Department of Bioengineering, University of California, Los Angeles" />

  <meta name="citation_author" content="Scott M. Carlson" />

  <meta name="citation_author_institution" content="Visterra, Inc." />

  <meta name="citation_author" content="Aaron S. Meyer" />

  <meta name="citation_author_institution" content="Department of Bioengineering, University of California, Los Angeles" />

  <meta name="citation_author_institution" content="Department of Bioinformatics, University of California, Los Angeles" />

  <meta name="citation_author_institution" content="Jonsson Comprehensive Cancer Center, University of California, Los Angeles" />

  <meta name="citation_author_institution" content="Eli and Edythe Broad Center of Regenerative Medicine and Stem Cell Research, University of California, Los Angeles" />

  <meta name="citation_author_orcid" content="0000-0003-4513-1840" />

  <meta name="twitter:creator" content="@aarmey" />

  <meta property="og:type" content="article" />

  <meta property="twitter:card" content="summary_large_image" />

  <link rel="icon" type="image/png" sizes="192x192" href="https://manubot.org/favicon-192x192.png" />

  <link rel="mask-icon" href="https://manubot.org/safari-pinned-tab.svg" color="#ad1457" />

  <meta name="theme-color" content="#ad1457" />

  <!-- end Manubot generated metadata -->'
keywords:
- cytokines
- IL-2
lang: en-US
title: Fitting a mechanistic model of IL-2/IL-15 response to single cell measurements using moment propagation
...






<small><em>
This manuscript
was automatically generated
on April 16, 2020.
</em></small>

## Authors



+ **Brian Orcutt-Jahns**<br><br>
  <small>
     Department of Bioengineering, University of California, Los Angeles
  </small>

+ **Zoe S. Kim**<br><br>
  <small>
     Department of Bioengineering, University of California, Los Angeles
  </small>

+ **Scott M. Carlson**<br><br>
  <small>
     Visterra, Inc.
  </small>

+ **Aaron S. Meyer**<br>
    ORCID
    [0000-0003-4513-1840](https://orcid.org/0000-0003-4513-1840)
    · Github
    [aarmey](https://github.com/aarmey)
    · twitter
    [aarmey](https://twitter.com/aarmey)<br>
  <small>
     Department of Bioengineering, University of California, Los Angeles; Department of Bioinformatics, University of California, Los Angeles; Jonsson Comprehensive Cancer Center, University of California, Los Angeles; Eli and Edythe Broad Center of Regenerative Medicine and Stem Cell Research, University of California, Los Angeles
  </small>



## Abstract {.page_break_before}

Text.

## Summary points

- Propogation of moments generally enables prediction and fitting between diverse models and single cell measurements.
- B
- C
- D

## Author Summary

Text.


## Introduction








New statistical approaches are needed to link single cell measurements to molecular pathway models. The nature of the statistical approach is dependent upon the underlying source of single cell variation. For example, stochastic differential equations effectively model processes wherein variation arises from intrinsic fluxuations (CITE). However, much of cell heterogeneity in cell response arises from variation that is extrinsic to the pathway being modeled, or even extrinsic to the cell due to variation in the extracellular environment [@KRuDS281]. In this situation, the pathway of interest can be treated as deterministic, but with varying inputs. One successful approach has been to separately fit a model to each cell independently [@bKxRf27N]. However, there are two limitations of this strategy: First, one usually wishes to fit both shared and varying rate parameters within a population of cells. For example, the binding affinity of a protein interaction should be shared among cells, while the level of protein expression might vary. Second, this approach requires sufficient measurements from the same cell. More often, one can easily collect many measurements in the same population of cells, but only make a limited set of measurements from an individual cell.






Here, we present a flexible statistical approach for fitting mechanistic or data-driven models to single cell measurements. This method, which relies on moment propagation, allows for integration of single cell measurements in the same populations and of extrinsic and intrinsic variability. We apply this technique to identify how receptor variation translates to variation in response to the common gamma chain cytokine receptor cytokines. These results enable rational cytokine engineering accounting for the widespread heterogeneity within cytokine-responsive cell populations. Simultaneously, our statistical approach can be applied to many molecular programs with prominent sources of cell variability.



### Visterra Notebooks Text

#### (1) Initial assumptions

To review, we've been assuming (1) that binding occurs identically in the endosome as compared to the surface, (2) constant rate of receptor expression, (3) active receptor complexes are endocytosed at a higher rate, and (4) active complexes in the endosome still have signaling capacity. Here, we've adjusted the model to make all binding 5-fold weaker in the endosome as compared to the surface. We're also using our fit expression levels of IL2Ra, IL2Rb, and gc for the YT-1 cell line (3.9, 0.7, 1.7 receptors/min/cell, respectively). The binding rate parameters used here are from fitting using our older endosome binding assumption, but we're in the process of running fitting again to update these. I doubt the results here will change meaningfully. We're also adding a sigmoidal relationship between number of active complexes and STAT5 activation like we discussed. However, we haven't yet run the fitting with this, and so the values here are all proportional to active receptor complexes.

#### Exchanging IL2Rb for IL2Ra affinity in CD25+ cells

A reminder this is using our inferred receptor expression levels from working with the CD25+ YT-1 cells. On the x-axis we're varying IL2Ra affinity (lower values mean tighter binding). The colors indicate differing IL2Rb affinities (again, higher values indicate weaker binding). The y-axis is the IL2 concentration at which you reach half-maximal activation in the wild-type cells (i.e. ~20 active receptor complexes). The black line shows the level that corresponds to wild type.

As expected, this shows there is indeed a trade-off, and one can reduce the affinity to IL2Rb after increasing it to IL2Ra. For example, the green line shows a 10X increase in IL2Ra affinity allows a 5X decrease in IL2Rb affinity to preserve the same threshold.

#### Exchange does change other IL2 response metrics

However, exactly how much you adjust each affinity is dependent upon how you quantify IL2 response. The example I gave is plotted below, where you can see you have the same half-maximal concentration, but level of maximal activation is now lower.

I know we also discussed looking at how fast ligand is consumed in cases such as this. That's very straightforward to calculate within the model but I've left anything about ligand consumption out for now.

#### IL2Rb affinity adjustment with variable CD25 expression

A second concern we discussed is separating CD25+/CD25- cells. Here, the lines indicate cells with reduced CD25 expression (relative to the cells above), then the x-axis is IL2Rb affinity. The y-axis is the same half-maximal activation threshold. Because the lines remain similarly spaced, this seems to indicate adjusting the IL2Rb affinity doesn't give you any better/worse discrimination between CD25+/CD25- cells. The difference in half maximal concentration for 10% CD25 and no CD25 is less than 2-fold. **Note this is for wild-type IL2Ra affinity.**

#### Same plot, with 10X higher IL2Ra affinity

Here is the same plot, but with a 10-fold higher IL2Ra affinity. Now there's much better separation of the cell populations—almost 100-fold difference in half-max for CD25+ vs CD25-, and just under 10-fold for 25% levels of CD25. Again IL2Rb affinity doesn't influence the level of specificity between cell populations. Note that modulating IL2Rb affinity could still have benefits in ligand consumption.

#### (2) Dose response when changing IL2Ra affinity alone

This plot shows how a CD25+ cell's dose-response behavior shifts when you increase the affinity of IL2 for IL2Ra. Lines with smaller Kd's correspond to higher binding affinity.

#### Dose response when changing IL2Rb affinity alone (for wt IL2Ra affinity)

This plot shows how a CD25+ cell's dose-response behavior shifts when you decrease the affinity of IL2 for IL2Rb. Lines with larger Kd's correspond to lower binding affinity.

#### Dose response when changing IL2Rb affinity alone (for 10X higher IL2Ra affinity)

#### (3) IL2 Degradation when changing IL2Ra affinity alone

This plot shows how a CD25+ cell's IL2-degradation behavior shifts when you increase the affinity of IL2 for IL2Ra. Lines with smaller Kd's correspond to higher binding affinity.

#### IL2 Degradation when changing IL2Rb affinity alone (for wt IL2Ra affinity)

This plot shows how a CD25+ cell's IL2-degradation behavior shifts when you decrease the affinity of IL2 for IL2Rb. Lines with larger Kd's correspond to lower binding affinity.

#### IL2 Degradation when changing IL2Rb affinity alone (for 10X higher IL2Ra affinity)


## Results

### Variation in receptor abundance drives variation in response

![**B1.** (A)](./output/figureB1.svg){#fig:B1}

### Quantifying receptor variability

![**B2.** (A)](./output/figureB2.svg){#fig:B2}

### Quantifying variation in pSTAT5 response

![**B3.** (A)](./output/figureB3.svg){#fig:B3}

### Moment propagation to link molecular and response variation

![**B4.** (A)](./output/figureB4.svg){#fig:B4}

![**B5.** (A)](./output/figureB5.svg){#fig:B5}

## Discussion



## Materials and Methods

All analysis was implemented in Python, and can be found at <https://github.com/meyer-lab/type-I-ckine-model>, release 1.0.

### Model

#### Base model

Cytokine (IL-2, -4, -7, -9, -15, & -21) binding to receptors was modeled using ordinary differential equations (ODEs). IL-2 and -15 each had two private receptors, one being a signaling-deficient α-chain (IL-2Rα & -15Rα) and the other being signaling-competent IL-2Rβ. The other four cytokines each had one signaling-competent private receptor (IL-7Rα, -9R, -4Rα, & -21Rα). JAK-STAT signaling is initiated when JAK-binding motifs are brought together. JAK binding sites are found on the intracellular regions of the γ~c~, IL-2Rβ, IL-4Rα, and IL-7Rα receptors; therefore all complexes which contained two signaling-competent receptors were deemed to be active species. Ligands were assumed to first bind a private receptor and then can dimerize with other private receptors or γ~c~ thereafter. Direct binding of ligand to γ~c~ was not included due to its very weak or absent binding [@De7G6qVz].

In addition to binding interactions, our model incorporated receptor-ligand trafficking. Receptor synthesis was assumed to occur at a constant rate. The endocytosis rate was defined separately for active (k~endo,a~) and inactive (k~endo~) receptors. f~sort~ fraction of species in the endosome were ultimately trafficked to the lysosome, and active species in the endosome had a sorting fraction of 1.0. All endosomal species not sent to lysosomes were recycled back to the cell surface. The lysosomal degradation and recycling rate constants were defined as k~deg~ and k~rec~, respectively. We assumed no autocrine ligand was produced by the cells. We assumed an endosomal volume of 10 fL and endosomal surface area half that of the plasma membrane [@vWk0iSTe]. All binding events were assumed to occur with 5-fold greater disassociation rate in the endosome due to its acidic pH [@WJgrtb3L].

Free receptors and complexes were measured in units of number per cell and soluble ligands were measured in units of concentration (nM). Due to these unit choices for our species, the rate constants for ligand binding to a free receptors had units of nM^-1^ min^-1^, rate constants for the forward dimerization of free receptor to complex had units of cell min^-1^ number^-1^. Dissociation rates had units of min^-1^. All ligand-receptor binding processes had an assumed forward rate (k~bnd~) of 10^7^ M^-1^ sec^-1^. All forward dimerization reaction rates were assumed to be identical, represented by k~fwd~. Reverse reaction rates were unique. The experimentally-derived affinities of 59 nM [@22AxiXSN] was used for IL-7 binding to its cognate private receptors, respectively. IL-2 and -15 were assumed to have affinities of 10 nM and 0.065 nM for their respective α-chains [@2F51YH7p; @1BgE7VeUq; @1EcAaOZOj], and affinities of 144 nM and 438 nM for their respective β-chains [@2F51YH7p]. Rates k~5~, k~10~, and k~11~ were set to their experimentally-determined dissassociation constants of 1.5, 12, and 63 nM [@2F51YH7p].

Initial values were calculated by assuming steady-state in the absence of ligand. Differential equation solving was performed using DifferentialEquations.jl in Julia (CITE). Model sensitivities were calculated by forward autodifferentiation using ForwardDiff.jl. A model solving tolerance of 10^-9^ was used throughout. We used unit tests for conservation of mass, equilibrium, and detailed balance to help ensure model correctness (CITE detailed balance ms).

#### Moment Propagation




#### Model fitting






### Experimental Methods

#### Receptor abundance quantitation

Cryopreserved PBMCs (ATCC, PCS-800-011, lot#81115172) were thawed to room temperature and slowly diluted with 9 mL pre-warmed RPMI-1640 medium (Gibco, 11875-093) supplemented with 10% fetal bovine serum (FBS, Seradigm, 1500-500, lot#322B15). Media was removed, and cells washed once more with 10 mL warm RPMI-1640 + 10% FBS. Cells were brought to 1.5x10^6^ cells/mL, distributed at 250,000 cells per well in a 96-well V-bottom plate, and allowed to recover 2 hrs at 37℃ in an incubator at 5% CO2. Cells were then washed twice with PBS + 0.1% BSA (PBSA, Gibco, 15260-037, Lot#2000843) and suspended in 50 µL PBSA + 10% FBS for 10 min on ice to reduce background binding to IgG.

Antibodies were diluted in PBSA + 10% FBS and cells were stained for 1 hr at 4℃ in darkness with a gating panel (Panel 1, Panel 2, Panel 3, or Panel 4) and one anti-receptor antibody, or an equal concentration of matched isotype/fluorochrome control antibody. Stain for CD25 was included in Panel 1 when CD122, CD132, CD127, or CD215 was being measured (CD25 is used to separate T~reg~s from other CD4+ T cells).

Compensation beads (Simply Cellular Compensation Standard, Bangs Labs, 550, lot#12970) and quantitation standards (Quantum Simply Cellular anti-Mouse IgG or anti-Rat IgG, Bangs Labs, 815, Lot#13895, 817, Lot#13294) were prepared for compensation and standard curve. One well was prepared for each fluorophore with 2 µL antibody in 50 µL PBSA and the corresponding beads. Bead standards were incubated for 1 hr at room temperature in the dark.

Both beads and cells were washed twice with PBSA. Cells were suspended in 120 µL per well PBSA, and beads to 50 µL, and analyzed using an IntelliCyt iQue Screener PLUS with VBR configuration (Sartorius) with a sip time of 35 and 30 secs for cells and beads, respectively. Antibody number was calculated from fluorescence intensity by subtracting isotype control values from matched receptor stains and calibrated using the two lowest binding quantitation standards. T~reg~ cells could not be gated in the absence of CD25, so CD4+ T cells were used as the isotype control to measure CD25 in T~reg~ populations. Cells were gated as shown in XXX. Measurements were performed using four independent staining procedures over two days. Separately, the analysis was performed with anti-receptor antibodies at 3x normal concentration to verify that receptor binding was saturated. Replicates were summarized by geometric mean.

#### pSTAT5 Measurement of IL-2 and -15 Signaling in PBMCs

Human PBMCs were thawed, distributed across a 96-well plate, and allowed to recover as described above. IL-2 (R&D Systems, 202-IL-010) or IL-15 (R&D Systems, 247-ILB-025) were diluted in RPMI-1640 without FBS and added to the indicated concentrations. To measure pSTAT5, media was removed, and cells fixed in 100 µL of 10% formalin (Fisher Scientific, SF100-4) for 15 mins at room temperature. Formalin was removed, cells were placed on ice, and cells were gently suspended in 50 µL of cold methanol (-30℃). Cells were stored overnight at -30℃. Cells were then washed twice with PBSA, split into two identical plates, and stained 1 hr at room temperature in darkness using antibody panels 4 and 5 with 50 µL per well. Cells were suspended in 100 µL PBSA per well, and beads to 50 µL, and analyzed on an IntelliCyt iQue Screener PLUS with VBR configuration (Sartorius) using a sip time of 35 seconds and beads 30 seconds. Compensation was performed as above. Populations were gated as shown in XXX, and the median pSTAT5 level extracted for each population in each well.

#### Recombinant proteins

IL-2/Fc fusion proteins were expressed using the Expi293 expression system according to manufacturer instructions (Thermo Scientific). Proteins were as human IgG1 Fc fused at the N- or C-terminus to human IL-2 through a (G4S)4 linker. C-terminal fusions omitted the C-terminal lysine residue of human IgG1. The AviTag sequence GLNDIFEAQKIEWHE was included on whichever terminus did not contain IL-2. Fc mutations to prevent dimerization were introduced into the Fc sequence [@cCd4xfrS]. Proteins were purified using MabSelect resin (GE Healthcare). Proteins were biotinylated using BirA enzyme (BPS Biosciences) according to manufacturer instructions, and extensively buffer-exchanged into phosphate buffered saline (PBS) using Amicon 10 kDa spin concentrators (EMD Millipore). The sequence of IL-2Rβ/γ Fc heterodimer was based on a reported active heterodimeric molecule (patent application US20150218260A1), with the addition of (G4S)2 linker between the Fc and each receptor ectodomain. The protein was expressed in the Expi293 system and purified on MabSelect resin as above. IL2-Rα ectodomain was produced with C-terminal 6xHis tag and purified on Nickel-NTA spin columns (Qiagen) according to manufacturer instructions. 

#### Octet binding assays

Binding affinity was measured on an OctetRED384 (ForteBio). Briefly, biotinylated monomeric IL-2/Fc fusion proteins were uniformly loaded to Streptavidin biosensors (ForteBio) at roughly 10% of saturation point and equilibrated for 10 minutes in PBS + 0.1% bovine serum albumin (BSA). Association time was up to 40 minutes in IL-2Rβ/γ titrated in 2x steps from 400 nM to 6.25 nM, or IL-2Rα from 25 nM to 20 pM, followed by dissociation in PBS + 0.1% BSA. A zero-concentration control sensor was included in each measurement and used as a reference signal. Assays were performed in quadruplicate across two days. Binding to IL-2Rα did not fit to a simple binding model so equilibrium binding was used to determine the K~D~ within each assay. Binding to IL-2Rβ/γ fit a 1:1 binding model so on-rate (k~on~), off-rate (k~off~) and K~D~ were determined by fitting to the entire binding curve. Kinetic parameters and K~D~ were calculated for each assay by averaging all concentrations with detectable binding signal (typically 12.5 nM and above).


## Acknowledgements

This work was supported by a research agreement with Visterra, Inc. **Competing financial interests:** S.M.C. and C.P. are employees of Visterra Inc.

## Author contributions statement

A.S.M. and S.M.C. conceived of the study. S.M.C. and C.P. performed the PBMC experiments and engineered the IL-2 fusion proteins. A.S.M, B.O.J., and Z.S.K. performed the computational analysis. All authors helped to design experiments and/or analyze the data.


## Author contributions statement

A.S.M. conceived of the study. A.C.W., A.M.F., A.S.M, and Z.S.K. performed the computational analysis. All authors helped to design experiments and/or analyze the data.


## References {.page_break_before}

<!-- Explicitly insert bibliography here -->
<div id="refs"></div>
