---
title: Tensor Factorization Maps Computation Performed by the Common Gamma Chain Receptors
author:
- name: Adam Weiner
  affilnum: a
- name: Ali Farhat
  affilnum: a
- name: Aaron S. Meyer
  affilnum: a,b
keywords: [IL2, IL15, IL21, IL4, IL7, IL9, common gamma chain, cytokines, receptors, immunology, T cells, NK cells]
affiliation:
- name: Department of Bioengineering, Jonsson Comprehensive Cancer Center, Eli and Edythe Broad Center of Regenerative Medicine and Stem Cell Research; University of California, Los Angeles
  key: a
- name: Contact info
  key: b
bibliography: ./Manuscript/References.bib
abstract: Many receptor families exhibit both pleitropy and redundancy in their regulation, with multiple ligands, receptors, and responding cell populations. The multivariate nature of these systems confounds intuition about therapeutic manipulation. The common γ-chain cytokine receptor dimerizes with complexes of the cytokines interleukin (IL)2, IL-4, IL-7, IL-9, IL-15, and IL-21 and their corresponding "private" receptors. These γ-chain cytokines have accrued broad interest as potential immune therapies because they potently modulate immune cell population types. Here, we build a reaction model for the diverse ligand-receptor interactions of common γ-chain cytokines enabling quantitative predictions of response. Using this quantitative model, we employ tensor factorization to build a quantitative map of regulation across the family. These results map the emergent behavior of the common γ-chain cytokines, and demonstrate an approach to generating interpretable guidelines for manipulation of complex receptor families.
link-citations: true
csl: ./Manuscript/Templates/nature.csl
---

# Summary points

- A mechanism-based binding model of the γ-chain cytokines successfully captures responses to IL-2, IL-15, IL-4, and IL-7
- Two
- Canonical polyadic decomposition maps responses across cell populations, receptors, cytokines, and dynamics

# Introduction

Cytokines are cell signaling proteins responsible for cellular communication within the immune system. The common γ chain (γc) cytokine receptors, with interleukins (IL)-2, 4, 7, 9, 15, and 21 as ligands, are integral for modulating both innate and adaptive immune responses and have accrued broad interest for their potential as immune therapies [@Rochman_2009]. γc ligands bind to private receptors specific to a given cytokine before interacting with the common γc receptor to induce signaling [@Walsh2010]. Studies showing the γc family’s involvement in immune system regulation have led researchers to believe that these ligands and receptors can be leveraged in the development of autoimmune, immunodeficient, and cancer therapies [@Pulliam2015]. Presence of the γc receptor has been implicated in spontaneous and induced lymphoproliferation through γc overexpression and knockout studies, indicating that immune cell populations could be selectively expanded or repressed [@Amorosi3304 @Vigliano2012]. Moreover, γc receptor mutations have been shown to cause Severe Combined Immunodeficiency (SCID) as a defective γc receptor prevents T and NK cells from maturing properly [@Wang9542]. 
In cancer immunotherapy, the antitumor effects of γc cytokines have been explored through clinical trials as these cytokines can show a synergistic increase in the effectiveness of combinatorial drugs.  The γc cytokines’ ability to regulate lymphocytes can impact both solid and hematological tumors, though there are many drawbacks to their therapeutic use due to the onset of severe side effects such as concomitant over-expansion of certain T cell sub-populations [@Pulliam2015]. Nonetheless, regulation of these cytokines is poorly understood due to the complex binding activation mechanism of their receptors, wherein a single common receptor (γc) must complex with a series of private receptors [@Walsh2010]. Further, any intervention imparts effects across multiple distinct cell populations, with each population having a unique response as a consequence of their broad receptor expression [@ring_mechanistic_2012] [@Cotarira17]. 
In addition to therapeutic use of the natural ligands, engineered ligands have been produced with potentially beneficial properties. The most common approach is to develop mutant ligands with altered binding affinities for specific receptors. Mutant IL-2 libraries have been created to see if certain mutant ligands have binding affinities that increase or decrease overall T cell activity [@Berndt1994] [@Collins7709]. Researchers addressing immunodeficiency and cancer immunotherapy issues have shown that engineering a mutant IL-2 with a higher binding affinity for one of its receptors, IL-2Rβ, induces greater cytotoxic T cell proliferation (thus improved antitumor responses) and proportionally less regulatory T cell expansion [@Levin2012]. This behavior can be understood considering IL-2’s typical mode of action, for which cells become more sensitive to upon the expression of its higher affinity receptor, IL-2Rα [@ring_mechanistic_2012]. Bypassing the functional requirement of IL-2 for IL-2Rα expression, researchers were able to differentially regulate cells according to their distinctive receptor expressions [@Levin2012]. Other researchers addressing autoimmunity issues have decreased ligand binding affinity by fusing IL-2 to immunoglobulin G1 (IgG1), which selectively expanded regulatory T cell populations. These IL-2-IgG1 ligands also have a longer half-life and expand fewer NK cells compared to native IL-2 [@Bell_2015] [@Peterson_2018].
In light of these observations, current interests in computationally modeling this receptor family have increased, and several studies have revealed insights into predicting its functions and behaviors. Early attempts at mathematically modeling the synergy between IL-2 and IL-4 in B and T cells found quantitative estimates of parameters characterizing proliferative responses to these two cytokines [@BURKE199742]. At the population level, computational models have explained how IL-2 uptake allows regulatory T cells to suppress effector T cell activation. This suppression was shown to rely on cell numbers, surface receptor densities (particularly that of IL-2Rα), available cytokine concentration, and timing [@Feinerman437]. Furthermore, computational models incorporating IL-4, IL-7, and IL-21 have revealed pathway cross-talk due to the hierarchy of affinities for receptor dimerization events involving the common γc receptor [@Gonnordeaal1253]. Nevertheless, these models lack a major component of intracellular signaling and endosomal trafficking. There is downregulation of IL-2Rα and IL-2Rβ when IL-2 is given to cells due to endocytosis of ligand-receptor complexes [@ring_mechanistic_2012] [@Duprez15091988]. This receptor-mediated endocytosis pathway can be extended to all receptors in the γC family [@Lamaze_2001]. Models accounting for trafficking have more accurately predicted ligand depletion and T cell proliferation in the case of IL-2 compared to their non-trafficking counterparts, further proving the importance of considering trafficking [@Fallon2000].
In this paper, we assemble a map of γc cytokine family regulation. We first built a family-wide mathematical model that incorporates both binding and trafficking kinetics. This more comprehensive model allows us to investigate emergent behavior, such as competition between cytokines. This cytokine family is inherently high dimensional—with multiple ligands, cognate receptors, and cells with distinct expression. Using tensor factorization methods, we then mapped the family-wide regulation. This map will prove useful in therapeutic design by showing how native or engineered ligands are addressed to specific immune cell populations based on their receptor expression levels.


- The common gamma chain cytokine family is central to regulation of the immune system
    - Giving IL-15 doses to monkeys caused expansion of NK and T-eff cell lines, making it a potential therapy in cancers and immunodeficiencies. [@Sneller6845]
    - Transgenic mice that overexpressed IL-15 revealed that IL-15 selectively propagates memory T cells by blocking IL-2-induced T cell death. [@Marks-Konczalik11445]
    - gamma chain cytokines have accrued interest as cancer immunotherapeutics due to their ability to regulate lymphocytes. [@Pulliam2015]
    - Mutation in gamma chain causes XSCID. [@NOGUCHI1993147]
    - Cells must express gamma chain in order to undergo lymphoproliferation. [@Amorosi3304]
    - Gamma chain is tied to hematopoietic cell proliferation. High expression correlates with malignancies and knockouts slow the proliferation of these malignancies. [@Vigliano2012]
- Many receptor families are composed of multiple ligands, receptors, and cell populations
	- [@Antebi_BMP]
    - hGH affinities do not correlate well with cell stimulation levels because of endosomal trafficking. [@Haugh2004]
    - EGFRs are being targeted in cancer therapies due to their increased expression in tumor cells. [@CIARDIELLO20031348]
- Efforts to engineer gc behavior and/or gc treatments
    - [@ring_mechanistic_2012]
    - Specific IL-2 amino acid substitutions render it incapable from binding to IL2Ra. [@Collins7709]
    - IL-2 mutant library created and tested to see if any mutations increased binding affinity to receptors and therefore greater T cell activity. [@Berndt1994]
    - "in vitro evolution" was used to create an IL-2 mutant with increased affinity for IL-2Rb. [@Levin2012]
- Efforts to model behavior of cytokine families
    - Model that explains the positive and negative feedback loops of TCR antigen recognition [@FrancoisE888]
    - Modeling that doesn't account for trafficking
        - Model of T cell proliferative response to IL-2 and IL-4. [@BURKE199742]
        - High level model of how IL-2 feedback loops regulate T cell populations. Focus on helper T and regulatory T cells. [@Feinerman437], [@Garcia2012]. 
    - Modeling that accounts for trafficking
        - IL-2 stimulates downregulation of IL2Ra due to endocytosis of IL2R complexes. [@Duprez15091988]
        - Accounting for trafficking allows for more accurate ligand depletion and cell proliferation model in case of IL-2. [@Fallon2000]
        - Trafficking and TCR inhibition of pSTAT5 activation can be used to improve reaction model of IL-2... potentially get rid of this [@Tkach2014TCT]
    - Focused on competition for common gamma chain
        - Showed how high IL2Ra abundnaces lowered the IL-2 EC50 but dampened the signaling responses of IL-7 and IL-15. Their work was mostly experimental (flow cytometry)  [@Cotarira17]
        - Abundance of gamma chain on surface played role in IL-2 mediated activation. [@Ponce2016]
- Transition to paper

