---
author-meta:
- Brian Orcutt-Jahns
- Zoe S. Kim
- Scott M. Carlson
- Aaron S. Meyer
date-meta: '2020-02-17'
header-includes: '<!--

  Manubot generated metadata rendered from header-includes-template.html.

  Suggest improvements at https://github.com/manubot/manubot/blob/master/manubot/process/header-includes-template.html

  -->

  <meta name="dc.format" content="text/html" />

  <meta name="dc.title" content="Coordinate Optimization of Cytokines" />

  <meta name="citation_title" content="Coordinate Optimization of Cytokines" />

  <meta property="og:title" content="Coordinate Optimization of Cytokines" />

  <meta property="twitter:title" content="Coordinate Optimization of Cytokines" />

  <meta name="dc.date" content="2020-02-17" />

  <meta name="citation_publication_date" content="2020-02-17" />

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
title: Coordinate Optimization of Cytokines
...






<small><em>
This manuscript
was automatically generated
on February 17, 2020.
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

- A
- B
- C
- D

## Author Summary

Text.


## Introduction

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

![**B1.** (A)](./output/figureB1.svg){#fig:B1}

![**B2.** (A)](./output/figureB2.svg){#fig:B2}

![**B3.** (A)](./output/figureB3.svg){#fig:B3}

![**B4.** (A)](./output/figureB4.svg){#fig:B4}

![**B5.** (A)](./output/figureB5.svg){#fig:B5}

## Discussion



## Materials and Methods

All analysis was implemented in Python, and can be found at <https://github.com/meyer-lab/type-I-ckine-model>, release 1.0.


## Acknowledgements

This work was supported by NIH DP5-OD019815 to A.S.M. **Competing financial interests:** The authors declare no competing financial interests.


## Author contributions statement

A.S.M. conceived of the study. A.C.W., A.M.F., A.S.M, and Z.S.K. performed the computational analysis. All authors helped to design experiments and/or analyze the data.


## References {.page_break_before}

<!-- Explicitly insert bibliography here -->
<div id="refs"></div>
