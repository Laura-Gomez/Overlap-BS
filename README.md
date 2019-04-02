# OVERLAP-BS

This program plots the TF binding sites that exist in a set of promoters

## DESCRIPTION

This program plots all binding sites annotated as part of a set of promoters. 

The promoters are grouped into two categories:
 - Promoters that contain at least one binding site that was recovered by an HT technology
 - Promoters that contain at least one  binding site that was not recovered by an HT technology.

## PLOT CONVETIONS

**Each row represents one promoter**

The **TRANSPARENCY LEVEL** indicates the leve of evidence for each BS:
 - Solid color: Strong evidence
 - Transparent color: Weak and no evidence

The **TYPE OF LINE** indicates the annotated effect of the BS:
 - Continuous line: Activator
 - Dashed line: Repressor
 - Dotted line: Dual

The **WIDTH OF LINE** indicates if the BS was found by the HT technology:
 - Thick line: Found BS
 - Thin line: Not found BS


### INPUT

For each input file, a list of columns is provided. All columns must exist. When the information is no available NA is allowed. The mandatory fields are shown in bold letters

**Binding sites from RegulonDB** [regulon.bs]
 - TF ID
 - **TF name**
 - TFBS ID
 - **TFBS left position**
 - **TFBS rigth position** 
 - TFBS strand
 - TFBS associated gene ID
 - TFBS associated TU 
 - **TFBS effect**
 - **TFBS associated promoter**
 - **TFBS distance to TSS**
 - TFBS sequence
 - TFBS evidence
 - **TFBS evidence level**

**HT DATA DIRECTORY** [HT/]
A BS is considered the same as one in RegulonDB if there is a match with the promoter name, the TF name and the TSS.DIST

Each TF must have a csv with the following columns:

 - found
 - ID1
 - TF.NAME
 - ID2
 - TFBS.CENTER
 - TFBS.LEFT
 - TFBS.RIGTH
 - TFBS.STRAND
 - ID3
 - TU
 - EFFECT
 - PROMOTER
 - TSS.DIST
 - SEQ
 - EVIDENCE.LIST
 - EVIDENCE
 - LAST

**TF LIST** [TF.LIST]
A list of TFs of interest.
One TF per line


## OUTPUT [ Plot.pdf ]

A PDF with four plots per TF:
 - Promoters that contain at least one binding site that was recovered by an HT technology:
  * Plot with only BSs with strong evidence
  * Plot with all BSs

 - Promoters that contain at least one  binding site that was not recovered by an HT technology:
  * Plot with only BSs with strong evidence
  * Plot with all BSs


### RUN

 Rscript --vanilla Promoters_Plot.R regulon.bs HT/ TF-LIST Plot.pdf
 
