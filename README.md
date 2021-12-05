# BC-CCS

## "Reclassifying Cancer: Defining tumour cell cycle activity in terms of its tissue of origin in over 13,000 samples" 
### Link to the manuscript Lundberg et al.

### Written by ARIAN LUNDBERG 

All R codes (R statistical software version 4.1.0) to reproduce main figures of 
"Reclassifying Cancer: Defining tumour cell cycle activity in terms of its tissue of origin in over 13,000 samples" 
Lundberg et al. are available in "~/RCodes/"

**The following files were used in the codes.**
 
### Gene.expression_data.RData 
This file is accessible on Dropbox https://www.dropbox.com/sh/2suqtrq6ore9hpy/AACDeYnE4Riz5Hcsva6gynIfa?dl=0

Gene.expression_data.RData is an ExpressionSet data (eSet.CCS: subset of data to only include CCS genes, eSet.Data: used to identify most variable genes, eSet.noTestis: data after exclusion of Testis samples, eSet.final: all the data after final exclusions, cc.Genes: cell-cycle score genes), for the analyses in this project; the data was obtained from 
https://xenabrowser.net/datapages/?hub=https://toil.xenahubs.net:443 

It also contains standard TCGA abbreviations and color codes for the samples and

annuploidy score (Aneuploidy scores and arm calls file - PANCAN_ArmCallsAndAneuploidyScore_092817.txt): https://gdc.cancer.gov/about-data/publications/pancan-aneuploidy

as well as tumour purity data: https://gdc.cancer.gov/about-data/publications/pancanatlas

also PAM50 subtypes: retrieved from cBioportal

### mutations_data.RData 
This file is accessible on Dropbox https://www.dropbox.com/sh/2suqtrq6ore9hpy/AACDeYnE4Riz5Hcsva6gynIfa?dl=0

Gene.expression_data.RData contains mutational information of the samples which were adjusted for number of tumours for each cancer type. 

### Contact

Please contact Arian Lundberg: **arian.lundberg@ucsf.edu** if you have any furthur questions. 

Sincerely,
Arian Lundberg
