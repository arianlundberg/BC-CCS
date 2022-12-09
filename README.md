# BC-CCS

## "Reclassifying tumour cell cycle activity in terms of its tissue of origin" 

### Written by ARIAN LUNDBERG 

The R codes for reproducing the main figures of "Reclassifying tumour cell cycle activity in terms of its tissue of origin" by Lundberg et al. are available in the "~/RCodes/" directory. 

The required Gene.expression_data.RData and mutations_data.RData files are available upon request by emailing nick.tobin@ki.se. 

These files contain the following data for the analyses in this project:

eSet.CCS: a subset of data that only includes CCS genes

eSet.Data: data used to identify the most variable genes

eSet.noTestis: data after exclusion of Testis samples

eSet.final: all data after final exclusions

cc.Genes: cell-cycle score genes, and mutational information of the samples that has been adjusted for the number of tumours for each cancer type

The data was obtained from https://xenabrowser.net/datapages/?hub=https://toil.xenahubs.net:443. 
It also contains standard TCGA abbreviations and color codes for the samples and aneuploidy score (Aneuploidy scores and arm calls file - PANCAN_ArmCallsAndAneuploidyScore_092817.txt): https://gdc.cancer.gov/about-data/publications/pancan-aneuploidy, as well as tumour purity data: https://gdc.cancer.gov/about-data/publications/pancanatlas. 

The PAM50 subtypes information was retrieved from cBioportal.

Please feel free to contact me via arian.lundberg@ucsf.edu or nick.tobin@ki.se if you have any further questions.

Sincerely,
Arian Lundberg
