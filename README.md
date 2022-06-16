# BC-CCS

## "Reclassifying tumour cell cycle activity in terms of its tissue of origin" 

### Written by ARIAN LUNDBERG 

All R codes (R statistical software version 4.1.0) to reproduce the main figures of "Reclassifying tumour cell cycle activity in terms of its tissue of origin" Lundberg et al. are stored in "~/RCodes/"

**The following files were used in the codes.**
 
### Gene.expression_data.RData and mutations_data.RData 
The files are available upon request due to file size limitations on Github. Please contact nick.tobin@ki.se

These files contain (eSet.CCS: a subset of data to only include CCS genes, eSet.Data: used to identify most variable genes, eSet.noTestis: data after exclusion of Testis samples, eSet.final: all the data after final exclusions, cc.Genes: cell-cycle score genes, mutational information of the samples which were adjusted for the number of tumours for each cancer type.), for the analyses in this project; the data was obtained from https://xenabrowser.net/datapages/?hub=https://toil.xenahubs.net:443

It also contains standard TCGA abbreviations and color codes for the samples and
aneuploidy score (Aneuploidy scores and arm calls file - PANCAN_ArmCallsAndAneuploidyScore_092817.txt): https://gdc.cancer.gov/about-data/publications/pancan-aneuploidy as well as tumour purity data: https://gdc.cancer.gov/about-data/publications/pancanatlas

The PAM50 subtypes information was retrieved from cBioportal.

### Contact info

Please feel free to contact me via **arian.lundberg@ucsf.edu** or **nick.tobin@ki.se** if you have any further questions. 

Sincerely,
Arian Lundberg
