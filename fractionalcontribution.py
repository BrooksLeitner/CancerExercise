#This function takes a numpy array of the MDV of any nutrient and returns the
#atom percent enrichment (APE), or the fraction of the molecule that is actually labelled
#The first value in the MDV is multiplied by 0, so the m0 contribution of the molecule
#does not count towards the APE

#This equation is from Text Box 4 of Fernandez-Garcia et al 2020, DOI: 10.1016/j.tibs.2019.12.002

#This package requires Numpy

import numpy as np

def FractionalAPE(NutrientMDV):
    carbons = np.arange(len(NutrientMDV))
    APE = sum(NutrientMDV*carbons)/(len(NutrientMDV)-1)
    return APE

def FractionalContribution(ProductAPE, SourceAPE):
    Sub_to_Product = ProductAPE/SourceAPE
    return Sub_to_Product
