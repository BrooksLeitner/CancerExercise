import numpy as np
import scipy

#Function to calculate turnover based on Joanne Kelleher's Linear Regression Modeling
#It takes the Tracer, Natural/Background (unenriched plasma), 
#and Observed/Experimental/Infused (plasma after infusion)
#in the format of mass distribution vectors (m0, m1, m2, m3, etc...) as fractions that add up to 1
#Then infusion rate in the units mg/kg/min. is taken
#These are taken in the format of Python lists

def Turnover_Ra(Tracer, Natural, Observed, infusion_rate):
    
    #Convert lists from raw MPEs into np arrays
    
    Tracer_ms = np.array(Tracer)
    Natural_ms = np.array(Natural)
    Observed_ms = np.array(Observed)
    
    #Perform subtractions to account for background MS signal
    observed_natural = Observed_ms - Natural_ms
    tracer_natural = Tracer_ms - Natural_ms
    
    #Obtain slope of linear regression based on simple format of y=mx+b
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(tracer_natural, observed_natural)
    
    #slope is our "p", the fraction of the total flux into blood due to tracer
    
    #use the slope to calculate Turnover (Rate of Appearance, Ra)
    Ra = (infusion_rate/slope) - infusion_rate
    
    return Ra
