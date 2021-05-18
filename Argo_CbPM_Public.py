#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 11:19:55 2020

@author: Lionel Alejandro Arteaga Quintero (lionel.arteagaquintero@nasa.gov, laaq@princeton.edu)

DISCLAIMER:
These code is provided as-is. There are no guarantees as to the presence of errors within the code. It is the user's responsibility to ensure that the code meets the user's needs.  However, please report any observed errors to the email listed above. 

CODE DESCRIPTION: Adapation of the Carbon-based Productivity Model (CbPM) by 
Westberry, et al. (2008) to estimate depth-resolved phytoplankton growth variables
forced by output from biogeochemical (bgc) profiling floats (adapted specifically for 
SOCCOM float data (https://soccom.princeton.edu)).
Code adapted from the satellite-based CbPM "updated" model at the OSU Ocean 
Productivity web site (http://sites.science.oregonstate.edu/ocean.productivity/cbpm2.code.php)

GENERAL MODEL (CbPM) DESCRIPTION: This is a spectrally resolved version of the cbpm, 
using nine separate wavelengths.  It is also depth resolved, integrating the effects 
from the surface down to a fixed depth of 200 m.

The CbPM algorithm estimates productivity using chl (m-1), bbp (m-1), surface 
irradiance (Einsteins m-2 d-1), k490 (m-1), mld (m), zno3 (m) and day length (hours).

Net primary productivity is phytoplankton carbon (Cphyto) \times growth rate, where 
carbon is  proportional to particulate backscatter (bbp). The orginal satellite-based 
CbPM converts satellite estimates of bbp to Cphyto. Here we used estimates
of Cphyto obtained from BGC-SOCCOM floats based on particulate 
organic carbon (POC) estimated from float measaurements of bbp.

    Cphyto (mg m^-3) = 0.19 * POC (mg m^-3) + 8.7 (Graff et al, 2015, Deep-Research I)
    (\mug l^-1) = (mg m^-3)

Phytoplankton growth rate is a function of nutrient and temperature stresst (f(nut,t))
and photoacclimation (f(Ig))

    growth rate (u) = umax * f(nut,T) * f(Ig)
    
where:    

    umax = 2

	f(nut,T) = ((Chl/C(z)) - (Chl/C)mu=0) / ((Chl/C)max - (Chl/C)mu=0)

	f(Ig) = 1 - exp (-5 * Ig)
    
and:

    (Chl/C(z)) = Ratio of float-basee chl and carbon at each depth (z)

	(Chl/C)max = 0.022 + (0.045-0.022) * exp (-3 * Ig) (max Chl/C given Ig)

    (Chl/C)mu=0 = (minimum) Chl:C at growth rate (mu) = 0    

The above items are analyzed for nine separate wavelengths, and is vertically 
resolved to a depth of 200 m.
    
For more details, please see the paper by Westberry, et al (2008).    
    
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
from AustinPetzold_1986 import AustinPetzold_1986
from daylength import daylength

def cbpm_bgcfloats(chl_z,Cphyto_z,irr,year,month,day,lat):
    # Spectral variables
    Lambda = [ 400, 412, 443, 490, 510, 555, 625, 670, 700];
    parFraction = [  0.0029, 0.0032, 0.0035, 0.0037, 0.0037, 0.0036, 0.0032, 0.0030, 0.0024];
    X = [ .11748, .122858, .107212, .07242, .05943, .03996, .04000, .05150, .03000];
    e = [ .64358, .653270, .673358, .68955, .68567, .64204, .64700, .69500, .60000];
    Kw = [ .01042, .007932, .009480, .01660, .03385, .06053, .28400, .43946, .62438];

    # Initialization variables for testing
    #year = 2015
    #month = 4
    #day = 24
    #lat = 30
    #chl_ml = 1 
    #irr = 30

    # Initialize necessary values for the model
    y0 = 0.0003 # min Chl:C ratio when mu = 0
    umax = 2.0 #after Banse (1991)
    chlC_z = np.zeros(200) # Phytoplankton Chl:C ratio (mg/mg)
    chlC_z[:] = np.nan
    nutTempFunc_z = np.zeros(200) # Nutrient dependent term
    nutTempFunc_z[:] = np.nan
    IgFunc_z = np.zeros(200) # Light dependent term
    IgFunc_z[:] = np.nan
    mu_z = np.zeros(200) # Phytoplankton growth rate
    mu_z[:] = np.nan
    prcnt_z = np.zeros(200) # Remaining light fraction
    prcnt_z[:] = np.nan
    pp_z = np.zeros(200) # Primary production
    pp_z[:] = np.nan
    Ezlambda = np.zeros([200,9]) # Fraction of light at nine wavelengths
    Ezlambda[:] = np.nan
    par_z = np.zeros(200) # Photosynthetically available radiation 
    par_z[:] = np.nan

    # Attenuation coefficient at 490nm based on Equation 8 of 
    # Morel, et al. (2007, Remote Sens. Environ.)
    chl_surf = chl_z[0]
    k490 = 0.0166 + 0.0773 * pow(chl_surf,0.6715);

    # Calculate length of day in hours
    Daylength = daylength(year,month,day,lat)

    # Multispectral comonent of the model
    klambda = np.zeros(len(Lambda))
    klambda[:] = np.nan
    E0 = np.zeros(len(Lambda))
    E0[:] = np.nan
    kbio = np.zeros(len(Lambda))
    kbio[:] = np.nan
    kd = np.zeros(len(Lambda))
    kd[:] = np.nan
    kdif = np.zeros(len(Lambda))
    kdif[:] = np.nan

    for n in range(len(Lambda)):
        klambda[n] = AustinPetzold_1986(Lambda[n],k490)
        E0[n] = irr * parFraction[n]
        #Ez_mld(i) = Eo(i) .* 0.975 .* exp(-klambda(i) .* mld/2.0);
    
    # Calculate Kd offset carry through to depth non-chl attenuation
    #for n in range(len(Lambda)):
    #   kbio = X[n] * pow(chl_ml,e[n])
    #  kd[n] = Kw[n] + kbio
    # kdif[n] = klambda[n] - kd[n]

    for z in range(len(chl_z)):
        if z == 0:        
            for n in range(len(Lambda)):
                Ezlambda[z,n] = E0[n] * 0.975 * np.exp(-klambda[n] * z)
        else:
            for n in range(len(Lambda)):
                kbio = X[n] * pow(chl_z[z-1],e[n]);   # after Morel and Maritorena (2001)
                kd[n] = Kw[n] + kbio
                #kdif[n] = klambda[n] - kd[n]
                #kd[n] = kdif[n] + kd[n];
                Ezlambda[z,n] = Ezlambda[z-1,n] * np.exp(-kd[n] * 1)
        #print(z)        
        par_z[z] = 0.0;
        for n in range(len(Lambda)-1):
            par_z[z] = ((Lambda[n+1] - Lambda[n]) * (Ezlambda[z,n+1] + Ezlambda[z,n]) /2) + par_z[z]
        #print(par_z[z])    
        
        chlC_z[z] = chl_z[z] / Cphyto_z[z]
        chlCarbonMax_z = 0.022 + (0.045-0.022) * np.exp(-3.0 * par_z[z] / Daylength)
        nutTempFunc_z[z] = (chlC_z[z] - y0) / (chlCarbonMax_z - y0)
        if nutTempFunc_z[z] > 1: 
            nutTempFunc_z[z] = 1
        IgFunc_z[z] = 1 - np.exp(-5.0 * par_z[z] / Daylength) 
        mu_z[z] = umax * nutTempFunc_z[z] * IgFunc_z[z]
        if mu_z[z] > umax:
            mu_z[z] = umax
        prcnt_z[z] = par_z[z] / irr * 0.975

        # Track 1% of surf. irradiance
        if prcnt_z[z] >= 0.01:
            mzeu = z
        
        # Depth-resolved primary production
        pp_z[z] = mu_z[z] * Cphyto_z[z]
        
    return pp_z,mu_z,par_z,prcnt_z   
       