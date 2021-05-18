#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 12:08:23 2020

@author: Lionel Alejandro Arteaga Quintero (laaq) (laaq@princeton.edu, lionel.arteagaquintero@nasa.gov, artionel@gmail.com)

Given a reference k490 vlaue, determine k(lambda) for a specified lambda.
Taken from the CbPM code of Westberry et al (2008) as 
in http://www.science.oregonstate.edu/ocean.productivity/cbpm2.code.php
Based on Austin and Petzold (1986).
"""
def AustinPetzold_1986(Lambda,k490):
    wave = [350, 360, 370, 380, 390, 400, 410, 420\
            , 430, 440, 450, 460, 470, 480, 490, 500, 510\
            , 520, 530, 540, 550, 560, 570, 580, 590, 600\
            , 610, 620, 630, 640, 650, 660, 670, 680, 690, 700];

    M = [2.1442, 2.0504, 1.9610, 1.8772, 1.8009, 1.7383, 1.7591,\
         1.6974, 1.6108, 1.5169, 1.4158, 1.3077, 1.1982, 1.0955,\
         1.0000, 0.9118, 0.8310, 0.7578, 0.6924, 0.6350, 0.5860,\
         0.5457, 0.5146, 0.4935, 0.4840, 0.4903, 0.5090, 0.5380,\
         0.6231, 0.7001, 0.7300, 0.7301, 0.7008, 0.6245, 0.4901, 0.2891];

    Kdw = [0.0510, 0.0405, 0.0331, 0.0278, 0.0242, 0.0217, 0.0200,\
           0.0189, 0.0182, 0.0178, 0.0176, 0.0176, 0.0179, 0.0193, 0.0224,\
           0.0280, 0.0369, 0.0498, 0.0526, 0.0577, 0.0640, 0.0723, 0.0842,\
           0.1065, 0.1578, 0.2409, 0.2892, 0.3124, 0.3296, 0.3290, 0.3559,\
           0.4105, 0.4278, 0.4521, 0.5116, 0.6514];

    # Interpolate to wavelength of interest
    for i in range(len(wave)):
        if wave[i] >= Lambda:
            l1 = wave[i]
            k1 = Kdw[i]
            m1 = M[i]
            l0 = wave[i-1]
            k0 = Kdw[i-1]
            m0 = M[i-1]
            break

    num = Lambda - l0
    den = l1 - l0
    frac = num / den
    
    kdiff = k1 - k0
    Kdw_l = k0 + frac * kdiff

    mdiff = m1 - m0
    M_l = m0 + frac*mdiff

    # Get reference wavelength (490 for now) and appply model
    ref = 14
    Kd = (M_l/M[ref]) * (k490 - Kdw[ref]) + Kdw_l
    return Kd





