#!/usr/bin/python

import argparse
from astropy import constants as c
from numpy import log10
import numpy as np
import os
import MR_spec as MR
from matplotlib import pyplot as plt

def _parse():
    '''Calculate the surface gravity from M and R'''
    p = argparse.ArgumentParser(description='Calculate logg from solar M and R')
    p.add_argument('M', help='Mass in solar units', type=float)
    p.add_argument('R', help='Radius in solar units', type=float)
    return p.parse_args()


def logg(M, R):
    """Mass and radius in solar units"""
    G = c.G.value * 1e3
    M *= c.M_sun.value * 1e3
    R *= c.R_sun.value * 1e2
    return log10(G*M/(R**2))

def BC_g_Andrae(Teff, Tsun = 5777):
    a  = [6e-2    , 6.731e-5, -6.647e-8,  2.859e-11, -7.197e-15]
    sa = [2.634e-2, 2.438e-5, -1.129e-9, -6.722e-12,  1.635e-15]
    if Teff < 4000:
        a  = [ 1.749,  1.977e-3, 3.737e-7, -8.966e-11, -4.183e-14]
        sa = [-2.487, -1.8762-3, 2.128e-7,  3.807e-10,  6-570e-14]

    bcg  = np.sum([a[i]*(Teff-Tsun)**i for i in range(len(a))])
    sbcg = np.sum([sa[i]*(Teff-Tsun)**i for i in range(len(sa))])
    return bcg, sbcg

def logg_gaia(mass, teff, gaia_gmag, gaia_paralax, Ag = 0, teff_sun=5777):

    MG = gaia_gmag + 5 - 5*log10(1000./gaia_paralax) - Ag
    bcg,sbcg = BC_g_Andrae(teff)
    mbol_star = MG + bcg
    mbol_sun = 4.74  #IAU Resolution 2015
    logg_sun =  logg(1,1)
    logg_gaia = log10(mass) + 4*log10(teff)-4*log10(teff_sun) + 0.4*mbol_star - 0.4*mbol_sun + logg_sun
    return logg_gaia

def logg_gaia_error(mass, emass, teff, eteff, gaia_gmag, egaia_gmag, gaia_paralax, egaia_paralax, Ag=0, teff_sun=5777, npoints = 10000):
    masss = np.random.normal(mass, emass, npoints)
    teffs = np.random.normal(teff, eteff, npoints)
    gaia_gmags = np.random.normal(gaia_gmag, egaia_gmag, npoints)
    gaia_paralaxs = np.random.normal(gaia_paralax, egaia_paralax, npoints)

    logg_gaia_dist = np.zeros(npoints)

    for i in range(npoints):
        logg_gaia_dist[i] = logg_gaia(masss[i], teffs[i], gaia_gmags[i], gaia_paralaxs[i])

    meanlogg = np.nanmean(logg_gaia_dist)
    stdlogg = np.nanstd(logg_gaia_dist)
    return meanlogg, stdlogg


def logg_mass_iteractive(teffs, loggs, fehs, gaia_gmag, gaia_paralax, Ag=0, teff_sun=5777):
    logg_ga = 0
    niter = 0
    logg_g = loggs
    while(abs(logg_g-logg_ga) > 0.01 and niter < 10):
        MT, Mcor, logM  = MR.mass_torres2010(teffs,logg_g,fehs)
        logg_ga = logg_g
        logg_g = logg_gaia(Mcor, teffs, gaia_gmag, gaia_paralax, Ag=Ag, teff_sun=teff_sun)
        #print(niter, Mcor, logg_ga, logg_g)
        niter +=1
    #print(logg_g, Mcor)
    return logg_g,Mcor


def logg_mass_iteractive_error(teffs, eteffs, loggs, eloggs, fehs, efehs, gaia_gmag, egaia_gmag, gaia_paralax, egaia_paralax, Ag=0, teff_sun=5777, npoints=10000):
    teffss = np.random.normal(teffs, eteffs, npoints)
    loggss = np.random.normal(loggs, eloggs, npoints)
    fehss = np.random.normal(fehs, efehs, npoints)
    gaia_gmags = np.random.normal(gaia_gmag, egaia_gmag, npoints)
    gaia_paralaxs = np.random.normal(gaia_paralax, egaia_paralax, npoints)
    logg_gaia_dist = np.zeros(npoints)
    mass_gaia_dist = np.zeros(npoints)
    for i in range(npoints):
        logg_gaia_dist[i], mass_gaia_dist[i] = logg_mass_iteractive(teffss[i], loggss[i], fehss[i], gaia_gmags[i], gaia_paralaxs[i])

    meanlogg = np.nanmean(logg_gaia_dist)
    stdlogg = np.nanstd(logg_gaia_dist)
    meanmass = np.nanmean(mass_gaia_dist)
    stdmass = np.nanstd(mass_gaia_dist)

#    plt.figure(1)
#    plt.hist(logg_gaia_dist)
#    plt.figure(2)
#    plt.hist(mass_gaia_dist)
#    plt.show()

    return meanlogg, stdlogg, meanmass, stdmass
       





def main():
    dir_script = os.path.dirname(__file__)
    if dir_script != "" :
        dir_script+='/'
    print(dir_script)
    args = _parse()
    print(('logg: %.2f' % logg(args.M, args.R)))

#    logg_mass_iteractive(5810, 21, 4.33, 0.03, 0.21, 0.02, 5.2832, 0.0028, 64.4048, 0.0771)
#    logg_mass_iteractive(6330, 130, 4.37, 0.2, 0.31, 0.15, 12.6441, 0.0028, 1.6613, 0.0126)
#    print(logg_mass_iteractive_error(6382, 22, 4.60, 0.04, 0.07, 0.02, 7.6937, 0.0028, 17.8725, 0.0206))



if __name__ == '__main__':
    main()
