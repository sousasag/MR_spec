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

def bc_flower(teff):
  """
  Get the bolometric correction using the relations from flower 1996
  """
  lteff=np.log10(teff)
  if (lteff<3.7):
    bcflow=-0.190537291496456e+05+0.155144866764412e+05*lteff-0.421278819301717e+04*(lteff*lteff)+0.381476328422343e+03*(lteff*lteff*lteff)
  if (lteff>=3.7 and lteff<3.9):
    bcflow=-0.370510203809015e+05+0.385672629965804e+05*lteff-0.150651486316025e+05*(lteff*lteff)+0.261724637119416e+04*(lteff*lteff*lteff)-0.170623810323864e+03*(lteff*lteff*lteff*lteff)
  if (lteff>=3.9):
    bcflow=-0.370510203809015e+05+0.385672629965804e+05*lteff-0.150651486316025e+05*(lteff*lteff)+0.261724637119416e+04*(lteff*lteff*lteff)-0.170623810323864e+03*(lteff*lteff*lteff*lteff)
  return bcflow

def get_logg_hip(mass, teff, paralax, vmag):
  bcflow = bc_flower(teff)
  logg_hip = np.log10(mass) + 4.*np.log10(teff/5777.) + 2.*np.log10(paralax/1000.) + 0.4 * (vmag + bcflow) + 0.11 + 4.44
  return logg_hip

def logg_mass_hip_interactive(teffs,loggs,fehs,vmag, paralax):
    logg_ga = -100000
    niter = 0
    logg_g = loggs
    while(abs(logg_g-logg_ga) > 0.01 and niter < 10):
        MT, Mcor, logM  = MR.mass_torres2010(teffs,logg_g,fehs)
        logg_ga = logg_g
        logg_g = get_logg_hip(Mcor, teffs, paralax, vmag)
        #print(niter, Mcor, logg_ga, logg_g)
        niter +=1
    if niter >= 10:
        print("Did not converge!")
    #print(logg_g, Mcor)
    return logg_g, Mcor

def logg_mass_hip_iteractive_error(teffs, eteffs, loggs, eloggs, fehs, efehs, vmag, evmag, paralax, eparalax, npoints=10000):
    teffss = np.random.normal(teffs, eteffs, npoints)
    loggss = np.random.normal(loggs, eloggs, npoints)
    fehss = np.random.normal(fehs, efehs, npoints)
    vmags = np.random.normal(vmag, evmag, npoints)
    paralaxs = np.random.normal(paralax, eparalax, npoints)
    logg_dist = np.zeros(npoints)
    mass_dist = np.zeros(npoints)
    for i in range(npoints):
#        print(i, teffss[i], loggss[i], fehss[i], gaia_gmags[i], gaia_paralaxs[i])
        logg_dist[i], mass_dist[i] = logg_mass_iteractive(teffss[i], loggss[i], fehss[i], vmags[i], paralaxs[i])

    meanlogg = np.nanmean(logg_dist)
    stdlogg = np.nanstd(logg_dist)
    meanmass = np.nanmean(mass_dist)
    stdmass = np.nanstd(mass_dist)

#    plt.figure(1)
#    plt.hist(logg_gaia_dist)
#    plt.figure(2)
#    plt.hist(mass_gaia_dist)
#    plt.show()

    return meanlogg, stdlogg, meanmass, stdmass





def logg_gaia(mass, teff, gaia_gmag, gaia_paralax, Ag = 0, teff_sun=5777, distance_gaia = 0):

    if distance_gaia == 0:
      MG = gaia_gmag + 5 - 5*log10(1000./gaia_paralax) - Ag
    else:
      MG = gaia_gmag + 5 - 5*log10(distance_gaia) - Ag  
    bcg,sbcg = BC_g_Andrae(teff)
    mbol_star = MG + bcg
    mbol_sun = 4.74  #IAU Resolution 2015
    logg_sun =  logg(1,1)
    logg_gaia = log10(mass) + 4*log10(teff)-4*log10(teff_sun) + 0.4*mbol_star - 0.4*mbol_sun + logg_sun
    return logg_gaia

def logg_gaia_error(mass, emass, teff, eteff, gaia_gmag, egaia_gmag, gaia_paralax, egaia_paralax, Ag=0, teff_sun=5777, distance_gaia = 0, e_distance_gaia=0, npoints = 10000):
    masss = np.random.normal(mass, emass, npoints)
    teffs = np.random.normal(teff, eteff, npoints)
    gaia_gmags = np.random.normal(gaia_gmag, egaia_gmag, npoints)
    gaia_paralaxs = np.random.normal(gaia_paralax, egaia_paralax, npoints)
    gaia_distance = np.random.normal(distance_gaia, e_distance_gaia, npoints)

    logg_gaia_dist = np.zeros(npoints)

    for i in range(npoints):
        logg_gaia_dist[i] = logg_gaia(masss[i], teffs[i], gaia_gmags[i], gaia_paralaxs[i], distance_gaia = gaia_distance[i])

    meanlogg = np.nanmean(logg_gaia_dist)
    stdlogg = np.nanstd(logg_gaia_dist)
    return meanlogg, stdlogg


def logg_mass_iteractive(teffs, loggs, fehs, gaia_gmag, gaia_paralax, Ag=0, teff_sun=5777, distance_gaia = 0, calib='torres'):
    logg_ga = -100000
    niter = 0
    logg_g = loggs
    while(abs(logg_g-logg_ga) > 0.01 and niter < 10):
        MT, Mcor, logM  = MR.mass_torres2010(teffs,logg_g,fehs, calib=calib)
        logg_ga = logg_g
        logg_g = logg_gaia(Mcor, teffs, gaia_gmag, gaia_paralax, Ag=Ag, teff_sun=teff_sun, distance_gaia=distance_gaia)
        #print(niter, Mcor, logg_ga, logg_g)
        niter +=1
    if niter >= 10:
        print("Did not converge!")
    #print(logg_g, Mcor)
    return logg_g,Mcor


def logg_mass_iteractive_error(teffs, eteffs, loggs, eloggs, fehs, efehs, gaia_gmag, egaia_gmag, gaia_paralax, egaia_paralax, Ag=0, teff_sun=5777, distance_gaia = 0, e_distance_gaia=0, npoints=10000, calib='torres'):
    teffss = np.random.normal(teffs, eteffs, npoints)
    loggss = np.random.normal(loggs, eloggs, npoints)
    fehss = np.random.normal(fehs, efehs, npoints)
    gaia_gmags = np.random.normal(gaia_gmag, egaia_gmag, npoints)
    gaia_paralaxs = np.random.normal(gaia_paralax, egaia_paralax, npoints)
    gaia_distance = np.random.normal(distance_gaia, e_distance_gaia, npoints)
    logg_gaia_dist = np.zeros(npoints)
    mass_gaia_dist = np.zeros(npoints)
    for i in range(npoints):
#        print(i, teffss[i], loggss[i], fehss[i], gaia_gmags[i], gaia_paralaxs[i])
        logg_gaia_dist[i], mass_gaia_dist[i] = logg_mass_iteractive(teffss[i], loggss[i], fehss[i], gaia_gmags[i], gaia_paralaxs[i], distance_gaia=gaia_distance[i], calib=calib)

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
       
def get_logg_mass_radius_gaia_torres(teffs, eteffs, loggs, eloggs, fehs, efehs, gaia_gmag, egaia_gmag, gaia_paralax, egaia_paralax, Ag=0, teff_sun=5777, distance_gaia = 0, e_distance_gaia=0, npoints=10000, calib='torres'):
    meanlogg, stdlogg, meanmass, stdmass = logg_mass_iteractive_error(teffs, eteffs, loggs, eloggs, fehs, efehs, gaia_gmag, egaia_gmag, gaia_paralax, egaia_paralax, Ag=Ag, teff_sun=teff_sun, distance_gaia=distance_gaia, e_distance_gaia=e_distance_gaia, npoints=npoints, calib=calib)
    r, er = MR.radius_torres2010_error(teffs, eteffs, meanlogg, stdlogg, fehs, efehs, npoints = npoints, calib=calib)
    return meanlogg, stdlogg, meanmass, stdmass, r, er


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
    teff = 5433; eteff = 62; loggs = 4.404; eloggs = 0.111; feh = -0.548; efeh=0.042;
    par = 5.3891; epar = 0.0263; gmag = 11.5837; egmag = 0.0028

#    print(logg_mass_iteractive_error(4684, 99, 4.16, 0.24, -0.203, 0.05, 9.389915, 0.002759, 25.5817, 0.0127))
#    print(MR.radius_torres2010_error(4684, 99, 4.533, 0.056, -0.203, 0.05))
    teff = 4773; eteff= 94; loggs= 4.59; eloggs= 1.46; feh = -0.22; efeh =  0.6; 
    par = 2.2893; epar =  0.0271; gmag = 14.5759; egmag = 0.0028


    teff = 4552; eteff= 154; loggs= 4.198; eloggs= 0.576; feh = -0.28; efeh =  0.067; 
    par = 36.3525; epar =  0.0161; gmag = 9.407; egmag = 0.00277

    loggh, eloggh, massh, emassh = logg_mass_iteractive_error(teff, eteff, loggs, eloggs, feh, efeh, gmag, egmag, par, epar)
    rh, erh = MR.radius_torres2010_error(teff, eteff, loggh, eloggh, feh, efeh)
    print(loggh, eloggh, massh, emassh, rh, erh)
    loggh, eloggh, massh, emassh = logg_mass_iteractive_error(teff, eteff, loggs, eloggs, feh, efeh, gmag, egmag, par, epar, calib='maxted')
    rh, erh = MR.radius_torres2010_error(teff, eteff, loggh, eloggh, feh, efeh, calib='maxted')
    print(loggh, eloggh, massh, emassh, rh, erh)


if __name__ == '__main__':
    main()
