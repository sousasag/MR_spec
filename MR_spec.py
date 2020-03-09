#!/usr/bin/python
## My first python code

##imports:

import numpy as np

## My functions:

def massTorres_online(teff, erteff, logg, erlogg, feh, erfeh):
    """Calculate a mass using the Torres calibration"""
    ntrials = 10000
    randomteff = teff + erteff * np.random.randn(ntrials)
    randomlogg = logg + erlogg * np.random.randn(ntrials)
    randomfeh = feh + erfeh * np.random.randn(ntrials)

    # Parameters for the Torres calibration:
    a1 = 1.5689
    a2 = 1.3787
    a3 = 0.4243
    a4 = 1.139
    a5 = -0.1425
    a6 = 0.01969
    a7 = 0.1010

    M = np.zeros(ntrials)
    logM = np.zeros(ntrials)
    for i in xrange(ntrials):
        X = np.log10(randomteff[i]) - 4.1
        logMass = a1 + a2 * X + a3 * X * X + a4 * X * X * X + a5 *\
            randomlogg[i] * randomlogg[i] + a6 * randomlogg[i] *\
            randomlogg[i] * randomlogg[i] + a7 * randomfeh[i]
        logM[i] = logMass
        M[i] = 10 ** logMass

    meanMasslog = np.mean(logM)
    sigMasslog = np.sqrt(np.sum([(logMi - meanMasslog)**2 for logMi in logM]) /
                         (ntrials - 1))
    sigMasslogTot = np.sqrt(0.027*0.027 + sigMasslog*sigMasslog)

    meanMass = 10**meanMasslog
    sigMass = 10**(meanMasslog + sigMasslogTot) - meanMass

    return meanMass, sigMass


def mass_torres2010(teff,logg,feh):
  """
  Get the mass of a star using the Torres calibrations from 2010
  """
  ai=[1.5689,1.3787,0.4243,1.139,-0.1425,0.01969,0.1010]
#  eai=[0.058,0.029,0.029,0.24,0.011,0.0019,0.014]
  X=np.log10(teff)-4.1
  logM=ai[0]+ai[1]*X+ai[2]*X**2.+ai[3]*X**3+ai[4]*logg**2+ai[5]*logg**3+ai[6]*feh
  MT=10**logM


  #limits of correction
  if MT >= 0.7 and MT <= 1.3:
    Mcor = 0.791*MT**2.-0.575*MT+0.701 # Santos et al 2013 correction
  else:
    Mcor = MT
  return (MT,Mcor, logM)



def radius_torres2010(teff,logg,feh):
  """
  Get the radius of a star using the Torres calibrations from 2010
  """
  bi=[2.4427,0.6679,0.1771,0.705,-0.21415, 0.02306,0.04173]
#  ebi=[0.038,0.016,0.027,0.13,0.0075,0.0013,0.0082]
  X=np.log10(teff)-4.1
  logR=bi[0]+bi[1]*X+bi[2]*X**2.+bi[3]*X**3+bi[4]*logg**2+bi[5]*logg**3+bi[6]*feh
  RT=10**logR
  return (RT, logR)


def radius_torres2010_error(teff, erteff, logg, erlogg, feh, erfeh, npoints = 10000):
  """
  Get the radius of a star using the Torres calibrations from 2010
  Compute the errors as in Santos et al. 2013
  """
  teffs = np.random.normal(teff, erteff, npoints)
  loggs = np.random.normal(logg, erlogg, npoints)
  fehs  = np.random.normal(feh, erfeh, npoints)

  radius_dist = np.zeros(npoints)
  logradius_dist = np.zeros(npoints)

  for i in range(npoints):
    RT, logR = radius_torres2010(teffs[i],loggs[i],fehs[i])
    radius_dist[i] = RT
    logradius_dist = logR

  meanlogR = np.mean(logradius_dist)
  stdlogR = np.std(logradius_dist)
  R = 10**meanlogR
  erlogR = np.sqrt(stdlogR**2. + 0.014**2.)
  erR =  10**(meanlogR + erlogR) - R
  return R, erR



def mass_torres2010_error(teff, erteff, logg, erlogg, feh, erfeh, npoints = 10000):
  """
  Get the mass of a star using the Torres calibrations from 2010
  Compute the errors as in Santos et al. 2013
  """
  teffs = np.random.normal(teff, erteff, npoints)
  loggs = np.random.normal(logg, erlogg, npoints)
  fehs  = np.random.normal(feh, erfeh, npoints)

  mass_dist = np.zeros(npoints)
  logmass_dist = np.zeros(npoints)
  for i in range(npoints):
    MT, Mcor, logM  = mass_torres2010(teffs[i],loggs[i],fehs[i])
    mass_dist[i] = Mcor
    logmass_dist[i] = logM
  meanlogM = np.mean(logmass_dist)
  stdlogM = np.std(logmass_dist)
  errorlogM = np.sqrt(stdlogM**2. + 0.027**2.)
  MT = 10**meanlogM
  #limits of correction
  if MT >= 0.7 and MT <= 1.3:
    Mcor = 0.791*MT**2.-0.575*MT+0.701 # Santos et al 2013 correction
  else:
    Mcor = MT
  errorM = 10**(meanlogM + errorlogM) - MT

  return MT, Mcor, errorM



def test_functions():
    #Solar
    teff = 5777
    er_teff = 50
    logg = 4.4
    er_logg = 0.1
    # vtur = 1.0
    # er_vtur = 0.1
    feh = 0.0
    er_feh = 0.05
    print (massTorres_online(teff, er_teff, logg, er_logg, feh, er_feh))
    print mass_torres2010_error(teff, er_teff, logg, er_logg, feh, er_feh)


### Main program:
def main():
  print "Hello"
  test_functions()


if __name__ == "__main__":
    main()
