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
  return (MT,Mcor)


def radius_torres2010_error(teff, erteff, logg, erlogg, feh, erfeh, npoints = 10000):
  """
  Get the radius of a star using the Torres calibrations from 2010
  Compute the errors as in Santos et al. 2013
  """
  teffs = np.random.normal(teff, erteff, npoints)
  loggs = np.random.normal(logg, erlogg, npoints)
  fehs  = np.random.normal(feh, erfeh, npoints)

  radius_dist = np.zeros(npoints)
  for i in range(npoints):
    RT = radius_torres2010(teffs[i],loggs[i],fehs[i])
    radius_dist[i] = RT
  return np.mean(radius_dist), np.std(radius_dist)



def mass_torres2010_error(teff, erteff, logg, erlogg, feh, erfeh, npoints = 10000):
  """
  Get the mass of a star using the Torres calibrations from 2010
  Compute the errors as in Santos et al. 2013
  """
  teffs = np.random.normal(teff, erteff, npoints)
  loggs = np.random.normal(logg, erlogg, npoints)
  fehs  = np.random.normal(feh, erfeh, npoints)

  mass_dist = np.zeros(npoints)
  for i in range(npoints):
    MT, Mcor = mass_torres2010(teffs[i],loggs[i],fehs[i])
    mass_dist[i] = Mcor
  MT, Mcor = mass_torres2010(teff,logg,feh)
  return MT, Mcor, np.mean(mass_dist), np.std(mass_dist), np.sqrt(np.std(mass_dist)**2. + 0.027**2.)



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
