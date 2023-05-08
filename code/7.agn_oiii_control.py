#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 13:29:18 2021

@author: joshua
"""
import os, sys
from astropy.io import fits
import numpy as np
import math as m
import matplotlib.pyplot as plt
import plotstyle as ps
from astropy.stats import bootstrap
#______________________________________________________________________________
def boot(data, cycles=1000, resample=0.68, func=np.mean):
    mask = ~np.isnan(data)
    if len(data[mask]) > 0:
        bootresult = bootstrap(data[mask], cycles, round(len(data[mask])*resample), bootfunc=func)
        med = np.nanmean(data)
        lower, upper = perc(bootresult)
        lower = med - lower
        upper = upper - med
    else:
        lower = 0
        upper = 0
    return lower, upper
#______________________________________________________________________________
def perc(dat, conflimit = 0.68):
    sdat = sorted(dat)
    lowindex = round(((1.0 - conflimit)/2) * len(dat))
    highindex = round(len(dat) - lowindex -  1)
    return sdat[lowindex], sdat[highindex]
#______________________________________________________________________________
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = '/0.data'

drpall = fits.getdata(datapath+'/0.drpall/drpall_clean_weight.fits')
pair = fits.getdata(filepath+'/1.pair.fits')

mr = drpall['NSA_ELPETRO_ABSMAG'][:,4]

mnsa = drpall['nsa_elpetro_mass']/0.7**2
z = drpall['nsa_z']

plateifu = drpall['plateifu']

dm = 10**pair['dlogm']
m1 = np.log10(mnsa/(1 + 1/dm))
m2 = np.log10(mnsa/(1 + dm))

m1 = np.log10(mnsa)

mnsa = np.log10(mnsa)
dm = pair['dlogm']

ctrl = (pair['ctrl']==1)&(pair['cagn']==1)
agnt = (pair['agnt']==1)
agnc = (pair['agnc']==1)
dagn = (pair['dagn']==1)  

bptt = pair['bptt']
bptc = pair['bptc']

lot = pair['loglo3t']
loc = pair['loglo3c']

lsun = m.log10(3.846*10**33)

lot = lot - lsun
loc = loc - lsun

lot[~np.isfinite(lot)|(lot<=4)|(lot>=9)] = np.nan
loc[~np.isfinite(loc)|(loc<=4)|(lot>=9)] = np.nan
#______________________________________________________________________________
zbs = 0.02
zbins = np.arange(0.02,0.15, zbs)

mbs = 0.2
mbins = np.arange(9.6, 11.4+mbs, mbs)
#______________________________________________________________________________
### The old 2-D plot ###
# Plot Parameters
fontsize = 18.0
fig, ax = plt.subplots(1,2, figsize=(12,6))

# Mass Panel
ax[0].scatter(mnsa[ctrl], lot[ctrl], c='lightgrey', label='Control Galaxy')
ax[0].scatter(m1[agnt], lot[agnt], c = 'deepskyblue', marker='D', label='Offset AGN')

ax[0].scatter(m1[agnt&dagn], lot[agnt&dagn], c = 'crimson', marker='D', label='Dual AGN')

# Z Panel
ax[1].scatter(z[ctrl], lot[ctrl], c='lightgrey', label='Control Galaxy')
ax[1].scatter(z[agnt], lot[agnt], c = 'deepskyblue', marker='D', label='Offset AGN')

ax[1].scatter(z[agnt&dagn], lot[agnt&dagn], c = 'crimson', marker='D', label='Dual AGN')

# Median luminosity in bin
mmed = []
mstd = []
for i in range(len(mbins)-1):
    ix = ctrl&(mnsa>mbins[i])&(mnsa<=mbins[i+1])
    mmed.append(np.log10(np.nanmedian(10**lot[ix])))
    mstd.append(boot(lot[ix]))
mmed = np.array(mmed)
mstd = np.array(mstd)

zmed = []
zstd = []
for i in range(len(zbins)-1):
    ix = ctrl&(z>zbins[i])&(z<=zbins[i+1])
    zmed.append(np.log10(np.nanmedian(10**lot[ix])))
    zstd.append(boot(lot[ix]))

zmed = np.array(zmed)
zstd = np.array(zstd)

ax[0].errorbar(mbins[1::]-mbs/2, mmed, xerr=mbs/2, yerr=mstd.T, c='k', \
         capsize=4, capthick=2, fmt='D')
ax[0].scatter(mbins[1::]-mbs/2, mmed, c='k', marker='D', label='Control Galaxy Median')
ax[1].errorbar(zbins[1::]-zbs/2, zmed, xerr=zbs/2, yerr=zstd.T, c='k', \
         capsize=4, capthick=2, fmt='D')

ax[0].set_xlim(9.4,11.6)
ax[0].set_ylim(4.8,9.2)
ax[1].set_xlim(-0.01, 0.16)
ax[1].set_ylim(4.8,9.2)

ax[0].set_xlabel(r'log(Stellar Mass/$M_\odot$)', fontsize=fontsize)
ax[0].set_ylabel(r'log(L[OIII] L/L$_{\odot}$)', fontsize=fontsize)
ax[1].set_xlabel('Redshift', fontsize=fontsize)
ax[1].set_ylabel(r'log(L[OIII] L/L$_{\odot}$)', fontsize=fontsize)

ps.legend(ax[0], fontsize=fontsize) 

ps.style(ax[0], fontsize=fontsize)
ps.style(ax[1], fontsize=fontsize, labelleft=False)
ps.ticks(ax[0], xmajor=0.5, ymajor=1, xminor=0.1, yminor=0.2)
ps.ticks(ax[1], xmajor=0.05, ymajor=1, xminor=0.01, yminor=0.2)

plt.savefig(filepath + '/22oiii_control.pdf', bbox_inches='tight')
plt.close('all')
#______________________________________________________________________________
lum=lot
mass=m1
mask=agnt
mlim=0.1
zlim=0.01
nctrl=20

ldiff = []
sep = []
d = []
dm2 =[]
m2 = []
z3 = []
for i in range(len(mass[mask])):
    mp = mass[mask][i]
    z2 = z[mask][i]
    
    # Find Pair's bins
    difm = np.abs(mnsa - mp) # |diff betw pair and bin centers|
    difz = np.abs(z - z2)
    
    ix = ctrl&(difm <= mlim)&(difz <= zlim)        
    if len(lum[ix]) >= nctrl:
        l = lum[ix]
        lmed = np.nanmedian(l)
        
        
        ldiff.append(lum[mask][i] - lmed)
        sep.append(pair['sepkpc'][mask][i])
        dm2.append(dm[mask][i])
        m2.append(mnsa[mask][i])
        z3.append(z[mask][i])
        
        if dagn[mask][i]==True:
            d.append(1)
        else:
            d.append(0)
ldiff = np.array(ldiff)
sep = np.array(sep)
dm2 = np.array(dm2)
m2 = np.array(m2)
z3 = np.array(z3)
#_______________________________________________________________________
# Plot enhancement
fontsize = 25

fig, ax = plt.subplots(1,1, figsize=(6,6))

ax.scatter(sep, ldiff, c='deepskyblue', marker='D', label='Offset AGN')
ax.scatter(sep[d==1], ldiff[d==1], c='crimson', marker='D', label='Dual AGN')
ax.plot([-1,26], [0,0], c='gray', linestyle='--')

rbs = 5
rbins = np.arange(0,25+rbs, rbs)

ldmed = np.zeros(rbins.shape) + np.nan
lstd = np.zeros((len(rbins), 2)) + np.nan
for i in range(len(rbins)-1):
    ix = (sep>rbins[i])&(sep<=rbins[i+1])
    ldmed[i] = np.nanmean(ldiff[ix])
    lstd[i] = boot(ldiff[ix])
    
ax.errorbar(rbins+rbs/2, ldmed, xerr=rbs/2, yerr=lstd.T, c='k', \
             capsize=4, capthick=2)
ax.scatter(rbins+rbs/2, ldmed, marker='D', c='k', label=r"Mean $\Delta$Log([OIII])")
    
ps.legend(ax, fontsize=fontsize)
ax.set_xlim(-1,26)
ax.set_ylim(-1.2,2.2)
ax.set_xlabel(r'r$_{\rm p}$ (kpc)', fontsize=fontsize)
ax.set_ylabel(r'$\Delta$log(L[OIII]) (dex)', fontsize=fontsize)

ps.style(ax, fontsize=fontsize)
ps.ticks(ax, xmajor=5, ymajor=1, xminor=1, yminor=0.2)
plt.savefig(filepath+'/7.agn_oiii_diff.pdf', bbox_inches='tight')
plt.close('all')
