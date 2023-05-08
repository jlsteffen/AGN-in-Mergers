#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 16:29:05 2021

@author: joshua
"""

import os, sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import plotstyle as ps
import math as m

def lumin(ma, alpha=-1.26, mstar=-20.71 + 5*m.log10(0.7), phi=0.0093*(0.7**3)*10**6):
    return np.log10(0.4*m.log(10)*phi*(10**(-0.4*(ma-mstar)))**(alpha + 1) * np.exp(-10**(-0.4*(ma-mstar))))

drpall = fits.getdata('/Volumes/Ext_drive/Research/final/0.drpall/drpall_main_weight.fits')

ew = drpall['esweight']
pw = drpall['eweight']
sw = drpall['sweight']

mag = drpall['NSA_ELPETRO_absmag'][:,4] + 5*m.log10(0.7)
mngtarg1 = drpall['mngtarg1']

# samples
pp = ((mngtarg1 & 2**10+2**12)!=0)
s = ((mngtarg1 & 2**11)!=0)
#plot all


# bin sizes
bs = 0.5
bins = np.arange(-24, -16 + bs, bs)

full = []
prim = []
sec  = []
for i in range(len(bins)-1):
    ix = (mag>=bins[i])&(mag<bins[i+1])
    if len(ew[ix]) > 0:
        full.append(ew[ix].sum())
    else:
        full.append(np.nan)
    
    if len(pw[ix&pp]) > 0:
        prim.append(pw[ix&pp].sum())
    else:
        prim.append(np.nan)
        
    if len(sw[ix&s]) > 0:
        sec.append(sw[ix&s].sum())
    else:
        sec.append(np.nan)

full = np.array(full)
prim = np.array(prim)
sec  = np.array(sec )
    
# plotting
fontsize = 18
fig, ax = plt.subplots(1,1, figsize=(6,6))
ax.scatter(mag, np.log10(ew), c='k', s=.1)
ax.scatter(mag[pp], np.log10(pw[pp]), c='deepskyblue', s=.1)
ax.scatter(mag[s], np.log10(sw[s]), c='crimson', s=.1)

ax.scatter(bins[0:-1], np.log10(full), c='k', marker='D', label='Main Sample', s=70)
ax.scatter(bins[0:-1], np.log10(prim), c='deepskyblue', marker='o', label='Primary+')
ax.scatter(bins[0:-1], np.log10(sec), c='crimson', marker='v', label='Secondary')

logn = lumin(bins)
ax.plot(bins, logn, 'k--', label='Montero-Dorta+09')
ax.fill_between(bins, lumin(bins, alpha=-1.26+0.05, mstar=-20.71 + 5*m.log10(0.7)+0.04, phi=0.0093*(0.7**3)*10**6 +0.07), \
                  lumin(bins, alpha=-1.26-0.05, mstar=-20.71 + 5*m.log10(0.7)-0.04, phi=0.0093*(0.7**3)*10**6 -0.07), \
                      facecolor='grey', edgecolor='k')

ax.set_xlabel('Absolute Magnitude (r-band)', fontsize=fontsize)
ax.set_ylabel(r'log(n) # / 10$^6$ MPC$^3$ / mag', fontsize=fontsize)

ps.legend(ax, fontsize=fontsize)

ax.set_xlim(-16.5, -24.5)

ps.style(ax, fontsize=fontsize)
ps.ticks(ax, xmajor=1, xminor=0.2, ymajor=1, yminor=0.2)

plt.savefig('lumfunc.pdf', bbox_inches='tight')
plt.close()