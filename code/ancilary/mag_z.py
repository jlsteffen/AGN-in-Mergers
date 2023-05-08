#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 15:14:18 2022

@author: joshua
"""

import sys, os
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import plotstyle as ps

drpall = fits.getdata('/Volumes/ssd2t/final_clean/0.drpall/drpall_clean_weight.fits')
z = drpall['nsa_z']
mag = drpall['NSA_ELPETRO_absmag'][:,5] + 5*np.log10(0.7)
mngtarg1 = drpall['mngtarg1']

p = ((mngtarg1 & 2**10)!=0)
s = ((mngtarg1 & 2**11)!=0)
c = ((mngtarg1 & 2**12)!=0)

fontsize = 20
fig, ax = plt.subplots(1,1, figsize=(6,6))
ax.scatter(z[c], mag[c], c='goldenrod', label='Color-Enhanced', s=5)
ax.scatter(z[p], mag[p], c='deepskyblue', label='Primary', s=5)
ax.scatter(z[s], mag[s], c='crimson', label='Secondary', s=5)

ax.set_xlim(-0.01, 0.16)
ax.set_ylim(-17, -25)
ax.set_xlabel('Redshift', fontsize=fontsize)
ax.set_ylabel('i-band Absolute Magnitude', fontsize=fontsize)

ps.legend(ax, fontsize=fontsize)

ps.ticks(ax, xmajor=0.05, ymajor=2, xminor=0.01, yminor=0.5)
ps.style(ax, fontsize=fontsize)

plt.savefig('mag_z.png', bbox_inches='tight')