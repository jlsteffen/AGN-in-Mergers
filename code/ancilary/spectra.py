#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 16:19:30 2022

@author: joshua
"""

import os, sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import plotstyle as ps
from PIL import Image
from astropy import wcs
import math as m
from copy import copy
#______________________________________________________________________________
plateifu = '7963-6102'

size = 300
scale = 40.0/size
def png(im, ax, plateifu):
    ind = np.where(drpall['plateifu']==plateifu)[0][0]
    IFU_size = int(drpall['ifudsgn'][ind][0:-2])
    IFU_arcsec = np.array(([19, 37, 61, 91, 127], [16.0, 21.0, 26.0, 31.0, 36.0]))
    ifu_scale = IFU_arcsec[1, np.where(IFU_size==IFU_arcsec[0])[0][0]]
    cent_pix = np.array([im.size[0]/2.0, im.size[1]/2.0])
    # Make a Astropy WCS
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [cent_pix[0], cent_pix[1]]
    w.wcs.cdelt = np.array([-scale/3600.0, -scale/3600.0])
    w.wcs.crval = [drpall['ifura'][ind], drpall['ifudec'][ind]]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ['deg', 'deg']
    w.wcs.set_pv([(2, 1, 45.0)])
    ax.imshow(im, aspect='auto')
    ax.set_title(plateifu, fontsize=fontsize)
    # IFU hexagons    
    hexagon = patches.RegularPolygon([size/2.0, size/2.0], 6,linestyle='--',
              radius = size/2.0 * ifu_scale/40.0, ec = 'lime', fc = 'none', orientation=m.pi/2)
    ax.add_patch(copy(hexagon))
    ax.text(0.1, 0.9, '(a)', fontsize=fontsize, color='w', ha='center', 
            va='center', transform=ax.transAxes)
    # Label the dAGN
    markers = ['o', 's', 'H']
    for j in range(len(data)):
        # Convert object ra and dec to by in pixel coordinates
        crd = np.array([[drpall['objra'][ind], drpall['objdec'][ind]]], np.float_)        
        pix = w.wcs_world2pix(crd, 1)
        obj_x = pix[0][0]
        obj_y = pix[0][1]
        ax.scatter(obj_x, obj_y, marker = markers[j], s=250, facecolors='none', edgecolors='lime')
        #ax.text(obj_x-5, obj_y - 17, str(j), fontsize=fontsize/1.5, color='w')
     
    loc1 = 75
    loc2 = 15
    ps.style(ax, fontsize=fontsize, labelbottom=False, labelleft=False, c='w')
    ps.ticks(ax, xmajor=loc1, ymajor=loc1, xminor=loc2, yminor=loc2)
#______________________________________________________________________________
# Paths and Files
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = os.path.abspath(os.path.join(filepath, '..', 'data'))

drpall = fits.getdata('/Volumes/ssd2t/final_clean/0.drpall/drpall_clean_weight.fits')

spec = fits.open('/Volumes/ssd2t/final_clean/1.spfit/r1arc/'+plateifu+'_spec.fits')
data = spec[1].data
#______________________________________________________________________________
# Initialize Plots
fontsize = 14.0
fig = plt.figure(figsize= (10,9))
#gs = plt.GridSpec(30, 16, width_ratios = [1]*16, height_ratios = [1]*30)
gs = plt.GridSpec(9, 10, width_ratios = [1]*10, height_ratios = [1]*9)
#______________________________________________________________________________
# Galaxy png
ax1 = plt.subplot(gs[0:4, 0:4])
im = Image.open(datapath +'/png/'+str(plateifu)+'.png')
png(im, ax1, plateifu)
#______________________________________________________________________________
# Full Spectra
ax6 = plt.subplot(gs[5:8, 0:10])
ax6.plot(10**data['log10lam'][0], data['galaxy'][0], c='k', linewidth=1, label='Observed')
ax6.plot(10**data['log10lam'][0], data['neat'][0], c='g', linewidth=0.75, linestyle='--', label='Continuum')
ax6.plot(10**data['log10lam'][0], data['emis'][0], c='deepskyblue', linewidth=0.75, linestyle='--', label='Emission')
ax6.plot(10**data['log10lam'][0], data['best'][0], c='crimson', linewidth=0.75, linestyle='--', label='Model')

ax6.plot(10**data['log10lam'][0], data['err'][0], c='grey', linewidth=0.75)

ax6.set_xlim(3450, 7050)
ax6.set_ylim(-1*data['err'][0].max()*10, data['galaxy'][0].max()*1.5)

ps.style(ax6, fontsize=fontsize, labelbottom=False)
ps.ticks(ax6, xmajor=500, ymajor=10, xminor=100, yminor=2)

ax6.set_ylabel(r'Flux [$10^{-17} erg/s/cm^2/\AA$]', fontsize=fontsize, labelpad=12)

ps.legend(ax6, fontsize=fontsize)
ax6.text(0.05, 0.9, '(d)', fontsize=fontsize, ha='center', va='center', transform=ax6.transAxes)
#______________________________________________________________________________
# Full Spectra residual
ax7 = plt.subplot(gs[8:9, 0:10])
y = (data['galaxy'][0] - data['best'])/data['err']
ax7.scatter(10**data['log10lam'][0], y, c='grey', s=0.1, marker='D')

ax7.set_xlim(3450, 7050)
ax7.set_ylim(-30,30)

ps.style(ax7, fontsize=fontsize)
ps.ticks(ax7, xmajor=500, ymajor=25, xminor=100, yminor=5)
ax7.set_xlabel(r'Wavelength [${\rm \AA}$]', fontsize=fontsize)

ax7.set_ylabel(r'$\Delta/\sigma$', fontsize=fontsize)
#______________________________________________________________________________
# H and K lines
ax2 = plt.subplot(gs[0:3, 4:7])
ax2.plot(10**data['log10lam'][0], data['galaxy'][0], c='k', linewidth=1)
ax2.plot(10**data['log10lam'][0], data['neat'][0], c='g', linewidth=0.75, linestyle='--')
#ax2.plot(10**data['log10lam'][0], data['emis'][0], c='deepskyblue', linewidth=0.75, linestyle='--')
ax2.plot(10**data['log10lam'][0], data['best'][0], c='crimson', linewidth=0.75, linestyle='--')

ax2.plot(10**data['log10lam'][0], data['err'][0], c='grey', linewidth=0.75)

x1, x2 = 3850, 4050
ax2.set_xlim(x1, x2)

ix = (10**data['log10lam'][0]>x1)&(10**data['log10lam'][0]<=x2)
mi = np.nanmin(data['galaxy'][0][ix])*0.75
ma = np.nanmax(data['galaxy'][0][ix])*1.25
ax2.set_ylim(mi, ma)

# Mark box on full spectra
ax6.plot([x1, x2, x2, x1, x1], [mi, mi, ma, ma, mi], c='k', linestyle=':')

ps.style(ax2, fontsize=fontsize, labelbottom=False)
ps.ticks(ax2, xmajor=100, ymajor=5, xminor=20, yminor=1)

ax2.set_title('Ca II H & K', fontsize=fontsize)
ax2.text(0.1, 0.9, '(b)', fontsize=fontsize, ha='center', va='center', transform=ax2.transAxes)
#______________________________________________________________________________
# H and K lines Residuals
ax3 = plt.subplot(gs[3:4, 4:7])
y = (data['galaxy'][0] - data['best'])/data['err']
ax3.scatter(10**data['log10lam'][0], y, c='grey', s=0.1, marker='D')

ax3.set_xlim(x1, x2)
ax3.set_ylim(-30,30)

ps.style(ax3, fontsize=fontsize)
ps.ticks(ax3, xmajor=100, ymajor=25, xminor=20, yminor=5)
#______________________________________________________________________________
# Halpha
ax4 = plt.subplot(gs[0:3, 7:10])
ax4.plot(10**data['log10lam'][0], data['galaxy'][0], c='k', linewidth=1)
ax4.plot(10**data['log10lam'][0], data['neat'][0], c='g', linewidth=0.75, linestyle='--')
#ax4.plot(10**data['log10lam'][0], data['emis'][0], c='deepskyblue', linewidth=0.75, linestyle='--')
ax4.plot(10**data['log10lam'][0], data['best'][0], c='crimson', linewidth=0.75, linestyle='--')

ax4.plot(10**data['log10lam'][0], data['err'][0], c='grey', linewidth=0.75)

x1, x2 = 6515, 6775
ax4.set_xlim(x1, x2)

ix = (10**data['log10lam'][0]>x1)&(10**data['log10lam'][0]<=x2)
mi = np.nanmin(data['galaxy'][0][ix])*0.75
ma = np.nanmax(data['galaxy'][0][ix])*1.25
ax4.set_ylim(mi, ma)

# Mark box on full spectra
ax6.plot([x1, x2, x2, x1, x1], [mi, mi, ma, ma, mi], c='k', linestyle=':')

ps.style(ax4, fontsize=fontsize, labelbottom=False)
ps.ticks(ax4, xmajor=100, ymajor=5, xminor=20, yminor=1)

ax4.set_title(r'H$\alpha$, NII, and SII', fontsize=fontsize)
ax4.text(0.1, 0.9, '(c)', fontsize=fontsize, ha='center', va='center', transform=ax4.transAxes)
#______________________________________________________________________________
# Halpha Residuals
ax5 = plt.subplot(gs[3:4, 7:10])
y = (data['galaxy'][0] - data['best'])/data['err']
ax5.scatter(10**data['log10lam'][0], y, c='grey', s=0.1, marker='D')

ax5.set_xlim(x1, x2)
ax5.set_ylim(-30,30)

ps.style(ax5, fontsize=fontsize)
ps.ticks(ax5, xmajor=100, ymajor=25, xminor=20, yminor=5)
#______________________________________________________________________________
plt.subplots_adjust(wspace=0.6, hspace=0.0)
plt.savefig('spectra.pdf', bbox_inches='tight')
