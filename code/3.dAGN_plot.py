#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 13:29:45 2021

Loop over each dAGN
print an image of the field w/ objects overlaid
print a bptmap of the field
print a bpt diagram of the field

@author: joshua
"""
print()
print('Initiallizing...')
print()

from astropy.io import fits
import numpy as np
import sys, os
import math as m
from PIL import Image
from astropy import wcs
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from copy import copy
import functions as f
import plotstyle as ps

import warnings
warnings.filterwarnings("ignore")

#______________________________________________________________________________
def kewley01(x):
    return ( 0.61 / (x - 0.47) ) + 1.19
def kauffmann03(x):
    return ( 0.61 / (x - 0.05) ) + 1.30
def schawinski07(x):
    return (1.05 * x) + 0.45
#______________________________________________________________________________
def png(im, ax, plateifu):
    cent_pix = np.array([im.size[0]/2.0, im.size[1]/2.0])
    # Make a Astropy WCS
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [cent_pix[0], cent_pix[1]]
    w.wcs.cdelt = np.array([-scale/3600.0, -scale/3600.0])
    w.wcs.crval = [ra, dec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ['deg', 'deg']
    w.wcs.set_pv([(2, 1, 45.0)])
    ax.imshow(im, aspect='auto')
    ax.text(s=plateifu, x=60, y=30, color='w', fontsize=30)
    # IFU hexagons    
    hexagon = patches.RegularPolygon([size/2.0, size/2.0], 6,linestyle='--',
              radius = size/2.0 * ifu_scale/40.0, ec = 'lime', fc = 'none', orientation=m.pi/2)
    ax.add_patch(copy(hexagon))
    
    # Label the dAGN
    markers = ['o', 's', 'H']
    for j in range(len(ix)):
        # Convert object ra and dec to by in pixel coordinates
        crd = np.array([[obj_ra[ix][j], obj_dec[ix][j]]], np.float_)        
        pix = w.wcs_world2pix(crd, 1)
        obj_x = pix[0][0]
        obj_y = pix[0][1]
        ax.scatter(obj_x, obj_y, marker = markers[j], s=250, facecolors='none', edgecolors='lime')
     
    loc1 = 75
    loc2 = 15
    ps.style(ax, fontsize=fontsize, labelbottom=False, labelleft=False, c='w')
    ps.ticks(ax, xmajor=loc1, ymajor=loc1, xminor=loc2, yminor=loc2)
#______________________________________________________________________________
def bptmap(bpt, ax, im):
    x1, x2 = (size-ifuscale[ind])/2, ifuscale[ind] + (size-ifuscale[ind])/2
    # BPT map
    col_dict={-1:"dimgray",
               0:"navy",
               1:"cyan",
               2:"g",
               3:"orange",
               4:"crimson"}

    labels = np.array(["Amb","Ret","SF","Comp","LINER", "Sey"])
    cm, fmt, tickz = ps.cbar_descrete(col_dict, labels)
    
    #ax2.set_title('BPT-WHAN Classification', fontsize=fontsize)
    
    im2 = ax.imshow(bpt, origin='lower', cmap=cm, extent=[x1,x2,x2,x1], vmin=-1.5, vmax=4.5, alpha=0.5, zorder=1)

    # Collage colorbar
    cbar_ax1 = fig.add_axes([0.4, 0.88, 0.225, 0.02])
    cb1 = fig.colorbar(im2, format=fmt, ticks=tickz, cax=cbar_ax1, orientation='horizontal')
    ps.cbar_style(cbar=cb1, fontsize=fontsize*1.25, side='top')        
    
    cbar_ax2 = fig.add_axes([0.67, 0.88, 0.225, 0.02])
    cb2 = fig.colorbar(im2, format=fmt, ticks=tickz, cax=cbar_ax2, orientation='horizontal')
    ps.cbar_style(cbar=cb2, fontsize=fontsize*1.25, side='top')

    cb1.ax.xaxis.set_label_position('top')
    cb1.ax.xaxis.set_ticks_position('top')
    
        
    loc1 = 75.0
    loc2 = 15.0
    
    ps.style(ax, fontsize=fontsize, labelbottom=False, labelleft=False, c='w')
    ps.ticks(ax, xmajor=loc1, ymajor=loc1, xminor=loc2, yminor=loc2)
    
    ax.imshow(im, aspect='auto', extent=[0,size,size,0], zorder=0)
    # IFU hexagons    
    hexagon = patches.RegularPolygon([size/2.0, size/2.0], 6,linestyle='--',
              radius = size/2.0 * ifu_scale/40.0, ec = 'lime', fc = 'none', orientation=m.pi/2)
    ax.add_patch(copy(hexagon))
    cent_pix = np.array([im.size[0]/2.0, im.size[1]/2.0])
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [cent_pix[0], cent_pix[1]]
    w.wcs.cdelt = np.array([-scale/3600.0, -scale/3600.0])
    w.wcs.crval = [ra, dec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ['deg', 'deg']
    w.wcs.set_pv([(2, 1, 45.0)])
    
    # Label the dAGN
    markers = ['o', 's', 'H']
    for j in range(len(ix)):
        # Convert object ra and dec to by in pixel coordinates
        crd = np.array([[obj_ra[ix][j], obj_dec[ix][j]]], np.float_)        
        pix = w.wcs_world2pix(crd, 1)
        obj_x = pix[0][0]
        obj_y = pix[0][1]
        ax.scatter(obj_x, obj_y, marker = markers[j], s=250, facecolors='none', edgecolors='lime', linewidth=1)
        #ax.text(obj_x-5, obj_y - 17, str(j), fontsize=fontsize/1.5, color='w')
    
#______________________________________________________________________________
def bptplot(n2ha, o3hb, bpt, ax):
    # BPT diagram
    sf = bpt==1
    comp = bpt==2
    sey = bpt==4
    lin = bpt==3
    ret = bpt==0
    amb = bpt < 0
    
    x = -1.25, 0.75
    y = -1.45, 1.45
        
    ax.scatter(n2ha[amb], o3hb[amb], c='dimgray', label='Ambiguous', marker='D', s=15)
    ax.scatter(n2ha[ret], o3hb[ret], c='navy', label='Retired', marker='D', s=15)
    ax.scatter(n2ha[sf], o3hb[sf], c='cyan', label='Starforming', marker='D', s=15)
    ax.scatter(n2ha[comp], o3hb[comp], c='g', label='Composite', marker='D', s=15)
    ax.scatter(n2ha[sey], o3hb[sey], c='crimson', label='Seyfert', marker='D', s=15)
    ax.scatter(n2ha[lin], o3hb[lin], c='orange', label='LINER', marker='D', s=15)
    
    # Target and companion marker on diagram
    markers = ['o', 's', 'H']
    for j in range(len(N2HA[sx][ix])):
        ax.scatter(N2HA[sx][ix][j], O3HB[sx][ix][j], c='w', edgecolor='k', marker=markers[j], s=150)
        ax.scatter(N2HA[sx][ix][j], O3HB[sx][ix][j], c='k', marker='+', s=150)
        
    xx = np.arange(x[0], 0.45, 0.01)
    yy = kewley01(xx)
    ax.plot(xx, yy, c='k', linewidth=2)
    xx = np.arange(x[0], 0.0, 0.01)
    yy = kauffmann03(xx)
    ax.plot(xx, yy, c='k', linewidth=2, linestyle='--')
    xx = np.arange(-0.18, x[1], 0.01)
    yy = schawinski07(xx)
    ax.plot(xx, yy, c='k', linewidth=2, linestyle='-.')
    
    ax.set_xlim(x)
    ax.set_ylim(y)
    
    ax.set_ylabel(r'log([O III]/H$\beta$)', fontsize=fontsize)
    ax.set_xlabel(r'log([N II]/H$\alpha$)', fontsize=fontsize)
    
    ps.ticks(ax, xmajor=0.5, ymajor=0.5, xminor=0.1, yminor=0.1)
    ps.style(ax, fontsize=fontsize, labelbottom=False, labelleft=False)  
#______________________________________________________________________________
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = filepath + '/0.data'

drpall = fits.getdata(datapath+'/0.drpall/drpall_clean_weight.fits')
plifu = drpall['plateifu']
ifuscale = drpall['ifurad']*2/40 * 300

spobj = fits.getdata(datapath+'/0.manga_specobj.fits')
index = spobj['index']
maintarg = spobj['maintarg']
spl = spobj['plateifu']

psam = fits.getdata(filepath+'/1.pair.fits')
dagn = (psam['dagn']==1)
cagn = (psam['cagn']==1)
plateifu = plifu[cagn]

bptfil = fits.getdata(datapath+'/3.match/bptw_r1kpc.fits')
N2HA = bptfil['N2HA']
O3HB = bptfil['O3HB']
#______________________________________________________________________________
# Information from the drpall
ras = drpall['ifura']
decs = drpall['ifudec']
ifus = drpall['ifudesignsize']

outpath = filepath + '/3.AGN_png'
if not os.path.exists(outpath):
    os.makedirs(outpath)
#______________________________________________________________________________
# pixel Size of the image and radius size of targeting circles
size = 300
scale = 40.0/size
fontsize=20
plateifu = ['7443-12703']
#______________________________________________________________________________

for i in range(len(plateifu)):
    if not os.path.exists(outpath + '/'+str(plateifu[i])+'.1pdf'):
        table = fits.getdata(datapath+'/1.spfit/r1kpc/'+plateifu[i]+'.fits')
        ind = np.where(drpall['plateifu']==plateifu[i])[0][0]
        obj_ra = table['ra']
        obj_dec = table['dec']
        sx = (spl==plateifu[i])
        ix = [np.where(maintarg[sx]==1)[0][0], index[psam['pindex'][ind]]]
#______________________________________________________________________________
        # Download images from skyserver    
        ra = ras[ind]
        dec = decs[ind]
        
        IFU_size = int(drpall['ifudsgn'][ind][0:-2])
        IFU_arcsec = np.array(([19, 37, 61, 91, 127], [16.0, 21.0, 26.0, 31.0, 36.0]))
        ifu_scale = IFU_arcsec[1, np.where(IFU_size==IFU_arcsec[0])[0][0]]
        
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=False, sharex=False, figsize=(20, 6))
#______________________________________________________________________________
        # Open the downloaded png for editting
        im = Image.open(datapath +'/0.thumbnails/png/'+str(plateifu[i])+'.png')
        png(im, ax1, plateifu[i])
        
        # BPT diagram
        fil = fits.getdata(datapath + '/6.bptmap/map/'+plateifu[i]+'.fits')
        bpt = fil['bptwhan']
        ha = fil['ha6565']
        hb = fil['hb4863']
        o3 = fil['oiii5008']
        n2 = fil['nii6585']

        bptmap(bpt, ax2, im)
        
        n2ha = np.log10(np.divide(n2, ha))
        o3hb = np.log10(np.divide(o3, hb))
        
        bptplot(n2ha, o3hb, bpt, ax3)
#______________________________________________________________________________
        plt.subplots_adjust(wspace = 0.15)
        plt.savefig(outpath + '/'+str(plateifu[i])+'.png', bbox_inches='tight')
        plt.close('all')
        
        f.update_progress((i+1.0)/np.float64(len(plateifu)))