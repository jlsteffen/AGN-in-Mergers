#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 11:20:07 2021

create maps of BPT and WHAN classifications with their respective line fluxes
and equivalent widths.

@author: joshua
"""

from astropy.io import fits
import numpy as np
import sys, os
import functions as f
from glob import glob
import warnings
warnings.filterwarnings("ignore")

linename = np.array(['OII3730', 'NeIII3870', 'NeIII3969', 'Hg4342', 'OIII4364',
            'HeII4687', 'Hb4863', 'OIII5008', 'NI5199', 'NaI5892',
            'NaI5898', 'OI6302', 'Ha6565', 'NII6585', 'SII6718',
            'SII6733', 'ArIII7138'])

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = os.path.abspath(os.path.join(filepath, '..'))


outpath = filepath + '/map/' 
if not os.path.exists(outpath):
    os.makedirs(outpath)

drpall = fits.getdata(datapath+'/0.drpall/drpall_clean_weight.fits')

files = glob(datapath + '/1.spfit/r1kpc/*[!_spec].fits')

for i in range(len(files)):
    
    fil = fits.open(files[i])
    plateifu = os.path.basename(files[i]).strip('.fits')    
    if not os.path.exists(outpath+'/'+plateifu+'.fits'):
        data = fil[1].data
        anc = fil[2].data
        
        linename = anc['linename'][0]
        binmap = anc['binmap'][0]
        bins = binmap.ravel()
        
        ha = data['flux'][:, 0, (linename=='Ha6565')]
        hb = data['flux'][:, 0, (linename=='Hb4863')]
        o3 = data['flux'][:, 0, (linename=='OIII5008')]
        n2 = data['flux'][:, 0, (linename=='NII6585')]
        ew = data['ew'][:, (linename=='Ha6565')]
        
        bpt = f.bpt(o3, hb, n2, ha)
        bptwhan = f.bpt_whan(o3, hb, n2, ha, ew)
        whan = f.whan(n2, ha)
        
        hamap = np.zeros(binmap.shape)*np.nan
        hbmap = np.zeros(binmap.shape)*np.nan
        o3map = np.zeros(binmap.shape)*np.nan
        n2map = np.zeros(binmap.shape)*np.nan
        ewmap = np.zeros(binmap.shape)*np.nan
        bptmap = np.zeros(binmap.shape)*np.nan
        bptwhanmap = np.zeros(binmap.shape)*np.nan
        whanmap = np.zeros(binmap.shape)*np.nan
        for j in range(bins.max() + 1):
            y, x = np.where(binmap == j)
            hamap[y,x] = ha[j]
            hbmap[y,x] = hb[j]
            o3map[y,x] = o3[j]
            n2map[y,x] = n2[j]
            ewmap[y,x] = ew[j]
            bptmap[y,x] = bpt[j]
            bptwhanmap[y,x] = bptwhan[j]
            whanmap[y,x] = whan[j]
            
        cols = []
        cols.append(fits.Column(name='Ha6565', format=str(len(hamap))+'D', array=hamap))
        
        cols.append(fits.Column(name='Hb4863', format=str(len(hamap))+'D', array=hbmap))
        cols.append(fits.Column(name='OIII5008', format=str(len(hamap))+'D', array=o3map))
        cols.append(fits.Column(name='NII6585', format=str(len(hamap))+'D', array=n2map))
        cols.append(fits.Column(name='EwHa', format=str(len(hamap))+'D', array=ewmap))
        
        cols.append(fits.Column(name='BPT', format=str(len(hamap))+'D', array=bptmap))
        cols.append(fits.Column(name='WHAN', format=str(len(hamap))+'D', array=whanmap))
        cols.append(fits.Column(name='BPTWHAN', format=str(len(hamap))+'D', array=bptwhanmap))
        
        cols = fits.ColDefs(cols)
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.writeto(outpath+'/'+plateifu+'.fits', overwrite=True)
        
    f.update_progress((i+1.0)/np.float64(len(files)))