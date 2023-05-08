#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 11:22:15 2021

Create the BPT and WHAN diagrams for the pair sample.

@author: joshua
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import plotstyle as ps

import warnings
warnings.filterwarnings("ignore")

def blue(x):
    b = 3.1682
    m = 0.16
    return b - (x+18) * m
def red(x):
    b = 4.7866
    m = 0.04
    return b - (x+18) * m

def kewley01(x):
    return ( 0.61 / (x - 0.47) ) + 1.19
def kauffmann03(x):
    return ( 0.61 / (x - 0.05) ) + 1.30
def schawinski07(x):
    return (1.05 * x) + 0.45


filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = filepath + '/0.data'

drpall = fits.getdata(datapath+'/3.match/drpall.fits')

ncomp = fits.getdata(filepath+'/1.pair.fits')
spec = fits.getdata(datapath+'/0.manga_specobj.fits')
spar = fits.getdata(datapath+'/3.match/r1kpc.fits')
#______________________________________________________________________________

BPT = np.concatenate((ncomp['bptt'], ncomp['bptc']), axis=0)
o3hb = np.concatenate((ncomp['o3hbt'], ncomp['o3hbc']), axis=0)
n2ha = np.concatenate((ncomp['n2hat'], ncomp['n2hac']), axis=0)
ewha = np.concatenate((ncomp['what'], ncomp['whac']), axis=0)
#______________________________________________________________________________
pairs = np.concatenate(((ncomp['npair']>=2), (ncomp['npair']>=2)), axis=0)
dagn = np.concatenate(((ncomp['dagn']==1), (ncomp['dagn']==1)), axis=0)

sf = (BPT==1)&(ewha>=6)
comp = (BPT==2)&(ewha>=6)
lin = (BPT==3)&(ewha>=6)
sey = (BPT==4)&(ewha>=6)
ret = (BPT==0)|(ewha<6)
amb = (BPT==-1)&(ewha>=6)

ewha = np.log10(ewha)
#______________________________________________________________________________
fig, (ax2, ax3) = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(20, 10))
fontsize=30.0

x = -1.25, 0.75
y = -1.45, 1.45
s = 15
"""BPT Diagram"""
ax2.scatter(n2ha[pairs&ret], o3hb[pairs&ret], c='gray', label='Retired (RG)', marker='D', s=s)
ax2.scatter(n2ha[pairs&sf], o3hb[pairs&sf], c='deepskyblue', label='Starforming (SF)', marker='D', s=s)
ax2.scatter(n2ha[pairs&comp], o3hb[pairs&comp], c='limegreen', label='Composite (Comp.)', marker='D', s=s)
ax2.scatter(n2ha[pairs&lin], o3hb[pairs&lin], c='orange', label='LINER', marker='D', s=s)
ax2.scatter(n2ha[pairs&sey], o3hb[pairs&sey], c='crimson', label='Seyfert', marker='D', s=s)

ax2.scatter(n2ha[dagn], o3hb[dagn], c='k', marker='s', s=75, label='dAGN')
ax2.scatter(n2ha[dagn&sf], o3hb[dagn&sf], c='deepskyblue', marker='+', s=55)
ax2.scatter(n2ha[dagn&comp], o3hb[dagn&comp], c='limegreen', marker='+', s=55)
ax2.scatter(n2ha[dagn&lin], o3hb[dagn&lin], c='orange', marker='+', s=55)
ax2.scatter(n2ha[dagn&sey], o3hb[dagn&sey], c='crimson', marker='+', s=55)

leg = ax2.legend(fontsize=fontsize/1.75, loc=3, fancybox=False, edgecolor='k')
leg.get_frame().set_linewidth(2.0)

xx = np.arange(x[0], 0.45, 0.01)
yy = kewley01(xx)
ax2.plot(xx, yy, c='k', linewidth=2)
xx = np.arange(x[0], 0.0, 0.01)
yy = kauffmann03(xx)
ax2.plot(xx, yy, c='k', linewidth=2, linestyle='--')
xx = np.arange(-0.18, x[1], 0.01)
yy = schawinski07(xx)
ax2.plot(xx, yy, c='k', linewidth=2, linestyle='-.')

ax2.set_xlim(x)
ax2.set_ylim(y)

ax2.set_ylabel(r'log([O III]/H$\beta$)', fontsize=fontsize)
ax2.set_xlabel(r'log([N II]/H$\alpha$)', fontsize=fontsize)

ps.ticks(ax2, xmajor=0.5, ymajor=0.5, xminor=0.1, yminor=0.1)
ps.style(ax2, fontsize=fontsize)

ax2.annotate('SF', (-1, -0.25), fontsize=fontsize/1.25, fontweight='bold')
ax2.annotate('Comp.', (-0.175, -0.9), fontsize=fontsize/1.25, fontweight='bold')
ax2.annotate('LINER', (0.25, -0.55), fontsize=fontsize/1.25, fontweight='bold')
ax2.annotate('Seyfert', (-0.5, 1.15), fontsize=fontsize/1.25, fontweight='bold')

#______________________________________________________________________________
"""WHAN Diagram"""
ax3.scatter(n2ha[pairs&ret], ewha[pairs&ret], c='gray', label='Retired', marker='D', s=s)
ax3.scatter(n2ha[pairs&sf], ewha[pairs&sf], c='deepskyblue', label='Starforming', marker='D', s=s)
ax3.scatter(n2ha[pairs&comp], ewha[pairs&comp], c='limegreen', label='Composite', marker='D', s=s)
ax3.scatter(n2ha[pairs&lin], ewha[pairs&lin], c='orange', label='AGN', marker='D', s=s)
ax3.scatter(n2ha[pairs&sey], ewha[pairs&sey], c='crimson', label='AGN', marker='D', s=s)

ax3.scatter(n2ha[dagn], ewha[dagn], c='k', marker='s', s=75)
ax3.scatter(n2ha[dagn&sf], ewha[dagn&sf], c='deepskyblue', marker='+', s=55)
ax3.scatter(n2ha[dagn&comp], ewha[dagn&comp], c='limegreen', marker='+', s=55)
ax3.scatter(n2ha[dagn&lin], ewha[dagn&lin], c='orange', marker='+', s=55)
ax3.scatter(n2ha[dagn&sey], ewha[dagn&sey], c='crimson', marker='+', s=55)

ax3.axhline(np.log10(3), c='k')
ax3.plot([-0.4, 0.75], [np.log10(6), np.log10(6)], c='k', linestyle=':')
ax3.plot([-0.4, -0.4], [np.log10(3), 2.75], c='k', linestyle='--')

ax3.set_xlim(x)
ax3.set_ylim(-0.75, 2.75)

ax3.set_ylabel(r'log(EW(H$\alpha$)', fontsize=fontsize)
ax3.set_xlabel(r'log([N II]/H$\alpha$)', fontsize=fontsize)

ps.ticks(ax3, xmajor=0.5, ymajor=0.5, xminor=0.1, yminor=0.1)
ps.style(ax3, fontsize=fontsize)

ax3.annotate('SF', (-1, 0.75), fontsize=fontsize/1.25, fontweight='bold')
ax3.annotate('sAGN', (0.2, 1.0), fontsize=fontsize/1.25, fontweight='bold')
ax3.annotate('wAGN', (0.35, 0.6), fontsize=fontsize/1.25, fontweight='bold')
ax3.annotate('RG', (-1, -0.25), fontsize=fontsize/1.25, fontweight='bold')
#______________________________________________________________________________
plt.savefig(filepath + '/2.bpt_whan.pdf', bbox_inches='tight', overwrite=True)

plt.close('all')