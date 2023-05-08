#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 14:38:01 2021

Make Figure 5 from Fu18

1. Plot 1/Vmax weights (eq3) for all pairs, highlight ones which are AGN
2. Sum AGN pair in sep (kpc) bins
3. Use eqs in Fu+18 to show the expected sums

For the error bars use bootstrap resampling for both the observed and expected
sums. In the expected sum errors I need to include the error from the model 
fit. 

Is there a dependence on mass ratio?

@author: joshua
"""

import os, sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import plotstyle as ps
import warnings
from astropy.stats import bootstrap

warnings.filterwarnings("ignore", message="divide by zero encountered in divide")
warnings.filterwarnings("ignore", message="divide by zero encountered") 
warnings.filterwarnings("ignore", message="invalid value encountered")

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = filepath + '/0.data'

# AGN fraction model parameters
file = open(filepath + '/20model.txt').read().split('\n')
f0, b, sig = [float(i.split(' = ')[1]) for i in file]
#______________________________________________________________________________
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
def make_error_boxes(ax, xdata, ydata, xerror, yerror, facecolor='deepskyblue',
                     edgecolor='None', alpha=0.5):
    # Loop over data points; create box from errors at each point
    errorboxes = [Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
                  for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T)]
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor)
    # Add collection to axes
    ax.add_collection(pc)
    # Plot errorbars
    #artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,
    #                      fmt='None', ecolor='b', label='Expected')
    return pc
#______________________________________________________________________________
def linfit(x, a, b):
    return a * x + b
def powfit(x, a, b, c):
    return a*x**-b + c
#______________________________________________________________________________
# AGN model fraction function
def mod(M, z, f0=f0, b=b, sigma=sig):
    return f0 * np.exp(-0.5*((M - b)/(sigma))**2) * (1+z)**4
#______________________________________________________________________________
def boot(data, cycles=1000, resample=0.68, func=np.mean):
    if len(data) > 0:
        mask = ~np.isnan(data)
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
def eb(ax, x, y, xerr, yerr, c, fmt='D', label=None, linestyle='-'):
    '''
    My own stylized errorbars
    '''
    er = ax.errorbar(x, y, xerr=xerr, fmt=fmt, c=c)
    er[-1][0].set_linestyle(linestyle)
    ax.errorbar(x, y, yerr=yerr, fmt=fmt, c=c, capsize=4, capthick=2, label=label)
def eb2(ax, x, y, xerr, yerr, c, fmt='D', label=None, linestyle='-'):
    '''
    My own stylized errorbars
    '''
    er = ax.errorbar(x, y, xerr=xerr, marker=fmt, c=c, ls=':')
    er[-1][0].set_linestyle(linestyle)
    ax.errorbar(x, y, yerr=yerr, marker=fmt, c=c, capsize=4, capthick=2, label=label, ls=':')
#______________________________________________________________________________
drpall = fits.getdata(datapath+'/0.drpall/drpall_clean_weight.fits')
pair = fits.getdata(filepath+'/1.pair.fits')
fil = fits.getdata(datapath+'/DR14/drpall-v2_1_2.fits')

sepkpc = pair['sepkpc']

ew = drpall['esweight']
mnsa = drpall['nsa_elpetro_mass']/0.7**2
z = drpall['nsa_z']

dm = 10**pair['dlogm']
m1 = np.log10(mnsa/(1 + 1/dm))
m2 = np.log10(mnsa - 10**m1)
#______________________________________________________________________________
# Sample Definition
seplim = 25
e = (pair['eagn']&(sepkpc<=seplim)).astype('int')
s = (pair['sagn']&(sepkpc<=seplim)).astype('int')
d = (pair['dagn']&(sepkpc<=seplim)).astype('int')
par = ((pair['npair']==2)&(sepkpc<=seplim)).astype('int')

# AGN model for primaries and secondaries
fagn = mod(m1, z)
fagns = mod(m2, z)
#______________________________________________________________________________
# Create adaptive bin sizes so that each bin contains the same number of pairs
btest = [0]
data = sorted(sepkpc[(par==1)&(~np.isnan(ew))&(ew!=-99)])

numbin = 5
num = len(data) // numbin + 0

sel = num
for i in range(numbin):
    if sel < len(data):
        btest.append(data[sel])
        sel+=num
    else:
        btest.append(data[-1])
rpbin = btest
bs = np.array([btest[i+1] - btest[i] for i in range(len(btest)-1)])
#______________________________________________________________________________
# equations
def eq_sagn(w, fagn, fagns, lam, xi):
    s1 = ((1-xi)*w*(1+lam)*fagn*(1-(1+lam)*fagns)).sum()
    s2 = ((1-xi)*w*(1+lam)*fagns*(1-(1+lam)*fagn)).sum()
    return s1+s2
def eq_dagn(w, fagn, fagns, lam, xi):
    d1 = (w*(1+lam)**2*fagn*fagns).sum()
    d2 = ((xi)*w*(1+lam)*fagn*(1-(1+lam)*fagns)).sum()
    d3 = ((xi)*w*(1+lam)*fagns*(1-(1+lam)*fagn)).sum()
    return d1+d2+d3
def eq_eagn(w, fagn, fagns, lam, xi):
    #return eq_sagn(w, fagn, fagns, lam, xi) + eq_dagn(w, fagn, fagns, lam, xi)
    e1 = (w*(1+lam)*fagn).sum()
    e2 = (w*(1+lam)*fagns).sum()
    e3 = (w*(1+lam)**2*fagn*fagns).sum()
    return e1+e2-e3
#______________________________________________________________________________   
va,vb,vc,vd = 0,0,0,0

# observed
obe = []
eobe = []
obs = []
eobs = []
obd = []
eobd = []
# expected
exe = []
exs = []
eexs = []
exd = []
for i in range(len(rpbin)-1):
    # observed
    ix = (par==1)&(sepkpc<rpbin[i+1])&(sepkpc>rpbin[i])&(~np.isnan(ew))&(ew!=-99)

    lam = powfit(sepkpc[ix], va, vb, 0) # this seems to work for S+D
    xi = powfit(sepkpc[ix], vc, vd, 0)            
        
    # observed
    dat = ew[ix&(d==1)]
    if len(dat)>1:
        obd.append((dat).sum())
    else:
        obd.append(np.nan)
    eobd.append(boot(np.log10(dat)))
    dat = ew[ix&(s==1)]
    if len(dat)>1:
        obs.append((dat).sum())
    else:
        obs.append(np.nan)   
    eobs.append(boot(np.log10(dat)))
    dat = ew[ix&(e==1)]
    if len(dat)>1:
        obe.append((dat).sum())
    else:
        obe.append(np.nan)
    eobe.append(boot(np.log10(dat)))
    
    # expected
    exs.append(eq_sagn(ew[ix], fagn[ix], fagns[ix], lam, xi))
    exd.append(eq_dagn(ew[ix], fagn[ix], fagns[ix], lam, xi))
    exe.append(eq_eagn(ew[ix], fagn[ix], fagns[ix], lam, xi))
    
    # Expected Errors
    dat = ew[ix]
    eexs.append(boot(np.log10(dat)))


obe = np.log10(np.array(obe))
eobe = np.array(eobe).T
exe = np.log10(np.array(exe))

obs = np.log10(np.array(obs))
obd = np.log10(np.array(obd))
eobs = np.array(eobs).T
eobd = np.array(eobd).T
exs = np.log10(np.array(exs))
exd = np.log10(np.array(exd))

eexs = np.array(eexs).T
eexs = (eexs**2 + 0.04**2)**0.5 # error from best-fit AGN fraction model
"""""""""""
Plotting
"""""""""""

fontsize = 22
fig, ax = plt.subplots(1,1, figsize=(6,6))

ellison13 = np.loadtxt(filepath+'/8.compare/Ellison13_fig10.csv', delimiter=',')
ellison13[0,0] = -5*0.7
ellison13_err = np.abs(np.log10(ellison13[:,2])-np.log10(ellison13[:,1]))

satyapal14 = np.loadtxt(filepath+'/8.compare/Satyapal14_fig5_blue.csv', delimiter=',')
satyapal14[0,0] = -5*0.7
satyapal14_err = np.abs(np.log10(satyapal14[:,2])-np.log10(satyapal14[:,1]))

shah20 = np.loadtxt(filepath+'/8.compare/Shah20_fig15.csv', delimiter=',')
shah20[0,0] = -5
shah20_lerr = np.abs(np.log10(shah20[:,2])-np.log10(shah20[:,1]))
shah20_uerr = np.abs(np.log10(shah20[:,3])-np.log10(shah20[:,1]))

mcalpine20 = np.loadtxt(filepath+'/8.compare/McAlpine20.csv', delimiter=',')
mcalpine20_err = np.abs(np.log10(mcalpine20[:,2]) - np.log10(mcalpine20[:,1]))

ax.fill_between([-10,0], y1=[-0.2,-0.2], y2=[1.25,1.25],facecolor='grey', alpha=0.5)

eb(ax, rpbin[1::]-bs/2, (obe-exe), xerr=bs/2, yerr=(eobe**2+eexs**2)**0.5, c='k', fmt='D', label='This Work (offset+dual AGN)')

eb(ax, ellison13[:,0]/0.7, np.log10(ellison13[:,1]), xerr=0, yerr=ellison13_err, c='deepskyblue', fmt='s')
eb(ax, satyapal14[:,0]/0.7, np.log10(satyapal14[:,1]), xerr=0, yerr=satyapal14_err, c='crimson', fmt='s')
eb(ax, shah20[:,0], np.log10(shah20[:,1]), xerr=0, yerr=[shah20_lerr, shah20_uerr], c='blueviolet', fmt='o')
eb(ax, mcalpine20[:,0], np.log10(mcalpine20[:,1]), xerr=0, yerr=mcalpine20_err, c='gold', fmt='P')

ax.scatter(ellison13[:,0]/0.7, np.log10(ellison13[:,1]), c='deepskyblue', marker='s', label='Ellison+13 (Optical)')
ax.scatter(satyapal14[:,0]/0.7, np.log10(satyapal14[:,1]), c='crimson', marker='s', label='Satyapal+14 (Optical+IR)')
ax.scatter(shah20[:,0], np.log10(shah20[:,1]), c='blueviolet', marker='o', label='Shah+20 (X-Ray)')
ax.scatter(mcalpine20[:,0], np.log10(mcalpine20[:,1]), c='gold', marker='P', label='McAlpine+20 (Eagle Simulation)')

ps.legend(ax, fontsize=fontsize)

ax.plot([-7,150], [0,0], c='k', linestyle='--')

ax.set_ylabel('AGN Excess (dex)', fontsize=fontsize)
ax.set_xlabel(r'r$_{\rm p}$ (kpc)', fontsize=fontsize)

ax.set_xlim(-10, 105)
ax.set_ylim(-0.2, 1.25)

ps.style(ax, fontsize)
ps.ticks(ax, xmajor=20, xminor=5, ymajor=0.5, yminor=0.1)

plt.savefig(filepath+'/8.compare.pdf', bbox_inches='tight')
plt.close('all')