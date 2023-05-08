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
from astropy.stats import bootstrap

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
filepath = os.path.abspath(os.path.join(filepath, '..'))
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

from scipy.optimize import minimize
method = 'powell'
def lam_fit(params):
    a, b = params
    exp = []
    obs = []
    for i in range(len(rpbin)-1):
        ix = (par==1)&(sepkpc<rpbin[i+1])&(sepkpc>rpbin[i])&(~np.isnan(ew))&(ew!=-99)
        lam = a * sepkpc[ix] ** (-b)
        exp.append(eq_eagn(ew[ix], fagn[ix], fagns[ix], lam, 0))
        
        dat = ew[ix&(e==1)]
        if len(dat)>1:
            obs.append((dat).sum())
        else:
            obs.append(np.nan)
    dif = np.log10(np.array(obs))-np.log10(np.array(exp))
    return np.sum(dif**2)

bnds = [(0,10),(0,2)]
x0 = (2,1)
resl = minimize(lam_fit, x0, bounds=bnds, method=method)
print(resl.x)
resl = minimize(lam_fit, resl.x, bounds=bnds, method=method)
print(resl.x)

def xi_fit(params):
    a, b = params
    exd = []
    obd = []
    exs = []
    obs = []
    for i in range(len(rpbin)-1):
        ix = (par==1)&(sepkpc<rpbin[i+1])&(sepkpc>rpbin[i])&(~np.isnan(ew))&(ew!=-99)
        #lam = resl.x[0] * sepkpc[ix] ** (-resl.x[1])
        xi = a * sepkpc[ix] ** (-b)
        exd.append(eq_dagn(ew[ix], fagn[ix], fagns[ix], 0, xi))
        exs.append(eq_sagn(ew[ix], fagn[ix], fagns[ix], 0, xi))
        
        dat = ew[ix&(d==1)]
        if len(dat)>1:
            obd.append((dat).sum())
        else:
            obd.append(np.nan)
        dat = ew[ix&(s==1)]
        if len(dat)>1:
            obs.append((dat).sum())
        else:
            obs.append(np.nan)

    difd = np.log10(np.array(obd))-np.log10(np.array(exd))
    difs = np.log10(np.array(obs))-np.log10(np.array(exs))
    return np.sum(difd**2+difs**2)

bnds = [(0,np.inf),(1,1)]
x0 = (1,1)
resx = minimize(xi_fit, x0, bounds=bnds, method=method)
print(resx.x)
resx = minimize(xi_fit, resx.x, bounds=bnds, method=method)
print(resx.x)

def comb_fit(params):
    va, vb, vc, vd = params
    exd = []
    obd = []
    exs = []
    obs = []
    exe = []
    obe = []
    for i in range(len(rpbin)-1):
        ix = (par==1)&(sepkpc<rpbin[i+1])&(sepkpc>rpbin[i])&(~np.isnan(ew))&(ew!=-99)
        
        lam = va * sepkpc[ix] ** (-vb)
        xi = vc * sepkpc[ix] ** (-vd)
        
        exd.append(eq_dagn(ew[ix], fagn[ix], fagns[ix], lam, xi))
        exs.append(eq_sagn(ew[ix], fagn[ix], fagns[ix], lam, xi))
        exe.append(eq_eagn(ew[ix], fagn[ix], fagns[ix], lam, 0))
        
        dat = ew[ix&(d==1)]
        if len(dat)>1:
            obd.append((dat).sum())
        else:
            obd.append(np.nan)
        dat = ew[ix&(s==1)]
        if len(dat)>1:
            obs.append((dat).sum())
        else:
            obs.append(np.nan)        
        dat = ew[ix&(e==1)]
        if len(dat)>1:
            obe.append((dat).sum())
        else:
            obe.append(np.nan)
        
    difd = np.log10(np.array(obd))-np.log10(np.array(exd))
    difs = np.log10(np.array(obs))-np.log10(np.array(exs))
    dife = np.log10(np.array(obe))-np.log10(np.array(exe))
    return np.sum(difd**2+difs**2+dife**2)

bnds = [(0,np.inf),(1,1),(0,np.inf),(1,1)]
x0 = (resl.x[0], resl.x[1], resx.x[0], resx.x[1])
x0 = (3.15,1,1.17,1)
resc = minimize(comb_fit, x0, bounds=bnds, method=method)
print(resc.x)
resc = minimize(comb_fit, resc.x, bounds=bnds, method=method)
print(resc.x)

va = np.round(resc.x[0], 2)
vb = np.round(resc.x[1], 2)
vc = np.round(resc.x[2], 2)
vd = np.round(resc.x[3], 2)

# reduced Chi squared minimization
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

difd = obd - exd
difs = obs - exs
dife = obe - exe
    
varrd = (eobd[0]**2+eobd[1]**2+eexs[0]**2+eexs[1]**2)
varrs = (eobs[0]**2+eobs[1]**2+eexs[0]**2+eexs[1]**2)
varre = (eobe[0]**2+eobe[1]**2+eexs[0]**2+eexs[1]**2)

chi2 = (difd**2/varrd).sum() + (difs**2/varrs).sum() + (dife**2/varre).sum()
rchi2 = chi2/(len(difd)+len(difs)+len(dife)-2)
#______________________________________________________________________________    
"""
Uncomment the next line to turn off the model during the plotting
"""
#va,vb,vc,vd = 0,0,0,0

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

#Plotting
fontsize = 18
fig = plt.figure(figsize = (18,7))
gs = plt.GridSpec(16, 48, width_ratios = [1]*48, height_ratios = [1]*16)

ax5 = plt.subplot(gs[0:10, 0:16])
ax6 = plt.subplot(gs[11:16, 0:16])
ax3 = plt.subplot(gs[0:10, 17:32])
ax4 = plt.subplot(gs[11:16, 17:32])
ax1 = plt.subplot(gs[0:10, 33:48])
ax2 = plt.subplot(gs[11:16, 33:48])
a = [ax1, ax2, ax3, ax4, ax5, ax6]

# Highlight regions between sets
x = np.repeat(rpbin, 2)[1:-1]
y1, y2 = np.repeat(obe, 2), np.repeat(exe, 2)
a[0].fill_between(x, y1, y2, where=y1-y2>0, step='post', facecolor='lightgreen', alpha=0.5)
a[0].fill_between(x, y1, y2, where=y1-y2<0, step='post', facecolor='red', alpha=0.5)
a[1].fill_between(x, y1-y2, y2-y2, where=y1-y2>0, step='post', facecolor='lightgreen', alpha=0.5)
a[1].fill_between(x, y1-y2, y2-y2, where=y1-y2<0, step='post', facecolor='red', alpha=0.5)

y1, y2 = np.repeat(obs, 2), np.repeat(exs, 2)
a[2].fill_between(x, y1, y2, where=y1-y2>0, step='post', facecolor='lightgreen', alpha=0.5)
a[2].fill_between(x, y1, y2, where=y1-y2<0, step='post', facecolor='red', alpha=0.5)
a[3].fill_between(x, y1-y2, y2-y2, where=y1-y2>0, step='post', facecolor='lightgreen', alpha=0.5)
a[3].fill_between(x, y1-y2, y2-y2, where=y1-y2<0, step='post', facecolor='red', alpha=0.5)

y1, y2 = np.repeat(obd, 2), np.repeat(exd, 2)
a[4].fill_between(x, y1, y2, where=y1-y2>0, step='post', facecolor='lightgreen', alpha=0.5)
a[4].fill_between(x, y1, y2, where=y1-y2<0, step='post', facecolor='red', alpha=0.5)
a[5].fill_between(x, y1-y2, y2-y2, where=y1-y2>0, step='post', facecolor='lightgreen', alpha=0.5)
a[5].fill_between(x, y1-y2, y2-y2, where=y1-y2<0, step='post', facecolor='red', alpha=0.5)

# Indiv Weights
a[0].scatter(sepkpc[par==1], np.log10(ew[par==1]), c='lightgrey', 
            label='All Pairs', s=10)
a[0].scatter(sepkpc[e==1], np.log10(ew[e==1]), c='lightskyblue', s=20, label='Pairs w/ AGN')

a[2].scatter(sepkpc[par==1], np.log10(ew[par==1]), c='lightgrey', 
            label='All Pairs', s=10)
a[2].scatter(sepkpc[s==1], np.log10(ew[s==1]), c='lightskyblue', s=20, label='Pairs w/ AGN')

a[4].scatter(sepkpc[par==1], np.log10(ew[par==1]), c='lightgrey', 
            label='All Pairs', s=10)
a[4].scatter(sepkpc[d==1], np.log10(ew[d==1]), c='lightskyblue', s=20, label='Pairs w/ AGN')

# Observed bins
a[0].errorbar(rpbin[1::]-bs/2, obe, xerr=bs/2, yerr=eobe, fmt='D', c='k', label='Total AGN Weights', capsize=4, capthick=2)
a[2].errorbar(rpbin[1::]-bs/2, obs, xerr=bs/2, yerr=eobs, fmt='D', c='k', label='Total AGN Weights', capsize=4, capthick=2)
a[4].errorbar(rpbin[1::]-bs/2, obd, xerr=bs/2, yerr=eobd, fmt='D', c='k',  label='Total AGN Weights', capsize=4, capthick=2)

# Expected bins
eb(a[0], rpbin[1::]-bs/2, exe, xerr=bs/2, yerr=eexs, c='dimgrey', fmt='.', linestyle='--', label='Stochastic')
eb(a[2], rpbin[1::]-bs/2, exs, xerr=bs/2, yerr=eexs, c='dimgrey', fmt='.', linestyle='--', label='Stochastic')
eb(a[4], rpbin[1::]-bs/2, exd, xerr=bs/2, yerr=eexs, c='dimgrey', fmt='.', linestyle='--', label='Stochastic')

# Excess plots
eb(a[1], rpbin[1::]-bs/2, obe-exe, xerr=bs/2, yerr=(eobe**2+eexs**2)**0.5, c='k', fmt='D')
eb(a[3], rpbin[1::]-bs/2, obs-exs, xerr=bs/2, yerr=(eobs**2+eexs**2)**0.5, c='k', fmt='D')
eb(a[5], rpbin[1::]-bs/2, obd-exd, xerr=bs/2, yerr=(eobd**2+eexs**2)**0.5, c='k', fmt='D')

a[1].plot([-2,47], [0,0], c='dimgrey', linestyle='--', linewidth=1.0)
a[3].plot([-2,47], [0,0], c='dimgrey', linestyle='--', linewidth=1.0)
a[5].plot([-2,47], [0,0], c='dimgrey', linestyle='--', linewidth=1.0)

a[4].set_ylabel(r'log n (# per $10^6$ Mpc$^3$)', fontsize=fontsize)
a[5].set_ylabel('Excess (dex)', fontsize=fontsize)
a[1].set_xlabel(r'r$_{\rm p}$ (kpc)', fontsize=fontsize)
a[3].set_xlabel(r'r$_{\rm p}$ (kpc)', fontsize=fontsize)
a[5].set_xlabel(r'r$_{\rm p}$ (kpc)', fontsize=fontsize)

ps.legend(a[4], fontsize=fontsize*1.5)

for i in range(6):
    a[i].set_xlim(-2,27)
    ps.ticks(a[i], xmajor=10, ymajor=1, xminor=2, yminor=0.2)
    
a[0].set_ylim(-1.8,1.5)
a[2].set_ylim(-1.8,1.5)
a[4].set_ylim(-1.8,1.5)

a[1].set_ylim(-0.75,2.3)
a[3].set_ylim(-0.75,2.3)
a[5].set_ylim(-0.75,2.3)

ps.style(a[0], fontsize=fontsize*1.5, labelbottom=False, labelleft=False)
ps.style(a[1], fontsize=fontsize*1.5, labelleft=False)
ps.style(a[2], fontsize=fontsize*1.5, labelleft=False, labelbottom=False)
ps.style(a[3], fontsize=fontsize*1.5, labelleft=False)
ps.style(a[4], fontsize=fontsize*1.5, labelbottom=False)
ps.style(a[5], fontsize=fontsize*1.5)

a[0].set_title('Offset + Dual AGN', fontsize=fontsize)
a[2].set_title('Offset AGN', fontsize=fontsize)
a[4].set_title('Dual AGN', fontsize=fontsize)

plt.subplots_adjust(wspace=-0.5, hspace=-0.5)
plt.savefig(filepath + '/1.vol_den_sep.pdf', bbox_inches='tight')
plt.close(fig)
#______________________________________________________________________________
# Compare with Fu18
fontsize = 18
fig, a = plt.subplots(1,3, figsize=(18, 6))

# Excess plots
eb(a[2], rpbin[1::]-bs/2, obe-exe, xerr=bs/2, yerr=(eobe**2+eexs**2)**0.5, c='k', fmt='D', label='This Work')
eb(a[1], rpbin[1::]-bs/2, obs-exs, xerr=bs/2, yerr=(eobs**2+eexs**2)**0.5, c='k', fmt='D', label='This Work')
eb(a[0], rpbin[1::]-bs/2, obd-exd, xerr=bs/2, yerr=(eobd**2+eexs**2)**0.5, c='k', fmt='D', label='This Work')

a[2].plot([-2,47], [0,0], c='dimgrey', linestyle='--', linewidth=1.0)
a[1].plot([-2,47], [0,0], c='dimgrey', linestyle='--', linewidth=1.0)
a[0].plot([-2,47], [0,0], c='dimgrey', linestyle='--', linewidth=1.0)

# Fu18 points
feagn = np.loadtxt('/Volumes/Ext_drive/Research/mpl11/compare/Fu18_eagn.csv', delimiter=',')
fsagn = np.loadtxt('/Volumes/Ext_drive/Research/mpl11/compare/Fu18_sagn.csv', delimiter=',')
fdagn = np.loadtxt('/Volumes/Ext_drive/Research/mpl11/compare/Fu18_dagn.csv', delimiter=',')
a[2].errorbar(feagn[:,0], feagn[:,1], xerr=[2.5,2.5,2.5,2.5,5], yerr=[np.abs(feagn[:,2]-feagn[:,1]),np.abs(feagn[:,3]-feagn[:,1])], fmt='D', c='crimson', label='Fu+18')
a[1].errorbar(fsagn[:,0], fsagn[:,1], xerr=[5,10], yerr=[np.abs(fsagn[:,2]-fsagn[:,1]),np.abs(fsagn[:,3]-fsagn[:,1])], fmt='D', c='crimson', label='Fu+18')
a[0].errorbar(fdagn[:,0], fdagn[:,1], xerr=[5,10], yerr=[np.abs(fdagn[:,2]-fdagn[:,1]),np.abs(fdagn[:,3]-fdagn[:,1])], fmt='D', c='crimson', label='Fu+18')

a[0].set_ylabel('Excess (dex)', fontsize=fontsize)
a[2].set_xlabel(r'r$_{\rm p}$ (kpc)', fontsize=fontsize)
a[1].set_xlabel(r'r$_{\rm p}$ (kpc)', fontsize=fontsize)
a[0].set_xlabel(r'r$_{\rm p}$ (kpc)', fontsize=fontsize)

ps.legend(a[0], fontsize=fontsize*1.5)
for i in range(3):
    a[i].set_xlim(-2,32)
    ps.ticks(a[i], xmajor=10, ymajor=1, xminor=2, yminor=0.2)
    a[i].set_ylim(-0.75,2.3)

ps.style(a[2], fontsize=fontsize*1.5, labelleft=False)
ps.style(a[1], fontsize=fontsize*1.5, labelleft=False)
ps.style(a[0], fontsize=fontsize*1.5)

a[2].set_title('Offset + Dual AGN', fontsize=fontsize)
a[1].set_title('Offset AGN', fontsize=fontsize)
a[0].set_title('Dual AGN', fontsize=fontsize)

plt.subplots_adjust(wspace=-0.0, hspace=-0.0)
plt.savefig(filepath + '/1.vol_den_sep_Fu18.pdf', bbox_inches='tight')
plt.close(fig)
