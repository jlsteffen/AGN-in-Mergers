#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 14:47:48 2021

Calculate the stochastic AGN feeding baseline.

Make redshift to stellar mass 2D histogram, like old pair fraction plot

@author: joshua
"""
import os, sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import plotstyle as ps
import functions as f
#______________________________________________________________________________
def frac1d(cand, al, bins):
    '''
    produce 1d histograms and fractions
    '''
    all_hist = np.histogram(al, bins = bins)[0].astype(float)
    cand_hist = np.histogram(cand, bins = bins)[0].astype(float)
    
    bet = f.binormial_er(a=cand_hist, b=all_hist)
    return bet[0], bet[1], bet[2]

def fit_par(cand_hist, all_hist, model, start, stop, spacing):
    '''
    Fit a model pair fraction
    '''
    cand_hist_poisson = np.sqrt(cand_hist) 
    all_hist_poisson = np.sqrt(all_hist)
    
    # get the variance histogram for the pair fraction histogram
    var_hist = np.sqrt(((f.divide(cand_hist_poisson**2, all_hist**2)) + \
                (f.divide(cand_hist**2 * all_hist_poisson**2, all_hist**4))**2))
    var_hist = var_hist * (pair_hist !=0) + 0.0 * (pair_hist == 0)
    # Get the best-fit parameter
    par_range = np.arange(start, stop, spacing)
    chi2_arr = []
    for i in range(len(par_range)):
        model2 = model*par_range[i]/100.0
        chi2 = np.nansum(f.divide((pair_hist-model2)**2, var_hist**2)) / len(np.where(model!=-0.05)[0])
        chi2_arr.append(chi2)
    chi2_arr = np.asarray(chi2_arr)
    return par_range[(chi2_arr == chi2_arr.min())]/100.0, chi2_arr[(chi2_arr == chi2_arr.min())]
#______________________________________________________________________________
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = filepath + '/0.data'

drpall = fits.getdata(datapath+'/0.drpall/drpall_clean_weight.fits')
plifu = drpall['plateifu']
mass = np.log10(drpall['nsa_elpetro_mass']/(0.7**2))
z = drpall['nsa_z']
#______________________________________________________________________________

#Control Selection
ncomp = fits.getdata(filepath+'/1.pair.fits')
ctrl = (ncomp['CTRL']==1)

pair = (ncomp['ngood']==2)

#AGN Selection
agn = ncomp['cagn']==1
pagn = ncomp['eagn']==1

z = drpall['nsa_z']
#______________________________________________________________________________
#Mass and IFU selection for pairs and full sample
all_mass = mass[ctrl]
cand_mass = mass[(agn)]

all_z = z[ctrl]
cand_z = z[(agn)]
#______________________________________________________________________________
# Histogram bin definitions
mbinsize = 0.2
zbinsize = 0.01
mass_bins = np.arange(8.4, 12.2 + mbinsize, mbinsize)
rad_bins = np.arange(-0.01, 0.17 + zbinsize, zbinsize)

mass_width = (mass_bins.max() - mass_bins.min())/(len(mass_bins)-1)
rad_width = (rad_bins.max() - rad_bins.min())/(len(rad_bins)-1)

mass_cent = mass_bins[1::] - mass_width/2.
rad_cent = rad_bins[1::] - rad_width/2.
#______________________________________________________________________________
# 2-D Histogram
cand_hist = np.histogram2d(cand_z, cand_mass, bins=[rad_bins, mass_bins])[0]
all_hist, xedges, yedges = np.histogram2d(all_z, all_mass, bins=[rad_bins, mass_bins])

# Empty regions in the total MaNGA parameter space are assigned the value NaN
min_size = 15
vmax = 0.2
mask = (all_hist<=min_size)
cand_hist[mask], all_hist[mask] = np.nan, np.nan

# 2D-Fraction plot
pair_hist = f.divide(cand_hist, all_hist)
pair_hist[mask] = np.nan

# 1D Masks
maski = (~mask).astype('int')

mc1d = []
for i in range(len(cand_mass)):
    im = np.where(np.abs(mass_cent-cand_mass[i]) == np.abs(mass_cent-cand_mass[i]).min())[0][0]
    ir = np.where(np.abs(rad_cent-cand_z[i]) == np.abs(rad_cent-cand_z[i]).min())[0][0]
    if maski[ir, im] == 1:
        mc1d.append(True)
    else:
        mc1d.append(False)
mc1d = np.array(mc1d)

cand_mass = cand_mass[mc1d]
cand_z = cand_z[mc1d]

ma1d = []
for i in range(len(all_mass)):
    im = np.where(np.abs(mass_cent-all_mass[i]) == np.abs(mass_cent-all_mass[i]).min())[0][0]
    ir = np.where(np.abs(rad_cent-all_z[i]) == np.abs(rad_cent-all_z[i]).min())[0][0]
    if maski[ir, im] == 1:
        ma1d.append(True)
    else:
        ma1d.append(False)
ma1d = np.array(ma1d)

all_mass = all_mass[ma1d]
all_z = all_z[ma1d]
#______________________________________________________________________________                        
# 1-D Histogram fractions
mass_fraction, mass_frac_lower, mass_frac_upper = frac1d(cand_mass, all_mass, mass_bins)
rad_fraction, rad_frac_lower, rad_frac_upper = frac1d(cand_z, all_z, rad_bins)
#______________________________________________________________________________
# Use stellar mass histogram to find b and sigma
from scipy.optimize import curve_fit
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

p0 = [0.12, 10.5, 0.5]
coeff, cov = curve_fit(gauss, mass_cent, mass_fraction, p0=p0)

from lmfit.models import GaussianModel
mod = GaussianModel()

pars = mod.guess(mass_fraction, x=mass_cent)
out = mod.fit(mass_fraction, pars, x=mass_cent)

b = out.params['center'].value
sigma = out.params['sigma'].value

# make model AGN fraction
model = np.empty((0,len(mass_cent)))
for i in range(len(rad_cent)):
    row = np.exp(-0.5*((mass_cent-b)**2 / (sigma)**2)) * (1+rad_cent)**4
    model = np.append(model, [row], axis=0)


fitting = fit_par(cand_hist, all_hist, model, 10, 30, 0.001)
print()
print('fitting pair = ' + str(fitting[0] *100.0) + ' %')
print('chi^2 = '+str(fitting[1]))
print()
model = fitting[0] * model * (pair_hist!=-0.05) + 0 * (pair_hist==-0.05)

model_n = model * all_hist
mass_model_fraction = f.divide(np.nansum(model_n, axis=0), np.nansum(all_hist, axis=0))
rad_model_fraction = f.divide(np.nansum(model_n, axis=1),np.nansum(all_hist, axis=1))

model[np.isnan(pair_hist)] = np.nan


# write a txt file for the model params
t = open("5.model.txt", "w+")
t.write("f = "+str(fitting[0][0])+"\nb = "+str(b)+"\ns = "+str(sigma))
t.close()

#______________________________________________________________________________
### The old 2-D plot ###
# Plot Parameters
fontsize = 20.0

fig = plt.figure(figsize= (12,12))
gs = plt.GridSpec(30, 16, width_ratios = [1]*16, height_ratios = [1]*30)

cmap = plt.cm.coolwarm
#______________________________________________________________________________
# 2-D Histogram
ax1 = plt.subplot(gs[2:20, 0:10])

# make empty bins NaNs
pair_hist[(all_hist<0)] = np.nan

# 2D Histogram show
hist = ax1.imshow(pair_hist.T, origin='lower', extent=[rad_bins.min(), rad_bins.max(),
                mass_bins.min(), mass_bins.max()], aspect='auto', cmap=cmap, vmin=0, vmax=vmax)

ax1.set_ylim((mass_bins.min(), mass_bins.max()))
ax1.set_xlim(rad_bins.min(), rad_bins.max())
ax1.set_ylabel(r'log(Stellar Mass/$M_\odot$)', fontsize=fontsize)
#______________________________________________________________________________
# Mass Histogram
ax2 = plt.subplot(gs[2:20, 10:16], sharey=ax1)
ax2.hist(cand_mass, bins=mass_bins, orientation='horizontal', facecolor=None, 
         weights=np.ones_like(cand_mass)/float(len(cand_mass)), alpha=0.5, label='AGN Dist.', hatch='///', edgecolor='blue', fill=None, linewidth=1.2)
ax2.hist(all_mass, bins=mass_bins, orientation='horizontal', color='grey', 
         weights=np.ones_like(all_mass)/float(len(all_mass)), alpha=0.5, label='All MaNGA', edgecolor='black', linewidth=1.2)
ax2.plot(mass_model_fraction[mass_fraction!=0], mass_cent[mass_fraction!=0], linestyle='--', c='green', label='Model Fraction')

ax2.errorbar(mass_fraction[mass_fraction!=0], mass_cent[mass_fraction!=0], yerr = mass_width/2., 
             xerr = [mass_frac_lower[mass_fraction!=0], mass_frac_upper[mass_fraction!=0]], fmt='s', c='r', 
             label='AGN Fraction', markerfacecolor='white', capsize=3.0, capthick=1.2)

ax2.xaxis.set_label_position('top')
ax2.set_ylim((mass_bins.min(), mass_bins.max()))
ax2.set_xlim(-0.025, 0.35)
ax2.set_xlabel('Fraction', fontsize=fontsize)
legend=ax2.legend(bbox_to_anchor=(0.025, -0.025), loc=2, borderaxespad=0., fontsize=fontsize, 
                  frameon=True, fancybox=False, framealpha=1, edgecolor='k')
legend.get_frame().set_linewidth(2)
#______________________________________________________________________________
# Ifu size Histogram
ax3 = plt.subplot(gs[20:30, 0:10])
ax3.hist(cand_z, bins=rad_bins, facecolor=None, weights=np.ones_like(cand_z)/float(len(cand_z)), alpha=0.5, 
         hatch='///', edgecolor='blue', fill=None, linewidth=1.2)
ax3.hist(all_z, bins=rad_bins, color='grey',weights=np.ones_like(all_z)/float(len(all_z)), alpha=0.5, edgecolor='black', linewidth=1.2)
ax3.plot(rad_cent[rad_fraction!=0], rad_model_fraction[rad_fraction!=0], linestyle='--', c='green')

ax3.errorbar(rad_cent[rad_fraction!=0], rad_fraction[rad_fraction!=0], xerr=rad_width/2., 
             yerr = [rad_frac_lower[rad_fraction!=0], rad_frac_upper[rad_fraction!=0]], 
             marker='s', c='r', fmt='s', markerfacecolor='white', capsize=3.0, capthick=1.2)

ax3.set_xlim(rad_bins.min(), rad_bins.max())
ax3.set_ylim(-0.025, 0.35)
ax3.set_xlabel('Redshift', fontsize=fontsize)
ax3.set_ylabel('Fraction', fontsize=fontsize)
#______________________________________________________________________________
# Model Subplot
inset_axes = inset_axes(ax1,
                        width="40%", # width = 30% of parent_bbox
                        height=2.5, # height : 1 inch
                        loc=4, borderpad=2)
model[all_hist<4] = np.nan
plt.imshow(model.T, cmap=cmap, origin='lower', extent=[rad_bins.min(), rad_bins.max(),
                mass_bins.min(), mass_bins.max()], aspect='auto', vmin=0, vmax=vmax)
plt.title('Model', fontsize=fontsize/1.5)

ps.ticks(inset_axes, xmajor=zbinsize*5, ymajor=mbinsize*4, xminor=zbinsize, yminor=mbinsize)
ps.style(inset_axes, fontsize=fontsize/1.5, labelbottom=True, labelleft=True) 
#______________________________________________________________________________
# Color bar
ax4 = plt.subplot(gs[1:2, 2:8])
cb = plt.colorbar(hist, cax=ax4, orientation="horizontal", ticks = np.linspace(0, vmax, int(vmax*10 + 1)), 
                  boundaries=np.linspace(0, vmax, 1000))
ax4.xaxis.set_ticks_position("top")
ps.cbar_style(cb, fontsize=fontsize)
#______________________________________________________________________________
# Axis parameters
plt.subplots_adjust(wspace=0)
plt.subplots_adjust(hspace=0)

ax4.set_title(' AGN Fraction', pad=30, fontsize=fontsize)

ax1.annotate(r'$\rm{(a)}$', (0, 11.92), fontsize=fontsize)
ax2.annotate(r'$\rm{(b)}$', (0.28, 11.92), fontsize=fontsize)
ax3.annotate(r'$\rm{(c)}$', (0, 0.3), fontsize=fontsize)

ps.ticks(ax1, xmajor=zbinsize*5, ymajor=mbinsize*4, xminor=zbinsize, yminor=mbinsize)
ps.style(ax1, fontsize=fontsize, labelbottom=False, labelleft=True) 
ps.ticks(ax2, xmajor=0.1, ymajor=mbinsize*4, xminor=0.05, yminor=mbinsize)
ps.style(ax2, fontsize=fontsize, labelbottom=False, labelleft=False, labeltop=True, labelright=True) 
ps.ticks(ax3, xmajor=zbinsize*5, ymajor=0.1, xminor=zbinsize, yminor=0.05)
ps.style(ax3, fontsize=fontsize, labelbottom=True, labelleft=True) 

plt.savefig(filepath + '/5.agn_fraction.pdf', bbox_inches='tight')
plt.close('all')
#______________________________________________________________________________
# 1-D Stellar Mass histogram
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(16, 8))
fontsize=30.0

hist = np.histogram(all_mass, bins=mass_bins, weights=np.ones_like(all_mass)/float(len(all_mass)))
ax1.fill_between(mass_cent, hist[0], np.array([0]*len(hist[0])), color='darkgrey', alpha=1.0, step='mid', linewidth=2.0, label='All MaNGA')

hist = np.histogram(cand_mass, bins=mass_bins, weights=np.ones_like(cand_mass)/float(len(cand_mass)))
ax1.fill_between(mass_cent, hist[0], np.array([0]*len(hist[0])), color='deepskyblue', alpha=0.5, step='mid', label='AGN')

ax1.plot(mass_cent[mass_fraction!=0], mass_model_fraction[mass_fraction!=0], linestyle='--', c='k', label='Model Fraction')

ax1.errorbar(mass_cent[mass_fraction!=0], mass_fraction[mass_fraction!=0], xerr = mass_width/2., 
             yerr = [mass_frac_lower[mass_fraction!=0], mass_frac_upper[mass_fraction!=0]], fmt='D', c='k', 
             label='AGN Fraction (This Work)', capsize=3.0, capthick=1.2, markerfacecolor='w')

fu18 = np.loadtxt('/Volumes/Ext_drive/Research/mpl11/20Fu18_fig3_mass_hist.csv', delimiter=',')

fulerr = np.abs(fu18[:,2] - fu18[:,0])
fuuerr = np.abs(fu18[:,3] - fu18[:,0])

fux = np.arange(9.375, 11.625 + 0.25, 0.25)

ax1.errorbar(fux, fu18[:,0], xerr=0.25/2, yerr=[fulerr, fuuerr], c='crimson', fmt='D',
             label='AGN Fraction (Paper I)', capsize=3.0, capthick=1.2, markerfacecolor='w')
#______________________________________________________________________________
# pair histogram

pair_mass = mass[pair]
pagn_mass = mass[pagn]

mbinsize = 0.5
mass_bins = np.arange(8.25, 12.25 + mbinsize, mbinsize)
mass_width = (mass_bins.max() - mass_bins.min())/(len(mass_bins)-1)
mass_cent = mass_bins[1::] - mass_width/2.

hist = np.histogram(pair_mass, bins=mass_bins, weights=np.ones_like(pair_mass)/float(len(pair_mass)))
ax2.fill_between(mass_cent, hist[0], np.array([0]*len(hist[0])), color='darkgrey', alpha=1.0, step='mid', linewidth=2.0, label='All Pairs')

hist = np.histogram(pagn_mass, bins=mass_bins, weights=np.ones_like(pagn_mass)/float(len(pagn_mass)))
ax2.fill_between(mass_cent, hist[0], np.array([0]*len(hist[0])), color='deepskyblue', alpha=0.5, step='mid', label='AGN')

mass_fraction, mass_frac_lower, mass_frac_upper = frac1d(pagn_mass, pair_mass, mass_bins)

ax2.errorbar(mass_cent[mass_fraction!=0], mass_fraction[mass_fraction!=0], xerr = mass_width/2., 
             yerr = [mass_frac_lower[mass_fraction!=0], mass_frac_upper[mass_fraction!=0]], fmt='D', c='k', 
             label='AGN Fraction (This Work)', capsize=3.0, capthick=1.2, markerfacecolor='w')
#______________________________________________________________________________
ax1.set_xlim(8.8, 12.2)
ax1.set_ylim(0, 0.42)
ax1.set_ylabel('Fraction', fontsize=fontsize)
ax1.set_xlabel(r'log(Stellar Mass/M$_\odot$)', fontsize=fontsize)

ax2.set_xlim(8.8, 12.2)
ax2.set_ylim(0, 0.52)
ax2.set_xlabel(r'log(Stellar Mass/M$_\odot$)', fontsize=fontsize)

ax1.set_title('Control Galaxies', fontsize=fontsize)
ax2.set_title('Paired Galaxies', fontsize=fontsize)

ps.ticks(ax1, ymajor=0.1, xmajor=1.0, yminor=0.02, xminor=0.2)
ps.style(ax1, fontsize=fontsize) 
ps.ticks(ax2, ymajor=0.1, xmajor=1.0, yminor=0.02, xminor=0.2)
ps.style(ax2, fontsize=fontsize) 

ps.legend(ax1, fontsize=fontsize, loc=2)
ps.legend(ax2, fontsize=fontsize, loc=2)

plt.savefig(filepath + '/5.agn_hist.pdf', bbox_inches='tight')
plt.close('all')
#______________________________________________________________________________