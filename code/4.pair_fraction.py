#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 15:16:54 2018

@author: joshua
"""
import os, sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import plotstyle as ps
import warnings
import functions as f
warnings.filterwarnings("ignore")
#______________________________________________________________________________
from scipy.special import betainc
def binormial_er(a, b, cl='cl'):
    frac = np.empty((a.shape))
    lower = np.empty((a.shape))
    upper = np.empty((a.shape))
    for i in range(len(a)):
        k = np.float64(a[i])
        n = np.float64(b[i])
        kn = np.divide(k, n, out=np.zeros(n.shape), where=n!=0)
        if cl == 'cl':
            cl = 0.683 # 1 sigma
        z = np.arange(10**4)*10**-4
        bet = betainc(k+1, n-k+1, z)
        il = np.where(abs(bet - (1.-cl)/2.) == min(abs(bet - (1.-cl)/2.)))[0]
        iu = np.where(abs(bet - (1. - (1-cl)/2.)) == min(abs(bet - (1. - (1.-cl)/2.))))[0]
        p_lower=z[il]
        p_upper=z[iu]
        el = kn - p_lower
        eu = p_upper - kn
        frac[i] = kn
        lower[i] = el
        upper[i] = eu
    return frac, lower, upper

def frac1d(cand, al, bins):
    '''
    produce 1d histograms and fractions
    '''
    all_hist = np.histogram(al, bins = bins)[0].astype(float)
    cand_hist = np.histogram(cand, bins = bins)[0].astype(float)
    
    bet = binormial_er(a=cand_hist, b=all_hist)
    return bet[0], bet[1], bet[2]

def fit_par(cand_hist, all_hist, model, start, stop, spacing):
    '''
    Fit a model pair fraction
    '''
    cand_hist_poisson = np.sqrt(cand_hist) 
    all_hist_poisson = np.sqrt(all_hist)
    
    # get the variance histogram for the pair fraction histogram
    var_hist = np.sqrt(((np.divide(cand_hist_poisson**2, all_hist**2, out=np.zeros(all_hist_poisson.shape), 
                                        where=all_hist!=0)) + \
                (np.divide(cand_hist**2 * all_hist_poisson**2, all_hist**4, out=np.zeros(cand_hist_poisson.shape), 
                                       where=cand_hist!=0))**2))
    var_hist = var_hist * (pair_hist !=0) + 0.0 * (pair_hist == 0)
    # Get the best-fit parameter
    par_range = np.arange(start, stop, spacing)
    chi2_arr = []
    for i in range(len(par_range)):
        model2 = model*par_range[i]/100.0
        chi2 = np.nansum(np.divide((pair_hist-model2)**2, var_hist**2, out=np.zeros(var_hist.shape), 
                                where=var_hist!=0)) / len(np.where(model!=-0.05)[0])
        chi2_arr.append(chi2)
    chi2_arr = np.asarray(chi2_arr)
    return par_range[(chi2_arr == chi2_arr.min())]/100.0, chi2_arr[(chi2_arr == chi2_arr.min())]
#______________________________________________________________________________
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = filepath + '/0.data'

drpall = fits.getdata(datapath+'/0.drpall/drpall_clean.fits')
plifu = drpall['plateifu']
mass = np.log10(drpall['nsa_elpetro_mass']/(0.7**2))
mngtarg1 = drpall['mngtarg1']
mngtarg3 = drpall['mngtarg3']
z = drpall['nsa_z']
#______________________________________________________________________________
'''
Pair Selection
'''
ncomp = fits.getdata(filepath + '/1.pair.fits')
pairs = (ncomp['npair']==2)
ctrl = ((mngtarg1 & 2**10+2**11+2**12)!=0)

kpc_DA = []
for i in z:
    kpc_DA.append(f.cosmo(i))
kpc_DA = np.array(kpc_DA)

ifurin = drpall['ifurin']
ifu_rad = ifurin * kpc_DA
#______________________________________________________________________________
'''
Mass and IFU selection for pairs and full sample
'''
all_mass = mass[ctrl]
cand_mass = mass[(pairs)]

all_rad = ifu_rad[ctrl]
cand_rad = ifu_rad[(pairs)]
#______________________________________________________________________________
# Histogram bin definitions
mbinsize = 0.175
rbinsize = 1.75
mass_bins = np.arange(8.225, 12.25 + mbinsize, mbinsize)
rad_bins = np.arange(-1.75, 42.0 + rbinsize, rbinsize)

mass_width = (mass_bins.max() - mass_bins.min())/(len(mass_bins)-1)
rad_width = (rad_bins.max() - rad_bins.min())/(len(rad_bins)-1)

mass_cent = mass_bins[1::] - mass_width/2.
rad_cent = rad_bins[1::] - rad_width/2.
#______________________________________________________________________________
# 2-D Histogram
cand_hist = np.histogram2d(cand_rad, cand_mass, bins=[rad_bins, mass_bins])[0]
all_hist, xedges, yedges = np.histogram2d(all_rad, all_mass, bins=[rad_bins, mass_bins])

pair_hist = np.divide(cand_hist, all_hist, out=np.zeros(cand_hist.shape), where=all_hist!=0)

# Empty regions in the total MaNGA parameter space are assigned the value -0.05
pair_hist = pair_hist * (all_hist != 0) + -0.05 * (all_hist == 0) 
#______________________________________________________________________________                        
# 1-D Histogram fractions
mass_fraction, mass_frac_lower, mass_frac_upper = frac1d(cand_mass, all_mass, mass_bins)
rad_fraction, rad_frac_lower, rad_frac_upper = frac1d(cand_rad, all_rad, rad_bins)
#______________________________________________________________________________
# make model pair fraction
model = np.empty((0,len(mass_cent)))
for i in range(len(rad_cent)):
    row = (rad_cent[i]/30.0)*(10**((mass_cent-11.0)))
    model = np.append(model, [row], axis=0)

fitting = fit_par(cand_hist, all_hist, model, 10, 20, 0.001)
print()
print('fitting pair = ' + str(fitting[0] *100.0) + ' %')
print('chi^2 = '+str(fitting[1]))
print()
model = fitting[0] * model * (pair_hist!=-0.05) + 0 * (pair_hist==-0.05)

model_n = model * all_hist
mass_model_fraction = np.divide(np.sum(model_n, axis=0),np.sum(all_hist, axis=0), 
                                out=np.zeros(np.sum(model_n, axis=0).shape), where=np.sum(all_hist, axis=0)!=0)
rad_model_fraction = np.divide(np.sum(model_n, axis=1),np.sum(all_hist, axis=1), 
                               out=np.zeros(np.sum(model_n, axis=1).shape), where=np.sum(all_hist, axis=1)!=0)

model = model * (pair_hist!=-0.05) + -0.05 * (pair_hist==-0.05)
#______________________________________________________________________________
# Cummulative pair fraction
def mod(par, r, m):
    return par * (r/30) * (10**m/10**11)
ew = drpall['esweight']

rmax = 20/0.7

Fpair = np.sum(ew[pairs]*mod(fitting[0], rmax, mass[pairs]))/np.sum(ew[pairs])
print('Fpair = '+str(Fpair*100)+'%')
#______________________________________________________________________________
# Plot Parameters
fontsize = 20.0

fig = plt.figure(figsize= (12,12))
gs = plt.GridSpec(30, 16, width_ratios = [1]*16, height_ratios = [1]*30)

cmap = plt.cm.coolwarm
#______________________________________________________________________________
# 2-D Histogram
ax1 = plt.subplot(gs[2:20, 0:10])

pair_hist[(all_hist<4)] = np.nan

hist = ax1.imshow(pair_hist.T, origin='lower', extent=[rad_bins.min(), rad_bins.max(),mass_bins.min(), mass_bins.max()], aspect='auto', cmap=cmap, vmin=0, vmax=0.5)

ax1.set_ylim((mass_bins.min(), mass_bins.max()))
ax1.set_xlim(rad_bins.min(), rad_bins.max())
ax1.set_ylabel(r'log(Stellar Mass/$M_\odot$)', fontsize=fontsize)
#______________________________________________________________________________
# Mass Histogram
ax2 = plt.subplot(gs[2:20, 10:16], sharey=ax1)
ax2.hist(cand_mass, bins=mass_bins, orientation='horizontal', facecolor=None, 
         weights=np.ones_like(cand_mass)/float(len(cand_mass)), alpha=0.5, label='Pairs', hatch='///', edgecolor='blue', fill=None, linewidth=1.2)
ax2.hist(all_mass, bins=mass_bins, orientation='horizontal', color='grey', 
         weights=np.ones_like(all_mass)/float(len(all_mass)), alpha=0.5, label='All MaNGA', edgecolor='black', linewidth=1.2)
ax2.plot(mass_model_fraction[mass_fraction!=0], mass_cent[mass_fraction!=0], linestyle='--', c='green', label='Model Fraction')

ax2.errorbar(mass_fraction[mass_fraction!=0], mass_cent[mass_fraction!=0], yerr = mass_width/2., 
             xerr = [mass_frac_lower[mass_fraction!=0], mass_frac_upper[mass_fraction!=0]], fmt='s', c='r', 
             label='Pair Fraction', markerfacecolor='white', capsize=3.0, capthick=1.2)

ax2.xaxis.set_label_position('top')
ax2.set_ylim((mass_bins.min(), mass_bins.max()))
ax2.set_xlim(-0.025, 0.45)
ax2.set_xlabel('Fraction', fontsize=fontsize)
legend=ax2.legend(bbox_to_anchor=(0.025, -0.025), loc=2, borderaxespad=0., fontsize=fontsize, 
                  frameon=True, fancybox=False, framealpha=1, edgecolor='k')
legend.get_frame().set_linewidth(2)
#______________________________________________________________________________
# Ifu size Histogram
ax3 = plt.subplot(gs[20:30, 0:10])
ax3.hist(cand_rad, bins=rad_bins, facecolor=None, weights=np.ones_like(cand_rad)/float(len(cand_rad)), alpha=0.5, 
         hatch='///', edgecolor='blue', fill=None, linewidth=1.2)
ax3.hist(all_rad, bins=rad_bins, color='grey',weights=np.ones_like(all_rad)/float(len(all_rad)), alpha=0.5, edgecolor='black', linewidth=1.2)
ax3.plot(rad_cent[rad_fraction!=0], rad_model_fraction[rad_fraction!=0], linestyle='--', c='green')

ax3.errorbar(rad_cent[rad_fraction!=0], rad_fraction[rad_fraction!=0], xerr=rad_width/2., 
             yerr = [rad_frac_lower[rad_fraction!=0], rad_frac_upper[rad_fraction!=0]], 
             marker='s', c='r', fmt='s', markerfacecolor='white', capsize=3.0, capthick=1.2)

ax3.set_xlim(rad_bins.min(), rad_bins.max())
ax3.set_ylim(-0.025, 0.65)
ax3.set_xlabel('IFU Radius (kpc)', fontsize=fontsize)
ax3.set_ylabel('Fraction', fontsize=fontsize)
#______________________________________________________________________________
# Model Subplot
inset_axes = inset_axes(ax1,
                        width="40%", # width = 30% of parent_bbox
                        height=2.5, # height : 1 inch
                        loc=4, borderpad=2)
model[all_hist<4] = np.nan
plt.imshow(model.T, cmap=cmap, origin='lower', extent=[rad_bins.min(), rad_bins.max(),mass_bins.min(), mass_bins.max()], aspect='auto', vmin=0, vmax=0.5)
plt.title('Model', fontsize=fontsize/1.5)

ps.ticks(inset_axes, xmajor=rbinsize*5, ymajor=mbinsize*5, xminor=rbinsize, yminor=mbinsize)
ps.style(inset_axes, fontsize=fontsize/1.5, labelbottom=True, labelleft=True) 
#______________________________________________________________________________
# Color bar
ax4 = plt.subplot(gs[1:2, 2:8])
cb = plt.colorbar(hist, cax=ax4, orientation="horizontal", ticks = np.linspace(0, 0.5, 6), boundaries=np.linspace(0, 0.5, 1000))
ax4.xaxis.set_ticks_position("top")
ps.cbar_style(cb, fontsize=fontsize)
#______________________________________________________________________________
# Axis parameters
plt.subplots_adjust(wspace=0)
plt.subplots_adjust(hspace=0)

ax4.set_title(' Pair Fraction', pad=30, fontsize=fontsize)

ax1.annotate(r'$\rm{(a)}$', (0, 12), fontsize=fontsize)
ax2.annotate(r'$\rm{(b)}$', (0.38, 12), fontsize=fontsize)
ax3.annotate(r'$\rm{(c)}$', (0, 0.53), fontsize=fontsize)

ps.ticks(ax1, xmajor=rbinsize*5, ymajor=mbinsize*5, xminor=rbinsize, yminor=mbinsize)
ps.style(ax1, fontsize=fontsize, labelbottom=False, labelleft=True) 
ps.ticks(ax2, xmajor=0.2, ymajor=mbinsize*5, xminor=0.05, yminor=mbinsize)
ps.style(ax2, fontsize=fontsize, labelbottom=False, labelleft=False, labeltop=True, labelright=True) 
ps.ticks(ax3, xmajor=rbinsize*5, ymajor=0.2, xminor=rbinsize, yminor=0.05)
ps.style(ax3, fontsize=fontsize, labelbottom=True, labelleft=True) 

plt.savefig(filepath + '/19pair_fraction.pdf', bbox_inches='tight')
plt.close('all')