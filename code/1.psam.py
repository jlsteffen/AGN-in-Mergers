#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 10:08:23 2021

Create a catalog of paired galaxies. Paired galaxies are selected by;

    1. The observation is in the Primary+ or Secondary samples
    2. The MaNGA target has an r-band absolute magnitude that is less than -19
    3. Both pairs are entirely contained within the IFU
    4. The Target galaxies is within the central 10% radius of the IFU
    5. the relative line-of-sight velocity between the galaxies be within 
       Δv < 500 km s−1 .
    6. the mass ratio between the galaxies be within 1:10

@author: joshua
"""
print('\nInitiallizing...\n')

from astropy.io import fits
import numpy as np
import sys, os
import functions as f
from astropy.coordinates import SkyCoord

import warnings
warnings.filterwarnings("ignore")

def dust_corr(loglum, halpha, hbeta, klam=3.5):
    # klam=3.5 is for oiii specifically
    return loglum + (0.79*klam*np.log10(halpha/(2.86*hbeta)))

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = filepath + '/0.data'

drpall = fits.getdata(datapath+'/0.drpall/drpall_clean_weight.fits')
plateifu = drpall['plateifu']
ifurin = drpall['ifurin']
ifurad = drpall['ifurad']
mngtarg1 = drpall['mngtarg1']
mag = drpall['NSA_ELPETRO_absmag'][:,4] + 5*np.log10(0.7)

drpall_all = fits.getdata(datapath+'/3.match/drpall.fits')

spobj = fits.getdata(datapath+'/0.manga_specobj.fits')
plifu = spobj['plateifu']

pars = fits.getdata(datapath+'/3.match/r1kpc.fits')
kin = pars['kinstar'][:,0,0]
sigma = pars['kinstar'][:,1,0]
vgas = pars['vel'][:,0,12] # Ha line vel
ewha = pars['ew'][:,12]  # Ha equivalent width
halpha = pars['flux'][:,0,12]
hbeta = pars['flux'][:,0,6]
mstar = np.log10(np.sum(pars['m_star'][:, 0], axis=1))

bptfil = fits.getdata(datapath+'/3.match/bptw_r1kpc.fits')
bpt = bptfil['bptwhan']
bpt[(bptfil['WHA']<6)&(bpt>=1)&(bpt!=5)] = 0

mag2 = drpall_all['NSA_ELPETRO_absmag'][:,4] + 5*np.log10(0.7)
valid = ((drpall_all['mngtarg1'] & 2**10+2**11+2**12)!=0)&(mag2<-19)

sim = fits.getdata(datapath+'/3.match/simard.fits')
rmag = sim['RGMAG']

lo3 = dust_corr(bptfil['LogLo3'], halpha, hbeta)
#______________________________________________________________________________
# other params
N500 = np.ones(len(drpall)) * -9999
NGOOD = np.ones(len(drpall)) * -9999
CRTL = np.zeros(len(drpall))
CRTLAGN = np.zeros(len(drpall))
SEP = np.ones(len(drpall)) * 9999
SEPKPC = np.ones(len(drpall)) * 9999
DM = np.ones(len(drpall)) * -9999 # subtraction of logged masses, target - companion
SIMDM = np.ones(len(drpall)) * -9999 # mass ratio from R-band Magnitudes
DV = np.ones(len(drpall)) * -9999

SIGMAt = np.ones(len(drpall)) * -9999 # stellar Velocity dispersion
SIGMAc = np.ones(len(drpall)) * -9999 # stellar Velocity dispersion

SEXDM = np.ones(len(drpall)) * -9999 # Mass ratios from SExtractor
PETRODM = np.ones(len(drpall)) * -9999
KRONDM = np.ones(len(drpall)) * -9999

N2HAt = np.ones(len(drpall)) * -9999
O3HBt = np.ones(len(drpall)) * -9999
N2HAc = np.ones(len(drpall)) * -9999
O3HBc = np.ones(len(drpall)) * -9999

WHAt = np.ones(len(drpall)) * -9999
WHAc = np.ones(len(drpall)) * -9999

LO3t = np.ones(len(drpall)) * -9999
LO3c = np.ones(len(drpall)) * -9999

BPTt = np.ones(len(drpall)) * -9999
BPTc = np.ones(len(drpall)) * -9999

MSTARt = np.ones(len(drpall)) * -9999
MSTARc = np.ones(len(drpall)) * -9999

SEMSTARt = np.ones(len(drpall)) * -9999
SEMSTARc = np.ones(len(drpall)) * -9999

EAGN = np.zeros(len(drpall)) # either AGN in pair
SAGN = np.zeros(len(drpall)) # single AGN in pair
DAGN = np.zeros(len(drpall)) # dual AGN in pair

AGNt = np.zeros(len(drpall)) # Is target AGN? 1 == True
AGNc = np.zeros(len(drpall)) # Is Companion AGN? 1 == True

PIND = np.ones(len(drpall)) * -9999 # index of the selected pair

KPC_DA = np.ones(len(drpall)) * -9999
for i in range(len(drpall)):    
    kpc_DA = f.cosmo(drpall['nsa_z'][i])
    KPC_DA[i] = kpc_DA
    
    ix = np.where(plateifu[i]==plifu)[0]
    spt = spobj['sptype'][ix]
    vel = kin[ix]
    vg = vgas[ix]
    
    b = bpt[ix]
    
    ind = spobj['index'][ix]
    
    p = ((spt==1)|(spt==4))&((mngtarg1[i] & 2**10+2**11+2**12)!=0)&(mag[i]<-19)
    NGOOD[i] = p.astype('int').sum()
    
    N2HAt[i] = bptfil['n2ha'][ix][ind==0]
    O3HBt[i] = bptfil['o3hb'][ix][ind==0]
    WHAt[i] = bptfil['WHA'][ix][ind==0]
    BPTt[i] = bptfil['bptwhan'][ix][ind==0]
    LO3t[i] = lo3[ix][ind==0]
    SIGMAt[i] = sigma[ix][ind==0]
    # Pairs
    if p.astype('int').sum() > 1 and p[ind==0] == True:
        c1 = SkyCoord(spobj['ra'][ix][ind==0], spobj['dec'][ix][ind==0], frame='icrs', unit='deg')
        c2 = SkyCoord(spobj['ra'][ix], spobj['dec'][ix], frame='icrs', unit='deg')
        
        sep = c1.separation(c2).arcsecond
        
        dm = mstar[ix][ind==0] - mstar[ix]
        
        vsel = np.abs(vel[ind==0] - vel)
        
        # SExtractor mass ratios
        if os.path.exists(datapath+'/5.segmap/rseg_spfit/'+plateifu[i]+'.fits'):
            fil = fits.open(datapath+'/5.segmap/rseg_spfit/'+plateifu[i]+'.fits')
            nspax = fil[2].data['nspaxel']*0.25
            cat = fits.getdata(datapath+'/5.segmap/rseg/'+plateifu[i]+'_cat.fits')
            semstar = np.log10(np.sum(fil[1].data['m_star'][:,0], axis=1)*nspax)[0]
            sedm = semstar[ind==0] - semstar
            SEMSTARt[i] = semstar[ind==0]
            # Magnitude Ratios
            petrodm =  -0.4*(cat['mag_petro'][ind==0] - cat['mag_petro'])
            krondm = -0.4*(cat['mag_auto'][ind==0] - cat['mag_auto'])
        dm_sel = dm
            
        MSTARt[i] = mstar[ix][ind==0]
        # Is the MaNGA target in the center of the IFU
        ci = SkyCoord(drpall['ifura'][i], drpall['ifudec'][i], frame='icrs', unit='deg')
        co = SkyCoord(drpall['objra'][i], drpall['objdec'][i], frame='icrs', unit='deg')
        dist = ci.separation(co).arcsecond
        
        # pairs
        p2 = p&(vsel<=500)&(np.abs(dm_sel)<=1.0)&(dist<=drpall['ifurad'][i]*0.1)&(sep<=ifurin[i])
        N500[i] = p2.astype('int').sum()
        # AGN samples
        agn = (((spt==1)&(b>=2))|(spt==4))#&(snr[ix]==True)

        # PAirs
        if p2.sum() == 2: 
            sep[ind==0] = 9999
            '''
            Pair Seleciton order
            1. Closest major Merger (dm<=0.5) with an AGN
            2. Closest Major Merger
            3. Closest merger with an AGN
            4. Closest merger
            '''
            if any(p2[(ind!=0)&(np.abs(dm_sel)<=0.5)]):
                sx = np.where(sep==sep[p2&(ind!=0)&(np.abs(dm_sel)<=0.5)].min())[0][0]
            else:
                sx = np.where(sep==sep[p2&(ind!=0)].min())[0][0]
            
            PIND[i] = sx
            
            DV[i] = vsel[sx]
            DM[i] = mstar[ix][ind==0] - mstar[ix][sx]
            if os.path.exists('/Volumes/ssd2t/final_clean/5.segmap/rseg_spfit/'+plateifu[i]+'.fits'):
                SEMSTARc[i] = semstar[sx]
                
                SEXDM[i] = semstar[ind==0] - semstar[sx]
                PETRODM[i] =  -0.4*(cat['mag_petro'][ind==0] - cat['mag_petro'][sx])
                KRONDM[i] = -0.4*(cat['mag_auto'][ind==0] - cat['mag_auto'][sx])
            # Simard mass ratios, only take the ratio if both sources have a 
            # recorded magnitude
            if rmag[ix][ind==0] != 0 and rmag[ix][sx] != 0:
                SIMDM[i] = -0.4*(rmag[ix][ind==0] - rmag[ix][sx])
            
            SEP[i] = sep[sx]
            SEPKPC[i] = sep[sx] * kpc_DA
            
            N2HAc[i] = bptfil['n2ha'][ix][sx]
            O3HBc[i] = bptfil['o3hb'][ix][sx]
            WHAc[i] = bptfil['WHA'][ix][sx]
            BPTc[i]  = bptfil['bptwhan'][ix][sx]
            LO3c[i] = lo3[ix][sx]
            MSTARc[i] = mstar[ix][sx]
            SIGMAc[i] = sigma[ix][sx]
            if (p2[sx]) & (agn[sx])|(agn[ind==0]):
                EAGN[i] = 1
                # Which comps are AGN?
                if agn[ind==0]:
                    AGNt[i] = 1
                if (p2[sx]) & (agn[sx]):
                    AGNc[i] = 1
                
            if (p2[sx]) & ((agn[sx]&~agn[ind==0])|(~agn[sx]&agn[ind==0])):
                SAGN[i] = 1
            if (p2[sx]) & (agn[sx])&(agn[ind==0]):
                DAGN[i] = 1
    
    # Controls
    if p.astype('int').sum() == 1 and p[ind==0] == True:
        CRTL[i] = 1
    cagn = (((spt==1)&(b>=2))|(spt==4))
    if p.astype('int').sum() == 1 and p[ind==0] == True and cagn[ind==0]:
        CRTLAGN[i] = 1
    
    f.update_progress((i+1.0)/np.float64(len(drpall)))
#______________________________________________________________________________
cols = []
cols.append(fits.Column(name='plateifu', format='11A', array=plateifu))
cols.append(fits.Column(name='NGOOD', format='1I', array=NGOOD))
cols.append(fits.Column(name='npair', format='1I', array=N500))
cols.append(fits.Column(name='CTRL', format='1I', array=CRTL))
cols.append(fits.Column(name='sep', format='1E', array=SEP))
cols.append(fits.Column(name='sepkpc', format='1E', array=SEPKPC))
cols.append(fits.Column(name='dlogm', format='1E', array=DM))
cols.append(fits.Column(name='rmag_dm', format='1E', array=SIMDM))
cols.append(fits.Column(name='SExDM', format='1E', array=SEXDM))
cols.append(fits.Column(name='Petrodrmag', format='1E', array=PETRODM))
cols.append(fits.Column(name='Krondrmag', format='1E', array=KRONDM))
cols.append(fits.Column(name='dv', format='1E', array=DV))

cols.append(fits.Column(name='EAGN', format='1I', array=EAGN))
cols.append(fits.Column(name='SAGN', format='1I', array=SAGN))
cols.append(fits.Column(name='DAGN', format='1I', array=DAGN))
cols.append(fits.Column(name='PINDEX', format='1I', array=PIND))
cols.append(fits.Column(name='CAGN', format='1I', array=CRTLAGN))

cols.append(fits.Column(name='AGNt', format='1I', array=AGNt))
cols.append(fits.Column(name='AGNc', format='1I', array=AGNc))

cols.append(fits.Column(name='N2HAt', format='1E', array=N2HAt))
cols.append(fits.Column(name='N2HAc', format='1E', array=N2HAc))
cols.append(fits.Column(name='O3HBt', format='1E', array=O3HBt))
cols.append(fits.Column(name='O3HBc', format='1E', array=O3HBc))

cols.append(fits.Column(name='WHAt', format='1E', array=WHAt))
cols.append(fits.Column(name='WHAc', format='1E', array=WHAc))

cols.append(fits.Column(name='BPTt', format='1I', array=BPTt))
cols.append(fits.Column(name='BPTc', format='1I', array=BPTc))

cols.append(fits.Column(name='LogLo3t', format='1E', array=LO3t))
cols.append(fits.Column(name='LogLo3c', format='1E', array=LO3c))

cols.append(fits.Column(name='Sigmat', format='1E', array=SIGMAt))
cols.append(fits.Column(name='Sigmac', format='1E', array=SIGMAc))

cols.append(fits.Column(name='mstart', format='1E', array=MSTARt))
cols.append(fits.Column(name='mstarc', format='1E', array=MSTARc))

cols.append(fits.Column(name='semstart', format='1E', array=SEMSTARt))
cols.append(fits.Column(name='semstarc', format='1E', array=SEMSTARc))

cols = fits.ColDefs(cols)
tbhdu = fits.BinTableHDU.from_columns(cols)
tbhdu.writeto(filepath+'/1.pair.fits', overwrite=True)
#______________________________________________________________________________
# Make a txt file with sample size stats

txt = 'SpecObj Subsamples \n'
txt+= 'MaNGA Target Catalog size = '+str(len(drpall))+'\n'
txt+= 'MaNGA Primary+ Sample  = '+str(len(drpall[((mngtarg1 & 2**10+2**11+2**12)!=0)]))+'\n'
txt+= 'MaNGA Primary+ Sample + Lum Func = '+str(len(drpall[((mngtarg1 & 2**10+2**11+2**12)!=0)&(mag<-19)]))+'\n'
txt+= 'Spectrocopically Confirmed Catalog Size = '+str(len(spobj))+'\n'
txt+= 'Valid spectroscopic objects (Primary+ + Lum Func) = '+str(len(spobj[valid]))+'\n'
txt+= 'Good Objects = '+str(len(spobj['RA'][(spobj['sptype']==1)&valid]))+'\n'
txt+= 'BLAGN Objects = '+str(len(spobj['RA'][(spobj['sptype']==4)&valid]))+'\n'
txt+= 'Borderline Objects = '+str(len(spobj['RA'][(spobj['sptype']==0)&valid]))+'\n'
txt+= 'Stellar Objects = '+str(len(spobj['RA'][(spobj['sptype']==-1)&valid]))+'\n'
txt+= 'z_off Objects = '+str(len(spobj['RA'][(spobj['sptype']==-2)&valid]))+'\n'
txt+= 'QSO Objects = '+str(len(spobj['RA'][(spobj['sptype']==-3)&valid]))+'\n'
txt+= 'Low SNR Objects = '+str(len(spobj['RA'][(spobj['sptype']==-4)&valid]))+'\n'
txt+= 'Other Objects = '+str(len(spobj['RA'][(spobj['sptype']==-9)&valid]))+'\n\n'



txt+= 'AGN in Whole Sample \n'
txt+= 'Ambiguous in Pairs = '+str(((spobj['sptype']==1)&(bpt<0)&valid).sum())+'\n'
txt+= 'Retired in Pairs = '+str(((spobj['sptype']==1)&(bpt==0)&valid).sum())+'\n'
txt+= 'SF in Pairs = '+str(((spobj['sptype']==1)&(bpt==1)&valid).sum())+'\n'
txt+= 'COMP in Pairs = '+str(((spobj['sptype']==1)&(bpt==2)&valid).sum())+'\n'
txt+= 'LINER in Pairs = '+str(((spobj['sptype']==1)&(bpt==3)&valid).sum())+'\n'
txt+= 'Seyfert in Pairs = '+str(((spobj['sptype']==1)&(bpt==4)&valid).sum())+'\n'
txt+= 'BLAGN in Pairs = '+str(((spobj['sptype']==4)&(bpt==5)&valid).sum())+'\n\n'

txt+= 'Pair Candidate Sets = '+str(len(NGOOD[NGOOD>1]))+'\n'
txt+= 'Pair Candidate Objects = '+str(NGOOD[NGOOD>1].sum())+'\n'
txt+= 'Confirmed Pair Sets = '+str(len(N500[N500>1]))+'\n'
txt+= 'Confirmed Pair Objects = '+str(N500[N500>1].sum())+'\n'
txt+= 'Pairs = '+str(len(N500[N500==2]))+'\n'
txt+= 'Triplets = '+str(len(N500[N500==3]))+'\n'
txt+= 'Quadruplets = '+str(len(N500[N500==4]))+'\n'
txt+= 'Quintuplets = '+str(len(N500[N500==5]))+'\n\n'

txt+= 'Ambiguous in Pairs = '+str(((N500==2)&(BPTt<0)).sum()+((N500==2)&(BPTc<0)).sum())+'\n'
txt+= 'Retired in Pairs = '+str(((N500==2)&(BPTt==0)).sum()+((N500==2)&(BPTc==0)).sum())+'\n'
txt+= 'SF in Pairs = '+str(((N500==2)&(BPTt==1)).sum()+((N500==2)&(BPTc==1)).sum())+'\n'
txt+= 'COMP in Pairs = '+str(((N500==2)&(BPTt==2)).sum()+((N500==2)&(BPTc==2)).sum())+'\n'
txt+= 'LINER in Pairs = '+str(((N500==2)&(BPTt==3)).sum()+((N500==2)&(BPTc==3)).sum())+'\n'
txt+= 'Seyfert in Pairs = '+str(((N500==2)&(BPTt==4)).sum()+((N500==2)&(BPTc==4)).sum())+'\n'
txt+= 'BLAGN in Pairs = '+str(((N500==2)&(BPTt==5)).sum()+((N500==2)&(BPTc==5)).sum())+'\n\n'

txt+= 'Controls = '+str(CRTL.sum())+'\n'
txt+= 'Ambiguous Controls = '+str(CRTL[BPTt==-9].sum())+'\n'
txt+= 'Retired Controls = '+str(CRTL[BPTt==0].sum())+'\n'
txt+= 'SF Controls = '+str(CRTL[BPTt==1].sum())+'\n'
txt+= 'Comp Controls = '+str(CRTL[BPTt==2].sum())+'\n'
txt+= 'LINER Controls = '+str(CRTL[BPTt==3].sum())+'\n'
txt+= 'Seyfert Controls = '+str(CRTL[BPTt==4].sum())+'\n'
txt+= 'BLAGN Controls = '+str(CRTL[BPTt==5].sum())+'\n\n'

txt+= 'Single AGN = '+str(SAGN.sum())+'\n'
txt+= 'Dual AGN = '+str(DAGN.sum())+'\n'
txt+= 'Single+Dual AGN = '+str(EAGN.sum())+'\n'

with open('1.stat.txt', 'w') as f:
    f.write(txt)