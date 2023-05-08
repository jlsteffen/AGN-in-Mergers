;+
; compare expected AGN volume densities with observed densities
; for full manga pair sample, a cut in r_p, and a cut in stellar mass
; 
; need to think about
; * binning issue
; * uncertainty estimates for small number statistics
;-

h = 0.7
minbclass = 2

; setup plot
loadct,0
psfile = 'figs/phi1agn.eps'
nx = 3 
ny = 1
setps,psfile,11*1.1*nx,11*1.3*ny,font='helvetica'
titsize = 1.5
labelsize = 1.2
multiplot,/default
multiplot,[nx,ny],xgap=0.015*1.0,/dox,/doy 

; load galaxy pair catalog
pcat = mrdfits('../7.sample/pair1kpc.fits',1)
pcat.plateifu = strtrim(pcat.plateifu) ; 176
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 149
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 119
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 113
	pcat.sep_arcsec le pcat.IFUrin)  ; 105
pcat = pcat[s]
help,pcat
print,minmax(pcat.sep_kpc)
print,minmax(alog10(pcat.NSA_ELPETRO_MASS/h^2))
; AGN flags 
flg_pagn = pcat.bclass[0] ge minbclass ; primary is AGN
flg_sagn = pcat.bclass[1] ge minbclass and pcat.bclass[0] lt minbclass  ; secondary is AGN
flg_bagn = min(pcat.bclass,dim=1) ge minbclass ; both comps are AGN

; load full target catalog w/ new weights
target = mrdfits('../../sample/target/DR14_newcosm.fits',1)
; match target catalog to pair catalog
match2,pcat.NSA_NSAID,target.NSA_NSAID,sa,sb
ptarg = target[sa]

; load NSA v1_0_1 catalog
nsa = mrdfits('../../sample/matched/nsa.fits',1)
match2,pcat.NSA_NSAID,nsa.NSAID,sa,sb
pnsa = nsa[sa]
; Wake 17: Petrosian magnitudes, with the Petrosian radius rP 
; defined as the major axis of the ellipse where the Petrosian 
; ratio η = 0.2, and with the aperture for the flux defined with 
; a major-axis radius 2x rP
; Elliptical SDSS-style Petrosian radius (r-band, no corrections)
r_ap = 2 * pnsa.ELPETRO_THETA_R ; arcsec
help,where(r_ap lt pcat.sep_arcsec) ; -1 

; compute f_AGN for all pairs using best-fit function
z = pcat.NSA_Z
mass_t = pcat.NSA_ELPETRO_MASS/h^2 ; total mass
mass_p = alog10(mass_t * 10.^pcat.mstar[0]/(10.^pcat.mstar[1]+$
	10.^pcat.mstar[0]))
mass_s = alog10(mass_t * 10.^pcat.mstar[1]/(10.^pcat.mstar[1]+$
	10.^pcat.mstar[0]))
;mass_p = alog10(pcat.NSA_ELPETRO_MASS/h^2)
;mass_s = mass_p + (pcat.mstar[1] - pcat.mstar[0])
; take best-fit AGN model
;save,p,perr,pcerr,bestnorm,dof,filename='pd2agn.idlsav'
;agnmod = p[0] * exp(-0.5d*(yy-p[1])^2/p[2]^2) * (1+xx)^4
restore,'pd2agn.idlsav'
fagn_p =  p[0] * exp(-0.5d*(mass_p-p[1])^2/p[2]^2) * (1+z)^4
fagn_s =  p[0] * exp(-0.5d*(mass_s-p[1])^2/p[2]^2) * (1+z)^4

; Panel 1 
; define X grid for binning
xtit = textoidl('r_p (kpc)')
ytit='log n'+textoidl(' (# per 10^6 Mpc^3)')
ytit2 = 'Excess (dex)'
x = pcat.sep_kpc 
w = ptarg.ESWEIGHT ; W = n*1e6 Mpc^3 => W is the # per 1e6 Mpc^3
;xgri = range(2.5,27.5,6)
;xbin = xgri*0 + 5.0
xgri = [3.5,8.5,13.5,18.5,25.5]
xbin = [5.0,5.0,5.0 ,5.0 ,9]
calc_phi,x,w,xgri,xbin,flg_bagn,flg_pagn,flg_sagn,fagn_p,fagn_s,$
	phi_1agn_exp=phi_1agn_exp,phi_1agn_obs=phi_1agn_obs,$
	phi_bagn_exp=phi_bagn_exp,phi_bagn_obs=phi_bagn_obs
; common arrays
xr = [-2.5,32.5]
y = alog10(w) 
idx = where(flg_pagn or flg_sagn)
blanks = replicate(' ',10)
yr = [-1.8,1.8]
yr2 = [-1.2,1.2]
phi_obs = phi_1agn_obs
phi_exp = phi_1agn_exp
mkexcess,xgri,xbin/2,phi_obs,phi_exp,xtit=xtit,ytit=ytit2,xticklen=0.04,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr2
al_legend,'(a)',/bottom,/left,box=0,chars=1.5,margin=0
mkplot,x,y,idx,xgri,xbin/2,phi_obs,phi_exp,ytit=ytit,xtickname=blanks,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr
xyouts,mean(xr),yr[1]-0.3,'All Pairs',align=0.5,/data,charsize=labelsize
al_legend,['Weights','Observed','Expected'],psym=[9,6,27],$ ;linsize=[0.5,0.5],$
	color=cgcolor(['dark gray','red','royal blue']),$
	/top,/right,background=cgcolor('white'),charsize=1.0
multiplot,/dox,/doy

; Panel 2 - W vs mass
xtit = 'log M (M'+sunsymbol()+')'
ytit = ''
ytit2= ''
rpmax = 10 ; kpc
s = where(pcat.sep_kpc le rpmax)
x = alog10(pcat[s].NSA_ELPETRO_MASS/h^2) 
w = ptarg[s].ESWEIGHT ; W = n*1e6 Mpc^3 => W is the # per 1e6 Mpc^3
print,minmax(x)
xgri2 = [9.5,10.25,10.75,11.25]
xbin2 = [1.0,0.5,0.5,0.5]
calc_phi,x,w,xgri2,xbin2,flg_bagn[s],flg_pagn[s],flg_sagn[s],fagn_p[s],fagn_s[s],$
	phi_1agn_exp=phi_1agn_exp,phi_1agn_obs=phi_1agn_obs,$
	phi_bagn_exp=phi_bagn_exp,phi_bagn_obs=phi_bagn_obs
xr = [8.5,12.0]
y = alog10(w) 
idx = where(flg_pagn[s] or flg_sagn[s]) 
blanks = replicate(' ',10)
yr = [-1.8,1.8]
yr2 = [-1.2,1.2]
phi_obs = phi_1agn_obs
phi_exp = phi_1agn_exp
mkexcess,xgri2,xbin2/2,phi_obs,phi_exp,xtit=xtit,ytit=ytit2,xticklen=0.04,$
	xtickint=1,xminor=5,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr2
al_legend,'(b)',/bottom,/left,box=0,chars=1.5,margin=0
mkplot,x,y,idx,xgri2,xbin2/2,phi_obs,phi_exp,ytit=ytit,xtickname=blanks,$
	xtickint=1,xminor=5,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr
xyouts,mean(xr),yr[1]-0.3,textoidl('r_p < ')+strc(rpmax)+' kpc',align=0.5,/data,charsize=labelsize
multiplot,/dox,/doy

; Panel 3 - W vs. r_p
xtit = textoidl('r_p (kpc)')
logM = alog10(pcat.NSA_ELPETRO_MASS/h^2)
mmin = 10.5
mmax = 12.0
s = where(logM gt mmin and logM lt mmax) 
x = pcat[s].sep_kpc 
w = ptarg[s].ESWEIGHT ; W = n*1e6 Mpc^3 => W is the # per 1e6 Mpc^3
calc_phi,x,w,xgri,xbin,flg_bagn[s],flg_pagn[s],flg_sagn[s],fagn_p[s],fagn_s[s],$
	phi_1agn_exp=phi_1agn_exp,phi_1agn_obs=phi_1agn_obs,$
	phi_bagn_exp=phi_bagn_exp,phi_bagn_obs=phi_bagn_obs
xr = [-2.5,32.5]
y = alog10(w) 
idx = where(flg_pagn[s] or flg_sagn[s]) 
phi_obs = phi_1agn_obs
phi_exp = phi_1agn_exp
mkexcess,xgri,xbin/2,phi_obs,phi_exp,xtit=xtit,ytit=ytit2,xticklen=0.04,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr2
al_legend,'(c)',/bottom,/left,box=0,chars=1.5,margin=0
mkplot,x,y,idx,xgri,xbin/2,phi_obs,phi_exp,ytit=ytit,xtickname=blanks,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr
xyouts,mean(xr),yr[1]-0.3,'log M/M'+sunsymbol()+' > '+strc(mmin),align=0.5,/data,charsize=labelsize
multiplot,/dox,/doy

theend:
multiplot,/default,/reset
device,/close
;spawn,'gv '+psfile+' &'

end

