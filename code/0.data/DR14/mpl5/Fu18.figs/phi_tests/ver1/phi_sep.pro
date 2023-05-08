; compare observed volume densities of AGNs in pairs vs. expected volume
; densities
h = 0.7
minbclass = 2

; setup plot
loadct,0
psfile = 'figs/phi_sep.eps'
nx = 2 
ny = 1
setps,psfile,11*1.1*nx,11*1.3*ny,font='helvetica'
titsize = 1.5
labelsize = 1.2
multiplot,/default
multiplot,[nx,ny],xgap=0.015*1.5,/dox,/doy 

;xtit = 'Pair Separation (kpc)'
;ytit='Volume Density ('+textoidl('log(# per 10^6 Mpc^3)')
xtit = textoidl('r_p (kpc)')
ytit='log n'+textoidl(' (# per 10^6 Mpc^3)')
;ytit2 = textoidl('log n_{obs} / n_{exp}')
ytit2 = 'AGN Excess'

; load galaxy pair catalog
pcat = mrdfits('../7.sample/pair1kpc.fits',1)
pcat.plateifu = strtrim(pcat.plateifu) ; 176
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 149
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 119
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 113
	pcat.sep_arcsec le pcat.IFUrin)  ; 105
pcat = pcat[s]
; AGN flag 
flg_pagn = pcat.bclass[0] ge minbclass ; primary is AGN
flg_sagn = pcat.bclass[1] ge minbclass ; secondary is AGN
flg_bagn = min(pcat.bclass,dim=1) ge minbclass ; both comps are AGN

; load full target catalog w/ new weights
target = mrdfits('../../sample/target/DR14_newcosm.fits',1)
; match target catalog to pair catalog
match2,pcat.NSA_NSAID,target.NSA_NSAID,sa,sb
ptarg = target[sa]

;; load f_AGN vs mass & Rifu
;restore,'figs/pd2d_Rifu_mass.sav'
;fagn = aphi/cphi
;mid = where(cden lt 1)
;fagn[mid] = 0.0
;; interpolate f_AGN map for all pairs
;Rifu = pcat.IFUrin*dangular(pcat.NSA_Z,/kpc)/(180.*3600./!pi) ; kpc
;mass = alog10(pcat.NSA_ELPETRO_MASS/h^2)
;xgri = xx[*,0]
;ygri = yy[0,*]
;xidx = interpol(findgen(n_elements(xgri)),xgri,Rifu)
;yidx = interpol(findgen(n_elements(ygri)),ygri,mass)
;fagn1 = interpolate(fagn,xidx,yidx)

; compute f_AGN for all pairs using best-fit function
z = pcat.NSA_Z
mass_p = alog10(pcat.NSA_ELPETRO_MASS/h^2)
mass_s = mass_p + (pcat.mstar[1] - pcat.mstar[0])
; take best-fit AGN model
;save,p,perr,pcerr,bestnorm,dof,filename='pd2agn.idlsav'
;agnmod = p[0] * exp(-0.5d*(yy-p[1])^2/p[2]^2) * (1+xx)^4
restore,'pd2agn.idlsav'
fagn_p =  p[0] * exp(-0.5d*(mass_p-p[1])^2/p[2]^2) * (1+z)^4
fagn_s =  p[0] * exp(-0.5d*(mass_s-p[1])^2/p[2]^2) * (1+z)^4

; define X grid for binning
x = pcat.sep_kpc 
xgri = [3.5,8.5,16,26]
xbin = [5.0,5.0,10,10]
;xgri = range(2.5,27.5,6)
;xbin = xgri*0+5.0
;xgri = [5.0,15,25]+1
;xbin = xgri*0 + 10
xerr = xbin/2
xnp = n_elements(xgri)
; statistical weights to sum up
w = ptarg.ESWEIGHT ; W = n*1e6 Mpc^3 => W is the # per 1e6 Mpc^3
phi_pagn_exp = fltarr(xnp,3)
phi_sagn_exp = fltarr(xnp,3)
phi_bagn_exp = fltarr(xnp,3)
phi_pair_obs = fltarr(xnp,3)
phi_pagn_obs = fltarr(xnp,3)
phi_sagn_obs = fltarr(xnp,3)
phi_bagn_obs = fltarr(xnp,3)
phi_1agn_exp = fltarr(xnp,3) ; either primary or secondary
phi_1agn_obs = fltarr(xnp,3) ; either primary or secondary

nboot = 1000
for i=0,n_elements(xgri)-1 do begin
	; expected volume densities
	insidebin = x ge xgri[i]-xbin[i]/2 and x lt xgri[i]+xbin[i]/2
	ind_pair = where(insidebin, npair)
	if npair gt 0 then begin
		; lower limit, mean, and upper limit
		ct = npair
		ind = ind_pair
		phi_pair_obs[i,*] = bootstrap_mean(w[ind],nboot=nboot)*ct
		phi_pagn_exp[i,*] = bootstrap_mean(w[ind]*fagn_p[ind],nboot=nboot)*ct
		phi_sagn_exp[i,*] = bootstrap_mean(w[ind]*fagn_s[ind],nboot=nboot)*ct
		phi_bagn_exp[i,*] = bootstrap_mean(w[ind]*fagn_p[ind]*fagn_s[ind],nboot=nboot)*ct
		phi_1agn_exp[i,*] = bootstrap_mean(w[ind]*fagn_p[ind]+$
			w[ind]*fagn_s[ind],nboot=nboot)*ct
	endif
	; observed volume densities
	ind = where(flg_pagn and insidebin)
	phi_pagn_obs[i,*] = calc_phi(w,ind,ind_pair,nboot)
	ind = where(flg_sagn and insidebin)
	phi_sagn_obs[i,*] = calc_phi(w,ind,ind_pair,nboot)
	ind = where(flg_bagn and insidebin)
	phi_bagn_obs[i,*] = calc_phi(w,ind,ind_pair,nboot)
	; either primary or secondary is AGN
	ind1 = where(flg_pagn and insidebin, ct1)
	ind2 = where(flg_sagn and insidebin, ct2)
	if ct1 gt 0 and ct2 gt 0 then ind = [ind1,ind2]
	if ct1 gt 0 and ct2 eq 0 then ind = ind1
	if ct1 eq 0 and ct2 gt 0 then ind = ind2
	phi_1agn_obs[i,*] = calc_phi(w,ind,ind_pair,nboot)
endfor

; common arrays
xr = [-2.5,32.5]
y = alog10(w) 
idxp = where(flg_pagn) ; primary is AGN
idxs = where(flg_sagn) ; secondary is AGN
idxb = where(flg_bagn) ; both are AGN
idx1 = where(flg_pagn or flg_sagn) ; either is AGN

blanks = replicate(' ',10)
; Panel: either as AGN
yr = [-1.8,1.8]
yr2 = [-1.2,1.2]
idx = idx1
phi_obs = phi_1agn_obs
phi_exp = phi_1agn_exp
mkexcess,xgri,xerr,phi_obs,phi_exp,xtit=xtit,ytit=ytit2,xticklen=0.04,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr2
al_legend,'(a)',/bottom,/left,box=0,chars=1.5,margin=0
mkplot,x,y,idx,xgri,xerr,phi_obs,phi_exp,ytit=ytit,xtickname=blanks,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr
xyouts,mean(xr),yr[1]-0.3,'Either as AGN',align=0.5,/data,charsize=labelsize
;al_legend,['Weights','Observed','Expected'],psym=[9,6,27],$ ;linsize=[0.5,0.5],$
;	color=cgcolor(['dark gray','red','royal blue']),$
;	/top,/right,background=cgcolor('white'),charsize=1.0
multiplot,/dox,/doy

; Panel: BAGN
ytit = ''
ytit2 = ''
idx = idxb
phi_obs = phi_bagn_obs
phi_exp = phi_bagn_exp
mkexcess,xgri,xerr,phi_obs,phi_exp,xtit=xtit,ytit=ytit2,xticklen=0.04,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr2
al_legend,'(b)',/bottom,/left,box=0,chars=1.5,margin=0
mkplot,x,y,idx,xgri,xerr,phi_obs,phi_exp,ytit=ytit,xtickname=blanks,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr
xyouts,mean(xr),yr[1]-0.3,'Binary AGN',align=0.6,/data,charsize=labelsize
al_legend,['Weights','Observed','Expected'],psym=[9,6,27],$ ;linsize=[0.5,0.5],$
	color=cgcolor(['dark gray','red','royal blue']),$
	/top,/right,background=cgcolor('white'),charsize=labelsize
multiplot,/dox,/doy


theend:
multiplot,/default,/reset
device,/close
;spawn,'gv '+psfile+' &'


end
