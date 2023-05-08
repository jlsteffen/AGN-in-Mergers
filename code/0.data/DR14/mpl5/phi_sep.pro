;+
; compare expected AGN volume densities with observed densities
; for full manga pair sample, a cut in r_p, and a cut in stellar mass
; 
; This version does not double count binary AGNs
;
; To run this, first define the functions by
;	.r phifun.pro
;-

h = 0.7
minbclass = 2

; setup plot
loadct,0
psfile = 'figs/phi_sep.eps'
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

;s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 149
;	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 142 (-7)
;	pcat.sep_arcsec le pcat.IFUrin and $ ; 131 (-11)
;	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0) ; 105 (-26)
	
; load full target catalog w/ new weights
target = mrdfits('../../sample/target/DR14_newcosm.fits',1)
; match target catalog to pair catalog
match2,pcat.NSA_NSAID,target.NSA_NSAID,sa,sb
ptarg = target[sa]

; compute f_AGN for all pairs using best-fit function
z = pcat.NSA_Z
mass_t = pcat.NSA_ELPETRO_MASS/h^2 ; total mass
mass_p = alog10(mass_t * 10.^pcat.mstar[0]/$
	(10.^pcat.mstar[0]+10.^pcat.mstar[1]))
mass_s = alog10(mass_t * 10.^pcat.mstar[1]/$
	(10.^pcat.mstar[0]+10.^pcat.mstar[1]))
; take best-fit AGN model
;save,p,perr,pcerr,bestnorm,dof,filename='pd2agn.idlsav'
;agnmod = p[0] * exp(-0.5d*(yy-p[1])^2/p[2]^2) * (1+xx)^4
restore,'pd2agn.idlsav'
fagn_p =  p[0] * exp(-0.5d*(mass_p-p[1])^2/p[2]^2) * (1+z)^4
fagn_s =  p[0] * exp(-0.5d*(mass_s-p[1])^2/p[2]^2) * (1+z)^4

; AGN flag 
flg_eagn = pcat.bclass[0] ge minbclass or pcat.bclass[1] ge minbclass 
flg_sagn = (pcat.bclass[0] ge minbclass or pcat.bclass[1] ge minbclass) and $
	   min(pcat.bclass,dim=1) lt minbclass
flg_bagn = min(pcat.bclass,dim=1) ge minbclass 
help,where(flg_eagn),where(flg_sagn),where(flg_bagn)

; for plot
xtit = textoidl('r_p (kpc)')
ytit='log n'+textoidl(' (# per 10^6 Mpc^3)')
ytit2 = 'Excess (dex)'
blanks = replicate(' ',10)
xr = [-2.5,32.5]

; raw data
x = pcat.sep_kpc 
w = ptarg.ESWEIGHT ; W = n*1e6 Mpc^3 => W is the # per 1e6 Mpc^3
y = alog10(w) 

; Panel 1 - Either AGN
xgri = [2.5,7.5,12.5,17.5,25]
xbin = [5.0,5.0, 5.0, 5.0,10]
xi = xbin*0
calc_phi,x,w,xgri,xbin,xi,flg_eagn,flg_sagn,flg_bagn,fagn_p,fagn_s,$
	phi_eagn_exp=phi_eagn_exp,phi_eagn_obs=phi_eagn_obs,$
	phi_sagn_exp=phi_sagn_exp,phi_sagn_obs=phi_sagn_obs,$
	phi_bagn_exp=phi_bagn_exp,phi_bagn_obs=phi_bagn_obs
; plotting
idx = where(flg_eagn)
yr = [-1.5,1.8]
yr2 = [-1.2,1.2]
phi_obs = phi_eagn_obs
phi_exp = phi_eagn_exp
mkexcess,xgri,xbin/2,phi_obs,phi_exp,xtit=xtit,ytit=ytit2,xticklen=0.04,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr2
al_legend,'(a)',/bottom,/left,box=0,chars=1.5,margin=0
mkplot,x,y,idx,xgri,xbin/2,phi_obs,phi_exp,ytit=ytit,xtickname=blanks,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr
xyouts,mean(xr),yr[1]-0.3,'Either AGNs',align=0.5,/data,charsize=labelsize
al_legend,['Weights','Observed','Expected'],psym=[9,6,27],$ ;linsize=[0.5,0.5],$
	color=cgcolor(['dark gray','red','royal blue']),$
	/top,/right,background=cgcolor('white'),charsize=1.0
multiplot,/dox,/doy
ytit =''
ytit2 = ''

; Panel 2 - Single AGN
;xgri = [2.5,7.5,12.5,17.5,25]
;xbin = [5.0,5.0, 5.0, 5.0,10]
xgri = [5.,20]
xbin = [10.,20]
xi = xbin*0
calc_phi,x,w,xgri,xbin,xi,flg_eagn,flg_sagn,flg_bagn,fagn_p,fagn_s,$
	phi_eagn_exp=phi_eagn_exp,phi_eagn_obs=phi_eagn_obs,$
	phi_sagn_exp=phi_sagn_exp,phi_sagn_obs=phi_sagn_obs,$
	phi_bagn_exp=phi_bagn_exp,phi_bagn_obs=phi_bagn_obs
; plotting
idx = where(flg_sagn) 
phi_obs = phi_sagn_obs
phi_exp = phi_sagn_exp
mkexcess,xgri,xbin/2,phi_obs,phi_exp,xtit=xtit,ytit=ytit2,xticklen=0.04,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr2
al_legend,'(b)',/bottom,/left,box=0,chars=1.5,margin=0
mkplot,x,y,idx,xgri,xbin/2,phi_obs,phi_exp,ytit=ytit,xtickname=blanks,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr
xyouts,mean(xr),yr[1]-0.3,'Single AGNs',align=0.5,/data,charsize=labelsize
multiplot,/dox,/doy

; Panel 3 - Binary AGN
xgri = [5.,20]
xbin = [10.,20]
xi = xbin*0
;xi = [0.40,0.15]
;xi = [0.40*1.4,0.15*1.5]
;xi = [0.40/1.4,0.15/2.5]
calc_phi,x,w,xgri,xbin,xi,flg_eagn,flg_sagn,flg_bagn,fagn_p,fagn_s,$
	phi_eagn_exp=phi_eagn_exp,phi_eagn_obs=phi_eagn_obs,$
	phi_sagn_exp=phi_sagn_exp,phi_sagn_obs=phi_sagn_obs,$
	phi_bagn_exp=phi_bagn_exp,phi_bagn_obs=phi_bagn_obs
; plotting
idx = where(flg_bagn) 
phi_obs = phi_bagn_obs
phi_exp = phi_bagn_exp
mkexcess,xgri,xbin/2,phi_obs,phi_exp,xtit=xtit,ytit=ytit2,xticklen=0.04,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr2
al_legend,'(c)',/bottom,/left,box=0,chars=1.5,margin=0
mkplot,x,y,idx,xgri,xbin/2,phi_obs,phi_exp,ytit=ytit,xtickname=blanks,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr
xyouts,mean(xr),yr[1]-0.3,'Binary AGNs',align=0.5,/data,charsize=labelsize
multiplot,/dox,/doy

theend:
multiplot,/default,/reset
device,/close
;spawn,'gv '+psfile+' &'


end
