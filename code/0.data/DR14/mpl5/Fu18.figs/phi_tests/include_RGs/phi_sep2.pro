; this version investigates whether the different results between
; Ellison+11 and ours is due to different AGN selection methods

pro mkplot,xx,yy,idx,xgri,xerr,phi_obs,phi_exp, _REF_EXTRA=extra
	; plot frame
	plot,[0,1],[0,1],/nodata,_extra=extra ; /xs,/ys,xr=xr,yr=yr,xtit=xtit,ytit=ytit
	syms = 0.8
	; pair sample
	plots,xx,yy,psym=cgsymcat(9),syms=syms,color=cgcolor('gray'),noclip=0
	; subsample
	plots,xx[idx],yy[idx],psym=cgsymcat(16),syms=syms,noclip=0,color=cgcolor('gray')
	
	; total Phi for expected AGNs in close pairs
	y = alog10(phi_exp[*,1])
	ylerr = (phi_exp[*,1]-phi_exp[*,0])/phi_exp[*,1]/alog(10)
	yuerr = (phi_exp[*,2]-phi_exp[*,1])/phi_exp[*,1]/alog(10)
	; add 0.04 dex uncertainty in f_agn model 
	moderr = 0.04
	ylerr = sqrt(ylerr^2+moderr^2)
	;oploterror,xgri,y,xerr,ylerr,/lobar,/nohat,psym=3,color=cgcolor('royal blue')
	;oploterror,xgri,y,xerr,yuerr,/hibar,/nohat,psym=3,color=cgcolor('royal blue')
	; draw boxes
	for i=0,n_elements(xgri)-1 do $
	tvbox,[xerr[i]*2,yuerr[i]+ylerr[i]],xgri[i],(y[i]-ylerr[i]+y[i]+yuerr[i])/2,$
		/data,color=cgcolor('royal blue'),noclip=0,lines=0,$
		/fill,/line_fill,orient=45,spacing=0.2

	; total Phi for observed AGNs in close pairs
	y = alog10(phi_obs[*,1])
	; d log(y) = 1/ln(10) dy/y
	ylerr = (phi_obs[*,1]-phi_obs[*,0])/phi_obs[*,1]/alog(10)
	yuerr = (phi_obs[*,2]-phi_obs[*,1])/phi_obs[*,1]/alog(10)
	oploterror,xgri,y,xerr,ylerr,/lobar,psym=3,color=cgcolor('red')
	oploterror,xgri,y,xerr,yuerr,/hibar,psym=3,color=cgcolor('red')
	plots,xgri,y,psym=cgsymcat(15),color=cgcolor('white'),noclip=0
	plots,xgri,y,psym=cgsymcat(6),color=cgcolor('red'),noclip=0

end

; compare observed volume densities of AGNs in pairs vs. expected volume
; densities
h = 0.7
minbclass = 2

; setup plot
loadct,0
psfile = 'figs/phi_sep2.eps'
nx = 3 
ny = 1
xtit='Transverse Separation (kpc)'
xtit = 'Pair Separation (kpc)'
ytit='log n'+textoidl(' (# per 10^6 Mpc^3)')
ytit='Volume Density ('+textoidl('log(# per 10^6 Mpc^3)')
setps,psfile,11*1.1*nx,11*ny,font='helvetica'
titsize = 1.5
labelsize = 1.2
multiplot,/default
multiplot,[nx,ny],xgap=0.015,/dox,/doy ;,mxtitle='',mytitle=ytit,$
	;mxtitsize=titsize,mytitsize=titsize,mxtitoffset=0.7

; load galaxy pair catalog
pcat = mrdfits('../7.sample/pair1kpc.fits',1)
pcat.plateifu = strtrim(pcat.plateifu) ; 176
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 149
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 119
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 113
	pcat.sep_arcsec le pcat.IFUrin)  ; 105
pcat = pcat[s]
; AGN flag - use NII/Ha BPT only, no WHAN 
flg_pagn = pcat.bpt[0] ge minbclass ; primary is AGN
flg_sagn = pcat.bpt[1] ge minbclass ; secondary is AGN
flg_bagn = min(pcat.bpt,dim=1) ge minbclass ; both comps are AGN

; load full target catalog w/ new weights
target = mrdfits('../../sample/target/DR14_newcosm.fits',1)
; match target catalog to pair catalog
match2,pcat.NSA_NSAID,target.NSA_NSAID,sa,sb
ptarg = target[sa]

; compute f_AGN for all pairs
z = pcat.NSA_Z
mass_p = alog10(pcat.NSA_ELPETRO_MASS/h^2)
mass_s = mass_p + (pcat.mstar[1] - pcat.mstar[0])
; take best-fit AGN model
;save,p,perr,pcerr,bestnorm,dof,filename='pd2agn2.idlsav'
;agnmod = p[0] * exp(-0.5d*(yy-p[1])^2/p[2]^2) * (1+xx)^4
restore,'pd2agn2.idlsav'
fagn_p =  p[0] * exp(-0.5d*(mass_p-p[1])^2/p[2]^2) * (1+z)^4
fagn_s =  p[0] * exp(-0.5d*(mass_s-p[1])^2/p[2]^2) * (1+z)^4

; define X grid for binning
x = pcat.sep_kpc 
;xbin =  5.0*2 ; kpc
;xmin =  5.0-xbin+1
;xmax = 25.0+xbin+1
;;xbin =  5.0*3 ; kpc
;;xmin =  7.5-xbin
;;xmax = 22.5+xbin
;xnp = fix(abs(xmax-xmin)/xbin)+1
;xgri = range(xmin,xmax,xnp)
;xerr = xgri*0 + xbin/2

xgri = [3.5,8.5,16,26]
xbin = [5.0,5.0,10,10]
;xgri = [5.0,15,25]
;xbin = [10,10,10]
;xgri = [3.5,8.5,13.5,18.5, 25.5]
;xbin = [5.0,5.0,5.0 ,5.0, 9]


xerr = xbin/2
xnp = n_elements(xgri)
; statistical weights to sum up
w = ptarg.ESWEIGHT ; per 1e6 Mpc^3
phi_pagn_exp = fltarr(xnp,3)
phi_sagn_exp = fltarr(xnp,3)
phi_bagn_exp = fltarr(xnp,3)
phi_pair_obs = fltarr(xnp,3)
phi_pagn_obs = fltarr(xnp,3)
phi_sagn_obs = fltarr(xnp,3)
phi_bagn_obs = fltarr(xnp,3)

nboot = 1000
for i=0,n_elements(xgri)-1 do begin
	; expected volume densities
	insidebin = x ge xgri[i]-xbin[i]/2 and x lt xgri[i]+xbin[i]/2
	ind = where(insidebin, npair)
	if npair gt 0 then begin
		ct = npair
		phi_pair_obs[i,*] = bootstrap_mean(w[ind],nboot=nboot)*ct
		phi_pagn_exp[i,*] = bootstrap_mean(w[ind]*fagn_p[ind],nboot=nboot)*ct
		phi_sagn_exp[i,*] = bootstrap_mean(w[ind]*fagn_s[ind],nboot=nboot)*ct
		phi_bagn_exp[i,*] = bootstrap_mean(w[ind]*fagn_p[ind]*fagn_s[ind],nboot=nboot)*ct
	endif
	; observed volume densities
	ind = where(flg_pagn and insidebin, ct)
	if ct gt 0 then begin
		if ct ge 5 then begin
			phi_pagn_obs[i,*] = bootstrap_mean(w[ind],nboot=nboot)*ct
		endif else begin
			phi_pagn_obs[i,1] = total(w[ind])
			err = binormial_ci(ct,npair)
			phi_pagn_obs[i,0] = total(w[ind])*(err[0]-err[4])/err[0]
			phi_pagn_obs[i,2] = total(w[ind])*(err[0]+err[4])/err[0]
		endelse
	endif
	ind = where(flg_sagn and insidebin, ct)
	if ct gt 0 then begin
		if ct ge 5 then begin
			phi_sagn_obs[i,*] = bootstrap_mean(w[ind],nboot=nboot)*ct
		endif else begin
			phi_sagn_obs[i,1] = total(w[ind])
			err = binormial_ci(ct,npair)
			phi_sagn_obs[i,0] = total(w[ind])*(err[0]-err[4])/err[0]
			phi_sagn_obs[i,2] = total(w[ind])*(err[0]+err[4])/err[0]
		endelse
	endif
	ind = where(flg_bagn and insidebin, ct)
	if ct gt 0 then begin
		print,xgri[i],ct,npair
		if ct ge 5 then begin
			phi_bagn_obs[i,*] = bootstrap_mean(w[ind],nboot=nboot)*ct
		endif else begin
			phi_bagn_obs[i,1] = total(w[ind])
			err = binormial_ci(ct,npair)
			phi_bagn_obs[i,0] = total(w[ind])*(err[0]-err[4])/err[0]
			phi_bagn_obs[i,2] = total(w[ind])*(err[0]+err[4])/err[0]
		endelse
	endif
endfor

; common arrays
xr = [-2.5,32.5]
;xr = [0.8,50]
yr = [-1.4,1.6]
;yr = [-3.0,1.7]
xx = pcat.sep_kpc
yy = alog10(ptarg.esweight) ; W = n*1e6 Mpc^3 => W is the # per 1e6 Mpc^3
idxp = where(flg_pagn) ; primary is AGN
idxs = where(flg_sagn) ; secondary is AGN
idxb = where(flg_bagn) ; both are AGN

; Panel 1: Primary as AGN
idx = idxp
phi_obs = phi_pagn_obs
phi_exp = phi_pagn_exp
mkplot,xx,yy,idx,xgri,xerr,phi_obs,phi_exp,ytit=ytit,$
	xtickint=10,xminor=5,ytickint=1.0,yminor=5,/xs,/ys,xr=xr,yr=yr,xtit=xtit
al_legend,'(a)',/top,/right,box=0,chars=1.5,margin=0
xyouts,mean(xr),yr[1]-0.3,'Primary as AGN',align=0.5,/data,charsize=labelsize
multiplot,/dox,/doy

; Panel 2: Secondary as AGN
idx = idxs
phi_obs = phi_sagn_obs
phi_exp = phi_sagn_exp
mkplot,xx,yy,idx,xgri,xerr,phi_obs,phi_exp,$
	xtickint=10,xminor=5,ytickint=1,yminor=5,/xs,/ys,xr=xr,yr=yr,xtit=xtit
al_legend,'(b)',/top,/right,box=0,chars=1.5,margin=0
xyouts,mean(xr),yr[1]-0.3,'Secondary as AGN',align=0.5,/data,charsize=labelsize
multiplot,/dox,/doy

; Panel 3: BAGN
idx = idxb
phi_obs = phi_bagn_obs
phi_exp = phi_bagn_exp
yr = yr-0.4
mkplot,xx,yy,idx,xgri,xerr,phi_obs,phi_exp,$
	xtickint=10,xminor=5,ytickint=1,yminor=5,/xs,/ys,xr=xr,yr=yr,xtit=xtit
; plot upper limit for last bin
i = 3
idx = where(x ge xgri[i]-xbin[i]/2 and x lt xgri[i]+xbin[i]/2)
plotsym,1,5.0,thick=4
oploterror,xgri[i],alog10(max(w[idx])),xerr[i],xbin[i]*0,psym=8,color=cgcolor('red')
; legend
al_legend,'(c)',/top,/right,box=0,chars=1.5,margin=0
xyouts,mean(xr),yr[1]-0.3,'Binary AGN',align=0.5,/data,charsize=labelsize
al_legend,['Weights','Observed','Expected'],psym=[9,6,27],$ ;linsize=[0.5,0.5],$
	color=cgcolor(['dark gray','red','royal blue']),$
	/bottom,/left,background=cgcolor('white'),charsize=labelsize
multiplot,/dox,/doy

theend:
multiplot,/default,/reset
device,/close
;spawn,'gv '+psfile+' &'


end
