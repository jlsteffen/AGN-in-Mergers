;+
; compute and plot 2D distributions of AGNs in pairs 
; in the plane of mass and projected separation
;-

opt1 = 'sep_kpc'
opt2 = 'mass'
h = 0.7

; for AGN selection
minbclass = 2

; setup plot
loadct,0
nx = 3
ny = 1
psfile = 'phi2d.eps'
print,psfile
setps,psfile,12*nx,12*1.2*ny,font='helvetica'
titsize = 1.5
labelsize = 1.2
multiplot,/default
pos=[0.05,0.1,0.99,1.0]
multiplot,[nx,ny],xgap=0,ygap=0,position=pos
;,mxtitle=xtit,mytitle=ytit

; load galaxy pair catalog
pcat = mrdfits('../7.sample/pair1kpc.fits',1)
pcat.plateifu = strtrim(pcat.plateifu) ; 176
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 149
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 119
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 113
	pcat.sep_arcsec le pcat.IFUrin)  ; 105
pcat = pcat[s]
; AGN flag 
flg_agn1 = pcat.bclass[0] ge minbclass
flg_agn2 = pcat.bclass[1] ge minbclass
; load full target catalog w/ new weights
target = mrdfits('../../sample/target/DR14_newcosm.fits',1)
; match target catalog to pair catalog
match2,pcat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
ptarg = target[sa]

; compute f_AGN for all pairs
z = pcat.NSA_Z
mass_p = alog10(pcat.NSA_ELPETRO_MASS/h^2)
mass_s = mass_p + (pcat.mstar[1] - pcat.mstar[0])
; take best-fit AGN model
;save,p,perr,pcerr,bestnorm,dof,filename='pd2agn.idlsav'
;agnmod = p[0] * exp(-0.5d*(yy-p[1])^2/p[2]^2) * (1+xx)^4
restore,'pd2agn.idlsav'
fagn_p =  p[0] * exp(-0.5d*(mass_p-p[1])^2/p[2]^2) * (1+z)^4
fagn_s =  p[0] * exp(-0.5d*(mass_s-p[1])^2/p[2]^2) * (1+z)^4

; count number density in each cell
h = 0.7 ; Hubble 

; define X grid
x = pcat.sep_kpc ; kpc
xtit = 'Pair Separation (kpc)'
xbin = 2.5*2 ; kpc
xmin =  1.25 -xbin
xmax = 28.75 +xbin
xnp = fix(abs(xmax-xmin)/xbin)+1
xgri = range(xmin,xmax,xnp)
xr = [xmin,xmax]+0.5*xbin*[-1,1]
xtickint=10
xminor=4

; define Y grid
y = alog10(pcat.NSA_ELPETRO_MASS/h^2) 
ytit = textoidl('log(Stellar Mass/M_{sun})')
ybin = 0.25*2 ; log(Msun)
ymin =  9.125 - 3*ybin
ymax = 11.625 + 2*ybin
ynp = fix(abs(ymax-ymin)/ybin)+1
ygri = range(ymin,ymax,ynp)
yr = [ymin,ymax]+0.5*ybin*[-1,1]
ytickint=1
yminor=4

; number of sources in each cell
pden = fltarr(xnp,ynp) ; pairs
aden1 = fltarr(xnp,ynp) ; AGNs in primary
aden2 = fltarr(xnp,ynp) ; AGNs in secondary
; sum of weights in each cell
pphi = fltarr(xnp,ynp)
pphi1 = fltarr(xnp,ynp)
pphi2 = fltarr(xnp,ynp)
aphi1 = fltarr(xnp,ynp)
aphi2 = fltarr(xnp,ynp) 

for i=0,n_elements(xgri)-1 do begin
	for j=0,n_elements(ygri)-1 do begin
		; pairs
		ind = where(x ge xgri[i]-xbin/2 and x lt xgri[i]+xbin/2 and $
	       		y ge ygri[j]-ybin/2 and y lt ygri[j]+ybin/2,ct)
		if ct gt 0 then begin
			pden[i,j] = ct
			pphi[i,j] = total(ptarg[ind].ESWEIGHT)
			pphi1[i,j] = total(ptarg[ind].ESWEIGHT*fagn_p[ind]) 
			pphi2[i,j] = total(ptarg[ind].ESWEIGHT*fagn_s[ind]) 
		endif
		; AGNs in Primaries
		ind = where(flg_agn1 and x ge xgri[i]-xbin/2 and x lt xgri[i]+xbin/2 and $
	       		y ge ygri[j]-ybin/2 and y lt ygri[j]+ybin/2,ct)
		if ct gt 0 then begin
			aden1[i,j] = ct
			aphi1[i,j] = total(ptarg[ind].ESWEIGHT)
		endif
		; AGNs in Secondaries
		ind = where(flg_agn2 and x ge xgri[i]-xbin/2 and x lt xgri[i]+xbin/2 and $
	       		y ge ygri[j]-ybin/2 and y lt ygri[j]+ybin/2,ct)
		if ct gt 0 then begin
			aden2[i,j] = ct
			aphi2[i,j] = total(ptarg[ind].ESWEIGHT)
		endif
	endfor
endfor

;; save key variables
;xx = xgri # (fltarr(ynp)+1.0)	
;yy = (fltarr(xnp)+1.0) # ygri	
;save,xx,yy,xr,yr, pden,aden1,aden2,pphi,pphi1,pphi2,aphi1,aphi2,$
;	filename=repstr(psfile,'.eps','.sav')

; 1D mask
mid = where(pden lt 1)

;;;;;;;;;
map = alog10(aphi1+aphi2)
zr = [-1,1] ; minmax(map)
zminor=4
ctable = 62
manga_mkplot,map,zr=zr,/colorbar,zsize=labelsize,$ 
	ztit=textoidl('\Phi_{agn} per cell (2.5 kpc \times 0.25 dex)'),$
	xr=xr,yr=yr,xtit='',ytit=ytit,$
	xtickint=xtickint,xminor=xminor,ytickint=ytickint,yminor=yminor,$
	ztickint=ztickint,zminor=zminor,$
	position=!p.position,sauron=ctable,charsize=titsize
	;c_map=c_map,c_levels=c_levels
multiplot

;;;;;;;;;
map = alog10(pphi1+pphi2)
map[mid] = !values.f_nan
manga_mkplot,map,zr=zr,/colorbar,zsize=labelsize,$ 
	ztit=textoidl('Expected \Phi_{agn} per cell'),$
	xr=xr,yr=yr,xtit=xtit,ytit='',$
	xtickint=xtickint,xminor=xminor,ytickint=ytickint,yminor=yminor,$
	ztickint=ztickint,zminor=zminor,$
	position=!p.position,sauron=ctable,charsize=titsize
oplot,x,y,psym=3,color=cgcolor('black')
multiplot

;;;;;;;;;
map = alog10(aphi1+aphi2)-alog10(pphi1+pphi2)
map[mid] = !values.f_nan
manga_mkplot,map,zr=zr,/colorbar,zsize=labelsize,$ 
	ztit=textoidl('Expected \Phi_{agn} per cell'),$
	xr=xr,yr=yr,xtit=xtit,ytit='',$
	xtickint=xtickint,xminor=xminor,ytickint=ytickint,yminor=yminor,$
	ztickint=ztickint,zminor=zminor,$
	position=!p.position,sauron=ctable,charsize=titsize
multiplot


theend:
multiplot,/default,/reset
device,/close
;spawn,'gv '+psfile+' &'


end
