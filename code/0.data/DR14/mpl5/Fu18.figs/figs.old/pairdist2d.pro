;+
; plot pair fraction in the plane of M* and Rifu
;-
pro mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin, _REF_EXTRA=extra
	;sauron_colormap
	loadct,0
	; reverse colormap
	TVLCT, r, g, b, /Get
	TVLCT, Reverse(r), Reverse(g), Reverse(b)
	; overwrite color table
	;mincolor = cgColor('White',0)
	maxcolor = cgColor('Black',!D.Table_Size-1)
	; aviod white color at minimum value
	s = where(map-zr[0] lt (zr[1]-zr[0])*1./256.,ct1)
	if ct1 gt 0 then map[s] = zr[0]+(zr[1]-zr[0])*1./256
	; aviod black color at maximum value
	s = where(map-zr[0] gt (zr[1]-zr[0])*255./256.,ct1)
	if ct1 gt 0 then map[s] = zr[0]+(zr[1]-zr[0])*255./256
	; use mincolor for NaNs
	ind = where(~finite(map),ct)
	if ct gt 0 then map[ind] = zr[0]
	; plotting ranges
	xr = [xmin,xmax]+1.0*xbin*[-1,1]
	yr = [ymin,ymax]+1.0*ybin*[-1,1]
	imgxr = [xmin,xmax]+0.5*xbin*[-1,1]
	imgyr = [ymin,ymax]+0.5*ybin*[-1,1]
	; plot frame
	plot,[0,1],[0,1],/nodata,xs=1+4,ys=1+4,xr=xr,yr=yr 
	; 2D histogram
	map2 = bytscl(map,min=zr[0],max=zr[1],/nan,top=255)
	oplotimage,map2,imgxrange=imgxr,imgyrange=imgyr
	; replot frame
	plot,[0,1],[0,1],/nodata,/noerase,/xs,/ys,xr=xr,yr=yr,color=maxcolor,_extra=extra 	
end

; setup plot
loadct,0
psfile = 'figs/pairdist2d.eps'
setps,psfile,15*1.1*5,15*3,font='helvetica'
titsize = 1.5
labelsize = 1.2
multiplot,/default
multiplot,[5,3],xgap=0.025,ygap=0.04,/dox,/doy

; load all DR14 catalog
ncat = mrdfits('../7.sample/nuc1kpc.fits',1) ; 2718 
s = where((ncat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0) ; 2618
ncat = ncat[s]
; load galaxy pair catalog
pcat = mrdfits('../7.sample/pair1kpc.fits',1)
pcat.plateifu = strtrim(pcat.plateifu) ; 176
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 149
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 119
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 113
	pcat.sep_arcsec le pcat.IFUrin)  ; 105
pcat = pcat[s]
; Pair index - PIND
match2,pcat.NSA_NSAID,ncat.NSA_NSAID,pind,sb
help,where(pcat.NSA_NSAID ne ncat[pind].NSA_NSAID)
;help,where(strtrim(ncat[pind].plateifu) ne strtrim(pcat.plateifu))

; load full target catalog w/ new weights
target = mrdfits('../../sample/target/DR14_newcosm.fits',1)
; match target catalog to DR14 catalog
match2,ncat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
ntarg = target[sa]
; match target catalog to pair catalog
match2,pcat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
ptarg = target[sa]

; define Rifu grid
xbin = 5.0/2 ; kpc
xmin = 2.5
xmax = 37.5
xnp = fix(abs(xmax-xmin)/xbin)+1
xgri = range(xmin,xmax,xnp)
; define stellar mass grid
ybin = 0.5/2 ; log(Msun)
ymin = 9.0
ymax = 11.75
ynp = fix(abs(ymax-ymin)/ybin)+1
ygri = range(ymin,ymax,ynp)

; count number density in each cell
h = 0.7 ; Hubble 
x = ncat.IFUrin*dangular(ncat.NSA_Z,/kpc)/(180.*3600./!pi) ; kpc
y = alog10(ncat.NSA_ELPETRO_MASS/h^2) ; log(Msun)

;; 2D histogram calculation
;nden = hist_2d(x,y,bin1=xbin,bin2=ybin,$
;	min1=xmin-xbin/2,max1=xmax+xbin/2,min2=ymin-ybin/2,max2=ymax+ybin/2)
;pden = hist_2d(x[pind],y[pind],bin1=xbin,bin2=ybin,$
;	min1=xmin-xbin/2,max1=xmax+xbin/2,min2=ymin-ybin/2,max2=ymax+ybin/2)
;; remove the extra column and row
;nx=(size(nden))[1]
;ny=(size(nden))[2]
;nden = nden[0:nx-2,0:ny-2]
;pden = pden[0:nx-2,0:ny-2]

; number of sources in each cell
nden = fltarr(xnp,ynp) ;  to nden from hist_2d
pden = fltarr(xnp,ynp)
aden = fltarr(xnp,ynp)
; sum of weights in each cell
nphi = fltarr(xnp,ynp)
pphi = fltarr(xnp,ynp)
aphi = fltarr(xnp,ynp)
for i=0,n_elements(xgri)-1 do begin
	for j=0,n_elements(ygri)-1 do begin
		; all targets
		ind = where(x ge xgri[i]-xbin/2 and x lt xgri[i]+xbin/2 and $
	       		y ge ygri[j]-ybin/2 and y lt ygri[j]+ybin/2,ct)
		if ct gt 0 then begin
			nden[i,j] = ct
			nphi[i,j] = total(ntarg[ind].ESWEIGHT)
		endif
		; pairs
		ind = where(x[pind] ge xgri[i]-xbin/2 and x[pind] lt xgri[i]+xbin/2 and $
	       		y[pind] ge ygri[j]-ybin/2 and y[pind] lt ygri[j]+ybin/2,ct)
		if ct gt 0 then begin
			pden[i,j] = ct
			pphi[i,j] = total(ntarg[pind[ind]].ESWEIGHT)
		endif
		; AGNs
		ind = where(x ge xgri[i]-xbin/2 and x lt xgri[i]+xbin/2 and $
	       		y ge ygri[j]-ybin/2 and y lt ygri[j]+ybin/2 and $
			ncat.bclass ge 2,ct)
		if ct gt 0 then begin
			aden[i,j] = ct
			aphi[i,j] = total(ntarg[ind].ESWEIGHT)
		endif
	endfor
endfor

; 1D mask
mid = where(nden lt 1)

; compute model pair fraction in each cell
fp_model = nden*0.0
for i=0,n_elements(xgri)-1 do begin
	for j=0,n_elements(ygri)-1 do begin
		z = median(ncat.nsa_z)
		m_min = ygri[j]-ybin/2 ;+0.28 ; log(Msun)
		m_max = ygri[j]+ybin/2 ;+0.28
		mu_min = 0.1 ; mass ratio
		mu_max = 1.0
		fgas_min = 0.0 ; gas fraction
		fgas_max = 1.0 
		; truncate at IFU radius
		sep_max = (xgri[i]*h) ; < 30*h) ; h^-1 kpc, maximum projected pair separation
		fp_model[i,j] = merger_rate_calculator(z,m_min,m_max,mu_min,mu_max,fgas_min,fgas_max,$
			RETURN_MERGER_FRACTION_PAIRS=sep_max,/USE_STELLAR_MASSES,$
			/USE_ALTERNATIVE_MASSFUNCTION,/quiet)
	endfor
endfor

; 2D pair fraction
; error of pair fraction
u = pden*1.0
du = sqrt(pden)
v = nden*1.0
dv = sqrt(nden)
ind = where(nden gt 0 and pden eq 0)
u[ind] = 1.0
du[ind] = 1.0
dfp2d = nden*0.0
dfp2d[*,*] = (eratio(u,du,v,dv))[*,1]
dfp2d[mid] = !values.f_nan
fp2d = nden*0.0
fp2d = (eratio(u,du,v,dv))[*,0]
fp2d[mid] = !values.f_nan
; compute best scaling for the model to fit the data 
data = fp2d
err  = dfp2d
model = fp_model
model[mid] = !values.f_nan
ind = where(finite(data))
for f=1.5,3.5,0.1 do begin
	chisq = total(((data[ind]-f*model[ind])/err[ind])^2)
	print,f,chisq
endfor
scale = 2.6

; plot
zr = [-0.1,1]*max(nden)
map = nden
map[mid] = !values.f_nan
mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin,$
	xtickint=10,xminor=5,ytickint=1.0,yminor=5,$
	xtit='IFU Radius (kpc)',$
	ytit=textoidl('log(Stellar Mass/M_{sun})'),charsize=titsize
oplot,x,y,psym=3,color=cgcolor('red') ;maxcolor
; contour
xx = xgri # (fltarr(ynp)+1.0)	
yy = (fltarr(xnp)+1.0) # ygri	
contour,nden,xx,yy,/data,/irregular,/overplot,levels=zr[1]*[0.1,0.2,0.5,0.8,0.9]
multiplot,/dox,/doy
;; legend
;;al_legend,'(a)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
;multiplot,/dox,/doy

; save key variables
save,xgri,ygri,xx,yy,fp_model,scale,nden,pden,nphi,pphi,filename='pairdist2d.sav'

zr = [-0.1, 0.4]
map = (pden*1.0)/nden
map[mid] = !values.f_nan
mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,$
	xtit='IFU Radius (kpc)',$
	ytit=textoidl('log(Stellar Mass/M_{sun})'),charsize=titsize
; legend
;al_legend,'(c)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

map = fp_model*scale
map[mid] = !values.f_nan
mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,$
	xtit='IFU Radius (kpc)',$
	ytit=textoidl('log(Stellar Mass/M_{sun})'),charsize=titsize
; legend
;al_legend,'(c)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

map = abs((pden*1.0)/nden-fp_model*scale)
map[mid] = !values.f_nan
mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,$
	xtit='IFU Radius (kpc)',$
	ytit=textoidl('log(Stellar Mass/M_{sun})'),charsize=titsize
; legend
;al_legend,'(c)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

zr = [0,3]
map = abs((pden*1.0)/nden-fp_model*scale)/dfp2d
map[mid] = !values.f_nan
mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,$
	xtit='IFU Radius (kpc)',$
	ytit=textoidl('log(Stellar Mass/M_{sun})'),charsize=titsize
; legend
;al_legend,'(c)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

;goto,theend

;;;;;;;;;;;;;;;
; plot summed weights
;;;;;;;;;;;;;;;
zr = minmax(nphi)
map = nphi
mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin,$
	xtickint=10,xminor=5,ytickint=1.0,yminor=5,$
	xtit='IFU Radius (kpc)',$
	ytit=textoidl('log(Stellar Mass/M_{sun})'),charsize=titsize
oplot,x,y,psym=3,color=cgcolor('red') ;maxcolor
; contour
xx = xgri # (fltarr(ynp)+1.0)	
yy = (fltarr(xnp)+1.0) # ygri	
contour,nphi,xx,yy,/data,/irregular,/overplot,levels=zr[1]*[0.1,0.2,0.5,0.8,0.9]
multiplot,/dox,/doy

zr = [-0.1, 0.4]
map = pphi/nphi
map[mid] = !values.f_nan
mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,$
	xtit='IFU Radius (kpc)',$
	ytit=textoidl('log(Stellar Mass/M_{sun})'),charsize=titsize
; legend
;oplot,x[pind],y[pind],psym=6,color=cgcolor('black') ;maxcolor
;al_legend,'(c)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

map = fp_model*scale
map[mid] = !values.f_nan
mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,$
	xtit='IFU Radius (kpc)',$
	ytit=textoidl('log(Stellar Mass/M_{sun})'),charsize=titsize
; legend
;al_legend,'(c)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

map = abs(pphi/nphi-fp_model*scale)
map[mid] = !values.f_nan
mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,$
	xtit='IFU Radius (kpc)',$
	ytit=textoidl('log(Stellar Mass/M_{sun})'),charsize=titsize
; legend
;al_legend,'(c)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

zr = [0,3]
map = abs(pphi/nphi-fp_model*scale)/dfp2d
map[mid] = !values.f_nan
mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,$
	xtit='IFU Radius (kpc)',$
	ytit=textoidl('log(Stellar Mass/M_{sun})'),charsize=titsize
; legend
;al_legend,'(c)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;
zr = [-0.1, 0.4]
map = aden/nden
;map[mid] = !values.f_nan
mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,$
	xtit='IFU Radius (kpc)',$
	ytit=textoidl('log(Stellar Mass/M_{sun})'),charsize=titsize
; legend
;oplot,x[pind],y[pind],psym=6,color=cgcolor('black') ;maxcolor
;al_legend,'(c)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

zr = [-0.1, 0.4]
map = aphi/nphi
;map[mid] = !values.f_nan
mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,$
	xtit='IFU Radius (kpc)',$
	ytit=textoidl('log(Stellar Mass/M_{sun})'),charsize=titsize
; legend
;oplot,x[pind],y[pind],psym=6,color=cgcolor('black') ;maxcolor
;al_legend,'(c)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

map = nden * (pden/nden) * (aden/nden)^2
map[mid] = !values.f_nan
zr = minmax(map,/nan)
mk2dhist,map,zr,xmin,xmax,xbin,ymin,ymax,ybin,$
	xtickint=10,xminor=4,ytickint=1.0,yminor=5,$
	xtit='IFU Radius (kpc)',$
	ytit=textoidl('log(Stellar Mass/M_{sun})'),charsize=titsize
; legend
;oplot,x[pind],y[pind],psym=6,color=cgcolor('black') ;maxcolor
;al_legend,'(c)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy


theend:
device,/close
;spawn,'gv '+psfile+' &'


end
