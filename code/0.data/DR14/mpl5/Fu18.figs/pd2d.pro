;+
; compute and plot pair and AGN fraction in the plane of 
; (1) mass/Mi and Rifu and (2) mass/Mi and redshift
; save figure in a figs/.esp file
; save 2D histograms in a figs/.sav file
;-

;opt1 = 'Rifu'
opt1 = 'z'

;opt2 = 'Mi'
opt2 = 'mass'

; for AGN selection
minbclass = 2
;minbclass = 3 

; setup plot
loadct,0
nx = 3
ny = 1
psfile = 'figs/pd2d_'+opt1+'_'+opt2+'.eps' ; '_'+strc(minbclass)+'.eps'
print,psfile
setps,psfile,12*nx,12*1.2*ny,font='helvetica'
titsize = 1.5
labelsize = 1.2
multiplot,/default
pos=[0.05,0.1,0.99,1.0]
multiplot,[nx,ny],xgap=0,ygap=0,position=pos
;,mxtitle=xtit,mytitle=ytit

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
; Pair flag
match2,pcat.NSA_NSAID,ncat.NSA_NSAID,sa,sb
flg_pair = sb ge 0 
; AGN flag 
flg_agn = ncat.bclass ge minbclass

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

; count number density in each cell
h = 0.7 ; Hubble 

if opt1 eq 'Rifu' then begin
; define X grid
x = ncat.IFUrin*dangular(ncat.NSA_Z,/kpc)/(180.*3600./!pi) ; kpc
xtit = 'IFU Radius (kpc)'
xbin = 2.5 ; kpc
xmin =  1.25 -xbin
xmax = 38.75 +xbin
xnp = fix(abs(xmax-xmin)/xbin)+1
xgri = range(xmin,xmax,xnp)
xr = [xmin,xmax]+0.5*xbin*[-1,1]
xtickint=10
xminor=4
endif

if opt1 eq 'z' then begin
; define X grid
x = ncat.NSA_Z
xtit = 'redshift'
xbin = 0.01 
xmin = 0.01 -1.5*xbin
xmax = 0.15 +1.5*xbin
xnp = fix(abs(xmax-xmin)/xbin)+1
xgri = range(xmin,xmax,xnp)
xr = [xmin,xmax]+0.5*xbin*[-1,1]
xtickint=0.05
xminor=5
endif

if opt2 eq 'Mi' then begin
; define Y grid
y = ncat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(h) 
ytit = textoidl('Absolute Magnitude (M_i)')
ybin = 0.5 
ymax = -24.0 - ybin
ymin = -17.5 + ybin
ynp = fix(abs((ymax-ymin)/ybin))+1
ygri = range(ymin,ymax,ynp)
yr = [ymin,ymax]-0.5*ybin*[-1,1]
ytickint=2
yminor=4
endif

if opt2 eq 'mass' then begin
; define Y grid
y = alog10(ncat.NSA_ELPETRO_MASS/h^2) 
ytit = textoidl('log(Stellar Mass/M_{sun})')
ybin = 0.25 ; log(Msun)
ymin =  9.125 - 3*ybin
ymax = 11.625 + 2*ybin
ynp = fix(abs(ymax-ymin)/ybin)+1
ygri = range(ymin,ymax,ynp)
yr = [ymin,ymax]+0.5*ybin*[-1,1]
ytickint=1
yminor=4
endif

; number of sources in each cell
nden = fltarr(xnp,ynp) ; all MaNGA
pden = fltarr(xnp,ynp) ; pairs
cden = fltarr(xnp,ynp) ; Control = All - Pairs
aden = fltarr(xnp,ynp) ; AGNs in Control Sample

; sum of weights in each cell
nphi = fltarr(xnp,ynp)
pphi = fltarr(xnp,ynp)
cphi = fltarr(xnp,ynp)
aphi = fltarr(xnp,ynp)

; Hopkins model pair fraction
h10mod = fltarr(xnp,ynp)

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
		ind = where(flg_pair and x ge xgri[i]-xbin/2 and x lt xgri[i]+xbin/2 and $
	       		y ge ygri[j]-ybin/2 and y lt ygri[j]+ybin/2,ct)
		if ct gt 0 then begin
			pden[i,j] = ct
			pphi[i,j] = total(ntarg[ind].ESWEIGHT)
		endif
		; Control sample 
		ind = where(~flg_pair and x ge xgri[i]-xbin/2 and x lt xgri[i]+xbin/2 and $
	       		y ge ygri[j]-ybin/2 and y lt ygri[j]+ybin/2,ct)
		if ct gt 0 then begin
			cden[i,j] = ct
			cphi[i,j] = total(ntarg[ind].ESWEIGHT)
		endif
		; AGNs in Control Sample
		ind = where(~flg_pair and flg_AGN and x ge xgri[i]-xbin/2 and x lt xgri[i]+xbin/2 and $
	       		y ge ygri[j]-ybin/2 and y lt ygri[j]+ybin/2,ct)
		if ct gt 0 then begin
			aden[i,j] = ct
			aphi[i,j] = total(ntarg[ind].ESWEIGHT)
		endif
		; Hopkins model pair fraction
		if opt1 eq 'Rifu' and opt2 eq 'mass' and xgri[i] gt 0 then begin
			z = 0.04
			m_min = ygri[j]-ybin/2 ;+0.28 ; log(Msun)
			m_max = ygri[j]+ybin/2 ;+0.28
			mu_min = 0.1 ; mass ratio
			mu_max = 1.0
			fgas_min = 0.0 ; gas fraction
			fgas_max = 1.0 
			; truncate at IFU radius
			sep_max = (xgri[i]*h) ; < 30*h) ; h^-1 kpc, maximum projected pair separation
			h10mod[i,j] = merger_rate_calculator(z,m_min,m_max,mu_min,mu_max,fgas_min,fgas_max,$
				RETURN_MERGER_FRACTION_PAIRS=sep_max,/USE_STELLAR_MASSES,$
				;/USE_ALTERNATIVE_MASSFUNCTION,$
				/quiet)
		endif
	endfor
endfor

; save key variables
xx = xgri # (fltarr(ynp)+1.0)	
yy = (fltarr(xnp)+1.0) # ygri	
save,xx,yy,xr,yr, nden,pden,cden,aden, nphi,pphi,cphi,aphi, h10mod,$
	filename=repstr(psfile,'.eps','.sav')

; 1D mask
mid = where(nden lt 1)

; contours
;c_map = nphi
;c_levels = max(c_map)*range(0.15,0.9,5)

;;;;;;;;;
map = nden
zr = [0,1]*max(map)
ztickint=20
zminor=4
ctable = 62
manga_mkplot,map,zr=zr,/colorbar,zsize=labelsize,$ 
	ztit=textoidl('N_{gal} per cell (2.5 kpc \times 0.25 dex)'),$
	xr=xr,yr=yr,xtit='',ytit=ytit,$
	xtickint=xtickint,xminor=xminor,ytickint=ytickint,yminor=yminor,$
	ztickint=ztickint,zminor=zminor,$
	position=!p.position,sauron=ctable,charsize=titsize
	;c_map=c_map,c_levels=c_levels
oplot,x,y,psym=3,color=cgcolor('black')
multiplot

;;;;;;;;;
map = pphi/nphi
map[mid] = !values.f_nan
zr = [-0.0,1.0]*0.4
ztickint=0.1
zminor=4
manga_mkplot,map,zr=zr,/colorbar,zsize=labelsize,$ 
	ztit=textoidl('Pair Fraction ( f_{pair}=\Phi_{Pair}/\Phi_{All} )'),$
	xr=xr,yr=yr,xtit=xtit,ytit='',$
	xtickint=xtickint,xminor=xminor,ytickint=ytickint,yminor=yminor,$
	ztickint=ztickint,zminor=zminor,$
	position=!p.position,sauron=ctable,charsize=titsize
multiplot

;;;;;;;;;
map = aphi/cphi
map[mid] = !values.f_nan
manga_mkplot,map,zr=zr,/colorbar,zsize=labelsize,$ 
	ztit=textoidl('AGN Fraction ( f_{agn}=\Phi_{AGN}/\Phi_{Control} )'),$
	xr=xr,yr=yr,xtit='',ytit='',$
	xtickint=xtickint,xminor=xminor,ytickint=ytickint,yminor=yminor,$
	ztickint=ztickint,zminor=zminor,$
	position=!p.position,sauron=ctable,charsize=titsize

theend:
multiplot,/default,/reset
device,/close
;spawn,'gv '+psfile+' &'


end
