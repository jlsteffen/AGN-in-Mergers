; plot distributions of IFUs that contain galaxy pairs 
; in Rifu/Re and Rifu in kpc

pro mkhist, x, pind, bin, _REF_EXTRA=extra
	; plot frame
	plot,[0,1],[0,1],/nodata,xs=1+4,ys=1+4,_extra=extra 
	; compute histogram values
	plothist,x,xh,yh,bin=bin,/noplot
	plothist,x[pind],xh1,yh1,bin=bin,/noplot
	; match the two histograms
	match2,strc(xh),strc(xh1),sa,sb
	xh = xh[sb]
	yh = yh[sb]
	; plot histograms
	plothist,x,bin=bin,peak=max(yh)*1.0/total(yh),/overplot,$
		/fill,fcolor=cgcolor('light gray')
	plothist,x[pind],bin=bin,peak=max(yh1)*1.0/total(yh1),/overplot,$
		/fill,/fline,forient=60,color=cgcolor('royal blue'),$
		fcolor=cgcolor('royal blue')
	; show Pair fraction vs. IFU size
	xerr = xh*0+bin/2
	tmp  = fltarr(5,n_elements(xh))
	for i=0,n_elements(xh)-1 do tmp[*,i] = binormial_ci(yh1[i],yh[i])
	y = tmp[0,*]
	ylerr = tmp[1,*]
	yuerr = tmp[2,*]
	oploterror,xh,y,xerr,ylerr,/lobar,psym=3,color=cgcolor('red')
	oploterror,xh,y,xerr,yuerr,/hibar,psym=3,color=cgcolor('red')
	plots,xh,y,psym=cgsymcat(15),color=cgcolor('white'),noclip=0
	plots,xh,y,psym=cgsymcat(6),color=cgcolor('red'),noclip=0
	;print,minmax(y),mean(y)
	; replot frame
	plot,[0,1],[0,1],/nodata,/noerase,/xs,/ys,_extra=extra 	
end

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
help,where(strtrim(ncat[pind].plateifu) ne strtrim(pcat.plateifu))

;; load full target catalog w/ new weights
;target = mrdfits('../../sample/target/DR14_newcosm.fits',1)
;; match target catalog to DR14 catalog
;match2,ncat.NSA_NSAID,target.NSA_NSAID,sa,sb
;help,where(sa lt 0) ; should be -1
;ntarg = target[sa]
;; match target catalog to pair catalog
;match2,pcat.NSA_NSAID,target.NSA_NSAID,sa,sb
;help,where(sa lt 0) ; should be -1
;ptarg = target[sa]

; setup plot
loadct,0
psfile = 'figs/pairdist.eps'
setps,psfile,15*1.1*4,15*1,font='helvetica'
titsize = 1.5
labelsize = 1.2
multiplot,/default
multiplot,[4,1],xgap=0.025,ygap=0.02,/dox,/doy

; angular scale
scale = dangular(ncat.NSA_Z,/kpc)/(180.*3600./!pi) ; kpc/arcsec
x = ncat.IFUrin*scale ; IFU inner radius in kpc
bin = 5.0
mkhist, x, pind, bin, xr=[-2,42], yr=[0,0.5],$
	xtickint=10,xminor=5,ytickint=0.1,yminor=5,$
	xtit='IFU Radius (kpc)',ytit='Fraction',charsize=titsize
; compute mean redshifts & stellar mass for each bin
xmax=37.5
xmin=2.5
np = fix(abs(xmax-xmin)/bin)+1
rifu = range(xmin,xmax,np)
meanz = fltarr(np)
meanms= fltarr(np)
fpair = fltarr(np)
fpair2 = fltarr(np)
fpair3 = fltarr(np)
h = 0.7 ; Hubble constant
for i=0,np-1 do begin
	s = where(x ge rifu[i]-bin/2 and x lt rifu[i]+bin/2,ct)
	if ct gt 0 then begin
		meanz[i] = mean(ncat[s].NSA_z)
		; NSA stellar masses are given in units of h^-2 Msun
		; Chabrier (2003) IMF (Blanton & Roweis 2007, K-correction paper)
		meanms[i] = alog10(mean(ncat[s].NSA_ELPETRO_MASS/h^2))
	endif
	; full Hopkins10 model
	redshift = meanz[i]
	m_min = meanms[i]-0.05 ; log(Msun)
	m_max = meanms[i]+0.05
	mu_min = 0.1 ; mass ratio
	mu_max = 1.0
	fgas_min = 0.0 ; gas fraction
	fgas_max = 1.0 
	; truncate at IFU radius
	sep_max = rifu[i]*h ; h^-1 kpc, maximum projected pair separation
	fpair[i] = merger_rate_calculator(redshift,m_min,m_max,mu_min,mu_max,fgas_min,fgas_max,$
		RETURN_MERGER_FRACTION_PAIRS=sep_max,/USE_STELLAR_MASSES,$
		/USE_ALTERNATIVE_MASSFUNCTION,/quiet)
	; truncate at 30kpc
	sep_max = (rifu[i]*h < 30*h) ; h^-1 kpc, maximum projected pair separation
	fpair2[i] = merger_rate_calculator(redshift,m_min,m_max,mu_min,mu_max,fgas_min,fgas_max,$
		RETURN_MERGER_FRACTION_PAIRS=sep_max,/USE_STELLAR_MASSES,$
		/USE_ALTERNATIVE_MASSFUNCTION,/quiet)
	; shift stellar mass
	fpair3[i] = $
	merger_rate_calculator(redshift,m_min+0.28,m_max+0.28,mu_min,mu_max,fgas_min,fgas_max,$
		RETURN_MERGER_FRACTION_PAIRS=sep_max,/USE_STELLAR_MASSES,$
		/USE_ALTERNATIVE_MASSFUNCTION,/quiet)
endfor
forprint,rifu,meanz,meanms,fpair,fpair2,fpair3
; compare w/ full Hopkins10 model
sep = range(0.0,42.0)
oplot,sep,interpol(fpair2,rifu,sep),lines=3
oplot,sep,interpol(fpair3,rifu,sep),lines=0

; t_fric ~ sep - dynamical friction timescale
maxfp  = 0.26 ; maximum pair fraction
maxsep = 30   ; kpc
t_fric = maxfp * (sep/maxsep)
oplot,sep,(t_fric<0.26),lines=2

;; projection effect, dN/di = sin(i) => <sin i> = !pi/4 = 0.7854
;maxsep = 30./(!pi/4) ; maximum intrinsic separation that can be detected
;maxfp  = 0.26 ; maximum pair fraction
;t_fric = maxfp * (sep/maxsep)
;oplot,sep*!pi/4,(t_fric<0.26),lines=2

; compare w/ best-fit Hopkins10 model
restore,'pairdist2d.sav'
;mid = where(nden lt 1)
;fp_model[mid] = !values.f_nan
;y = mean(scale*fp_model,dim=2)
y = total(fp_model*nden*scale,2)/total(nden,2)
oplot,xgri,y,psym=10,lines=2,color=cgcolor('dark green')
;oplot,sep,interpol(y,xgri,sep),lines=2,color=cgcolor('red')

; legend
al_legend,'(a)',/top,/left,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

; plot histogram
x = ncat.IFUrin/ncat.NSA_ELPETRO_TH50_R 
bin = 0.4
mkhist, x, pind, bin, xr=[0.2,4.3], yr=[0,0.4],$
	xtickint=1,xminor=5,ytickint=0.1,yminor=5,$
	xtit='IFU Radius / Effective Radius',ytit='',charsize=titsize
hor,0.04,lines=2
oplot,[1.5,1.5],[0,0.36],lines=2
oplot,[2.5,2.5],[0,0.36],lines=2
xyouts,1.5,0.37,'Primary+',/data,align=0.5
xyouts,2.5,0.37,'Secondary',/data,align=0.5
; legend
al_legend,'(b)',/top,/left,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

; Mi - absolute magnitudes
h = 0.7 ; hubble constant
x = ncat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(h) 
xr = [-25,-17]
xtit = textoidl('Absolute Magnitude (M_i)')
bin = 0.8
mkhist, x, pind, bin, xr=xr, yr=[0,0.5],$
	xtickint=2,xminor=4,ytickint=0.1,yminor=5,$
	xtit=xtit,ytit='',charsize=titsize
; legend
al_legend,'(c)',/top,/left,box=0,chars=2,margin=0; ,background=cgcolor('white')
al_legend,['All MaNGA','Close Pairs','Pair Fraction'],$ ; ,'AGN Fraction'],$
	textcolors=cgcolor(['Black','Blue','Red']),$ ;,'Dark Green'])
	charsize=titsize,/top,/right

multiplot,/dox,/doy

; Stellar mass
x = alog10(ncat.NSA_ELPETRO_MASS/h^2) 
xr = [8.5,12.5]
xtit = textoidl('Stellar Mass (log(M_{sun}))')
bin = 0.25
mkhist, x, pind, bin, xr=xr, yr=[0,0.5],$
	xtit=xtit,ytit='',charsize=titsize
; compare w/ best-fit Hopkins10 model
restore,'pairdist2d.sav'
y = total(fp_model*nden*scale,1)/total(nden,1)
oplot,ygri,y,psym=10,lines=2,color=cgcolor('dark green')
; legend
al_legend,'(d)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

device,/close
;spawn,'gv '+psfile+' &'

end

