; plot sample distributions in Mi and redshift

PRO buildhist, x, y, xmin, xmax, bin, xh=xh, yh=yh ; bincen=bincen, phi=phi
	; build histogram based on X and weight (Y)
	np = fix(abs(xmax-xmin)/bin)+1
	bincen = range(xmin,xmax,np)
	phi = fltarr(np)
	for i=0,np-1 do begin
		s = where(x ge bincen[i]-bin/2 and x lt bincen[i]+bin/2,ct)
		if ct gt 0 then phi[i] = total(y[s])
	endfor
	; rearrange the array for histogram plots
	xh = fltarr(np*2)
	xh[0:*:2] = bincen-bin/2
	xh[1:*:2] = bincen+bin/2
	yh = fltarr(np*2)
	yh[0:*:2] = phi
	yh[1:*:2] = phi
END

PRO mkhist, x, y, xmin, xmax, bin, bc, pind, _REF_EXTRA=extra
	; plot frame
	plot,[0,1],[0,1],/nodata,xs=1+4,ys=1+4,_extra=extra
	; all DR14 - provide denominator
	buildhist, x, y, xmin, xmax, bin, xh=xh, yh=yh0
	np = n_elements(xh)
	; all DR14 AGNs - provide numerator
	s = where(bc ge 2)
	buildhist, x[s], y[s], xmin, xmax, bin, xh=xh, yh=yh
	polyfill,[xh[0],xh,xh[np-1]],[0,yh/yh0,0],/fill,/data,color=cgcolor('light gray')
	oplot,[xh[0],xh,xh[np-1]],[0,yh/yh0,0],lines=0,color=cgcolor('black')
	;; all pairs
	;s = pind
	;buildhist, x[s], y[s], xmin, xmax, bin, xh=xh, yh=yh
	;oplot,[xh[0],xh,xh[np-1]],[0,yh/yh0,0],lines=2,color=cgcolor('black')
	; break down different AGN populations
	s = where(bc eq 2) ; Composites
	buildhist, x[s], y[s], xmin, xmax, bin, xh=xh, yh=yh
	oplot,[xh[0],xh,xh[np-1]],[0,yh/yh0,0],lines=0,color=cgcolor('royal blue')
	polyfill,[xh[0],xh,xh[np-1]],[0,yh/yh0,0],/line_fill,/data,orient=0,color=cgcolor('royal blue')
	s = where(bc eq 3) ; LINER
	buildhist, x[s], y[s], xmin, xmax, bin, xh=xh, yh=yh	
	oplot,[xh[0],xh,xh[np-1]],[0,yh/yh0,0],lines=0,color=cgcolor('dark green')
	polyfill,[xh[0],xh,xh[np-1]],[0,yh/yh0,0],/line_fill,orient=-45,$
		/data,color=cgcolor('forest green')
	s = where(bc ge 4) ; Seyfert 1 & 2
	buildhist, x[s], y[s], xmin, xmax, bin, xh=xh, yh=yh	
	oplot,[xh[0],xh,xh[np-1]],[0,yh/yh0,0],lines=0,color=cgcolor('red')
	polyfill,[xh[0],xh,xh[np-1]],[0,yh/yh0,0],/line_fill,orient=45,/data,color=cgcolor('red')
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
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 168
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 128
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 122
	pcat.sep_arcsec le pcat.IFURin)  ; 114
pcat = pcat[s]
; Pair index
match2,pcat.NSA_NSAID,ncat.NSA_NSAID,pind,sb
help,where(pcat.NSA_NSAID ne ncat[pind].NSA_NSAID)

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

; setup plot
loadct,0
psfile = 'figs/agnfrac.eps'
setps,psfile,15*1.1*2,15*1,font='helvetica'
titsize = 1.5
labelsize = 1.2
multiplot,/default
multiplot,[2,1],xgap=0.025,ygap=0.02,/dox,/doy

; input data
x = ncat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y = ntarg.ESWEIGHT
bc = ncat.bclass
pind = pind
xtit = textoidl('Absolute Magnitude (M_i)')
ytit = textoidl('AGN Fraction')
xr = [-17,-25]
yr = [0,0.25]
xmin = -23.5
xmax = -18.5
bin = 1.0
mkhist, x, y, xmin, xmax, bin, bc, pind, $
	xtickint=2.0,xminor=4,ytickint=0.1,yminor=10,$	
	xtit=xtit,ytit=ytit,xr=xr,yr=yr,charsize=titsize
; legend
al_legend,'(a)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
al_legend,['All AGNs','Composites','LINERs','Seyferts'],$ ; ,'AGN Fraction'],$
	textcolors=cgcolor(['Black','Blue','dark green','Red']),$ ;,'Dark Green'])
	/top,/left

multiplot,/dox,/doy

; input data
h = 0.7 ; hubble constant
x = alog10(ncat.NSA_ELPETRO_MASS/h^2) 
y = ntarg.ESWEIGHT
bc = ncat.bclass
pind = pind
xtit = textoidl('Stellar Mass (log(M_{sun}))')
ytit = '' ; textoidl('Number of sources')
xr = [8.5,12.5]
yr = [0,0.25]
xmin = 8.75
xmax = 11.75
bin =  0.5
mkhist, x, y, xmin, xmax, bin, bc, pind, $
	xtickint=1.0,xminor=10,ytickint=0.1,yminor=10,$
	xtit=xtit,ytit=ytit,xr=xr,yr=yr,charsize=titsize
; legend
al_legend,'(b)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

;; input data
;x = ncat.NSA_Z
;y = ntarg.ESWEIGHT
;bc = ncat.bclass
;pind = pind
;xtit = 'Redshift'
;ytit = '' ; textoidl('Number of sources')
;xr = [0.0,0.16]
;yr = [0,0.25]
;xmin = 0.02
;xmax = 0.14
;bin =  0.02
;mkhist, x, y, xmin, xmax, bin, bc, pind, xtickint=0.05,xtit=xtit,ytit=ytit,xr=xr,yr=yr,charsize=titsize
;; legend
;al_legend,'(b)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
;multiplot,/dox,/doy

device,/close
spawn,'gv '+psfile+' &'

end
