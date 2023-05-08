; plot distributions of IFUs that contain AGNs 
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
	; show fraction
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
pind = where(ncat.bclass ge 2)

; setup plot
loadct,0
psfile = 'figs/agndist.eps'
setps,psfile,15*1.1*3,15*1,font='helvetica'
titsize = 1.5
labelsize = 1.2
multiplot,/default
multiplot,[3,1],xgap=0.025,ygap=0.02,/dox,/doy

; angular scale
scale = dangular(ncat.NSA_Z,/kpc)/(180.*3600./!pi) ; kpc/arcsec
x = ncat.IFUrin*scale ; IFU inner radius in kpc
bin = 5.0
mkhist, x, pind, bin, xr=[-2,42], yr=[0,0.5],$
	xtickint=10,xminor=5,ytickint=0.1,yminor=5,$
	xtit='IFU Radius (kpc)',ytit='Fraction',charsize=titsize
; legend
al_legend,'(a)',/top,/left,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

; plot histogram
x = ncat.IFUrin/ncat.NSA_ELPETRO_TH50_R 
bin = 0.4
mkhist, x, pind, bin, xr=[0.2,4.3], yr=[0,0.4],$
	xtickint=1,xminor=5,ytickint=0.1,yminor=5,$
	xtit='IFU Radius / Effective Radius',ytit='',charsize=titsize
hor,0.16,lines=2
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
al_legend,['All MaNGA','All AGNs','AGN Fraction'],$ ; ,'AGN Fraction'],$
	textcolors=cgcolor(['Black','Blue','Red']),$ ;,'Dark Green'])
	charsize=titsize,/top,/right

multiplot,/dox,/doy

device,/close
spawn,'gv '+psfile+' &'

end


