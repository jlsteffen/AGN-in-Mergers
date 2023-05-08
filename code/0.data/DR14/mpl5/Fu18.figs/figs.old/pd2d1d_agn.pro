;+
; plot AGN fraction in the plane of Mass and Rifu
; use manga_mkplot for 2D
; use mkhist for 1D
; This version combines 2D and 1D distributions in one plot
;-

pro mkhist, x, pind, bin, r90=r90, _REF_EXTRA=extra
	; compute histogram values
	plothist,x,xh,yh,bin=bin,/noplot
	plothist,x[pind],xh1,yh1,bin=bin,/noplot
	; match the two histograms
	match2,strc(xh),strc(xh1),sa,sb
	xh = xh[sb]
	yh = yh[sb]
	; compute fractions 
	xerr = xh*0+bin/2
	tmp  = fltarr(5,n_elements(xh))
	for i=0,n_elements(xh)-1 do tmp[*,i] = binormial_ci(yh1[i],yh[i])
	y = tmp[0,*]
	ylerr = tmp[1,*]
	yuerr = tmp[2,*]
	; plot histograms
	if ~keyword_set(r90) then begin
		plothist,x,bin=bin,peak=max(yh)*1.0/total(yh),$
			/fill,fcolor=cgcolor('light gray'),xs=1+4,ys=1+4,_extra=extra
		plothist,x[pind],bin=bin,peak=max(yh1)*1.0/total(yh1),/overplot,$
			/fill,/fline,forient=60,color=cgcolor('royal blue'),$
			fcolor=cgcolor('royal blue')
		; show fraction
		oploterror,xh,y,xerr,ylerr,/lobar,psym=3,color=cgcolor('red')
		oploterror,xh,y,xerr,yuerr,/hibar,psym=3,color=cgcolor('red')
		plots,xh,y,psym=cgsymcat(15),color=cgcolor('white'),noclip=0
		plots,xh,y,psym=cgsymcat(6),color=cgcolor('red'),noclip=0
		plot,[0,1],[0,1],/nodata,/noerase,/xs,/ys,_extra=extra 		
	endif else begin
		plothist,x,bin=bin,peak=max(yh)*1.0/total(yh),/rotate,$
			/fill,fcolor=cgcolor('light gray'),xs=1+4,ys=1+4,_extra=extra
		plothist,x[pind],bin=bin,peak=max(yh1)*1.0/total(yh1),/rotate,/overplot,$
			/fill,/fline,forient=60,color=cgcolor('royal blue'),$
			fcolor=cgcolor('royal blue')
		; show fraction
		oploterror,y,xh,ylerr,xerr,/lobar,psym=3,color=cgcolor('red')
		oploterror,y,xh,yuerr,xerr,/hibar,psym=3,color=cgcolor('red')
		plots,y,xh,psym=cgsymcat(15),color=cgcolor('white'),noclip=0
		plots,y,xh,psym=cgsymcat(6),color=cgcolor('red'),noclip=0
		plot,[0,1],[0,1],/nodata,/noerase,xs=1+8,/ys,_extra=extra 		
	endelse
end

; setup plot
nx = 3
ny = 2
psfile = 'figs/pd2d1d_agn.eps'
setps,psfile,12*1.4,12*1.5,font='helvetica'
titsize = 1.5
labelsize = 1.2

pos1 = [0.12,0.4,0.65,1.0]
pos2 = [0.12,0.1,0.65,0.4]
pos3 = [0.65,0.4,0.98,0.88]
pos4 = [0.65,0.1,0.98,0.4]

h = 0.7 ; Hubble

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
flg_pair = sb ge 0
; Control Sample
s = where(~flg_pair)
ncat = ncat[s]
; AGN index
pind = where(ncat.bclass ge 2)

blanks = replicate(' ',10)
; restore previously saved stuff
restore,'figs/pd2d_Rifu_mass.sav'
map = aden/cden
mid = where(cden lt 1)
zr = [0,0.4]
map[mid] = !values.f_nan
ctable = 62
c_map = cden
c_levels = max(c_map)*range(0.05,0.8,4)
;c2_map = pden
;c2_levels = max(c2_map)*range(0.1,0.8,3)
manga_mkplot,map,zr=zr,/colorbar,zsize=labelsize,$ 
	ztit=textoidl('AGN Fraction'),$
	xr=xr,yr=yr,xtit='',ytit='Stellar Mass',$
	xtickint=10,xminor=4,ytickint=1,yminor=4,$
	xtickname=blanks,$
	position=pos1,sauron=ctable,charsize=titsize,$
	c_map=c_map,c_levels=c_levels,/noerase ;,$
	;c2_map=c_map2,c2_levels=c_levels2
al_legend,'(a)',/top,/left,box=0,chars=labelsize,margin=0,background=cgcolor('white')

; fraction range to plot
fr = [-0.05,0.45]
; histogram 1
loadct,0
; IFU size in kpc
x = ncat.IFUrin*dangular(ncat.NSA_Z,/kpc)/(180.*3600./!pi) ; IFU inner radius in kpc
bin = 5.0
mkhist, x, pind, bin, /noerase, pos=pos2, xr=xr, yr=fr,$
	xtickint=10,xminor=4,ytickint=0.2,yminor=4,$
	xtit='IFU Radius (kpc)',ytit='Fraction',charsize=titsize
al_legend,'(c)',/top,/left,box=0,chars=labelsize,margin=0,background=cgcolor('white')
hor,0.0,lines=2

; Stellar mass
x = alog10(ncat.NSA_ELPETRO_MASS/h^2) 
xtit = textoidl('Stellar Mass (log(M_{sun}))')
bin = 0.25
mkhist, x, pind, bin, /r90, /noerase, pos=pos3, yr=yr, xr=fr,$
	xtickint=0.2,xminor=4,ytickint=1.0,yminor=4,$
	xtit='',ytit='',charsize=titsize,xtickname=blanks,ytickname=blanks
ver,0.0,lines=2
axis,xaxis=1,xr=fr,/xs,xtit='Fraction',xtickint=0.2,xminor=4,charsize=titsize
; legend
al_legend,'(b)',/top,/right,box=0,chars=labelsize,margin=0,background=cgcolor('white')

; labels
plot,[0,1],/nodata,xs=4,ys=4,pos=pos4,/noerase
al_legend,['MaNGA Control','All AGNs','AGN Fraction'],$ ; ,'AGN Fraction'],$
	textcolors=cgcolor(['Black','Blue','Red']),$ ;,'Dark Green'])
	charsize=labelsize,/top,/right

theend:
device,/close
spawn,'gv '+psfile+' &'

end
