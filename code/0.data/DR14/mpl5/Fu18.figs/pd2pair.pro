;+
; plot Close Pair fraction in the plane of Mass and Rifu
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
psfile = 'figs/pd2pair.eps'
setps,psfile,12*1.4,12*1.5,font='helvetica'
titsize = 1.5
labelsize = 1.5

pos1 = [0.12,0.4,0.65,1.0]
pos1a =[0.40,0.45,0.63,0.65]
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
help,where(pcat.NSA_NSAID ne ncat[pind].NSA_NSAID)
help,where(strtrim(ncat[pind].plateifu) ne strtrim(pcat.plateifu))

blanks = replicate(' ',10)
; restore previously computed 2D histograms from pd2d.pro
restore,'figs/pd2d_Rifu_mass.sav'
mid = where(nden lt 4) ; pixels to mask out
cid = where(pden gt 0 and nden ge 4) ; pixels to include in chisq calculation
data = pden/nden
data[mid] = !values.f_nan
; error of pair fraction - Poisson noise
u = pden*1.0
u[where(nden gt 0 and pden eq 0)] = 1.0
du = sqrt(pden)
v = nden*1.0
dv = sqrt(nden)
err = data*0
err[*,*] = (eratio(u,du,v,dv))[*,1]
err[mid] = !values.f_nan

; plot main panel
zr = [-0.04,0.4]
map = pden/nden
map[mid] = !values.f_nan
ctable = 62
c_map = nden
c_levels = max(c_map)*range(0.05,0.8,4)
print,c_levels
;c2_map = pden
;c2_levels = max(c2_map)*range(0.1,0.8,3)
manga_mkplot,map,zr=zr,/colorbar,zsize=labelsize,$ 
	ztit=textoidl('Pair Fraction'),$
	xr=xr,yr=yr,xtit='',$
	ytit = 'log(Stellar Mass/M'+sunsymbol()+')',$
	xtickint=10,xminor=4,ytickint=1,yminor=4,ztickint=0.1,zminor=4,$
	xtickname=blanks,$
	position=pos1,sauron=ctable,charsize=titsize,$
	c_map=c_map,c_levels=c_levels,/noerase ;,$
	;c2_map=c_map2,c2_levels=c_levels2
al_legend,'(a)',/top,/left,box=0,chars=labelsize,margin=0,background=cgcolor('white')

; find best-fit scaling factor for Hopkins model
for f = 2.45, 3.20, 0.01 do begin
	chisq = total(((data[cid]-f*h10mod[cid])/err[cid])^2)
	print,f,chisq
endfor
;     2.45000      31.8765
;     2.83000      30.8277
;     3.20000      31.8194
; best-fit Hopkins model 
;h10mod *= 2.83

; single powerlaw model
y0 = 11.0
pl = 1.0
plmod = xx/30d * 10d^((yy-y0)*pl)
for amp = 0.08, 0.107, 0.001 do begin
	chisq = total(((data[cid]-amp*plmod[cid])/err[cid])^2)
	print,amp,pl,chisq
	;plmod[mid] = !values.f_nan
	;atv,[data,amp*plmod,abs(data-amp*plmod)],/linear,min=-0.1,max=0.5
endfor
;    0.0810000      1.00000       30.491822
;    0.0930000      1.00000       29.476770
;     0.106000      1.00000       30.570368
amp = 0.093
h10mod = amp*plmod

;; mass-independent single powerlaw model
;y0 = 11.0
;plmod = (xx < 20.)/30d
;for amp = 0.05, 0.30, 0.01 do begin
;	chisq = total(((data[cid]-amp*plmod[cid])/err[cid])^2)
;	print,amp,chisq
;	;plmod[mid] = !values.f_nan
;	;atv,[data,amp*plmod,abs(data-amp*plmod)],/linear,min=-0.1,max=0.5
;endfor
;;     0.100000       34.671201 
;amp = 0.1
;h10mod = amp*plmod

;; double powerlaw model
;y0 = 11.0
;p1 = 0.2
;p2 = 1.7
;plmod = xx/30d * (10d^((yy-y0)*p1) + 10d^((yy-y0)*p2))
;  amp = total(data[s])/total(plmod[s])
;print,amp
;h10mod = amp*plmod

; show model as inset
map = h10mod
map[mid] = !values.f_nan
manga_mkplot,map,zr=zr,$ 
	xr=xr,yr=yr,xtit='',ytit='',$
	xtickint=10,xminor=4,ytickint=1,yminor=4,$
	position=pos1a,sauron=ctable,charsize=1.0,$
	;c_map=c_map,c_levels=c_levels,$
	/noerase,tit='Model'

; fraction range to plot
fr = [-0.02,0.43]
; histogram 1
loadct,0
; IFU size in kpc
x = ncat.IFUrin*dangular(ncat.NSA_Z,/kpc)/(180.*3600./!pi) ; IFU inner radius in kpc
bin = 2.5
mkhist, x, pind, bin, /noerase, pos=pos2, xr=xr, yr=fr,$
	xtickint=10,xminor=4,ytickint=0.2,yminor=4,$
	xtit='IFU Radius (kpc)',ytit='Fraction',charsize=titsize
hor,0.0,lines=2
; compare w/ best-fit Hopkins10 model
a = xx[*,0]
b = total(h10mod*nden,2,/nan)/total(nden,2,/nan)
oplot,a,b,lines=2,color=cgcolor('dark green')
; legend 
al_legend,'(c)',/top,/left,box=0,chars=labelsize,margin=0,background=cgcolor('white')

; Stellar mass
x = alog10(ncat.NSA_ELPETRO_MASS/h^2) 
bin = 0.25
mkhist, x, pind, bin, /r90, /noerase, pos=pos3, yr=yr, xr=fr,$
	xtickint=0.2,xminor=4,ytickint=1.0,yminor=4,$
	xtit='',ytit='',charsize=titsize,xtickname=blanks,ytickname=blanks
ver,0.0,lines=2
axis,xaxis=1,xr=fr,/xs,xtit='Fraction',xtickint=0.2,xminor=4,charsize=titsize
; compare w/ best-fit Hopkins10 model
frac = total(h10mod*nden,1,/nan)/total(nden,1,/nan)
s = where(finite(frac))
oplot,frac[s],yy[0,s],lines=2,color=cgcolor('dark green')
; legend
al_legend,'(b)',/top,/right,box=0,chars=labelsize,margin=0,background=cgcolor('white')

; labels
plot,[0,1],/nodata,xs=4,ys=4,pos=pos4,/noerase
al_legend,['All MaNGA','Close Pairs','Pair Fraction','Model Fraction'],$
	textcolors=cgcolor(['Black','Blue','Red','Dark Green']),$
	colors=cgcolor(['Black','Blue','Red','Dark Green']),$
	lines=[0,0,0,2],linsize=0.5,$
	psym=[-3,-3,-6,-3],$
	charsize=1.0,/top,/left

theend:
device,/close
;spawn,'gv '+psfile+' &'

end
