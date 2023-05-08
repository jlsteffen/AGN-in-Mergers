;+
; plot AGN fraction in the plane of Mass and redshift
; use manga_mkplot for 2D
; use mkhist for 1D
; This version combines 2D and 1D distributions in one plot
;-

function mymodel, x, p
	; X - dummy parameter to be compatible with MPFITFUN
	COMMON chi2_block, xx, yy, cid
	model = p[0]*exp(-0.5d*(yy-p[1])^2/p[2]^2)*(1.0+xx)^4
	return,model[cid]
end

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

; common block
;undefine,chi2_block,xx,yy,cid
COMMON chi2_block, xx, yy, cid

; setup plot
nx = 3
ny = 2
psfile = 'figs/pd2agn.eps'
setps,psfile,12*1.4,12*1.5,font='helvetica'
titsize = 1.5
labelsize = 1.5

pos1 = [0.12,0.4,0.65,1.0]
pos1a =[0.40,0.45,0.63,0.65]
pos2 = [0.12,0.1,0.65,0.4]
pos3 = [0.65,0.4,0.98,0.88]
pos4 = [0.65,0.1,0.98,0.4]

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
print,'Control Sample:',n_elements(ncat),n_elements(pind)
; 2513, 365

; original datapoints
h = 0.7 ; Hubble
;x = ncat.IFUrin*dangular(ncat.NSA_Z,/kpc)/(180.*3600./!pi) 
x = ncat.NSA_Z 
y = alog10(ncat.NSA_ELPETRO_MASS/h^2) 

blanks = replicate(' ',10)
; restore previously saved stuff
restore,'figs/pd2d_z_mass.sav'
mid = where(cden lt 10) ; pixels to mask out

; plot main panel
map = aden/cden
map[mid] = !values.f_nan
ctable = 62
c_map = cden
c_levels = max(c_map)*range(0.05,0.8,4)
print,c_levels
;c2_map = pden
;c2_levels = max(c2_map)*range(0.1,0.8,3)
zr = [-0.04,0.4]
manga_mkplot,map,zr=zr,/colorbar,zsize=labelsize,$ 
	ztit=textoidl('AGN Fraction'),$
	xr=xr,yr=yr,xtit='',$
	ytit='log(Stellar Mass/M'+sunsymbol()+')',$
	xtickint=0.05,xminor=5,ytickint=1,yminor=4,ztickint=0.1,zminor=4,$
	xtickname=blanks,$
	position=pos1,sauron=ctable,charsize=titsize,$
	c_map=c_map,c_levels=c_levels,/noerase ;,$
	;c2_map=c_map2,c2_levels=c_levels2
;oplot,x,y,psym=3,color=cgcolor('black')
al_legend,'(a)',/top,/left,box=0,chars=labelsize,margin=0,background=cgcolor('white')

; build model for AGN fraction
data = aden/cden
data[mid] = !values.f_nan
; try to get the best fit
u = aden*1.0
du = sqrt(aden)
v = cden*1.0
dv = sqrt(cden)
; compute ratio and its error
err = nden*0.0
err[*,*] = (eratio(u,du,v,dv))[*,1]
err[mid] = !values.f_nan

;; start fitting
;ind = where(aden gt 0 and cden ge 10)
;help,ind
;result = [0.,0,0,0,0]
;; fixed parameters
;z0  = 0.05
;pl  = 4.0
;for amp = 0.2,0.3,0.01 do begin
;	for sig = 0.4, 0.6, 0.01 do begin
;		for cen = 10.0, 11.0, 0.01 do begin
;	model = amp * exp(-0.5d*(yy-cen)^2/sig^2) * ((1+xx)/(1+z0))^pl
;	chisq = total(((data[ind]-model[ind])/err[ind])^2)
;	chisq2 = total((data[ind]-model[ind])^2) ; treat every pixel equally
;	result = [[result],[[amp,sig,cen,chisq,chisq2]]]
;		endfor
;	endfor
;endfor
;result = result[*,1:*]
;print,min(result[3,*],id1),result[*,id1]
;print,min(result[4,*],id2),result[*,id2]

; run MPFIT
cid = where(aden gt 0 and cden ge 10)
P0 = [0.15,10.5,0.5]
y = data[cid]
yerr = cid*0+mean(err,/nan)
P = mpfitfun('mymodel',x,y,yerr,P0,perror=perr,$
	status=status,bestnorm=bestnorm,dof=dof) 
PCERR = PERR * SQRT(BESTNORM / DOF)
print,p,perr,pcerr,bestnorm,dof
save,p,perr,pcerr,bestnorm,dof,filename='pd2agn.idlsav'
agnmod = p[0] * exp(-0.5d*(yy-p[1])^2/p[2]^2) * (1+xx)^4

; show model as inset
map = agnmod
map[mid] = !values.f_nan
manga_mkplot,map,zr=zr,xr=xr,yr=yr,xtit='',ytit='',$
	xtickint=0.05,xminor=5,ytickint=1,yminor=4,$
	position=pos1a,sauron=ctable,charsize=1.0,$
	;c_map=c_map,c_levels=c_levels,$
	/noerase,tit='Model'

; fraction range to plot
fr = [-0.02,0.35]
; histogram 1
loadct,0
; redshift
x = ncat.NSA_Z 
bin = 0.02/2
mkhist, x, pind, bin, /noerase, pos=pos2, xr=xr, yr=fr,$
	xtickint=0.05,xminor=5,ytickint=0.1,yminor=2,$
	xtit='Redshift',ytit='Fraction',charsize=titsize
al_legend,'(c)',/top,/left,box=0,chars=labelsize,margin=0 ;,background=cgcolor('white')
hor,0.0,lines=2
; compare w/ best-fit model
npx = (size(xx))[1]
a = xx[1:npx-2,0]
b = (total(agnmod*cden,2)/total(cden,2))[1:npx-2]
oplot,a,b,lines=2,color=cgcolor('dark green')
;oplot,rebin(a,n_elements(a)/2),rebin(b,n_elements(b)/2),lines=2,color=cgcolor('dark green')

; Stellar mass
x = alog10(ncat.NSA_ELPETRO_MASS/h^2) 
xtit = textoidl('Stellar Mass (log(M_{sun}))')
bin = 0.25
mkhist, x, pind, bin, /r90, /noerase, pos=pos3, yr=yr, xr=fr,$
	xtickint=0.1,xminor=2,ytickint=1.0,yminor=4,$
	xtit='',ytit='',charsize=titsize,xtickname=blanks,ytickname=blanks
ver,0.0,lines=2
axis,xaxis=1,xr=fr,/xs,xtit='Fraction',xtickint=0.1,xminor=2,charsize=titsize
; compare w/ best-fit model
frac = total(agnmod*cden,1,/nan)/total(cden,1,/nan)
s = where(finite(frac))
oplot,frac[s],yy[0,s],lines=2,color=cgcolor('dark green')
; legend
al_legend,'(b)',/top,/right,box=0,chars=labelsize,margin=0,background=cgcolor('white')

; labels
plot,[0,1],/nodata,xs=4,ys=4,pos=pos4,/noerase
al_legend,['Control Dist.','AGN Dist.','AGN Fraction','Model Fraction'],$
	textcolors=cgcolor(['Black','Blue','Red','Dark Green']),$
	colors=cgcolor(['Black','Blue','Red','Dark Green']),$
	lines=[0,0,0,2],linsize=0.5,$
	psym=[-3,-3,-6,-3],$
	charsize=1.0,/top,/left

theend:
undefine,chi2_block

device,/close
;spawn,'gv '+psfile+' &'

end
