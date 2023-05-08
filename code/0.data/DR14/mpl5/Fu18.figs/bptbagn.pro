;+
; NII/Ha BPT diagram 
; WHAN diagram
;-

PRO mkplot, x,y,bc, px,py,dpx,dpy,pbc, _REF_EXTRA=extra
	
	; plot frame
	plot,[0,1],[0,1],/nodata,/xs,/ys,_extra=extra

	; plot all DR14
	;s = where(bc eq 0) ; retired galaxies
	;plots,x[s],y[s],psym=cgSymCat(16),$
	;	color=cgcolor('gray'),noclip=0,syms=0.5
	s = where(bc ge 1)
	plots,x[s],y[s],psym=cgSymCat(16),$
		color=cgcolor('gray'),noclip=0,syms=0.5
		
	;; BAGNs
	;syms = 1.5	
	;s = where(min(pbc,dim=1) ge 2)
	;;for i=0,n_elements(s)-1 do $ ; connecting lines
	;;	oplot,px[*,s[i]],py[*,s[i]],lines=0,color=cgcolor('green')
	;plots,px[0,s],py[0,s],psym=cgSymCat(15),noclip=0,syms=syms,color=cgcolor('red')
	;plots,px[0,s],py[0,s],psym=cgSymCat(6),noclip=0,syms=syms
	;plots,px[1,s],py[1,s],psym=cgSymCat(16),noclip=0,syms=syms,color=cgcolor('sky blue')
	;plots,px[1,s],py[1,s],psym=cgSymCat(9),noclip=0,syms=syms
	;oploterror,px[*,s],py[*,s],dpx[*,s],dpy[*,s],psym=3,/nohat

	; all other pairs
	syms = 1.0
	s = where(pbc ge 1)
	px2 = px[s]
	py2 = py[s]
	plots,px2,py2,psym=cgSymCat(16),noclip=0,syms=syms,color=cgcolor('dark gray')
	plots,px2,py2,psym=cgSymCat(9),noclip=0,syms=syms,color=cgcolor('black')
	; add error bars
	;oploterror,px[*,s],py[*,s],dpx[*,s],dpy[*,s],psym=3,/nohat
	
	; BAGNs
	syms = 1.5
	s = where(min(pbc,dim=1) ge 2)
	for i=0,n_elements(s)-1 do begin ; connecting lines
		oplot,px[*,s[i]],py[*,s[i]],lines=0,color=cgcolor('red')
		;print,i,px[0,s[i]],py[0,s[i]],px[1,s[i]],py[1,s[i]]
	endfor
	px2 = px[*,s]
	py2 = py[*,s]
	pbc2 = pbc[*,s]
	dpx2 = dpx[*,s]
	dpy2 = dpy[*,s]
	i = where(pbc2 eq 2)
	plots,px2[i],py2[i],psym=cgSymCat(15),noclip=0,syms=syms,color=cgcolor('sky blue')
	plots,px2[i],py2[i],psym=cgSymCat(6),noclip=0,syms=syms
	i = where(pbc2 eq 3)
	plots,px2[i],py2[i],psym=cgSymCat(15),noclip=0,syms=syms,color=cgcolor('green')
	plots,px2[i],py2[i],psym=cgSymCat(6),noclip=0,syms=syms
	i = where(pbc2 ge 4)
	plots,px2[i],py2[i],psym=cgSymCat(15),noclip=0,syms=syms,color=cgcolor('orange')
	plots,px2[i],py2[i],psym=cgSymCat(6),noclip=0,syms=syms
	; add error bars
	oploterror,px[*,s],py[*,s],dpx[*,s],dpy[*,s],psym=3,/nohat

end

; load galaxy pair catalog
pcat = mrdfits('../7.sample/pair1kpc.fits',1)
pcat.plateifu = strtrim(pcat.plateifu) ; 176
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 149
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 119
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 113
	pcat.sep_arcsec le pcat.IFUrin) ; 105
pcat = pcat[s]
; load all DR14 catalog
ncat = mrdfits('../7.sample/nuc1kpc.fits',1) ; 2718 
s = where((ncat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0) ; 2618
ncat = ncat[s]

; setup plot
loadct,0
psfile = 'figs/bptbagn.eps'
nx = 3
ny = 1
setps,psfile,15*1.1*nx,15*ny,font='helvetica'
multiplot,/default
multiplot,[nx,ny],xgap=0.03,/dox,/doy
titsize = 1.5
labelsize = 1.2

; [N II] BPT diagram
xr=[-1.3,0.6]
yr=[-1.2,1.4]
xtit='log([N II]/H'+cgSymbol('alpha')+')'
ytit='log([O III]/H'+cgSymbol('beta')+')'
x = ncat.n2ha
y = ncat.o3hb
bc = ncat.bclass
px = pcat.n2ha
py = pcat.o3hb
dpx = pcat.dn2ha
dpy = pcat.do3hb
pbc = pcat.bclass
mkplot,x,y,bc, px,py,dpx,dpy,pbc,xr=xr,yr=yr,$
	xtit=xtit,ytit=ytit,charsize=titsize
; Kewley06 extreme SB line
x = range(-2.5,1.5,40)
oplot,x,1.19+0.61/((x<0.47)-0.47),lines=2 
; Kauffmann03 pure SF line
oplot,x,1.30+0.61/((x<0.05)-0.05),lines=0 
; Schawinski07 Seyfert/Liner line
x = range(-0.18,1.5,40)
oplot,x,1.05*x+0.45,lines=3
; labels
xyouts,-0.02,-1.0,'Composite',align=0.5,/data,charsize=labelsize
xyouts,-0.30,1.05,'Seyfert',align=1.0,/data,charsize=labelsize
xyouts,0.53,0.0,'LINER',align=1.0,/data,charsize=labelsize
xyouts,-1.2,-0.5,'Starforming',align=0,/data,charsize=labelsize
; Reddening arrow
; for A(V)=5 mag => E(B-V)=A(V)/R, R=3.1
lambda = [6585.28,6564.61,5008.240,4862.683] 	
ccm_UNRED,lambda,[1.0,1.0,1.0,1.0],5/3.1,funred,R_V=3.1
px = -1.1 ; xr[1]-0.70*(xr[1]-xr[0])-0.05 
py = -1.0 ; yr[1]-0.9*(yr[1]-yr[0])-0.1
ax = alog10(funred[0]/funred[1])  ; arrow x 
ay = alog10(funred[2]/funred[3])  ; arrow y
l = sqrt(ax^2+ay^2)               ; length of the arrow
x1 = px+ax & y1 = py+ay
rsin = 0.5*sin(22.5*!dtor)
rcos = 0.5*cos(22.5*!dtor)
; 0.5 : Length of arrowhead legs relative to vectorlength.
; 22.5 deg = 1/2 times the angle between the arrowhead legs.
plots,[px,x1,x1-(ax*rcos+ay*rsin), x1,x1-(ax*rcos-ay*rsin)], $
      [py,y1,y1-(ay*rcos-ax*rsin), y1,y1-(ay*rcos+ax*rsin)], THICK=4.0
xyouts,noclip=0,px-0.00,py+0.02,'A!DV!N = 5',align=0.5,chars=1.5,chart=5
; legend
al_legend,'(a)',/top,/right,box=0,chars=2,margin=0
al_legend,['EW(H'+cgSymbol('alpha')+')> 3'+cgSymbol('angstrom'),$
	'All MaNGA','Close Pairs','Binary AGNs'],$
	psym=[3,16,16,6],color=cgcolor(['white','gray','dark gray','black']),$
	/top,/left,background=cgcolor('white')
multiplot,/dox,/doy

; [S II] BPT diagram
xr=[-1.1,0.6]
yr=[-1.2,1.4]
xtit='log([S II]/H'+cgSymbol('alpha')+')'
ytit='log([O III]/H'+cgSymbol('beta')+')'
x = ncat.s2ha
y = ncat.o3hb
bc = ncat.bclass
px = pcat.s2ha
py = pcat.o3hb
dpx = pcat.ds2ha
dpy = pcat.do3hb
pbc = pcat.bclass
mkplot,x,y,bc, px,py,dpx,dpy,pbc,xr=xr,yr=yr,$
	xtit=xtit,ytit=ytit,charsize=titsize
; Kewley06 [S II] classification lines
x = range(xr[0],xr[1],40)
oplot,x,1.30+0.72/((x<0.27)-0.32),lines=2 
oplot,x>(-0.4),1.89*(x>(-0.4))+0.76,lines=3 
; labels
xyouts,-0.6,0.95,'Seyfert',align=1.0,/data,charsize=labelsize
xyouts,0.53,0.6,'LINER',align=1.0,/data,charsize=labelsize
xyouts,-0.1,-1.00+0.05,'Starforming',align=1,/data,charsize=labelsize
xyouts,-0.1,-1.13+0.05,'& Composite',align=1,/data,charsize=labelsize
; Reddening arrow
; for A(V)=5 mag => E(B-V)=A(V)/R, R=3.1
lambda = [6726.0,6564.61,5008.240,4862.683] 	
ccm_UNRED,lambda,[1.0,1.0,1.0,1.0],5/3.1,funred,R_V=3.1
px = -0.9 ; xr[1]-0.70*(xr[1]-xr[0])-0.05 
py = -1.0 ; yr[1]-0.9*(yr[1]-yr[0])-0.1
ax = alog10(funred[0]/funred[1])  ; arrow x 
ay = alog10(funred[2]/funred[3])  ; arrow y
l = sqrt(ax^2+ay^2)               ; length of the arrow
x1 = px+ax & y1 = py+ay
rsin = 0.5*sin(22.5*!dtor)
rcos = 0.5*cos(22.5*!dtor)
; 0.5 : Length of arrowhead legs relative to vectorlength.
; 22.5 deg = 1/2 times the angle between the arrowhead legs.
plots,[px,x1,x1-(ax*rcos+ay*rsin), x1,x1-(ax*rcos-ay*rsin)], $
      [py,y1,y1-(ay*rcos-ax*rsin), y1,y1-(ay*rcos+ax*rsin)], THICK=4.0
xyouts,noclip=0,px-0.00,py+0.02,'A!DV!N = 5',align=0.5,chars=1.5,chart=5
; legend
al_legend,'(b)',/top,/right,box=0,chars=2,margin=0
multiplot,/dox,/doy

; [O I] diagram
xr = [-2.5,0.0]
yr = [-1.2,1.4]
xtit='log([O I]/H'+cgSymbol('alpha')+')'
ytit='log([O III]/H'+cgSymbol('beta')+')'
x = ncat.o1ha
y = ncat.o3hb
bc = ncat.bclass
px = pcat.o1ha
py = pcat.o3hb
dpx = pcat.do1ha
dpy = pcat.do3hb
pbc = pcat.bclass
mkplot,x,y,bc, px,py,dpx,dpy,pbc,xr=xr,yr=yr,$
	xtit=xtit,ytit=ytit,charsize=titsize
; Kewley06 classification lines
x = range(xr[0],xr[1],40)
oplot,x,1.33+0.73/((x<(-0.64))+0.59),lines=2
oplot,x>(-1.1),1.18*(x>(-1.1))+1.30,lines=3
; labels
xyouts,-1.6,0.95,'Seyfert',align=1.0,/data,charsize=labelsize
xyouts,-0.5,-0.1,'LINER',align=0.0,/data,charsize=labelsize
xyouts,-1.0,-1.00+0.05,'Starforming',align=1,/data,charsize=labelsize
xyouts,-1.0,-1.13+0.05,'& Composite',align=1,/data,charsize=labelsize
; Reddening arrow
; for A(V)=5 mag => E(B-V)=A(V)/R, R=3.1
lambda = [6300.0,6564.61,5008.240,4862.683] 	
ccm_UNRED,lambda,[1.0,1.0,1.0,1.0],5/3.1,funred,R_V=3.1
px = -2.0  
py = -1.0 
ax = alog10(funred[0]/funred[1])  ; arrow x 
ay = alog10(funred[2]/funred[3])  ; arrow y
l = sqrt(ax^2+ay^2)               ; length of the arrow
x1 = px+ax & y1 = py+ay
rsin = 0.5*sin(22.5*!dtor)
rcos = 0.5*cos(22.5*!dtor)
; 0.5 : Length of arrowhead legs relative to vectorlength.
; 22.5 deg = 1/2 times the angle between the arrowhead legs.
plots,[px,x1,x1-(ax*rcos+ay*rsin), x1,x1-(ax*rcos-ay*rsin)], $
      [py,y1,y1-(ay*rcos-ax*rsin), y1,y1-(ay*rcos+ax*rsin)], THICK=4.0
xyouts,noclip=0,px-0.00,py+0.02,'A!DV!N = 5',align=0.5,chars=1.5,chart=5
; legend
al_legend,'(c)',/top,/right,box=0,chars=2,margin=0

device,/close
spawn,'gv '+psfile+' &'

end
