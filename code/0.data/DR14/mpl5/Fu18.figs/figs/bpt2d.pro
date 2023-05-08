;+
; Generate diagnostic diagrams for BAGN candidates.
; This version makes images of different IFUs the same size
; It also plots the polygon-extracted spectra' measurements
;-

; load pair sample catalog
pcat = mrdfits('../7.pairsample.fits',1)
pcat.plateifu = strtrim(pcat.plateifu)
Rin = (pcat.IFUdiam+3.5-2.31)/2*sqrt(3)/2 ; inner circle radius

; choose a subclass
Mi_cut = -18.0
s = where(min(pcat.bclass,dim=1) ge 2 and $ ; 18
	(pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 18
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 16 (8146-12705 8978-6101) 
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 16
	pcat.sep_arcsec le Rin and $ ; 15 (8553-9102)
	pcat.NSA_ELPETRO_ABSMAG[5]+5*alog10(0.7) lt Mi_cut) ; 15
pcat = pcat[s]

outdir = 'figs/bpt2d/'
plateifus = pcat.plateifu 
infile = '$MANGA_DIR/hfdap/MPL-5/spfit/'+plateifus+'.fits'

; latex plotting commands
pdir = '../9.bagn/'
if n_elements(plateifus) le 8 then begin
	forprint,'\plotone{'+pdir+outdir+plateifus+'.eps}',$
	text=outdir+'part1.tex',/nocomment 
endif else begin
	forprint,'\plotone{'+pdir+outdir+plateifus[0:7]+'.eps}',$
	text=outdir+'part1.tex',/nocomment 
	forprint,'\plotone{'+pdir+outdir+plateifus[8:*]+'.eps}',$
	text=outdir+'part2.tex',/nocomment 
endelse
;stop

; set up output directory
if ~file_test(outdir) then spawn,'mkdir '+outdir $
	else spawn,'rm -f '+outdir+'*.eps'
; set up plotting window
nx = 7
ny = 1

; manga IFUs sizes
; https://trac.sdss.org/wiki/MANGA/Information
; 3.5" increased due to 3-point dithering
bundles = hash()
bundles[19] = 12.5 + 3.5 ; 19 fibers -> 12.5" diameter
bundles[37] = 17.5 + 3.5 
bundles[61] = 22.5 + 3.5
bundles[91] = 27.5 + 3.5
bundles[127] = 32.5 + 3.5

blank = replicate(' ',10)
for i=0,n_elements(infile)-1 do begin
;for i=11,11 do begin
	if i eq 0 or i eq 8 then cbar = 1 else cbar = 0

	; get plate-ifudsgn from filename
	fname = repstr(exfilename(infile[i]),'.fits','')
	; set up plot
	psfile = outdir+fname+'.eps'
	exwidth = 2.0
	if cbar eq 1 then $
		setps,psfile,8.5*nx+exwidth,8.5*1.2*ny,font='helvetica' $
		else $
		setps,psfile,8.5*nx+exwidth,8.5*ny,font='helvetica' 
	multiplot,/default
	pos = [0.001d,0.001,8.5*nx/(8.5*nx+exwidth),0.99]
	multiplot,[nx,ny],/dox,/doy,position=pos ;,/nomargin
	
	; load best-fit parameters, binmap, drp
	pars = mrdfits(infile[i],1,/silent)
	vbin = mrdfits(infile[i],2,/silent)
	drp = mrdfits(infile[i],3,/silent)

	; BPT classification
	linename = strtrim(pars[0].name,2)+strc(round(pars[0].lambda))
	io3 = where(linename eq 'OIII5008')
	ihb = where(linename eq 'Hb4863')
	in2 = where(linename eq 'NII6585')
	iha = where(linename eq 'Ha6565')
	io1 = where(linename eq 'OI6302')
	is2a = where(linename eq 'SII6718')
	is2b = where(linename eq 'SII6733')
	;o1ha = alog10(pars.flux[io1,0]/pars.flux[iha,0])
	n2ha = alog10(pars.flux[in2,0]/pars.flux[iha,0])
	;s2ha = alog10((pars.flux[is2a,0]+pars.flux[is2b,0])/pars.flux[iha,0])
	o3hb = alog10(pars.flux[io3,0]/pars.flux[ihb,0])

	; H-alpha EW, rest-frame
	; dV = z*c ~= ln(lambda1/lambda0)*c = log10(lambda1/lambda0)*ln(10)*c
	c = 299792.4580d ; Speed of light in km/s
	dV = abs(pars[0].log10lam-alog10(6564.632))*c*alog(10.)
	idx_con = where(dV gt 200. and dV lt 400.)
	c_ha = median(pars.best[idx_con]-pars.emis[idx_con],dim=1)
	w_ha = pars.flux[iha,0]/c_ha 
	;w_ha_err = w_ha * emlerr(pars.aon[iha],pars.sigma_obs[iha])
	whanclass = whan(w_ha,n2ha,code=whancode)

	binmap = vbin.binmap
	xdim = (size(binmap))[1]
	ydim = (size(binmap))[2]
	tmp = fltarr(xdim,ydim)
		
	; SDSS finder chart
	basename = strc(drp.plate)+'-'+drp.ifudsgn
	; matches the values used by sdss.pro
	scale = 0.25 
	npix = fix(40./scale)
	xr = [-0.5,0.5]*scale*npix
	; SDSS image w/o spec label
	jpgfile = '~/work/manga/cats/sdss/'+basename+'b.png'
	read_jpeg,jpgfile,sdss_img
	; make WCS header
	mkhdr,h_sdss,fltarr(npix,npix)
	putast,h_sdss,[[-1,0],[0,1]]*scale/3600d, [npix/2+0.5,npix/2-1.5], $
		[drp.ifura,drp.ifudec], ['RA---TAN','DEC--TAN']	
	; compute position
	position=!p.position
	r_ypos = position[3]-position[1]
	if cbar eq 1 then pos = position-[0,0,0,0.2*r_ypos] else pos=position
	; display image
	tvimage,sdss_img,true=1,position=pos
	; invisible axes
	plot,[1,1],[1,1],/nodata,xr=[0,npix],yr=[0,npix],xs=1+4,ys=1+4,pos=pos 
	; show NSA ELPETRO Ellipse
	adxy,h_sdss,pcat[i].objra,pcat[i].objdec,xc,yc
	;plots,xc,yc,psym=1,color=cgcolor('red'),syms=5,thick=2
	rmax = pcat[i].NSA_ELPETRO_TH50_R/scale ; Re - arcsec -> pixel
	rmin = rmax*pcat[i].NSA_ELPETRO_BA ; b/a ratio
	pos_ang = pcat[i].NSA_ELPETRO_PHI+90 ; PA (east of N)
	tvellipse,rmax,rmin,xc,yc,pos_ang,/data,color=cgcolor('red'),lines=2,noclip=0
	; show SDSS components
	adxy,h_sdss,pcat[i].ras,pcat[i].decs,xc,yc
	plots,xc[0]+0.5,yc[0]+0.5,psym=cgsymcat(6),syms=3,color=cgcolor('green')
	plots,xc[1]+0.5,yc[1]+0.5,psym=cgsymcat(9),syms=3,color=cgcolor('green')
	; plot arcsec-interval tickmarks
	plot,[1,1],[1,1],/nodata,xr=xr,yr=xr,/xs,/ys,/noerase,$
		xtickint=10,xminor=5,ytickint=10,yminor=5,$		
		position=pos,xtickname=blank,color=cgcolor('white')
	; draw IFU around target
	nfiber = drp.IFUDESIGNSIZE
	a = findgen(7)*(!pi*2/6.)
	oplot,bundles[nfiber]*0.5*cos(a),bundles[nfiber]*0.5*sin(a),$
		color=cgcolor('sky blue'),lines=2
	; label the plot
	xyouts,0,xr[1]*0.82,drp.plateifu+' z='+string(drp.nsa_z,f='(f5.3)'),$
		chars=1.7,color=cgcolor('white'),align=0.5	
	xyouts,xr[0]*0.9,xr[0]*0.93,string(pcat[i].sep_arcsec,f='(f4.1)')+$
		'"='+string(pcat[i].sep_kpc,f='(f4.1,"kpc")'),$
		chars=1.5,color=cgcolor('white'),align=0
	; plot black outline
	plot,[0,1],[0,1],/nodata,/xs,/ys,/noerase,$
		xtickint=10,xminor=-1,ytickint=10,yminor=-1,$		
		position=pos,xtickname=blank,color=cgcolor('black')
	multiplot,/dox

	; MaNGA binning map
	binmap = vbin.binmap
	; make WCS header
	mkhdr,h_manga,binmap
	putast,h_manga,[[sxpar(vbin.hdr,'CD1_1'),0],[0,sxpar(vbin.hdr,'CD2_2')]], $
		[sxpar(vbin.hdr,'CRPIX1'),sxpar(vbin.hdr,'CRPIX2')], $
		[sxpar(vbin.hdr,'CRVAL1'),sxpar(vbin.hdr,'CRVAL2')], $
		['RA---TAN','DEC--TAN']	
	; match to SDSS image
	hastrom,binmap,h_manga,binmap2,hdr,h_sdss
	; get axes ranges
	xr = [-0.5,0.5]*sxpar(hdr,'NAXIS1')*sxpar(hdr,'cd2_2')*3600
	yr = [-0.5,0.5]*sxpar(hdr,'NAXIS2')*sxpar(hdr,'cd2_2')*3600
	pxsize = (sxpar(hdr,'cd2_2')*3600)^2

	; match SDSS_img w/ cube WCS for contour plot
	sdss2 = total(sdss_img,1)
	c_map = sdss2
	c_levels = max(sdss2)*range(0.15,0.9,5)

	; H-alpha A/N & EQW mask
	ha_aon = manga_mkmap(binmap,pars.aon[iha])
	ha_ew  = manga_mkmap(binmap,w_ha)
	; 2D mask
	mask = ha_aon*0+!values.f_nan
	;mask[where(ha_aon gt 3.0 and ha_ew ge 3)] = 1.0
	mask[where(ha_aon gt 3.0)] = 1.0
	; 1D mask
	;m1d = where(pars.aon[iha] gt 3.0 and w_ha ge 3.0)
	m1d = where(pars.aon[iha] gt 3.0)
	; disable mask
	;mask[*,*] = 1.0
	;m1d = indgen(n_elements(pars.aon[iha]))

	; Ha velocity map
	ctable = 1 
	map = manga_mkmap(binmap,pars.vel[iha,0]) 
	hastrom,map*mask,h_manga,map2,h2,hdr,interp=0
	; aviod black color at maximum velocity
	s = where(map2 gt 395,ct1)
	if ct1 gt 0 then map2[s] = 395
	manga_mkplot,map2,zr=[-400,400],colorbar=cbar,$ 
		ztit='V (H'+cgSymbol('alpha')+') (km/s)',$
		xr=xr,yr=yr,xtit='',ytit='',$
		xtickint=10,xminor=5,ytickint=10,yminor=5,$
		position=!p.position,sauron=ctable,xtickname=blank,$
		ztickint=200,zminor=5,$		
		c_map=c_map,c_levels=c_levels
	oplot,bundles[nfiber]*0.5*cos(a),bundles[nfiber]*0.5*sin(a),$
		color=cgcolor('dark gray'),lines=2
	multiplot,/dox
	
	; H-alpha lineflux map
	ctable = 33
	map = manga_mkmap(binmap,alog10(pars.flux[iha,0]))
	hastrom,map*mask,h_manga,map2,h2,hdr,interp=0
	manga_mkplot,map2,zr=[-0.5,2.5],colorbar=cbar,$ 
		ztit='log(H'+cgSymbol('alpha')+textoidl(') (10^{-17} erg s^{-1} cm^{-2})'),$
		xr=xr,yr=yr,xtit='',ytit='',$
		xtickint=10,xminor=5,ytickint=10,yminor=5,$
		position=!p.position,sauron=ctable,xtickname=blank,$
		ztickint=1.0,zminor=5,$		
		c_map=c_map,c_levels=c_levels 
	oplot,bundles[nfiber]*0.5*cos(a),bundles[nfiber]*0.5*sin(a),$
		color=cgcolor('dark gray'),lines=2
	multiplot,/dox

	;; [O III] lineflux map
	;ctable = 10	
	;map = manga_mkmap(binmap,alog10(pars.flux[io3,0]))
	;hastrom,map*mask,h_manga,map2,h2,hdr,interp=0
	;manga_mkplot,map2,zr=[-2.0,2.0],colorbar=cbar,$
	;	ztit=textoidl('log([O III]) (10^{-17} erg s^{-1} cm^{-2})'),$
	;	xr=xr,yr=yr,xtit='',ytit='',$
	;	xtickint=10,xminor=5,ytickint=10,yminor=5,$
	;	position=!p.position,sauron=ctable,xtickname=blank,$
	;	ztickint=1.0,zminor=5,$		
	;	c_map=c_map,c_levels=c_levels 
	;oplot,bundles[nfiber]*0.5*cos(a),bundles[nfiber]*0.5*sin(a),$
	;	color=cgcolor('dark gray'),lines=2
	;multiplot,/dox

	; EW (Ha) map
	ctable = 33
	map = manga_mkmap(binmap,alog10(w_ha))
	hastrom,map*mask,h_manga,map2,h2,hdr,interp=0
	manga_mkplot,map2,zr=[-0.5,2.5],colorbar=cbar,$ 
		ztit='log(EW H'+cgSymbol('alpha')+') ('+cgSymbol('Angstrom')+')',$
		xr=xr,yr=yr,xtit='',ytit='',$
		xtickint=10,xminor=5,ytickint=10,yminor=5,$
		position=!p.position,sauron=ctable,xtickname=blank,$
		ztickint=1.0,zminor=5,$		
		c_map=c_map,c_levels=c_levels 		
	oplot,bundles[nfiber]*0.5*cos(a),bundles[nfiber]*0.5*sin(a),$
		color=cgcolor('dark gray'),lines=2
	multiplot,/dox

	; OIII/Hb map
	ctable = 38 ; ;38 ; 10
	map = manga_mkmap(binmap,o3hb) 
	hastrom,map*mask,h_manga,map2,h2,hdr,interp=0
	manga_mkplot,map2,zr=[-1.0,1.0],colorbar=cbar,$ 
		ztit='log([O III]/H'+cgsymbol('beta')+')',$
		xr=xr,yr=yr,xtit='',ytit='',$
		xtickint=10,xminor=5,ytickint=10,yminor=5,$
		position=!p.position,sauron=ctable,xtickname=blank,$
		ztickint=0.5,zminor=5,$		
		c_map=c_map,c_levels=c_levels 		
	oplot,bundles[nfiber]*0.5*cos(a),bundles[nfiber]*0.5*sin(a),$
		color=cgcolor('dark gray'),lines=2
	multiplot,/dox

	; NII/Ha BPT map
	ctable = 38
	zr = [1,4]+[-0.5,3.5]/4
	bptclass = bpt_k06_n2(o3hb,n2ha,code=bptcode) ; code = 1-4
	;ind = where(w_ha lt 3,ct)
	;if ct gt 0 then bptcode[ind] = 0 ; set to ambiguous
	map = manga_mkmap(binmap,bptcode) ;,hdr
	hastrom,map*mask,h_manga,map2,h2,hdr,interp=0
	manga_mkplot,map2,zr=zr,colorbar=cbar,$
		ztit=textoidl('[NII]/H\alpha BPT Class'),$
		xr=xr,yr=yr,xtit='',ytit='',$
		position=!p.position,sauron=ctable,$
		xtickint=10,xminor=5,ytickint=10,yminor=5,$
		ztickname=['SF','Comp','LINER','Sey'],$
		zminor=-1,ztickint=1,$
		xtickname=blank,c_map=c_map,c_levels=c_levels
	oplot,bundles[nfiber]*0.5*cos(a),bundles[nfiber]*0.5*sin(a),$
		color=cgcolor('dark gray'),lines=2
	multiplot,/dox,/doy
	; NII/Ha BPT diagram
	position=!p.position
	r_xpos = position[2]-position[0]
	r_ypos = position[3]-position[1]
	if cbar eq 1 then pos = position-[0,0,0,0.2*r_ypos] else pos=position
	xr=[-1.4,0.5]
	yr=[-1.2,1.4]
	if keyword_set(cbar) then begin
		plot,n2ha,o3hb,/nodata,psym=1,xr=xr,yr=yr,pos=pos,xs=1+8,ys=1+8,$
			xtickname=blank,ytickname=blank,color=255 
		axis,xaxis=1,xr=xr,/xs,xtit='log([N II]/H'+cgsymbol('alpha')+')',color=255
		axis,yaxis=1,yr=yr,/ys,ytit='log([O III]/H'+cgsymbol('beta')+')',color=255
	endif else begin
		plot,n2ha,o3hb,/nodata,psym=1,xr=xr,yr=yr,pos=pos,xs=1,ys=1,$
			xtickname=blank,ytickname=blank,color=255
	endelse
	for code=1,4 do begin
		; high EW points
		ind = where(bptcode[m1d] eq code and w_ha[m1d] ge 3,ct2)
		if ct2 gt 0 then begin
			color=bytscl(code*1.0,min=zr[0],max=zr[1],top=255)
			plotsym,0,/fill
			plots,n2ha[m1d[ind]],o3hb[m1d[ind]],psym=8,syms=1.0,color=color,noclip=0
		endif
		; low EW points
		ind = where(bptcode[m1d] eq code and w_ha[m1d] lt 3,ct2)
		if ct2 gt 0 then begin
			color=bytscl(code*1.0,min=zr[0],max=zr[1],top=255)
			plotsym,0
			plots,n2ha[m1d[ind]],o3hb[m1d[ind]],psym=8,syms=1.0,color=color,noclip=0
		endif
	endfor
	; load polygon components fit parameters
	x = pcat[i].n2ha
	y = pcat[i].o3hb
	dx = pcat[i].dn2ha
	dy = pcat[i].do3hb
	; show components line ratios & errors
	plots,x[0],y[0],psym=cgSymCat(15),syms=2,color=cgcolor('red')
	plots,x[1],y[1],psym=cgSymCat(16),syms=2,color=cgcolor('sky blue')
	plots,x[0],y[0],psym=cgSymCat(6),syms=2,color=255
	plots,x[1],y[1],psym=cgSymCat(9),syms=2,color=255
	oploterror,x,y,dx,dy,psym=3,/nohat

	; Kewley06 extreme SB line
	x = range(-2.5,1.5,40)
	oplot,x,1.19+0.61/((x<0.47)-0.47),lines=2,color=255 
	; Kauffmann03 pure SF line
	oplot,x,1.30+0.61/((x<0.05)-0.05),lines=0,color=255 
	;; Kauffmann03 Seyfert/Liner line 
	;x = range(-0.2,1.5,40)
	;oplot,x,2.145*x+0.465,color=255,lines=2
	; Schawinski07 Seyfert/Liner line
	x = range(-0.18,1.5,40)
	oplot,x,1.05*x+0.45,lines=3,color=255

	multiplot,/dox

	device,/close
	;spawn,'gv '+psfile+' &'
endfor

end
