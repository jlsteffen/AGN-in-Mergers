;+
; In 4.spfit we extracted and fit the nuclear spectra for all
; merging components in the merger sample, here we assemble useful observed 
; and derived parameters for all spec-confirmed mergers
;-

FOR JJ=2,2 DO BEGIN

CASE jj OF
	1: begin
		indir = 'pair1arc'
		apradi = 1.0
		apunit = 'arc'
		end
	2: begin
		indir = 'pair1kpc'
		apradi = 1.3 
		apunit = 'kpc'
		end
	3: begin
		indir = 'pair2arc'
		apradi = 2.0 
		apunit = 'arc'
		end
ENDCASE

; load DRPALL
drpall = mrdfits('../../sample/matched/drpall.fits',1) ; 2718
drpall.plateifu = strtrim(drpall.plateifu,2)
; limit to spec-confirmed mergers
tmp = mrdfits('../6.specclass/pair.fits',1) 
ind = where(strtrim(tmp.class,2) eq 'Merger',nobj)
drpall = drpall[ind]

; Pair-specific tags to be combined with DRPALL 
tags = ['sep_arcsec','sep_kpc','aprad',$
	'BPT','WHAN','BCLASS',$ ; classification info
	'Lambda','AoN','Sigma_Obs','AoN2','AoN3',$ ; line info
	'RAs','Decs','Rmag','Rpmag','Rpmag_2',$ ; model/PSF/Petro mags
	'FO3','dFO3','Lo3',$ ; 1d-17 erg/s/cm2, df/f, erg/s
	'O3Hb','dO3Hb','N2Ha','dN2Ha','WHa','dWHa',$
	'O1Ha','dO1Ha','S2Ha','dS2Ha',$
	'Vgas','Sgas','Vstar','Sstar','Mstar',$
	'dVgas','dSgas','dVstar','dSstar'] 
vals = ['0.','0.','0.',$
	replicate('[-9,-9]',3),$
	'fltarr(7)',replicate('fltarr(7,2)',4),$
	replicate('dblarr(2)',n_elements(tags)-11)] 
str = mrd_struct(tags,vals,nobj)

fmt = '(" Progress:",f6.2,"% done",$,%"\r")'
maxsteps = nobj
for i=0,nobj-1 do begin
	; components SDSS catalog	
	comps0 = mrdfits('../3.comps/'+drpall[i].plateifu+'.fits',1,/silent)
	
	; original spfit parameters - take BLR results if exist
	infile = '../4.spfit/'+indir+'/spfit/'+drpall[i].plateifu+'.fits'
	brfile = '../4.spfit/'+indir+'/spfit_blr/'+drpall[i].plateifu+'-br.fits'
	if file_test(brfile) then infile = brfile
	pars0 = mrdfits(infile,1,/silent)

	; sort comps & pars base on distance to manga target
	; primary = manga target = index 0
	dis = sphdist(comps0.raj2000,comps0.dej2000,drpall[i].objra,drpall[i].objdec,/deg)*3600
	idx = sort(dis)
	pars = pars0[idx]
	comps = comps0[idx]

	; get emission line indices 
	if i eq 0 then begin
		linename = strtrim(pars[0].name,2)+strc(round(pars[0].lambda))
		io3 = where(linename eq 'OIII5008')
		ihb = where(linename eq 'Hb4863')
		in2 = where(linename eq 'NII6585')
		iha = where(linename eq 'Ha6565')
		io1 = where(linename eq 'OI6302')
		is2a = where(linename eq 'SII6718')
		is2b = where(linename eq 'SII6733')
	endif

	; note that fit_manga use mean to combine spectra in each bin
	; find out the number of spaxels inside the aperture
	if apunit eq 'arc' then pixscale = 1.0 ; arcsec/pixel
	if apunit eq 'kpc' then $
		pixscale = dangular(drpall[i].NSA_Z,/kpc)/(180.*3600./!pi) ; kpc/pixel
	rad = apradi/pixscale ; pixel 
	npix = !pi*rad^2   ; aperture area in pixels
	str[i].aprad = rad ; pixel

	; emission line fluxes
	fn2 = pars.flux[in2,0]
	fha = pars.flux[iha,0]
	fo3 = pars.flux[io3,0]
	fhb = pars.flux[ihb,0]
	fo1 = pars.flux[io1,0]
	fs2a= pars.flux[is2a,0]
	fs2b= pars.flux[is2b,0]
	; fractional errors, df/f
	dfn2 = emlerr(pars.aon[in2],pars.sigma_obs[in2])
	dfha = emlerr(pars.aon[iha],pars.sigma_obs[iha])
	dfo3 = emlerr(pars.aon[io3],pars.sigma_obs[io3])
	dfhb = emlerr(pars.aon[ihb],pars.sigma_obs[ihb])
	dfo1 = emlerr(pars.aon[io1],pars.sigma_obs[io1])
	dfs2a= emlerr(pars.aon[is2a],pars.sigma_obs[is2a])
	dfs2b= emlerr(pars.aon[is2b],pars.sigma_obs[is2b])
	dfs2 = sqrt((dfs2a*fs2a)^2+(dfs2b*fs2b)^2)/(fs2a+fs2b)

	; Amplitude-to-Noise ratio 1
	; amplitude from best-fit Gaussian
	; noise from robust_sigma of residual (obs - best-fit stellar&emission line)
	idx = [in2,iha,io3,ihb,io1,is2a,is2b]
	str[i].lambda = pars[0].lambda[idx]
	str[i].aon = pars.aon[idx] 
	str[i].sigma_obs = pars.sigma_obs[idx]
	c_kms = 299792.46d ; speed-of-light
	for k=0,1 do begin
		; Amplitude-to-Noise ratio 2
		; amplitude from best-fit Gaussian
		; noise from input error array
		s = where(pars[k].good eq 1)
		; dV = z*c ~= ln(lambda1/lambda0)*c = log10(lambda1/lambda0)*ln(10)*c
		str[i].aon2[*,k] = $
			interpol(pars[k].emis[s]/pars[k].err[s],pars[k].log10lam[s],$
			alog10(pars[k].lambda[idx])+pars[k].vel[idx]/c_kms/alog(10))
		; Amplitude-to-Noise ratio 3
		; amplitude from residual spectrum (obs - best-fit stellar)
		; noise from input error array
		str[i].aon3[*,k] = $
			interpol((pars[k].galaxy-(pars[k].best-pars[k].emis))[s]/pars[k].err[s],$
			pars[k].log10lam[s],alog10(pars[k].lambda[idx])+pars[k].vel[idx]/c_kms/alog(10))
	endfor
	; line ratios
	str[i].n2ha = alog10(fn2/fha)
	str[i].o3hb = alog10(fo3/fhb)
	str[i].o1ha = alog10(fo1/fha)
	str[i].s2ha = alog10((fs2a+fs2b)/fha)
	; error propagation: d(log a) = (1/ln10) (da/a)
	str[i].dn2ha = sqrt(dfn2^2+dfha^2)/alog(10)
	str[i].do3hb = sqrt(dfo3^2+dfhb^2)/alog(10)
	str[i].do1ha = sqrt(dfo1^2+dfha^2)/alog(10)
	str[i].ds2ha = sqrt(dfs2^2+dfha^2)/alog(10)
	; NII/Ha BPT classification
	bptclass = bpt_k06_n2(str[i].o3hb,str[i].n2ha,code=bptcode)
	str[i].bpt = bptcode

	; WHAN classification
	;	CODE - integer, sAGN=4, wAGN=3, SF=1, RG=0, Invalid=-1
	; H-alpha EW, rest-frame
	; dV = z*c ~= ln(lambda1/lambda0)*c = log10(lambda1/lambda0)*ln(10)*c
	c = 299792.4580d ; Speed of light in km/s
	dV = abs(pars[0].log10lam-alog10(6564.632))*c*alog(10.)
	idx_con = where(dV gt 200. and dV lt 400.)
	c_ha = median(pars.best[idx_con]-pars.emis[idx_con],dim=1)
	str[i].wha = pars.flux[iha,0]/c_ha 
	str[i].dwha = str[i].wha * emlerr(pars.aon[iha],pars.sigma_obs[iha])
	whanclass = whan(str[i].wha,str[i].n2ha,code=whancode)
	str[i].whan = whancode

	; stellar mass
	str[i].Mstar = alog10(total(pars.m_star[*,0],1)*npix)
	; [O III] fluxes - 1d-17 erg/s/cm2
	str[i].fo3 = fo3*npix
	str[i].dfo3 = dfo3
	; 1d-17 erg/s/cm2 * Mpc^2 -> erg/s
	; 1d-17 * 4 * !PI * (3.086d24)^2 = 1.1967453e+33
	str[i].Lo3 = fo3*npix*(lumdist(drpall[i].nsa_z,/silent))^2*1.197d33
	
	; kinematics
	str[i].Vgas = pars.vel[iha,0]
	str[i].Sgas = pars.sigma[iha,0]
	str[i].Vstar = pars.kinstar[0,0]
	str[i].Sstar = pars.kinstar[1,0]
	; errors
	str[i].dVgas = pars.vel[0,1]
	str[i].dVstar = pars.kinstar[0,1]
	str[i].dSstar = pars.kinstar[1,1]
	; sig_obs d sig_obs = sig_int d sig_int
	; => d sig_int = (sig_obs^2/sig_int) * (d sig_obs / sig_obs)
	velscale = 69.0 ; km/s per pixel
	str[i].dSgas = (pars.sigma_obs[iha]*velscale)^2/pars.sigma[iha,0] * $
			emlerr(pars.aon[iha],pars.sigma_obs[iha],/sigma)

	; projected separations
	str[i].sep_arcsec = sphdist(comps[0].raj2000,comps[0].dej2000,$
		comps[1].raj2000,comps[1].dej2000,/deg)*3600d ; arcsec
	angscale = dangular(drpall[i].nsa_z,/kpc)/206264.81d ; kpc/arcsec
	str[i].sep_kpc = angscale * str[i].sep_arcsec ; kpc
	str[i].RAs = comps.raj2000
	str[i].Decs = comps.dej2000
	str[i].rmag = comps.rmag
	str[i].rpmag = comps.rpmag
	str[i].rpmag_2 = comps.rpmag_2

	; print progress
	if (i mod round(maxsteps/10d0) eq 0) then print,f=fmt,100.*i/maxsteps
endfor
print, '' ;; don't overwrite the final line

; BPT	CODE - Seyfert=4, LINER=3, Comp=2, SF=1, Ambi=0, Invalid=-1
; WHAN 	CODE - sAGN=4, wAGN=3, SF=1, RG=0, Invalid=-1 (WHa < 3 = RG, WHa > 6 = sAGN)
;;;;;;;;;;;;
; BCLASS CODES:
; -3/-2/-1: unreliable or invalid classifications
; 0/1/2/3/4/5: RG/SF/COMP/LINER/SEYFERT/BLAGN
;;;;;;;;;;;;

; Known BLAGNs: 16 cases
; these all have poor fits because SPFIT cannot properly fit BLAGNs 
; give those w/ bclass = 5
readcol,'../../sample/blagn/t1agn.txt',plates,ifus,f='l,l',comment='#'
basename = strc(plates)+'-'+strc(ifus)
match2,basename,drpall.plateifu,sa,sb
idx = where(sa ne -1,ct)
if ct gt 0 then str[sa[idx]].bclass[0] = 5

; Additional BLAGN
i = where(drpall.plateifu eq '8711-12701')
str[i].bclass[1] = 5

for i=0,1 do begin
	; W(Ha) < 3 A
	; Retired Galaxies: bclass = 0
	s = where(str.wha[i] lt 3 and str.bclass[i] ne 5,ct) 
	if ct gt 0 then str[s].bclass[i] = 0
	; W(Ha) > 3 A
	; (1) retain BPT classification
	s = where(str.wha[i] ge 3 and str.bclass[i] ne 5,ct)
	if ct gt 0 then str[s].bclass[i] = str[s].bpt[i]
	; (2) Ha AoN < 4 -> set bclass = -2 for those 
	s = where(str.wha[i] ge 3 and str.bclass[i] ne 5 and $
		str.bpt[i] ne -1 and str.aon[1,i] lt 4,ct)
	if ct gt 0 then str[s].bclass[i] = -2
	; (3) when the target fell outside of IFU 
	offset = sphdist(drpall.OBJRA,drpall.OBJDEC,drpall.IFURA,drpall.IFUDEC,/deg)*3600
	s = where(offset gt drpall.IFURAD,ct)
	if ct gt 0 then	str[s].bclass = -3
	; (4) use WHAN class for the two sources at z ~ 0.114
	s = where(drpall.nsa_z ge 0.113 and drpall.nsa_z lt 0.1155 $
		and str.whan[i] gt 0,ct)
	if ct gt 0 then begin
		forprint,drpall[s].plateifu,str[s].whan[i]
		str[s].bclass[i] = str[s].whan[i]-1
	endif
endfor

; save combined structure
mwrfits,struct_combine(str,drpall),indir+'.fits',/create

; save # of each class in a log file
pcat = mrdfits(indir+'.fits',1) ; 156 
;print,minmax(abs(pcat.vstar[0]-pcat.vstar[1]))
;print,minmax(pcat.sep_kpc) 
;print,minmax(pcat.sep_arcsec)
Mi_cut = -18.0
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 149
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 119
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 113
	pcat.sep_arcsec le pcat.IFUrin and $ ; 105
	pcat.NSA_ELPETRO_ABSMAG[5]+5*alog10(0.7) lt Mi_cut) ; 105
print,minmax(abs(pcat[s].vstar[0]-pcat[s].vstar[1]))
print,minmax(pcat[s].sep_kpc) 
print,minmax(pcat[s].sep_arcsec)
; link PNG files for mergers
manga_regroup,pcat[s].plateifu,'merger_final/',indir='../../5.specpng/'
; print the log file
jk = where(pcat[s].bclass gt -9,cta) ; all
jk = where(pcat[s].bclass lt  0,cti) ; invalid - 8454-6102 [OIII] in sky line
jk = where(pcat[s].bclass eq  0,ct0)
jk = where(pcat[s].bclass eq  1,ct1)
jk = where(pcat[s].bclass eq  2,ct2)
jk = where(pcat[s].bclass eq  3,ct3)
jk = where(pcat[s].bclass eq  4,ct4)
jk = where(pcat[s].bclass eq  5,ct5)
forprint,['Main Sample','Unclass','RG', 'SF', 'Comp', 'LINER', 'Sey2', 'Sey1'],$
       [cta,cti,ct0,ct1,ct2,ct3,ct4,ct5],f='(a10,i6)',text=indir+'.txt'

endfor

end
