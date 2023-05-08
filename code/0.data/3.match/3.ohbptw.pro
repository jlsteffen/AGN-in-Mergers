; Joint BPT and WHAN classification using SPFIT results
; 
outdir = './'

; MaNGA SpecObj catalog
spobj = mrdfits('../0.manga_specobj.fits',1)
spobj.plateifu = strtrim(spobj.plateifu)

for j=0,1 do begin
  if j eq 0 then indir='r1arc' else indir='r1kpc'

;;;;;;;;;;;;;;;;
; 12+log(O/H) cat 
;---------------
; SPFIT best-fit parameters
spftpars = mrdfits(indir+'.fits',1)
; SPFIT line names and wavelengths
vbin = mrdfits('../1.spfit/'+indir+'/7443-1901.fits',2)
; calculate log(O/H) using various methods
ohstr = caloh(spftpars, vbin, outfile=outdir+'oh_'+indir+'.fits')

;;;;;;;;;;;;;;;;
; BPT/WHAN cat
;---------------
; get SSP parameters - Age & MH
ssp = mrdfits('/Volumes/scr/manga/spfit/MPL-11/miuscat-thin/7443-1901.fits',1)
; create output BPT-WHAN structure
tags = ['BPT','WHAN','BPTWHAN',$ ; classification flags 
	'AoN','z',$ ; em-line A/N ratio, abs-line redshift
	'logM','dlogM','Age','MH',$  ; stellar mass in log(Msun)
	'FO3','dFO3','logLo3',$	; [OIII]5007 line flux: 1d-17 erg/s/cm2, df/f, erg/s
	'WHa','dWHa',$     	; de-redshifted H-alpha EW in Angstrom
	'O3Hb' ,'N2Ha' ,'O1Ha' ,'S2Ha' ,$  ; log(line ratios)
	'dO3Hb','dN2Ha','dO1Ha','dS2Ha',$  ; errors
	'Vgas' ,'Sgas' ,'Vstar' ,'Sstar',$ ; gas & stellar kinematics
	'dVgas','dSgas','dVstar','dSstar'] ; errors
vals = ['-9','-9','-9','fltarr(7)',$
	replicate('-9.9',n_elements(tags)-4)] 
nobj = n_elements(spobj)
str = mrd_struct(tags,vals,nobj)
; get emission line indices 
linename = strtrim(vbin.linename,2)
linewave = vbin.linewave
io2 = where(linename eq 'OII3730')  ; 3727 tied to 3730 @ 1:1
io3 = where(linename eq 'OIII5008') ; 4960 tied to 5008 @ 0.35:1
ihb = where(linename eq 'Hb4863')
in2 = where(linename eq 'NII6585')  ; 6549 tied to 6585 @ 0.34:1
iha = where(linename eq 'Ha6565')
io1 = where(linename eq 'OI6302')
is2a= where(linename eq 'SII6718') 
is2b= where(linename eq 'SII6733')

; loop through all objects
fmt = '(" Progress:",f6.2,"% done",$,%"\r")'
for i=0,nobj-1 do begin
	; skip non-good fits
	if spobj[i].sptype ne 1 then continue 
	pars = spftpars[i]

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

	; Amplitude to noise ratio
	str[i].aon = pars.aon[[in2,iha,io3,ihb,io1,is2a,is2b]]

	; WHAN classification
	;	CODE - integer, sAGN=4, wAGN=3, SF=1, RG=0, Invalid=-1
	; H-alpha EW, rest-frame
	str[i].wha = pars.EW[iha]
	; this error doesn't include continuum uncertainty
	str[i].dwha = str[i].wha * emlerr(pars.aon[iha],pars.sigma_obs[iha])
	whanclass = whan(str[i].wha,str[i].n2ha,code=whancode)
	str[i].whan = whancode

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
	; redshift - dz = dv/c
	str[i].z = pars.redshift + pars.kinstar[0,0]/c_kms

	; because fit_manga use mean to combine spectra in each bin
	; we need to find out the number of spaxels inside the aperture
	; to recover the total line flux and stellar mass
   if j eq 0 then begin
      rad = 1.0 ; arcsec
	endif else begin
      objz = pars.redshift
      apradi = 1.3 ; kpc
		pltscale = dangular(objz,/kpc)/(180.*3600./!pi) ; kpc/arcsec
		rad = apradi/pltscale ; kpc->arcsec
   endelse
	area = !pi*rad^2 ; arcsec^2, total aperture size
	
   ; [O III] fluxes - 1d-17 erg/s/cm2/arcsec^2 -> erg/s/cm2
	str[i].fo3 = fo3*area
	; fractional error - df/f
	str[i].dfo3 = dfo3
	; 1d-17 erg/s/cm2 * Mpc^2 -> erg/s
	; 1d-17 * 4 * !PI * (3.086d24)^2 = 1.1967453e+33
	str[i].logLo3 = alog10(str[i].fo3*(lumdist(pars.redshift,/silent))^2*1.197d33)
	; stellar mass (logMsun), age (Gyr), and [M/H]
	str[i].logM = alog10(total(pars.m_star[*,0],1)*area)
	str[i].dlogM = sqrt(total(pars.m_star[*,1]^2))/total(pars.m_star[*,0],1)/alog(10)
	str[i].age = total(ssp.age*pars.m_star[*,0])/total(pars.m_star[*,0])
	str[i].MH = alog10(total(10.^ssp.mh*pars.m_star[*,0])/total(pars.m_star[*,0]))
	
	; print progress
	if (i mod round(nobj/100d0) eq 0) then print,f=fmt,100.*i/nobj
endfor
print, '' ;; don't overwrite the final line

;;;;;;;;;;;;
; Joint BPT+WHAN classification codes:
; 0/1/2/3/4/5: RG/SF/COMP/LINER/SEYFERT/BLAGN
bptclass = bpt_k06_n2(str.o3hb,str.n2ha,EW=str.WHa,code=bptwhan)
str.bptwhan = bptwhan
; BLAGN identified w/ MaNGA spectra
str[where(spobj.sptype eq 4)].bptwhan = 5

mwrfits,str,outdir+'bptw_'+indir+'.fits',/create

endfor

end
