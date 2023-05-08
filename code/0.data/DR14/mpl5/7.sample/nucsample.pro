;+
; In 4.spfit_all we extracted and fit all of the nuclear
; spectra from DR14 with a range of aperture sizes, here 
; we assemble basic observed and derived parameters 
; for all galaxy's nuclear spectra  
;-

FOR JJ=1,2 DO BEGIN

CASE jj OF
	1: begin
		indir = 'nuc1arc'
		apradi = 1.0
		apunit = 'arc'
		end
	2: begin
		indir = 'nuc1kpc'
		apradi = 1.3 
		apunit = 'kpc'
		end
	3: begin
		indir = 'nuc2arc'
		apradi = 2.0 
		apunit = 'arc'
		end
ENDCASE

; output: assembled catalog w/ key pars extracted from SPFIT output
outfile = indir+'0.fits'

; load DRPALL
drpall = mrdfits('../../sample/matched/drpall.fits',1) ; 2718
drpall.plateifu = strtrim(drpall.plateifu,2)

; tags to be combined with DRPALL 
tags = ['BPT','WHAN','BCLASS','APrad',$
	'Lambda','AoN','Sigma_Obs','AoN2','AoN3',$ ; line info
	'FO3','dFO3','Lo3',$ ; 1d-17 erg/s/cm2, df/f, erg/s
	'O3Hb','dO3Hb','N2Ha','dN2Ha','WHa','dWHa',$
	'O1Ha','dO1Ha','S2Ha','dS2Ha',$
	'Vgas','Sgas','Vstar','Sstar','Mstar',$
	'dVgas','dSgas','dVstar','dSstar'] 	
vals = ['-9','-9','-9','-9.9',$
	replicate('fltarr(7)',5),$
	replicate('0.d',n_elements(tags)-9)] 
str = mrd_struct(tags,vals,n_elements(drpall))

fmt = '(" Progress:",f6.2,"% done",$,%"\r")'
maxsteps = n_elements(drpall)
for i=0,n_elements(drpall)-1 do begin
	; spfit parameters 
	infile = '../4.spfit_all/'+indir+'/'+drpall[i].plateifu+'.fits'
	if ~file_test(infile) then continue
	pars = mrdfits(infile,1,/silent)

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

	; save emission line fluxes
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
	str[i].lambda = pars.lambda[idx]
	str[i].aon = pars.aon[idx] 
	str[i].sigma_obs = pars.sigma_obs[idx]
	; Amplitude-to-Noise ratio 2
	; amplitude from best-fit Gaussian
	; noise from input error array
	s = where(pars.good eq 1)
	c_kms = 299792.46d ; speed-of-light
	; dV = z*c ~= ln(lambda1/lambda0)*c = log10(lambda1/lambda0)*ln(10)*c
	str[i].aon2= interpol(pars.emis[s]/pars.err[s],pars.log10lam[s],$
		alog10(pars.lambda[idx])+pars.vel[idx]/c_kms/alog(10))
	; Amplitude-to-Noise ratio 3
	; amplitude from residual spectrum (obs - best-fit stellar)
	; noise from input error array
	str[i].aon3= interpol((pars.galaxy-(pars.best-pars.emis))[s]/pars.err[s],pars.log10lam[s],$
		alog10(pars.lambda[idx])+pars.vel[idx]/c_kms/alog(10))
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
	;	CODE - Seyfert=4, LINER=3, Comp=2, SF=1, Ambi=0, Invalid=-1	
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
	str[i].dfo3 = dfo3 ; fractional error
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

	; print progress
	if (i mod round(maxsteps/100d0) eq 0) then print,f=fmt,100.*i/maxsteps
endfor
print, '' ;; don't overwrite the final line

mwrfits,str,outfile,/create

ENDFOR

end
