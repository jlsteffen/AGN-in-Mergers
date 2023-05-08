;+
; define functions to be used in 
; phi1agn.pro & phi2agn.pro
;-
pro mkexcess,xgri,xerr,phi_obs,phi_exp, _REF_EXTRA=extra
; plot excess panel
	; Lower panel
	pos = !p.position
	pos[3] = 0.4
	plot,[0,1],[0,1],/nodata,pos=pos,_extra=extra 
	hor,0.0,lines=2,color=cgcolor('royal blue')

	; total Phi for expected AGNs in close pairs
	y1 = alog10(phi_exp[*,1])
	; d log(y) = 1/ln(10) dy/y
	;y1lerr = (phi_exp[*,1]-phi_exp[*,0])/phi_exp[*,1]/alog(10)
	;y1uerr = (phi_exp[*,2]-phi_exp[*,1])/phi_exp[*,1]/alog(10)
	; directly convert range to log scale
	y1lerr = alog10(phi_exp[*,1])-alog10(phi_exp[*,0])
	y1uerr = alog10(phi_exp[*,2])-alog10(phi_exp[*,1])
	; add 0.04 dex uncertainty in f_agn model 
	moderr = 0.04
	y1lerr = sqrt(y1lerr^2+moderr^2)
	y1uerr = sqrt(y1uerr^2+moderr^2)
	; total Phi for observed AGNs in close pairs
	y2 = alog10(phi_obs[*,1])
	y2lerr = alog10(phi_obs[*,1])-alog10(phi_obs[*,0])
	y2uerr = alog10(phi_obs[*,2])-alog10(phi_obs[*,1])
	; compute AGN excess
	y = y2-y1
	ylerr = sqrt(y2lerr^2+y1lerr^2)
	yuerr = sqrt(y2uerr^2+y1uerr^2)

	; show excess
	oploterror,xgri,y,xerr,ylerr,/lobar,psym=3,color=cgcolor('red')
	oploterror,xgri,y,xerr,yuerr,/hibar,psym=3,color=cgcolor('red')
	plots,xgri,y,psym=cgsymcat(15),color=cgcolor('white'),noclip=0
	plots,xgri,y,psym=cgsymcat(6),color=cgcolor('red'),noclip=0

	; upper limits
	s = where(phi_obs[*,1] eq 0,ct)
	if ct gt 0 then begin
		y = alog10(phi_obs[s,0])-alog10(phi_exp[s,1])
		plotsym,1,5.0,thick=4
		s1 = where(y lt 1.2,ct2)
		if ct2 gt 0 then $
		oploterror,xgri[s],y,xerr[s],s*0,psym=8,color=cgcolor('red')
		s2 = where(y gt 1.2,ct2)
		if ct2 gt 0 then $
		oplot,xgri[s[s2]],s2*0+1.2,psym=8,color=cgcolor('red')
	endif

end

pro mkplot,xx,yy,idx,xgri,xerr,phi_obs,phi_exp, _REF_EXTRA=extra
; plot main panel
	; phi are given in lower limit, mean, and upper limit
	; Upper panel
	pos = !p.position
	pos[1] = 0.4 
	plot,[0,1],[0,1],/nodata,pos=pos,_extra=extra 
	syms = 0.8
	; show 1/Vmax weights
	plots,xx,yy,psym=cgsymcat(9),syms=syms,color=cgcolor('dark gray'),noclip=0
	plots,xx[idx],yy[idx],psym=cgsymcat(16),syms=syms,noclip=0,color=cgcolor('dark gray')
	; total Phi for expected AGNs in close pairs
	y = alog10(phi_exp[*,1])
	ylerr = alog10(phi_exp[*,1])-alog10(phi_exp[*,0])
	yuerr = alog10(phi_exp[*,2])-alog10(phi_exp[*,1])
	; add 0.04 dex uncertainty in f_agn model 
	moderr = 0.04
	ylerr = sqrt(ylerr^2+moderr^2)
	yuerr = sqrt(yuerr^2+moderr^2)
	; draw boxes
	for i=0,n_elements(xgri)-1 do $
	tvbox,[xerr[i]*2,yuerr[i]+ylerr[i]],xgri[i],(y[i]-ylerr[i]+y[i]+yuerr[i])/2,$
		/data,color=cgcolor('royal blue'),noclip=0,lines=0,$
		/fill,/line_fill,orient=45,spacing=0.2

	; total Phi for observed AGNs in close pairs
	y = alog10(phi_obs[*,1])
	ylerr = alog10(phi_obs[*,1])-alog10(phi_obs[*,0])
	yuerr = alog10(phi_obs[*,2])-alog10(phi_obs[*,1])
	oploterror,xgri,y,xerr,ylerr,/lobar,psym=3,color=cgcolor('red')
	oploterror,xgri,y,xerr,yuerr,/hibar,psym=3,color=cgcolor('red')
	plots,xgri,y,psym=cgsymcat(15),color=cgcolor('white'),noclip=0
	plots,xgri,y,psym=cgsymcat(6),color=cgcolor('red'),noclip=0	
	; upper limits
	s = where(phi_obs[*,1] eq 0,ct)
	if ct gt 0 then begin
		y = alog10(phi_obs[s,0])
		plotsym,1,5.0,thick=4
		oploterror,xgri[s],y,xerr[s],s*0,psym=8,color=cgcolor('red')
	endif
end

function calc_phi_obs, w, ind_agn, ind_pair, nboot
; calculate observed volume densities of AGNs
; bootstrap_mean:
;    Returns a 3-element vector containing the lower limit, mean, and
;    upper limit.
	phi_obs = fltarr(3) ; lower limit, mean, and upper limit
	jk = where(ind_agn ge 0,nagn)
	jk = where(ind_pair ge 0,npair)
	if nagn gt 0 then begin
		; use bootstrap when there are more than 4 AGNs
		if nagn ge 3 then begin
			phi_obs = bootstrap_mean(w[ind_agn],nboot=nboot)*nagn
		endif else begin
			phi_obs[1] = total(w[ind_agn])
			; convert binormial confidence intervals
			; implicit assumption is that all objects have the same
			; 1/Vmax weights
			err = binormial_ci(nagn,npair)
			; [k/n,lower_err,upper_err,eratio(k,sqrt(k),n,sqrt(n))]
			meanw = mean(w[ind_pair])
			phi_obs[0] = (err[0]-err[1]) * npair * meanw
			phi_obs[2] = (err[0]+err[2]) * npair * meanw
		endelse
	endif else begin
		; 3-sigma upper limit
		nmax = npair * (binormial_ci(0,npair,cl=0.997))[2]
		phi_obs[0] = mean(w[ind_pair]) * nmax
	endelse
	return,phi_obs
end

pro calc_phi,x,w,xgri,xbin,xi,flg_bagn,flg_pagn,flg_sagn,fagn_p,fagn_s,$
	phi_1agn_exp=phi_1agn_exp,phi_1agn_obs=phi_1agn_obs,$
	phi_bagn_exp=phi_bagn_exp,phi_bagn_obs=phi_bagn_obs
; calculate expected and observed volume densities of AGNs
; bootstrap_mean:
;    Returns a 3-element vector containing the lower limit, mean, and
;    upper limit.
	xnp = n_elements(xgri)
	phi_bagn_exp = fltarr(xnp,3)
	phi_bagn_obs = fltarr(xnp,3)
	phi_1agn_exp = fltarr(xnp,3) ; either primary or secondary
	phi_1agn_obs = fltarr(xnp,3) ; either primary or secondary
	nboot = 1000
	for i=0,n_elements(xgri)-1 do begin
		; expected volume densities
		insidebin = x ge xgri[i]-xbin[i]/2 and x lt xgri[i]+xbin[i]/2
		ind_pair = where(insidebin, npair)
		if npair gt 0 then begin
			ind = ind_pair
			ct = npair
			wbagn = w[ind]*fagn_p[ind]*fagn_s[ind] + $
				w[ind]*fagn_p[ind]*xi[i]
			phi_bagn_exp[i,*] = bootstrap_mean(wbagn,nboot=nboot)*ct
			w1agn = w[ind]*fagn_p[ind] + w[ind]*fagn_s[ind] - $
				w[ind]*fagn_p[ind]*fagn_s[ind]
			phi_1agn_exp[i,*] = bootstrap_mean(w1agn,nboot=nboot)*ct
		endif else continue
		
		; observed volume densities
		; binary AGNs
		ind = where(flg_bagn and insidebin,nbagn)
		phi_bagn_obs[i,*] = calc_phi_obs(w,ind,ind_pair,nboot)
		; either primary or secondary is AGN
		ind1 = where(flg_pagn and insidebin, ct1)
		ind2 = where(flg_sagn and insidebin, ct2)
		if ct1 gt 0 and ct2 gt 0 then ind = [ind1,ind2]
		if ct1 gt 0 and ct2 eq 0 then ind = ind1
		if ct1 eq 0 and ct2 gt 0 then ind = ind2
		phi_1agn_obs[i,*] = calc_phi_obs(w,ind,ind_pair,nboot)

		; print results
		print,xgri[i],nbagn,ct1+ct2,npair		
	endfor

	return
end

