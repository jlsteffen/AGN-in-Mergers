; use Hopkins10's merger_rate_calculator.pro to explore the merger
; fraction vs. stellar mass

goto,plot2

outfile = 'fpair_vs_mass.png'
np = 20
mass = range(9.0,12.0,np) 
fpair1 = fltarr(np)
fpair2 = fltarr(np)

for i=0,np-1 do begin
	redshift = 0.05
	m_min = mass[i]-0.1 ; log(Msun)
	m_max = mass[i]+0.1
	mu_min = 0.1 ; mass ratio
	mu_max = 1.0
	fgas_min = 0.0 ; gas fraction
	fgas_max = 1.0
	h = 0.7 ; Hubble constant, h=H0/100
	sep_max = 40*h  ; h^-1 kpc, projected maximum pair separation
	fpair1[i] = merger_rate_calculator(redshift,m_min,m_max,mu_min,mu_max,fgas_min,fgas_max,$
		RETURN_MERGER_FRACTION_PAIRS=sep_max,/USE_STELLAR_MASSES,/quiet)
	fpair2[i] = merger_rate_calculator(redshift,m_min,m_max,mu_min,mu_max,fgas_min,fgas_max,$
		RETURN_MERGER_FRACTION_PAIRS=sep_max,/USE_STELLAR_MASSES,/quiet,/USE_ALTERNATIVE_MASSFUNCTION)
endfor
plot,mass,fpair1,/ylog,xtit='Stellar Mass',ytit='Pair Fraction'
oplot,mass,fpair2,lines=2
ver,11.4,lines=2
ver,10.0,lines=2
;oplot,mass,0.006*(exp(((mass-9)/1.6)^2.5)),lines=2,color=cgcolor('green')
;oplot,mass,0.01*exp((mass-10-0.4)*2.8),lines=2,color=cgcolor('red')
; double powerlaw
logx = mass
logy1 = -2.04 + 0.23*(logx-10)
oplot,logx,10d^logy1,lines=3,color=cgcolor('red')
logy2 = -1.05 + 1.7*(logx-11.4)
oplot,logx,10d^logy2,lines=3,color=cgcolor('red')
oplot,logx,10d^logy1+10d^logy2,lines=2,color=cgcolor('red')

; evolving powerlaw
logy = -2.04 + (logx-10)*(0.23+0.4*((logx-10)>0))
oplot,logx,10d^logy,lines=2,color=cgcolor('green')

save_screen,outfile
;stop

plot2:
;;;;;;;;;;;;;;;;;;;;;;;;;;;
outfile = 'fpair_vs_z.png'
np = 10
z = range(0.0,0.5,np) 
fpair1 = fltarr(np)
fpair2 = fltarr(np)
for i=0,np-1 do begin
	redshift = z[i]
	m_min =  9.5 ; log(Msun)
	m_max = 11.5
	mu_min = 0.1 ; mass ratio
	mu_max = 1.0
	fgas_min = 0.0 ; gas fraction
	fgas_max = 1.0
	h = 0.7 ; Hubble constant, h=H0/100
	sep_max = 40*h  ; h^-1 kpc, projected maximum pair separation
	fpair1[i] = merger_rate_calculator(redshift,m_min,m_max,mu_min,mu_max,fgas_min,fgas_max,$
		RETURN_MERGER_FRACTION_PAIRS=sep_max,/USE_STELLAR_MASSES,/quiet)
	fpair2[i] = merger_rate_calculator(redshift,m_min,m_max,mu_min,mu_max,fgas_min,fgas_max,$
		RETURN_MERGER_FRACTION_PAIRS=sep_max,/USE_STELLAR_MASSES,/quiet,/USE_ALTERNATIVE_MASSFUNCTION)
endfor
; compute slope of evolution
; assume f_pair ~ (1+z)^alpha, solve alpha
z1 = 0.01
z2 = 0.15
alpha = alog(interpol(fpair1,z,z1)/interpol(fpair1,z,z2))/alog((1+z1)/(1+z2))
plot,z,fpair1,xtit='Redshift',ytit='Pair Fraction',/xs,/ys,$
	tit=textoidl('\alpha = '+strc(alpha)) ;yr=[0.0,0.1]
oplot,z,fpair2,lines=2
ver,[z1,z2],lines=1
oplot,z,fpair1[0]*(1+z)^alpha,color=cgcolor('red'),lines=2

save_screen,outfile


end
