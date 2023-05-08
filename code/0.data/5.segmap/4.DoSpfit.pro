;+
; Run fit_manga on all plateifus that contain kinematic pairs or
; multiples, using sextractor segmentation maps to define the 
; extraction regions
;
; % Time elapsed: 9569.5381 seconds.
;-

; pair sample
ncomps = mrdfits('../4.pairs/ncomp.fits',1)
pair = ncomps[where(ncomps.ngood ge 2,npair)]
plateifu = strtrim(pair.plateifu)

; string -> long
plates   = long(((strsplit(plateifu,'-',/ex)).ToArray())[*,0])
ifus     = long(((strsplit(plateifu,'-',/ex)).ToArray())[*,1])

;; test run
;fit_manga,plates[0],ifus[0],tag='MPL-11',mdegree=6,ds9reg='rseg/',outdir='rseg_spfit/',$
;	/ignore_drp3qual,/overwrite,/verb,/saveplot,prebin=1
;stop

; run fit_manga in parallel mode
cmd = 'fit_manga'
inpdir = 'rseg/'
outdir = 'rseg_spfit/'
extra = ',tag=''MPL-11'',mdegree=6,'+$
	'ds9reg='''+inpdir+''','+$
	'outdir='''+outdir+''','+$
	'/ignore_drp3qual,/overwrite,prebin=1,/quiet'
outputs=inpdir+strc(plates)+'-'+strc(ifus)+'.fits' 
manga_parallel,plates,ifus,cmd,extra,ncpu=6,outputs=outputs,/skip

end
