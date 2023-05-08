;;;;;;;;;;;;
; figure out the resampling rate of bootstrap_median.pro 
;;;;;;;;;;;;
n = 1000 ; # of elements in value array
nboot = 500
mix = long(randomu(seed,nboot,n)*n) ; (nboot * n) uniform [0-n] distribution array
; count number of uniq indices
nuniq = fltarr(nboot)
for i=0,nboot-1 do nuniq[i] = n_elements(rem_dup(mix[i,*]))
; 63% unique indices, meaning replacing 1/e = 36.8% of the sample
print,minmax(nuniq)/n,mean(nuniq)/n,median(nuniq)/n

;;;;;;;;;;;
; How can we test if bootstrap makes sense?
;;;;;;;;;;;
; Gaussian dist. centered on zero w/ stddev of unity
values = randomn(seed,n,/normal) 
vals = values[0:50]
xr = [0,1200]
yr=[0,2]
plot,[0,1],[0,1],/nodata,xr=xr,yr=yr
for nboot=100,1000,100 do begin
	tmp = bootstrap_median(vals,nboot=nboot)
	plots,nboot,(tmp[2]-tmp[0])/2,psym=6
endfor

end
