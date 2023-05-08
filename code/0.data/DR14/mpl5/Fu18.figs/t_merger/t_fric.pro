; trying to figure out t_fric vs. separation (D) 
; by comparing with the D vs. t plot from a 
; merger simulation
; 
; Here t_fric(r) is the time a merger spends at 
; separations less than r, and we assume t_fric ~ r^alpha
; 
; t(r < r_i) = t_m - t(r) = t_m * (r_i/r0)^a
; so that t(r < r0) = t_m, t_m is the total merging time
; => r_i/r0 = (1 - t/t_m)^(1/a)
; 
; The comparison shows alpha = 1 gives the best-fit to the simulations
; This is in agreement with the galaxy cannibalism inspiral time formula 
; derived in Binney & Tremaine Eq. 8.17.

read_png,'Volonteri15Fig1.png',im
wind,xsize=(size(im))[2],ysize=(size(im))[3]
tvimage,im[0:2,*,*]

pos = [0.137,0.7473,0.98,0.9643]
plot,[0,1],/nodata,xr=[0,1.35],yr=[5e-3,600],/xs,/ys,$
	color=cgcolor('red'),pos=pos,/noerase,/yl

tm = 1.03 ; Gyr
r0 = 80 ; kpc
t = range(0,tm,1e3)
alpha = 1 ; t(r < r_i) ~ (r_i/r0)^alpha
r = r0 * (1.0 - t/tm)^(1./alpha)
oplot,t,r,lines=2,color=cgcolor('red'),thick=2

alpha = 2 ; t(r < r_i) ~ (r_i/r0)^alpha
r = r0 * (1.0 - t/tm)^(1./alpha)
oplot,t,r,lines=2,color=cgcolor('blue'),thick=2

ver,tm,lines=2,color=cgcolor('red')

save_screen,'t_fric.png'

end



