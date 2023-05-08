str = mrdfits('nuc1kpc.fits',1)
str.plateifu = strtrim(str.plateifu)

; BPT	CODE - Seyfert=4, LINER=3, Comp=2, SF=1, Ambi=0, Invalid=-1
; WHAN 	CODE - sAGN=4, wAGN=3, SF=1, RG=0, Invalid=-1 (WHa < 3 = RG, WHa > 6 = sAGN)

; make combined classes - bclass
; -1/0/1/2/3/4/5: unreliable/RG/SF/COMP/LINER/SEYFERT/BLAGN
bclass = str.bpt*0-9 ; default class -9 

; Total 2718 galaxies in MPL-5=DR14

; 16 Known BLAGNs
; these all have poor fits because SPFIT cannot properly fit BLAGNs 
; give those bclass = 5
readcol,'../../sample/blagn/t1agn.txt',plates,ifus,f='l,l',comment='#'
basename = strc(plates)+'-'+strc(ifus)
match2,basename,str.plateifu,sa,sb
help,where(basename ne str[sa].plateifu)
str[sa].bpt = 5

; W(Ha) > 3 A: 1446 cases
; 1. 2 have Ha AoN < 4 -> set bpt = -1 for those 
plot,str.wha,str.aon[1],psym=3,/xlog,/ylog,xr=[0.1,1e3],xtit='Ha EW (A)',ytit='Ha AoN'
ver,3,lines=3
hor,4,lines=2
s = where(str.bpt ne -1 and str.bpt ne 5 and str.wha ge 3 and str.aon[1] lt 4)
str[s].bpt = -1
oplot,str[s].wha,str[s].aon[1],psym=6
save_screen,'pngs/AoN_vs_WHa.png'
forprint,';'+str[s].plateifu,str[s].nsa_z,str[s].wha,str[s].aon,str[s].bpt
;8551-9101       0.040500500       5.5132057      3.68336 1 wrong redshift
;8618-3702        0.12988100       4.7004711      3.05386 3 poor cont fit

; 2. retain BPT classification
ind = where(str.wha ge 3)
bclass[ind] = (str.bpt)[ind]

; 5 out of 1446 have BPT = -1 (invalid classification because o3hb = -Infinity)
; all of those should be excluded from emission-line galaxy sample
s = where(str.bpt eq -1 and str.wha ge 3) 
forprint,';'+str[s].plateifu,str[s].nsa_z,str[s].o3hb,str[s].n2ha
;8155-12702       0.14868400        Infinity     -0.58422363 not a galaxy
;8244-12704       0.11513700       -Infinity     -0.32778659 [OIII] in sky line 
;8454-6102        0.11359500       -Infinity      0.18659249 [OIII] in sky line
;8461-6101       0.094323500        Infinity     -0.51468921 wrong redshift z=0.027
;8623-12703      0.052177900       -Infinity     -0.42028224 f/g star on nucleus

; 19 out of 1446 have weak lines (AoN2 < 3)
; those are all fine to keep their original BPT classification 
s = where(str.bpt ne -1 and str.bpt ne 5 and str.wha ge 3 and min(str.aon2[0:3,*],dim=1) lt 3)
s1 = where(str.bpt ne -1 and str.bpt ne 5 and str.wha ge 3)
plot,str[s1].n2ha,str[s1].o3hb,psym=3,xtit='[N II]/Ha',ytit='[O III]/Hb'
oplot,str[s].n2ha,str[s].o3hb,psym=6 
save_screen,'pngs/BPT_N2.png'

end
