
FOR JJ=2,2 DO BEGIN

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

; load assembled SPFIT measurements
str = mrdfits(indir+'0.fits',1)
outfir = indir+'.fits'

; load DRPALL
drpall = mrdfits('../../sample/matched/drpall.fits',1) ; 2718
drpall.plateifu = strtrim(drpall.plateifu,2)

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
help,where(basename ne drpall[sa].plateifu)
str[sa].bclass = 5

; W(Ha) < 3 A: 1205 cases
; Retired Galaxies: bclass = 0
s = where(str.wha lt 3 and str.bclass ne 5) 
str[s].bclass = 0

; W(Ha) > 3 A: 1497 cases
; (1) retain BPT classification
s = where(str.wha ge 3 and str.bclass ne 5)
str[s].bclass = str[s].bpt
; 6 out of 1446 have BPT = -1 (invalid classification because o3hb = -Infinity)
; all of those should be excluded from emission-line galaxy sample
s = where(str.wha ge 3 and str.bclass ne 5 and str.bpt eq -1) 
forprint,';'+drpall[s].plateifu,drpall[s].nsa_z,str[s].o3hb,str[s].n2ha
;7443-3703       0.028599800        Infinity      -2.4503431 wrong redshift z=0.277 (not in main sample)
;8155-12702       0.14868400        Infinity     -0.58422363 this is a star
;8244-12704       0.11513700       -Infinity     -0.32778659 [OIII] in sky line 
;8454-6102        0.11359500       -Infinity      0.18659249 [OIII] in sky line
;8461-6101       0.094323500        Infinity     -0.51468921 wrong redshift z=0.027
;8623-12703      0.052177900       -Infinity     -0.42028224 f/g star on nucleus
; use WHAN class for sources at z ~ 0.114
s = where(drpall.nsa_z ge 0.113 and drpall.nsa_z lt 0.116 and str.whan gt 0) 
str[s].bclass = str[s].whan-1
forprint,';'+drpall[s].plateifu,drpall[s].nsa_z,str[s].bpt,str[s].whan,str[s].bclass
; (2) 2 have Ha AoN < 4 -> set bpt = -1 for those 
s = where(str.wha ge 3 and str.bclass ne 5 and str.bpt ne -1 and str.aon[1] lt 4)
str[s].bclass = -2
forprint,';'+drpall[s].plateifu,drpall[s].nsa_z,str[s].wha,str[s].aon[1],str[s].bpt
;8551-9101       0.040500500       5.5132057      3.68336 1 wrong redshift
;8618-3702        0.12988100       4.7004711      3.05386 3 broad H-alpha (added to t1AGN catalog)
; (3) 1 case where the target fell outside of IFU 
; this is also a known BLAGN, bclass=5 originally
;8239-3701     0.018076700       0.0000000      0.00000          5
offset = sphdist(drpall.OBJRA,drpall.OBJDEC,drpall.IFURA,drpall.IFUDEC,/deg)*3600
s = where(offset gt drpall.IFURAD,ct)
str[s].bclass = -3
forprint,';'+drpall[s].plateifu,drpall[s].nsa_z,str[s].wha,str[s].aon[1],str[s].bclass

; 20 W(Ha)>3 cases have weak [OIII] or Hb lines (AoN2 < 3)
; keep their original BPT classification 
s1 = where(str.wha ge 3 and str.bclass ne 5 and str.bclass ge 0)
s2 = where(str.wha ge 3 and str.bclass ne 5 and str.bclass ge 0 and min(str.aon2[0:3,*],dim=1) lt 3)
plot,str[s1].n2ha,str[s1].o3hb,psym=3,xtit='[N II]/Ha',ytit='[O III]/Hb'
oplot,str[s2].n2ha,str[s2].o3hb,psym=6 
save_screen,'pngs/BPT_N2'+indir+'.png'

; save # of each class in a log file
s = where((drpall.mngtarg1 and 2L^10+2L^11+2L^12) ne 0) ; 2618
forprint,['Main Sample','Unclass','RG', 'SF', 'Comp', 'LINER', 'Sey2', 'Sey1'],$
	[n_elements(where(str[s].bclass gt -9)),$
	n_elements(where(str[s].bclass lt 0)),$
	n_elements(where(str[s].bclass eq 0)),$
	n_elements(where(str[s].bclass eq 1)),$
	n_elements(where(str[s].bclass eq 2)),$
	n_elements(where(str[s].bclass eq 3)),$
	n_elements(where(str[s].bclass eq 4)),$
	n_elements(where(str[s].bclass eq 5))],f='(a10,i6)',text=indir+'.txt'

; save combined structure
mwrfits,struct_combine(str,drpall),indir+'.fits',/create

endfor

end
