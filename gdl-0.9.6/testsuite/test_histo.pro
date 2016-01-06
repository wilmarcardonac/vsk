;
; AC 01-Jun-2007
; SA 30-Aug-2009 (TEST_HISTO_BASIC)
; AC 06-Dec-2011 (adding TEST_HISTO_NAN)
; AC 20-Feb-2013 (adding TEST_HISTO_UNITY_BIN)
;
pro TEST_HISTO_RANDOMU, nbp=nbp, nan=nan
;
if (N_ELEMENTS(nbp) EQ 0) then nbp=1e2
a=randomu(seed,nbp)
;
if KEYWORD_SET(nan) then begin
    j=[round(nbp/3.), round(nbp/2.), round(nbp*2/3.)]
    a[j]=!values.f_nan
    print, j
endif
plot, a, psym=10
end
;
; based on a IDL example
;
pro TEST_HISTO_GAUSS, test=test
;
; Two-hundred values ranging from -5 to 4.95:  
X = FINDGEN(200) / 20. - 5.  
; Theoretical normal distribution, scale so integral is one:  
Y = 1/SQRT(2.*!PI) * EXP(-X^2/2) * (10./200)  
; Approximate normal distribution with RANDOMN,  
; then form the histogram.  
H = HISTOGRAM(RANDOMN(SEED, 2000), $  
  BINSIZE = 0.4, MIN = -5., MAX = 5.)/2000.  
;
h_x=FINDGEN(26) * 0.4 - 4.8
; Plot the approximation using "histogram mode."  
PLOT,h_x, H, PSYM = 10  
; Overplot the actual distribution:  
OPLOT, X, Y * 8.  
;
if KEYWORD_SET(test) then stop
;
end

; SA: intended for checking basic histogram functionality
pro TEST_HISTO_BASIC

  ; for any input if MAX/MIN kw. value is the max/min element of input
  ; it shoud be counted in the last/first bins
  message, 'TEST 01', /continue
  for e = -1023, 1023 do begin
    input = [-2d^e, 2d^e]
    if ( $
      ~array_equal(histogram(input, max=input[1], min=input[0], nbins=2, $
        reverse=ri), [1,1]) $
    ) then begin
      print, histogram(input, max=input[1], min=input[0], nbins=2)
      message, 'FAILED: ' + string(e)
    endif
  endfor
  ignored = histogram([0.], min=0, max=0, reverse=ri) 

  ; test if binsize=(max-min)/(nbins-1) when nbins is set and binsize is not set
  message, 'TEST 02', /continue
  for type = 1, 15 do if type lt 6 or type gt 11 then begin ; data-type loop 
    data = make_array(100, type=type, index = type ne 7)
    for nbins = 2, 1000 do begin
      a = histogram(data, nbins=nbins, loc=l)
      if total(l[0:1] * [-1, 1]) ne (max(data)-min(data))/(nbins-1) then message, 'FAILED'
      a = histogram(data, nbins=nbins, max=max(data), loc=l)
      if total(l[0:1] * [-1, 1]) ne (max(data)-min(data))/(nbins-1) then message, 'FAILED'
      a = histogram(data, nbins=nbins, min=min(data), loc=l)
      if total(l[0:1] * [-1, 1]) ne (max(data)-min(data))/(nbins-1) then message, 'FAILED'
      a = histogram(data, nbins=nbins, min=min(data), max=max(data), loc=l)
      if total(l[0:1] * [-1, 1]) ne (max(data)-min(data))/(nbins-1) then message, 'FAILED'
    endfor
  endif

  ; TODO: test other possible keyword/input combinations...

end
;
; ----------------------------------
;
; array "b" did not contain +/- Inf, it is OK
; array "c" did contain +/- Inf, it is not OK without /nan
;
pro TEST_HISTO_NAN, all=all
;
a=FINDGEN(10)
;
b=a
b[5]=!values.f_nan
;
c=b
c[7]=!values.f_infinity
;
print, HISTOGRAM(a)
print, HISTOGRAM(b)
print, HISTOGRAM(b,/nan)
if KEYWORD_SET(all) then begin
   print, HISTOGRAM(c)
   print, HISTOGRAM(c,/nan)
endif
;
print, HISTOGRAM(a, bin=2)
print, HISTOGRAM(b, bin=2)
print, HISTOGRAM(b, bin=2,/nan)
if KEYWORD_SET(all) then begin
   print, HISTOGRAM(c, bin=2)
   print, HISTOGRAM(c, bin=2,/nan)
endif
;
print, HISTOGRAM(a, nbin=4)
print, HISTOGRAM(b, nbin=4)
print, HISTOGRAM(b, nbin=4,/nan)
if KEYWORD_SET(all) then begin
   print, HISTOGRAM(c, nbin=4)
   print, HISTOGRAM(c, nbin=4,/nan)
endif
;
end
;
; see bug report 3602623
; http://sourceforge.net/tracker/?func=detail&aid=3602623&group_id=97659&atid=618683
;
; TBC: the effect seems to be different on 32b and 64b machines ...
;
pro TEST_HISTO_UNITY_BIN, nbp, display=display, test=test, help=help
;
if KEYWORD_SET(help) then begin
   print, 'pro TEST_HISTO_UNITY_BIN, nbp, display=display, test=test, help=help'
   return
endif 
;
if (N_PARAMS()) EQ 0 then nbp=13000
;
; if 13000 points, we create a shawtooth with 10 points in each unity bin ...
ramp=LINDGEN(nbp) mod 1300
;
h1 = HISTOGRAM(ramp, bin=1)
h2 = HISTOGRAM(ramp)
;
diff=TOTAL(ABS(h2 - h1))
;
if (diff GT 0.0) then begin
   MESSAGE, 'error !', /continue
endif
;
if KEYWORD_SET(display) then begin
   plot, h1, yrange=[-1, 21], /ystyle
   oplot, h2, psym=2
endif 
;
if KEYWORD_SET(test) then STOP
;
end
;
; ------------------------------------------------------------------
;
pro TEST_HISTO
;
TEST_HISTO_UNITY_BIN
;
end

