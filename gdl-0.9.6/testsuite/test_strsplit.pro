;
; some tests for STRPLIT
;
; Lea Noreskal, June 2010, under GNU GPL V2 or later
;
; Alain Coulais, October 2010: 
; -- few extra cases. Also adding tests on count= and lenght= keywords
;
; Alain Coulais, August 2011: 
; -- inaccurate tests found, corrected (ARRAY_EQUAL(a,b,/NO_TYPECONC))
; -- change from /quiet to /verbose
; -- no exit on error if /debug or /test (suppose to be interactive ;-)
; -- adding test for bug 3286746 (STR_SEP)
; -- adding basic test for TexToIDL
;
; A bug 554 about "length" was corrected in July 2013.
;
pro IPRINT, texte, indice
print, texte, indice
indice=indice+1
end
;
pro TEST_STRSPLIT, no_exit=no_exit, verbose=verbose, test=test, debug=debug, help=help
;
if KEYWORD_SET(help) then begin
   print, 'pro TEST_STRSPLIT, no_exit=no_exit, verbose=verbose, test=test, debug=debug, help=help'
   return
endif
;
nb_pbs=0
indice=0
;
str='$$$$$$$$$$' ; LONG 0
str1='$chops &up str*ings.' ; [1]
str2='ch$ops &up str*ings.' ; [0,3]
str3='$ch$ops &up str*ings.' ; [1,4]
str4='../foo.txt' ; [2,7]
str5=' t e s t '
str6='qwerty'
str7='aa;aaaaaaaaa;aaaaaaaaaa'  ; cf bug 554
;
; When Search Pattern is not in Input String (default Search pattern
; is white space ' ')
;
tab=STRSPLIT(str6)
if (tab ne 0) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str basic 0', /continue
   nb_pbs=nb_pbs+1
endif
tab=STRSPLIT(str6,/extract)
if (tab ne str6) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str basic 0 bis', /continue
   nb_pbs=nb_pbs+1
endif
;
if KEYWORD_SET(debug) then IPRINT, 'ici', indice
;
; When Search Pattern is a stupid void input ''
;
tab=STRSPLIT(str5,'')
if ((SIZE(tab,/type) NE 3) and (tab ne 0)) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str basic 1', /continue
   nb_pbs=nb_pbs+1
endif
tab=STRSPLIT(str5,'',/extract)
if (tab ne '') then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str basic 1 bis', /continue
   nb_pbs=nb_pbs+1
endif
;
if KEYWORD_SET(debug) then IPRINT, 'ici', indice
;
; When Search Pattern is not in Input String
;
tab=STRSPLIT(str5,'a')
if ((SIZE(tab,/type) NE 3) and (tab ne 0)) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str basic 2', /continue
   nb_pbs=nb_pbs+1
endif
tab=STRSPLIT(str5,'a',/extract)
if (tab ne str5) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str basic 2 bis', /continue
   nb_pbs=nb_pbs+1
endif
;
;
if KEYWORD_SET(debug) then IPRINT, 'ici', indice
;
tab=STRSPLIT(str,'$')
if ((SIZE(tab,/type) NE 3) and (tab ne 0)) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str', /continue
   nb_pbs=nb_pbs+1
endif

ext=STRSPLIT(str,'$',/ext)
;print , str , ' Fext> ' , strjoin(ext,'-') , N_ELEMENTS(ext)

if(N_ELEMENTS(ext) ne 1) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str extract', /continue
   nb_pbs=nb_pbs+1
endif

tab=STRSPLIT(str,'$',/preserve)
res=LINDGEN(11)
if (ARRAY_EQUAL(tab,res, /NO_TYPECONV) eq 0) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str preserve_null', /continue
   nb_pbs=nb_pbs+1
endif

ext=STRSPLIT(str,'$',/ext,/preserve)
;print , str , ' Fext> ' , strjoin(ext,'-') , N_ELEMENTS(ext)
if (N_ELEMENTS(ext) ne 11) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str extract', /continue
   nb_pbs=nb_pbs+1
endif

ext=STRSPLIT(str,'$',/preserve,length=length)
if ARRAY_EQUAL(length,replicate(0,STRLEN(str)+1)) eq 0 then begin 
   if KEYWORD_SET(verbose) then MESSAGE, 'error str len+preserve', /continue
   nb_pbs=nb_pbs+1
endif

ext=STRSPLIT(str,'$',length=length)
if (length ne 0) then begin 
   if KEYWORD_SET(verbose) then MESSAGE, 'error str len', /continue
   nb_pbs=nb_pbs+1
endif 
;
if KEYWORD_SET(debug) then IPRINT, 'ici', indice
;
tab=STRSPLIT(str1,'$')
res=[1L]
if (ARRAY_EQUAL(tab,res, /NO_TYPECONV) eq 0) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str1', /continue
   nb_pbs=nb_pbs+1
endif

ext=STRSPLIT(str1,'$',/ext)
;print , str1 , ' Fext> ' , strjoin(ext,'-') , N_ELEMENTS(ext)

res=['chops &up str*ings.']
if( (N_ELEMENTS(ext) ne 1) AND (ARRAY_EQUAL(tab,res, /NO_TYPECONV) eq 0 )) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str extract', /continue
   nb_pbs=nb_pbs+1
endif

ext=STRSPLIT(str1,'$',/ext,/preserve)
;print , str1 , ' Fext> ' , strjoin(ext,'-') , N_ELEMENTS(ext)
if((N_ELEMENTS(ext) ne 2)) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str1 extract preserves', /continue
   nb_pbs=nb_pbs+1
endif
;
if KEYWORD_SET(debug) then IPRINT, 'ici', indice

tab=STRSPLIT(str2,'$')
res=LONG([0,3])
if (ARRAY_EQUAL(tab,res, /NO_TYPECONV) eq 0) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str2', /continue
   nb_pbs=nb_pbs+1
endif

ext=STRSPLIT(str2,'$',/ext)
;print , str2 , ' Fext> ' , strjoin(ext,'-') , N_ELEMENTS(ext)

res=['ch' , 'ops &up str*ings.' ]
if ((N_ELEMENTS(ext) ne 2) AND (ARRAY_EQUAL(tab,res, /NO_TYPECONV) eq 0 )) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str extract', /continue
   nb_pbs=nb_pbs+1
endif

ext=STRSPLIT(str2,'$',/ext,/preserve)
;print , str2 , ' Fext> ' , strjoin(ext,'-') , N_ELEMENTS(ext)
if ((N_ELEMENTS(ext) ne 2)) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str2 extract preserve', /continue
   nb_pbs=nb_pbs+1
endif
;
;
tab=STRSPLIT(str2,'$',/preserve)
res=LONG([0,3])
if (ARRAY_EQUAL(tab,res,/NO_TYPECONV) eq 0) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str2 preserve_null', /continue
   nb_pbs=nb_pbs+1
endif

tab=STRSPLIT(str3,'$')
res=LONG([1,4])
if (ARRAY_EQUAL(tab,res, /NO_TYPECONV) eq 0) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str3', /continue
   nb_pbs=nb_pbs+1
endif

ext=STRSPLIT(str3,'$',/ext)
;print , str3 , ' Fext> ' , strjoin(ext,'-') , N_ELEMENTS(ext)

res=['ch' , 'ops &up str*ings.' ]
if ((N_ELEMENTS(ext) ne 2) AND (ARRAY_EQUAL(tab,res, /NO_TYPECONV) eq 0 )) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str extract', /continue
   nb_pbs=nb_pbs+1
endif
;

ext=STRSPLIT(str3,'$',/ext,/preserve)
;print , str3 , ' Fext> ' , strjoin(ext,'-') , N_ELEMENTS(ext)
if ((N_ELEMENTS(ext) ne 3)) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str3 extract preserve', /continue
   nb_pbs=nb_pbs+1
endif
;
if KEYWORD_SET(debug) then IPRINT, 'ici', indice
;
tab=STRSPLIT(str4,'.')
res=LONG([2,7])
if (ARRAY_EQUAL(tab,res, /NO_TYPECONV) eq 0) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str4', /continue
   nb_pbs=nb_pbs+1
endif

tab=STRSPLIT(str4,'.',/preserve)
res=LONG([0,1,2,7])
if (ARRAY_EQUAL(tab,res, /NO_TYPECONV) eq 0) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str4 preserve_null', /continue
   nb_pbs=nb_pbs+1
endif

ext=STRSPLIT(str4,'.',/ext,/preserve)
;print , str4 , ' Fext> ' , strjoin(ext,'-') , N_ELEMENTS(ext)
if ((N_ELEMENTS(ext) ne 4)) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str4 extract , preserve', /continue
   nb_pbs=nb_pbs+1
endif
tab=STRSPLIT(str4,'.',length=length)
if ARRAY_EQUAL(length,[4,3]) eq 0 then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str4 length', /continue
   nb_pbs=nb_pbs+1
endif
tab=STRSPLIT(str4,'.',length=length,/preserve)
if ARRAY_EQUAL(length,[0,0,4,3]) eq 0 then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str4 length+preserve', /continue
   nb_pbs=nb_pbs+1
endif
;
; few tests on str7
;
pos=STRSPLIT(str7,';',length=length)
if ARRAY_EQUAL(length,[2,9,10]) eq 0 then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str7 length', /continue
   nb_pbs=nb_pbs+1
endif
if ARRAY_EQUAL(pos,[0,3,13]) eq 0 then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str7 pos', /continue
   nb_pbs=nb_pbs+1
endif
;
pos=STRSPLIT(str7,'.',length=length)
if ARRAY_EQUAL(length,23) eq 0 then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str7 length no sep', /continue
   nb_pbs=nb_pbs+1
endif
if (pos NE 0) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str7 pos no sep', /continue
   nb_pbs=nb_pbs+1
endif
;
pos=STRSPLIT(str7,'a',length=length)
if ARRAY_EQUAL(length,[1,1]) eq 0 then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str7 length <<a>>', /continue
   nb_pbs=nb_pbs+1
endif
if ARRAY_EQUAL(pos,[2,12]) eq 0 then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str7 pos <<a>>', /continue
   nb_pbs=nb_pbs+1
endif
;
pos=STRSPLIT(str7,'a',length=length,/preserve)
exp_pos=LONG([INDGEN(3), INDGEN(9)+4,INDGEN(10)+14])
exp_len=LONARR(STRLEN(str7)-1)
exp_len[2]=1
exp_len[2+9]=1
;
if ARRAY_EQUAL(pos,exp_pos) eq 0 then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str7 pos <<a>> /preserve', /continue
   nb_pbs=nb_pbs+1
endif
if ARRAY_EQUAL(length,exp_len) eq 0 then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str7 length <<a>> /preserve', /continue
   nb_pbs=nb_pbs+1
endif
;
; tests on str5
;
tab=STRSPLIT(str5,/preserve)
res=LONG([0,1,3,5,7,9])
if (ARRAY_EQUAL(tab,res, /NO_TYPECONV) eq 0) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str5 preserve_null', /continue
   nb_pbs=nb_pbs+1
endif

ext=STRSPLIT(str5,/preserve,/extract)
if((N_ELEMENTS(ext) ne 6)) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error str5 extract , preserve', /continue
   nb_pbs=nb_pbs+1
endif
;
; bug found via STR_SEP 3286746 in the Patch section
;
tab=STR_SEP('ahasadfasdf','dfa')
res=['ahasa','sdf']
if (ARRAY_EQUAL(tab,res, /NO_TYPECONV) eq 0) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error STR_SEP 3286746', /continue
   nb_pbs=nb_pbs+1
endif
tab=STR_SEP('ahasadfasdfa','dfa')
res=['ahasa','s','']
if (ARRAY_EQUAL(tab,res, /NO_TYPECONV) eq 0) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'error STR_SEP 3286746, /preserve_null', /continue
   nb_pbs=nb_pbs+1
endif
;
; we must return a STRARR except for ''
;
tab=STRSPLIT('eeee','e',/extract)
if (SIZE(tab, /n_dim) NE 0) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'bad size for void string', /continue
   nb_pbs=nb_pbs+1
endif
tab=STRSPLIT('eeABCee','e',/extract)
if (SIZE(tab, /n_dim) NE 1) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'bad size for 1D STRARR', /continue
   nb_pbs=nb_pbs+1
endif
if (tab NE 'ABC') then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'bad value in STRARR', /continue
   nb_pbs=nb_pbs+1
endif
tab=STRSPLIT('eeABCeeABCee','e',/extract)
if (SIZE(tab, /n_dim) NE 1) then begin
   if KEYWORD_SET(verbose) then MESSAGE, 'bad size for 2 elements STRARR', /continue
   nb_pbs=nb_pbs+1
endif
;
vide=''
if (EXECUTE('res=TEXTOIDL(vide)') EQ 0) then begin
   MESSAGE, /CONTINUE, "Missing TexToIDL in your GDL_PATH or IDL_PATH"
endif else begin
   res='!7l!X!U2!N'
   tab=TEXTOIDL('\mu^2')
   if (res NE tab) then begin
      if KEYWORD_SET(verbose) then MESSAGE, 'error when using TexToIDL', /continue
      nb_pbs=nb_pbs+1
   endif
endelse
;
line="======================================="
MESSAGE, /Continue, line
MESSAGE, /Continue, " "
mess=' errors encoutered during STRSPLIT tests'
if (nb_pbs GT 0) then mess=STRING(nb_pbs)+mess else mess='NO'+mess
MESSAGE, /Continue, mess
MESSAGE, /Continue, " "
MESSAGE, /Continue, line
;
; if /debug OR /test nodes, we don't want to exit
if (nb_pbs GT 0) then begin
    if ~(KEYWORD_SET(debug) or KEYWORD_SET(test) or KEYWORD_SET(no_exit)) then EXIT, status=1
endif
;
if KEYWORD_SET(test) then STOP
;
end



