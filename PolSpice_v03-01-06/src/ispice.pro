pro paramfile2thincscript, paramfile, module, thinc_script, omp=omp

jobname = 'myJob'
nickname = module

openw, lunit2, thinc_script, /get_lun
printf, lunit2, 'from thinc import *'
printf, lunit2, ' '
printf, lunit2, 'BeginPipe()'
printf, lunit2, jobname+' = NewJob("'+module+'", label = "'+nickname+'", cast = True)'

; nl = numlines(paramfile) for IDL before 5.6
nl = file_lines(paramfile)
openr, lunit1, paramfile,    /get_lun
line_in = ''
line_out = ''
for i=0, nl-1 do begin
    readf, lunit1, line_in
    line_in = strtrim(line_in,2)
    if (line_in ne '' && strmid(line_in,0,1) ne '#' ) then begin
        parse = str_sep(line_in,'=')
        kwd = strtrim(parse[0],2)
        val = strtrim(parse[1],2)

        line_out = 'Set('+jobname+', "'+kwd+'", "'+val+'")'
        ;;;;;;print, thinc_script, line_out
        printf, lunit2, line_out
    endif
endfor

extra = ''
if keyword_set(omp) then begin
    extra += ', exportEnviron={"OMP_NUM_THREADS":"'+strtrim(omp,2)+'"}'
endif
printf,lunit2, 'Submit('+jobname+extra+')'
printf,lunit2, 'EndPipe()'
free_lun, lunit2

return
end

function spice_yes_no, kwd
; kwd not set -> NO
; kwd="NO"    -> NO
; kwd="YES"   -> YES
if keyword_set(kwd) then begin ; kwd set, test for NO
    yes_no = (strupcase(kwd) eq 'NO')? 'NO' : 'YES'
endif else begin ; kwd not set -> NO
    yes_no = 'NO'
endelse
return, yes_no
end


pro ispice, mapfile1, clfile $
            , apodizesigma	 =  apodizesigma $
            , apodizetype	 =  apodizetype $
            , binpath            = binpath $
            , fwhm1		 =  fwhm1 $
            , beam_file1	 =  beam_file1 $
            , fwhm2		 =  fwhm2 $
            , beam_file2	 =  beam_file2 $
            , cl_outmap_file	 =  cl_outmap_file $
            , cl_inmap_file	 =  cl_inmap_file $
            , cl_outmask_file	 =  cl_outmask_file $
            , cl_inmask_file	 =  cl_inmask_file $
            , corfile		 =  corfile $
            , covfileout	 =  covfileout $
            , decouple		 =  decou_usr $
;            , dry		 =  dry $
            , fits_out		 =  fits_out_usr $
            , dmc                = dmc $
            , help               = help $
            , keep_tmp_files     = keep_tmp_files $
            , kernelsfileout	 =  kernelsfileout $
            , lmw1               = lmw1 $  
            , lmw2               = lmw2 $  
            , mapfile2		 =  mapfile2 $
            , maskfile1		 =  maskfile1 $
            , maskfile2		 =  maskfile2 $
            , maskfilep1	 =  maskfilep1 $
            , maskfilep2	 =  maskfilep2 $
            , nlmax	         =  nlmax $
            , normfac		 =  normfac $
            , npairsthreshold	 =  npairsthreshold $
            , noisecorfile	 =  noisecorfile $
            , noiseclfile	 =  noiseclfile $
;            , overwrite		 =  overwrite $
            , pixelfile		 =  pixelfile_usr $
            , polarization	 =  polar_usr $
            , show_cl            = show_cl $
            , silent             = silent $
            , subav		 =  subav_usr $
            , subdipole		 =  subdipole_usr $
            , symmetric_cl       =  symm_usr $
            , tf_file		 =  tf_file $
            , thetamax		 =  thetamax $
            , tolerance		 =  tolerance $
;;            , verbosity		 =  verbosity $
            , weightfile1	 =  weightfile1 $
            , weightfilep1	 =  weightfilep1 $
            , weightfile2	 =  weightfile2 $
            , weightfilep2	 =  weightfilep2 $
            , weightpower1	 =  weightpower1 $
            , weightpowerp1	 =  weightpowerp1 $
            , weightpower2	 =  weightpower2 $
            , weightpowerp2	 =  weightpowerp2 $
            , windowfilein	 =  windowfilein $
            , windowfileout	 =  windowfileout $
            , xtramapfile1       =  xtramapfile1 $
            , xtramapfile2       =  xtramapfile2 
;+
; NAME:
;        ISPICE
;
; PURPOSE:
;        Interface to Spice C(l) (angular power spectrum) estimator
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;        ispice, map [ clfile , 
;    apodizesigma=, apodizetype=, binpath=, dmc=, fwhm1=,
;    beam_file1=, fwhm2=, beam_file2=, 
;    cl_outmap_file=, cl_inmap_file=, cl_outmask_file=, cl_inmask_file=, corfile=, covfileout=
;    decouple=, dmc=, fits_out =, help=, keep_tmp_files =, kernelsfileout =,    lmw1=, lmw2=, 
;    mapfile2=, maskfile1=, maskfile2=, maskfilep1=, maskfilep2=, 
;    nlmax=, normfac=, npairsthreshold=, noisecorfile=, noiseclfile=, 
;    pixelfile=, polarization=, show_cl=, silent=, subav=, subdipole=, symmetric_cl=,
;    tf_file=, thetamax=, tolerance=, weightfile1=, weightfilep1=, weightfile2=, weightfilep2=, 
;    weightpower1=, weightpowerp1=, weightpower2=, weightpowerp2=, windowfilein=, windowfileout=,
;    xtramapfile1=, xtramapfile2= ]
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;  creation: 2008-03-02, EH, IAP
;  2009-01-27: outputs kernels
;  2009-07-23: pixelfile='NO' will not correct for pixel window 
;  2010-01-08: bug correction of the above
;  2010-01-11: added subdipole, 
;              can process Planck-HFI data in DMC format (DMC= keyword)
;  2010-05: added covfileout
;  2012-03: added xtramapfile1 and xtramapfile2
;  2014-04-01: support linear combination of input maps on the fly:
;     allow map1_in and mapfile2 to be vector of FITS files or objects
;     added lmw1 and lmw2 (wrappers for listmapweights1_* and listmapweights2_)
; 2014-04-16: keywords DECOUPLE, POLARIZATION, SUBAV, SUBDIPOLE, SYMMETRIC_CL  now accept
;            "YES" and "NO" values (as well as 1 and 0)
; 2015-03-02: bug correction on lmw1 and lmw2
; 2015-08-26: added TOLERANCE; set explicitely FITS_OUT=0 to prevent FITS output
;
;-

do_dmc = keyword_set(dmc)
options = do_dmc ? '' : '-optinfile'
local = {routine: 'ispice', exe: 'spice', options: options}
syntax = [local.routine+', map1_in, [clfile ,' $
,'    apodizesigma=, apodizetype=, binpath=, dmc=, fwhm1=,' $
,'    beam_file1=, fwhm2=, beam_file2=, ' $
,'    cl_outmap_file=, cl_inmap_file=, cl_outmask_file=, cl_inmask_file=, corfile=, covfileout=,' $
,'    decouple=, fits_out=, help=, keep_tmp_files =, kernelsfileout =,       lmw1 =, lmw2 =, ' $
,'    mapfile2=, maskfile1=, maskfile2=, maskfilep1=, maskfilep2=, ' $
,'    nlmax=, normfac=, npairsthreshold=, noisecorfile=, noiseclfile=, ' $
,'    pixelfile=, polarization=, show_cl=, silent=, subav=, subdipole=, symmetric_cl=' $
,'    tf_file=, thetamax=, tolerance=, weightfile1=, weightfilep1=, weightfile2=, weightfilep2=, ' $
,'    weightpower1=, weightpowerp1=, weightpower2=, weightpowerp2=, windowfilein=, windowfileout=' $
,'    xtramapfile1=, xtramapfile2= ]' ]

if keyword_set(help) then begin
    doc_library,local.routine
    return
endif

if (n_params() eq 0) then begin
    print,syntax,form='(a)'
    print
    print,local.routine+',/help     for extended help'
    return
endif
if (n_params() gt 2) then begin
    print,syntax,form='(a)'
    message,'wrong number of arguments'
endif

if (size(mapfile1,/tname) ne 'STRING') then begin
    message,'First argument (input map) must be a FITS file',/info
    print
    print,syntax
    return
endif
multimap1 = (n_elements(mapfile1) gt 1)
multimap2 = (n_elements(mapfile2) gt 1)

; polarization = keyword_set(polar_usr) ? 'YES' : 'NO'
; decouple     = keyword_set(decou_usr) ? 'YES' : 'NO'
; symmetric_cl = keyword_set(symm_usr)  ? 'YES' : 'NO'
; subav        = keyword_set(subav_usr) ? 'YES' : 'NO'
; subdipole    = keyword_set(subdipole_usr) ? 'YES' : 'NO'
polarization = spice_yes_no(polar_usr) 
decouple     = spice_yes_no(decou_usr) 
symmetric_cl = spice_yes_no(symm_usr)  
subav        = spice_yes_no(subav_usr) 
subdipole    = spice_yes_no(subdipole_usr) 
pixelfile    = (keyword_set(pixelfile_usr) || ~(keyword_set(noisecorfile) || keyword_set(noiseclfile)) ) ? 'YES' : 'NO'
if (defined(pixelfile_usr) && strupcase(pixelfile_usr) eq 'NO') then pixelfile = 'NO'
fits_out     = defined(fits_out_usr) ? spice_yes_no(fits_out_usr) : "YES"

;-------------------
hpx_xface_generic, fullpath, tmp_par_file, binpath, init=local, cxx=cxx
NoFile = keyword_set(cxx) ? " " : " '' "

; deal with online data
; tmp_clfile   = hpx_mem2file((arg_present(clfile) || defined(clfile)) ? clfile : NoFile, /out)
tmp_clfile   = hpx_mem2file((arg_present(clfile) || defined(clfile)) ? clfile : 'NO', /out)

; writes parameter file
openw,lunit,tmp_par_file, /get_lun
printf,lunit,'# parameter file for IDL interface to '+fullpath
printf,lunit,'# written: '+systime()+' by '+local.routine
printf,lunit,' '

printf,lunit,hpx_add_parameter('apodizesigma',	apodizesigma,	skip_if_not_set=~do_dmc, default="NO") ; let code choose default
printf,lunit,hpx_add_parameter('apodizetype' ,	apodizetype,	skip_if_not_set=~do_dmc, default=0)
printf,lunit,hpx_add_parameter('beam',          fwhm1   ,	skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('beam_file',	beam_file1,	/skip_if_not_set) 
printf,lunit,hpx_add_parameter('beam2',         fwhm2,   	skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('beam_file2',    beam_file2,	/skip_if_not_set)
printf,lunit,hpx_add_parameter('clfile',        tmp_clfile,  	/skip_if_not_set)
printf,lunit,hpx_add_parameter('cl_outmap_file',  cl_outmap_file,	      /skip_if_not_set)
printf,lunit,hpx_add_parameter('cl_inmap_file',   cl_inmap_file,	      /skip_if_not_set)
printf,lunit,hpx_add_parameter('cl_outmask_file', cl_outmask_file,	      /skip_if_not_set)
printf,lunit,hpx_add_parameter('cl_inmask_file',  cl_inmask_file,	      /skip_if_not_set)
printf,lunit,hpx_add_parameter('corfile',	  corfile,  	              /skip_if_not_set)
printf,lunit,hpx_add_parameter('covfileout',      covfileout,	              /skip_if_not_set)
printf,lunit,hpx_add_parameter('decouple',	  decouple,     skip_if_not_set=~do_dmc, default="NO")
if (~do_dmc) then begin
    printf,lunit,hpx_add_parameter('dry',   	  'NO')
    ;printf,lunit,hpx_add_parameter('fits_out',	  'YES')
    printf,lunit,hpx_add_parameter('fits_out',	  fits_out)
endif
printf,lunit,hpx_add_parameter('kernelsfileout',  kernelsfileout,	     /skip_if_not_set)
;
if (multimap1) then begin
    if (do_dmc && keyword_set(polar_usr)) then begin
        if (n_elements(mapfile1) mod 3 ne 0) then message,'Expected (I,Q,U) triplet in map'
        for i=0,n_elements(mapfile1)-1,3 do begin
            k = i/3+1
            printf,lunit,hpx_add_parameter('listmapfiles1_' +strtrim(k,2), mapfile1[i],   /expand)
            printf,lunit,hpx_add_parameter('listmapfilesQ1_'+strtrim(k,2), mapfile1[i+1], /expand)
            printf,lunit,hpx_add_parameter('listmapfilesU1_'+strtrim(k,2), mapfile1[i+2], /expand)
        endfor        
    endif else begin
        for i=0,n_elements(mapfile1)-1 do begin
            ;print,'listmapfiles1_'+strtrim(i+1,2)+' = '+mapfile1[i]
            printf,lunit,hpx_add_parameter('listmapfiles1_'+strtrim(i+1,2), mapfile1[i], /expand)
        endfor
    endelse
endif else begin
    printf,lunit,hpx_add_parameter('mapfile',         mapfile1,      /expand)
    printf,lunit,hpx_add_parameter('extramapfile',    xtramapfile1,          /skip_if_not_set)
endelse
if (multimap2) then begin
    if (do_dmc && keyword_set(polar_usr)) then begin
        if (n_elements(mapfile2) mod 3 ne 0) then message,'Expected (I,Q,U) triplet in mapfile2'
        for i=0,n_elements(mapfile2)-1,3 do begin
            k = i/3+1
            printf,lunit,hpx_add_parameter('listmapfiles2_' +strtrim(k,2), mapfile2[i],   /expand)
            printf,lunit,hpx_add_parameter('listmapfilesQ2_'+strtrim(k,2), mapfile2[i+1], /expand)
            printf,lunit,hpx_add_parameter('listmapfilesU2_'+strtrim(k,2), mapfile2[i+2], /expand)
        endfor        
    endif else begin
        for i=0,n_elements(mapfile2)-1 do begin
            ;print,'listmapfiles2_'+strtrim(i+1,2)+' = '+mapfile2[i]
            printf,lunit,hpx_add_parameter('listmapfiles2_'+strtrim(i+1,2), mapfile2[i], /expand)
        endfor
    endelse
endif else begin
    printf,lunit,hpx_add_parameter('mapfile2',	      mapfile2,		     /skip_if_not_set)
    printf,lunit,hpx_add_parameter('extramapfile2',   xtramapfile2,	     /skip_if_not_set)
endelse
if defined(lmw1) then begin
;    for i=0, n_elements(lmw1)-1 do printf,lunit,hpx_add_parameter('listmapweights1_'+strtrim(i+1,2), lmw1[i+1],   /skip_if_not_set, default=1.0)
    for i=0, n_elements(lmw1)-1 do printf,lunit,hpx_add_parameter('listmapweights1_'+strtrim(i+1,2), lmw1[i],   /skip_if_not_set, default=1.0)
endif
if defined(lmw2) then begin
;    for i=0, n_elements(lmw2)-1 do printf,lunit,hpx_add_parameter('listmapweights2_'+strtrim(i+1,2), lmw2[i+1],   /skip_if_not_set, default=1.0)
    for i=0, n_elements(lmw2)-1 do printf,lunit,hpx_add_parameter('listmapweights2_'+strtrim(i+1,2), lmw2[i],   /skip_if_not_set, default=1.0)
endif
;
printf,lunit,hpx_add_parameter('maskfile',	  maskfile1,		     /skip_if_not_set)
printf,lunit,hpx_add_parameter('maskfile2',	  maskfile2,		     /skip_if_not_set)
printf,lunit,hpx_add_parameter('maskfilep',	  maskfilep1,	             /skip_if_not_set)
printf,lunit,hpx_add_parameter('maskfilep2',	  maskfilep2,	             /skip_if_not_set)
printf,lunit,hpx_add_parameter('nlmax',		  nlmax,	   skip_if_not_set=~do_dmc, default=-1)
printf,lunit,hpx_add_parameter('normfac',	  normfac,         skip_if_not_set=~do_dmc, default=1.0)
;printf,lunit,hpx_add_parameter('npairsthreshold', npairsthreshold, skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('npairsthreshold', npairsthreshold, skip_if_not_set=~do_dmc, default=0.0)
printf,lunit,hpx_add_parameter('noisecorfile',	  noisecorfile,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('noiseclfile',	  noiseclfile,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('overwrite',	  'YES')
printf,lunit,hpx_add_parameter('pixelfile',	  pixelfile,     skip_if_not_set=~do_dmc, default="YES")
printf,lunit,hpx_add_parameter('polarization',	  polarization,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('subav',		  subav,         skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('subdipole',	  subdipole,	 skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('symmetric_cl',	  symmetric_cl,  skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('tf_file',	  tf_file,		     /skip_if_not_set)
printf,lunit,hpx_add_parameter('thetamax',	  thetamax,      skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('tolerance',	  tolerance,     skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('verbosity',	  keyword_set(silent)?0:2)
printf,lunit,hpx_add_parameter('weightfile',	  weightfile1,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('weightfile2',	  weightfile2,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('weightfilep',	  weightfilep1,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('weightfilep2',	  weightfilep2,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('weightpower',	  weightpower1,   skip_if_not_set=~do_dmc, default=1.0)
printf,lunit,hpx_add_parameter('weightpower2',	  weightpower2,   skip_if_not_set=~do_dmc, default=1.0)
printf,lunit,hpx_add_parameter('weightpowerp',	  weightpowerp1,  skip_if_not_set=~do_dmc, default=1.0)
printf,lunit,hpx_add_parameter('weightpowerp2',   weightpowerp2,  skip_if_not_set=~do_dmc, default=1.0)
printf,lunit,hpx_add_parameter('windowfilein',	  windowfilein,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('windowfileout',   windowfileout,	     /skip_if_not_set)
free_lun, lunit

; execute command
if (do_dmc) then begin
    ;thinc_file = '/tmp/hivon_thinc.py'
    thinc_file = getenv('IDL_TMPDIR')+'spice_from_idl_'+string(long(systime(1)*100 mod 1.e8), form='(i8.8)')+'.py'
    paramfile2thincscript, tmp_par_file, 'spice', thinc_file, omp=8
    spawn, 'python '+thinc_file
endif else begin
    hpx_xface_generic, /run, fullpath, tmp_par_file, silent=silent
endelse

; deal with online data
if (arg_present(clfile) || defined(clfile)) then hpx_file2mem, tmp_clfile, clfile, /cl, show_cl = show_cl

; to_remove
hpx_xface_generic, clean = ~keyword_set(keep_tmp_files)


return
end
