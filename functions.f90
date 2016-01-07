Module functions

  Implicit none
     
contains

  Subroutine generate_gaussian_cmb_map()

    use healpix_types
    use pix_tools, only: nside2npix
    use fitstools, only: getsize_fits, input_map, output_map
    use arrays
    use fiducial

    Implicit none

    Integer*4 :: p
    Character(len=8),parameter :: fmt = '(I4.4)'
    Character(len=4) :: x

    Do p=1,number_of_cmb_simulations

       write(x,fmt) p

       call write_parameter_file_synfast(p)

       call system('synfast -d '//trim(PATH_TO_SYNFAST_PARAMETER_FILE)//''//trim(x)//'.par')

    End Do

  End subroutine generate_gaussian_cmb_map

  Subroutine write_parameter_file_synfast(iseed)

    use fiducial
    Implicit none

    Integer*4 :: iseed
    Character(len=8),parameter :: fmt = '(I4.4)'
    Character(len=4) :: x

    write(x,fmt) iseed

    open(UNIT_SYNFAST_PAR_FILE,file=PATH_TO_SYNFAST_PARAMETER_FILE//trim(x)//'.par')

    write(UNIT_SYNFAST_PAR_FILE,*) 'simul_type = ', simul_type

    write(UNIT_SYNFAST_PAR_FILE,*) 'nsmax = ', nsmax

    write(UNIT_SYNFAST_PAR_FILE,*) 'nlmax = ', nlmax

    write(UNIT_SYNFAST_PAR_FILE,*) 'infile = ', infile

    write(UNIT_SYNFAST_PAR_FILE,*) 'iseed = ', iseed

    write(UNIT_SYNFAST_PAR_FILE,*) 'fwhm_arcmin = ', fwhm_arcmin

    write(UNIT_SYNFAST_PAR_FILE,*) 'beam_file = ', beam_file

    write(UNIT_SYNFAST_PAR_FILE,*) 'almsfile = ', almsfile

    write(UNIT_SYNFAST_PAR_FILE,*) 'plmfile = ', plmfile

    write(UNIT_SYNFAST_PAR_FILE,*) 'outfile = ', PATH_TO_CMB_MAPS//trim(x)//'.fits'

    write(UNIT_SYNFAST_PAR_FILE,*) 'outfile_alms = ', outfile_alms

    write(UNIT_SYNFAST_PAR_FILE,*) 'winfiledir = ', PATH_TO_HEALPIX_DATA

    close(UNIT_SYNFAST_PAR_FILE)

  End subroutine write_parameter_file_synfast

End module functions
