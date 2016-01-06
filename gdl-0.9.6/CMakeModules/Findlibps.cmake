find_path(LIBPSPKG libps.pc PATH_SUFFIXES lib lib/pkgconfig lib64/pkgconfig)
include(FindPackageHandleStandardArgs)
if(LIBPSPKG)
	set(ENV{PKG_CONFIG_PATH} ${LIBPSPKG}) # pkg search path
	include(FindPkgConfig)
	pkg_check_modules(LIBPS libps)
	if(LIBPS_FOUND)
		find_package_handle_standard_args(LIBPS DEFAULT_MSG LIBPS_LIBRARIES)
	endif(LIBPS_FOUND)
else(LIBPSPKG) # no libps.pc file
	find_library(LIBPS_LIBRARIES NAMES libps pslib)
	find_path(LIBPS_INCLUDE_DIRS NAMES libps/pslib.h)	
	find_package_handle_standard_args(LIBPS DEFAULT_MSG LIBPS_LIBRARIES)
endif(LIBPSPKG)

mark_as_advanced(
LIBPSPKG
LIBPS
LIBPS_LIBRARIES
)
