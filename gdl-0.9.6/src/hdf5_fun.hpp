/***************************************************************************
                          hdf5_fun.hpp  -  HDF5 GDL library function
                             -------------------
    begin                : Aug 02 2004
    copyright            : (C) 2004 by Peter Messmer
    email                : messmer@users.sourceforge.net
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef HDF5_FUN_HPP_
#define HDF5_FUN_HPP_

#include "datatypes.hpp"
#include "envt.hpp"

/** 
 * needs the hdf5 library, obtainable from http://hdf.ncsa.uiuc.edu/HDF5/
 */
#define H5_USE_16_API 1
#include "hdf5.h"

namespace lib {

  /**
   * wrapper routines for hdf5 function calls. Only minimal set currently
   * implemented in order to get things going. The rest is bare-bones work
   */
  BaseGDL* h5f_open_fun( EnvT* e);
  BaseGDL* h5d_open_fun( EnvT* e);
  BaseGDL* h5d_read_fun( EnvT* e);
  BaseGDL* h5d_get_space_fun( EnvT* e);
  BaseGDL* h5s_get_simple_extent_dims_fun( EnvT* e);
  void h5f_close_pro( EnvT* e);
  void h5d_close_pro( EnvT* e);
  void h5s_close_pro( EnvT* e);

  // SA:
  BaseGDL* h5f_is_hdf5_fun( EnvT* e);
  BaseGDL* h5_get_libversion_fun( EnvT* e);
  BaseGDL* h5d_get_type_fun( EnvT* e);
  BaseGDL* h5t_get_size_fun( EnvT* e);
  BaseGDL* h5a_open_name_fun( EnvT* e);
  BaseGDL* h5a_open_idx_fun( EnvT* e);
  BaseGDL* h5a_get_name_fun( EnvT* e);
  BaseGDL* h5a_get_space_fun( EnvT* e);
  BaseGDL* h5a_get_type_fun( EnvT* e);
  BaseGDL* h5a_get_num_attrs_fun( EnvT* e);
  BaseGDL* h5a_read_fun( EnvT* e);
  BaseGDL* h5g_open_fun( EnvT* e);
  void h5a_close_pro( EnvT* e );
  void h5t_close_pro( EnvT* e );
  void h5g_close_pro( EnvT* e );

} // namespace

#endif

