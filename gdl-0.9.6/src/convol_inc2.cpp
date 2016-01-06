/***************************************************************************
                          convol_inc2.cpp  -  convol function EDGE_TRUNCATE
                             -------------------
    begin                : Sep 19 2004
    copyright            : (C) 2004 by Marc Schellens
    email                : m_schellens@users.sourceforge.net
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

// to be included from convol.cpp

#ifdef INCLUDE_CONVOL_INC_CPP

// for all result elements
for (SizeT a = 0; a < nA; ++aInitIx[1]) {
  bool regular = true;
  for (SizeT aSp = 1; aSp < nDim;) {
    if (aInitIx[ aSp] < this->dim[ aSp]) {
	  regArr[ aSp] = 
	    aInitIx[aSp] >= aBeg[aSp] && aInitIx[aSp] < aEnd[ aSp];

	  if( regular)
	    for(; aSp<nDim; ++aSp)
          if (!regArr[ aSp]) {
		  regular = false; 
		  break;
		}
	  break;
	}

      aInitIx[ aSp] = 0;
      regArr[ aSp] = !aBeg[ aSp];
      if( aBeg[ aSp]) regular = false; 
	  
      ++aInitIx[ ++aSp];
    }

  if (regular) {
      a += aBeg0;

      if( center)
      for (long aInitIx0 = aBeg0; aInitIx0 < aEnd0; ++aInitIx0, ++a) {
#if defined(CONVOL_BYTE__) || defined (CONVOL_UINT__) || defined (CONVOL_INT__)
	    DLong res_a = 0;
#else
	    Ty& res_a = (*res)[ a];
#endif

	    long* kIx = kIxArr;
        for (SizeT k = 0; k < nK; k += kDim0) {
		SizeT aLonIx = aInitIx0 + kIx[0];
          for (SizeT rSp = 1; rSp < nDim; ++rSp) aLonIx += (aInitIx[ rSp] + kIx[ rSp]) * aStride[ rSp];

		for( SizeT k0=0; k0<kDim0; ++k0) res_a += ddP[ aLonIx+k0] * ker[ k+k0]; 

		kIx += kDim0_nDim;
	    }

        res_a /= scale;
        res_a += bias;

#if defined(CONVOL_BYTE__) || defined (CONVOL_UINT__) || defined (CONVOL_INT__)
            CONVERT_CONVOL_TO_ORIG;
#endif
      } else
      for (long aInitIx0 = aBeg0; aInitIx0 < aEnd0; ++aInitIx0, ++a) {
#if defined(CONVOL_BYTE__) || defined (CONVOL_UINT__) || defined (CONVOL_INT__)
	    DLong res_a = 0;
#else
	    Ty& res_a = (*res)[ a];
#endif

	    long* kIx = kIxArr;
        for (SizeT k = 0; k < nK; k += kDim0) {
		SizeT aLonIx = aInitIx0 + kIx[0];
          for (SizeT rSp = 1; rSp < nDim; ++rSp) aLonIx += (aInitIx[ rSp] + kIx[ rSp]) * aStride[ rSp];

		for( SizeT k0=0; k0<kDim0; ++k0) res_a += ddP[ aLonIx-k0] * ker[ k+k0]; 

		kIx += kDim0_nDim;
	    }

        res_a /= scale;
        res_a += bias;
        
#if defined(CONVOL_BYTE__) || defined (CONVOL_UINT__) || defined (CONVOL_INT__)
            CONVERT_CONVOL_TO_ORIG;
#endif
	  }

      a += dim0_aEnd0;
    } // if( regular) // in dim 1-n
  else {
      a += dim0; // update a
    } 
} // for(...)

#endif
