/***************************************************************************
                          list.hpp  - for LIST objects
                             -------------------
    begin                : July 22 2013
    copyright            : (C) 2013 by M. Schellens et al.
    email                : m_schellens@users.sf.net

***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef LIST_HPP_
#define LIST_HPP_

// #include <envt.hpp>
void LIST__ToStream( DStructGDL* oStructGDL, std::ostream& o, SizeT w, SizeT* actPosPtr);

DStructGDL*GetSELF( BaseGDL* selfP, EnvUDT* e);

namespace lib {

   BaseGDL* list_fun( EnvT* e);

   void list__cleanup( EnvUDT* e);

   BaseGDL* LIST___OverloadIsTrue( EnvUDT* e);

   BaseGDL* LIST___OverloadBracketsRightSide( EnvUDT* e);
   void LIST___OverloadBracketsLeftSide( EnvUDT* e);
   BaseGDL* LIST___OverloadPlus( EnvUDT* e);
   BaseGDL* LIST___OverloadEQOp( EnvUDT* e);
   BaseGDL* LIST___OverloadNEOp( EnvUDT* e);

   void list__add( EnvUDT* e);
   
   void list__remove_pro( EnvUDT* e);

   void list__reverse( EnvUDT* e);
   
   BaseGDL* list__remove_fun( EnvUDT* e);

   BaseGDL* list__toarray( EnvUDT* e);

   BaseGDL* list__isempty( EnvUDT* e);
   BaseGDL* list__count( EnvUDT* e);
   BaseGDL* list__where( EnvUDT* e);
}

#endif
