// Copyright (C) 2008  Werner Smekal
//
// This file is part of PLplot.
//
// PLplot is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published
// by the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// PLplot is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with PLplot; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
//

#ifndef __WXWIDGETS_H__
#define __WXWIDGETS_H__

#include <vector>
#include <memory>

// plplot headers
#include "plplotP.h"
#include "wxwidgets_comms.h"

// some special wxWidgets headers
#include <wx/wx.h>
#include <wx/spinctrl.h>
#ifdef wxUSE_GRAPHICS_CONTEXT
#include <wx/dcgraph.h>
#endif

class wxPLplotFrame;

class FontGrabber
{
public:
    FontGrabber();
    wxFont GetFont( PLUNICODE fci, PLFLT scaledFontSize, bool underlined );
    bool lastWasCached( ){ return m_lastWasCached; }
private:
    wxFont    m_prevFont;
    PLUNICODE m_prevFci;
    PLFLT     m_prevScaledFontSize;
    bool      m_prevUnderlined;
    bool      m_lastWasCached;
};

// base device class
class wxPLDevice
{
public:
    wxPLDevice( PLStream *pls, char * mfo, PLINT text, PLINT hrshsym );
    virtual ~wxPLDevice( void );

    void DrawLine( short x1a, short y1a, short x2a, short y2a );
    void DrawPolyline( short *xa, short *ya, PLINT npts );
    void ClearBackground( PLStream* pls, PLINT x1 = -1, PLINT y1 = -1, PLINT x2 = -1, PLINT y2 = -1 );
    void FillPolygon( PLStream *pls );
    void SetWidth( PLStream *pls );
    void SetColor( PLStream *pls );
    void SetDC( PLStream *pls, wxDC* dc );
    void EndPage( PLStream* pls );
    void BeginPage( PLStream* pls );
    void SetSize( PLStream* pls, int width, int height );
    void ProcessString( PLStream* pls, EscText* args );
    void FixAspectRatio( bool fix );
    void Locate( PLStream* pls, PLGraphicsIn *graphicsIn );
    void Flush( PLStream* pls );
    void PreDestructorTidy( PLStream *pls );

private:
    void DrawTextLine( PLUNICODE* ucs4, int ucs4Len, PLFLT baseFontSize, bool drawText, PLINT &superscriptLevel, bool &underlined );
    void DrawTextSection( char* utf8_string, PLFLT baseFontSize, PLFLT yOffset, bool underlined, bool drawText );
    void TransmitBuffer( PLStream* pls, unsigned char transmissionType );
    void SetupMemoryMap();
    wxRegion GetClipRegion();

    //The DC we will draw on if given by the user
    wxDC *m_dc;
    bool m_useDcTextTransform;
    //for the gcdc case we may need to store the graphics context for use
    // with text transformations
    wxGraphicsContext *m_gc;
    wxPen             m_pen;
    wxBrush           m_brush;

    //Size and Scale
    //As far as plplot is concerned the size of the window is SHRT_MAX by
    //SHRT_MAX which gives us the best resolution.
    const PLFLT m_plplotEdgeLength;
    PLFLT       m_width;   //native width
    PLFLT       m_height;  //native height
    PLFLT       m_xScale;  //conversion from native width to plplotEdgeLength
    PLFLT       m_yScale;  //conversion from native height to plplotEdgeLength
    PLFLT       m_xAspect; //values which when multiplied by m_plplotEdgeLength give an aspect
    PLFLT       m_yAspect; //ratio equal to the native aspect ratio, the biggest of which is 1.0
    PLFLT       m_scale;   //MAX(m_scalex, m_scaley)
    bool        m_fixedAspect;

    // font variables
    static const int m_max_string_length = 500;
    //bool m_underlined;
    FontGrabber      m_fontGrabber;
    wxCoord          m_textWidth, m_textHeight, m_textDescent, m_textLeading;
    PLUNICODE        m_fci;

    //memory of previous single character string (likely plot point)
    PLUNICODE m_prevSingleCharString;
    PLINT     m_prevSingleCharStringWidth;
    PLINT     m_prevSingleCharStringHeight;

    //Text positioning related variables
    wxCoord m_superscriptHeight;          //distance between superscript top and baseline
    wxCoord m_subscriptDepth;             //distance between subscript base and baseline
    PLFLT   m_lineSpacing;
    //PLFLT   m_yOffset;
    PLINT   m_posX;
    PLINT   m_posY;
    PLFLT   m_rotation;

    //variables for dealing with sending/receiving commands
    //via a memory map
    char         m_mfo[PLPLOT_MAX_PATH];
    PLNamedMutex m_mutex;
    size_t       m_localBufferPosition;
    PLMemoryMap  m_outputMemoryMap;
};


struct dev_entry
{
    wxString dev_name;
    wxString dev_menu_short;
    wxString dev_menu_long;
    wxString dev_file_app;
    bool     pixelDevice;
};

//--------------------------------------------------------------------------
//  Declarations for the device.
//--------------------------------------------------------------------------

void plD_init_wxwidgets( PLStream * );
void plD_init_wxpng( PLStream * );
void plD_line_wxwidgets( PLStream *, short, short, short, short );
void plD_polyline_wxwidgets( PLStream *, short *, short *, PLINT );
void plD_eop_wxwidgets( PLStream * );
void plD_bop_wxwidgets( PLStream * );
void plD_tidy_wxwidgets( PLStream * );
void plD_state_wxwidgets( PLStream *, PLINT );
void plD_esc_wxwidgets( PLStream *, PLINT, void * );

void wx_set_dc( PLStream* pls, wxDC* dc );
void wx_set_buffer( PLStream* pls, wxImage* buffer );

//--------------------------------------------------------------------------
//  Debug functions
//--------------------------------------------------------------------------

// define if you want debug output
// #define _DEBUG //
// #define _DEBUG_VERBOSE //
void Log_Verbose( const char *fmt, ... );
void Log_Debug( const char *fmt, ... );


//--------------------------------------------------------------------------
// Font style and weight lookup tables
//--------------------------------------------------------------------------
const wxFontFamily fontFamilyLookup[5] = {
    wxFONTFAMILY_SWISS,      // sans-serif
    wxFONTFAMILY_ROMAN,      // serif
    wxFONTFAMILY_TELETYPE,   // monospace
    wxFONTFAMILY_SCRIPT,     // script
    wxFONTFAMILY_SWISS       // symbol
};

const int          fontStyleLookup[3] = {
    wxFONTSTYLE_NORMAL,      // upright
    wxFONTSTYLE_ITALIC,      // italic
    wxFONTSTYLE_SLANT        // oblique
};

const int          fontWeightLookup[2] = {
    wxFONTWEIGHT_NORMAL,     // medium
    wxFONTWEIGHT_BOLD        // bold
};

#endif // __WXWIDGETS_H__
