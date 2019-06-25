/*------------------------------------------------------------------------*
 * File: map.h 
 * $Author: hpc $
 * $Id: globalmap.h 17 2006-12-11 21:50:24Z hpc $ *
 *------------------------------------------------------------------------*/

#ifndef Map_included
#define Map_included

#include <qframe.h>
#include <qpixmap.h>
#include <qpainter.h>
#include <realtype.h>

class SelectGlobalDataDlgImpl;
class Manager;

class Map : public QFrame
{
    friend class SelectGlobalDataDlgImpl;

public:
               Map( SelectGlobalDataDlgImpl* globalDlg, const char* name = NULL );
    void       DrawContents( QPainter* );
    void       mousePressEvent( QMouseEvent* );
    void       DrawLatLon( QPainter* painter );
    void       DrawHighlightClmn( bool highlight = TRUE );
    void       SetHighlightClmn( real_t lat, real_t lon );
    void       SetNumLatsLons( int nlat, int nlon ) { numLats = nlat,
                                                          numLons = nlon; }
private:
    
    void       paintEvent( QPaintEvent * );
    SelectGlobalDataDlgImpl* theGlobalDlg;
    QPixmap*   theMapPM;
    int        xHCol;
    int        yHCol;
    int        numLats, numLons;
    int        vertIncr;
    int        horizIncr;
    int        frameBorder;   // thickness of the frame border
    void       DrawHighlightClmn( QPainter* painter );
    void       SetIncrements();
    
};
    
#endif // Map_included

