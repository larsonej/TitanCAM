/*------------------------------------------------------------------------*
 * File: numlined.h 
 * $Author: hpc $
 * $Id: numlined.h 19 2007-02-16 19:32:47Z hpc $ *
 *------------------------------------------------------------------------*/

#ifndef NUMLINED_H
#define NUMLINED_H

#include <qkeycode.h>
#include <qlined.h>
#include <ctype.h>

#include "manager.h"
#include "timeconvert.h"

//
//  NumLineEdit is a subclass of QLineEdit that ignores non
//  numeric entries, i.e., you can't type alphabetic characters
//  into it.
//
class NumLineEdit : public QLineEdit
{
public:
    NumLineEdit( QWidget *parent=0, const char *name=0 ):
            QLineEdit( parent, name ) 
    { 
        QFont font( "helvetica", 18, 63, 0 );
        setFont( font );
    }
private:
    void keyPressEvent( QKeyEvent *e ) {
          // ignore alphabetic key entries
        if ( isalpha( e->ascii() ) ) {
            e->ignore();
            clearFocus();
            return;
        }
        else
            QLineEdit::keyPressEvent( e );
    }
      // set end step when the mouse cursor leaves the window
    void leaveEvent( QEvent* e  ) {
        int endstep;
        TimeConverter tc( text().ascii() );
        tc.Steps() > MANAGER.MaxStep() ? endstep = MANAGER.MaxStep() : endstep =  tc.Steps();
        MANAGER.SetEndStep( endstep); 
    }
};

   
#endif // NUMLINED_H
