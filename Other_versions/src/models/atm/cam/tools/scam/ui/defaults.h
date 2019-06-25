/*------------------------------------------------------------------------*
 * File: defaults.h 
 * $Author: cam_titan $
 * $Id: defaults.h 62 2008-04-23 22:59:18Z cam_titan $ *
 *------------------------------------------------------------------------*/

#ifndef DEFAULTS_H
#define DEFAULTS_H

#include <string>
#include <vector>
#include <fstream>
#include "realtype.h"
#include "max.h"
#include "ioerr.h"

class Defaults
{
public:
    enum OpenMode { READ, WRITE };

    Defaults( const string& fileName, OpenMode mode = READ ) throw ( DfltsErr );
    ~Defaults(); 
    
    string  GetStringDefault(  const string& defaultname ) throw ( DfltsErr );
    int     GetIntDefault(  const string& defaultname ) throw ( DfltsErr );
    real_t   GetRealDefault(  const string& defaultname ) throw ( DfltsErr );

    void  WriteDefault( const string& defaultname, int    defaultvalue ) throw ( DfltsErr );
    void  WriteDefault( const string& defaultname, real_t  defaultvalue ) throw ( DfltsErr );
    void  WriteDefault( const string& defaultname, const string&  defaultvalue ) throw ( DfltsErr );

private:

    string  _filename;
    OpenMode _mode;
    vector< string > defaultsText;
    ofstream *output;
};



#endif /* DEFAULTS_H */
 
