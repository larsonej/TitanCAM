/*------------------------------------------------------------------------*
 * File: sicfile.h 
 * $Author: hpc $
 * $Id: sicfile.h 19 2007-02-16 19:32:47Z hpc $
 *------------------------------------------------------------------------*/

#ifndef SICFILE_H
#define SICFILE_H

#include "realtype.h"
#include "ncfile.h"

class Sicfile: public NcFile
{
public:
    Sicfile( const string& name, OpenMode mode, bool clobber );
    ~Sicfile();
};

#endif // SICFILE_H

