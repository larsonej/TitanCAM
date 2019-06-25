/*------------------------------------------------------------------------*
 * File: dataset.h 
 * $Author: cam_titan $
 * $Id: dataset.h 62 2008-04-23 22:59:18Z cam_titan $ *
 *------------------------------------------------------------------------*/
#ifndef DATASET_H
#define DATASET_H

#include "realtype.h"
#include "ncfile.h"


const int MAX_DATASET_VARS = 500;

class Dataset : public ncfile::NcFile
{

public:
    Dataset( string filename, int type )
        : NcFile( filename, READ ), filetype( type ) {}
    Dataset::Dataset( int type )
        : NcFile(), filetype( type ) {}
    int FindIdx( string varname, real_t value );
    bool CheckLat( real_t lat, int latIdx );
    bool CheckLon( real_t lon, int lonIdx );
    bool CheckVars( string& missingVars );
    static string TypeDescription( int type );
    static string TypeName( int type );
    string lsmpftfile;
    
private:
    int   numvars;
    int   filetype;
};


#endif // DATASET_H



