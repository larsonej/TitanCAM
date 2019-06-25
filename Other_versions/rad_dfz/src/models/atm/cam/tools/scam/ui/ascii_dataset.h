/*------------------------------------------------------------------------*
 * File: ascii_dataset.h 
 * $Author: hpc $
 * $Id: ascii_dataset.h 17 2006-12-11 21:50:24Z hpc $ *
 *------------------------------------------------------------------------*/
#ifndef ASCII_DATASET_H
#define ASCII_DATASET_H

#include "realtype.h"
#include "ncfile.h"

class Ascii_Dataset 
{

public:
    Ascii_Dataset( string filename, int type )
        : _name(filename),_type(type) {}
    Ascii_Dataset::Ascii_Dataset( int type )
        : _type(type) {}
    const string& name() const {return _name;}
    int type(){return _type;}
    
protected:
    string   _name;
    int   _type;
};

#endif // ASCII_DATASET_H



