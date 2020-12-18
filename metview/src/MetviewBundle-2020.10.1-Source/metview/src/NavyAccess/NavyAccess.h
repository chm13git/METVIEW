/***************************** LICENSE START ***********************************

 Copyright 2012 ECMWF and INPE. This software is distributed under the terms
 of the Apache License version 2.0. In applying this license, ECMWF does not
 waive the privileges and immunities granted to it by virtue of its status as
 an Intergovernmental Organization or submit itself to any jurisdiction.

 ***************************** LICENSE END *************************************/

#ifndef NAVYACCESS_H
#define NAVYACCESS_H

//*****************************************************************
//  Application NavyAccess
//
//  Access GRIB/BUFR/NETCDF data from the Navy's storage directories
//*****************************************************************

#include <Metview.h>

class NavyAccess : public MvService
{
public:
    NavyAccess(char* name) :
        MvService(name){};
    void serve(MvRequest&, MvRequest&);

    // Build the path
    bool buildPath(MvRequest&, string&);

    // Build filenames
    bool build_filenames(MvRequest&, vector<string>&);
    bool build_filenames_model(MvRequest&, vector<string>&);
    bool build_filenames_observation(MvRequest&, vector<string>&);
    bool build_filenames_ww3(MvRequest&, vector<string>&, vector<string>&, vector<string>&);
    bool build_filenames_icon(vector<string>&, vector<string>&, vector<string>&,
                              vector<string>&, vector<string>&, string&, vector<string>&);
    bool build_filenames_icon_single(vector<string>&, vector<string>&, vector<string>&,
                                     vector<string>&, vector<string>&, vector<string>&);

    // Build output data file
    string grib_icon_model(string&, vector<string>&);
    string grib_icon_single_model(string&, vector<string>&);

    // Add messages to a grib file
    bool add_messages(string& ffile, vector<string>& vparams,
                    vector<long>& vsteps, vector<long>& vlevels,
                    string& outname);

    // Split string according to a single caracter
    vector<string> split(const string&, char);

};

#endif
