/***************************** LICENSE START ***********************************

 Copyright 2012 ECMWF and INPE. This software is distributed under the terms
 of the Apache License version 2.0. In applying this license, ECMWF does not
 waive the privileges and immunities granted to it by virtue of its status as
 an Intergovernmental Organization or submit itself to any jurisdiction.

 ***************************** LICENSE END *************************************/

#include "NavyAccess.h"

#include <iostream>
#include <sstream>

#define NAVY_ACCESS_PATH_PREFIX_DEFAULT "/Users/fernandoii/metview/Navy/"

void NavyAccess::serve(MvRequest& in, MvRequest& out)
{
    cout << "Navy Access::serve in..." << endl;
    in.print();

    int i;

    // Build the path
    string spath;
    if (!buildPath(in, spath))
        return;

    // Build the filename(s)
    vector<string> vfiles;
    if (!build_filenames(in, vfiles))
        return;

    // Handling Observation data
    string stype = (const char*)in("TYPE");
    if (stype == "OBSERVATION") {
        // Check if file exists
        string ffile = spath + vfiles[0];
        FILE* f  = fopen(ffile.c_str(), "r");
        if (!f) {
            string error = "NavyAccess-> FILE NOT FOUND: ";
            error += ffile.c_str();
            setError(1, error.c_str());
            return;
        }
        fclose(f);

        // Create output request
        out.setVerb("BUFR");
        out("PATH") = ffile.c_str();
    }

    // Handling model data
    else if (stype == "MODEL") {
        string outname;
        string smodel = (const char*)in("MODEL");
        if (smodel == "WW3_ICON" || smodel == "WW3_COSMO" || smodel == "WW3_GFS") {
            // FAMII20200812: at the moment, returns only one netcdf file
            // Check if file exists
            string ffile = spath + vfiles[0];
            FILE* f  = fopen(ffile.c_str(), "r");
            if (!f) {
                string error = "NavyAccess-> FILE NOT FOUND: ";
                error += ffile.c_str();
                setError(1, error.c_str());
                return;
            }
            fclose(f);

            // Create output request
            out.setVerb("NETCDF");
            out("PATH") = ffile.c_str();
        }
        else {
            if (smodel == "ICON")
                outname =  this->grib_icon_model(spath, vfiles);
            else if (smodel == "ICON_SINGLE")
                outname = this->grib_icon_single_model(spath, vfiles);
            else {
                string error = "NavyAccess-> ERROR: MODEL GEONETCAST not implemented yet";
                setError(1, error.c_str());
                return;
            }

            // Create output request
            out.setVerb("GRIB");
            out("PATH") = outname.c_str();
        }
    }

    // Handling GeonetCast data
    else if (stype == "GEONETCAST") {
        string error = "NavyAccess-> ERROR: MODEL GEONETCAST not implemented yet";
        setError(1, error.c_str());
        return;
    }

    cout << "Navy Access::serve out..." << endl;
    out.print();

    return;
}

bool NavyAccess::buildPath(MvRequest& in, string& spath)
{
    // 1. Fixed part of the path. Try to get it from the environment variable:
    // METVIEW_NAVY_ACCESS_PATH. If not defined, use a default value.
    const char* fpath = getenv("METVIEW_NAVY_ACCESS_PATH");
    spath = fpath ? fpath : NAVY_ACCESS_PATH_PREFIX_DEFAULT;

    // 2. Part of the path related to the TYPE
    string stype = (const char*)in("TYPE");
    spath = spath + stype + "/";
    if ( stype == "MODEL" ) {
        string smodel = (const char*)in("MODEL");
        spath = spath + smodel + "/";
    }
    else if ( stype == "OBSERVATION" ) {
        string sobs = (const char*)in("OBSTYPE");
        spath = spath + sobs + "/";
    }
    else if ( stype == "GEONETCAST" ) {
        //spath = spath + ;
    }

    return true;
}

bool NavyAccess::build_filenames(MvRequest& in, vector<string>& vfn)
{
    // Part of the filename related to the model
    string stype = (const char*)in("TYPE");
    if (stype == "MODEL") {
        if (!build_filenames_model(in, vfn))
            return false;
    }
    else if (stype == "GEONETCAST") {
        string error = "NavyAccess-> ERROR: GEONETCAST not implemented yet ";
        setError(1, error.c_str());
        return false;
    }
    else if (stype == "OBSERVATION") {
        if (!build_filenames_observation(in, vfn))
            return false;
    }
    else {
        string error = "NavyAccess-> ERROR: TYPE not recognized ";
        setError(1, error.c_str());
        return false;
    }

for (int i = 0; i < vfn.size(); i++)
    cout << "FILENAMES: " << vfn[i] << endl;

    return true;
}

bool NavyAccess::build_filenames_model(MvRequest& in, vector<string>& vfn)
{
    int i, j, k, l, m;

    // Get all values related to parameter DATE
    int ndates = in.countValues("DATE");
    vector<string> vdates;
    vdates.reserve(ndates);
    for (i = 0; i < ndates; ++i)
        vdates.push_back((const char*)in("DATE",i));

    // Get all values related to parameter TIME
    // OBS: this version only allowed on time value, but
    // the algorithm below accepts many values
    int ntimes = in.countValues("TIME");
    vector<string> vtimes;
    vtimes.reserve(ntimes);
    std::stringstream ss;
    int iaux;
    for (i = 0; i < ntimes; ++i) {
        iaux = (int)in("TIME",i) / 100;  //hhmm
        ss << std::setw(2) << std::setfill('0') << iaux;
        vtimes.push_back(ss.str());
        ss.str("");
    }

    // Handle model WW3
    string smodel = (const char*)in("MODEL");
    if (smodel == "WW3_ICON" || smodel == "WW3_COSMO" || smodel == "WW3_GFS")
        return build_filenames_ww3(in, vdates, vtimes, vfn);

    // Get STEP values
    int nsteps = in.countValues("STEP");
    vector<string> vsteps;
    vsteps.reserve(nsteps);
    for (i = 0; i < nsteps; ++i) {
        iaux = (int)in("STEP",i);
        ss << std::setw(3) << std::setfill('0') << iaux;
        vsteps.push_back(ss.str());
        ss.str("");
    }

    // Get LEVEL_TYPE value
    string slevel_type = (const char*)in("LEVEL_TYPE");
    string sltype;
    if (slevel_type == "PRESSURE_LEVEL")
        sltype = "pressure-level";
    else if (slevel_type == "SURFACE" || slevel_type == "10_METER")
        sltype = "single-level";

    // Get LEVEL_LIST values
    vector<string> vlevels;
    int nlevels = 0;
    if ( (const char*)in("LEVEL_LIST") ) {
        nlevels = in.countValues("LEVEL_LIST");
        vlevels.reserve(nlevels);
        for (i = 0; i < nlevels; ++i)
            vlevels.push_back((const char*)in("LEVEL_LIST",i));
    }

    // Get PARAMETER values
    int npars = in.countValues("PARAMETER");
    vector<string> vparams;
    vparams.reserve(npars);
    for (i = 0; i < npars; ++i)
        vparams.push_back((const char*)in("PARAMETER",i));

    // Build filenames
    string slevel, sstep, stime, sdate, sdate_dir, sparam;
    if (smodel == "ICON") {
        if (!build_filenames_icon(vdates, vtimes, vsteps, vlevels, vparams,
                                  sltype, vfn))
            return false;

/*
        string fn_prefix = "icon_global_icosahedral_" + sltype + "_";
        string fn_suffix = "_regulargrid.grib2";
        for (i = 0; i < ndates; i++) {
            sdate_dir = vdates[i] + "/";
            sdate = fn_prefix + vdates[i];
            for (j = 0; j < ntimes; j++) {
                stime = sdate + vtimes[j] + "_";
                for (k = 0; k < nsteps; k++) {
                    sstep = stime + vsteps[k] + "_";
                    if (nlevels) {
                        for (l = 0; l < nlevels; l++) {
                            slevel = sstep + vlevels[l] + "_";
                            for (m = 0; m < npars; m++) {
                                vfn.push_back(sdate_dir + slevel + vparams[m] + fn_suffix);
                            }
                        }
                    }
                    else {
                        for (m = 0; m < npars; m++)
                            vfn.push_back(sdate_dir + sstep + vparams[m] + fn_suffix);
                    }
                }
            }
        }*/
        
    }
    else if (smodel == "ICON_SINGLE") {  // TEST TO READ FILE ICON TOTAL: icon13km_00.grib2
        if (!build_filenames_icon_single(vdates, vtimes, vsteps, vlevels, vparams, vfn))
                return false;


        /*        string fn_prefix = "icon13km_";
        string fn_suffix = ".grib2";
        for (i = 0; i < ndates; i++) {
            sdate_dir = vdates[i] + "/";
            for (j = 0; j < ntimes; j++) {
                for (k = 0; k < nsteps; k++) {
                    sstep = "," + vsteps[k];
                    if (nlevels) {
                        for (l = 0; l < nlevels; l++) {
                            slevel = "," + vlevels[l];
                            for (m = 0; m < npars; m++) {
                                sparam = "," + vparams[m];
                                vfn.push_back(sdate_dir + fn_prefix + vtimes[j] + fn_suffix + sstep + sparam + slevel);
                            }
                        }
                    }
                    else {
                        for (m = 0; m < npars; m++) {
                            sparam = "," + vparams[m];
                            vfn.push_back(sdate_dir + fn_prefix + vtimes[j] + fn_suffix + sstep + sparam);
                        }
                    }
                }
            }
        }*/
    }
    else {
        string error = "NavyAccess-> ERROR: MODEL not recognized ";
        setError(1, error.c_str());
        return false;
    }

    return true;
}

bool NavyAccess::build_filenames_observation(MvRequest& in, vector<string>& vfn)
{
    int i, j;

    // Get all values related to parameters DATE/TIME
    int ndates = in.countValues("DATE");
    vector<string> vdates;
    vdates.reserve(ndates);
    for (i = 0; i < ndates; ++i)
        vdates.push_back((const char*)in("DATE",i));

    int ntimes = in.countValues("TIME");
    vector<string> vtimes;
    vtimes.reserve(ntimes);
    std::stringstream ss;
    int iaux;
    for (i = 0; i < ntimes; ++i) {
        iaux = (int)in("TIME",i) / 100;  //hhmm
        ss << std::setw(2) << std::setfill('0') << iaux;
        vtimes.push_back(ss.str());
        ss.str("");
    }

    // Get "UNKNOWN" values
    string ssat = (const char*)in("UNKNOWN1");
    string sunk = (const char*)in("UNKNOWN2");

    // Build filenames
    string stime, sdate, sfn;
    string sobs = (const char*)in("OBSTYPE");
    if (sobs == "ASCAT") {
        string fn_prefix = "ascat_";
        string fn_suffix = ".l2_bufr";
        for (i = 0; i < ndates; i++) {
            sdate = fn_prefix + vdates[i];
        //for (j = 0; j < ntimes; j++) {
        //    stime = sdate + vtimes[j];
            if (ssat == "metopb")
                sfn = sdate + "_123900_" + ssat + "_40885_eps_o_" + sunk + fn_suffix;
            else
                sfn = sdate + "_123900_" + ssat + "_09043_eps_o_" + sunk + fn_suffix;

            vfn.push_back(sfn);
        }
    }
    else {
        string error = "NavyAccess-> ERROR: OBSTYPE not recognized ";
        setError(1, error.c_str());
        return false;
    }

    return true;
}

bool NavyAccess::build_filenames_ww3(MvRequest& in, vector<string>& vdates, vector<string>& vtimes, vector<string>& vfn)
{
    // Build filename prefix
    string fn_prefix;
    string smodel = (const char*)in("MODEL");
    if (smodel == "WW3_ICON")
        fn_prefix = "ww3ico13_met_";
    else if (smodel == "WW3_COSMO")
        fn_prefix = "ww3icosmo_met_";
    else if (smodel == "WW3_GFS")
        fn_prefix= "ww3gfs_met_";

    // Build filenames
    string stime, sdate;
    string fn_suffix = ".nc";
    for (int i = 0; i < vdates.size(); i++) {
        sdate = fn_prefix + vdates[i];
        for (int j = 0; j < vtimes.size(); j++) {
            stime = sdate + vtimes[j];
            vfn.push_back(stime + fn_suffix);
        }
    }

    return true;
}

bool NavyAccess::build_filenames_icon(vector<string>& vdates, vector<string>& vtimes,
                                      vector<string>& vsteps, vector<string>& vlevels,
                                      vector<string>& vparams, string& sltype,
                                      vector<string>& vfn)
{
    string slevel, sstep, stime, sdate, sdate_dir;
    string fn_prefix = "icon_global_icosahedral_" + sltype + "_";
    string fn_suffix = "_regulargrid.grib2";
    for (int i = 0; i < vdates.size(); i++) {
        sdate_dir = vdates[i] + "/";
        sdate = fn_prefix + vdates[i];
        for (int j = 0; j < vtimes.size(); j++) {
            stime = sdate + vtimes[j] + "_";
            for (int k = 0; k < vsteps.size(); k++) {
                sstep = stime + vsteps[k] + "_";
                if (vlevels.size()) {
                    for (int l = 0; l < vlevels.size(); l++) {
                        slevel = sstep + vlevels[l] + "_";
                        for (int m = 0; m < vparams.size(); m++)
                            vfn.push_back(sdate_dir + slevel + vparams[m] + fn_suffix);
                    }
                }
                else {
                    for (int m = 0; m < vparams.size(); m++)
                        vfn.push_back(sdate_dir + sstep + vparams[m] + fn_suffix);
                }
            }
        }
    }

    return true;
}

bool NavyAccess::build_filenames_icon_single(vector<string>& vdates, vector<string>& vtimes,
                                      vector<string>& vsteps, vector<string>& vlevels,
                                      vector<string>& vparams, vector<string>& vfn)
{
    string slevel, sstep, sdate_dir, sparam;
    string fn_prefix = "icon13km_";
    string fn_suffix = ".grib2";
    for (int i = 0; i < vdates.size(); i++) {
        sdate_dir = vdates[i] + "/";
        for (int j = 0; j < vtimes.size(); j++) {
            for (int k = 0; k < vsteps.size(); k++) {
                sstep = "," + vsteps[k];
                if (vlevels.size()) {
                    for (int l = 0; l < vlevels.size(); l++) {
                        slevel = "," + vlevels[l];
                        for (int m = 0; m < vparams.size(); m++) {
                            sparam = "," + vparams[m];
                            vfn.push_back(sdate_dir + fn_prefix + vtimes[j] + fn_suffix + sstep + sparam + slevel);
                        }
                    }
                }
                else {
                    for (int m = 0; m < vparams.size(); m++) {
                        sparam = "," + vparams[m];
                        vfn.push_back(sdate_dir + fn_prefix + vtimes[j] + fn_suffix + sstep + sparam);
                    }
                }
            }
        }
    }

    return true;
}

string NavyAccess::grib_icon_model(string& spath, vector<string>& vfiles)
{
    // Build the requested fieldset
    // Open each input grib file and add them to the output GRIB file

    // Auxilliary variables for GribApi
    char shortName[20];
    size_t len;
    int error       = 0;
    grib_handle* h  = NULL;
    grib_context* c = grib_context_get_default();

    // Create the output file name
    string outname(marstmp());
    for (int i = 0; i < (signed)vfiles.size(); ++i) {
        // Full filename
        string ffile = spath + vfiles[i];

        // Open the input file
        FILE* f = fopen(ffile.c_str(), "r");
        if (!f) {
            string error = "NavyAccess-> FILE NOT FOUND: ";
            error += ffile.c_str();
            setError(1, error.c_str());
            return std::string();
        }

        // Loop on all GRIB messages in file
        long level;
        while ((h = grib_handle_new_from_file(c, f, &error)) != NULL) {
            grib_write_message(h, outname.c_str(), "a");
            grib_handle_delete(h);
        }
        fclose(f);
    }

    return outname;
}

string NavyAccess::grib_icon_single_model(string& spath, vector<string>& vfiles)
{
    // Build the requested fieldset
    // For each input file filter the requested fields by PARAM/LEVEL/STEP
    // and add them to the output GRIB file

    // Create the output file name
    string sdate = "";
    vector<string> vparams;
    vector<long> vlevels;
    vector<long> vsteps;
    string outname(marstmp());
    vector<string> vf;
    for (int i = 0; i < (signed)vfiles.size(); ++i) {
        // Split the string to retrieve the filename, step, parameter and level
        vf = split(vfiles[i], ',');

        // Check if this is the same file as the previous one
        vector<string> fn = split(vf[0], '/');
        if (fn[0] == sdate) {
            vsteps.push_back(stol(vf[1]) * 60); // convert to minutes
            vparams.push_back(vf[2]);
            if (vf.size() > 3)
               vlevels.push_back(stol(vf[3]));

            continue;
        }
        
        // Add requested messages to the output file
        if ( i != 0 ) {
            string ffile = spath + vf[0];
            if (!this->add_messages(ffile, vparams, vsteps, vlevels, outname))
                return std::string();
        }

        // Update variables
        vsteps.clear();
        vparams.clear();
        vlevels.clear();
        vsteps.push_back(stoi(vf[1]) * 60); // convert to minutes
        vparams.push_back(vf[2]);
        if (vf.size() > 3)
           vlevels.push_back(stoi(vf[3]));
    }

    // Check if needs to save the last set of messages
    if (vparams.size())
    {
        string ffile = spath + vf[0];
        if (!this->add_messages(ffile, vparams, vsteps, vlevels, outname))
            return std::string();
    }

    return outname;
}

bool NavyAccess::add_messages(string& ffile, vector<string>& vparams,
                                 vector<long>& vsteps, vector<long>& vlevels,
                                 string& outname)
{
    // Auxilliary variables for GribApi
    char shortName[20];
    size_t len;
    int error       = 0;
    grib_handle* h  = NULL;
    grib_context* c = grib_context_get_default();

    // Open input grib file
    FILE* f = fopen(ffile.c_str(), "r");
    if (!f) {
        string error = "NavyAccess-> FILE NOT FOUND: ";
        error += ffile.c_str();
        setError(1, error.c_str());
        return false;
    }

    // Loop on all GRIB messages in file
    long ilevel;
    long istep;
    string sparam;
    while ((h = grib_handle_new_from_file(c, f, &error)) != NULL) {
        len = 20;
        grib_get_string(h, "shortName", shortName, &len);
        sparam = string(shortName);
        grib_get_long(h, "level", &ilevel);
        grib_get_long(h, "endStep", &istep);
        if(std::find(vparams.begin(), vparams.end(), sparam) != vparams.end())
            if(std::find(vsteps.begin(), vsteps.end(), istep) != vsteps.end())
                if (vlevels.empty() || (std::find(vlevels.begin(), vlevels.end(), ilevel) != vlevels.end()))
                    grib_write_message(h, outname.c_str(), "a");

        grib_handle_delete(h);
    }
    return true;
}

vector<string> NavyAccess::split(const string &s, char delim)
{
    vector<string> result;
    std::stringstream ss(s);
    std::string item;

    while (std::getline (ss, item, delim))
        result.push_back (item);

    return result;
}

int main(int argc, char** argv)
{
    MvApplication theApp(argc, argv);
    NavyAccess data((char*)"NAVY_ACCESS");

    theApp.run();
}
