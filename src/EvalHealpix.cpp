/***************************************************************************

                             -------------------
    begin                : Thu
    copyright            : (C) 2016 by Diamanti Alessio, Baietta Alessia, Magenta Letizia
    email                : alessio.diama@gmail.com
/***************************************************************************/

/***************************************************************************

 *                                                                         *
 *   This program is free software for non commercial purpose              *
 *   and for public research institutes; you can redistribute it and/or    *
 *   modify it under the terms of the GNU General Public License.          *
 *   For commercial purpose see appropriate license terms                  *
 *                                                                         *
 ***************************************************************************/



#include <math.h>

#include <FitsUtils.h>
#include <Selection.h>
#include <Eval.h>
#include <PilParams.h>
#include <iostream>

#include <math.h>
#include "EvalHealpix.h"
#include <healpix_map.h>
#include <healpix_base.h>
#include <fitshandle.h>
#include <healpix_map_fitsio.h>
#include <arr.h>


#include <opencv/cv.h>
#include <opencv/highgui.h>

#define DEBUG 1
#define MAX_NEIGHBOORS 20
using namespace std;


int EvalCountsHealpix(const char *outfile, double mres, double tmin,
               double tmax, double mdim, double la, double ba,
               double lonpole, double emin, double emax, double fovradmax,
               double fovradmin, double albrad, int phasecode, int filtercode,
               const char *selectionFilename,  const char *templateFilename,
               Intervals &intervals, vector< vector<int> > &counts, bool saveMaps)
{
    int status = 0;

    long mxdim = long(mdim / mres + 0.1); // dimension (in pixels) of the map
#ifdef DEBUG
    cout << "mdim: " << mdim << " mres: " << mres << " mxdim: " << mxdim << endl;
#endif
    long npixels = mxdim * mxdim;
    counts.resize(intervals.Count());

    //crea l'array di pixels
    for (int i=0; i < intervals.Count(); i++) {
        counts[i].resize(npixels);
        for (int j=0; j < npixels; j++)
            counts[i][j] = 0;
    }

    int hdutype;
    fitsfile* selectionFits;
    if (fits_open_file(&selectionFits, selectionFilename, READONLY, &status)) {
        cerr << "ERROR opening selection file " << selectionFilename << endl;
        return -1;
    }
    fits_movabs_hdu(selectionFits, 2, &hdutype, &status);

    fitsfile* templateFits;
    if (fits_open_file(&templateFits, templateFilename, READWRITE, &status)) {
        cerr << "ERROR opening template file " << templateFilename << endl;
        return -1;
    }
    fits_movabs_hdu(templateFits, 2, &hdutype, &status);
    long oldnrows;
    fits_get_num_rows(templateFits, &oldnrows, &status);
    fits_delete_rows(templateFits, 1, oldnrows, &status);

#ifdef DEBUG
    cout << "Evaluating counts.." << endl;
#endif

    int totalCounts = 0;
    cout<<"mres= "<<mres<<endl;
    int Nside = pow(2,mres);
    Healpix_Map<int> map((int)mres,NEST); // NEST is chosen for seek of efficency
 	for(int i=0;i<map.Npix();i++){ //initialize the healpix map to all zeros
		map[i]=0;
	}
    for (int intvIndex = 0; intvIndex < intervals.Count(); intvIndex++) {
#ifdef DEBUG
        cout << "Interval #" << intvIndex << endl;
#endif
        Intervals sIntv;
        sIntv.Add(intervals[intvIndex]);
        string selExpr = selection::TimesExprString(sIntv);
#ifdef DEBUG
        cout << "selExpr: " << selExpr << endl;
#endif
        fits_select_rows(selectionFits, templateFits, (char*)selExpr.c_str(), &status);
#ifdef DEBUG
      //  cout << "Rows from " << tempFilename << " selected" << endl;
#endif
        long nrows;
        fits_get_num_rows(templateFits, &nrows, &status);
#ifdef DEBUG
        cout << "Reading all " << nrows << " rows" << endl;
#endif
        int raColumn, decColumn;
        fits_get_colnum(templateFits, 1, (char*)"RA", &raColumn, &status);
        fits_get_colnum(templateFits, 1, (char*)"DEC", &decColumn, &status);

        double ra, dec, l, b, the, x, y, i = 0, ii = 0;
        double baa = ba * DEG2RAD;
        double laa = la * DEG2RAD;
        for (long k = 0; k<nrows; k++) {
            fits_read_col(templateFits, TDOUBLE, raColumn, k+1, 1, 1, NULL, &ra, NULL, &status);
            fits_read_col(templateFits, TDOUBLE, decColumn, k+1, 1, 1, NULL, &dec, NULL, &status);
            Euler(ra, dec, &l, &b, 1);
            l *= DEG2RAD;
            b *= DEG2RAD;

			pointing point = pointing((M_PI/2)-b,l); // Encodes an angular position on unitary sphere as colatitude and longitude
			int index = map.ang2pix(point);
			map[index]++;

        }
        if (nrows > 0)
            fits_delete_rows(templateFits, 1, nrows, &status);
    }
#ifdef DEBUG
    cout << "Ending evaluation" << endl;
#endif


	if(saveMaps){
		fitshandle handle = fitshandle() ;
		handle.create("./healpix_map.FITS");
		write_Healpix_map_to_fits(handle,map,PLANCK_INT64);

		//fitshandle handleC = fitshandle() ;
		//handleC.create("./convolved_&_thresholded_healpix_map.FITS");

	//	write_Healpix_map_to_fits(handleC,convolved_map,PLANCK_INT64);
	//	handleC.set_key("RADECSYS",string("FK5"),"World Coordinate System ");
	//	handleC.set_key("TCTYP2",string("GLON-HPX"),"X coordinate type ");
	//	handleC.set_key("TCTYP3",string("GLAT-HPX"),"Y coordinate type ");
		fits_close_file(selectionFits, &status);
		fits_close_file(templateFits, &status);
    }

    return status;
}