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

    /*struct prjprm *prj = new(prjprm);
    prjini(prj);
    ProjectionType proj;
    if (!LoadProjection(projection, proj)) {
        cerr << "Error loading projection '" << projection << "'" << endl;
        return false;
    }*/

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
    //int mapArray[786432];
 	for(int i=0;i<12*pow(Nside,2);i++){ //initialize the healpix map to all zeros
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
        //cout << "Rows from " << tempFilename << " selected" << endl;
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
				//std::cout<<"l= "<<l<<std::endl;
				//	std::cout<<"b= "<<b<<std::endl;

			//std::cout<<"map_order= "<<map.Order()<<std::endl;
		//	std::cout<<"nPix= "<<map.Npix()<<std::endl;

			/*for(int i=0;i<13;i++){
				std::cout<<"map["<<i<<"](prev)=\t"<<map[i]<<std::endl;
			}*/
			pointing point = pointing((M_PI/2)-b,l); // Encodes an angular position on unitary sphere as colatitude and longitude
			int index = map.ang2pix(point);
			//std::cout<<"i: \t"<<index<<std::endl;
			map[index]++;
			//mapArray[index]++;
		/*	for(int i=0;i<12;i++){
				std::cout<<"map["<<i<<"](next)=\t"<<map[i]<<std::endl;
			}*/
            //abbiamo (l,b) in coordinate galattiche
            //da qui va a determinare quale è il pixel a cui l'evento gamma appartiene (e.g. il fotone)

          /*  the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
            if (the < -1.0)
                the = M_PI;
            else if (the > 1.0)
                the = 0.0;
            else
                the = acos(the);

            if (proj == ARC) {
                x = RAD2DEG/Sinaa(the) * cos(b)*sin(l-laa);
                y = RAD2DEG/Sinaa(the) * (sin(b)*cos(baa) - cos(b)*sin(baa)*cos(l-laa));

                i = (int)floor(((-x+(mdim/2.))/mres));
                ii = (int)floor(((y+(mdim/2.))/mres));
            }
            else if (proj == AIT) {
                l=l-laa;

                if (l < M_PI)
                l=-l;
                else
                l=2*M_PI -l;

                x=RAD2DEG*(sqrt(2.0)*2.0*cos(b)*sin(l/2.0))/sqrt(1.0 + cos(b)*cos(l/2.0) );
                y=RAD2DEG*(sqrt(2.0)*sin(b))/sqrt(1.0 + cos(b)*cos(l/2.0) );

                i=(int)floor(((x+(mdim/2.))/mres));
                ii=(int)floor(((y+(mdim/2.))/mres));
            }*/

           /* if (inmap(i, ii, mxdim)) {
                counts[intvIndex][ii*mxdim+i]++;
                totalCounts++;
            }*/
        }
        if (nrows > 0)
            fits_delete_rows(templateFits, 1, nrows, &status);
    }
#ifdef DEBUG
    cout << "Ending evaluation" << endl;
#endif
	/*Create a grey-scale image from the counting just done*/
	int max = -1;
	int nPix = map.Npix();
	for(int i=0;i<nPix;i++){
		if(map[i]>= max){
			max = (int)map[i];
		}
	}
#ifdef DEBUG
	cout<<"max = "<< max<<endl;
#endif
	float data[nPix];
	for(int i=0;i<nPix;i++){
		data[i] = ((float)map[i]*255.0)/(float)max;
	}
	float convolved_data[nPix];
	Healpix_Map<float> convolved_map((int)mres,NEST);
	float kernel_side = 9.0;

	float k_1= 0.150342;
	/*float kernel[8] = { 0.077847,	0.123317,	0.077847,
				0.123317,	0.123317,
				0.077847,	0.123317,	0.077847};
				*/
	/*float kernel_reordered[8]={	0.059912,0.094907,0.059912,0.094907,0.059912,0.094907,0.059912,	0.094907};
	float kernel_reordered_2[16]={0.023528,0.033969,0.038393,0.033969,0.023528,0.033969,0.038393,0.033969,0.023528,0.033969,0.038393,0.033969,0.023528,0.033969,0.038393,0.033969};
	*/


	/*sigma = 1.5 side=9*/
	float kernel_reordered[8]={0.044672,0.055338,0.044672,0.055338,0.044672,0.055338,0.044672,0.055338};
	float kernel_reordered_2[16]={0.012358,	0.023496,0.029106,0.023496,0.012358,0.023496,0.029106,0.023496,0.012358,0.023496,0.029106,0.023496,0.012358,0.023496,0.029106,0.023496};
	float kernel_reordered_3[24]={0.00145,0.004233,0.008048,0.00997,0.008048,0.004233,0.00145,0.004233,0.008048,0.00997,0.008048,0.004233,0.00145,0.004233,0.008048,0.00997,0.008048,0.004233,0.00145,0.004233,0.008048,0.00997,0.008048,0.004233};
	float kernel_reordered_4[32]={0.000072,0.000323,0.000944,0.001794,0.002222,0.001794,0.000944,0.000323,0.000072,0.000323,0.000944,0.001794,0.002222,0.001794,0.000944,0.000323,0.,.000323,0.000944,0.001794,0.002222,0.001794,0.000944,0.000323,0.000072,0.000323,0.000944,0.001794,0.002222,0.001794,0.000944,0.000323};
	/*float kernel_reordered2[24]={};*/


	fix_arr<int,81> arr[10];
	fix_arr<int,8> tempA;
	int temp;
	for(int i=0;i<nPix;i++){
		temp=0;
		temp +=map[i]*k_1;
		map.neighbors(i,tempA);//extract the 8-neighborood of the i-th pixel
		for(int j=0;j<8;j++){
			arr[0][j]=tempA[j];
		}
		int z = 0;
		int separator;
		int lim = (kernel_side/2.0)-1;
		//cout<<"lim="<<lim<<endl;
		for(int times=0;times<lim;times++){ //if necessary extract the other neighbors
			//cout<<"times="<<times<<endl;
			separator = 2*(times+1);
			//cout<<"separator="<<separator<<endl;
			for(int s=0;s<8+(8*times);s++){
				if(s==0){
					if(arr[times][s]!= -1){ //bottom-left corner
						map.neighbors(arr[times][s],tempA); // 8-neighborhood of arr[s]}
						arr[times+1][s+z]= tempA[z];
						z++;
						arr[times+1][s+z]= tempA[z];
						arr[times+1][s+(15+(8*times))]= tempA[7];
					}else{
						arr[times+1][s+z]= -1;
						z++;
						arr[times+1][s+z]= -1;
						arr[times+1][s+(15+(8*times))]= -1;
					}
				}else if((s!=0) && (s<separator)){//left side
					if(arr[times][s]!=0){
						map.neighbors(arr[times][s],tempA); // 8-neighborood of arr[s]
						arr[times+1][s+1]=tempA[z];
					}else{
						arr[times+1][s+1]=-1;
					}
				}else if(s==separator){ // top-left corner
					if(arr[times][s]!=0){
						map.neighbors(arr[times][s],tempA); // 8-neighborood of arr[s]
						arr[times+1][s+z]= tempA[z];
						z++;
						arr[times+1][s+z]= tempA[z];
						z++;
						arr[times+1][s+z]= tempA[z];
					}else{
						arr[times+1][s+z]= -1;
						z++;
						arr[times+1][s+z]= -1;
						z++;
						arr[times+1][s+z]= -1;
					}
				}else if((s>separator)&& (s<separator*2)){//up-side
					if(arr[times][s]!=0){
						map.neighbors(arr[times][s],tempA); // 8-neighborood of arr[s]
						arr[times+1][s+3]=tempA[z];
					}else{
						arr[times+1][s+(times+1)]=-1;
					}
				}else if(s==separator*2){ // top-right corner
					if(arr[times][s]!=0){
						map.neighbors(arr[times][s],tempA); // 8-neighborood of arr[s]
						arr[times+1][s+z]= tempA[z];
						z++;
						arr[times+1][s+z]= tempA[z];
						z++;
						arr[times+1][s+z]= tempA[z];
					}else{
						arr[times+1][s+z]= -1;
						z++;
						arr[times+1][s+z]= -1;
						z++;
						arr[times+1][s+z]= -1;
					}
				}else if((s>separator*2)&& (s<separator*3)){// right side
					if(arr[times][s]!=0){
						map.neighbors(arr[times][s],tempA); // 8-neighborood of arr[s]
						arr[times+1][s+5]=tempA[z];
					}else{
						arr[times+1][s+5]=-1;
					}
				}else if(s==separator*3){ //bottom-right corner
					if(arr[times][s]!=0){
						map.neighbors(arr[times][s],tempA); // 8-neighborood of arr[s]
						arr[times+1][s+z]= tempA[z];
						z++;
						arr[times+1][s+z]= tempA[z];
						z++;//z=7
						arr[times+1][s+z]= tempA[z];
					}else{
						map.neighbors(arr[times][s],tempA); // 8-neighborood of arr[s]
						arr[times+1][s+z]= -1;
						z++;
						arr[times+1][s+z]= -1;
						z++;//z=7
						arr[times+1][s+z]= -1;
					}
				}else if((s>separator*3)&& (s<separator*4)){ // bottom side
					if(arr[times][s]!=0){
						map.neighbors(arr[times][s],tempA); // 8-neighborood of arr[s]
						arr[times+1][s+7]=tempA[z];
					}else{
						arr[times+1][s+7]=-1;
					}
				}
			}
			z=0;
		}

		for(int k=0;k<8;k++){
			if(arr[0][k]!=-1){
				temp += data[arr[0][k]]*kernel_reordered[(k+4)%8];
			}
		}
		for(int k=0;k<16;k++){
			if(arr[1][k]!=-1){
				temp += data[arr[1][k]]*kernel_reordered_2[(k+8)%16];
			}
		}

		for(int k=0;k<24;k++){
			if(arr[2][k]!=-1){
				temp += data[arr[2][k]]*kernel_reordered_3[(k+12)%24];
			}

		}
		for(int k=0;k<32;k++){
			if(arr[3][k]!=-1){
				temp += data[arr[3][k]]*kernel_reordered_4[(k+16)%32];
			}
		}
		convolved_data[i] = temp;
	}
	for(int i=0;i<nPix;i++){
		convolved_map[i]= convolved_data[i] ;
	}




	fitshandle handle = fitshandle() ;
	handle.create("./healpix_map.FITS");
	write_Healpix_map_to_fits(handle,map,PLANCK_INT64);

	fitshandle handleC = fitshandle() ;
	handleC.create("./convolved_healpix_map.FITS");
	write_Healpix_map_to_fits(handleC,convolved_map,PLANCK_INT64);
   /* if (saveMaps) {
        vector<unsigned short> sum;
        sum.resize(npixels);
        for (unsigned int i=0; i<mxdim; i++)
            for (unsigned int j=0; j<mxdim; j++)
                sum[i*mxdim+j] = 0;
        for (int intvIndex=0; intvIndex<intervals.Count(); intvIndex++)
            for (unsigned int i=0; i<mxdim; i++)
                for (unsigned int j=0; j<mxdim; j++)
                    sum[i*mxdim+j] += counts[intvIndex][i*mxdim+j];

        int bitpix = USHORT_IMG;
        long naxis = 2;
        long naxes[2] = { mxdim, mxdim };
        fitsfile * mapFits;
        cout << "Creating file " << outfile << endl;
        if (fits_create_file(&mapFits, outfile, &status) != 0) {
            cerr << "Error opening file " << outfile << endl;
            return status;
        }
        cout << "creating Counts Map...................................." << endl;
        fits_create_img(mapFits, bitpix, naxis, naxes, &status);
        cout << "writing Counts Map with " << totalCounts << " events" << endl;
        fits_write_img(mapFits, bitpix, 1, npixels, &sum[0], &status);
        cout << "writing header........................................" << endl << endl;
        fits_update_key(mapFits, TDOUBLE, "CRVAL1", &la, NULL, &status);
        fits_update_key(mapFits, TDOUBLE, "CRVAL2", &ba, NULL, &status);
        char projstr1[FLEN_FILENAME];
        char projstr2[FLEN_FILENAME];
        if (proj == AIT) {
            strcpy(projstr1,"GLON-AIT");
            strcpy(projstr2,"GLAT-AIT");
        } else {
            strcpy(projstr1,"GLON-ARC");
            strcpy(projstr2,"GLAT-ARC");
        }
        fits_update_key(mapFits, TSTRING, "CTYPE1", projstr1, NULL, &status);
        fits_update_key(mapFits, TSTRING, "CTYPE2", projstr2, NULL, &status);
        double xx = mdim/mres/2+0.5;
        fits_update_key(mapFits, TDOUBLE, "CRPIX1", &xx, NULL, &status);
        fits_update_key(mapFits, TDOUBLE, "CRPIX2", &xx, NULL, &status);
        xx = -mres;
        fits_update_key(mapFits, TDOUBLE, "CDELT1", &xx, NULL, &status);
        xx = mres;
        fits_update_key(mapFits, TDOUBLE, "CDELT2", &xx, NULL, &status);
        char unit[] = "deg";
        fits_update_key(mapFits, TSTRING, "CUNIT1", unit, NULL, &status);
        fits_update_key(mapFits, TSTRING, "CUNIT2", unit, NULL, &status);
        char str3[] = "FK5";
        fits_update_key(mapFits, TSTRING,  "RADESYS", str3, NULL, &status);
        xx = 2000.0;
        fits_update_key(mapFits, TDOUBLE,  "EQUINOX", &xx, NULL, &status);
        fits_update_key(mapFits, TDOUBLE,  "LONPOLE", &lonpole, NULL, &status);
        WriteTime(mapFits, tmin, tmax);
        fits_update_key(mapFits, TDOUBLE,  "MINENG", &emin, NULL, &status);
        fits_update_key(mapFits, TDOUBLE,  "MAXENG", &emax, NULL, &status);
        char str1[] = "AGILE";
        fits_update_key(mapFits, TSTRING,  "TELESCOP", str1, NULL, &status);
        char str2[] = "GRID";
        fits_update_key(mapFits, TSTRING,  "INSTRUME", str2, NULL, &status);
        char str6[] = "T";
        fits_update_key(mapFits, TSTRING,  "PIXCENT", str6, NULL, &status);
        char str7[FLEN_FILENAME] = "";
        fits_update_key(mapFits, TSTRING,  "BUNIT", str7, NULL, &status);
        fits_update_key(mapFits, TDOUBLE,  "FOV", &fovradmax, "Maximum off-axis angle (deg)", &status);
        fits_update_key(mapFits, TDOUBLE,  "FOVMIN", &fovradmin, "Minimum off-axis angle (deg)", &status);
        fits_update_key(mapFits, TDOUBLE,  "ALBEDO", &albrad, "Earth zenith angle (deg)", &status);
        fits_update_key(mapFits, TINT,  "PHASECOD", &phasecode, "Orbital phase code", &status);
        fits_update_key(mapFits, TINT,  "FILTERCO", &filtercode, "Event filter code", &status);
        cout << "Closing " << outfile << endl;
        fits_close_file(mapFits, &status);
    }*/
    fits_close_file(selectionFits, &status);
    fits_close_file(templateFits, &status);

    return status;
}

