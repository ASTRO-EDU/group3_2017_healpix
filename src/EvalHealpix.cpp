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
#define MAX_NEIGHBOORS 10
using namespace std;


int EvalCountsHealpix(const char *outfile, double mres, double tmin,
               double tmax, double mdim, double la, double ba,
               double lonpole, double emin, double emax, double fovradmax,
               double fovradmin, double albrad, int phasecode, int filtercode,
               const char *selectionFilename,  const char *templateFilename,
               Intervals &intervals, vector< vector<int> > &counts, bool saveMaps,bool smooth, bool binarize)
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
			//std::cout<<"i: \t"<<index<<std::endl;
			map[index]++;

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
Healpix_Map<float> convolved_map((int)mres,NEST);
float convolved_data[nPix];
	if(smooth){
		float data[nPix];
		for(int i=0;i<nPix;i++){
			data[i] = ((float)map[i]*255.0)/(float)max;
		}
		for(int i=0;i<nPix;i++){
				convolved_data[i]=0;
		}
		cout<<"max = "<< max<<endl;
		float kernel_side = 11.0;

         /*sigma = 5*/
		/*float kernel[11][11]={{0.004411,	0.005278,	0.006068,	0.006704,	0.007117,	0.00726,	0.007117,	0.006704,	0.006068,	0.005278,	0.004411},
							{0.005278,	0.006315,	0.00726,	0.008021,	0.008515,	0.008687,	0.008515,	0.008021,	0.00726,	0.006315,	0.005278},
							{0.006068,	0.00726,	0.008347,	0.009222,	0.00979,	0.009988,	0.00979,	0.009222,	0.008347,	0.00726,	0.006068},
							{0.006704,	0.008021,	0.009222,	0.010189,	0.010817,	0.011034,	0.010817,	0.010189,	0.009222,	0.008021,	0.006704},
							{0.007117,	0.008515,	0.00979,	0.010817,	0.011483,	0.011714,	0.011483,	0.010817,	0.00979,	0.008515,	0.007117},
							{0.00726,	0.008687,	0.009988,	0.011034,	0.011714,	0.01195,	0.011714,	0.011034,	0.009988,	0.008687,	0.00726},
							{0.007117,	0.008515,	0.00979,	0.010817,	0.011483,	0.011714,	0.011483,	0.010817,	0.00979,	0.008515,	0.007117},
							{0.006704,	0.008021,	0.009222,	0.010189,	0.010817,	0.011034,	0.010817,	0.010189,	0.009222,	0.008021,	0.006704},
							{0.006068,	0.00726,	0.008347,	0.009222,	0.00979,	0.009988,	0.00979,	0.009222,	0.008347,	0.00726,	0.006068},
							{0.005278,	0.006315,	0.00726,	0.008021,	0.008515,	0.008687,	0.008515,	0.008021,	0.00726,	0.006315,	0.005278},
							{0.004411,	0.005278,	0.006068,	0.006704,	0.007117,	0.00726,	0.007117,	0.006704,	0.006068,	0.005278,	0.004411}
							};*/


		float kernel[11][11] = {{
		0.000002,	0.00001,	0.000047,	0.000136,	0.000259,	0.00032,	0.000259,	0.000136,	0.000047,	0.00001,	0.000002},
		{0.00001,	0.000072,	0.000322,	0.000939,	0.001785,	0.002212,	0.001785,	0.000939,	0.000322,	0.000072,	0.00001},
		{.000047,	0.000322,	0.001443,	0.004212,	0.008008,	0.009921,	0.008008,	0.004212,	0.001443,	0.000322,	0.000047},
		{0.000136,	0.000939,	0.004212,	0.012297,	0.02338,	0.028963,	0.02338,	0.012297,	0.004212,	0.000939,	0.000136},
		{0.000259,	0.001785,	0.008008,	0.02338,	0.044453,	0.055067,	0.044453,	0.02338,	0.008008,	0.001785,	0.000259},
		{0.00032,	0.002212,	0.009921,	0.028963,	0.055067,	0.068216,	0.055067,	0.028963,	0.009921,	0.002212,	0.00032},
		{0.000259,	0.001785,	0.008008,	0.02338,	0.044453,	0.055067,	0.044453,	0.02338,	0.008008,	0.001785,	0.000259},
		{0.000136,	0.000939,	0.004212,	0.012297,	0.02338,	0.028963,	0.02338,	0.012297,	0.004212,	0.000939,	0.000136},
		{0.000047,	0.000322,	0.001443,	0.004212,	0.008008,	0.009921,	0.008008,	0.004212,	0.001443,	0.000322,	0.000047},
		{0.00001,	0.000072,	0.000322,	0.000939,	0.001785,	0.002212,	0.001785,	0.000939,	0.000322,	0.000072,	0.00001},
		{0.000002,	0.00001,	0.000047,	0.000136,	0.000259,	0.00032,	0.000259,	0.000136,	0.000047,	0.00001,	0.000002}
		}; // kernel for sigma=1.5 e kernel_side=11
		float kernel_reordered_Array[MAX_NEIGHBOORS][8*MAX_NEIGHBOORS];
		int lim = (kernel_side/2.0)-1; //counts how many neighborhoods must be extracted
		/*Reorder the kernel according to neighborhoods representation */
		int centre = kernel_side/2.0;
		float k_1 = kernel[centre][centre];
		for(int times=0;times<lim;times++){
			int i=0;
			int count=0;
			int cursor_row = centre+(1*times);
			cout<<"cursor_row"<<cursor_row<<endl;

			int cursor_column = centre-(times+1);
			cout<<"cursor_column"<<cursor_column<<endl;

			while(count<3+(2*times)){
				cursor_row--;
				kernel_reordered_Array[times][i]= kernel[cursor_row][cursor_column];
				i++;
				count++;
			}
			count=0;
			while(count<2+(2*times)){
				cursor_column++;
				count++;
				kernel_reordered_Array[times][i]= kernel[cursor_row][cursor_column];
				i++;
			}
			count=0;
			while(count<2+(2*times)){
				cursor_row++;
				count++;
				kernel_reordered_Array[times][i]= kernel[cursor_row][cursor_column];
				i++;
			}
			count=0;
			while(count<1+(2*times)){
				cursor_column--;
				count++;
				kernel_reordered_Array[times][i]= kernel[cursor_row][cursor_column];
				i++;
			}

		}

		/*Gaussian kernel for sigma = 1.5 side=11*/
		/*TODO : optimal kernel size for sigma=1.5 is 11=> enlarge the side for better approssimation of Gaussian Beam*/
		/*float k_1= 	0.068552; // central value of the kernel matrix
		float kernel_reordered[8]={0.044672,0.055338,0.044672,0.055338,0.044672,0.055338,0.044672,0.055338};
		float kernel_reordered_2[16]={0.012358,	0.023496,0.029106,0.023496,0.012358,0.023496,0.029106,0.023496,0.012358,0.023496,0.029106,0.023496,0.012358,0.023496,0.029106,0.023496};
		float kernel_reordered_3[24]={0.00145,0.004233,0.008048,0.00997,0.008048,0.004233,0.00145,0.004233,0.008048,0.00997,0.008048,0.004233,0.00145,0.004233,0.008048,0.00997,0.008048,0.004233,0.00145,0.004233,0.008048,0.00997,0.008048,0.004233};
		float kernel_reordered_4[32]={0.000072,0.000323,0.000944,0.001794,0.002222,0.001794,0.000944,0.000323,0.000072,0.000323,0.000944,0.001794,0.002222,0.001794,0.000944,0.000323,0.,.000323,0.000944,0.001794,0.002222,0.001794,0.000944,0.000323,0.000072,0.000323,0.000944,0.001794,0.002222,0.001794,0.000944,0.000323};
		float kernel_reordered_5[40]={0.000002,0.00001,0.000047,0.000136,0.000259,0.00032,0.000259,0.000136,0.000047,0.00001,0.000002,0.00001,0.000047,0.000136,0.000259,0.00032,0.000259,0.000136,0.000047,0.00001,0.000002,0.00001,0.000047,0.000136,0.000259,0.00032,0.000259,0.000136,0.000047,0.00001,0.000002,0.00001,0.000047,0.000136,0.000259,0.00032,0.000259,0.000136,0.000047,0.00001};
		/*float kernel_reordered2[24]={};*/
		/*float* kernel_reordered_Array[5];
		kernel_reordered_Array[0]=kernel_reordered;
		kernel_reordered_Array[1]=kernel_reordered_2;
		kernel_reordered_Array[2]=kernel_reordered_3;
		kernel_reordered_Array[3]=kernel_reordered_4;
		kernel_reordered_Array[4]=kernel_reordered_5;*/


		fix_arr<int,81> arr[MAX_NEIGHBOORS]; // array for holding the (i+1)-neighborood
		fix_arr<int,8> tempA; //holds temporary neighbors
		int temp; // accumulator for convolution sums
		for(int i=0;i<nPix;i++){
			temp=0;
			temp +=data[i]*k_1;
			map.neighbors(i,tempA);//extract the 8-neighborood of the i-th pixel
			for(int j=0;j<8;j++){
				arr[0][j]=tempA[j];
			}
			int z = 0; //used to select neighbors
			int separator; //indexing dividing arr[times][s] to reflect sides of neighborhood matrix

			/* if necessary extract the other neighbors */
			for(int times=0;times<lim;times++){ //if kernel side is 3 there is no need to extract other neighbors
				separator = 2*(times+1);
				for(int s=0;s<8+(8*times);s++){
					if(s==0){//bottom-left corner
						if(arr[times][s]!= -1){
							map.neighbors(arr[times][s],tempA); // 8-neighborhood of arr[s]}
							arr[times+1][s+z]= tempA[z];
							z++;
							arr[times+1][s+z]= tempA[z];
							arr[times+1][s+(15+(8*times))]= tempA[7]; // last neighbor pixel
						}else{ //if arr[times][s] does not exit also his neighbor does not exist
							arr[times+1][s+z]= -1;
							z++;
							arr[times+1][s+z]= -1;
							arr[times+1][s+(15+(8*times))]= -1;
						}
					}else if((s!=0) && (s<separator)){//left side
						if(arr[times][s]!=-1){
							map.neighbors(arr[times][s],tempA); // 8-neighborood of arr[s]
							arr[times+1][s+1]=tempA[z];
						}else{
							arr[times+1][s+1]=-1;
						}
					}else if(s==separator){ // top-left corner
						if(arr[times][s]!=-1){
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
						if(arr[times][s]!=0-1){
							map.neighbors(arr[times][s],tempA); // 8-neighborood of arr[s]
							arr[times+1][s+3]=tempA[z];
						}else{
							arr[times+1][s+(times+1)]=-1;
						}
					}else if(s==separator*2){ // top-right corner
						if(arr[times][s]!=-1){
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
						if(arr[times][s]!=-1){
							map.neighbors(arr[times][s],tempA); // 8-neighborood of arr[s]
							arr[times+1][s+5]=tempA[z];
						}else{
							arr[times+1][s+5]=-1;
						}
					}else if(s==separator*3){ //bottom-right corner
						if(arr[times][s]!=-1){
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
						if(arr[times][s]!=-1){
							map.neighbors(arr[times][s],tempA); // 8-neighborood of arr[s]
							arr[times+1][s+7]=tempA[z];
						}else{
							arr[times+1][s+7]=-1;
						}
					}
				}
				z=0;
			}
			/*Finally calculates the convolution sum */
			for(int q=0;q<lim+1;q++){
				for(int k=0;k<8*(q+1);k++){
					if(arr[0][k]!=-1){
						temp += data[arr[q][k]]*kernel_reordered_Array[q][(k+4*(q+1))%(8*(q+1))];
					}
				}
			}
			convolved_data[i] = temp;
		}
	}

	/*thresholding*/

	/*first create the histogram*/
	int histogram[255];
	for(int i=0;i<256;i++){
		histogram[i]=0;;
	}
	for(int i=0;i<nPix;i++){
		histogram[(int)convolved_data[i]]++;
	}
	float mass_prob;

	float percent_acc = 0;
	int index=-1;
	float treshold = -1;
	cout<<"nPix "<<nPix<<endl;
	 /* calcolate the treshold as 99% of all pixels*/
	for(int i=0;i<256;i++){
			cout<<"histogram[i]= "<<histogram[i]<<endl;
		mass_prob= (float)histogram[i]/(float)nPix;

		percent_acc+=mass_prob;
		cout<<"percent_acc= "<<percent_acc<<endl;
		if(percent_acc>=0.998){
			treshold= i;
			break;
		}
	}
	cout<<"treshold= "<<treshold<<endl;

	for(int i=0;i<nPix;i++){
		if(convolved_data[i]<=treshold){
			convolved_data[i]=0;
		}else{
			convolved_data[i]=255;
		}
	}

	for(int i=0;i<nPix;i++){
		convolved_map[i]= convolved_data[i] ;
	}

	if(saveMaps){
		fitshandle handle = fitshandle() ;
		handle.create("./healpix_map.FITS");
		write_Healpix_map_to_fits(handle,map,PLANCK_INT64);

		fitshandle handleC = fitshandle() ;
		handleC.create("./convolved_&_thresholded_healpix_map.FITS");
		write_Healpix_map_to_fits(handleC,convolved_map,PLANCK_INT64);
		fits_close_file(selectionFits, &status);
		fits_close_file(templateFits, &status);
    }

    return status;
}