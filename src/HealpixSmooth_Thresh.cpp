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


int SmoothAndThresh( double mres, bool saveMaps, float thresh)
{
   int status = 0;
	/*Create a grey-scale image from the counting just done*/
	Healpix_Map<int> map((int)mres,NEST); // NEST is chosen for seek of efficency
	read_Healpix_map_from_fits ("./healpix_map.FITS", map);
	int max = -1;
	long int nPix = map.Npix();
	for(int i=0;i<nPix;i++){
		if(map[i]>= max){
			max = (int)map[i];
		}
	}
#ifdef DEBUG
	cout<<"max = "<< max<<endl;
#endif

 float convolved_data[nPix];
 bool smooth = true;
	if(smooth){


	 float data[nPix];
		for(int i=0;i<nPix;i++){
			data[i] = ((float)map[i]*255.0)/(float)max;
		}
		for(int i=0;i<nPix;i++){
				convolved_data[i]=0;
		}
		float kernel_side = 19.0;

		float kernel[19][19]{{0.000002,	0.000006,	0.000014,	0.000028,	0.000051,	0.000084,	0.000124,	0.000163,	0.000193,	0.000203,	0.000193,	0.000163,	0.000124,	0.000084,	0.000051,	0.000028,	0.000014,	0.000006,	0.000002},
								{0.000006,	0.000015,	0.000035,	0.000071,	0.000131,	0.000215,	0.000316,	0.000416,	0.000491,	0.000519,	0.000491,	0.000416,	0.000316,	0.000215,	0.000131,	0.000071,	0.000035,	0.000015,	0.000006},
								{0.000014,	0.000035,	0.00008,	0.000163,	0.000299,	0.000491,	0.000722,	0.00095,	0.001121,	0.001184,	0.001121,	0.00095,	0.000722,	0.000491,	0.000299,	0.000163,	0.00008,	0.000035,	0.000014},
								{0.000028,	0.000071,	0.000163,	0.000334,	0.000612,	0.001004,	0.001476,	0.001944,	0.002293,	0.002423,	0.002293,	0.001944,	0.001476,	0.001004,	0.000612,	0.000334,	0.000163,	0.000071,	0.000028},
								{0.000051,	0.000131,	0.000299,	0.000612,	0.001121,	0.00184,	0.002705,	0.003562,	0.004201,	0.004439,	0.004201,	0.003562,	0.002705,	0.00184,	0.001121,	0.000612,	0.000299,	0.000131,	0.000051},
								{0.000084,	0.000215,	0.000491,	0.001004,	0.00184,	0.00302,	0.004439,	0.005845,	0.006895,	0.007285,	0.006895,	0.005845,	0.004439,	0.00302,	0.00184,	0.001004,	0.000491,	0.000215,	0.000084},
								{0.000124,	0.000316,	0.000722,	0.001476,	0.002705,	0.004439,	0.006525,	0.008593,	0.010136,	0.010709,	0.010136,	0.008593,	0.006525,	0.004439,	0.002705,	0.001476,	0.000722,	0.000316,	0.000124},
								{0.000163,	0.000416,	0.00095,	0.001944,	0.003562,	0.005845,	0.008593,	0.011315,	0.013347,	0.014102,	0.013347,	0.011315,	0.008593,	0.005845,	0.003562,	0.001944,	0.00095,	0.000416,	0.000163},
								{0.000193,	0.000491,	0.001121,	0.002293,	0.004201,	0.006895,	0.010136,	0.013347,	0.015743,	0.016634,	0.015743,	0.013347,	0.010136,	0.006895,	0.004201,	0.002293,	0.0011201,	0.000491,	0.000193},
								{0.000203,	0.000519,	0.001184,	0.002423,	0.004439,	0.007285,	0.010709,	0.014102,	0.016634,	0.017575,	0.016634,	0.014102,	0.010709,	0.007285,	0.004439,	0.002423,	0.001184,	0.000519,	0.000203},
								{0.000193,	0.000491,	0.001121,	0.002293,	0.004201,	0.006895,	0.010136,	0.013347,	0.015743,	0.016634,	0.015743,	0.013347,	0.010136,	0.006895,	0.004201,	0.002293,	0.001121,	0.000491,	0.000193},
								{0.000163,	0.000416,	0.00095,	0.001944,	0.003562,	0.005845,	0.008593,	0.011315,	0.013347,	0.014102,	0.013347,	0.011315,	0.008593,	0.005845,	0.003562,	0.001944,	0.00095,	0.000416,	0.000163},
								{0.000124,	0.000316,	0.000722,	0.001476,	0.002705,	0.004439,	0.006525,	0.008593,	0.010136,	0.010709,	0.010136,	0.008593,	0.006525,	0.004439,	0.002705,	0.001476,	0.000722,	0.000316,	0.000124},
								{0.000084,	0.000215,	0.000491,	0.001004,	0.00184,	0.00302,	0.004439,	0.005845,	0.006895,	0.007285,	0.006895,	0.005845,	0.004439,	0.00302,	0.00184,	0.001004,	0.000491,	0.000215,	0.000084},
								{0.000051,	0.000131,	0.000299,	0.000612,	0.001121,	0.00184,	0.002705,	0.003562,	0.004201,	0.004439,	0.004201,	0.003562,	0.002705,	0.00184,	0.001121,	0.000612,	0.000299,	0.000131,	0.000051},
								{0.000028,	0.000071,	0.000163,	0.000334,	0.000612,	0.001004,	0.001476,	0.001944,	0.002293,	0.002423,	0.002293,	0.001944,	0.001476,	0.001004,	0.000612,	0.000334,	0.000163,	0.000071,	0.000028},
								{0.000014,	0.000035,	0.00008,	0.000163,	0.000299,	0.000491,	0.000722,	0.00095,	0.001121,	0.001184,	0.001121,	0.00095,	0.000722,	0.000491,	0.000299,	0.000163,	0.00008,	0.000035,	0.000014},
								{0.000006,	0.000015,	0.000035,	0.000071,	0.000131,	0.000215,	0.000316,	0.000416,	0.000491,	0.000519,	0.000491,	0.000416,	0.000316,	0.000215,	0.000131,	0.000071,	0.000035,	0.000015,	0.000006},
								{0.0000020,	0.000006,	0.000014,	0.000028,	0.000051,	0.000084,	0.000124,	0.000163,	0.000193,	0.000203,	0.000193,	0.000163,	0.000124,	0.000084,	0.000051,	0.000028,	0.000014,	0.000006,	0.000002}
								};//sigma = 3
		/*float kernel[11][11] = {{0.000002,	0.00001,	0.000047,	0.000136,	0.000259,	0.00032,	0.000259,	0.000136,	0.000047,	0.00001,	0.000002},
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
								}; // kernel for sigma=1.5 e kernel_side=11*/

		/*float kernel[13][13]={{0.000006,	0.000022,	0.000067,	0.000158,	0.000291,	0.000421,	0.000476,	0.000421,	0.000291,	0.000158,	0.000067,	0.000022,	0.000006},
							{0.000022,	0.000086,	0.000258,	0.000608,	0.001121,	0.001618,	0.001829,	0.001618,	0.001121,	0.000608,	0.000258,	0.000086,	0.000022},
							{0.000067,	0.000258,	0.000777,	0.00183,	0.003375,	0.004873,	0.005508,	0.004873,	0.003375,	0.00183,	0.000777,	0.000258,	0.000067},
							{0.000158,	0.000608,	0.00183,	0.004312,	0.007953,	0.011483,	0.012978,	0.011483,	0.007953,	0.004312,	0.00183,	0.000608,	0.000158},
							{0.000291,	0.001121,	0.003375,	0.007953,	0.014669,	0.021179,	0.023938,	0.021179,	0.014669,	0.007953,	0.003375,	0.001121,	0.000291},
							{0.000421,	0.001618,	0.004873,	0.011483,	0.021179,	0.030579,	0.034561,	0.030579,	0.021179,	0.011483,	0.004873,	0.001618,	0.000421},
							{0.000476,	0.001829,	0.005508,	0.012978,	0.023938,	0.034561,	0.039062,	0.034561,	0.023938,	0.012978,	0.005508,	0.001829,	0.000476},
							{0.000421,	0.001618,	0.004873,	0.011483,	0.021179,	0.030579,	0.034561,	0.030579,	0.021179,	0.011483,	0.004873,	0.001618,	0.000421},
							{0.000291,	0.001121,	0.003375,	0.007953,	0.014669,	0.021179,	0.023938,	0.021179,	0.014669,	0.007953,	0.003375,	0.001121,	0.000291},
							{0.000158,	0.000608,	0.00183,	0.004312,	0.007953,	0.011483,	0.012978,	0.011483,	0.007953,	0.004312,	0.00183,	0.000608,	0.000158},
							{0.000067,	0.000258,	0.000777,	0.00183,	0.003375,	0.004873,	0.005508,	0.004873,	0.003375,	0.00183,	0.000777,	0.000258,	0.000067},
							{0.000022,	0.000086,	0.000258,	0.000608,	0.001121,	0.001618,	0.001829,	0.001618,	0.001121,	0.000608,	0.000258,	0.000086,	0.000022},
							{0.000006,	0.000022,	0.000067,	0.000158,	0.000291,	0.000421,	0.000476,	0.000421,	0.000291,	0.000158,	0.000067,	0.000022,	0.000006}
							};//sigma=2*/
		float kernel_reordered_Array[MAX_NEIGHBOORS][8*MAX_NEIGHBOORS];
		int lim = (kernel_side/2.0)-1; //counts how many neighborhoods must be extracted

		/*Reorder the kernel according to neighborhoods representation */
		int centre = kernel_side/2.0;
		float k_1 = kernel[centre][centre]; // central element
		for(int times=0;times<centre;times++){
			int i=0;
			int count=0;
			int cursor_row = centre+(1*times)+2;
			int cursor_column = centre-(times+1);
			while(count<3+(2*times)){ // left side
				cursor_row--;
				kernel_reordered_Array[times][i]= kernel[cursor_row][cursor_column];
				i++;
				count++;
			}
			count=0;
			while(count<2+(2*times)){ //up side
				cursor_column++;
				count++;
				kernel_reordered_Array[times][i]= kernel[cursor_row][cursor_column];
				i++;
			}
			count=0;
			while(count<2+(2*times)){ // right side
				cursor_row++;
				count++;
				kernel_reordered_Array[times][i]= kernel[cursor_row][cursor_column];
				i++;
			}
			count=0;
			while(count<1+(2*times)){ // down side
				cursor_column--;
				count++;
				kernel_reordered_Array[times][i]= kernel[cursor_row][cursor_column];
				i++;
			}

		}


		fix_arr<int,81> arr[MAX_NEIGHBOORS]; // array for holding the (i+1)-neighborood
		fix_arr<int,8> tempA; //holds temporary neighbors
		float temp; // accumulator for convolution sums
		for(int i=0;i<nPix;i++){
			temp=0;
			temp +=data[i]*k_1; // product of central elements

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
						}else{//if arr[times][s] does not exit also his neighbor does not exist
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
						}else{//if arr[times][s] does not exit also his neighbor does not exist
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
						}else{//if arr[times][s] does not exit also his neighbor does not exist
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
						}else{//if arr[times][s] does not exit also his neighbor does not exist
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
						}else{//if arr[times][s] does not exit also his neighbor does not exist
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
						}else{//if arr[times][s] does not exit also his neighbor does not exist
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
						}else{//if arr[times][s] does not exit also his neighbor does not exist
							arr[times+1][s+7]=-1;
						}
					}
				}
				z=0;
			}

			/*Finally calculates the convolution sum */
			for(int q=0;q<lim+1;q++){
				for(int k=0;k<8*(q+1);k++){
					if(arr[q][k]!=-1){
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
		histogram[i]=0;
	}
	for(int i=0;i<nPix;i++){
		histogram[(int)convolved_data[i]]++;
	}
	float mass_prob;
	float percent_acc = 0;
	int index=-1;
	float treshold = -1;
	 /* calcolate the treshold as 99% of all pixels*/
	for(int i=0;i<256;i++){
		mass_prob= (float)histogram[i]/(float)nPix;

		percent_acc+=mass_prob;
		if(percent_acc>=thresh){
			treshold= i;
			break;
		}
	}

#ifdef DEBUG
	    cout << "Ending evaluation" << endl;
#endif

	for(int i=0;i<nPix;i++){ /* Suppress the pixel how's intensities are under the threshold  */
		if(convolved_data[i]<=treshold){
			convolved_data[i]=0;
		}
	}
	Healpix_Map<float> convolved_map((int)mres,NEST); /* copies the data into the Healpix map */
	for(int i=0;i<nPix;i++){
		convolved_map[i]= convolved_data[i] ;
	}

	if(saveMaps){
		fitshandle handleC = fitshandle() ;
		handleC.create("./convolved_&_thresholded_healpix_map.FITS");
		write_Healpix_map_to_fits(handleC,convolved_map,PLANCK_INT64);
		handleC.set_key("RADECSYS",string("FK5"),"World Coordinate System ");
		handleC.set_key("TCTYP2",string("GLON-HPX"),"X coordinate type ");
		handleC.set_key("TCTYP3",string("GLAT-HPX"),"Y coordinate type ");

    }

    return status;
}