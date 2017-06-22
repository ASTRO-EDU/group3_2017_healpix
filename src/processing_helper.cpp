#include <iostream>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include "EvalHealpix.h"


//using namespace std;


void convolution(const cv::Mat & input, cv::Mat & output, std::vector<float> & kernel, int kernel_side){
	//=====================EXERCISE 4==========================
	/*for (int i = 0; i < input.rows; i++) {
		for (int j = 0; j < input.cols; j++) {
			((uchar*)(output.data))[output.step*i + j] = ((uchar*)(input.data))[input.step*i + j];

		}
	}*/
	int k = (kernel_side - 1) / 2;
	//std::cout << "k=" << k << std::endl;
	for (int i = k; i < input.rows - k; i++) {
		for (int j = k; j < input.cols - k; j++) {
			float temp = 0;
			for (int m = -k; m <= k; m++) {
				for (int n = -k; n <= k; n++) {
					temp += (((uchar*)(input.data))[input.step*(i - m) + j - n])*kernel[3 * (m + k) + n + k];
				}
			}
			((uchar*)(output.data))[output.step*i + j] = temp;
		}
	}
	//apply the convolutional operation on each pixel of input using the passed kernel, save the result in output

}