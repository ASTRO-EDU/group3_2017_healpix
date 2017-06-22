#include <FitsUtils.h>
#include <Selection.h>
#include <Eval.h>
#include <PilParams.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>

int EvalCountsHealpix(const char *outfile, double mres, double tmin,
               double tmax, double mdim, double la, double ba,
               double lonpole, double emin, double emax, double fovradmax,
               double fovradmin, double albrad, int phasecode, int filtercode,
               const char *selectionFilename,  const char *templateFilename,
               Intervals &intervals, std::vector< std::vector<int> > &counts,
               bool saveMaps);

void convolution(const cv::Mat & input, cv::Mat & output, std::vector<float> & kernel, int kernel_side);