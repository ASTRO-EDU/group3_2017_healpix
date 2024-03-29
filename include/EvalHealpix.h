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

int SmoothAndThresh(double mres, bool saveMaps, float thresh);