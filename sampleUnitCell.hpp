#include <cstdint>

extern "C" int32_t sample_rate;				//controls the binning of counts
extern "C" const int32_t *img;  //the image, readonly
extern "C" const uint8_t *mask; //array holding our mask, good == 1, readonly
extern "C" int32_t img_width, img_height;
extern "C" int32_t *uc_sum;					//the average of all complete unitcells, write access
extern "C" int32_t uc_width, uc_height;
extern "C" double a1x,a1y,a2x,a2y; //basis vectors, a1 runs along width of unitcell
extern "C" double offset_XS,offset_YS; //a Lattice Point at a1=0 and a2=0 
extern "C" int32_t num_threads;
extern "C" int32_t *view;
extern "C" int32_t view_width, view_height;
extern "C" int32_t uc_sample_rate;
extern "C" double view_zoom;

/*public functions*/



//set Parameters and does a run
extern "C" int32_t sampleUnitCell(							
						const int32_t *frame,//1D array
						uint8_t *mk,//1D array
						int32_t w,//img_width
						int32_t h,//img_height
						int32_t *ucell,//1D array
						int32_t ucw,//uc_width
						int32_t uch,//uc_height
						double v1x, double v1y,//a1x,a1y 
						double v2x, double v2y,//a2x,a2y
						double orx, double ory,//offset_XS,offset_YS 
						int32_t srate //sample rate
);

//set Parameters for corrected periodic view of unitcell and does arun
extern "C" int32_t viewUnitCell(							
						int32_t *vv,//1D array
						int32_t vw,//view_width
						int32_t vh,//view_height
						int32_t vs,//sample rate for viewer
						double vz
						
);


//Extracts a statistical Moment of the undelying Histogramm in each pixel
//Always refers to the last run of sampling
extern "C" int32_t viewMoment(							
						int32_t order,
						float *vv//1D array						
);




// does a single run and returns number of complete unitcells
//should be checked for >=1
extern "C" int32_t transsum( void );

//re - creates the undistorted periodic view with the last seetings
extern "C" int32_t viewuc( void );



