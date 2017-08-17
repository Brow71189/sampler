#include <cstdlib>
#include <math.h>
#include <climits>
#include <string.h>
#include <algorithm>
#include <pthread.h>
#include <unistd.h>
#include <assert.h>
#include "sampleUnitCell.hpp"

/*Macros*/

#define UPDATE_IDET (idet = 1.0 / (  a1x *  a2y -  a1y *  a2x ))
#define xd_a(a1,a2) ( a1x*(a1) + a2x*(a2) )
#define yd_a(a1,a2) ( a1y*(a1) + a2y*(a2) )
#define a1d(x,y) ( ( a2y*(x) - a2x*(y) ) * idet )
#define a2d(x,y) ( ( a1x*(y) - a1y*(x) ) * idet )

/*global variables*/
const int32_t *img = nullptr;
const uint8_t *mask = nullptr;
int32_t img_width = 2048, img_height=2048;
int32_t *uc_sum = nullptr;
int32_t uc_width = 32, uc_height = 32;
int32_t *view = nullptr;
int32_t view_width = 256, view_height = 256;
 
int32_t sample_rate= 8;
int32_t uc_sample_rate = 8;
double view_zoom = 1.0;
double scale=1.0;

double a1x=32.0,a1y=0.0,a2x=0.0,a2y=32.0;
double offset_XS = 1024.0;
double offset_YS = 1024.0;
int32_t num_threads = 1;//actual use
pthread_t threads[8];//max number of of threads
int32_t ids[8];

int32_t  **unitcells = nullptr;	
bool *dirty = nullptr;			
double idet = 1.0;
bool first_call = true;
int32_t uc_gridX; 
int32_t uc_gridY; 
int32_t uc_rangeX;
int32_t uc_rangeY;
int32_t uc_len;
int32_t uc_area;
int32_t img_area;
int32_t view_area;
int32_t area;
int32_t task = 1; //1..sample, 2..re-sample
int32_t good_cells = 1;


/*private Prototypes*/

void init_sampling(void);
void set_uc_spacings(void);
int32_t ind00_a1_a2(double a1, double a2);
int32_t ind0_a1_a2(double a1, double a2);
int32_t ind1_a1_a2(double a1len, double a2len);
void resample(int32_t start, int32_t end);
int32_t stretch(int32_t start, int32_t end);
void init_xorshift1024star(void);
double rand_d(void);
void* run(void*);

/*implementations*/

extern "C" int32_t sampleUnitCell(	
						const int32_t *frame,
						uint8_t *mk,
						int32_t w,
						int32_t h,
						int32_t *ucell,
						int32_t ucw,
						int32_t uch,
						double v1x, double v1y, 
						double v2x, double v2y,
						double orx, double ory, 
						int32_t srate)
{
	img = frame;
	mask = mk;
	img_width = w;
	img_height = h;
	uc_sum = ucell;
	uc_width = ucw;
	uc_height = uch;
	a1x = v1x; a1y = v1y;
	a2x = v2x; a2y = v2y;
	offset_XS = orx;
	offset_YS = ory;	
	sample_rate = srate;
	
	uc_area = uc_width*uc_height;
	img_area = img_width * img_height;
	
	
	return transsum();
}

extern "C" int32_t viewUnitCell(							
						int32_t *vv,//1D array
						int32_t vw,//view_width
						int32_t vh,//view_height
						int32_t vs,//sampling for view
						double vz//zoom for view
)
{
	view = vv;
	view_width = vw;
	view_height = vh;
	view_area = view_width * view_height;
	uc_sample_rate = vs;
	view_zoom = vz;
	return viewuc();
}

extern "C" int32_t viewMoment(	int32_t order, float *vv)
{
	if(first_call)
	{	return -1;}
	const double inv_order = 1.0/order;
	for(int32_t uc_pix=0; uc_pix < uc_area; ++uc_pix)
	{
		int64_t sum = 0.0;
		for(int32_t uc_ind = 0; uc_ind < uc_len; ++uc_ind)
		{
			if(unitcells[uc_ind] != nullptr)
			{	sum+=unitcells[uc_ind][uc_pix];}
		}
		const double avg = ((double)sum/good_cells);
		switch(order)
		{
			case 0:
				vv[uc_pix]=(float)sum;
			break;
			case 1:
				vv[uc_pix]=(float)avg;	
			break;
			default:
			{	
				double sumO = 0.0;
				for(int32_t uc_ind = 0; uc_ind < uc_len; ++uc_ind)
				{
					if(unitcells[uc_ind] != nullptr)
					{	sumO += pow(unitcells[uc_ind][uc_pix]-avg,order);}
				}
				vv[uc_pix]=(float)pow(sumO,inv_order);
			}
			break;
		}
	}
	return good_cells;
}


void set_uc_spacings(void)
{
	uc_gridX =  (int32_t)(0.5*( fabs(a1x) + fabs(a2x) ));
	uc_gridY =  (int32_t)(0.5*( fabs(a1y) + fabs(a2y) ));
	if(uc_gridX == 0) {	uc_gridX = 1;}
	if(uc_gridY == 0) { uc_gridY = 1;}
	uc_rangeX = img_width / uc_gridX + 1; //We have to account for incomplete  uc at the edges
	uc_rangeY = img_height / uc_gridY + 1; //in both directions
	uc_len = uc_rangeX * uc_rangeY;		
}



void init_sampling(void)
{
	
	
	if( unitcells != nullptr )
	{
		int32_t **uc(unitcells);
		for(int32_t raw_ind = 0; raw_ind < uc_len ; ++raw_ind)
		{	delete[] *(uc++);}
		delete[] unitcells;
		unitcells = nullptr;
		delete[] dirty;
		dirty = nullptr;	
	}
	
	set_uc_spacings();
	dirty = new bool[uc_len];
	unitcells = new int32_t*[uc_len];
	std::fill_n(dirty,uc_len,false);
	std::fill_n(unitcells,uc_len,nullptr);
	return;
}

int32_t ind00_a1_a2(double a1, double a2)
{
	const int32_t x( (int32_t) ( floor(xd_a(a1, a2) + offset_XS) ) );
	const int32_t y( (int32_t) ( floor(yd_a(a1, a2) + offset_YS) ) );
	return ( ( (x<0) || (y<0) || (x>=img_width) || (y>=img_height) ) ? -1 : ( (x/uc_gridX) + (y/uc_gridY) * uc_rangeX ) );
}

int32_t ind0_a1_a2(double a1, double a2)
{
	//define the four corners
	const double a1l( floor(a1) );
	const double a1h( a1l + 1.0 );
	const double a2l( floor(a2) );
	const double a2h( a2l + 1.0 ); 
	//check if any of them lies outside the frame
	return ( ( ind00_a1_a2(a1l, a2h) < 0 ) ||
			 ( ind00_a1_a2(a1h, a2l) < 0 ) ||
			 ( ind00_a1_a2(a1h, a2h) < 0 ) ) 
			? -1 : ind00_a1_a2(a1l, a2l);
}

int32_t ind1_a1_a2(double a1len, double a2len)
{
	const double fa1( a1len - floor(a1len) );
	const double fa2( a2len - floor(a2len) );
	const int32_t ia1( ( (int32_t)floor(fa1 * uc_width) + uc_width ) % uc_width );
	const int32_t ia2( ( (int32_t)floor(fa2 * uc_height) + uc_height ) % uc_height );
	return (ia1 + uc_width * ia2);
}


unsigned long states[16];
const double range_H( 1.0/((double)ULONG_MAX) ); //prceision loss of ~ 10bit approx 1000 identical values in a period
int32_t p_xor;

void init_xorshift1024star(void)
{		
		for(int32_t i = 0; i < 16; ++i)
		{	//we only use the lower 16 bit of standard rand();
			unsigned long p0 = ( (unsigned long)(rand() & 0x0000FFFF) );
			unsigned long p1 = ( (unsigned long)(rand() & 0x0000FFFF) ) << 16;
			unsigned long p2 = ( (unsigned long)(rand() & 0x0000FFFF) ) << 32;
			unsigned long p3 = ( (unsigned long)(rand() & 0x0000FFFF) ) << 48;
			
			states[i] = ( p0 | p1 | p2 | p3 );
		}
		while( 0 != (p_xor = rand()%16) ){};
}

double rand_d(void)
{
	long s0( states[ p_xor ] );
	long s1( states[ p_xor = ( (p_xor+1) & 15) ] );
	s1 ^= (s1 << 31);
	s1 ^= (s1 >> 11);
	s0 ^= (s0 >> 30);
	return ( (double)  ( ( states[p_xor] = (s0 ^ s1) ) * (unsigned long)(1181783497276652981)  ) ) * range_H;
}



extern "C" int32_t transsum(void)
{
	if(img==nullptr)
	{	return -1;}
	UPDATE_IDET;
	init_sampling();
	if(first_call)
	{
		first_call = false;
		init_xorshift1024star();
	}
	task = 1;// 1 means sampling
	area = img_area; 
	for(int32_t i = 1; i < num_threads; ++i )
	{
		ids[i]=i;
		pthread_create(&threads[i],nullptr,run,&ids[i]);	
	}
	const int32_t start = 0;
	const int32_t end = img_area/num_threads;
	resample(start,end);
	for(int32_t done = 1;done < num_threads; ++done)
	{
		pthread_join(threads[done],nullptr);
	}
	good_cells = 0;
	if(uc_sum != nullptr)
	{
		memset(uc_sum,0,uc_area*sizeof(int32_t));
		int32_t **uc( unitcells );
		for(int32_t ind(0); ind < uc_len ; ++ind)
		{
			if( *uc != nullptr ) 
			{
				++good_cells;
				int32_t *s( uc_sum );
				int32_t *a( *uc );
				
				for(int32_t ind0(0); ind0 < uc_area; ++ind0)
				{	
					const int32_t val( *(a++) );
					*(s++) += val;		
				}	
			}
			++uc;
		}
	}
	return good_cells; 
}

void resample(const int32_t start, const int32_t end)
{
	const int32_t *fr(img+start);
	for(int32_t raw_ind(start); raw_ind < end ; ++raw_ind)
	{
		const int32_t val( *(fr++) );
		const int32_t dv0( 1 + (val-sample_rate/2)/sample_rate);
		const int32_t x( raw_ind % img_width);
		const int32_t y( raw_ind / img_width);
		const int32_t flg( (mask!=nullptr) ? mask[raw_ind] : 1 );
		const double cx( (double)x - offset_XS );
		const double cy( (double)y - offset_YS );
		//ignore flg at initial test so the dirty markings can affect the tentative unitcell
		for(int32_t v(0); v < val; v+=dv0) //very clean sampling where a pixel may even be shared among adjacent unitcells
		{
			int32_t dv( ( (v+dv0) > val ) ? (val-v) : dv0  );
			if(dv==0) {	break;}
			const double subx( rand_d() );
			const double suby( rand_d() ); 
			
			{
				const double  a1len(  a1d(cx+subx,cy+suby) );
				const double  a2len(  a2d(cx+subx,cy+suby) );
				const int32_t uc_key(  ind0_a1_a2( a1len,  a2len) );
				
				if( (uc_key != -1) && (!dirty[uc_key]) ) 
				//the corrsponding unitcell is enitrely inside the frame and not yet known to be dirty
				{
					//while(busy[uc_key]);//,std::memory_order_relaxed
					//busy[uc_key] = true;
					if( (unitcells[uc_key] == nullptr) && (flg==1) ) //in an interesting area
					{
						unitcells[uc_key] = new int32_t[uc_area];
						memset(unitcells[uc_key], 0, uc_area*sizeof(int32_t));
					}
					 
					if( !(flg==1) )  //anything we are not interested in
					{
						dirty[uc_key] = true;
						if(unitcells[uc_key] != nullptr) //could spare this check
						{
							delete[] unitcells[uc_key];
							unitcells[uc_key] = nullptr;
						}	 
					} 
					if( !dirty[uc_key] )
					{	
						const int32_t loc_ind( ind1_a1_a2( a1len, a2len) );
						unitcells[uc_key][loc_ind] += dv;
					}
					//busy[uc_key]=false;//,std::memory_order_relaxed
				}
			}
			if(!(flg==1))//anything else than graphene
			{	break;} //just mark the place as dirty and move on
		}	
	}
}

int32_t viewuc(void)
{
	
	if(first_call)
	{	return -1;}
	scale = 1.0/view_zoom;
	task = 2;// 2 means streching
	area = img_area; 
	for(int32_t i = 1; i < num_threads; ++i )
	{
		ids[i]=i;
		pthread_create(&threads[i],nullptr,run,&ids[i]);	
	}
	const int32_t start = 0;
	const int32_t end = view_area/num_threads;
	
	stretch(start,end);
	for(int32_t done = 1;done < num_threads; ++done)
	{
		pthread_join(threads[done],nullptr);
	}
	return view_area;
}




void* run(void* id)
{
	const int32_t myid = *((int32_t*)(id)); 
	assert(myid >= 0);
	assert(myid <= num_threads);
	
	
	const int32_t start = myid*area/num_threads;
	const int32_t end = (myid+1)*area/num_threads;
	switch(task)
	{
		case 1: resample(start,end);
		break;
		case 2: stretch(start,end);
		break;	
	}
	pthread_exit(id);	
}


int32_t stretch(int32_t start, int32_t end)
{
	
	int32_t *fr(view+start);
	for(int32_t raw_ind = start; raw_ind < end ; ++raw_ind)
	{
		const double x( raw_ind % view_width);
		const double y( raw_ind / view_width);
		int32_t sum = 0;
		
		for(int32_t i=0; i < uc_sample_rate; ++i)
		{
			const double subx( rand_d() );
			const double suby( rand_d() ); 
			const double  a1len(  scale * a1d(x + subx,y + suby) );
			const double  a2len(  scale * a2d(x + subx,y + suby) );
			const int32_t loc_ind( ind1_a1_a2( a1len, a2len) );	
			sum += uc_sum[loc_ind];
		}
		*(fr++) = sum/uc_sample_rate;			
	}
	return view_area;
}
