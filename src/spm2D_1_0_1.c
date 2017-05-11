/*
 *  SPM2D
 *  Version 1.0.1 - 22/11/2012
 *  Jop van Rooij
 *  j.a.vanrooij@uu.nl 
 *  Jop.VanRooij@jic.ac.uk
 */

/*
 * This code includes:
 * - Segmentation Potts Model
 * - Cell size estimation algorithm
 * - A background subtraction image filter
 * - A Gaussian blur image filter
 */

/*
 * Features:
 * - Segment multiple images sequentially
 * - Set parameters in command line
 * - Save different stages of segmentation, including intermediate segmentation results
 * - Include parameters and/or time in output directory
 * - Use filtered image, or one of the original color channels for cell seed
 */

/* Changes in 1.0.1
 * - Changed a bug that cuased the gaussian filter to potentially use the wrong target
 * - Corrected spelling mistakes
 */

//libraries
#include "excal2.h"
#include <ctype.h>

//extra data for cell structures
typedef struct extras {
	int color;
	
	//for persistence vector
	FLOAT vectorx;
	FLOAT vectory;
	FLOAT originx;
	FLOAT originy;	
	FLOAT nu;
} Extras;

//definitions
#define CELLWALL      0
#define CELL          1
#define MAXTYPES      2
#define EXTRAS ((Extras*)(cells.extras))

//cellstructures
static Cell cells;

//planes
static TYPE **state;
static TYPE **image;
static TYPE **template;
static TYPE **tempstate;
static TYPE **tempstate2;

static FLOAT **energy;
static FLOAT **ftempstate;

//pngs
static Png pngCellstateInter;
static Png pngCellstateInit;
static Png pngCellstateRaw;
static Png pngCellstateDet;
static Png pngCellstateFinal;

static Png pngCellcolors;
static Png pngTemplate;
static Png pngEnergy;
static Png pngCombined;

//strings
static char a_string[STRING];
static char displayname[STRING];
static char imagedirname[STRING];
static char outputdirname[STRING];
static char parfilename[STRING];

//internal parameters
static double mu;
static double areaconstraint;
static int kappa;
static double nu;
static int totalj;
static int npixels;
static double pixelsize;
static double coverage;
static double tcoverage;
static int totalarea;
static int totaltargetarea;
static double noiselevel;
static int currentz;
static int firstimage;
static int tijd;
static int ncells;

//parameters that will be read from parfile
static int hamiltonianneigh;
static int imageinteraction;
static int targetcoverage;
static int estimationoffset;
static int temperaturescale;
static int persistencestrength;
static int backgroundsubtraction;
static int gaussianblur;
static int runtime;

static int seedgaussian;
static int mincellarea;
static int persistencelength;
static int persistencescale;
static int maxcells;
static int commandline;

static int cellmarker;
static int seedmarker;
static int includepars;
static int includetime;
static int displayinterval;
static int infointerval;
static int saveinit;
static int saveinter;
static int saveraw;
static int savedet;
static int savefinal;
static int savecellcolors;
static int savecombined;
static int savetemplate;
static int saveenergy;

//external variables
extern double temperature;
extern double dissipation;
extern double lambda;
extern int boundary;
extern int ncol;
extern int nrow;
extern int dh;
extern int sigma;
extern int sigmaneigh;
extern int iposition;
extern int jposition;
extern int ineigh;
extern int jneigh;
extern int specialcelltype;
extern int targetarea;
extern int seed;
extern int graphics;
extern int scale;
extern int linewidth;

//setup functions
void SetUp(int ac,char *av[]);
void Allocate();
void AllocateCells();
void AllocateOutput(char *av[]);
void CloseOutput();
void MyInDat(char filename[],char *format,char *text,void *value);
void ReadPars(char name[]);
void SetOutput(char name[]);
void SetPars(char *av[]);

//SPM functions
void SetSegmentation();
void SetCellseed();
void SetTemplate();
void Extra();
void Aftermath();
void UpdateSegmentation();
void UpdatePersistence();

//IO functions
void PrintInfo();
void SaveGraphics();

//Support fucntions
void RemoveStochasticity2D(TYPE **target, FLOAT **template);
void Watershed2D(TYPE **target, TYPE **template, TYPE **tempstate);
void RemoveProtrutions2D(TYPE **target, TYPE **tempstate);
void RenameCells2D(TYPE **target);

double EstimateTargetCoverage2D(TYPE **target, TYPE **tempstate);
void FFindLocalMinima2D(TYPE **target, FLOAT **template, TYPE **tempstate);
void BackgroundSubtraction2D(TYPE **target, int neighsize, TYPE **tempstate);
void FGaussianFilter2D(FLOAT **target, int kernelsize); 
void FStretchHistogram2D(FLOAT **target);
int CountClusters2D(TYPE **target);

//color functions
int Code2Color(int state, int which);
int Color2Code(int intensity, int which);
int RGB2Code(int r,int g,int b);
int MixColors(int color1, int color2, double opacity);

int nlay = 100;


int main(int argc,char *argv[])
{
	SetUp(argc,argv);
	
	for(currentz=1;currentz<=nlay;currentz++) {  
		
		ReadSizePNG(&nrow,&ncol,argv[currentz+firstimage]);
		ReadPNG(image,1,1,argv[currentz+firstimage]);
		
		printf("\nLoading image\n");
		printf("\tReading: %s\n",argv[currentz+firstimage]);
		printf("\tImage size is %d x %d pixels\n",ncol,nrow);
		printf("\tcomplete\n");
		
		
		SetTemplate();
		SetCellseed();
		SetSegmentation();
		
		//run segmentation
		printf("\nStarting segmentation\n");
		for(tijd=0;tijd<runtime;tijd++) {
			
			UpdateSegmentation();
			if(persistencelength && !(tijd%persistencelength)) UpdatePersistence();
			
			//if(displayinterval && !(tijd%displayinterval)) DisplayGraphics();
			if(infointerval && !(tijd%infointerval)) PrintInfo();
			
			if(!hamiltonianneigh) Hamiltonian(state,&cells,&Extra,&Aftermath);
			else WideHamiltonian(state,&cells,hamiltonianneigh,&Extra,&Aftermath);
			
			if(saveinter && !(tijd%saveinter)) PlanePNG(state,&pngCellstateInter,0);
			
		}
		
		printf("Segmentation complete\n");
		
		//save output
		if(saveraw) PlanePNG(state,&pngCellstateRaw,0);
		RemoveStochasticity2D(state,energy);
		if(savedet) PlanePNG(state,&pngCellstateDet,0);
		Watershed2D(state,template,tempstate);
		RemoveProtrutions2D(state,tempstate);
		if(savefinal) PlanePNG(state,&pngCellstateFinal,0);
		SaveGraphics();
		printf("\n");
	}
	
	CloseOutput();
	
	printf("All images segmentated, exiting now\n");
	return 0;
}

//setup functions
void SetUp(int ac,char *av[]) 
{  
	int maxncol,maxnrow;
	
	printf("Seed: %d\n",seed);
	
	// check # of arguments
	if(ac<4){ 
		fprintf(stdout,"Usage: %s outputdirectory parfile image(s)\n",av[0]);
		exit(EXIT_FAILURE);
	}
	
	// read parameters
	fprintf(stdout,"Reading parameters from %s\n",av[2]);
	ReadPars(av[2]); 
	
	// check # of arguments (again)
	if(commandline && ac<13) {
		fprintf(stdout,"Usage: %s outputdirectory parfile HN II TC EO TS PS BS GB RT image(s)\n",av[0]);
		fprintf(stdout,"Usage: or turn of \'commandline\' in parfile: %s\n",av[2]);
		
		fprintf(stdout,"\n");
		fprintf(stdout,"HN -> hamiltonianneigh\n");
		fprintf(stdout,"II -> imageinteraction\n");
		fprintf(stdout,"TC -> targetcoverage\n");
		fprintf(stdout,"EO -> estimationoffset\n");
		fprintf(stdout,"TS -> temperaturescale\n");
		fprintf(stdout,"PS -> persistencestrength\n");
		fprintf(stdout,"BS -> backgroundsubtraction\n");
		fprintf(stdout,"GB -> gaussianblur\n");
		fprintf(stdout,"RT -> runtime\n");
		fprintf(stdout,"\n");
		
		exit(EXIT_FAILURE);
	}
	
	//set first image
	if(commandline) firstimage=11;
	else firstimage=2;
	
	// check number and size of images
	maxncol=0;
	maxnrow=0;
	nlay=ac-(firstimage+1);
	for(currentz=1;currentz<=nlay;currentz++) {
		ReadSizePNG(&nrow,&ncol,av[currentz+firstimage]);
		maxnrow=max(maxnrow,nrow);
		maxncol=max(maxncol,ncol);
	}
	nrow=maxnrow;
	ncol=maxncol; 
	
	//Allocate memory
	Allocate();
	
	//Allocate cells
	AllocateCells();
	
	// set parameters
	SetPars(av);
	
	// set output
	SetOutput(av[1]);
	
	// allocate images
	AllocateOutput(av);
	
	// open display
	//OpenDisplay(displayname,maxnrow,3*maxncol);
	
	printf("\nSetup\n");
	printf("\tnumber of images: %d\n",nlay);
	printf("\tmax image size is %d x %d pixels\n",ncol,nrow);
	printf("\tcomplete\n");
}

void Allocate()
{
	// allocate planes
	state=New();
	
	image=New();
	template=New();
	
	tempstate=New();
	tempstate2=New();
	
	energy=FNew();
	
	ftempstate=FNew();
	
	printf("\tplanes allocated\n");
	
	// allocate colors
	ColorFullRGB();
	printf("\tcolors allocated\n");
	printf("\tcomplete\n");
}

void AllocateCells()
{
	double c1,c2,c3;
	int r,g,b;
	
	printf("\nAllocateCells\n");
	
	CNew(&cells,maxcells,MAXTYPES);
	if((cells.extras=(void *)calloc((size_t)maxcells,sizeof(Extras)))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(EXIT_FAILURE);
	}
	InitCellPosition(&cells);
	
	// assign bright colors to cells
	EXTRAS[0].color=0;
	CELLS(cells,
			if(c>specialcelltype) {
				
				c1=0.5+0.5*RANDOM();
				c2=RANDOM();
				c3=0.5*RANDOM();
				
				if(RANDOM()<0.33) {
					r=rint(c1*255.);
					if(RANDOM()<0.5) {
						g=rint(c2*255.);
						b=rint(c3*255.);
					}
					else {
						g=rint(c3*255.);
						b=rint(c2*255.);
					}
				}
				else if(RANDOM()<0.5) {
					g=rint(c1*255.);
					if(RANDOM()<0.5) {
						r=rint(c2*255.);
						b=rint(c3*255.);
					}
					else {
						r=rint(c3*255.);
						b=rint(c2*255.);
					}
				}
				else {
					b=rint(c1*255.);
					if(RANDOM()<0.5) {
						r=rint(c2*255.);
						g=rint(c3*255.);
					}
					else {
						r=rint(c3*255.);
						g=rint(c2*255.);
					}
				}
				
				EXTRAS[c].color=RGB2Code(r,g,b);
			}
			);
	printf("\tallocated memory for %d cells\n",maxcells);
	printf("\tcomplete\n");
}

void AllocateOutput(char *av[])
{
	printf("\nAllocateOutput\n");
	// make output dir
	snprintf(a_string,STRING,"mkdir -p %s",outputdirname);
	system(a_string);
	
	//save parfile
	if(snprintf(a_string,STRING,"%s/%s",outputdirname,parfilename)>=STRING) { 
		fprintf(stderr,"warning: outputdirname too long: %s\n",a_string);
	}
	SaveOptions(av[2],a_string);
	
	// open pngs
	if(saveinit) {
		if(snprintf(imagedirname,STRING,"%s/CellstateInitial",outputdirname)>=STRING) { 
			fprintf(stderr,"warning: outputdirname too long: %s\n",imagedirname);
		}
		OpenPNG(&pngCellstateInit,imagedirname);	
	}
	
	if(saveinter) {
		if(nlay>1) {
			if(snprintf(imagedirname,STRING,"%s/CellstateIntermediate/%.5d",outputdirname,currentz-1)>=STRING) { 
				fprintf(stderr,"warning: outputdirname too long: %s\n",imagedirname);
			}
			OpenPNG(&pngCellstateInter,imagedirname);	
		}
		else {
			if(snprintf(imagedirname,STRING,"%s/CellstateIntermediate",outputdirname)>=STRING) { 
				fprintf(stderr,"warning: outputdirname too long: %s\n",imagedirname);
			}
			OpenPNG(&pngCellstateInter,imagedirname);	
		}
	}
	
	if(saveraw) {
		if(snprintf(imagedirname,STRING,"%s/CellstateRaw",outputdirname)>=STRING) { 
			fprintf(stderr,"warning: outputdirname too long: %s\n",imagedirname);
		}
		OpenPNG(&pngCellstateRaw,imagedirname);	
	}
	
	if(savedet) {
		if(snprintf(imagedirname,STRING,"%s/CellstateDeterministic",outputdirname)>=STRING) { 
			fprintf(stderr,"warning: outputdirname too long: %s\n",imagedirname);
		}
		OpenPNG(&pngCellstateDet,imagedirname);	
	}
	
	if(savefinal) {
		if(snprintf(imagedirname,STRING,"%s/CellstateFinal",outputdirname)>=STRING) { 
			fprintf(stderr,"warning: outputdirname too long: %s\n",imagedirname);
		}
		OpenPNG(&pngCellstateFinal,imagedirname);	
	}
	
	if(savecellcolors) {
		if(snprintf(imagedirname,STRING,"%s/Cellcolors",outputdirname)>=STRING) { 
			fprintf(stderr,"warning: outputdirname too long: %s\n",imagedirname);
		}
		OpenPNG(&pngCellcolors,imagedirname);	
	}
	
	if(savetemplate) {
		if(snprintf(imagedirname,STRING,"%s/Template",outputdirname)>=STRING) { 
			fprintf(stderr,"warning: outputdirname too long: %s\n",imagedirname);
		}
		OpenPNG(&pngTemplate,imagedirname);	
	}
	
	if(saveenergy) {
		if(snprintf(imagedirname,STRING,"%s/Energy",outputdirname)>=STRING) { 
			fprintf(stderr,"warning: outputdirname too long: %s\n",imagedirname);
		}
		OpenPNG(&pngEnergy,imagedirname);	
	}
	
	if(savecombined) {
		if(snprintf(imagedirname,STRING,"%s/Combined",outputdirname)>=STRING) { 
			fprintf(stderr,"warning: outputdirname too long: %s\n",imagedirname);
		}
		OpenPNG(&pngCombined,imagedirname);	
	}
	
	printf("\tcomplete\n");
}

void CloseOutput()
{
	if(saveinit) ClosePNG(&pngCellstateInit);	
	if(saveinter) ClosePNG(&pngCellstateInter);	
	if(saveraw) ClosePNG(&pngCellstateRaw);	
	if(savedet) ClosePNG(&pngCellstateDet);	
	if(savefinal) ClosePNG(&pngCellstateFinal);	
	if(savecellcolors) ClosePNG(&pngCellcolors);	
	if(savetemplate) ClosePNG(&pngTemplate);	
	if(saveenergy) ClosePNG(&pngEnergy);	
	if(savecombined) ClosePNG(&pngCombined);	
	
	//CloseDisplay();
}

void MyInDat(char filename[],char *format,char *text,void *value)
{
	if(!InDat(filename,format,text,value)) {
		printf("cannot find \"%s\" in parameter file; exitting now!\n",text);
		exit(EXIT_FAILURE);
	}
}

void ReadPars(char name[]) 
{
	char conversion[10];
	
#ifdef _SMALL
	snprintf(conversion,10,"%%f");
#else
	snprintf(conversion,10,"%%lf");
#endif
	
	printf("\nReadPars\n");
	printf("\tReading parameters from: %s\n",name);
	
	ReadOptions(name);
	//MyInDat(name,"%d","",&);
	//MyInDat(name,conversion,"",&);
	
	MyInDat(name,"%d","hamiltonianneigh",&hamiltonianneigh);
	MyInDat(name,"%d","imageinteraction",&imageinteraction);
	MyInDat(name,"%d","targetcoverage",&targetcoverage);
	MyInDat(name,"%d","estimationoffset",&estimationoffset);
	MyInDat(name,"%d","temperaturescale",&temperaturescale);
	MyInDat(name,"%d","persistencestrength",&persistencestrength);
	
	MyInDat(name,"%d","backgroundsubtraction",&backgroundsubtraction);
	MyInDat(name,"%d","gaussianblur",&gaussianblur);
	
	MyInDat(name,"%d","runtime",&runtime);
	
	MyInDat(name,"%d","seedgaussian",&seedgaussian);
	
	MyInDat(name,"%d","mincellarea",&mincellarea);
	
	MyInDat(name,"%d","persistencelength",&persistencelength);
	MyInDat(name,"%d","persistencescale",&persistencescale);
	
	MyInDat(name,"%d","maxcells",&maxcells);
	MyInDat(name,"%d","commandline",&commandline);
	
	MyInDat(name,"%d","cellmarker",&cellmarker);
	MyInDat(name,"%d","seedmarker",&seedmarker);
	
	MyInDat(name,"%d","includepars",&includepars);
	MyInDat(name,"%d","includetime",&includetime);
	
	MyInDat(name,"%d","displayinterval",&displayinterval);
	MyInDat(name,"%d","infointerval",&infointerval);
	
	MyInDat(name,"%d","saveinit",&saveinit);
	MyInDat(name,"%d","saveinter",&saveinter);
	MyInDat(name,"%d","saveraw",&saveraw);
	MyInDat(name,"%d","savedet",&savedet);
	MyInDat(name,"%d","savefinal",&savefinal);
	
	MyInDat(name,"%d","savecellcolors",&savecellcolors);
	MyInDat(name,"%d","savetemplate",&savetemplate);
	MyInDat(name,"%d","saveenergy",&saveenergy);
	MyInDat(name,"%d","savecombined",&savecombined);
	
	printf("\tcomplete\n");
}

void SetOutput(char name[])
{
	int pathlength;
	
	time_t sec;
	
	sec = time(NULL);
	
	printf("\nSetOutput\n");
	
	pathlength=strlen(name);
	if(name[pathlength-1]=='/') name[pathlength-1]='\0';
	name[0]=toupper(name[0]);
	
	snprintf(parfilename,STRING,"settings.par");
	
	if(includepars) {
		snprintf(a_string,STRING,"%s/SPM2D_HN%dII%dTC%dEO%dTS%dPS%dBS%dGB%dRT%d",name,hamiltonianneigh,imageinteraction,targetcoverage,estimationoffset,temperaturescale,persistencestrength,backgroundsubtraction,gaussianblur,runtime);
	}
	else {
		snprintf(a_string,STRING,"%s/SPM2D",name);
	}
	
	
	if(includetime) {
		if(snprintf(outputdirname,STRING,"%s_%.12ld",a_string,sec)>=STRING) { 
			fprintf(stderr,"warning: outputdirname too long: %s\n",outputdirname);
		}
	}
	else {
		if(snprintf(outputdirname,STRING,"%s",a_string)>=STRING) { 
			fprintf(stderr,"warning: outputdirname too long: %s\n",outputdirname);
		}
	}
	
	printf("\toutputdirname: %s\n",outputdirname);
	
	snprintf(displayname,STRING,"SPM2D :: %s",outputdirname);
	printf("\tcomplete\n");
}

void SetPars(char *av[]) 
{
	int neighsize;
	int targetchance0;
	
	printf("\nSetPars\n");
	
	// set unused CPM pars to zero
	lambda=0.0;
	dissipation=0.0;
	targetarea=0;
	
	// read pars from commandline
	if(commandline) { 
		
		hamiltonianneigh=atoi(av[3]);
		imageinteraction=atoi(av[4]);
		targetcoverage=atoi(av[5]);
		estimationoffset=atoi(av[6]);
		temperaturescale=atoi(av[7]);
		persistencestrength=atoi(av[8]);
		
		backgroundsubtraction=atoi(av[9]);
		gaussianblur=atoi(av[10]);
		
		runtime=atoi(av[11]);
	}
	
	// set misc pars
	if(!graphics) displayinterval=0;
	if(!persistencestrength) persistencelength=0;
	else if (!persistencelength || !persistencescale) { 
		fprintf(stdout,"ERROR: persistence is used, but \'persistencelength\' and/or \'persistencescale\' is specified\n");
		exit(EXIT_FAILURE);
	}
	boundary=FIXED;
	scale=2;
	
	//set J-values
	cells.J[CELLWALL][CELLWALL]=0;
	cells.J[CELLWALL][CELL]=1000;
	cells.J[CELL][CELLWALL]=1000;
	cells.J[CELL][CELL]=2000;
	
	// calculate relative pars
	if(hamiltonianneigh) neighsize=(double)NeighbourNumber(hamiltonianneigh);
	else neighsize=(double)NeighbourNumber(2);
	
	totalj=neighsize*cells.J[CELL][CELL];
	
	mu=(double)imageinteraction*(double)totalj*0.01; 
	nu=(double)persistencestrength*(double)totalj*0.01;
	
	printf("\timageinteraction: %d -> mu: %f\n",imageinteraction,mu);
	printf("\tpersistencestrength: %d -> nu: %f\n",persistencestrength,nu);
	
	//temperature is adjusted such that chance0 for specified percentage of worst possible change
	temperature=500.;
	targetchance0=rint((double)temperaturescale*(double)totalj*0.01);
	while(abs(cells.chance0-targetchance0) && temperature>0.0) {
		temperature+=1.0;
		temperature/=((double)cells.chance0/(double)targetchance0);
		temperature-=1.0;
		temperature=max(0.0,temperature);
		ResetBoltzmann(&cells);
	}
	printf("\ttemperaturescale: %d -> temperature: %f\n",temperaturescale,temperature);
	printf("\tseed: %d\n",seed);
	printf("\tcomplete\n");
}

//SPM functions
void SetSegmentation()
{
	printf("\nSetSegmentation\n");
	totaltargetarea=(int)(tcoverage*(double)npixels);
	printf("\ttotaltargetarea: %d\n",totaltargetarea);
	
	areaconstraint=0.0; //kappa is lambda for total targetarea
	coverage=0.;
	
	MemoryFree();
	
	//open new directory for saving intermediate images
	if(currentz>1 && nlay>1 && saveinter) {
		ClosePNG(&pngCellstateInter);	
		if(snprintf(imagedirname,STRING,"%s/CellstateIntermediate/%.5d",outputdirname,currentz-1)>=STRING) { 
			fprintf(stderr,"warning: outputdirname too long: %s\n",imagedirname);
		}
		OpenPNG(&pngCellstateInter,imagedirname);	
	}
	
	//set neigh for Hamiltonian
	if(hamiltonianneigh) NeighbourNumber(hamiltonianneigh);
	else NeighbourNumber(2);	
	printf("\tcomplete\n");
}

void SetCellseed()
{	
	int i;
	int j;
	int signal;
	
	signal=0;
	
	if(seedmarker==-1) {
		PLANE(
				ftempstate[i][j]=(double)(template[i][j]);
				if(ftempstate[i][j]>0.0) signal=1;					  
				);
	}
	
	else if(seedmarker==cellmarker) {
		PLANE(
				ftempstate[i][j]=(double)(Code2Color(image[i][j],cellmarker));
				if(ftempstate[i][j]>0.0) signal=1;
				);
	}
	
	else {
		PLANE(
				ftempstate[i][j]=(double)(255.-Code2Color(image[i][j],seedmarker));
				if(ftempstate[i][j]>0.0) signal=1;
				);
	}
	
	if(!signal) { 
		fprintf(stdout,"ERROR: specified colorchannel for seed marker is empty, set correct color channel in parfile\n");
		exit(EXIT_FAILURE);
	}
	
	FGaussianFilter2D(ftempstate,seedgaussian);	
	FFindLocalMinima2D(state,ftempstate,tempstate);
	RenameCells2D(state);
	
	printf("\nSetCellseed\n");
	
	ncells=1;
	PLANE(ncells=max(ncells,state[i][j]););
	
	if(ncells>=maxcells) {
		fprintf(stdout,"ERROR: not enough cells allocated, please increase maxcells and try again. Allocated: %d, required: %d\n",maxcells,ncells+1);
		exit(EXIT_FAILURE);
	}
	
	// set cells
	UpdateCFill(state,&cells);
	UpdateCellPosition(state,&cells);
	UpdateCellShape(&cells);
	
	CELLS(cells,
			if(c>specialcelltype && cells.area[c]) {
				
				cells.celltype[c]=CELL;
				
				// make sure cells have at least size 9
				i=rint(cells.shape[c].meany);
				j=rint(cells.shape[c].meanx);
				
				NEIGHBOURS( 
							  if(!(i+y<1 || i+y>nrow || j+x<1 || j+x>ncol)) {
								  if(!state[i+y][j+x]) {
									  state[i+y][j+x]=c;
								  }
							  }
							  );
			}
			);
	
	UpdateCFill(state,&cells);
	UpdateCellPosition(state,&cells);
	UpdateCellShape(&cells);
	
	CELLS(cells,
			if(c>specialcelltype && cells.area[c]) {
				
				EXTRAS[c].vectorx = 0;
				EXTRAS[c].vectory = 0; 
				
				EXTRAS[c].originx = cells.shape[c].meanx;
				EXTRAS[c].originy = cells.shape[c].meany;
				
				EXTRAS[c].nu=0.0;
			}
			);
	
	if(saveinit) PlanePNG(state,&pngCellstateInit,0);
	
	printf("\t%d cells initiated\n",ncells);
	printf("\tcomplete\n");}

void SetTemplate()
{
	int i,threshold;
	double histogram[256];
	double energyfunction[256];
	int signal;
	double minval;
	
	npixels=nrow*ncol;
	pixelsize=1./(double)(nrow*ncol);
	
	signal=0;
	PLANE(
			template[i][j]=Code2Color(image[i][j],cellmarker);
			if(template[i][j]) signal=1;
			);
	
	if(!signal) { 
		fprintf(stdout,"ERROR: specified colorchannel for cell marker is empty, set correct color channel in parfile\n");
		exit(EXIT_FAILURE);
	}
	
	if(backgroundsubtraction) BackgroundSubtraction2D(template,backgroundsubtraction,tempstate);
	
	PLANE(ftempstate[i][j]=(double)template[i][j];);
	
	if(gaussianblur) FGaussianFilter2D(ftempstate,gaussianblur);
	
	FStretchHistogram2D(ftempstate);
	
	PLANE(
			template[i][j]=rint(ftempstate[i][j]);
			template[i][j]=min(max(0,template[i][j]),255);
			);
	
	if(!targetcoverage) {
		tcoverage=EstimateTargetCoverage2D(template,tempstate);
		tcoverage+=(double)estimationoffset*0.01;
	}
	else tcoverage=(double)targetcoverage*0.01;
	
	printf("\nSetTemplate\n");
	//update offset
	for(i=0;i<=255;i++) {
		histogram[i]=0.;
	}
	
	PLANE(histogram[template[i][j]]+=pixelsize;);
	
	for(i=1;i<=255;i++) {
		histogram[i]+=histogram[i-1];
	}
	
	threshold=0;
	minval=0.0;
	for(i=0;i<=255;i++) {
		if(histogram[255-i]>0.0) minval=histogram[255-i];
		if(histogram[i]<tcoverage) { 
			threshold=i;
		}
	}
	tcoverage=histogram[threshold];
	printf("\tOffset intensity: %d, Coverage: %f, Min Coverage: %f\n",threshold,histogram[threshold],minval);
	
	//energy function -> goes from -1 (signal=0) via 0 (signal=offset) to 1 (signal=255)
	for(i=0;i<=255;i++) {
		if(histogram[i]<histogram[threshold]) {
			energyfunction[i]=histogram[i]-minval;
			energyfunction[i]/=(histogram[threshold]-minval);
			energyfunction[i]-=1;
		}
		else {
			energyfunction[i]=histogram[i]-histogram[threshold];
			energyfunction[i]/=(1.-histogram[threshold]);
		}
	}
	
	PLANE(
			signal=template[i][j];
			energy[i][j]=energyfunction[signal];
			);
	
	printf("\tcomplete\n");
}

void Extra()
{
	double xshift,yshift;
	double hypotxy;
	double signal;
	
	signal=energy[iposition][jposition]; //intensity in normalised images
	
	if(sigma) { 
		dh-=(int)(mu*signal);
		if(totalarea<=totaltargetarea)
			dh+=kappa;
		else dh-=kappa;
	}
	
	if(sigmaneigh) { 
		dh+=(int)(mu*signal);
		if(totalarea>=totaltargetarea)
			dh+=kappa;
		else dh-=kappa;
	}
	
	if(sigma && !sigmaneigh) { 
		if(signal<0. && cells.area[sigma]<=mincellarea) { //don't shrink beyond minsize if signal is negative
			dh=cells.chance0;
		}
	}
	
	if(sigma && sigmaneigh) {
		
		if(persistencestrength && signal<0.0 && (tijd*10)<(runtime*9)) { //persistence
			//persistence sigma
			xshift = (double)jposition-cells.shape[sigma].meanx;
			yshift = (double)iposition-cells.shape[sigma].meany;
			// cos alpha = x.y/vector length 
			hypotxy=hypot(xshift,yshift);
			// check if hypotxyz is not too close to 0, if so the cell is probably dying 
			if (!(hypotxy < EPSILON)) {
				dh+=(int)(EXTRAS[sigma].nu*-signal*(xshift*EXTRAS[sigma].vectorx+yshift*EXTRAS[sigma].vectory)/hypotxy);
			}
			//persistence sigmaneigh
			xshift = (double)jposition-cells.shape[sigmaneigh].meanx;
			yshift = (double)iposition-cells.shape[sigmaneigh].meany;
			// cos alpha = x.y/vector length */
			hypotxy=hypot(xshift,yshift);
			// check if hypotxyz is not too close to 0, if so the cell may be dying */
			if (!(hypotxy < EPSILON)) {
				dh-=(int)(EXTRAS[sigmaneigh].nu*-signal*(xshift*EXTRAS[sigmaneigh].vectorx+yshift*EXTRAS[sigmaneigh].vectory)/hypotxy);
			} 
		}
	}
	
}

void Aftermath()
{
	if(sigma) totalarea--;
	if(sigmaneigh) totalarea++;
	if(persistencestrength) CellPosition(&cells,iposition,jposition,sigma,sigmaneigh);
}

void UpdateSegmentation()
{	
	double dcoverage,dtcoverage,tempcoverage;
	
	ncells=0;
	CELLS(cells,
			if(c>specialcelltype && cells.area[c]) {
				ncells++;
			}
			);
	
	tempcoverage=coverage;
	coverage=0.;
	noiselevel=0.;
	totalarea=0;
	
	// update coverage and noiselevel
	PLANE(
			if(state[i][j]) {
				coverage+=pixelsize;
				totalarea++;
				if(energy[i][j]>0.0 && state[i][j]) noiselevel+=pixelsize;
			}
			);
	
	dcoverage=coverage-tcoverage; // difference between coverage and targetcoverage
	dtcoverage=coverage-tempcoverage; // change in coverage
	
	// update kappa -> energy contribution of total targetarea
	if(dcoverage<-0.01 && dtcoverage<0.01) areaconstraint+=0.001; // increase kappa when coverage too small and decreasing
	else if(dcoverage>0.01 && dtcoverage>-0.01) areaconstraint+=0.001; // increase kappa when coverage too big and increasing
	else kappa=max(0.0,kappa-0.001); // decrease kappa otherwise
	
	kappa=rint((double)totalj*areaconstraint);
}

void UpdatePersistence()
{
	double hypotxy;
	double vectorlength;
	
	//update persistence
	CELLS(cells,
			if(c>specialcelltype && cells.area[c]) {
				hypotxy = hypot(cells.shape[c].meanx-EXTRAS[c].originx,cells.shape[c].meany-EXTRAS[c].originy);
				if (hypotxy > EPSILON) {
					EXTRAS[c].vectorx = (cells.shape[c].meanx-EXTRAS[c].originx)/hypotxy;
					EXTRAS[c].vectory = (cells.shape[c].meany-EXTRAS[c].originy)/hypotxy;
					EXTRAS[c].originx = cells.shape[c].meanx;
					EXTRAS[c].originy = cells.shape[c].meany;
				}
				vectorlength=min(1.0,hypotxy/(double)persistencescale);
				EXTRAS[c].nu=vectorlength*nu;
			}
			);
}

//IO functions

void PrintInfo()
{	
	printf("\n");
	printf("\tTijd: %d/%d\t\n",tijd,runtime);
	printf("\tCoverage: %f/%f\t\n",coverage,tcoverage);
	printf("\tTotalarea: %d, totaltargetarea: %d\n",totalarea,totaltargetarea);
	printf("\tNumber of cells: %d\n",ncells);
	printf("\tAreaconstraint: %f, lambda: %d\n",areaconstraint,kappa);
	printf("\tNoiselevel: %f\n",noiselevel);
	printf("\n");
}

void SaveGraphics()
{
	int x,y;
	
	printf("\nSaveGraphics\n");
	
	PLANE(
			tempstate[i][j]=EXTRAS[state[i][j]].color;
			);
	PLANE(
			for(x=0;x<=1;x++) {
				for(y=0;y<=1;y++) {
					if(!(i+y>nrow || j+x>ncol)) {
						if(state[i][j]!=state[i+y][j+x]) {
							tempstate[i][j]=RGB2Code(0,0,255);
						}
					}
				}
			}
			);
	
	if(savecellcolors) PlanePNG(tempstate,&pngCellcolors,0);
	if(savecombined) {
		PLANE(
				tempstate2[i][j]=MixColors(image[i][j],tempstate[i][j],(double)savecombined/100.);
				);
		PlanePNG(tempstate2,&pngCombined,0); 
	}
	
	if(savetemplate) {
		PLANE(
				tempstate[i][j]=Color2Code(template[i][j],0);
				);
		PlanePNG(tempstate,&pngTemplate,0); 
	}
	
	if(saveenergy) {
		PLANE(
				if(energy[i][j]<0.0) tempstate[i][j]=RGB2Code(0,-(rint)(255.*energy[i][j]),0);
				else tempstate[i][j]=RGB2Code((rint)(255.*energy[i][j]),0,(rint)(255.*energy[i][j]));
				);
		PlanePNG(tempstate,&pngEnergy,0);
	}
	printf("\tcomplete\n");
}

//Support functions
void RemoveStochasticity2D(TYPE **target, FLOAT **template)
{
	int nremoved,delta;
	int c,x,y;
	
	double distance,mindistance;
	
	UpdateCellPosition(state,&cells);
	UpdateCellShape(&cells);
	
	printf("\nRemoveStochasticity2D\n");
	
	//remove cells from high intensity pixels
	nremoved=0;
	PLANE(
			if(state[i][j]) {
				if(template[i][j]>0.0) {
					state[i][j]=0;
					nremoved++;
				}
			}
			);
	printf("\t%d pixels at high intensity removed\n",nremoved);
	
	//remove loose ends
	UpdateCFill(state,&cells);
	UpdateCellPosition(state,&cells);
	UpdateCellShape(&cells);
	
	PLANE(tempstate[i][j]=0;);
	
	CELLS(cells,
			if(c>specialcelltype && cells.area[c]) {
				x=rint(cells.shape[c].meanx);
				y=rint(cells.shape[c].meany);
				if(state[y][x]==c) {
					tempstate[y][x]=c;
				}
				else {
					mindistance=hypot((double)nrow,(double)ncol);
					PLANE(
							if(state[i][j]==c) {
								distance=hypot(cells.shape[c].meanx-(double)j,cells.shape[c].meany-(double)i);
								if(distance<mindistance) {
									mindistance=distance;
									x=j;
									y=i;
								}
							}
							);
					tempstate[y][x]=c;
				}
			}	 
			);
	
	delta=1;
	while(delta) {
		delta=0;
		PLANE(
				if(state[i][j] && !tempstate[i][j]) {
					c=state[i][j];
					NEIGHBOURS(
								  if(x*y==0 && (x || y)) {
									  if(!(i+y<1 || i+y>nrow || j+x<1 || j+x>ncol)) {
										  if(tempstate[i+y][j+x]==c) {
											  delta++;
											  tempstate[i][j]=c;
										  }
									  }
								  }
								  );
				}
				);
	}
	
	nremoved=0;
	PLANE(
			if(state[i][j] && !tempstate[i][j]) nremoved++;
			state[i][j]=tempstate[i][j];
			);
	
	printf("\t%d disconnected pixels removed\n",nremoved);
	printf("\tcomplete\n");
}

void Watershed2D(TYPE **target, TYPE **template, TYPE **tempstate) 
{
	int start;
	int delta;
	int maxval;
	
	printf("\nWatershed2D\n");
	
	//find starting intensity for watershed
	start=0;
	PLANE(
			if(!target[i][j]) delta=max(start,template[i][j]);
			);
	printf("\tstarting water level: %d\n",start);
	
	for(maxval=start;maxval<=255;maxval++) {
		delta=1;
		while(delta) {
			
			PLANE(
					tempstate[i][j]=target[i][j];
					);
			
			PLANE(
					if(!target[i][j]) {
						
						//random check of first order neighbours
						ASYNNEIGHBOURS(
											if(x*y==0) {
												//printf("Code check -> x: %d, y: %d, z: %d\n",x,y,z);
												if(!(i+y<1 || i+y>nrow || j+x<1 || j+x>ncol)) {
													if(target[i+y][j+x] && template[i+y][j+x]<=maxval) {
														tempstate[i][j]=target[i+y][j+x];
													}
												}
											}
											)
					}
					);
			
			//check for changes
			delta=0;
			PLANE(
					if(target[i][j]!=tempstate[i][j]) delta++;
					target[i][j]=tempstate[i][j];
					);
		}
	}
	printf("\tcomplete\n");
}

void RemoveProtrutions2D(TYPE **target, TYPE **tempstate) 
{
	int nswitched;
	int n,nneigh,delta,self;
	
	int neighid[9];
	int neighsize[9];
	
	printf("\nRemoveProtrutions2D\n");
	nswitched=0;
	
	delta=1;
	while(delta) {
		
		PLANE(tempstate[i][j]=target[i][j];);
		
		PLANE(
				self=target[i][j];
				nneigh=0;
				NEIGHBOURS(
							  if(!(i+y<1 || i+y>nrow || j+x<1 || j+x>ncol)) {
								  if(target[i+y][j+x]!=self) nneigh++;
							  }
							  );
				
				if(nneigh>=5) {
					
					for(n=1;n<=8;n++) {
						neighid[n]=0;
						neighsize[n]=0;
					}
					
					NEIGHBOURS(
								  if(!(i+y<1 || i+y>nrow || j+x<1 || j+x>ncol)) {
									  if(target[i+y][j+x]!=self) {
										  delta=0;
										  for(n=1;n<=8;n++) {
											  if(delta==0 && (neighid[n]==0 || neighid[n]==target[i+y][j+x])) {
												  delta=1;
												  neighid[n]=target[i+y][j+x];
												  neighsize[n]++;
											  }
										  }
									  }
								  }
								  );
					
					for(n=1;n<=8;n++) {
						if(neighsize[n]>=5) {
							tempstate[i][j]=neighid[n];
						}
					}
				}
				);
		
		delta=0;
		PLANE(
				if(target[i][j]!=tempstate[i][j]) delta++;
				target[i][j]=tempstate[i][j];
				);
		nswitched+=delta;
	}
	printf("\t%d pixels switched\n",nswitched);
	printf("\tcomplete\n");
}

void RenameCells2D(TYPE **target)
{
	// renames the cells from 1 to ncells
	
	int c;
	int cellnumber;
	int maxc;
	
	printf("\nRenameCells2D\n");
	
	maxc=1;
	PLANE(
			maxc=max(maxc,target[i][j]);
			);
	
	//int cellid[maxc+1];
    int *cellid = (int *) malloc(sizeof(int) * (maxc+1));
	
	for(c=1;c<=maxc;c++) {
		cellid[c]=0;
	}
	
	cellnumber=0;
	PLANE(
			c=target[i][j];
			if(c && cellid[c]==0) {
				cellnumber++;
				cellid[c]=cellnumber;
			}
			);
	
	printf("\tnumber of cells: %d\n",cellnumber);
	
	PLANE(
			c=target[i][j];
			if(c) {
				target[i][j]=cellid[c];
			}
			);
	
	printf("\tcomplete\n");
}

double EstimateTargetCoverage2D(TYPE **target, TYPE **tempstate)
{
	printf("\nFindCoverage2D\n");
	
	int threshold;
	int checkpoint1;
	int checkpoint2;
	int checkpoint3;
	int c,neighc;
	int delta;
	
	int nclusters[256];
	double tcoverage[256];
	
	checkpoint1=0;
	checkpoint2=0;
	checkpoint3=0;
	
	nclusters[0]=0;
	tcoverage[0]=1.0;
	
	threshold=1;
	
	while(!checkpoint3 && threshold<=256) {
		
		c=0;
		tcoverage[threshold]=0.0;
		
		PLANE(tempstate[i][j]=0;);
		
		PLANE(
				if(target[i][j]>=threshold) {
					neighc=0;
					NEIGHBOURS(
								  if(!(i+y<1 || i+y>nrow || j+x<1 || j+x>ncol)) {
									  if(tempstate[i+y][j+x]) {
										  neighc=tempstate[i+y][j+x];
									  }
								  }
								  );
					if(neighc) tempstate[i][j]=neighc;
					else {
						c++;
						tempstate[i][j]=c;
					}
				}
				else {
					tcoverage[threshold]+=pixelsize;
				}
				);
		
		delta=1;
		while(delta) {
			delta=0;
			PLANE(
					if(tempstate[i][j]) {
						NEIGHBOURS(
									  if(!(i+y<1 || i+y>nrow || j+x<1 || j+x>ncol)) {
										  if(tempstate[i+y][j+x]) {
											  if(tempstate[i+y][j+x]!=tempstate[i][j]) {
												  delta++;
												  tempstate[i][j]=min(tempstate[i][j],tempstate[i+y][j+x]);
												  tempstate[i+y][j+x]=min(tempstate[i][j],tempstate[i+y][j+x]);
											  }
										  }
									  }
									  );
					}
					);
		}
		
		nclusters[threshold]=CountClusters2D(tempstate);
		
		if(!checkpoint1) {
			if(nclusters[threshold]<nclusters[threshold-1]) {
				checkpoint1=threshold-1;
			}
		}
		
		else if(!checkpoint2) {
			if(nclusters[threshold]>nclusters[threshold-1]) {
				checkpoint2=threshold-1;
			}
		}
		
		else {
			if(nclusters[threshold]<nclusters[threshold-1]) {
				checkpoint3=threshold-1;
			}
		}
		printf("\tthreshold:  %d, nclusters: %d, coverage: %f\n",threshold,nclusters[threshold],tcoverage[threshold]);
		threshold++;
	}
	
	
	if(!(checkpoint1 && checkpoint2 && checkpoint3)) {
		fprintf(stdout,"ERROR: image not suitable for estimating target coverage, set target coverage manually\n");
		exit(EXIT_FAILURE);
	}
	
	else printf("\testimated target coverage: %f (threshold:  %d, number of clusters: %d\n",tcoverage[checkpoint2],checkpoint2,nclusters[checkpoint2]);
	printf("\tcomplete\n");
	return tcoverage[checkpoint2];
}

void FFindLocalMinima2D(TYPE **target, FLOAT **template, TYPE **tempstate)
{
	int c;
	int ncells;
	int delta;
	int localmin;
	
	printf("\nFFindLocalMinima2D\n");
	
	PLANE(target[i][j]=0;);
	ncells=1;
	
	PLANE(
			localmin=1;
			NEIGHBOURS( 
						  if(!(i+y<1 || i+y>nrow || j+x<1 || j+x>ncol)) {
							  if(template[i][j]>template[i+y][j+x]) {
								  localmin=0;
							  }
						  }
						  );
			
			if(localmin) {
				NEIGHBOURS( 
							  if(!(i+y<1 || i+y>nrow || j+x<1 || j+x>ncol)) {
								  if(!target[i][j]) target[i][j]=target[i+y][j+x];
								  else if(target[i+y][j+x]) {
									  target[i][j]=min(target[i][j],target[i+y][j+x]);
									  target[i+y][j+x]=min(target[i][j],target[i+y][j+x]);
								  }
							  }
							  );
				if(!target[i][j]) {
					target[i][j]=ncells;
					ncells++;
				}
			}
			);
	
	delta=1;
	while(delta) {
		PLANE(tempstate[i][j]=target[i][j];);
		PLANE(
				if(target[i][j]) {
					NEIGHBOURS( 
								  if(!(i+y<1 || i+y>nrow || j+x<1 || j+x>ncol)) {
									  if(target[i+y][j+x]) {
										  target[i][j]=min(target[i][j],target[i+y][j+x]);
										  target[i+y][j+x]=min(target[i][j],target[i+y][j+x]);
									  }
								  }
								  );
				}
				);
		delta=0;
		PLANE(
				if(target[i][j]!=tempstate[i][j]) {
					delta++;
				}
				);
	}	
	
	//filter for false minima
	//int nomin[ncells+1];
	int *nomin = (int *) malloc(sizeof(int) * (ncells+1));
	
	for(c=1;c<=ncells;c++) {
		nomin[c]=0;
	}
	PLANE(
			if(target[i][j]) {
				NEIGHBOURS( 
							  if(!(i+y<1 || i+y>nrow || j+x<1 || j+x>ncol)) {
								  if(!target[i+y][j+x] && template[i][j]>=template[i+y][j+x]) {
									  nomin[target[i][j]]=1;
								  }
							  }
							  );
			}
			);
	
	PLANE(
			if(target[i][j] && nomin[target[i][j]]) {
				target[i][j]=0;
			}
			);
	
	printf("\t%d local minima found\n",ncells);
	printf("\tcomplete\n");
}

void BackgroundSubtraction2D(TYPE **target, int neighsize, TYPE **tempstate)
{
	double histogram[256];
	int c,npixels,median;
	double delta,delta2;
	
	printf("\nBackgroundSubtraction2D\n");
	
	NeighbourNumber(neighsize);
	PLANE(
			for(c=0;c<=255;c++) {
				histogram[c]=0.0;
			}
			histogram[target[i][j]]=1.0;
			npixels=1;
			
			WIDENEIGHBOURS( 
								if(!(i+y<1 || i+y>nrow || j+x<1 || j+x>ncol)) {
									histogram[target[i+y][j+x]]+=1.0;
									npixels++;
								}
								);
			
			for(c=0;c<=255;c++) {
				histogram[c]/=(double)npixels;
			}
			
			for(c=1;c<=255;c++) {
				histogram[c]+=histogram[c-1];
			}
			
			median=0;
			for(c=0;c<=254;c++) {
				delta=fabs(histogram[c]-0.5);
				delta2=fabs(histogram[c+1]-0.5);
				if(delta<delta2 && median==0) median=c;
			}
			
			tempstate[i][j]=target[i][j]-median;
			tempstate[i][j]=max(0,tempstate[i][j]);
			);
	
	PLANE(
			target[i][j]=tempstate[i][j];
			);
	printf("\tcomplete\n");
}

void FGaussianFilter2D(FLOAT **target, int kernelsize)
{
	double gauss_weight;
	double stdev;
	//double tempresult[max(nrow,ncol)+1];
	double totalweight;
	//double weight[kernelsize+1];
	
    double *tempresult = (double *) malloc(sizeof(double) * (max(nrow, ncol)+1));
    double *weight = (double *) malloc(sizeof(double) * (kernelsize+1));
	
	int i;
	int j;
	int x;
	int y;
	
	printf("\nFGaussianFilter2D\n");
	
	stdev=(2*(double)kernelsize+1.)/6.;
	
	totalweight=0;
	
	for(i=0;i<=kernelsize;i++) {
		weight[i]=(1./(sqrt(2.*M_PI)*stdev))*pow(M_E,-((i*i)/(2.*stdev*stdev)));
		totalweight+=weight[i];
	}
	
	totalweight*=2.;
	totalweight-=weight[0];
	
	printf("\tkernelsize: %d, stdev: %f, totalweight: %f\n",kernelsize*2+1,stdev,totalweight);
	
	// Filter in x direction
	for(i=1;i<=nrow;i++) {
		for(j=1;j<=ncol;j++) {
			tempresult[j]=0.0;
			gauss_weight=0.0;
			for(x=-kernelsize;x<=kernelsize;x++) {
				if(j+x>=1 && j+x<=ncol) {
					if(!(target[i][j+x]<0.0)) {
						tempresult[j]+=target[i][j+x]*weight[abs(x)];
						gauss_weight+=weight[abs(x)]; 
					}
				}
			}
			tempresult[j]/=gauss_weight;
		}
		for(j=1;j<=ncol;j++) target[i][j]=tempresult[j];
	}
	printf("\tfiltered in x direction\n");
	
	// Filter in y direction
	for(j=1;j<=ncol;j++) {
		for(i=1;i<=nrow;i++) {
			tempresult[i]=0.0;
			gauss_weight=0.0;
			for(y=-kernelsize;y<=kernelsize;y++) {
				if(i+y>=1 && i+y<=nrow) {
					if(!(target[i+y][j]<0.0)) {
						tempresult[i]+=target[i+y][j]*weight[abs(y)];
						gauss_weight+=weight[abs(y)]; 
					}
				}
			}
			tempresult[i]/=gauss_weight;
		}
		for(i=1;i<=nrow;i++) target[i][j]=tempresult[i];
	}
	printf("\tfiltered in y direction\n");
	printf("\tcomplete\n");
}

void FStretchHistogram2D(FLOAT **target)
{
	double maxval;
	double minval;
	
	printf("\nFStretchHistogram2D\n");
	
	minval=255.;
	maxval=0.;
	
	PLANE(
			minval=min(minval,target[i][j]);
			maxval=max(maxval,target[i][j]);
			);
	
	if(maxval-minval>EPSILON) { 
		printf("\tminval: %f, maxval: %f\n",minval,maxval);
		PLANE(
				target[i][j]=255.*((target[i][j]-minval)/(maxval-minval));
				);
	}
	else {
		printf("\tno variation in histogram\n");
	}
	
	printf("\tcomplete\n");
}

int CountClusters2D(TYPE **target)
{
	int c,nc,maxc;
	
	maxc=0;
	PLANE(maxc=max(maxc,target[i][j]););
	
	//int size[maxc+1];
	int *size = (int *) malloc(sizeof(int) * (maxc+1));
	
	for(c=0;c<=maxc;c++) {
		size[c]=0;
	}
	
	PLANE(
			c=target[i][j];
			size[c]++;
			);
	
	nc=0;
	for(c=1;c<=maxc;c++) {
		if(size[c]) nc++;
	}
	
	return nc;
}

//color functions
int Code2Color(int state, int which)
{
	int red,green,blue,white;
	
	state=max(0,min(16777215,state));
	which=max(0,min(3,which));
	
	blue=state%256;
	if(which==3) return blue;
	else {
		green=((state-blue)/256)%256;
		if(which==2)return green;
		else {
			red=(state-blue-256*green)/65536;
			if(which==1)return red;
			else {
				white=(int)((double)(red+green+blue)/3.);
				return white;
			}
		}
	}
}

int Color2Code(int intensity, int which)
{
	intensity=max(0,min(255,intensity));
	which=max(0,min(3,which));
	
	if(which==0) return intensity+256*intensity+65536*intensity;
	else if(which==1) return 65536*intensity;
	else if(which==2) return 256*intensity;
	else return intensity;
}

int RGB2Code(int r,int g,int b)
{
	r=max(0,min(255,r));
	g=max(0,min(255,g));
	b=max(0,min(255,b));
	return b+256*g+65536*r;
}

int MixColors(int color1, int color2, double opacity)
{
	int red1,green1,blue1;
	int red2,green2,blue2;
	int r,g,b;
	
	blue1=color1%256;
	green1=((color1-blue1)/256)%256;
	red1=(color1-blue1-256*green1)/65536;
	
	blue2=color2%256;
	green2=((color2-blue2)/256)%256;
	red2=(color2-blue2-256*green2)/65536;
	
	r=(int)((double)red1*(1.-opacity)+(double)red2*opacity);
	g=(int)((double)green1*(1.-opacity)+(double)green2*opacity);
	b=(int)((double)blue1*(1.-opacity)+(double)blue2*opacity);
	
	r=max(0,min(255,r));
	g=max(0,min(255,g));
	b=max(0,min(255,b));
	
	return b+256*g+65536*r;
}
