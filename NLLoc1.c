/*
 * Copyright (C) 1999-2000 Anthony Lomax <lomax@faille.unice.fr>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


/*   NLLoc.c

	Program to do non-linear earthquake location in 3-D grid models

*/

/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */


/*
	history:

	ver 01    26SEP1997  AJL  Original version
	ver 02    08JUN1998  AJL  Metropolis added


.........1.........2.........3.........4.........5.........6.........7.........8

*/


/* References */
/*
	TV82	Tarantola and Valette,  (1982)
		"Inverse Problems = Quest for Information",
		J Geophys 50, 159-170.
	MEN92	Moser, van Eck and Nolet,  (1992)
		"Hypocenter Determination ... Shortest Path Method",
		JGR 97, B5, 6563-6572.
*/




#define PNAME  "NLLoc"

#include "GridLib.h"
#include "ran1.h"
#include "octtree.h"


/* defines */

/* local error codes (-55000) */
#define OBS_FILE_SKIP_INPUT_LINE 		-55011
#define OBS_FILE_ARRIVALS_CROSS_YEAR_BOUNDARY 	-55022
#define OBS_FILE_INVALID_PHASE			-55033
#define OBS_FILE_INVALID_DATE			-55044
#define OBS_FILE_END_OF_EVENT			-55055
#define OBS_FILE_END_OF_INPUT			-55066
#define GRID_NOT_INSIDE				-55111

#define SMALLEST_EVENT_YEAR	1800
#define LARGEST_EVENT_YEAR	2100




/*------------------------------------------------------------*/
/* structures */

/* gaussian errors location paramters */
/*	see (TV82, eq. 10-14; MEN92, eq. 22) */
typedef struct
{
 	double SigmaT; 		/* theoretical error coeff for travel times */
 	double CorrLen; 	/* model corellation length */
 	DMatrix WtMtrx; 	/* weight matrix */
 	double WtMtrxSum; 	/* sum of elements of weight matrix */
 	long double meanObs; 	/* weighted mean of obs arrival times */
 	double meanPred; 	/* weighted mean of predicted travel times */
}
GaussLocParams;


/* scatter paramters */
typedef struct
{
 	int npts; 		/* number of scatter points */
}
ScatterParams;


/* station/inst/component parameters */
typedef struct
{
 	char label[ARRIVAL_LABEL_LEN]; 	/* label (i.e. station name) */
 	char inst[INST_LABEL_LEN]; 	/* instrument */
 	char comp[COMP_LABEL_LEN]; 	/* component */
	double amp_fact_ml_hb;		/* amplitude scale factor */
	double sta_corr_ml_hb;		/* station correction */
	double sta_corr_md_fmag;		/* station correction MD_FMAG */
}
CompDesc;


/* station label alias */
typedef struct
{
 	char name[ARRIVAL_LABEL_LEN]; 	/* original label (i.e. station name) */
 	char alias[ARRIVAL_LABEL_LEN]; 	/* alias, i.e. new label */
	int byr, bmo, bday;		/* begin year, month, day of validity */
	int eyr, emo, eday;		/* end year, month, day of validity */
}
AliasDesc;


/* excluded phase desc */
typedef struct
{
 	char label[ARRIVAL_LABEL_LEN]; 	/* label (i.e. station name) */
 	char phase[ARRIVAL_LABEL_LEN]; 	/* phase */
}
ExcludeDesc;


/* time delays */
typedef struct
{
 	char label[ARRIVAL_LABEL_LEN];
 	char phase[ARRIVAL_LABEL_LEN];
	int n_residuals;
	double delay;	/* time delay (sec) */
}
TimeDelayDesc;


/* Phase Identification */
typedef struct
{
	char phase[ARRIVAL_LABEL_LEN];
	char id_string[MAXLINE];
}
PhaseIdent;


/* Event time information extracted from filename or hypocenter line in obs file */
typedef struct
{
	int year, month, day;
	int hour, min;
	double sec;
}
EventTimeExtract;


/* Octtree */

#define OCTREE_UNDEF_VALUE -VERY_SMALL_DOUBLE

typedef struct
{
	int init_num_cells_x, init_num_cells_y, init_num_cells_z;
		// num nodes for each side of initial Tree3D
	double min_node_size;	// size of smallest side of smallest allowed node
	int max_num_nodes;	// maximum number of nodes to evaluate
	int num_scatter;	// number of scatter points to output
}
OcttreeParams;


/* Metropolis */
typedef struct
{
	double x, y, z;
	double dx;
	double likelihood;
}
WalkParams;


/* Magnitude */
typedef struct
{
	int type;		/* magnitude calculation method */

		/* ML - Hutton & Boore, BSSA, v77, n6, Dec 1987 */
	double amp_fact_ml_hb;	/* amplitude scale factor */
	double hb_n;	/* "n" in H&B eq. (2) */
	double hb_K;	/* "K" in H&B eq. (2) */
	double hb_Ro;	/* reference distance "100" in H&B eq. (2) */
	double hb_Mo;	/* reference magnitude "3.0" in H&B eq. (2) */

		/* coda duration (FMAG) - HYPOELLIPSE users manual chap 4;
			Lee et al., 1972; Lahr et al., 1975; Bakun and Lindh, 1977 */
	double fmag_c1, fmag_c2, fmag_c3, fmag_c4, fmag_c5;
}
MagDesc;



/*------------------------------------------------------------*/
/* globals  */

/* Gaussian error parameters */
EXTERN_TXT GaussLocParams Gauss;

/* Scatter parameters */
EXTERN_TXT ScatterParams Scatter;


/* events */
EXTERN_TXT int NumEvents;
EXTERN_TXT int NumEventsLocated;
EXTERN_TXT int NumLocationsCompleted;

#define MAX_NUM_OBS_FILES 1000
EXTERN_TXT int NumObsFiles;

/* number of arrivals used for location */
EXTERN_TXT int NumArrivalsLocation;

/* observations filenames */
char fn_loc_obs[MAX_NUM_OBS_FILES][FILENAME_MAX_SMALL];
/* filetype */
char ftype_obs[MAXLINE];

/* filenames */
char fn_loc_grids[FILENAME_MAX], fn_path_output[FILENAME_MAX];

/* location search type (grid, simulated annealing, Metropolis, etc) */
#define SEARCH_GRID  	0
#define SEARCH_MET  	1
#define SEARCH_OCTTREE  2
EXTERN_TXT int SearchType;

/* location method (misfit, etc) */
#define METH_UNDEF  		0
#define METH_GAU_ANALYTIC  	1
#define METH_GAU_TEST  		2
EXTERN_TXT int LocMethod;
EXTERN_TXT double DistStaGridMax;
EXTERN_TXT int MinNumArrLoc;
EXTERN_TXT int MaxNumArrLoc;
EXTERN_TXT int MinNumSArrLoc;
EXTERN_TXT double VpVsRatio;

/* location signature */
EXTERN_TXT char LocSignature[MAXLINE_LONG];

/* location grids */
#define MAX_NUM_LOCATION_GRIDS 10
EXTERN_TXT GridDesc  LocGrid[MAX_NUM_LOCATION_GRIDS];
EXTERN_TXT int NumLocGrids;
EXTERN_TXT int LocGridSave[MAX_NUM_LOCATION_GRIDS]; /* !should be in GridDesc */
EXTERN_TXT int Num3DGridReadToMemory, MaxNum3DGridMemory;

/* related hypocenter file pointers */
FILE *pSumFileHypNLLoc[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileHypo71[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileHypoEll[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileHypoInv[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileAlberto4[MAX_NUM_LOCATION_GRIDS];

/* related flags */
EXTERN_TXT int iWriteHypHeader[MAX_NUM_LOCATION_GRIDS];
/* hypocenter filetype saving flags */
EXTERN_TXT int iSaveNLLocEvent, iSaveNLLocSum,
		iSaveHypo71Event, iSaveHypo71Sum,
		iSaveHypoEllEvent, iSaveHypoEllSum,
		iSaveHypoInvSum, iSaveAlberto4Sum;


/* phase identification */
#define MAX_NUM_PHASE_ID 50
EXTERN_TXT PhaseIdent PhaseID[MAX_NUM_PHASE_ID];
EXTERN_TXT int NumPhaseID;

/* Extracted Filename Information */
EXTERN_TXT EventTimeExtract EventTime;

/* magnitude calculation */
#define MAG_UNDEF  	0
#define MAG_ML_HB  	1
#define MAG_MD_FMAG  	2
EXTERN_TXT int NumMagnitudeMethods;
#define MAX_NUM_MAG_METHODS  	2
EXTERN_TXT MagDesc Magnitude[MAX_NUM_MAG_METHODS];

/* station/inst/component parameters */
#define MAX_NUM_COMP_DESC 1000
EXTERN_TXT CompDesc Component[MAX_NUM_COMP_DESC];
EXTERN_TXT int NumCompDesc;

/* arrival label alias */
#define MAX_NUM_LOC_ALIAS 1000
#define MAX_NUM_LOC_ALIAS_CHECKS 2*MAX_NUM_LOC_ALIAS
EXTERN_TXT AliasDesc LocAlias[MAX_NUM_LOC_ALIAS];
EXTERN_TXT int NumLocAlias;

/* exclude arrivals */
#define MAX_NUM_LOC_EXCLUDE 1000
EXTERN_TXT ExcludeDesc LocExclude[MAX_NUM_LOC_EXCLUDE];
EXTERN_TXT int NumLocExclude;

/* station delays */
#define WRITE_RESIDUALS 0
#define WRITE_RES_DELAYS 1
#define WRITE_PDF_RESIDUALS 2
#define WRITE_PDF_DELAYS 3
#define MAX_NUM_STA_DELAYS 1000
EXTERN_TXT TimeDelayDesc TimeDelay[MAX_NUM_STA_DELAYS];
EXTERN_TXT int NumTimeDelays;

/* station list */
int NumStations;
SourceDesc StationList[MAX_NUM_ARRIVALS];

/* fixed origin time parameters */
EXTERN_TXT int FixOriginTimeFlag = 0;

/* Metropolis */
EXTERN_TXT WalkParams Metrop;	/* walk parameters */
EXTERN_TXT int MetNumSamples;	/* number of samples to evaluate */
EXTERN_TXT int MetLearn;	/* learning length in number of samples for
					calculation of sample statistics */
EXTERN_TXT int MetEquil;	/* number of samples to equil before using */
EXTERN_TXT int MetStartSave;	/* number of sample to begin saving */
EXTERN_TXT int MetSkip;		/* number of samples to wait between saves */
EXTERN_TXT double MetStepInit;	/* initial step size (km) (< 0.0 for auto) */
EXTERN_TXT double MetStepMin;	/* minimum step size (km) */
EXTERN_TXT double MetStepFact;	/* step size factor */
EXTERN_TXT double MetProbMin;	/* minimum likelihood necessary after learn */
EXTERN_TXT int MetUse;		/* number of samples to use
					= MetNumSamples - MetEquil */


/* Octtree */
EXTERN_TXT OcttreeParams octtreeParams;	/* Octtree parameters */
EXTERN_TXT Tree3D* octTree;		/* Octtree */
EXTERN_TXT ResultTreeNode* resultTreeRoot;	/* Octtree results tree root node*/


/* take-off angles */
int angleMode;		/* angle mode - ANGLE_MODE_NO, ANGLE_MODE_YES */
int iAngleQualityMin;	/* minimum quality for angles to be used */
#define ANGLE_MODE_NO	0
#define ANGLE_MODE_YES	1
#define ANGLE_MODE_UNDEF	-1



/*------------------------------------------------------------*/
/** hashtable routines for accumulating station statistics */

/* from Kernigham and Ritchie, C prog lang, 2nd ed, 1988, sec 6.6 */

struct staStatNode {	/* station statistics node for linked list */

	struct staStatNode *next;	/* next entry in list */
	char label[ARRIVAL_LABEL_LEN];  /* arrival label (station name) */
	char phase[ARRIVAL_LABEL_LEN];  /* arrival phase id */
	int flag_ignore;  		/* ignore flag  = 1 if phase not used for misfit calc */
		/* standard (max like point) residuals */
	double residual_min;		/* minimum residual */
	double residual_max;		/* maximum residual */
	double residual_sum;		/* accumulated residuals */
	double residual_square_sum;		/* accumulated square of residuals */
	double weight_sum;		/* accumulated weights */
	int num_residuals;		/* number of residuals accumulated */
		/* statistical (PDF weighted) residuals */
	double pdf_residual_sum;	/* accumulated pdf residuals */
	double pdf_residual_square_sum;	/* accumulated square of pdf residuals */
	int num_pdf_residuals;		/* number of residuals accumulated */
	double delay;			/* input time delay */

};
typedef struct staStatNode StaStatNode;

#define HASHSIZE 46
			/* table of pointers to list of StaStatNode */
EXTERN_TXT StaStatNode *hashtab[MAX_NUM_LOCATION_GRIDS][HASHSIZE];
			/* maxumum residual values to include in statistics */
EXTERN_TXT int NRdgs_Min;
EXTERN_TXT double RMS_Max, Gap_Max;
EXTERN_TXT double P_ResidualMax;
EXTERN_TXT double S_ResidualMax;

/* hashtable function declarations */
unsigned hash(char* , char* );
StaStatNode *lookup(int , char* , char* );
StaStatNode *InstallStaStatInTable(int, char* , char* , int , double ,
				double , double , double , double);
int WriteStaStatTable(int , FILE *, double , int , double ,
			double , double , int);
void UpdateStaStat(int , ArrivalDesc *, int , double , double );

/** end of hashtable routines */
/*------------------------------------------------------------*/



/*------------------------------------------------------------*/
/** grid memory management routines */

typedef struct gridMem {	/* 3D grid data in memory */

	GridDesc* pgrid;	/* pointer to copy of grid description structure */
	float* buffer;		/* corresponding buffer (contiguous floats) */
	float*** array;		/* corresponding array access to buffer */
	int grid_read;		/* gread read flag  = 1 if grid has been read from disk */
	int active;		/* active flag  = 1 if grid is being used in current location */


} GridMemStruct;

/* array of gridMem structure pointers for storing list of grids in memory */
GridMemStruct** GridMemList = NULL;
int GridMemListSize = 0;
int GridMemListNumElements = 0;

/* GridLib wrapper functions */
float* NLL_AllocateGrid(GridDesc* pgrid);
void NLL_FreeGrid(GridDesc* pgrid);
float*** NLL_CreateGridArray(GridDesc* pgrid);
void NLL_DestroyGridArray(GridDesc* pgrid);
int NLL_ReadGrid3dBuf(GridDesc* pgrid, FILE* fpio);
GridMemStruct* GridMemList_AddGridDesc(GridDesc* pgrid);
void GridMemList_AddElement(GridMemStruct* pnewGridMemStruct);
void GridMemList_RemoveElementAt(int index);
GridMemStruct* GridMemList_ElementAt(int index);
int GridMemList_IndexOfGridDesc(int verbose, GridDesc* pgrid);
int GridMemList_NumElements();


/** end of grid memory management routines */
/*------------------------------------------------------------*/



/*------------------------------------------------------------*/
/* function declarations */

int ReadNLLoc_Input(FILE* );
int GetNLLoc_Files(char* );
int GetNLLoc_Method(char* );
int GetNLLoc_SearchType(char* );
int GetNLLoc_FixOriginTime(char* );
int GetObservations(FILE* , char* , char* , ArrivalDesc*  , int* , int* , int, HypoDesc* , int* , int* );
int GetNextObs(FILE* , ArrivalDesc *, char* , int );
void removeSpace(char *str);
int IsPhaseID(char *, char *);
int IsGoodDate(int , int , int );
int EvalPhaseID(char *);
int ReadArrivalSheets(int , ArrivalDesc* , double );
int IsDuplicateArrival(ArrivalDesc *, int , int , char *);

int WriteHypo71(FILE *, HypoDesc* , ArrivalDesc* , char* , int , int );
int WriteHypoEll(FILE *, HypoDesc* , ArrivalDesc* , char* , int , int );
int WriteHypoInverseArchive(FILE *, HypoDesc* , ArrivalDesc* , int ,
		char* , int , int );
int WriteHypoAlberto4(FILE *, HypoDesc* , ArrivalDesc* , char* );
int OpenSummaryFiles(char *);
int CloseSummaryFiles();

int GetCompDesc(char* );
int GetLocAlias(char* );
int GetLocExclude(char* line1);
int GetTimeDelays(char* );

int LocGridSearch(int , int , int , ArrivalDesc* ,  GridDesc* ,
		GaussLocParams* , HypoDesc* );
int LocMetropolis(int , int , int , ArrivalDesc *,
	GridDesc* , GaussLocParams* , HypoDesc* , WalkParams* , float* );
int SaveBestLocation(int , int , ArrivalDesc *,  GridDesc* , GaussLocParams* , HypoDesc* ,
		double , int );
int ConstWeightMatrix(int , ArrivalDesc* , GaussLocParams* );
void CalcCenteredTimesObs(int , ArrivalDesc* , GaussLocParams* , HypoDesc* );
inline void CalcCenteredTimesPred(int , ArrivalDesc* , GaussLocParams* );
inline double CalcSolutionQuality(int , ArrivalDesc* , GaussLocParams* , int, double* );
inline double CalcSolutionQuality_GAU_ANALYTIC(int , ArrivalDesc* , GaussLocParams* , int, double* );
inline double CalcSolutionQuality_GAU_TEST(int , ArrivalDesc* , GaussLocParams* , int, double* );
long double CalcMaxLikeOriginTime(int , ArrivalDesc* , GaussLocParams* );
inline void UpdateProbabilisticResiduals(int , ArrivalDesc *, double);
int CalcConfidenceIntrvl(GridDesc* , HypoDesc* , char* );
int HomogDateTime(ArrivalDesc* , int , HypoDesc* );
int StdDateTime(ArrivalDesc* , int , HypoDesc* );
int SetOutName(ArrivalDesc* , char* , char* , int );
int Locate(int , char *);
int SaveLocation(int , char *);
int GenEventScatterGrid(GridDesc* , HypoDesc* , ScatterParams* , char* );
void InitializeArrivalFields(ArrivalDesc *);
int isExcluded(char *label, char *phase);
int EvaluateArrivalAlias(ArrivalDesc *);
int ApplyTimeDelays(ArrivalDesc *);

double CalcAzimuthGap(ArrivalDesc *, int );

void InitializeMetropolisWalk(GridDesc* ,  ArrivalDesc* , int ,
		WalkParams* , int , double );
inline int GetNextMetropolisSample(WalkParams* , double , double ,
	double , double , double , double , double* , double* , double* );
inline int MetropolisTest(WalkParams* , double );

double CalculateVpVsEstimate(HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals);

int CalculateMagnitude(HypoDesc* , ArrivalDesc* , int , CompDesc* , int, MagDesc* );
int findStaInstComp(ArrivalDesc* parr, CompDesc* pcomp, int nCompDesc);
double Calc_ML_HuttonBoore(double amplitude, double dist, double depth, double sta_corr,
			double hb_n, double hb_K, double hb_Ro, double hb_Mo);
double Calc_MD_FMAG(double coda_dur, double dist, double depth, double sta_corr,
	double fmag_c1, double fmag_c2, double fmag_c3, double fmag_c4, double fmag_c5);

int addToStationList(SourceDesc *stations, int numStations, ArrivalDesc *arrival, int nArrivals);

int getTravelTimes(ArrivalDesc *arrival, int num_arr_loc, double xval, double yval, double zval);

Tree3D*  InitializeOcttree(GridDesc* ptgrid, ArrivalDesc* parrivals,
		int numArrLoc,  OcttreeParams* pParams);
int LocOctree(int ngrid, int num_arr_total, int num_arr_loc,
		ArrivalDesc *arrival,
		GridDesc* ptgrid, GaussLocParams* gauss_par, HypoDesc* phypo,
		OcttreeParams* pParams, Tree3D* pOctTree, float* fdata, double* poct_node_value_max);
int GenEventScatterOcttree(OcttreeParams* pParams, double oct_node_value_max, float* fscatterdata);
double integrateResultTree(ResultTreeNode* prtree, double sum, double oct_node_value_max);
int getScatterSampleResultTree(ResultTreeNode* prtree, OcttreeParams* pParams,
		double integral, float* fdata, int npoints, int* pfdata_index,
		double oct_node_value_max);



/*** program to generate synthetic travel times to stations */

#define NARGS 2

main(int argc, char *argv[])
{

	int istat, n;
	int i_end_of_input, iLocated;
	int narr, nsta, ngrid, nObsFile;
	int numArrivalsIgnore, numSArrivalsLocation;
	int maxArrExceeded = 0;
	char fn_root_out[FILENAME_MAX], fname[FILENAME_MAX];
	FILE *fp_obs, *fpio;



	/* set program name */
	strcpy(prog_name, PNAME);

	/* check command line for correct usage */

	if (argc != NARGS) {
		disp_usage(prog_name, "<control file>");
		exit(EXIT_ERROR_USAGE);
	}



	/* set constants */

	SetConstants();
	NumLocGrids = 0;
	NumEvents = NumEventsLocated = NumLocationsCompleted = 0;
	NumCompDesc = 0;
	NumLocAlias = 0;
	NumLocExclude = 0;
	NumTimeDelays = 0;
	NumPhaseID = 0;
	DistStaGridMax = 0.0;
	MinNumArrLoc = 0;
	MinNumSArrLoc = 0;
	MaxNumArrLoc = MAX_NUM_ARRIVALS;
	FixOriginTimeFlag = 0;
	Scatter.npts = -1;
	for (n = 0; n < MAX_NUM_MAG_METHODS; n++)
		Magnitude[n].type = MAG_UNDEF;
	NumMagnitudeMethods = 0;
	MaxNum3DGridMemory = -1;

	/* read control file */

	strcpy(fn_control, argv[1]);
	if ((fp_control = fopen(fn_control, "r")) == NULL) {
		puterr("FATAL ERROR: opening control file.");
		exit(EXIT_ERROR_FILEIO);
	} else {
		NumFilesOpen++;
	}


	/* read control file */

	if ((istat = ReadNLLoc_Input(fp_control)) < 0) {
		puterr("FATAL ERROR: reading control file.");
		exit(EXIT_ERROR_FILEIO);
	}
	fclose(fp_control);
	NumFilesOpen--;


	/* initialize random number generator */

	SRAND_FUNC(RandomNumSeed);
	if (message_flag >= 4)
		test_rand_int();


	/* open summary output file */

	if ((istat = OpenSummaryFiles(fn_path_output)) < 0) {
		puterr("FATAL ERROR: opening hypocenter summary files.");
		exit(EXIT_ERROR_FILEIO);
	}


	/* perform location for each observation file */

	for (nObsFile = 0; nObsFile < NumObsFiles; nObsFile++) {

	    i_end_of_input = 0;

	    putmsg(2, "");
	    sprintf(MsgStr, "... Reading observation file %s",
		fn_loc_obs[nObsFile]);
	    putmsg(1, MsgStr);

	    /* open observation file */

    	    if ((fp_obs = fopen(fn_loc_obs[nObsFile], "r")) == NULL) {
		puterr2("ERROR: opening observations file",
			fn_loc_obs[nObsFile]);
		continue;
	    } else {
		NumFilesOpen++;
	    }

	    /* extract info from filename */
	    if ((istat = ExtractFilenameInfo(fn_loc_obs[nObsFile], ftype_obs))
			< 0)
		puterr("WARNING: error extractng information from filename.");


	    /* read arrivals and locate event for each  */
	    /*		event (set of observations) in file */

	    NumArrivals = 0;
	    while (1)
	    {

		iLocated = 0;

		if (i_end_of_input)
			break;

		if (NumArrivals != OBS_FILE_SKIP_INPUT_LINE) {
			putmsg(2, "");
	    		sprintf(MsgStr,
"Reading next set of observations (Files open: Tot:%d Buf:%d Hdr:%d  Alloc: %d) ...",
				NumFilesOpen, NumGridBufFilesOpen, NumGridHdrFilesOpen, NumAllocations);
	   		 putmsg(1, MsgStr);
		}


		/* read next set of observations */

		NumArrivalsLocation = 0;
		if ((NumArrivals = GetObservations(fp_obs,
				ftype_obs, fn_loc_grids, Arrival,
				&i_end_of_input, &numArrivalsIgnore,
				MaxNumArrLoc, &Hypocenter,
				&maxArrExceeded, &numSArrivalsLocation)) == 0)
			break;

		if (NumArrivals < 0)
			goto cleanup;


		/* set number of arrivals to be used in location */

		NumArrivalsLocation = NumArrivals - numArrivalsIgnore;

		putmsg(2, "");
		SetOutName(Arrival + 0, fn_path_output, fn_root_out, 1);
		sprintf(MsgStr,
"... %d observations read, %d will be used for location (%s).",
			NumArrivals, NumArrivalsLocation, fn_root_out);
		putmsg(1, MsgStr);


		/* check for minimum number of arrivals */

		if (NumArrivalsLocation < MinNumArrLoc) {
			sprintf(MsgStr,
"WARNING: too few observations to locate, skipping event.");
			putmsg(1, MsgStr);
			sprintf(MsgStr,
"INFO: %d observations needed (specified in control file entry LOCMETH).",
			MinNumArrLoc);
			putmsg(2, MsgStr);
			goto cleanup;
		}


		/* check for minimum number of S arrivals */

		if (numSArrivalsLocation < MinNumSArrLoc) {
			sprintf(MsgStr,
"WARNING: too few S observations to locate, skipping event.");
			putmsg(1, MsgStr);
			sprintf(MsgStr,
"INFO: %d S observations needed (specified in control file entry LOCMETH).",
			MinNumSArrLoc);
			putmsg(2, MsgStr);
			goto cleanup;
		}


		/* process arrivals */

		/* add stations to station list */

		if (iSaveNLLocSum)
			NumStations += addToStationList(StationList, NumStations, Arrival, NumArrivals);

		/* sort to get location arrivals in time order */

		if ((istat =
			SortArrivalsIgnore(Arrival, NumArrivals)) < 0) {
				puterr(
				"ERROR: sorting arrivals by ignore flag.");
			goto cleanup;
		}
		if ((istat =
			SortArrivalsTime(Arrival, NumArrivalsLocation)) < 0) {
				puterr("ERROR: sorting arrivals by time.");
			goto cleanup;
		}


		/* construct weight matrix (TV82, eq. 10-9; MEN92, eq. 12) */

		if ((istat = ConstWeightMatrix(NumArrivalsLocation, Arrival,
				&Gauss)) < 0) {
			puterr("ERROR: constructing weight matrix - NLLoc requires non-zero observation or modelisation errors.");
			/* close time grid files and continue */
			goto cleanup;
		}


		/* calculate weighted mean of obs arrival times   */
		/*	(TV82, eq. A-38) */

		CalcCenteredTimesObs(NumArrivalsLocation, Arrival, &Gauss,
						&Hypocenter);


		/* preform location for each grid */

	    	sprintf(MsgStr,
"Locating... (Files open: Tot:%d Buf:%d Hdr:%d  Alloc: %d  3DMem: %d) ...",
		NumFilesOpen, NumGridBufFilesOpen, NumGridHdrFilesOpen, NumAllocations, Num3DGridReadToMemory);
	   		 putmsg(1, MsgStr);

		for (ngrid = 0; ngrid < NumLocGrids; ngrid++)
			if ((istat = Locate(ngrid, fn_root_out)) < 0) {
				if (istat == GRID_NOT_INSIDE)
					break;
				else {
					puterr(
					  "ERROR: location failed.");
					goto cleanup;
				}
			}

		NumEventsLocated++;
		if (istat == 0 && ngrid == NumLocGrids)
			NumLocationsCompleted++;
		iLocated = 1;


		cleanup: ;

		NumEvents++;


		/* release grid buffer or sheet storage */

		for (narr = 0; narr < NumArrivalsLocation; narr++) {
			DestroyGridArray(&(Arrival[narr].sheetdesc));
			FreeGrid(&(Arrival[narr].sheetdesc));
			NLL_DestroyGridArray(&(Arrival[narr].gdesc));
			NLL_FreeGrid(&(Arrival[narr].gdesc));
		}

		/* close time grid files (opened in function GetObservations) */

		for (narr = 0; narr < NumArrivalsLocation; narr++)
			CloseGrid3dFile(&(Arrival[narr].fpgrid),
				&(Arrival[narr].fphdr));

		if (iLocated) {
			putmsg(2, "");
	    		sprintf(MsgStr,
				"Finished event location, output files: %s.*",
				fn_root_out);
			putmsg(1, MsgStr);
		} else
			putmsg(2, "");

	    }  /* next event */

	    putmsg(2, "");
	    sprintf(MsgStr, "...end of observation file detected.");
	    putmsg(1, MsgStr);
	    fclose(fp_obs);
	    NumFilesOpen--;

	}  /* next observation file */

	putmsg(2, "");
	sprintf(MsgStr,
"No more observation files.  %d events read,  %d events located,  %d locations completed.",
		NumEvents, NumEventsLocated, NumLocationsCompleted);
	putmsg(1, MsgStr);
	putmsg(2, "");


	/* write cumulative arrival statistics */
	for (ngrid = 0; ngrid < NumLocGrids; ngrid++)
		if (LocGridSave[ngrid]) {
			sprintf(fname, "%s.sum.grid%d.loc.stat",
						fn_path_output, ngrid);
			if ((fpio = fopen(fname, "w")) == NULL) {
				puterr2(
"ERROR: opening cumulative phase statistics output file", fname);
				return(-1);
			} else {
				NumFilesOpen++;
			}
			WriteStaStatTable(ngrid, fpio,
				RMS_Max, NRdgs_Min, Gap_Max,
				P_ResidualMax, S_ResidualMax, WRITE_RESIDUALS);
			WriteStaStatTable(ngrid, fpio,
				RMS_Max, NRdgs_Min, Gap_Max,
				P_ResidualMax, S_ResidualMax, WRITE_RES_DELAYS);
			WriteStaStatTable(ngrid, fpio,
				RMS_Max, NRdgs_Min, Gap_Max,
				P_ResidualMax, S_ResidualMax,
					WRITE_PDF_RESIDUALS);
			WriteStaStatTable(ngrid, fpio,
				RMS_Max, NRdgs_Min, Gap_Max,
				P_ResidualMax, S_ResidualMax, WRITE_PDF_DELAYS);
			close(fpio);
		}

	/* write station list */
	for (ngrid = 0; ngrid < NumLocGrids; ngrid++)
		if (LocGridSave[ngrid]) {
			sprintf(fname, "%s.sum.grid%d.loc.stations",
						fn_path_output, ngrid);
			if ((fpio = fopen(fname, "w")) == NULL) {
				puterr2(
"ERROR: opening cumulative phase statistics output file", fname);
				return(-1);
			} else {
				NumFilesOpen++;
			}
			WriteStationList(fpio, StationList, NumStations);
			close(fpio);
		}

	CloseSummaryFiles();

	exit(EXIT_NORMAL);

}



/*** function to perform grid search location */

int Locate(int ngrid, char* fn_root_out)
{

	int istat, n, narr;
	char fnout[FILENAME_MAX];

	FILE *fpio;
	char fname[FILENAME_MAX];
	int nScatterSaved = -1;
	float *fdata, ftemp;
	int iSizeOfFdata;
	double alpha_2;
	double oct_node_value_max;



	/* write message */

	putmsg(2, "");
	if (SearchType == SEARCH_GRID)
		sprintf(MsgStr, "Searching Grid %d:", ngrid);
	else if (SearchType == SEARCH_MET)
		sprintf(MsgStr, "Applying Metropolis within Grid %d:", ngrid);
	else if (SearchType == SEARCH_OCTTREE)
		sprintf(MsgStr, "Applying Octtree search within Grid %d:", ngrid);
	putmsg(2, MsgStr);
	if (message_flag >= 3)
		display_grid_param(LocGrid + ngrid);



	/* set output name */
	sprintf(fnout, "%s.grid%d", fn_root_out, ngrid);
	strcpy(Hypocenter.fileroot, fnout);

	/* initialize hypocenter fields */
	sprintf(Hypocenter.locStat, "LOCATED");
	sprintf(Hypocenter.locStatComm, "Location completed.");
	Hypocenter.x = Hypocenter.y = Hypocenter.z = 0.0;
	Hypocenter.ix = Hypocenter.iy = Hypocenter.iz = -1;


	/* search type dependent initializations */

	if (SearchType == SEARCH_GRID) {

		/* check that current grid is contained within first grid */

		if (!IsGridInside(LocGrid + ngrid, LocGrid, 0)) {
			puterr(
	"WARNING: this grid not entirely contained inside 0th grid, ending search for this event.");
			return(GRID_NOT_INSIDE);
		}

		/* initialize 3D location grid */

		/* allocate location grid */
		LocGrid[ngrid].buffer = AllocateGrid(LocGrid + ngrid);
		if (LocGrid[ngrid].buffer == NULL) {
			puterr(
	"ERROR: allocating memory for 3D location grid buffer.");
			return(EXIT_ERROR_MEMORY);
		}
		/* create array access pointers */
		LocGrid[ngrid].array = CreateGridArray(LocGrid + ngrid);
		if (LocGrid[ngrid].array == NULL) {
			puterr(
	"ERROR: creating array for accessing 3D location grid buffer.");
			return(EXIT_ERROR_MEMORY);
		}
		LocGrid[ngrid].sum = 0.0;


		/* reset y-z dual-sheet grids (3D time grids) */

		for (narr = 0; narr < NumArrivalsLocation; narr++) {
			if (Arrival[narr].sheetdesc.type == GRID_TIME)
				Arrival[narr].sheetdesc.origx =
						VERY_LARGE_DOUBLE;
		}



	} else if (SearchType == SEARCH_MET) {

/* test change 17JAN2000 AJL */
/*		InitializeMetropolisWalk(LocGrid + ngrid ,
			Arrival, NumArrivalsLocation, &Metrop,
			MetNumSamples, MetStepInit);
*/

		InitializeMetropolisWalk(LocGrid + ngrid ,
			Arrival, NumArrivalsLocation, &Metrop,
			MetLearn + MetEquil, MetStepInit);

		/* allocate scatter array for saved samples */
		iSizeOfFdata = (1 + MetUse / MetSkip) * 4 * sizeof(float);
		if ((fdata = (float *)
				malloc(iSizeOfFdata)) == NULL)
			return(EXIT_ERROR_LOCATE);
		NumAllocations++;

	} else if (SearchType == SEARCH_OCTTREE) {

		octTree = InitializeOcttree(LocGrid + ngrid ,
			Arrival, NumArrivalsLocation, &octtreeParams);
		NumAllocations++;

		/* allocate scatter array for saved samples */
		iSizeOfFdata = octtreeParams.num_scatter * 4 * sizeof(float);
		iSizeOfFdata = (12 * iSizeOfFdata) / 10; // sample may be slightly larger than requested
		if ((fdata = (float *) malloc(iSizeOfFdata)) == NULL)
			return(EXIT_ERROR_LOCATE);
		NumAllocations++;

	}


	/* since sorted, reset companion indices */
	if (VpVsRatio > 0.0) {
		//for (narr = 0; narr < NumArrivalsLocation; narr++) {
		for (narr = 0; narr < NumArrivals; narr++) {
			if (Arrival[narr].n_companion < 0)
				continue;
			if (IsPhaseID(Arrival[narr].phase, "S") &&
					(Arrival[narr].n_companion =
					IsDuplicateArrival(Arrival, narr, narr, "P")) < 0) {
				puterr("ERROR: cannot find companion arrival.");
				return(EXIT_ERROR_LOCATE);
			}
		}
	}


	/* do search */

	if (SearchType == SEARCH_GRID) {

		/* grid-search location (fill location grid) */
		if ((istat =
			LocGridSearch(ngrid, NumArrivals, NumArrivalsLocation,
				Arrival, LocGrid + ngrid,
				&Gauss, &Hypocenter)) < 0) {
			puterr("ERROR: in grid search location.");
			return(EXIT_ERROR_LOCATE);
		}

	} else if (SearchType == SEARCH_MET) {

		/* Metropolis location (random walk) */
		if ((nScatterSaved =
			LocMetropolis(ngrid, NumArrivals, NumArrivalsLocation,
				Arrival, LocGrid + ngrid,
				&Gauss, &Hypocenter, &Metrop, fdata)) < 0) {
			puterr("ERROR: in Metropolis location.");
			return(EXIT_ERROR_LOCATE);
		}

	} else if (SearchType == SEARCH_OCTTREE) {

		/* Octree location (importance sampling) */
		if ((nScatterSaved =
			LocOctree(ngrid, NumArrivals, NumArrivalsLocation,
				Arrival, LocGrid + ngrid,
				&Gauss, &Hypocenter, &octtreeParams,
				octTree, fdata, &oct_node_value_max)) < 0) {
			puterr("ERROR: in Octree location.");
			return(EXIT_ERROR_LOCATE);
		}

	}



	/* clean up dates */
	StdDateTime(Arrival, NumArrivals, &Hypocenter);

	/* determine azimuth gap */
	Hypocenter.gap = CalcAzimuthGap(Arrival, NumArrivalsLocation);

	/* re-sort arrivals by distance */
	if ((istat = SortArrivalsDist(Arrival, NumArrivals)) < 0) {
		puterr("ERROR: sorting arrivals by distance.");
		return(EXIT_ERROR_LOCATE);
	}

	/* save distance to closest station */
	Hypocenter.dist = Arrival->dist;


	/* search type dependent processing */

	alpha_2 = 3.53;  /* value for 68% conf
				(see Num Rec, 2nd ed, sec 15.6) */
	if (SearchType == SEARCH_GRID && LocGridSave[ngrid]) {

		/* calculate confidence intervals and save to disk */

		if (LocGrid[ngrid].type == GRID_PROB_DENSITY)  {
			if ((istat = CalcConfidenceIntrvl(LocGrid + ngrid,
					&Hypocenter, fnout)) < 0) {
				puterr("ERROR: calculating confidence intervals.");
				return(EXIT_ERROR_LOCATE);
			}

			/* generate probabilistic scatter of events */

			if ((istat = GenEventScatterGrid(LocGrid + ngrid,
					&Hypocenter,
					&Scatter, fnout)) < 0) {
				puterr("ERROR: calculating event scatter.");
			}

			/* calculate "traditional" statistics */
			Hypocenter.expect = CalcExpectation(
						LocGrid + ngrid, NULL);
			istat = rect2latlon(0, Hypocenter.expect.x, Hypocenter.expect.y,
				&(Hypocenter.expect_dlat), &(Hypocenter.expect_dlong));
			Hypocenter.cov = CalcCovariance(LocGrid + ngrid,
				&Hypocenter.expect, NULL);
			Hypocenter.ellipsoid = CalcErrorEllipsoid(
					&Hypocenter.cov, alpha_2);

		} else {
			Hypocenter.probmax = -1.0;
		}

	} else if ((SearchType == SEARCH_MET || SearchType == SEARCH_OCTTREE ) && LocGridSave[ngrid]) {

		if (SearchType == SEARCH_OCTTREE && nScatterSaved == 0) // not saved during search
			nScatterSaved = GenEventScatterOcttree(
				&octtreeParams, oct_node_value_max, fdata);

		/* write scatter file */
		sprintf(fname, "%s.loc.scat", fnout);
		if ((fpio = fopen(fname, "w")) != NULL) {
			/* write scatter file header informaion */
			fseek(fpio, 0, SEEK_SET);
			fwrite(&nScatterSaved, sizeof(int), 1, fpio);
			ftemp = (float) Hypocenter.probmax;
			fwrite(&ftemp, sizeof(float), 1, fpio);
			/* skip header record */
			fseek(fpio, 4 * sizeof(float), SEEK_SET);
			/* write scatter samples */
			fwrite(fdata, 4 * sizeof(float), nScatterSaved, fpio);
			fclose(fpio);
		} else {
			puterr("ERROR: opening scatter output file.");
			return(EXIT_ERROR_IO);
		}

		/* calculate "traditional" statistics */
		Hypocenter.expect =
			CalcExpectationSamples(fdata, nScatterSaved);
		istat = rect2latlon(0, Hypocenter.expect.x, Hypocenter.expect.y,
			&(Hypocenter.expect_dlat), &(Hypocenter.expect_dlong));
		Hypocenter.cov = CalcCovarianceSamples(fdata, nScatterSaved,
			&Hypocenter.expect);
		if (nScatterSaved) {
			Hypocenter.ellipsoid = CalcErrorEllipsoid(
					&Hypocenter.cov, alpha_2);
		}
	}





	/* search type dependent results saving */

	if (SearchType == SEARCH_GRID) {

		/* save location grid to disk */

		if (LocGridSave[ngrid])
			if ((istat = WriteGrid3dBuf(LocGrid + ngrid, NULL,
					fnout, "loc")) < 0) {
				puterr("ERROR: writing location grid to disk.");
				return(EXIT_ERROR_IO);
			}
	} else if (SearchType == SEARCH_MET || SearchType == SEARCH_OCTTREE) {

		/* save location grid header to disk */

		if (LocGridSave[ngrid])
			if ((istat = WriteGrid3dHdr(LocGrid + ngrid, NULL,
					fnout, "loc")) < 0) {
				puterr(
				"ERROR: writing grid header to disk.");
				return(EXIT_ERROR_IO);
			}
	}



	/* display and save minimum misfit location to file */

	if (LocGridSave[ngrid]) {
		/* calculate magnitudes */
		Hypocenter.amp_mag = MAG_NULL;
		Hypocenter.num_amp_mag = 0;
		Hypocenter.dur_mag = MAG_NULL;
		Hypocenter.num_dur_mag = 0;
		for (n = 0; n < MAX_NUM_MAG_METHODS; n++)
			CalculateMagnitude(&Hypocenter, Arrival, NumArrivals,
				Component, NumCompDesc, Magnitude + n);
		/* calculate estimated VpVs ratio */
		CalculateVpVsEstimate(&Hypocenter, Arrival, NumArrivals);
		/* save location */
		if ((istat = SaveLocation(ngrid, fnout)) < 0)
			return(istat);
		/* update station statistics table */
		if (strncmp(Hypocenter.locStat, "LOCATED", 7) == 0
				&& Hypocenter.rms <= RMS_Max
				&& Hypocenter.nreadings >= NRdgs_Min
				&& Hypocenter.gap <= Gap_Max)
			UpdateStaStat(ngrid, Arrival, NumArrivals,
				P_ResidualMax, S_ResidualMax);
	}



	/* search type dependent cleanup */

	if (SearchType == SEARCH_GRID) {

		/* free grid memory */

		DestroyGridArray(LocGrid + ngrid);
		FreeGrid(LocGrid + ngrid);

		/* intialize next grid origin location */

		if (ngrid < NumLocGrids - 1) {
			if (LocGrid[ngrid + 1].autox)
				LocGrid[ngrid + 1].origx = Hypocenter.x
				- 0.5 * (double) (LocGrid[ngrid + 1].numx - 1)
						* LocGrid[ngrid + 1].dx;
			if (LocGrid[ngrid + 1].autoy)
				LocGrid[ngrid + 1].origy = Hypocenter.y
				- 0.5 * (double) (LocGrid[ngrid + 1].numy - 1)
						* LocGrid[ngrid + 1].dy;
			if (LocGrid[ngrid + 1].autoz)
				LocGrid[ngrid + 1].origz = Hypocenter.z
				- 0.5 * (double) (LocGrid[ngrid + 1].numz - 1)
						* LocGrid[ngrid + 1].dz;

			/* try to make sure new grid is inside initial grid */
			if (!IsGridInside(LocGrid + ngrid + 1, LocGrid, 1))
				puterr(
"WARNING: cannot get next grid entirely contained inside 0th grid.");
		}

	} else if (SearchType == SEARCH_MET) {

		/* free saved samples memory */
		free(fdata);
		NumAllocations--;

	} else if (SearchType == SEARCH_OCTTREE) {

		// free results tree - IMPORTANT!
		freeResultTree(resultTreeRoot);

		/* free octree memory */
		freeTree3D(octTree);
		NumAllocations--;

		/* free saved samples memory */
		free(fdata);
		NumAllocations--;

	}



	/* re-sort to get location arrivals in time order */

	if ((istat =
		SortArrivalsIgnore(Arrival, NumArrivals)) < 0) {
			puterr(
			"ERROR: sorting arrivals by ignore flag.");
		return(EXIT_ERROR_LOCATE);
	}
	if ((istat = SortArrivalsTime(Arrival, NumArrivalsLocation)) < 0) {
		puterr("ERROR: sorting arrivals by time.");
		return(EXIT_ERROR_LOCATE);
	}



	/* search type dependent return */

	if (SearchType == SEARCH_GRID) {

		return(0);

	} else if (SearchType == SEARCH_MET || SearchType == SEARCH_OCTTREE) {

		if (nScatterSaved == 0)
			return(1);

		return(0);
	}



}



/** function to initialize Metropolis walk */

void InitializeMetropolisWalk(GridDesc* ptgrid, ArrivalDesc* parrivals, int
	numArrLoc, WalkParams* pMetrop, int numSamples, double initStep)
{
	int narr;
	double xlen, ylen, zlen, dminlen;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	SourceDesc* pstation;


	/* set walk limits equal to grid limits */
	xmin = ptgrid->origx;
	xmax = xmin + (double) (ptgrid->numx - 1) * ptgrid->dx;
	ymin = ptgrid->origy;
	ymax = ymin + (double) (ptgrid->numy - 1) * ptgrid->dy;

	/* find station with earliest arrival and non-zero weight */
	narr = 0;
	while(narr < numArrLoc && parrivals[narr].weight < 0.001)
		narr++;

	/* initialize walk location */
	if (narr < numArrLoc)
		pstation = &(parrivals[narr].station);
	if (narr < numArrLoc &&
			pstation->x >= xmin && pstation->x <= xmax
			&& pstation->y >= ymin && pstation->y <= ymax) {
		/* start walk at location of station with earliest arrival */
		pMetrop->x = pstation->x;
		pMetrop->y = pstation->y;
	} else {
		/* start walk at grid center */
		pMetrop->x = ptgrid->origx
			+ (double) (ptgrid->numx - 1) * ptgrid->dx / 2.0;
		pMetrop->y = ptgrid->origy
			+ (double) (ptgrid->numy - 1) * ptgrid->dy / 2.0;
	}
	/* start walk at grid center depth */
	pMetrop->z = ptgrid->origz
		+ (double) (ptgrid->numz - 1) * ptgrid->dz / 2.0;

	/* calculate initial step size */
	if (initStep < 0.0) {
		xlen = (double) ptgrid->numx * ptgrid->dx / 2.0;
		ylen = (double) ptgrid->numy * ptgrid->dy / 2.0;
		zlen = (double) ptgrid->numz * ptgrid->dz / 2.0;
		dminlen = xlen < ylen ? xlen : ylen;
		dminlen = dminlen < zlen ? dminlen : zlen;
		/* step is size that tiles plane parallel to max len sides */
		pMetrop->dx = sqrt((xlen * ylen * zlen / dminlen) / (double) numSamples);
		/* step is size that tiles search volume */
		/*pMetrop->dx = pow(xlen * ylen * zlen / (double) numSamples, 1.0/3.0);*/
	} else {
		pMetrop->dx = initStep;
	}

	sprintf(MsgStr,
		"INFO: Metropolis initial step size: %lf", pMetrop->dx);
	putmsg(3, MsgStr);

	/* set likelihood */
	pMetrop->likelihood = -1.0;

}



/** function to display and save minimum misfit location to file */

int SaveLocation(int ngrid, char *fnout)
{
	int istat;
	char *pchr;
	char sys_command[MAXLINE_LONG];
	char fname[FILENAME_MAX], frootname[FILENAME_MAX], fpathname[FILENAME_MAX];
	char *ppath;
	FILE *fp_tmp;

		/* set signature string */
		sprintf(Hypocenter.signature, "%s   %s:v%s %s",
			LocSignature, PNAME, PVER, CurrTimeStr());
		while(pchr = strchr(Hypocenter.signature, '\n'))
			*pchr = ' ';

		/* display hypocenter to std out */
		if (message_flag >= 3)
			WriteLocation(stdout, &Hypocenter, Arrival,
				NumArrivals, fnout, 1, 1, 0, LocGrid + ngrid, 0);

		// get path to output files
		strcpy(fpathname, fnout);
		if ((ppath = strrchr(fpathname, '/')) != NULL
				|| (ppath = strrchr(fpathname, '\\')) != NULL)
			*(ppath + 1) = '\0';
		else
			strcpy(fpathname, "");

		/*  save requested hypocenter/phase formats */

		if (iSaveNLLocEvent) {
			/* write NLLoc hypocenter to event file */
			sprintf(frootname, "%s.loc", fnout);
			sprintf(fname, "%s.hyp", frootname);
			if ((istat = WriteLocation(NULL, &Hypocenter, Arrival,
					NumArrivals, fname, 1, 1, 0,
					LocGrid + ngrid, 0)) < 0) {
				puterr(
"ERROR: writing location to event file.");
				return(EXIT_ERROR_IO);
			}
			/* copy event file to last.hyp */
			sprintf(sys_command,
				"cp %s %slast.hyp", fname, fpathname);
			system(sys_command);
			sprintf(fname, "%s.hdr", frootname);
			sprintf(sys_command,
				"cp %s %slast.hdr", fname, fpathname);
			system(sys_command);
			sprintf(fname, "%s.scat", frootname);
			sprintf(sys_command,
				"cp %s %slast.scat", fname, fpathname);
			system(sys_command);
		}
		if (iSaveNLLocSum) {
			/* write NLLoc hypocenter to summary file */
			if ((istat = WriteLocation(pSumFileHypNLLoc[ngrid],
					&Hypocenter,
					Arrival, NumArrivals, fnout, 0, 1, 0,
					LocGrid + ngrid, 0)) < 0) {
				puterr(
"ERROR: writing location to summary file.");
				return(EXIT_ERROR_IO);
			}
			/* copy event grid header to .sum header */
			sprintf(sys_command,
				"cp %s.loc.hdr %s.sum.grid%d.loc.hdr",
				fnout, fn_path_output, ngrid);
			system(sys_command);
		}

		if (iSaveHypo71Event) {
			/* write HYPO71 hypocenter to event file */
			WriteHypo71(NULL, &Hypocenter, Arrival, fnout, 1, 1);
		}
		if (iSaveHypo71Sum) {
			/* write HYPO71 hypocenter to summary file */
			WriteHypo71(pSumFileHypo71[ngrid], &Hypocenter,
				Arrival, fnout, iWriteHypHeader[ngrid], 0);
		}

		if (iSaveHypoEllEvent) {
			/* write pseudo-HypoEllipse hypo to event file */
			WriteHypoEll(NULL, &Hypocenter, Arrival, fnout, 1, 1);
		}
		if (iSaveHypoEllSum) {
			/* write pseudo-HypoEllipse hypo to summary file */
			WriteHypoEll(pSumFileHypoEll[ngrid], &Hypocenter,
				Arrival,
				fnout, iWriteHypHeader[ngrid], 0);
		}

		if (iSaveHypoInvSum) {
			/* write HypoInverseArchive (FPFIT) hypocenter to summary file */
			iWriteHypHeader[ngrid] = 0;
			WriteHypoInverseArchive(pSumFileHypoInv[ngrid], &Hypocenter,
				Arrival, NumArrivals, fnout, iWriteHypHeader[ngrid], 1);
			/* also write to last.hypo_inv */
			sprintf(fname, "%slast.hypo_inv", fpathname);
//printf(">>%s\n", fname);
			if ((fp_tmp = fopen(fname, "w")) != NULL) {
				WriteHypoInverseArchive(fp_tmp, &Hypocenter,
					Arrival, NumArrivals, fnout, iWriteHypHeader[ngrid], 1);
				fclose(fp_tmp);
			}
		}

		if (iSaveAlberto4Sum) {
			/* write Alberto 4 SIMULPS format */
			WriteHypoAlberto4(pSumFileAlberto4[ngrid], &Hypocenter,
				Arrival, fnout);
		}

		iWriteHypHeader[ngrid] = 0;

	return(0);
}



/*** function read observation file and open station grid files */

#define ARRIVAL_NULL_STR "?"

int GetObservations(FILE* fp_obs, char* ftype_obs, char* fn_grids,
	ArrivalDesc *arrival, int *pi_end_of_input,
	int* pnignore, int maxNumArrivals, HypoDesc* phypo,
	int* pMaxArrExceeded, int *pnumSArrivals)
{

	int nobs, nobs_read, istat, ntry, nLocate, n_compan;
	char filename[FILENAME_MAX];


	*pnumSArrivals = 0;


	/** read observations to arrival array */

	nobs = 0;
	*pnignore = 0;
	nLocate = 0;
	ntry = 0;
	while ((istat = GetNextObs(fp_obs, arrival + nobs,
				ftype_obs, ntry++ == 0)) != EOF) {

		if (istat == OBS_FILE_END_OF_EVENT)
			if (nobs == 0)
				continue;
			else
				break;
		if (istat == OBS_FILE_END_OF_INPUT) {
				*pi_end_of_input = 1;
				break;
		}
		if (istat == OBS_FILE_INVALID_PHASE)
			continue;
		if (istat == OBS_FILE_INVALID_DATE)
			continue;
		if (istat == OBS_FILE_SKIP_INPUT_LINE)
			continue;
		if (nobs == MAX_NUM_ARRIVALS - 1) {
			if (!*pMaxArrExceeded) {
				*pMaxArrExceeded = 1;
				putmsg(1,
"WARNING: maximum number of arrivals exceeded.");
			}
			continue;
		}


		/* check for aliased arrival label */
		if ((istat = EvaluateArrivalAlias(arrival + nobs)) < 0)
			;

		/* check for time delays */
		if ((istat = ApplyTimeDelays(arrival + nobs)) < 0)
			;

		/* set some fields */
		InitializeArrivalFields(arrival + nobs);

		/* check some fields */
		if (!isgraph(arrival[nobs].comp[0]))
			strcpy(arrival[nobs].comp, ARRIVAL_NULL_STR);
		if (!isgraph(arrival[nobs].onset[0]))
			strcpy(arrival[nobs].onset, ARRIVAL_NULL_STR);
		if (!isgraph(arrival[nobs].first_mot[0]))
			strcpy(arrival[nobs].first_mot, ARRIVAL_NULL_STR);
 		if (arrival[nobs].coda_dur < VERY_SMALL_DOUBLE)
 			arrival[nobs].coda_dur = CODA_DUR_NULL;
 		if (arrival[nobs].amplitude < VERY_SMALL_DOUBLE)
 			arrival[nobs].amplitude = AMPLITUDE_NULL;
 		if (arrival[nobs].period < VERY_SMALL_DOUBLE)
 			arrival[nobs].period = PERIOD_NULL;


		/* display arrival parameters */
		sprintf(MsgStr,
"Arrival %d:  %s (%s)  %s %s %s %d  %4.4d %2.2d %2.2d   %2.2d %2.2d %lf  Unc: %s %lf  Amp: %lf  Dur: %lf  Per: %lf",
			nobs,
			arrival[nobs].label,
			arrival[nobs].time_grid_label,
			arrival[nobs].onset,
			arrival[nobs].phase,
			arrival[nobs].first_mot,
			arrival[nobs].quality,
			arrival[nobs].year,
			arrival[nobs].month,
			arrival[nobs].day,
			arrival[nobs].hour,
			arrival[nobs].min,
			arrival[nobs].sec,
			arrival[nobs].error_type,
			arrival[nobs].error,
			arrival[nobs].amplitude,
			arrival[nobs].coda_dur,
			arrival[nobs].period);
		putmsg(2, MsgStr);


		/* remove blanks/whitespace from phase string */
		removeSpace(arrival[nobs].phase);
		/* check for valid P or S phase code */
		if (!IsPhaseID(arrival[nobs].phase, "P")
				&& !IsPhaseID(arrival[nobs].phase, "S")) {
			sprintf(MsgStr,
"WARNING: phase code not in P or S phase id list, rejecting observation: %s %s", arrival[nobs].label, arrival[nobs].phase);
			putmsg(2, MsgStr);
			continue;
		}

		/* check for duplicate arrival */
		if (IsDuplicateArrival(arrival, nobs + 1, nobs, NULL) >= 0) {
			sprintf(MsgStr,
"WARNING: duplicate arrival, rejecting observation: %s %s", arrival[nobs].label, arrival[nobs].phase);
			putmsg(2, MsgStr);
			continue;
		}

		nobs++;
	}

	nobs_read = nobs;



	/** check for minimum number of arrivals */

//printf("*pi_end_of_input %d  nobs_read %d  MinNumArrLoc %d\n", *pi_end_of_input, nobs_read, MinNumArrLoc);

	if (!(*pi_end_of_input) && (nobs_read > 0 && nobs_read < MinNumArrLoc)) {
		sprintf(MsgStr,
"WARNING: too few observations to locate, skipping event.");
		putmsg(1, MsgStr);
		sprintf(MsgStr,
"INFO: %d observations needed (specified in control file entry LOCMETH).",
		MinNumArrLoc);
		putmsg(2, MsgStr);
		return(-1);
	}



	/** process arrivals, determine if reject or ignore */


	/* homogenize date/time */
	if ((istat = HomogDateTime(arrival, nobs_read, phypo)) < 0) {
		puterr("ERROR: in arrival date/times.");
		if (istat == OBS_FILE_ARRIVALS_CROSS_YEAR_BOUNDARY)
			puterr("ERROR: arrivals cross year boundary.");
		return(-1);
	}


	/* sort to get arrivals in time order */
	if ((istat = SortArrivalsTime(arrival, nobs_read)) < 0) {
		puterr("ERROR: sorting arrivals by time.");
		return(-1);
	}


	/* check each arrival and initialize */

	Num3DGridReadToMemory = 0;
	for (nobs = 0; nobs < nobs_read; nobs++) {

		sprintf(MsgStr, "Checking Arrival %d:  %s (%s)  %s %s %s %d",
			nobs,
			arrival[nobs].label,
			arrival[nobs].time_grid_label,
			arrival[nobs].onset,
			arrival[nobs].phase,
			arrival[nobs].first_mot,
			arrival[nobs].quality);
		putmsg(2, MsgStr);

		/* if Vp/Vs > 0, check for companion phase,
				its time grid will be used for times */
		strcpy(arrival[nobs].fileroot, "\0");
		arrival[nobs].n_companion = -1;
		arrival[nobs].tfact = 1.0;
		if (VpVsRatio > 0.0) {
			/* try finding previously initialized companion phase */
			if (IsPhaseID(arrival[nobs].phase, "S") &&
					(n_compan =
					IsDuplicateArrival(arrival, nobs, nobs, "P")) >= 0 &&
					arrival[n_compan].flag_ignore == 0) {
				arrival[nobs].tfact = VpVsRatio;
				arrival[nobs].gdesc.type = arrival[n_compan].gdesc.type;
				arrival[nobs].station = arrival[n_compan].station;
				arrival[nobs].n_companion = n_compan;
				sprintf(MsgStr,
"INFO: using companion phase %d travel time grids.", arrival[nobs].n_companion);
					putmsg(2, MsgStr);
				// save filename as grid identifier (needed for GridMemList)
				strcpy(arrival[nobs].gdesc.title, arrival[n_compan].gdesc.title);
			}
		}

		/* no previously initialized companion phase */
		if (arrival[nobs].n_companion < 0) {

			/* try to open corresponding time grid file */
			if (IsPhaseID(arrival[nobs].phase, "P"))
				sprintf(arrival[nobs].fileroot, "%s.%s.%s", fn_grids,
					"P", arrival[nobs].time_grid_label);
			else if (IsPhaseID(arrival[nobs].phase, "S"))
				sprintf(arrival[nobs].fileroot, "%s.%s.%s", fn_grids,
					"S", arrival[nobs].time_grid_label);
			sprintf(filename, "%s.time", arrival[nobs].fileroot);
			/* try opening time grid file for this phase */
			istat = OpenGrid3dFile(filename,
					&(arrival[nobs].fpgrid),
					&(arrival[nobs].fphdr),
					&(arrival[nobs].gdesc), "time",
					&(arrival[nobs].station));

			/* try opening time grid file for uninitialized P companion phase */
			if (istat < 0 && VpVsRatio > 0.0 &&
					IsPhaseID(arrival[nobs].phase, "S")) {
				arrival[nobs].tfact = VpVsRatio;
				sprintf(arrival[nobs].fileroot, "%s.%s.%s", fn_grids,
						"P", arrival[nobs].time_grid_label);
				sprintf(filename, "%s.time", arrival[nobs].fileroot);
				istat = OpenGrid3dFile(filename,
					&(arrival[nobs].fpgrid),
					&(arrival[nobs].fphdr),
					&(arrival[nobs].gdesc), "time",
					&(arrival[nobs].station));
				sprintf(MsgStr,
"INFO: using P phase travel time grid file: %s", filename);
					putmsg(2, MsgStr);
			}
			// save filename as grid identifier (needed for GridMemList)
			strcpy(arrival[nobs].gdesc.title, filename);

			if (istat < 0) {
				sprintf(MsgStr,
"WARNING: cannot open time grid file: %s: rejecting observation: %s %s",
filename, arrival[nobs].label, arrival[nobs].phase);
				putmsg(2, MsgStr);
				CloseGrid3dFile(&(arrival[nobs].fpgrid),
					&(arrival[nobs].fphdr));
				strcpy(arrival[nobs].fileroot, "\0");
				goto RejectArrival;
			}


			/* check that search grid is inside time grid (3D grids) */

			if (arrival[nobs].gdesc.type == GRID_TIME &&
				!IsGridInside(LocGrid + 0, &(arrival[nobs].gdesc), 0)) {
				sprintf(MsgStr,
"WARNING: initial location search grid not contained inside arrival time grid, rejecting observation: %s %s", arrival[nobs].label, arrival[nobs].phase);
				putmsg(1, MsgStr);
				CloseGrid3dFile(&(arrival[nobs].fpgrid),
					&(arrival[nobs].fphdr));
				goto RejectArrival;
			}


			/* (2D grids) check that
				(1) greatest distance from search grid to station
					is inside time grid, and
				(2) distance from center of search grid to station
					is less than DistStaGridMax
							(if DistStaGridMax > 0) */

			if (arrival[nobs].gdesc.type == GRID_TIME_2D &&
					(istat = IsGrid2DBigEnough(LocGrid + 0,
					&(arrival[nobs].gdesc),
					&(arrival[nobs].station),
					DistStaGridMax,
					LocGrid[0].origx + (LocGrid[0].dx
						* (double) (LocGrid[0].numx - 1)) / 2.0,
					LocGrid[0].origy + (LocGrid[0].dy
						* (double) (LocGrid[0].numy - 1)) / 2.0)
				) != 1)
			{
				CloseGrid3dFile(&(arrival[nobs].fpgrid),
					&(arrival[nobs].fphdr));
				if (istat == -1) {
					sprintf(MsgStr,
"WARNING: greatest distance from initial 3D location search grid to station \n\texceeds 2D time grid size, rejecting observation: %s %s", arrival[nobs].label, arrival[nobs].phase);
					putmsg(2, MsgStr);
					goto RejectArrival;
				} else if (istat == -2) {
					sprintf(MsgStr,
"WARNING: distance from grid center to station \n\texceeds maximum station distance, ignoring observation in misfit calculation: %s %s", arrival[nobs].label, arrival[nobs].phase);
					putmsg(2, MsgStr);
					arrival[nobs].flag_ignore = 1;
					goto IgnoreArrival;
				}
				if (istat == -3) {
					sprintf(MsgStr,
"WARNING: depth range of initial 3D location search grid exceeds that of 2D time grid size, rejecting observation: %s %s", arrival[nobs].label, arrival[nobs].phase);
					putmsg(2, MsgStr);
					goto RejectArrival;
				}
			}

		} /* if (arrival[nobs].n_companion < 0) */



		/* check if arrival is excluded */

		if (isExcluded(arrival[nobs].label, arrival[nobs].phase)) {
			sprintf(MsgStr,
"INFO: arrival is excluded, ignoring observation in misfit calculation: %s %s",
				arrival[nobs].label, arrival[nobs].phase);
			arrival[nobs].flag_ignore = 1;
			CloseGrid3dFile(&(arrival[nobs].fpgrid),
				&(arrival[nobs].fphdr));
			putmsg(1, MsgStr);
			goto IgnoreArrival;
		}



		/* check if maximum number of arrivals for location exceeded */

		if (nLocate == maxNumArrivals) {
			putmsg(2,
"WARNING: maximum number of arrivals for location exceeded, \n\tignoring observation in misfit calculation.");
			arrival[nobs].flag_ignore = 1;
			CloseGrid3dFile(&(arrival[nobs].fpgrid),
				&(arrival[nobs].fphdr));
			goto IgnoreArrival;
		}



		/** arrival to be used for location */

		/* no further processing for arrivals with companion time grids */
		if (arrival[nobs].n_companion >= 0)
			goto AcceptArrival;

		/** prepare time grids access in memory or on disk */

		/* construct dual-sheet description
			(2D grids for all search types,
				and 3D grids for grid-search) */

		if (SearchType == SEARCH_GRID
				|| arrival[nobs].gdesc.type == GRID_TIME_2D) {

			arrival[nobs].sheetdesc = arrival[nobs].gdesc;
			arrival[nobs].sheetdesc.numx = 2;
			/* allocate grid */
			arrival[nobs].sheetdesc.buffer =
				AllocateGrid(&(arrival[nobs].sheetdesc));
			if (arrival[nobs].sheetdesc.buffer == NULL) {
				puterr(
"ERROR: allocating memory for arrival sheet buffer.");
				goto RejectArrival;
				//return(EXIT_ERROR_MEMORY);
			}
			/* create array access pointers */
			arrival[nobs].sheetdesc.array =
				CreateGridArray(&(arrival[nobs].sheetdesc));
			if (arrival[nobs].sheetdesc.array == NULL) {
				puterr(
"ERROR: creating array for accessing arrival sheet buffer.");
				goto RejectArrival;
				//return(EXIT_ERROR_MEMORY);
			}
			arrival[nobs].sheetdesc.origx = VERY_LARGE_DOUBLE;

		}


		/* read 3D grid into memory (3D grids for Metropolis or Octtree search) */

		if ((SearchType == SEARCH_MET || SearchType == SEARCH_OCTTREE)
				&& arrival[nobs].gdesc.type == GRID_TIME
				&& (MaxNum3DGridMemory < 0 ||
				Num3DGridReadToMemory < MaxNum3DGridMemory)) {

			/* allocate grid */
			arrival[nobs].gdesc.buffer =
				NLL_AllocateGrid(&(arrival[nobs].gdesc));
			if (arrival[nobs].gdesc.buffer == NULL) {
				puterr(
"ERROR: allocating memory for arrival time grid buffer; grid will be read from disk.");
			} else {
				/* create array access pointers */
				arrival[nobs].gdesc.array =
					NLL_CreateGridArray(&(arrival[nobs].gdesc));
				if (arrival[nobs].gdesc.array == NULL) {
					puterr(
"ERROR: creating array for accessing arrival time grid buffer.");
					goto RejectArrival;
					//return(EXIT_ERROR_MEMORY);
				}
				/* read time grid */
				if (NLL_ReadGrid3dBuf(&(arrival[nobs].gdesc),
						arrival[nobs].fpgrid) < 0) {
					puterr(
"ERROR: reading arrival time grid buffer.");
					goto RejectArrival;
					//return(EXIT_ERROR_MEMORY);
				}
				CloseGrid3dFile(&(arrival[nobs].fpgrid),
					&(arrival[nobs].fphdr));
				Num3DGridReadToMemory++;
			}
		}


		/* read time grid and close file (2D grids)*/

		if (arrival[nobs].gdesc.type == GRID_TIME_2D) {
			istat = ReadArrivalSheets(1, &(arrival[nobs]), 0.0);
			CloseGrid3dFile(&(arrival[nobs].fpgrid),
				&(arrival[nobs].fphdr));
			if (istat  < 0) {
				sprintf(MsgStr,
"ERROR: reading arrival travel time sheets (2D grid), rejecting observation: %s %s", arrival[nobs].label, arrival[nobs].phase);
				puterr(MsgStr);
				goto RejectArrival;
			}
		}

		/* arrival accepted, use for location */
		AcceptArrival:
		arrival[nobs].flag_ignore = 0;
		nLocate++;
		if (IsPhaseID(arrival[nobs].phase, "S"))
			(*pnumSArrivals)++;


		/* arrival accepted, ignore for location */
		IgnoreArrival:

		*pnignore += arrival[nobs].flag_ignore;

		continue;


		/* arrival rejected */
		RejectArrival:
		arrival[nobs].flag_ignore = 999;
		sprintf(MsgStr, "   Rejected Arrival %d:  %s (%s)  %s %s %s %d",
			nobs,
			arrival[nobs].label,
			arrival[nobs].time_grid_label,
			arrival[nobs].onset,
			arrival[nobs].phase,
			arrival[nobs].first_mot,
			arrival[nobs].quality);
		putmsg(2, MsgStr);

	}


	/* sort to get rejected arrivals at end of arrivals array */

	if ((istat = SortArrivalsIgnore(arrival, nobs_read)) < 0) {
		puterr("ERROR: sorting arrivals by ignore flag.");
			return(-1);
	}


	/* avoid returning 0 if arrivals were read, return 0 indicates end of file */
	if (nLocate + *pnignore == 0 && nobs_read > 0) {
		sprintf(MsgStr,
			"WARNING: %d arrivals read, but none accepted for location.", nobs_read);
		putmsg(2, MsgStr);
		return(-1);
	}

	return(nLocate + *pnignore);
}



/*** function to initialize arrival fields */

void InitializeArrivalFields(ArrivalDesc *arrival)
{

 	arrival->obs_centered = 0.0;
 	arrival->pred_travel_time = 0.0;
 	arrival->pred_travel_time_best = -LARGE_DOUBLE;
	arrival->pred_centered = 0.0;
 	arrival->cent_resid = 0.0;
 	arrival->obs_travel_time = 0.0;
 	arrival->residual = 0.0;
	arrival->weight = 0.0;
 	arrival->dist = 0.0;
	arrival->azim = 0.0;
 	arrival->ray_azim = 0.0;
	arrival->ray_dip = 0.0;
	arrival->ray_qual = 0;
 	arrival->pdf_residual_sum = 0.0;
 	arrival->pdf_weight_sum = 0.0;

 	arrival->fpgrid = NULL;
 	arrival->fphdr = NULL;
	arrival->gdesc.buffer = NULL;
}


/*** function to test if arrival is excluded */

int isExcluded(char *label, char *phase) {

	int nexclude;

	for (nexclude = 0; nexclude < NumLocExclude; nexclude++) {
		if (strcmp(label, LocExclude[nexclude].label) == 0
				&& strcmp(phase, LocExclude[nexclude].phase) == 0)
			return(1);
	}

	return(0);

}



/*** function to evaluate arrival label aliases */

int EvaluateArrivalAlias(ArrivalDesc *arrival)
{
	int nAlias;
	int checkAgain = 1, icount = 0, aliasApplied = 0;
	char *pchr, tmpLabel[MAXLINE];


	strcpy(tmpLabel, arrival->label);

	sprintf(MsgStr, "Checking for station name alias: %s", tmpLabel);
	putmsg(3, MsgStr);


	/* evaluate aliases until no replacement done */

	while(checkAgain && icount < MAX_NUM_LOC_ALIAS_CHECKS) {

		checkAgain = 0;
		icount++;

		for (nAlias = 0; nAlias < NumLocAlias; nAlias++) {

			/* check if alias can be rejected */

			if (strcmp(LocAlias[nAlias].name, tmpLabel) != 0)
				continue;
			if (LocAlias[nAlias].byr > arrival->year)
				continue;
			if (LocAlias[nAlias].byr == arrival->year) {
				if (LocAlias[nAlias].bmo > arrival->month)
					continue;
				if (LocAlias[nAlias].bmo == arrival->month
					&& LocAlias[nAlias].bday > arrival->day)
				continue;
			}
			if (LocAlias[nAlias].eyr < arrival->year)
				continue;
			if (LocAlias[nAlias].eyr == arrival->year) {
				if (LocAlias[nAlias].emo < arrival->month)
					continue;
				if (LocAlias[nAlias].emo == arrival->month
					&& LocAlias[nAlias].eday < arrival->day)
				continue;
			}

			/* apply alias */

			aliasApplied = 1;

			strcpy(tmpLabel, LocAlias[nAlias].alias);
			sprintf(MsgStr, " -> %s", tmpLabel);
			putmsg(3, MsgStr);

			/* if alias label is same as arrival label,
				 end check to avoid infinite recurssion */
			if (strcmp(tmpLabel, arrival->label) == 0)
				checkAgain = 0;
			else
				checkAgain = 1;

			break;

		}
	}


	/* check for possible recursion in alias specifications */

	if (icount >= MAX_NUM_LOC_ALIAS_CHECKS) {
		putmsg(3, "");
		puterr(
"ERROR: possible infinite recursion in station name alias.");
		return(-1);
	}


	/* update arrival label */

	strcpy(arrival->time_grid_label, tmpLabel);
	if ((pchr = strrchr(tmpLabel, '_')) != NULL)
		*pchr = '\0';
//18JUL2002 AJL (IRSN)
// to prevent station label change in output file
//	strcpy(arrival->label, tmpLabel);


	/* return if no alias applied */

	if (!aliasApplied) {
		putmsg(3, "");
		return(0);
	}


	putmsg(3, "");


	return(0);
}



/*** function to apply time delays */

int ApplyTimeDelays(ArrivalDesc *arrival)
{
	int nDelay;
	double tmp_delay;


	sprintf(MsgStr, "Checking for time delay: %s %s",
		arrival->label, arrival->phase);
	putmsg(3, MsgStr);


	/* check time delays for match to label/phase */

	for (nDelay = 0; nDelay < NumTimeDelays; nDelay++) {

		if (strcmp(TimeDelay[nDelay].label, arrival->label) == 0
			&& strcmp(TimeDelay[nDelay].phase, arrival->phase) == 0)
		{
		    /* apply station delays */

		    tmp_delay = TimeDelay[nDelay].delay;
		    arrival->delay = 0.0;
		    if (IsPhaseID(arrival->phase, "P")
				&& fabs(tmp_delay) > VERY_SMALL_DOUBLE) {
// DELAY_CORR			arrival->sec -= tmp_delay;
			arrival->delay = tmp_delay;
			sprintf(MsgStr,
			    "   P delay of %lf sec subtracted from obs time.",
				tmp_delay);
			putmsg(3, MsgStr);
		    } else if (IsPhaseID(arrival->phase, "S")
				&& fabs(tmp_delay) > VERY_SMALL_DOUBLE) {
// DELAY_CORR			arrival->sec -= tmp_delay;
			arrival->delay = tmp_delay;
			sprintf(MsgStr,
			    "    S delay of %lf sec subtracted from obs time.",
				tmp_delay);
			putmsg(3, MsgStr);
		    }
		    break;
		}
	}


	putmsg(3, "");


	return(0);
}




/*** function to extract arrival information observation file name */

/* returns: 0 if nothing done, 1 if info extracted, -1 if unexpected error */

int ExtractFilenameInfo(char *filename, char *type_obs)
{
	int istat;
	char *filepos, *extpos;

	if (strcmp(ftype_obs, "RENASS_DEP") == 0) {
		/* find beginning of filename */
		if ((filepos = strrchr(filename, '/')) == NULL)
				return(-1);
		/* get date/time from filename */
		/* try long format (i.e. NICE199801311202.dep) */
		if ((extpos = strstr(filepos, ".dep")) != NULL &&
				extpos - filepos - 12 >= 0)
		{
			if ((istat = sscanf(extpos - 12, "%4d%2d%2d%2d%2d",
				&EventTime.year, &EventTime.month, &EventTime.day,
				&EventTime.hour, &EventTime.min))
					!= 5)
				return(-1);
		}
		/* try short format (i.e. g504210802.dep) */
		else if ((extpos = strstr(filepos, ".dep")) != NULL &&
				extpos - filepos - 9 >= 0)
		{
			if ((istat = sscanf(extpos - 9, "%1d%2d%2d%2d%2d",
				&EventTime.year, &EventTime.month, &EventTime.day,
				&EventTime.hour, &EventTime.min))
					!= 5)
				return(-1);
			EventTime.year += 1990;
		}

		return(1);
	}

	return(0);

}


/*** function to read arrival from observation file */

int GetNextObs(FILE* fp_obs, ArrivalDesc *arrival, char* ftype_obs, int nfirst)
{
	int istat, n, iloop;
	char *cstat;
	char chr, instruction[MAXSTRING];
	char chrtmp[MAXSTRING];

	double psec;
	double weight, ttime;

	int ioff;

	static char line[MAXLINE_LONG];

	static int date_saved, year_save, month_save, day_save;
	static int hour_save, min_save;
	static int check_for_S_arrival;

	static int in_hypocenter_event;

	int ifound;
	
	// ETH LOC format
	char eth_line_key[MAXSTRING];
	int eth_use_loc, eth_use_mag;


	/* if no obs read for this event, set date saved flag to 0 */
	if (nfirst) {
		date_saved = 0;
		check_for_S_arrival = 0;
		in_hypocenter_event = 0;
	}


	/* check for special control instructions */
	iloop = 1;
	while (iloop && !check_for_S_arrival) {
	    iloop = 0;
	    if ((chr = fgetc(fp_obs)) == '!'/* || chr == '#'*/) {
		ungetc(chr, fp_obs);
		fgets(line, MAXLINE_LONG, fp_obs);
		/* read instruction */
		istat = sscanf(line + 1, "%s", instruction);
		if (istat == EOF)
			return(OBS_FILE_END_OF_INPUT);
		if (istat == 1) {
			if (strcmp(instruction, "END_EVENT") == 0) {
				return(OBS_FILE_END_OF_EVENT);
			} else if (strcmp(instruction, "END_FILE") == 0) {
				return(OBS_FILE_END_OF_INPUT);
			} else if (strcmp(instruction, "SKIP_NEXT_LINE") == 0) {
				fgets(line, MAXLINE_LONG, fp_obs);
				iloop = 1;
				continue;
			} else {		// skip this line
				iloop = 1;
				continue;
			}
		}
		if (istat != 1) {
			puterr2("WARNING: unrecognized control statement",
				line);
		}
	    } else {
		ungetc(chr, fp_obs);
	    }
	}


	/* set field defaults */
 	strcpy(arrival->label, "?");
 	strcpy(arrival->inst, "?");
 	strcpy(arrival->comp, "?");
 	strcpy(arrival->onset, "?");
 	strcpy(arrival->phase, "?");
 	strcpy(arrival->first_mot, "?");
 	arrival->quality = 99;
 	arrival->coda_dur = CODA_DUR_NULL;
 	arrival->amplitude = AMPLITUDE_NULL;
 	arrival->period = PERIOD_NULL;
 	arrival->amp_mag = MAG_NULL;
 	arrival->dur_mag = MAG_NULL;


	/* attempt to read obs based on obs file type */

	if (strcmp(ftype_obs, "NLLOC_OBS") == 0)
	{

		/* read next line */
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
 		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);

		istat = ReadArrival(line, arrival, IO_ARRIVAL_OBS);

		if (istat < 1)
			return(OBS_FILE_END_OF_EVENT);

		/* convert error to quality */
 		if ((arrival->quality = Err2Qual(arrival)) < 0)
			 arrival->quality = 99;
		return(istat);
	}
	
	

	else if (strcmp(ftype_obs, "HYPO71") == 0 ||
		strcmp(ftype_obs, "HYPO71_OV") == 0 ||
		strcmp(ftype_obs, "HYPOELLIPSE") == 0)
	{

		if (check_for_S_arrival) {
			/* check for S phase input in last input line read */

			/* set S read offset to allow correction of incorrect hypo71 format */
			ioff = 0;
			if (strcmp(ftype_obs, "HYPO71_OV") == 0) {
				istat = ReadFortranString(line, 31, 1, chrtmp);
				if  (strcmp(chrtmp, " ") != 0)
					ioff = -1;
			}

			/* check for zero or blank S phase time */
			istat = ReadFortranString(line, 32 + ioff, 5, chrtmp);
			if (istat > 0
				&& strncmp(chrtmp, "     ", 5) != 0
				&& strncmp(chrtmp, "    0", 5) != 0)
			{
				/* read S phase input in last input line read */
				istat = ReadFortranString(line, 1, 4, arrival->label);
				TrimString(arrival->label);
				istat += ReadFortranInt(line, 10, 2, &arrival->year);
		// temporary fix for HYPO71 Y2K problem
		if (arrival->year < 20)
			arrival->year += 100;
				arrival->year += 1900;
				istat += ReadFortranInt(line, 12, 2, &arrival->month);
				istat += ReadFortranInt(line, 14, 2, &arrival->day);
				istat += ReadFortranInt(line, 16, 2, &arrival->hour);
				istat += ReadFortranInt(line, 18, 2, &arrival->min);
				istat += ReadFortranReal(line, 32 + ioff, 5, &arrival->sec);
				/* check for integer sec format */
				if (arrival->sec > 99.999)
				arrival->sec /= 100.0;
				istat += ReadFortranReal(line, 20, 5, &psec);
				/* check for P second >= 60.0 */
				if (psec >= 60.0 && arrival->sec < 60.0)
					arrival->sec += 60.0;
				arrival->phase[0] = '\0';
				istat += ReadFortranString(line, 38 + ioff, 1, arrival->phase);
				TrimString(arrival->phase);
				istat += ReadFortranInt(line, 40 + ioff, 1, &arrival->quality);
			}
		}

		/* check for S arrival input found */
		if (check_for_S_arrival && istat == 10
				&& IsPhaseID(arrival->phase, "S")
				&& IsGoodDate(arrival->year,
						arrival->month, arrival->day)) {

			//strcpy(arrival->phase, "S");

			/* set error fields */
			strcpy(arrival->error_type, "GAU");
			if (arrival->quality >= 0 &&
			    arrival->quality < NumQuality2ErrorLevels) {
				arrival->error =
					Quality2Error[arrival->quality];
			} else {
				arrival->error =
				   Quality2Error[NumQuality2ErrorLevels - 1];
				puterr("WARNING: invalid arrival weight.");
			}

			line[0] = '\0';
			check_for_S_arrival = 0;
			return(istat);
		} else if (check_for_S_arrival) {
			check_for_S_arrival = 0;
			return(OBS_FILE_SKIP_INPUT_LINE);
		}




		/* read next line */
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);

		/* read formatted P (or S) arrival input */
		istat = ReadFortranString(line, 1, 4, arrival->label);
		TrimString(arrival->label);
		istat += ReadFortranString(line, 5, 1, arrival->onset);
		TrimString(arrival->onset);
		istat += ReadFortranString(line, 6, 1, arrival->phase);
		TrimString(arrival->phase);
		istat += ReadFortranString(line, 7, 1, arrival->first_mot);
		TrimString(arrival->first_mot);
		istat += ReadFortranInt(line, 8, 1, &arrival->quality);
		istat += ReadFortranInt(line, 10, 2, &arrival->year);
// temporary fix for HYPO71 Y2K problem
if (arrival->year < 20)
	arrival->year += 100;
		arrival->year += 1900;
		istat += ReadFortranInt(line, 12, 2, &arrival->month);
		istat += ReadFortranInt(line, 14, 2, &arrival->day);
		istat += ReadFortranInt(line, 16, 2, &arrival->hour);
		istat += ReadFortranInt(line, 18, 2, &arrival->min);
		istat += ReadFortranReal(line, 20, 5, &arrival->sec);

		if (istat != 11) {
			line[0] = '\0';
			return(OBS_FILE_END_OF_EVENT);
		}

		/* read optional amplitude/period fields */
		istat += ReadFortranReal(line, 44, 4, &arrival->amplitude);
//printf("AMPLITUDE: %lf (%s)\n", arrival->amplitude, line);
		istat += ReadFortranReal(line, 48, 3, &arrival->period);
// 20030826 AJL coda_dur added
		istat += ReadFortranReal(line, 71, 5, &arrival->coda_dur);


		/* check for valid phase code */
		if (IsPhaseID(arrival->phase, "P")) {
			//strcpy(arrival->phase, "P");
			check_for_S_arrival = 1;

		} else if (IsPhaseID(arrival->phase, "S")) {
			//strcpy(arrival->phase, "S");
			check_for_S_arrival = 0;
		} else
			return(OBS_FILE_END_OF_EVENT);

		if (!IsGoodDate(arrival->year, arrival->month, arrival->day))
			return(OBS_FILE_END_OF_EVENT);


		/* check for integer sec format */
		if (arrival->sec > 99.999)
			arrival->sec /= 100.0;

		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}

	
	
	else if (strcmp(ftype_obs, "NCSN_Y2K_5") == 0)
	{

		if (check_for_S_arrival) {
			/* check for S phase input in last input line read */

			/* check for zero or blank S phase time */
			istat = ReadFortranString(line, 48, 1, chrtmp);
			if (istat > 0
				&& (strncmp(chrtmp, "S", 1) == 0
				|| strncmp(chrtmp, "s", 1) == 0))
			{
				/* read S phase input in last input line read */
				istat = ReadFortranString(line, 1, 5, arrival->label);
				TrimString(arrival->label);
				istat += ReadFortranInt(line, 18, 4, &arrival->year);
				istat += ReadFortranInt(line, 22, 2, &arrival->month);
				istat += ReadFortranInt(line, 24, 2, &arrival->day);
				istat += ReadFortranInt(line, 26, 2, &arrival->hour);
				istat += ReadFortranInt(line, 28, 2, &arrival->min);
				istat += ReadFortranReal(line, 42, 5, &arrival->sec);
				/* check for integer sec format */
//				if (arrival->sec > 99.999)
					arrival->sec /= 100.0;
				istat += ReadFortranReal(line, 30, 5, &psec);
					psec /= 100.0;
				/* check for P second >= 60.0 */
				if (psec >= 60.0 && arrival->sec < 60.0)
					arrival->sec += 60.0;
				arrival->phase[0] = '\0';
				istat += ReadFortranString(line, 48, 1, arrival->phase);
				TrimString(arrival->phase);
				istat += ReadFortranString(line, 47, 1, arrival->onset);
				TrimString(arrival->onset);
				istat += ReadFortranInt(line, 50, 1, &arrival->quality);
			}
		}

		/* check for S arrival input found */
		if (check_for_S_arrival && istat == 11
				&& IsPhaseID(arrival->phase, "S")
				&& IsGoodDate(arrival->year,
						arrival->month, arrival->day)) {

			/* set error fields */
			strcpy(arrival->error_type, "GAU");
			if (arrival->quality >= 0 &&
			    arrival->quality < NumQuality2ErrorLevels) {
				arrival->error =
					Quality2Error[arrival->quality];
			} else {
				arrival->error =
				   Quality2Error[NumQuality2ErrorLevels - 1];
				puterr("WARNING: invalid arrival weight.");
			}

			line[0] = '\0';
			check_for_S_arrival = 0;
			return(istat);
		} else if (check_for_S_arrival) {
			check_for_S_arrival = 0;
			return(OBS_FILE_SKIP_INPUT_LINE);
		}




		/* read next line */
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);

		/* read formatted P (or S) arrival input */
		istat = ReadFortranString(line, 1, 5, arrival->label);
		TrimString(arrival->label);
		istat += ReadFortranString(line, 9, 1, arrival->comp);
		TrimString(arrival->comp);
		istat += ReadFortranString(line, 10, 3, arrival->inst);
		TrimString(arrival->inst);
		istat += ReadFortranString(line, 14, 1, arrival->onset);
		TrimString(arrival->onset);
		if (strpbrk(arrival->onset, "XYZ") != NULL) {	// skip RTP (phases with X,Y,Z onset)
			return(OBS_FILE_SKIP_INPUT_LINE);
		}
		istat += ReadFortranString(line, 15, 1, arrival->phase);
		TrimString(arrival->phase);
		istat += ReadFortranString(line, 16, 1, arrival->first_mot);
		TrimString(arrival->first_mot);
		istat += ReadFortranInt(line, 17, 1, &arrival->quality);
/*
		if (arrival->quality > 3) {			// skip phases with quality > 4
			return(OBS_FILE_SKIP_INPUT_LINE);
		}
*/
		istat += ReadFortranInt(line, 18, 4, &arrival->year);
		istat += ReadFortranInt(line, 22, 2, &arrival->month);
		istat += ReadFortranInt(line, 24, 2, &arrival->day);
		istat += ReadFortranInt(line, 26, 2, &arrival->hour);
		istat += ReadFortranInt(line, 28, 2, &arrival->min);
		istat += ReadFortranReal(line, 30, 5, &arrival->sec);
		/* check for integer sec format */
//		if (arrival->sec > 99.999)
			arrival->sec /= 100.0;

		if (istat != 13) {
			line[0] = '\0';
			return(OBS_FILE_END_OF_EVENT);
		}

		/* read optional amplitude/period fields */
		istat += ReadFortranReal(line, 55, 7, &arrival->amplitude);
			arrival->amplitude /= 100.0;
		istat += ReadFortranReal(line, 84, 3, &arrival->period);
			arrival->period /= 100.0;


		/* check for valid phase code */
		if (IsPhaseID(arrival->phase, "P")) {
			//strcpy(arrival->phase, "P");
			check_for_S_arrival = 1;

		} else if (IsPhaseID(arrival->phase, "S")) {
			//strcpy(arrival->phase, "S");
			check_for_S_arrival = 0;
		} else
			return(OBS_FILE_END_OF_EVENT);

		if (!IsGoodDate(arrival->year, arrival->month, arrival->day))
			return(OBS_FILE_END_OF_EVENT);

		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}
	
	

	else if (strcmp(ftype_obs, "PAOLO_OV") == 0 ||
		strcmp(ftype_obs, "ALBERTO_3D") == 0 ||
		strcmp(ftype_obs, "ALBERTO_3D_4") == 0 ||
		strcmp(ftype_obs, "SIMULPS") == 0)
	{
/*
"PAOLO_OV"
87 119 23 5 45.00 40 49.17  14 25.80   0.43   0.00
 BKEP0 45.51 BKWP0 45.61 BKNP0 45.66 BAFP0 45.74 BKES1 45.86 BKWS1 46.04
 BKNS1 46.06 BAFS1 46.20

"ALBERTO_3D"
891018  0 4 15.34 37  2.66 121 52.54  16.92   5.80
JTGVIPU018.3000JBZVIPU018.7600JECVIPD018.5700JRGVIPU018.4500JHLVIPD018.5500
HLTVIPU024.9600BSLVIPU025.5800BPCVEPD025.4900HJSVIPU025.4900

"SIMULPS"
 from Stephan Husen   email:  stephan@tomo.ig.erdw.ethz.ch
summary line containing event origin time, location, and magnitude
 format: a4,a2,1x,a2,i2,1x,f5.2,i3,a1,f5.2,1x,i3,a1,f5.2,f7.2,f7.2)
           iyrmo, iday, ihr, mino, seco, ltde, cns, eltm, lnde, cew, elnm, dep, mag
observed traveltimes to stations, six sets per line
 each set is in the format: a4,a4,f6.2       sta, rmki, tt        rmki -> (I,E)(P,S,SP)(U,D,-)(0-4)
blank terminates event's data
example:
84 3 6  6 7 58.58 46N19.78   7E26.85   3.77   2.00
SIE_IP-1  1.45DIX2IP-1  4.94DIX_IP-1  4.95EMS_IP-1  8.14EMV_IP-1  8.54MMK2IP-1  8.58
MMK_IP-1  8.67STG_IP-1 16.48SLE_IP-1 27.54WIL_IP-1 28.14SAX_IP-1 30.78

84 3 7  4 7 56.92 46N20.89   7E26.59   6.09   1.50
ZZB_IP-1  1.03ZZE_IP-1  1.06ZZA_IP-1  1.07ZZC_IP-1  1.11ZZD_IP-1  1.11ZZF_IP-1  1.11


*/
		/* check for end of event (assumes no blanks after last phase) */
		chr = fgetc(fp_obs);
		if (chr == EOF) {
			return(OBS_FILE_END_OF_INPUT);
		} else if (chr == '\n') {
			if ((chr = fgetc(fp_obs)) != EOF && !isdigit(chr)) {
				// another phase follows
				ungetc(chr, fp_obs);
			} else if (chr == '0') {
				// end of event "0" line
				// read to end of line
				while ((chr = fgetc(fp_obs)) != EOF && chr != '\n')
					;
				return(OBS_FILE_END_OF_EVENT);
			} else {
				// end of event
				ungetc(chr, fp_obs);
				return(OBS_FILE_END_OF_EVENT);
			}
		} else {
			ungetc(chr, fp_obs);
		}


		/* check for event hypocenter line */
		if ((chr = fgetc(fp_obs)) != EOF && isdigit(chr)) {
			// read hypocenter time
			ungetc(chr, fp_obs);
			/* read next line */
			cstat = fgets(line, MAXLINE_LONG, fp_obs);
			if (cstat == NULL)
				return(OBS_FILE_END_OF_INPUT);
			istat = ReadFortranInt(line, 1, 2, &EventTime.year);
			istat += ReadFortranInt(line, 3, 2, &EventTime.month);
			istat += ReadFortranInt(line, 5, 2, &EventTime.day);
			istat += ReadFortranInt(line, 8, 2, &EventTime.hour);
			istat += ReadFortranInt(line, 10, 2, &EventTime.min);
			if (istat == 5) {
	// temporary fix for HYPO71 Y2K problem
	if (EventTime.year < 20)
		EventTime.year += 100;
				EventTime.year += 1900;
				// read to end of line
				//while ((chr = fgetc(fp_obs)) != EOF && chr != '\n')
				//	;
				if (chr == EOF)
					return(OBS_FILE_END_OF_INPUT);
			}
		} else {
			if (chr == EOF)
				return(OBS_FILE_END_OF_INPUT);
			ungetc(chr, fp_obs);
		}


		if (strcmp(ftype_obs, "PAOLO_OV") == 0) {
		// try to read phase arrival input
			istat = fscanf(fp_obs, " %s %lf", line, &arrival->sec);
			if (istat == EOF)
				return(OBS_FILE_END_OF_INPUT);
			istat += ReadFortranString(line, 1, 3, arrival->label);
			istat += ReadFortranString(line, 4, 1, arrival->phase);
			istat += ReadFortranInt(line, 5, 1, &arrival->quality);
			if (istat != 5)
				return(OBS_FILE_END_OF_EVENT);

		} else if (strcmp(ftype_obs, "SIMULPS") == 0) {
			// try to read phase arrival input
			istat = fscanf(fp_obs, "%14c", line);
			line[14]= '\0';
			if (istat == EOF)
				return(OBS_FILE_END_OF_INPUT);
// printf("<%s>\n", line);
			istat += ReadFortranString(line, 1, 4, arrival->label);
			istat += ReadFortranString(line, 5, 1, arrival->onset);
			istat += ReadFortranString(line, 6, 1, arrival->phase);
			istat += ReadFortranString(line, 7, 1, arrival->first_mot);
			istat += ReadFortranInt(line, 8, 1, &arrival->quality);
			istat += ReadFortranReal(line, 9, 6, &arrival->sec);

			if (istat != 7)
				return(OBS_FILE_END_OF_EVENT);

		} else if (strcmp(ftype_obs, "ALBERTO_3D") == 0
				|| strcmp(ftype_obs, "ALBERTO_3D_4") == 0) {
			// try to read phase arrival input
			istat = fscanf(fp_obs, "%15c", line);
			line[15]= '\0';
			if (istat == EOF)
				return(OBS_FILE_END_OF_INPUT);
//printf("<%s>\n", line);
			if (strcmp(ftype_obs, "ALBERTO_3D_4") == 0) {
				istat += ReadFortranString(line, 1, 4, arrival->label);
				istat++;
			} else {
				istat += ReadFortranString(line, 1, 3, arrival->label);
				istat += ReadFortranString(line, 4, 1, arrival->comp);
			}
			istat += ReadFortranString(line, 5, 1, arrival->onset);
			istat += ReadFortranString(line, 6, 1, arrival->phase);
			istat += ReadFortranString(line, 7, 1, arrival->first_mot);
			istat += ReadFortranInt(line, 8, 1, &arrival->quality);
			istat += ReadFortranReal(line, 9, 7, &arrival->sec);

			if (istat != 8)
				return(OBS_FILE_END_OF_EVENT);
       }


		arrival->year = EventTime.year;
		arrival->month = EventTime.month;
		arrival->day = EventTime.day;
		arrival->hour = EventTime.hour;
		arrival->min = EventTime.min;


		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}



	else if (strcmp(ftype_obs, "ETH_LOC") == 0)
	{

/* example:
 1 1   200.0   300.0    1.71 1                        INST
   0.000     0.000  10.00 0000 00 00 00 00  0.00 0    TRIAL
maraini                                               AUTOR
2001/12/04 01:29                                      Local           20010748
FUORN   P       E   14.910 1 0.14      448 NWA       0
FUORN   S       E   16.140 1 0.14      448 NWA       0
OSS2    P       ID  16.428 1 0.09     2554 AHP       0
OSS     P       ED  16.524 1 0.11     1672 AHP       0
BERNI   P       E   17.930 1 0.18       46 NWA       0
BERNI   S       E   21.139 1 0.18       46 NWA       0
CHDAW   P       Q   20.016 1 0.56       10 NWA       0
VDL2    P       IU  23.704 1 0.20      246 AHP       0
VDL     P       E   23.863 1 0.12      188 AHP       0
VDL     P       E   23.822 1 0.28        4 NWA       1
VDL     S       E   31.177 1 0.28        4 NWA       1
                                                      SKIP
20011204012913146568N010226E00614Ml1239814161001001000SEDG024196006511         B
KP200112040128                                        ARCHIVE
*/

		// read line
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);

		// check for end of event (assumes empty lines between event)
		if (LineIsBlank(line)) {
			// end of event
			return(OBS_FILE_END_OF_EVENT);
		}


		// read event hypocenter line
		if (!in_hypocenter_event) {
			// find origin time line
			ifound = 0;
			while((istat = ReadFortranString(line, 55, 4, eth_line_key)) > 0) {
				if (strcmp(eth_line_key, "Loca") == 0 
						|| strcmp(eth_line_key, "Regi") == 0
						|| strcmp(eth_line_key, "Tele") == 0) {
					ifound = 1;
					break;
				}
				// read next line
				cstat = fgets(line, MAXLINE_LONG, fp_obs);
				if (cstat == NULL)
					return(OBS_FILE_END_OF_INPUT);
			}
			if (!ifound)
				return(OBS_FILE_END_OF_EVENT);
				
			// read hypocenter time
			//2001/12/04 01:29                                      Local           20010748
			istat = ReadFortranInt(line, 1, 4, &EventTime.year);
			istat += ReadFortranInt(line, 6, 2, &EventTime.month);
			istat += ReadFortranInt(line, 9, 2, &EventTime.day);
			istat += ReadFortranInt(line, 12, 2, &EventTime.hour);
			istat += ReadFortranInt(line, 15, 2, &EventTime.min);
//printf("EventTime: %d/%d/%d %d:%d\n",
//EventTime.year, EventTime.month, EventTime.day,
//EventTime.hour, EventTime.min);
			if (istat == 5) {
				in_hypocenter_event = 1;
				// read until phase line reached
				while((istat = ReadFortranString(line, 55, 4, eth_line_key)) > 0
						&& strcmp(line, "    ") != 0) {
					cstat = fgets(line, MAXLINE_LONG, fp_obs);
					if (cstat == NULL)
						return(OBS_FILE_END_OF_INPUT);
				}
			} else {
				return(OBS_FILE_END_OF_EVENT);
			}
		}

		// check if valid phase line: key is blank or strlen < 55
		if ((istat = ReadFortranString(line, 55, 4, eth_line_key)) > 0
				&& strcmp(line, "    ") != 0) {
			return(OBS_FILE_SKIP_INPUT_LINE);
		}

		// read phase arrival input
		//CHDAW   P       Q   20.016 1 0.56       10 NWA       0
		//OSS     P       ED  16.524 1 0.11     1672 AHP       0
		//BERNI   P       E   17.930 1 0.18       46 NWA       0

		istat = ReadFortranString(line, 1, 5, arrival->label);
		istat += ReadFortranString(line, 9, 6, arrival->phase);
		istat += ReadFortranString(line, 17, 1, arrival->onset);
		istat += ReadFortranString(line, 18, 1, arrival->first_mot);
                istat += ReadFortranReal(line, 20, 7, &arrival->sec);
//              istat += ReadFortranReal(line, 21, 6, &arrival->sec);
		istat += ReadFortranInt(line, 28, 1, &eth_use_loc);
		istat += ReadFortranReal(line, 29, 5, &arrival->period);
		istat += ReadFortranReal(line, 34, 9, &arrival->amplitude);
		istat += ReadFortranString(line, 44,9, arrival->inst);
		istat += ReadFortranInt(line, 54, 1, &eth_use_mag);


		if (istat != 10) {
			return(OBS_FILE_END_OF_EVENT);
		}

		/* remove blanks/whitespace */
		removeSpace(arrival->label);
		removeSpace(arrival->phase);
		removeSpace(arrival->inst);
		
		if (EvalPhaseID(arrival->phase) < 0) {
			puterr2("WARNING: phase ID not found", arrival->phase);
			return(OBS_FILE_INVALID_PHASE);
		}

		arrival->year = EventTime.year;
		arrival->month = EventTime.month;
		arrival->day = EventTime.day;
		arrival->hour = EventTime.hour;
		arrival->min = EventTime.min;

		/* convert onset to error */
//Uncertainties of phase readings:
//                              I             E            Q
//for Key 'Loca':            < +/- 0.05   < +/- 0.2     > +/- 0.2
//for Key 'Regi' or 'Tele'   < +/- 0.20   < +/- 1.0     > +/- 1.0	
		strcpy(arrival->error_type, "GAU");
		if (eth_use_loc == 0) {
			arrival->error = 9999.9;
		} else if (strcmp(arrival->onset, "I") == 0) {
			arrival->error = 0.05;
		} else if (strcmp(arrival->onset, "E") == 0) {
			arrival->error = 0.2;
		} else if (strcmp(arrival->onset, "Q") == 0) {
			arrival->error = 0.2;
		} else {
			arrival->error = 9999.9;
		}
		arrival->quality = 0;
		
		// check if should use amplitude
		if (eth_use_mag != 0)
			arrival->amplitude == AMPLITUDE_NULL;

		return(istat);
	}

	

	else if (strcmp(ftype_obs, "HYPOCENTER") == 0)
	{

/* example:
   67  827 1632 12.4 L
  ZAK0SZ  PG      1632  33.30
  ZAK0SZ  SG      1632  50.80
  IRK0SZ  PG      1632  29.50
  IRK0SZ  SG      1632  43.90
  KHT0SZ  PG      1632  34.00
*/

		/* read line */
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);

		/* check for end of event (assumes empty lines between event) */
		if (LineIsBlank(line)) {
			/* end of event */
			return(OBS_FILE_END_OF_EVENT);
		}


		/* read event hypocenter line */
		if (!in_hypocenter_event) {
			/* read hypocenter time */
			istat = ReadFortranInt(line, 1, 5, &EventTime.year);
			istat += ReadFortranInt(line, 7, 2, &EventTime.month);
			istat += ReadFortranInt(line, 9, 2, &EventTime.day);
			if (istat == 3) {
				if (EventTime.year < 2000)
					EventTime.year += 1900;
						/* !! assume any 2 digit year is 1900-1999 */
				in_hypocenter_event = 1;
				/* read next line */
				cstat = fgets(line, MAXLINE_LONG, fp_obs);
				if (cstat == NULL)
					return(OBS_FILE_END_OF_INPUT);
			} else {
				return(OBS_FILE_END_OF_EVENT);
			}
		}


		/* read phase arrival input */

		istat = ReadFortranString(line, 3, 3, arrival->label);
		istat += ReadFortranInt(line, 6, 1, &arrival->quality);	// !!! is this field quality???
		istat += ReadFortranString(line, 7, 2, arrival->comp);
		istat += ReadFortranString(line, 11, 6, arrival->phase);
		istat += ReadFortranInt(line, 19, 2, &arrival->hour);
		istat += ReadFortranInt(line, 21, 2, &arrival->min);
		// format of seconds varies, so cannot use fixed fortran read
		//istat += ReadFortranReal(line, 25, 5, &arrival->sec);
		istat += sscanf(line, "%*s %*s %*s %lf", &arrival->sec);

		if (istat != 7) {
			return(OBS_FILE_END_OF_EVENT);
		}

		if (EvalPhaseID(arrival->phase) < 0) {
			puterr2("WARNING: phase ID not found", arrival->phase);
			return(OBS_FILE_INVALID_PHASE);
		}


		arrival->year = EventTime.year;
		arrival->month = EventTime.month;
		arrival->day = EventTime.day;


		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}

	
	
	else if (strcmp(ftype_obs, "HYPODD") == 0)
	{

/* example:
# 1985  1 24  2 19 58.71  37.8832 -122.2415    9.80 1.40  0.15  0.51  0.02      38542
NCCSP       2.850  -1.000   P
NCCBW       3.430  -1.000   P
NCCMC       2.920   0.200   P
NCCAI       3.440  -1.000   P
NCCBR       3.940   1.000   P
NCCLC       4.610   0.500   P
NCJPR       4.470   0.100   P
NCNLH       6.000   0.200   P
NCCSH       6.130   0.100   P
NCNTA       5.630   0.100   P
NCJMG       6.240   0.200   P
*/

		if (in_hypocenter_event) {
			/* check for end of event (assumes no blanks after last phase) */
			chr = fgetc(fp_obs);
			if (chr == EOF) {
				return(OBS_FILE_END_OF_INPUT);
			} else if (chr == '#') {
				// end of event
				ungetc(chr, fp_obs);
				return(OBS_FILE_END_OF_EVENT);
			} else {
				ungetc(chr, fp_obs);
			}
		}
		
		/* read line */
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);


		/* read event hypocenter line */
		if (!in_hypocenter_event) {
			/* read hypocenter time */
			istat = sscanf(line, "# %d %d %d %d %d %lf",
				&EventTime.year, &EventTime.month, &EventTime.day,
				&EventTime.hour, &EventTime.min, &EventTime.sec);
			if (istat == 6) {
				in_hypocenter_event = 1;
				/* read next line */
				cstat = fgets(line, MAXLINE_LONG, fp_obs);
				if (cstat == NULL)
					return(OBS_FILE_END_OF_INPUT);
			} else {
				return(OBS_FILE_END_OF_EVENT);
			}
		}


		// read phase arrival input

		istat = sscanf(line, "%s %lf %lf %s", chrtmp, &ttime, &weight, arrival->phase);
		if (istat != 4)
			return(OBS_FILE_END_OF_EVENT);
			
		// trim label (to agree with Alberto's 3 char convention)
		if (strlen(chrtmp) > 4 && strncmp(chrtmp, "NC", 2) == 0)
			strcpy(arrival->label, chrtmp + 2);
		else
			strcpy(arrival->label, chrtmp);
		
		arrival->hour = EventTime.hour;
		arrival->min = EventTime.min;
		arrival->sec = EventTime.sec + ttime;

		if (EvalPhaseID(arrival->phase) < 0) {
			puterr2("WARNING: phase ID not found", arrival->phase);
			return(OBS_FILE_INVALID_PHASE);
		}

		arrival->year = EventTime.year;
		arrival->month = EventTime.month;
		arrival->day = EventTime.day;

		// convert weight to quality (following Appendix C in hypoDD.pdf)
		weight = fabs(weight);
		if (weight > 0.99)
			arrival->quality = 0;	
		else if (weight > 0.49)
			arrival->quality = 1;	
		else if (weight > 0.19)
			arrival->quality = 2;	
		else if (weight > 0.09)
			arrival->quality = 3;	
		else
			arrival->quality = 4;	
		// convert quality to error
		Qual2Err(arrival);

		return(istat);
	}

	
	
	else if (strcmp(ftype_obs, "RENASS_WWW") == 0)
	{

		/* read next line */
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
 		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);


		/* read formatted P arrival input */
		istat = sscanf(line, "%s %s %dh %dmn %lfsec",
			arrival->label, arrival->phase, &arrival->hour,
			&arrival->min, &arrival->sec);

		if (istat != 5) {
			line[0] = '\0';
			return(OBS_FILE_SKIP_INPUT_LINE);
		}

		arrival->quality = 0;
		arrival->year = 1900;
		arrival->month = 01;
		arrival->day = 01;


		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}

	else if (strcmp(ftype_obs, "RENASS_DEP") == 0)
	{

		/* read next line */
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
 		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);


		/* read formatted arrival input */
		istat = ReadFortranString(line, 1, 4, arrival->label);
		TrimString(arrival->label);
		istat += ReadFortranString(line, 20, 4, arrival->phase);
		TrimString(arrival->phase);
		istat += ReadFortranInt(line, 25, 2, &arrival->hour);
		istat += ReadFortranInt(line, 28, 2, &arrival->min);
		istat += ReadFortranReal(line, 30, 5, &arrival->sec);

		if (istat != 5) {
			line[0] = '\0';
			return(OBS_FILE_SKIP_INPUT_LINE);
		}

		if (EvalPhaseID(arrival->phase) < 0) {
			puterr2("WARNING: phase ID not found", arrival->phase);
			return(OBS_FILE_INVALID_PHASE);
		}


		arrival->quality = 0;
		arrival->year = EventTime.year;
		arrival->month = EventTime.month;
		arrival->day = EventTime.day;
		if (arrival->hour != EventTime.hour)
			puterr(
"WARNING: filename and arrival hours do not match.");
		if (arrival->min != EventTime.min)
			puterr(
"WARNING: filename and arrival minutes do not match.");


		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}

	else if (strcmp(ftype_obs, "SEISAN") == 0 ||
			strcmp(ftype_obs, "NORDIC") == 0)
	{
		/* read next line */
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
 		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);

		/* check for blank line */
		if (LineIsBlank(line))
			return(OBS_FILE_END_OF_EVENT);

		/* check for date info */
		if (!date_saved) {
			istat = ReadFortranInt(line, 2, 4, &year_save);
			if (year_save < 500)
				year_save += 1900;
			istat += ReadFortranInt(line, 7, 2, &month_save);
			istat += ReadFortranInt(line, 9, 2, &day_save);
/*printf(
"SEISAN_DATE Read: istat %d   -  %d %d %d\n",
istat, year_save, month_save, day_save);*/
			if (istat == 3 && IsGoodDate(year_save,
						month_save, day_save)) {
				date_saved = 1;
				return(OBS_FILE_SKIP_INPUT_LINE);
			}
		}


		/* read formatted P arrival input */
		istat = ReadFortranString(line, 2, 4, arrival->label);
		TrimString(arrival->label);
		istat += ReadFortranString(line, 10, 1, arrival->onset);
		TrimString(arrival->onset);
		istat += ReadFortranString(line, 11, 4, arrival->phase);
		TrimString(arrival->phase);
		istat += ReadFortranInt(line, 15, 1, &arrival->quality);
		istat += ReadFortranString(line, 17, 1, arrival->first_mot);
		TrimString(arrival->first_mot);
		istat += ReadFortranInt(line, 19, 2, &arrival->hour);
		istat += ReadFortranInt(line, 21, 2, &arrival->min);
		istat += ReadFortranReal(line, 24, 5, &arrival->sec);
		istat += ReadFortranReal(line, 30, 4, &arrival->coda_dur);
		istat += ReadFortranReal(line, 34, 7, &arrival->amplitude);
		istat += ReadFortranReal(line, 42, 4, &arrival->period);
		if (date_saved) {
			arrival->year = year_save;
			arrival->month = month_save;
			arrival->day = day_save;
		} else {
			arrival->year = 1900;
			arrival->month = 01;
			arrival->day = 01;
		}

/*printf(
"SEISAN Read: istat %d   -  %s %s %s %d %s %d %d %d %d %d %lf %lf\n",
istat, arrival->label, arrival->onset, arrival->phase, arrival->quality, arrival->first_mot, arrival->year, arrival->month, arrival->day, arrival->hour, arrival->min, arrival->sec, arrival->coda_dur);*/

		if (istat != 11)
			return(OBS_FILE_SKIP_INPUT_LINE);

		if (EvalPhaseID(arrival->phase) < 0) {
			puterr2("WARNING: phase ID not found", arrival->phase);
			return(OBS_FILE_INVALID_PHASE);
		}

		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}
	else if (strcmp(ftype_obs, "ETH_LOC") == 0)
	{
/*
 1 1   200.0   300.0    1.71 1                        INST
   0.000     0.000  10.00 0000 00 00 00 00  0.00 0    TRIAL
maraini                                               AUTOR
2001/12/04 01:29                                      Local           20010748
FUORN   P       E   14.910 1 0.14      448 NWA       0
FUORN   S       E   16.140 1 0.14      448 NWA       0
OSS2    P       ID  16.428 1 0.09     2554 AHP       0
OSS     P       ED  16.524 1 0.11     1672 AHP       0
BERNI   P       E   17.930 1 0.18       46 NWA       0
BERNI   S       E   21.139 1 0.18       46 NWA       0
CHDAW   P       Q   20.016 1 0.56       10 NWA       0
VDL2    P       IU  23.704 1 0.20      246 AHP       0
VDL     P       E   23.863 1 0.12      188 AHP       0
VDL     P       E   23.822 1 0.28        4 NWA       1
VDL     S       E   31.177 1 0.28        4 NWA       1
                                                      SKIP
20011204012913146568N010226E00614Ml1239814161001001000SEDG024196006511         B
KP200112040128                                        ARCHIVE
*/
		// read next line
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
 		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);

		// check for blank line
		if (LineIsBlank(line))
			return(OBS_FILE_END_OF_EVENT);
			
		// check for end of phase list
		istat = ReadFortranString(line, 55, 4, chrtmp);
		if (strcmp(chrtmp, "SKIP") == 0 || strcmp(chrtmp, "SEDG") == 0
				|| strcmp(chrtmp, "ARCH") == 0)
			return(OBS_FILE_END_OF_EVENT);
			
		if (!date_saved) {
			// skip header lines
 			while (strcmp(chrtmp, "Loca") != 0) {
				cstat = fgets(line, MAXLINE_LONG, fp_obs);
				istat = ReadFortranString(line, 55, 4, chrtmp);
 			}
			// get date / time info
			istat = ReadFortranInt(line, 1, 4, &year_save);
			istat += ReadFortranInt(line, 6, 2, &month_save);
			istat += ReadFortranInt(line, 9, 2, &day_save);
			istat += ReadFortranInt(line, 12, 2, &hour_save);
			istat += ReadFortranInt(line, 15, 2, &min_save);
printf(
"ETH_LOC_DATE <%s>\n", line);
printf(
"ETH_LOC_DATE Read: istat %d   -  %d/%d/%d %d:%d\n",
istat, year_save, month_save, day_save, hour_save, min_save);
			if (istat == 5 && IsGoodDate(year_save,
						month_save, day_save)) {
				date_saved = 1;
			}
			return(OBS_FILE_SKIP_INPUT_LINE);
		}


		/* read formatted P arrival input */
/*
OSS2    P       ID  16.428 1 0.09     2554 AHP       0
*/
		istat = ReadFortranString(line, 1, 4, arrival->label);
		TrimString(arrival->label);
		istat += ReadFortranString(line, 5, 1, arrival->comp);
		TrimString(arrival->comp);
		istat += ReadFortranString(line, 9, 7, arrival->phase);
		TrimString(arrival->phase);
		istat += ReadFortranString(line, 17, 1, arrival->onset);
		TrimString(arrival->onset);
		istat += ReadFortranString(line, 18, 1, arrival->first_mot);
		TrimString(arrival->first_mot);
		istat += ReadFortranReal(line, 21, 6, &arrival->sec);
		//istat += ReadFortranInt(line, 15, 1, &arrival->quality);
arrival->quality = 0;
		//istat += ReadFortranReal(line, 30, 4, &arrival->coda_dur);
		//istat += ReadFortranReal(line, 34, 7, &arrival->amplitude);
		//istat += ReadFortranReal(line, 42, 4, &arrival->period);
		if (date_saved) {
			arrival->year = year_save;
			arrival->month = month_save;
			arrival->day = day_save;
			arrival->hour = hour_save;
			arrival->min = min_save;
		} else {
			arrival->year = 1900;
			arrival->month = 01;
			arrival->day = 01;
			arrival->hour = 00;
			arrival->min = 00;
		}

printf(
"ETH_LOC Read: istat %d   -  %s %s %s %d %s %d %d %d %d %d %lf %lf\n",
istat, arrival->label, arrival->onset, arrival->phase, arrival->quality, arrival->first_mot,
arrival->year, arrival->month, arrival->day, arrival->hour, arrival->min, arrival->sec,
arrival->coda_dur);


		if (istat != 6)
			return(OBS_FILE_SKIP_INPUT_LINE);

		if (EvalPhaseID(arrival->phase) < 0) {
			puterr2("WARNING: phase ID not found", arrival->phase);
			return(OBS_FILE_INVALID_PHASE);
		}

		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}
	else
	{
		puterr2("ERROR: unrecognized observation file type:", ftype_obs);
		return(OBS_FILE_END_OF_INPUT);
	}
}


/*** function remove whitespace from string */

void removeSpace(char *str) {

	int n, m;

	n = 0;
	while (str[n] != '\0' && n < 1000000)
		if (isspace(str[n])) {
			m = n;
			while (str[m] != '\0' && m < 1000000) {
				str[m] = str[m + 1];
				m++;
			}
		} else {
			n++;
		}

}


/*** function to check if phase ID is in phase ID list for given phase */

int IsPhaseID(char *phase_in, char *phase_check)
{
	int npha_id;
	char test_str[PHASE_LABEL_LEN + 4];

//printf("phase_in %s, phase_check %s, ", phase_in, phase_check);
	/* check for blank phase_in */
	if (strstr("              ", phase_in) != NULL)
//{printf(" NO - BLANK\n");
		return(0);
//}

	removeSpace(phase_in);

	sprintf(test_str, " %s ", phase_in);
	for (npha_id = 0; npha_id < NumPhaseID; npha_id++) {
//printf(" compare (<%s> <%s>)\n", PhaseID[npha_id].id_string, test_str);
		if (strcmp(PhaseID[npha_id].phase, phase_check) == 0) {
			if (strstr(PhaseID[npha_id].id_string, test_str) != NULL)
//{printf(" YES (%s %s)\n", PhaseID[npha_id].id_string,  phase_in);
				return(1);
//}
		}
	}

//printf(" NO\n");
	return(0);
}



/*** function to evaluate phase ID from phase ID lists */

int EvalPhaseID(char *phase_in)
{
	int npha_id;

	for (npha_id = 0; npha_id < NumPhaseID; npha_id++) {
		if (IsPhaseID(phase_in, PhaseID[npha_id].phase)) {
			//strcpy(phase_in, PhaseID[npha_id].phase);
			return(0);
		}
	}

	return(-1);
}


/*** function to check if a date is reasonable */

int IsGoodDate(int iyear, int imonth, int iday)
{
	if (iyear >= SMALLEST_EVENT_YEAR && iyear <= LARGEST_EVENT_YEAR
			&& imonth > 0 && imonth < 13
			&& iday > 0 && iday < 32)
		return(1);

	return(0);
}



/*** function to homogenize date / time of arrivals */

int HomogDateTime(ArrivalDesc *arrival, int num_arrivals, HypoDesc* phypo)
{
	int narr;
	int dofymin = 10000, yearmin = 10000;
	int test_month, test_day;

	for (narr = 0; narr < num_arrivals; narr++) {
		arrival[narr].day_of_year =
			DayOfYear(arrival[narr].year, arrival[narr].month,
					arrival[narr].day);
		if (arrival[narr].day_of_year < dofymin)
			dofymin = arrival[narr].day_of_year;
		if (arrival[narr].year < yearmin)
			yearmin = arrival[narr].year;
		if (arrival[narr].year != yearmin) {
//			puterr(
//"ERROR: arrivals cross year boundary, ignoring observation set.");
			return(OBS_FILE_ARRIVALS_CROSS_YEAR_BOUNDARY);
		}
	}

	for (narr = 0; narr < num_arrivals; narr++) {
		if (arrival[narr].day_of_year > dofymin) {
			arrival[narr].day_of_year--;
			arrival[narr].day--;
			arrival[narr].hour += 24;
		}
	}

	for (narr = 0; narr < num_arrivals; narr++)
		arrival[narr].obs_time = (long double) arrival[narr].sec
// DELAY_CORR
			- (long double) arrival[narr].delay
			+ 60.0L * ((long double) arrival[narr].min
			+ 60.0L * (long double) arrival[narr].hour);

	if (!FixOriginTimeFlag) {
		/* initialize hypocenter year/month/day
						if origin time not fixed */
		phypo->year = yearmin;
		MonthDay(yearmin, dofymin, &(phypo->month), &(phypo->day));
	} else {
		/* homogenize hypocenter otime if origin time fixed */
		MonthDay(yearmin, dofymin, &test_month, &test_day);
		if (phypo->year != yearmin || test_month != phypo->month
						|| test_day != phypo->day) {
			puterr(
"ERROR: earliest arrivals year/month/day does not match fixed origin time year/month/day, ignoring observation set.");
			return(OBS_FILE_ARRIVALS_CROSS_YEAR_BOUNDARY);
		}
		phypo->time = (long double) phypo->sec
			+ 60.0L * ((long double) phypo->min
			+ 60.0L * (long double) phypo->hour);
		phypo->min = 0;
		phypo->hour = 0;
	}

	return(0);


}




/*** function to standardize date / time of arrivals */

int StdDateTime(ArrivalDesc *arrival, int num_arrivals, HypoDesc* phypo)
{
	int narr;
	int dofymin = 10000, yearmin = 10000;
	double rms_resid = 0.0, weight_sum = 0.0;
	long double sec_tmp, hyp_time_tmp;


	for (narr = 0; narr < num_arrivals; narr++) {
		/* calc obs travel time and residual */
		arrival[narr].obs_travel_time =
			arrival[narr].obs_time - phypo->time;
		arrival[narr].residual = arrival[narr].obs_travel_time -
			arrival[narr].pred_travel_time;
		rms_resid += arrival[narr].weight *
			arrival[narr].residual * arrival[narr].residual;
		weight_sum += arrival[narr].weight;
		/* convert time to year/month/day/hour/min */
// DELAY_CORR		sec_tmp = arrival[narr].obs_time;
		sec_tmp = arrival[narr].obs_time + arrival[narr].delay;
		arrival[narr].hour = (int) (sec_tmp / 3600.0L);
		sec_tmp -= (long double) arrival[narr].hour * 3600.0L;
		arrival[narr].min = (int) (sec_tmp / 60.0L);
		sec_tmp -= (long double) arrival[narr].min * 60.0L;
		arrival[narr].sec = (double) sec_tmp;
		MonthDay(arrival[narr].year, arrival[narr].day_of_year,
			&(arrival[narr].month), &(arrival[narr].day));
	}

	phypo->rms = sqrt(rms_resid / weight_sum);

	hyp_time_tmp = phypo->time;
	phypo->hour = (int) (hyp_time_tmp / 3600.0L);
	hyp_time_tmp -= (long double) phypo->hour * 3600.0L;
	phypo->min = (int) (hyp_time_tmp / 60.0L);
	hyp_time_tmp -= (long double) phypo->min * 60.0L;
	phypo->sec = (double) hyp_time_tmp;

}




/*** function to set output file root name using arrival time */

int SetOutName(ArrivalDesc *arrival, char* out_file_root, char* out_file,
	int isec)
{

	if (isec)
		sprintf(out_file, "%s.%4.4d%2.2d%2.2d.%2.2d%2.2d%2.2d",
		out_file_root, arrival->year, arrival->month, arrival->day,
		arrival->hour, arrival->min, (int) arrival->sec);
	else
		sprintf(out_file, "%s.%4.4d%2.2d%2.2d.%2.2d%2.2d",
		out_file_root, arrival->year, arrival->month, arrival->day,
		arrival->hour, arrival->min);

}


/*** function to check for duplicate label and phase in arrival */

int IsDuplicateArrival(ArrivalDesc *arrival, int num_arrivals, int ntest,
			char *phase_test)
{
	int narr;

	if (phase_test == NULL) {
		for (narr = 0; narr < num_arrivals; narr++) {
			if (narr != ntest
					&& (IsPhaseID(arrival[narr].phase, "P") &&
							IsPhaseID(arrival[ntest].phase, "P")
						|| IsPhaseID(arrival[narr].phase, "S") &&
							IsPhaseID(arrival[ntest].phase, "S"))
					&& !strcmp(arrival[narr].time_grid_label,
						arrival[ntest].time_grid_label))
				return(narr);
		}
	} else {
		for (narr = 0; narr < num_arrivals; narr++) {
			if (narr != ntest
					&& !strcmp(arrival[narr].time_grid_label,
						arrival[ntest].time_grid_label)
					&& IsPhaseID(arrival[narr].phase, phase_test))
				return(narr);
		}
	}

	return(-1);

}



/*** function to perform grid search location */

int LocGridSearch(int ngrid, int num_arr_total, int num_arr_loc,
		ArrivalDesc *arrival,
		GridDesc* ptgrid, GaussLocParams* gauss_par, HypoDesc* phypo)
{

	int istat;
	int ix, iy, iz, narr, n_compan;
	int iGridType;
	int nReject, numGridReject = 0, numStaReject = 0;
	double xval, yval, zval;
	double yval_grid;
	/*double travel_time;*/
	double value;
	double misfit;
	double misfit_min = VERY_LARGE_DOUBLE, misfit_max = -VERY_LARGE_DOUBLE;
	double dlike;





	/* get solution quality at each grid point */

	putmsg(3, "");
	putmsg(3, "Calculating solution over grid...");

	iGridType = ptgrid->type;

	xval = ptgrid->origx;

	/* loop over grid points */

	for (ix = 0; ix <  ptgrid->numx; ix++) {

		/* read y-z sheets for arrival travel-times (3D grids) */
		if ((istat = ReadArrivalSheets(num_arr_loc, arrival,
					xval)) < 0)
			puterr("ERROR: reading arrival travel time sheets.");

		yval = ptgrid->origy;
		for (iy = 0; iy <  ptgrid->numy; iy++) {
			zval = ptgrid->origz;
			for (iz = 0; iz <  ptgrid->numz; iz++) {


			    /* loop over observed arrivals */

			    nReject = 0;
			    for (narr = 0; narr < num_arr_loc; narr++) {
				/* check for companion */
				if ((n_compan = arrival[narr].n_companion) >= 0) {
					if ((arrival[narr].pred_travel_time =
						arrival[n_compan].pred_travel_time) < 0.0)
					    nReject++;
				/* else check grid type */
				} else if (arrival[narr].sheetdesc.type
							== GRID_TIME) {
					/* 3D grid */
					if ((arrival[narr].pred_travel_time =
					    (double) ReadAbsInterpGrid3d(
						NULL,
						&(arrival[narr].sheetdesc),
						xval, yval, zval)) < 0.0)
					nReject++;
				} else {
					/* 2D grid (1D model) */
					yval_grid = GetEpiDist(
						&(arrival[narr].station),
						xval, yval);
					if ((arrival[narr].pred_travel_time =
					    ReadAbsInterpGrid2d(
						NULL,
						&(arrival[narr].sheetdesc),
						yval_grid, zval)) < 0.0)
					nReject++;
				}
				arrival[narr].pred_travel_time *= arrival[narr].tfact;
			    }

			    if (nReject) {
					numGridReject++;
					numStaReject += nReject;
					misfit = -1.0;
					if (iGridType == GRID_MISFIT)
						value = -1.0;
					else if (iGridType == GRID_PROB_DENSITY)
						value = 0.0;
			    } else {

				/* calc misfit or prob density */

				value = CalcSolutionQuality(num_arr_loc,
						arrival, gauss_par,
						iGridType, &misfit);
				ptgrid->array[ix][iy][iz] = value;
				if (iGridType == GRID_MISFIT) {
					ptgrid->sum += value;
				} else if (iGridType == GRID_PROB_DENSITY) {
					dlike = exp(value);
					ptgrid->sum += dlike;
					/* update  probabilitic residuals */
					UpdateProbabilisticResiduals(
						num_arr_loc, arrival, dlike);
				}

				/* check for minimum misfit */
				if (misfit < misfit_min) {
					misfit_min = misfit;
					phypo->misfit = misfit;
					phypo->ix = ix;
					phypo->iy = iy;
					phypo->iz = iz;
					phypo->x = xval;
					phypo->y = yval;
					phypo->z = zval;
			 	 	for (narr = 0; narr < num_arr_loc; narr++)
						arrival[narr].pred_travel_time_best =
							arrival[narr].pred_travel_time;
				}
				if (misfit > misfit_max)
					misfit_max = misfit;

			    }

				zval += ptgrid->dz;
			}
			yval += ptgrid->dy;
		}
		xval += ptgrid->dx;
	}


	/* give warning if grid points rejected */

	if (numGridReject > 0) {
		sprintf(MsgStr, "WARNING: %d grid locations rejected; travel times for an average of %.2lf arrival observations were not valid.",
			numGridReject, (double) numStaReject / numGridReject);
		putmsg(1, MsgStr);
	}


	/* construct search information string */
	sprintf(phypo->searchInfo,
		"GRID nPts %d\0", ix * iy *iz);
	/* write message */
	/*putmsg(1, phypo->searchInfo);*/


	/* re-calculate solution and arrival statistics for best location */

	SaveBestLocation(num_arr_total, num_arr_loc, arrival,  ptgrid,
		gauss_par, phypo,
		misfit_max, iGridType);


	return(0);

}



/*** function to perform Metropolis location */

#define MAX_NUM_MET_TRIES 1000

int LocMetropolis(int ngrid, int num_arr_total, int num_arr_loc,
		ArrivalDesc *arrival,
		GridDesc* ptgrid, GaussLocParams* gauss_par, HypoDesc* phypo,
		WalkParams* pMetrop, float* fdata)
{

	int istat;
	int ntry, nSamples, nSampStat, narr, ipos;
	long ngenerated;
	int maxNumTries;
	int writeMessage = 0;
	int iGridType;
	int nReject, numClipped = 0, numGridReject = 0, numStaReject = 0;
	int iAbort = 0, iReject = 0;
	int iAccept, numAcceptDeepMinima = 0;
	double xval, yval, zval;
	double currentMetStepFact;

	double value, dlike, dlike_max = -VERY_LARGE_DOUBLE;

	double misfit;
	double misfit_min = VERY_LARGE_DOUBLE, misfit_max = -VERY_LARGE_DOUBLE;

	double xmin, xmax, ymin, ymax, zmin, zmax;
	double dx_init, dx_test;

	int nScatterSaved;

	Vect3D expect = {0.0, 0.0, 0.0};

	double xmean_sum = 0.0, ymean_sum = 0.0, zmean_sum = 0.0;
	double xvar_sum = 0.0, yvar_sum = 0.0, zvar_sum = 0.0;
	double xvar = 0.0, yvar = 0.0, zvar = 0.0;
	double dsamp, dsamp2;



	/* get solution quality at each sample on random walk */

	putmsg(3, "");
	putmsg(3, "Calculating solution along Metropolis walk...");

	iGridType = GRID_PROB_DENSITY;

	/* set walk limits equal to grid limits */
	xmin = ptgrid->origx;
	xmax = xmin + (double) (ptgrid->numx - 1) * ptgrid->dx;
	ymin = ptgrid->origy;
	ymax = ymin + (double) (ptgrid->numy - 1) * ptgrid->dy;
	zmin = ptgrid->origz;
	zmax = zmin + (double) (ptgrid->numz - 1) * ptgrid->dz;

	/* save intiial values */
	currentMetStepFact = MetStepFact;
	dx_init = pMetrop->dx;


	/* loop over walk samples */

	nSamples = 0;
	nSampStat = 0;
	nScatterSaved = 0;
	ipos = 0;
	ntry = 0;
	ngenerated = 0;
 	maxNumTries = MAX_NUM_MET_TRIES;
	while (nSamples < MetNumSamples
			&& (nSamples <= MetLearn || ntry < maxNumTries)) {

		ntry++;
		ngenerated++;
		istat = GetNextMetropolisSample(pMetrop,
				xmin, xmax, ymin, ymax,
				zmin, zmax, &xval, &yval, &zval);
		if (nSamples > MetEquil && istat > 0)
			numClipped += istat;

 		/* get travel times for observed arrivals */

		nReject = getTravelTimes(arrival, num_arr_loc, xval, yval, zval);


		if (nReject) {
			numGridReject++;
			numStaReject += nReject;
			misfit = -1.0;
			dlike = 0.0;
		} else {

			/* calc misfit or prob density */
			value = CalcSolutionQuality(num_arr_loc,
					arrival, gauss_par,
					iGridType, &misfit);
			dlike = gauss_par->WtMtrxSum * exp(value);

			/* apply Metropolis test */
			iAccept = MetropolisTest(pMetrop, dlike);

			/* if not accepted, but at maxNumTries... */
			if (!iAccept && ntry == maxNumTries) {
				/* if not learning, accept anyway since
						may be stuck in a deep minima */
				if (nSamples >= MetLearn
						&& numAcceptDeepMinima++ < 5) {
					iAccept = 1;
//printf("Max Num Tries: accept deep minima\n");

					/* try reducing step size */
					currentMetStepFact /= 2.0;
					ntry = 0;
//printf("            +: step ch: was %lf\n", pMetrop->dx);

					/*if (pMetrop->dx > MetStepMin) {
						pMetrop->dx /= 2.0;
//printf("            +: step ch: %lf -> %lf\n", 2.0 * pMetrop->dx, pMetrop->dx);
						ntry = 0;
					}*/

				/* if learning, try reducing step size */
				} else if (nSamples < MetLearn
						&& pMetrop->dx > MetStepMin) {
					pMetrop->dx /= 2.0;
//printf("Max Num Tries: step ch: %lf -> %lf\n", 2.0 * pMetrop->dx, pMetrop->dx);
					ntry = 0;
				}
			}


			if (iAccept) {

				ntry = 0;
				nSamples++;

				/* check for minimum misfit */
				if (misfit < misfit_min) {
					misfit_min = misfit;
					dlike_max = dlike;
					phypo->misfit = misfit;
					phypo->x = xval;
					phypo->y = yval;
					phypo->z = zval;
				  	for (narr = 0; narr < num_arr_loc; narr++)
						arrival[narr].pred_travel_time_best =
							arrival[narr].pred_travel_time;
				}
				if (misfit > misfit_max)
					misfit_max = misfit;

				/* update sample location */
				pMetrop->x = xval;
				pMetrop->y = yval;
				pMetrop->z = zval;
				pMetrop->likelihood = dlike;

				/* if learning, update sample statistics */
				if (nSamples > MetLearn / 2
					&& nSamples <= MetLearn + MetEquil) {

					xmean_sum += xval;
					ymean_sum += yval;
					zmean_sum += zval;
					xvar_sum += xval * xval;
					yvar_sum += yval * yval;
					zvar_sum += zval * zval;
					nSampStat++;
				}

				/* if equilibrating, update Met step */
				if (nSamples > MetLearn
					&& nSamples <= MetLearn + MetEquil) {

					/* update Met step */
					dsamp = (double) nSampStat;
					dsamp2 = dsamp * dsamp;
					xvar = xvar_sum  / dsamp -
					    xmean_sum * xmean_sum / dsamp2;
					yvar = yvar_sum  / dsamp -
					    ymean_sum * ymean_sum / dsamp2;
					zvar = zvar_sum  / dsamp -
					    zmean_sum * zmean_sum / dsamp2;
					dx_test = currentMetStepFact * pow(
					   sqrt(xvar) * sqrt(yvar) * sqrt(zvar)
		/*/ (double) (MetUse / MetSkip),*/
						/ (double) MetUse,
 						1.0/3.0);
//if (pMetrop->dx != dx_test) printf("equil step ch: %lf -> %lf\n", pMetrop->dx, dx_test);

					if (dx_test > MetStepMin)
						pMetrop->dx = dx_test;
					else
						pMetrop->dx = MetStepMin;
				}

				/* if saving samples */
				if (nSamples > MetStartSave
						&& nSamples % MetSkip == 0) {

					/* save sample to scatter file */
					fdata[ipos++] = xval;
					fdata[ipos++] = yval;
					fdata[ipos++] = zval;
					fdata[ipos++] = dlike;

					/* update  probabilitic residuals */
					if (1)
						UpdateProbabilisticResiduals(
							num_arr_loc, arrival, 1.0);


					nScatterSaved++;
				}

				if (nSamples % 1000 == 1
					|| nSamples == MetLearn / 2)
					writeMessage = 1;

			}


if (writeMessage || ntry == maxNumTries - 1) {
		sprintf(MsgStr,
"Metropolis: n %d x %.2lf y %.2lf z %.2lf  xm %.2lf ym %.2lf zm %.2lf  xdv %.2lf ydv %.2lf zdv %.2lf  dx %.2lf  li %.2le", nSamples, pMetrop->x, pMetrop->y, pMetrop->z, xmean_sum / dsamp, ymean_sum / dsamp, zmean_sum / dsamp, sqrt(xvar), sqrt(yvar), sqrt(zvar), pMetrop->dx, pMetrop->likelihood);
		putmsg(3, MsgStr);
		writeMessage = 0;
}

		}


		/* check abort search conditions */

		/* failure to accept sample after maxNumTries */
		if (nSamples > MetLearn && ntry >= maxNumTries) {
			sprintf(MsgStr,
"ERROR: failed to accept new Metropolis sample after %d tries, aborting location.", ntry);
			puterr(MsgStr);
			sprintf(phypo->locStatComm, "%s", MsgStr);
			iAbort = 1;
			break;
		}

		/* maximum likelihood too low after learning stage */
		if (nSamples == MetLearn && dlike_max < MetProbMin) {
			sprintf(MsgStr,
"ERROR: after learning stage (%d samples), best probability = %.2le is less than ProbMin = %.2le, aborting location.",
				MetLearn, dlike_max, MetProbMin);
			puterr(MsgStr);
			sprintf(phypo->locStatComm, "%s", MsgStr);
			iAbort = 1;
			break;
		}

	}


	/* give warning if sample points clipped */

	if (numClipped > 0) {
		sprintf(MsgStr,
"WARNING: %d Metropolis samples clipped at search grid boundary.",
			numClipped);
		putmsg(1, MsgStr);
	}


	/* give warning if grid points rejected */

	if (numGridReject > 0) {
		sprintf(MsgStr,
"WARNING: %d Metropolis samples rejected; travel times for an average of %.2lf arrival observations were not valid.",
			numGridReject, (double) numStaReject / numGridReject);
		putmsg(1, MsgStr);
	}


	/* check reject location conditions */

	/* maximum like hypo on edge of grid */
	if (isOnGridBoundary(phypo->x, phypo->y, phypo->z,
			ptgrid, pMetrop->dx)) {
		sprintf(MsgStr,
"WARNING: max prob location on grid boundary, rejecting location.");
		putmsg(1, MsgStr);
		sprintf(phypo->locStatComm, "%s", MsgStr);
		iReject = 1;
	}

	/* construct search information string */
	sprintf(phypo->searchInfo,
"METROPOLIS nSamp %ld nAcc %d nSave %d nClip %d Dstep0 %lf Dstep %lf\0",
		ngenerated, nSamples, nScatterSaved, numClipped, dx_init, pMetrop->dx);
	/* write message */
	putmsg(1, phypo->searchInfo);


	/* check for termination */
	if (iAbort) {
		sprintf(Hypocenter.locStat, "ABORTED");
	} else if (iReject) {
		sprintf(Hypocenter.locStat, "REJECTED");
	}


	/* re-calculate solution and arrival statistics for best location */

	SaveBestLocation(num_arr_total, num_arr_loc, arrival,  ptgrid,
		gauss_par, phypo, misfit_max, iGridType);


	return(nScatterSaved);

}




/*** function to create next metropolis sample */
/* move sample random distance and direction */

#define NTRY_MAX 10

inline int GetNextMetropolisSample(WalkParams* pMetrop, double xmin, double xmax,
	double ymin, double ymax, double zmin, double zmax,
	double* pxval, double* pyval, double* pzval)
{

	int iClip = 0;
	double valx, valy, valz, valsum, norm;
	double x, y, z;


	/* get unit vector in random direction */

	do {
		valx = get_rand_double(-1.0, 1.0);
		valy = get_rand_double(-1.0, 1.0);
		valz = get_rand_double(-1.0, 1.0);
		valsum = valx * valx + valy * valy + valz * valz;
	} while (valsum < SMALL_DOUBLE);

	norm = pMetrop->dx / sqrt(valsum);

	/* add step to last sample location */

	x = pMetrop->x + norm * valx;
	y = pMetrop->y + norm * valy;
	z = pMetrop->z + norm * valz;


	/* crude clip against grid boundary */
	/* clip needed because travel time lookup requires that
			location is within initial search grid */
	if (x < xmin) {
		x = xmin;
		iClip = 1;
	} else if (x > xmax) {
		x = xmax;
		iClip = 1;
	}
	if (y < ymin) {
		y = ymin;
		iClip = 1;
	} else if (y > ymax) {
		y = ymax;
		iClip = 1;
	}
	if (z < zmin) {
		z = zmin;
		iClip = 1;
	} else if (z > zmax) {
		z = zmax;
		iClip = 1;
	}


	/* update sample location */

	*pxval = x;
	*pyval = y;
	*pzval = z;

	return(iClip);

}


/*** function to test new metropolis string */

inline int MetropolisTest(WalkParams* pMetrop, double likelihood_new)
{

	int accept;
	double prob;

	/* compare with last sample using Mosegaard & Tarantola eq (17) */

	if (likelihood_new >= pMetrop->likelihood)
		return(1);
	else if ((prob = get_rand_double(0.0, 1.0))
			< likelihood_new / pMetrop->likelihood)
		return(1);
	else
		return(0);

}


/** function to re-calculate solution and arrival statistics for best location */
/* some quantities are calculated only for arrivals used in location
		(num_arr_loc) others for all arrivals (num_arr_total) */

int SaveBestLocation(int num_arr_total, int num_arr_loc, ArrivalDesc *arrival,
		GridDesc* ptgrid, GaussLocParams* gauss_par, HypoDesc* phypo,
		double misfit_max, int iGridType)
{
	int istat, narr, n_compan, iopened;
	char filename[FILENAME_MAX];
	double value, misfit;


//printf("SaveBestLocation num_arr_total %d num_arr_loc %d\n", num_arr_total, num_arr_loc);

	/* loop over observed arrivals */
	for (narr = 0; narr < num_arr_total; narr++) {
		arrival[narr].dist = GetEpiDist(
				&(arrival[narr].station), phypo->x, phypo->y);
		arrival[narr].azim = GetEpiAzim(
				&(arrival[narr].station), phypo->x, phypo->y);

		/* get best travel time */

		arrival[narr].pred_travel_time = 0.0;
		n_compan = arrival[narr].n_companion;
		iopened = 0;
		/* check for stored best travel time */
		if (arrival[narr].pred_travel_time_best > 0.0) {
			/* load stored best travel time */
			arrival[narr].pred_travel_time = arrival[narr].pred_travel_time_best;
		/* check for companion travel time */
		} else if (n_compan >= 0 && arrival[n_compan].pred_travel_time_best > 0.0) {
			/* load companion stored best travel time */
			arrival[narr].pred_travel_time =
					arrival[n_compan].pred_travel_time_best;
			arrival[narr].pred_travel_time *= arrival[narr].tfact;
		} else {
			/* temporarily open time grid file
					and read time for ignored arrivals */
			sprintf(filename, "%s.time", arrival[narr].fileroot);
			if ((istat = OpenGrid3dFile(filename,
					&(arrival[narr].fpgrid),
					&(arrival[narr].fphdr),
					&(arrival[narr].gdesc), "time",
					&(arrival[narr].station))) < 0)
				continue;
			//iopened = 1;
			/* check grid type, read travel time */
			if (arrival[narr].gdesc.type == GRID_TIME) {
				/* 3D grid */
				if (arrival[narr].fpgrid != NULL)
				    arrival[narr].pred_travel_time =
					(double) ReadAbsInterpGrid3d(
					arrival[narr].fpgrid, &(arrival[narr].gdesc),
					phypo->x, phypo->y, phypo->z);
				if (arrival[narr].pred_travel_time < -LARGE_DOUBLE)
					arrival[narr].pred_travel_time = 0.0;
			} else {
				/* 2D grid (1D model) */
				if (arrival[narr].fpgrid != NULL)
				    arrival[narr].pred_travel_time =
					ReadAbsInterpGrid2d(
					arrival[narr].fpgrid, &(arrival[narr].gdesc),
					arrival[narr].dist, phypo->z);
				if (arrival[narr].pred_travel_time < -LARGE_DOUBLE)
					arrival[narr].pred_travel_time = 0.0;
			}
			CloseGrid3dFile(&(arrival[narr].fpgrid),
				&(arrival[narr].fphdr));
			arrival[narr].pred_travel_time *= arrival[narr].tfact;

		}

		/* read angles */
		/* angle grid file name */
		if (n_compan >= 0)
			sprintf(filename, "%s.angle", arrival[n_compan].fileroot);
		else
			sprintf(filename, "%s.angle", arrival[narr].fileroot);
		if (angleMode == ANGLE_MODE_YES) {
			if (arrival[narr].gdesc.type == GRID_TIME) {
				/* 3D grid */
				ReadTakeOffAnglesFile(filename,
					phypo->x, phypo->y, phypo->z,
					&(arrival[narr].ray_azim),
					&(arrival[narr].ray_dip),
					&(arrival[narr].ray_qual), -1.0);
			} else {
				/* 2D grid (1D model) */
				ReadTakeOffAnglesFile(filename,
					0.0, arrival[narr].dist, phypo->z,
					&(arrival[narr].ray_azim),
					&(arrival[narr].ray_dip),
					&(arrival[narr].ray_qual), arrival[narr].azim);
			}
		}

//		/* close time grid file for ignored arrivals */
//		if (iopened)
//			CloseGrid3dFile(&(arrival[narr].fpgrid),
//				&(arrival[narr].fphdr));
	}

	/* calc misfit or prob density */
	value = CalcSolutionQuality(num_arr_loc, arrival, gauss_par,
			iGridType, &misfit);

	/* calculate maximum likelihood origin time */
	if (!FixOriginTimeFlag)
		phypo->time =
			CalcMaxLikeOriginTime(num_arr_loc, arrival, gauss_par);
	if (iGridType == GRID_PROB_DENSITY) {
		phypo->probmax = exp(value);
	}

	/* set misc hypo fields */
	phypo->grid_misfit_max = misfit_max;
	istat = rect2latlon(0, phypo->x, phypo->y,
		&(phypo->dlat), &(phypo->dlong));
	phypo->depth = phypo->z;
	phypo->nreadings = num_arr_loc;

	return(0);

}



/*** function to read y-z travel time sheet from disk for each arrival */

int ReadArrivalSheets(int num_arrivals, ArrivalDesc *arrival, double xsheet)
{

	int istat, narr, ixsheet;
	float **array_tmp;
	double sheet_origx, sheet_dx;


	/* loop over arrivals */

	for (narr = 0; narr < num_arrivals; narr++) {

		/* skip sheet read if arrival has companion */
		if (arrival[narr].n_companion >= 0)
				continue;

		/* skip sheet read or set xsheet to zero for 2D grid */
		if (arrival[narr].gdesc.type == GRID_TIME_2D) {
			if (arrival[narr].sheetdesc.origx < LARGE_DOUBLE)
				continue;
			xsheet = 0.0;
		}

		sheet_origx = arrival[narr].sheetdesc.origx;
		sheet_dx = arrival[narr].sheetdesc.dx;


		/* check which sheets are required from disc */

		/* both required sheets already read */
		if (sheet_origx <= xsheet && xsheet < sheet_origx + sheet_dx)
			continue;

		/* find x index in disk grid of lower plane of dual sheet */
		ixsheet = (int) ((xsheet - arrival[narr].gdesc.origx)
						/ arrival[narr].gdesc.dx);
		if (ixsheet < 0 || ixsheet > arrival[narr].gdesc.numx - 1) {
			puterr("WARNING: invalid ixsheet value:");
			sprintf(MsgStr, "  Arr: %d  ixsheet: %d", narr, ixsheet);
			puterr(MsgStr);
		}

		/* one required sheet already read */
		if (sheet_origx + sheet_dx <= xsheet &&
				xsheet < sheet_origx + 2.0 * sheet_dx)
		{
			/* exchange sheet pointers */
			array_tmp =  arrival[narr].sheetdesc.array[0];
			arrival[narr].sheetdesc.array[0] =
				arrival[narr].sheetdesc.array[1];
			arrival[narr].sheetdesc.array[1] = array_tmp;

			/* read next sheet if xsheet not exactly on last sheet */
/*			if (fabs(xsheet - (sheet_origx + sheet_dx)) */
/*					> VERY_SMALL_DOUBLE) { */
				/* read new sheet */
				if ((istat =
					ReadGrid3dBufSheet(
					arrival[narr].sheetdesc.array[1][0],
					&(arrival[narr].gdesc),
					arrival[narr].fpgrid, ixsheet + 1)) < 0)
						puterr(
"ERROR: reading new arrival travel time sheet.");
/*			} */

			/* set dual-sheet origin */
			arrival[narr].sheetdesc.origx += sheet_dx;
		}

		/* no required sheets already read */
		else
		{

			/* read lower sheet */
			if ((istat =
				ReadGrid3dBufSheet(
				arrival[narr].sheetdesc.array[0][0],
				&(arrival[narr].gdesc),
				arrival[narr].fpgrid, ixsheet)) < 0)
					puterr(
"ERROR: reading lower arrival travel time sheet.");

			/* read upper sheet if not at last sheet */
			if (ixsheet + 1 < arrival[narr].gdesc.numx) {
				if ((istat =
					ReadGrid3dBufSheet(
					arrival[narr].sheetdesc.array[1][0],
					&(arrival[narr].gdesc),
					arrival[narr].fpgrid, ixsheet + 1)) < 0)
						puterr(
"ERROR: reading upper arrival travel time sheet.");
			}

			/* set dual-sheet origin */
			arrival[narr].sheetdesc.origx =
				(double) ixsheet * sheet_dx
					+ arrival[narr].gdesc.origx;
		}

/*Narr %d O %lf %lf %lf  N %d %d %d  dx %lf %lf %lf\n", narr, arrival[narr].sheetdesc.origx, arrival[narr].sheetdesc.origy, arrival[narr].sheetdesc.origz, arrival[narr].sheetdesc.numx, arrival[narr].sheetdesc.numy, arrival[narr].sheetdesc.numz, arrival[narr].sheetdesc.dx, arrival[narr].sheetdesc.dy, arrival[narr].sheetdesc.dz);*/
	}

}




/*** function to consruct weight matrix (inverse of covariance matrix) */

int ConstWeightMatrix(int num_arrivals, ArrivalDesc *arrival,
		GaussLocParams* gauss_par)
{
	int istat, nrow, ncol;
	double sigmaT2, corr_len2;
	double dx, dy, dz, dist2;
	double weight_sum;
	DMatrix cov_matrix, null_mtrx = NULL;
	SourceDesc *sta1, *sta2;


	/* allocate square matrix  */

	cov_matrix = dmatrix(0, num_arrivals - 1, 0, num_arrivals - 1);


	/* set constants */

	sigmaT2 = gauss_par->SigmaT * gauss_par->SigmaT;
	corr_len2 = gauss_par->CorrLen * gauss_par->CorrLen;
	if (corr_len2 < VERY_SMALL_DOUBLE)
		corr_len2 = 1.0;


	/* load covariances */

	for (nrow = 0; nrow < num_arrivals; nrow++) {
		sta1 = &(arrival[nrow].station);
		for (ncol = 0; ncol <= nrow; ncol++) {
		    sta2 = &(arrival[ncol].station);

		    /* travel time error (TV82, eq. 10-14; MEN92, eq. 22) */
		    if (strcmp(arrival[nrow].phase, arrival[ncol].phase) == 0) {
			dx = sta1->x - sta2->x;
			dy = sta1->y - sta2->y;
			dz = sta1->z - sta2->z;
			dist2 = dx * dx + dy * dy + dz * dz;
			cov_matrix[nrow][ncol] = cov_matrix[ncol][nrow] =
				sigmaT2 * exp(-0.5 * dist2 / corr_len2);
		    } else {
			/* different phase types, assumed no correlation  */
			cov_matrix[nrow][ncol] = cov_matrix[ncol][nrow] = 0.0;
		    }

		    /* obs time error */

		    if (ncol == nrow)
			cov_matrix[nrow][ncol] +=
			    arrival[nrow].error * arrival[nrow].error;

		}
	}

	if (message_flag >= 3)
		DisplayDMatrix("Covariance", cov_matrix, num_arrivals,
			num_arrivals);


	/* invert covariance matrix to obtain weight matrix */

	if ((istat = dgaussj(cov_matrix, num_arrivals, null_mtrx, 0)) < 0) {
		puterr("ERROR: inverting covariance matrix.");
		return(-1);
	}

	if (message_flag >= 3)
		DisplayDMatrix("Weight", cov_matrix, num_arrivals,
			num_arrivals);

	/* get row weights & sum of weights */

	weight_sum = 0.0;
	for (nrow = 0; nrow < num_arrivals; nrow++) {
		arrival[nrow].weight = 0.0;
		for (ncol = 0; ncol < num_arrivals; ncol++) {
			arrival[nrow].weight += cov_matrix[nrow][ncol];
			weight_sum += cov_matrix[nrow][ncol];
		}
	}
	for (nrow = 0; nrow < num_arrivals; nrow++) {
		arrival[nrow].weight = (double) num_arrivals
			* arrival[nrow].weight / weight_sum;
		if (arrival[nrow].weight < 0.0) {
			sprintf(MsgStr,
"ERROR: negative observation weight: %s %s %s weight: %lf",
				arrival[nrow].label, arrival[nrow].inst,
				arrival[nrow].comp, arrival[nrow].weight);
			puterr(MsgStr);
			puterr("   Gaussian model error (see LOCGAU) may be too large relative to obs uncertainty (see LOCQUAL2ERR).");
		}
	}
	sprintf(MsgStr, "Weight Matrix sum: %lf", weight_sum);
	putmsg(3, MsgStr);

	/* set global variables */

	gauss_par->WtMtrx = cov_matrix;
	gauss_par->WtMtrxSum = weight_sum;

	return(0);

}



/*** function to calculate weighted mean of observed arrival times */
/*		(TV82, eq. A-38) */

void CalcCenteredTimesObs(int num_arrivals, ArrivalDesc *arrival,
				GaussLocParams* gauss_par, HypoDesc* phypo)
{

	int nrow, ncol, narr;
	long double sum = 0.0L, weighted_mean;
	DMatrix wtmtx;


	if (!FixOriginTimeFlag) {

		/* calculate weighted mean of observed times */

		wtmtx = gauss_par->WtMtrx;

		for (nrow = 0; nrow < num_arrivals; nrow++)
			for (ncol = 0; ncol < num_arrivals; ncol++)
				sum += (long double) wtmtx[nrow][ncol] *
					arrival[ncol].obs_time;
		weighted_mean = sum / (long double) gauss_par->WtMtrxSum;

	} else {

		/* use fixed origin time as reference */

		weighted_mean = phypo->time;
	}


	/* calculate centered observed times */

	putmsg(2, "");
	putmsg(2, "Delayed, Sorted, Centered Observations:");
	for (narr = 0; narr < num_arrivals; narr++) {
		arrival[narr].obs_centered =
			(double) (arrival[narr].obs_time - weighted_mean);
		sprintf(MsgStr,
"  %3d  %-6s %-6s %2.2d:%2.2d:%7.4lf - %7.4lfs -> %8.4lf (%10.4lf)",
				narr, arrival[narr].label, arrival[narr].phase,
				arrival[narr].hour, arrival[narr].min,
				arrival[narr].sec, arrival[narr].delay, arrival[narr].obs_centered,
				((double) arrival[narr].obs_time));
		putmsg(2, MsgStr);
	}

	gauss_par->meanObs = weighted_mean;

}


/*** function to calculate weighted mean of predicted travel times */
/*		(TV82, eq. A-38) */

inline void CalcCenteredTimesPred(int num_arrivals, ArrivalDesc *arrival,
				GaussLocParams* gauss_par)
{

	int nrow, ncol, narr;
	double sum = 0.0, weighted_mean, pred_time_row;
	DMatrix wtmtx;
	double *wtmtxrow;


	if (!FixOriginTimeFlag) {

		wtmtx = gauss_par->WtMtrx;

		for (nrow = 0; nrow < num_arrivals; nrow++) {
			wtmtxrow =  wtmtx[nrow];
			pred_time_row = arrival[nrow].pred_travel_time;
			for (ncol = 0; ncol < num_arrivals; ncol++)
				sum += (double) *(wtmtxrow + ncol)
					* pred_time_row;
		}

		weighted_mean = sum / gauss_par->WtMtrxSum;

	} else {

		/* for fixed origin time use travel time directly */
		weighted_mean = 0.0;

	}



	/* calculate centered predicted times */

	for (narr = 0; narr < num_arrivals; narr++)
		arrival[narr].pred_centered =
			arrival[narr].pred_travel_time - weighted_mean;


	gauss_par->meanPred = (double) weighted_mean;

}



/*** function to calculate probability density */
/*		(MEN92, eq. 14) */

inline double CalcSolutionQuality(int num_arrivals, ArrivalDesc *arrival,
			GaussLocParams* gauss_par, int itype, double* pmisfit)
{

	if (LocMethod == METH_GAU_ANALYTIC)
		return (CalcSolutionQuality_GAU_ANALYTIC(num_arrivals, arrival,
			gauss_par, itype, pmisfit));
	else if (LocMethod == METH_GAU_TEST)
		return (CalcSolutionQuality_GAU_TEST(num_arrivals, arrival,
			gauss_par, itype, pmisfit));
	else
		return(-1.0);

}



/*** function to calculate probability density */
/*	sum of individual L2 residual probablities */

inline double CalcSolutionQuality_GAU_TEST(int num_arrivals, ArrivalDesc *arrival,
			GaussLocParams* gauss_par, int itype, double* pmisfit)
{

	int nrow, ncol, narr;

	double misfit;
	double ln_prob_density, rms_misfit;

	double arr_row_res;
	DMatrix wtmtx;
	double *wtmtxrow;

	double pred_min = VERY_LARGE_DOUBLE, pred_max = -VERY_LARGE_DOUBLE;
	double prob_max = -VERY_LARGE_DOUBLE;
	double tshift, tshift_at_prob_max;
	double tstep, tstart, tstop;
	double misfit_min = VERY_LARGE_DOUBLE;
	double misfit_tmp, prob;


	wtmtx = gauss_par->WtMtrx;


	/* calculate weighted mean of predicted travel times  */
	/*		(TV82, eq. A-38) */

	CalcCenteredTimesPred(num_arrivals, arrival, gauss_par);


	/* calculate residuals (TV82, eqs. 10-12, 10-13; MEN92, eq. 15) */

	for (narr = 0; narr < num_arrivals; narr++) {
		arrival[narr].cent_resid = arrival[narr].obs_centered - arrival[narr].pred_centered;
		if (arrival[narr].pred_centered < pred_min)
			pred_min = arrival[narr].pred_centered;
		if (arrival[narr].pred_centered > pred_max)
			pred_max = arrival[narr].pred_centered;
	}
//printf("pred_min %lf  pred_max %lf\n", pred_min, pred_max);


	// find maximum prob time shift

	tstep = (pred_max - pred_min) / 10.0;
	tstart = pred_min;
	tstop = pred_max;
	while (tstep > (pred_max - pred_min) / 1000000.0) {
		for (tshift = tstart; tshift <= tstop; tshift += tstep) {
			misfit = 0.0;
			prob = 0.0;
			for (nrow = 0; nrow < num_arrivals; nrow++) {
				wtmtxrow = wtmtx[nrow];
				arr_row_res = arrival[nrow].cent_resid + tshift;
				for (ncol = 0; ncol <= nrow; ncol++) {
					if (ncol != nrow) {
						misfit_tmp = (double) *(wtmtxrow + ncol) *
							arr_row_res * (arrival[ncol].cent_resid + tshift);
						//prob += 2.0 * exp(-misfit_tmp);
						//misfit += 2.0 * misfit_tmp;
					} else {
						misfit_tmp = (double) *(wtmtxrow + ncol) *
							arr_row_res * (arrival[ncol].cent_resid + tshift);
						prob += exp(-0.5 * misfit_tmp);
						misfit += misfit_tmp;
					}
				}
			}
			prob /= num_arrivals;
			if (prob > prob_max) {
			//if (misfit < misfit_min) {
				misfit_min = misfit;
				//misfit_min = -log(prob / num_arrivals);
				prob_max = prob;
				tshift_at_prob_max = tshift;
			}
		}
		tstart = tshift_at_prob_max - tstep;
		tstop = tshift_at_prob_max + tstep;
		tstep /= 10.0;
	}
	misfit = misfit_min;
//printf("tshift_at_prob_max %lf  prob_max %lf\n", tshift_at_prob_max, prob_max);


	/* return misfit or ln(prob density) */

	if (itype == GRID_MISFIT) {
		/* convert misfit to rms misfit */
		rms_misfit = sqrt(misfit / num_arrivals);
		*pmisfit = rms_misfit;
		return(rms_misfit);
	} else if (itype == GRID_PROB_DENSITY) {
		//ln_prob_density = -0.5 * misfit;
		ln_prob_density = log(prob_max) * num_arrivals * num_arrivals;
		//rms_misfit = sqrt(misfit / num_arrivals);
		rms_misfit = sqrt(misfit);
		*pmisfit = rms_misfit;
		return(ln_prob_density);
	} else {
		return(-1.0);
	}



}



/*** function to calculate probability density */
/*		(MEN92, eq. 14) */

inline double CalcSolutionQuality_GAU_ANALYTIC(int num_arrivals, ArrivalDesc *arrival,
			GaussLocParams* gauss_par, int itype, double* pmisfit)
{

	int nrow, ncol, narr;

	double misfit = 0.0;
	double ln_prob_density, rms_misfit;

	double arr_row_res;
	DMatrix wtmtx;
	double *wtmtxrow;

	wtmtx = gauss_par->WtMtrx;


	/* calculate weighted mean of predicted travel times  */
	/*		(TV82, eq. A-38) */

	CalcCenteredTimesPred(num_arrivals, arrival, gauss_par);


	/* calculate residuals (TV82, eqs. 10-12, 10-13; MEN92, eq. 15) */

	for (narr = 0; narr < num_arrivals; narr++)
		arrival[narr].cent_resid = arrival[narr].obs_centered -
						arrival[narr].pred_centered;

	for (nrow = 0; nrow < num_arrivals; nrow++) {
		wtmtxrow = wtmtx[nrow];
		arr_row_res = arrival[nrow].cent_resid;
		for (ncol = 0; ncol <= nrow; ncol++)
			if (ncol != nrow)
			    misfit += 2.0 * (double) *(wtmtxrow + ncol) *
				arr_row_res * arrival[ncol].cent_resid;
			else
			    misfit += (double) *(wtmtxrow + ncol) *
				arr_row_res * arrival[ncol].cent_resid;
	}


	/* return misfit or ln(prob density) */

	if (itype == GRID_MISFIT) {
		/* convert misfit to rms misfit */
		rms_misfit = sqrt(misfit / num_arrivals);
		*pmisfit = rms_misfit;
		return(rms_misfit);
	} else if (itype == GRID_PROB_DENSITY) {
		ln_prob_density = -0.5 * misfit;
		rms_misfit = sqrt(misfit / num_arrivals);
		*pmisfit = rms_misfit;
		return(ln_prob_density);
	} else {
		return(-1.0);
	}



}



/*** function to calculate maximum likelihood estimate for origin time */
/*		(MEN92, eq. 19) */

long double CalcMaxLikeOriginTime(int num_arrivals, ArrivalDesc *arrival,
			GaussLocParams* gauss_par)
{

	int nrow, ncol, narr;
	long double time_est = 0.0L;

	DMatrix wtmtx;

	wtmtx = gauss_par->WtMtrx;



/* NOTE: the following (MEN92, eq. 19) is unnecessary

	for (nrow = 0; nrow < num_arrivals; nrow++)
		for (ncol = 0; ncol <= nrow; ncol++) {
			if (ncol != nrow)
			    time_est += 2.0L * (long double) wtmtx[nrow][ncol]
				* (arrival[nrow].obs_time -
				(long double) arrival[nrow].pred_travel_time);
			else
			    time_est += (long double) wtmtx[nrow][ncol] *
				(arrival[nrow].obs_time -
				(long double) arrival[nrow].pred_travel_time);
		}
	return(time_est / (long double) gauss_par->WtMtrxSum);
*/

	return( gauss_par->meanObs - (long double) gauss_par->meanPred);

}



/*** function to update arrival probabilistic residuals */

inline void UpdateProbabilisticResiduals(int num_arrivals,
		ArrivalDesc *arrival, double prob)
{

	int narr;


	/* update probabilistc residuals */

	for (narr = 0; narr < num_arrivals; narr++) {
		arrival[narr].pdf_residual_sum += prob *
			arrival[narr].cent_resid;
 		arrival[narr].pdf_weight_sum += prob;
	}

}





/*** function to calculate confidence intervals and save to file */
/*		(MEN92, eq. 25ff) */

#define N_STEPS_SRCH	101
#define N_STEPS_CONF	11

int CalcConfidenceIntrvl(GridDesc* ptgrid, HypoDesc* phypo,
				char* filename)
{
	FILE *fpio;
	char fname[FILENAME_MAX];

	int ix, iy, iz;
	int iconf, isrch;
	double srch_level, srch_incr, conf_level, conf_incr, prob_den;
	double srch_sum[N_STEPS_SRCH];
	double contour[N_STEPS_CONF];
	double sum_volume;


	/* write message */
	putmsg(2, "");
	putmsg(2, "Calculating confidence intervals over grid...");


	for (isrch = 0; isrch < N_STEPS_SRCH; isrch++)
		srch_sum[isrch] = 0.0;


	/* accumulate approx integral of probability density in search bins */
	/*  and normalize sum of bin values * dx*dy*dz */

	sum_volume = ptgrid->sum * ptgrid->dx * ptgrid->dy * ptgrid->dz;
	phypo->probmax /= sum_volume;
	srch_incr = phypo->probmax / (N_STEPS_SRCH - 1);
	for (ix = 0; ix <  ptgrid->numx; ix++) {
	    for (iy = 0; iy <  ptgrid->numy; iy++) {
		for (iz = 0; iz <  ptgrid->numz; iz++) {

			ptgrid->array[ix][iy][iz] = (float) (
				exp((double) ptgrid->array[ix][iy][iz])
					/ sum_volume );
			prob_den = ptgrid->array[ix][iy][iz];
			srch_level = 0.0;
			for (isrch = 0; isrch < N_STEPS_SRCH; isrch++) {
				if (prob_den >= srch_level)
					srch_sum[isrch] += prob_den;
				srch_level += srch_incr;
			}

		}
	    }
	}
	ptgrid->sum = 1.0;

	/* normalize by 100% confidence level sum */

	for (isrch = 1; isrch < N_STEPS_SRCH; isrch++)
		srch_sum[isrch] /= srch_sum[0];
	srch_sum[0] = 1.0;


	/* open confidence interval file */

	sprintf(fname, "%s.loc.conf", filename);
	if ((fpio = fopen(fname, "w")) == NULL) {
		puterr("ERROR: opening confidence interval output file.");
		return(-1);
	} else {
		NumFilesOpen++;
	}

	/* find confidence levels and write to file */

	conf_incr = 1.0 / (N_STEPS_CONF - 1);
	conf_level = 1.0;
	iconf = N_STEPS_CONF - 1;
	for (isrch = 0; isrch < N_STEPS_SRCH; isrch++) {
		if (srch_sum[isrch] <= conf_level) {
			contour[iconf] = (double) isrch * srch_incr;
			fprintf(fpio, "%lf C %.2lf\n",
				contour[iconf], conf_level);
			if (--iconf < 0)
				break;
			conf_level -= conf_incr;
		}
	}


	fclose(fpio);
	NumFilesOpen--;

	return(0);

}



/*** function to generate sample (scatter) of location PDF */

int GenEventScatterGrid(GridDesc* ptgrid, HypoDesc* phypo, ScatterParams* pscat,
	char* filename)
{
	FILE *fpio;
	char fname[FILENAME_MAX];

	int ix, iy, iz;
	double origx, origy, origz;
	double dx, dy, dz, dvol;
	double xval, yval, zval;
	float fdata[4];
	int tot_npoints = 0;
	double xnpt, xnpoints;
	double probmax, prob_den, prob;



	/* return if no scatter samples requested */
	if (pscat->npts < 1)
		return(0);

	/* write message */
	putmsg(2, "");
	putmsg(2, "Generating event scatter file...");

	/* open scatter file */

	sprintf(fname, "%s.loc.scat", filename);
	if ((fpio = fopen(fname, "w")) == NULL) {
		puterr("ERROR: opening scatter output file.");
		return(-1);
	} else {
		NumFilesOpen++;
	}
	/* skip header record (used later to store number of samples taken) */
	fseek(fpio, 4 * sizeof(float), SEEK_SET);


	/* generate N=Scatter->npts events with prob P=prob_den/probmax at  */
	/*	uniformly-randomly chosen grid locations */

	origx = ptgrid->origx;
	origy = ptgrid->origy;
	origz = ptgrid->origz;
	dx = ptgrid->dx;
	dy = ptgrid->dy;
	dz = ptgrid->dz;
	dvol = dx * dy * dz;
	probmax = phypo->probmax;

	for (ix = 0; ix <  ptgrid->numx; ix++) {
	    for (iy = 0; iy <  ptgrid->numy; iy++) {
		for (iz = 0; iz <  ptgrid->numz; iz++) {

			prob_den =  ptgrid->array[ix][iy][iz];

			xnpoints = (double) pscat->npts * dvol * prob_den;

			xval = origx + (double) ix * dx;
			yval = origy + (double) iy * dy;
			zval = origz + (double) iz * dz;

			while (xnpoints > 0.0) {

			    if (xnpoints > 1.0 ||
					xnpoints - (double) ((int) xnpoints)
						> get_rand_double(0.0, 1.0)) {
				fdata[0] =
				    xval + get_rand_double(-dx / 2.0, dx / 2.0);
				fdata[1] =
				    yval + get_rand_double(-dy / 2.0, dy / 2.0);
				fdata[2] =
				    zval + get_rand_double(-dz / 2.0, dz / 2.0);
				fdata[3] = prob_den;
				fwrite(fdata, sizeof(float), 4, fpio);

				tot_npoints++;
			    }

			    xnpoints -= 1.0;

			}

		}
	    }
	}


	/* write header informaion */
	fseek(fpio, 0, SEEK_SET);
	fwrite(&tot_npoints, sizeof(int), 1, fpio);
	fdata[0] = (float) probmax;
	fwrite(fdata, sizeof(float), 1, fpio);

	fclose(fpio);
	NumFilesOpen--;

	/* write message */
	sprintf(MsgStr, "  %d points generated.", tot_npoints);
	putmsg(2, MsgStr);
	sprintf(MsgStr, "  (%d points requested, dvol= %lf, probmax=%lf)",
			pscat->npts, dvol, probmax);
	putmsg(2, MsgStr);

	return(0);

}



/*** function to read input file */

int ReadNLLoc_Input(FILE* fp_input)
{
	int istat, iscan;
	char param[MAXLINE] = "\0", *pchr;
	char line[MAXLINE_LONG], *fgets_return;

	int flag_control = 0, flag_outfile = 0, flag_grid = 0, flag_search = 0,
		flag_method = 0, flag_comment = 0, flag_signature = 0,
		flag_hyptype = 0, flag_gauss = 0, flag_trans = 0, flag_comp = 0,
		flag_phstat = 0, flag_phase_id = 0, flag_qual2err = 0,
		flag_mag = 0, flag_alias = 0, flag_exclude = 0, flag_time_delay = 0,
		flag_otime = 0, flag_angles = 0;
	int flag_include = 1;



	/* read each input line */

	/* read each input line */

	while ((fgets_return = fgets(line, 4*MAXLINE, fp_input)) != NULL
			|| fp_include != NULL) {


		/* check for end of include file */

		if (fgets_return == NULL && fp_include != NULL) {
			SwapBackIncludeFP(&fp_input);
			continue;
		}


		istat = -1;

		/*read parameter line */

		if ((iscan = sscanf(line, "%s", param)) < 0 )
			continue;

		/* skip comment line or white space */

		if (strncmp(param, "#", 1) == 0 || isspace(param[0]))
			istat = 0;


		/* read include file params and set input to include file */

		if (strcmp(param, "INCLUDE") == 0)
			if ((istat = GetIncludeFile(strchr(line, ' '),
							&fp_input)) < 0) {
				puterr("ERROR: processing include file.");
				flag_include = 0;
			}


		/* read control params */

		if (strcmp(param, "CONTROL") == 0)
			if ((istat = get_control(strchr(line, ' '))) < 0)
				puterr("ERROR: readingng control params.");
			else
				flag_control = 1;


		/* read grid params */

		if (strcmp(param, "LOCGRID") == 0)
    			if ((istat = GetNLLoc_Grid(strchr(line, ' '))) < 0)
				puterr("ERROR: reading grid parameters.");
			else
				flag_grid = 1;


		/* read file names */

		if (strcmp(param, "LOCFILES") == 0)
			if ((istat = GetNLLoc_Files(strchr(line, ' '))) < 0)
			  puterr("ERROR: reading NLLoc output file name.");
			else
				flag_outfile = 1;


		/* read output file types names */

		if (strcmp(param, "LOCHYPOUT") == 0)
			if ((istat = GetNLLoc_HypOutTypes(strchr(line, ' '))) < 0)
				puterr("ERROR: reading NLLoc hyp output file types.");
			else
				flag_hyptype = 1;


		/* read search type */

		if (strcmp(param, "LOCSEARCH") == 0)
			if ((istat = GetNLLoc_SearchType(strchr(line, ' '))) < 0)
				puterr("ERROR: reading NLLoc search type.");
			else
				flag_search = 1;


		/* read method */

		if (strcmp(param, "LOCMETH") == 0)
			if ((istat = GetNLLoc_Method(strchr(line, ' '))) < 0)
				puterr("ERROR: reading NLLoc method.");
			else
				flag_method = 1;


		/* read fixed origin time parameters */

		if (strcmp(param, "LOCFIXOTIME") == 0)
			if ((istat = GetNLLoc_FixOriginTime(
						strchr(line, ' '))) < 0)
				puterr("ERROR: reading NLLoc fixed origin time params.");
			else
				flag_otime = 1;


		/* read phase identifier values */

		if (strcmp(param, "LOCPHASEID") == 0)
			if ((istat = GetPhaseID(strchr(line, ' '))) < 0)
				puterr("ERROR: reading phase identifier values.");
			else
				flag_phase_id = 1;


		/* read quality2error values */

		if (strcmp(param, "LOCQUAL2ERR") == 0)
			if ((istat = GetQuality2Err(strchr(line, ' '))) < 0)
				puterr("ERROR: reading quality2error values.");
			else
				flag_qual2err = 1;


		/* read comment */

		if (strcmp(param, "LOCCOM") == 0) {
			strcpy(Hypocenter.comment, strchr(line, ' ') + 1);
			*(strchr(Hypocenter.comment, '\n')) = '\0';
			sprintf(MsgStr, "LOCCOMMENT:  %s\n",
				Hypocenter.comment);
			putmsg(2, MsgStr);
			flag_comment = 1;
		}


		/* read signature */

		if (strcmp(param, "LOCSIG") == 0) {
			strcpy(LocSignature, strchr(line, ' ') + 1);
			*(strchr(LocSignature, '\n')) = '\0';
			sprintf(MsgStr, "LOCSIGNATURE:  %s\n",
				LocSignature);
			putmsg(2, MsgStr);
			flag_signature = 1;
		}


		/* read gauss params */

		if (strcmp(param, "LOCGAU") == 0)
			if ((istat = GetNLLoc_Gaussian(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Gaussian parameters.");
			else
				flag_gauss = 1;



		/* read phase statistics params */

		if (strcmp(param, "LOCPHSTAT") == 0)
			if ((istat = GetNLLoc_PhaseStats(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Phase Statistics parameters.");
			else
				flag_phstat = 1;


		/* read take-off angles params */

		if (strcmp(param, "LOCANGLES") == 0)
			if ((istat = GetNLLoc_Angles(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Take-off Angles parameters.");
			else
				flag_angles = 1;


		/* read magnitude calculation params */

		if (strcmp(param, "LOCMAG") == 0)
			if ((istat = GetNLLoc_Magnitude(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Magnitude Calculation parameters.");
			else
				flag_mag = 1;


		/* read component params */

		if (strcmp(param, "LOCCMP") == 0)
			if ((istat = GetCompDesc(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Component description parameters.");
			else
				flag_comp = 1;


		/* read alias params */

		if (strcmp(param, "LOCALIAS") == 0)
			if ((istat = GetLocAlias(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Alias parameters.");
			else
				flag_alias = 1;


		/* read exclude params */

		if (strcmp(param, "LOCEXCLUDE") == 0)
			if ((istat = GetLocExclude(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Exclude parameters.");
			else
				flag_exclude = 1;


		/* read station time delay params */

		if (strcmp(param, "LOCDELAY") == 0)
			if ((istat = GetTimeDelays(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Time Delay parameters.");
			else
				flag_time_delay = 1;


		/*read transform params */

		if (strcmp(param, "TRANS") == 0)
    			if ((istat = get_transform(0, strchr(line, ' '))) < 0)
			    puterr("ERROR: reading transformation parameters.");
			else
				flag_trans = 1;


		/* unrecognized input */

		if (istat < 0) {
			if ((pchr = strchr(line, '\n')) != NULL)
				*pchr = '\0';
			sprintf(MsgStr, "Skipping input: %s", line);
			putmsg(4, MsgStr);
		}

	}


	/* check for missing required input */

	if (!flag_control)
		puterr("ERROR: reading control (CONTROL) params.");
	if (!flag_outfile)
		puterr("ERROR: reading i/o file (LOCFILES) params.");
	if (!flag_trans)
		puterr("ERROR: reading transformation (TRANS) params.");
	if (!flag_grid)
		puterr("ERROR: reading grid (LOCGRID) params.");
	if (!flag_search)
		puterr("ERROR: reading search type (LOCSEARCH) params.");
	if (!flag_method)
		puterr("ERROR: reading method (LOCMETH) params.");
	if (!flag_gauss)
		puterr("ERROR: reading Gaussian (LOCGAU) params.");
	if (!flag_qual2err)
		puterr("ERROR: reading Quality2Error (LOCQUAL2ERR) params.");


	/* check for missing optional input */

	if (!flag_comment) {
		sprintf(MsgStr, "INFO: no comment (LOCCOM) params read.");
		putmsg(1, MsgStr);
		Hypocenter.comment[0] = '\0';
	}
	if (!flag_signature) {
		sprintf(MsgStr, "INFO: no signature (LOCSIG) params read.");
		putmsg(1, MsgStr);
		LocSignature[0] = '\0';
	}
	if (!flag_hyptype) {
		sprintf(MsgStr,
"INFO: no hypocenter output file type (LOCHYPOUT) params read.");
		putmsg(1, MsgStr);
		sprintf(MsgStr,
"INFO: DEFAULT: \"LOCHYPOUT SAVE_NLLOC_ALL SAVE_HYPOINV_SUM\"");
		putmsg(1, MsgStr);
		iSaveNLLocEvent = iSaveNLLocSum = 1;
		iSaveHypo71Sum = iSaveHypoEllSum = 0;
		iSaveHypo71Event = iSaveHypoEllEvent = 0;
		iSaveHypoInvSum = 1;
		iSaveAlberto4Sum = 1;
	}
	if (!flag_phase_id) {
		sprintf(MsgStr,
"INFO: no phase identifier (LOCPHASEID) values read.");
		putmsg(1, MsgStr);
	}
	if (!flag_mag) {
		sprintf(MsgStr, "INFO: no Magnitude Calculation (LOCMAG) params read.");
		putmsg(1, MsgStr);
	}
	if (!flag_phstat) {
		sprintf(MsgStr,
"INFO: no PhaseStatistics (LOCPHSTAT) params read.");
		putmsg(1, MsgStr);
		RMS_Max = VERY_LARGE_DOUBLE;
		NRdgs_Min = -1;
		Gap_Max = VERY_LARGE_DOUBLE;
		P_ResidualMax = VERY_LARGE_DOUBLE;
		S_ResidualMax = VERY_LARGE_DOUBLE;
	}
	if (!flag_comp) {
		sprintf(MsgStr, "INFO: no Component Descirption (LOCCMP) params read.");
		putmsg(1, MsgStr);
	}
	if (!flag_alias) {
		sprintf(MsgStr, "INFO: no Alias (LOCALIAS) params read.");
		putmsg(1, MsgStr);
	}
	if (!flag_exclude) {
		sprintf(MsgStr, "INFO: no Exclude (LOCEXCLUDE) params read.");
		putmsg(1, MsgStr);
	}
	if (!flag_time_delay) {
		sprintf(MsgStr, "INFO: no Time Delay (LOCDELAY) params read.");
		putmsg(1, MsgStr);
	}
	if (!flag_otime) {
		sprintf(MsgStr, "INFO: no Fixed Origin Time (LOCFIXOTIME) params read.");
		putmsg(1, MsgStr);
	}
	if (!flag_angles) {
		sprintf(MsgStr, "INFO: no Take-off Angles (LOCANGLES) params read,");
		putmsg(1, MsgStr);
		sprintf(MsgStr, "      default is angleMode=ANGLES_NO, qualtiyMin=5 .");
		putmsg(1, MsgStr);
		angleMode = ANGLE_MODE_NO;
		iAngleQualityMin = 5;
	}


	return (flag_include * flag_control * flag_outfile * flag_grid *
		flag_search *
		flag_method * flag_gauss * flag_qual2err * flag_trans - 1);
}



/*** function to read output file name ***/

int GetNLLoc_Files(char* line1)
{
	int nObsFile;
	char fnobs[FILENAME_MAX];

	sscanf(line1, "%s %s %s %s", fnobs, ftype_obs, fn_loc_grids,
			fn_path_output);
//printf("FILENAME_MAX %d  fn_loc_grids <%s>\n", FILENAME_MAX, fn_loc_grids);

	/* check for wildcards in observation file name */
	NumObsFiles = ExpandWildCards(fnobs, fn_loc_obs, MAX_NUM_OBS_FILES);

	sprintf(MsgStr,
		"LOCFILES:  ObsType: %s  InGrids: %s.*  OutPut: %s.*",
		 ftype_obs, fn_loc_grids, fn_path_output);
	putmsg(2, MsgStr);
	for (nObsFile = 0; nObsFile < NumObsFiles; nObsFile++) {
		sprintf(MsgStr,
			"   Obs File: %3d  %s", nObsFile,
			fn_loc_obs[nObsFile]);
			putmsg(2, MsgStr);
		}
	if (NumObsFiles == MAX_NUM_OBS_FILES)
		putmsg(1,
"LOCFILES: WARNING: maximum number of files/events reached");

	return(0);
}




/*** function to read search type ***/

int GetNLLoc_SearchType(char* line1)
{
	int istat, ierr;

	char search_type[MAXLINE];


	istat = sscanf(line1, "%s", search_type);

	if (istat != 1)
		return(-1);

	if (strcmp(search_type, "GRID") == 0) {

		SearchType = SEARCH_GRID;
		istat = sscanf(line1, "%s %d",
			search_type, &(Scatter.npts));
		if (istat != 2)
			return(-1);

		sprintf(MsgStr, "LOCSEARCH:  Type: %s NumScatter %d",
			search_type, Scatter.npts);
		putmsg(2, MsgStr);

	} else if (strcmp(search_type, "MET") == 0) {

		SearchType = SEARCH_MET;
		istat = sscanf(line1, "%s %d %d %d %d %d %lf %lf %lf %lf",
			search_type, &MetNumSamples, &MetLearn, &MetEquil,
			&MetStartSave, &MetSkip,
			&MetStepInit, &MetStepMin, &MetStepFact, &MetProbMin);
		ierr = 0;

		sprintf(MsgStr,
"LOCSEARCH:  Type: %s  numSamples %d  numLearn %d  numEquilibrate %d  startSave %d  numSkip %d  stepInit %lf  stepMin %lf  stepFact %lf  probMin %lf",
			search_type, MetNumSamples, MetLearn, MetEquil,
			MetStartSave, MetSkip,
			MetStepInit, MetStepMin, MetStepFact, MetProbMin);
		putmsg(2, MsgStr);

		if (checkRangeInt("LOCSEARCH", "numSamples", MetNumSamples, 1, 0, 0, 0) != 0)
			ierr = -1;
		if (checkRangeInt("LOCSEARCH", "numLearn", MetLearn, 1, 0, 0, 0) != 0)
			ierr = -1;
		if (checkRangeInt("LOCSEARCH", "numEquilibrate", MetEquil, 1, 0, 0, 0) != 0)
			ierr = -1;
		if (checkRangeInt("LOCSEARCH", "startSave", MetStartSave, 1, 0, 0, 0) != 0)
			ierr = -1;
		if (checkRangeInt("LOCSEARCH", "numSkip", MetSkip, 1, 1, 0, 0) != 0)
			ierr = -1;
		if (checkRangeDouble("LOCSEARCH", "stepMin", MetStepMin, 1, 0.0, 0, 0.0) != 0)
			ierr = -1;
		if (ierr < 0)
			return(-1);
		if (istat != 10)
			return(-1);

		//?? AJL 17JAN2000 MetUse = MetNumSamples - MetEquil;
		MetUse = MetNumSamples - MetStartSave;

		/* check for "normal" StartSave value*/
		if (MetStartSave < MetLearn + MetEquil) {
			sprintf(MsgStr,
"LOCSEARCH:  WARNING: Metropolis StartSave < NumLearn + NumEquilibrate.");
			putmsg(1, MsgStr);
		}

	} else if (strcmp(search_type, "OCT") == 0) {

		SearchType = SEARCH_OCTTREE;
		istat = sscanf(line1, "%s %d %d %d %lf %d %d",
			search_type, &octtreeParams.init_num_cells_x,
			&octtreeParams.init_num_cells_y, &octtreeParams.init_num_cells_z,
			&octtreeParams.min_node_size, &octtreeParams.max_num_nodes,
			&octtreeParams.num_scatter);
		ierr = 0;

		sprintf(MsgStr,
"LOCSEARCH:  Type: %s  init_num_cells_x %d  init_num_cells_y %d  init_num_cells_z %d  min_node_size %lf  max_num_nodes %d  num_scatter %lf",
			search_type, octtreeParams.init_num_cells_x, octtreeParams.init_num_cells_y,
			octtreeParams.init_num_cells_z,
			octtreeParams.min_node_size, octtreeParams.max_num_nodes,
			octtreeParams.num_scatter);
		putmsg(2, MsgStr);

		if (checkRangeInt("LOCSEARCH", "init_num_cells_x",
				octtreeParams.init_num_cells_x, 1, 0, 0, 0) != 0)
			ierr = -1;
		if (checkRangeInt("LOCSEARCH", "init_num_cells_y",
				octtreeParams.init_num_cells_y, 1, 0, 0, 0) != 0)
			ierr = -1;
		if (checkRangeInt("LOCSEARCH", "init_num_cells_z",
				octtreeParams.init_num_cells_z, 1, 0, 0, 0) != 0)
			ierr = -1;
		if (checkRangeDouble("LOCSEARCH", "min_node_size",
				octtreeParams.min_node_size, 1, 0.0, 0, 0.0) != 0)
			ierr = -1;
		if (checkRangeInt("LOCSEARCH", "max_num_nodes",
				octtreeParams.max_num_nodes, 1, 0, 0, 0) != 0)
			ierr = -1;
		if (checkRangeInt("LOCSEARCH", "num_scatter",
				octtreeParams.num_scatter, 1, 0, 0, 0) != 0)
			ierr = -1;
		if (ierr < 0)
			return(-1);
		if (istat != 7)
			return(-1);

	}


	return(0);
}




/*** function to read requested hypocenter output file types ***/

int GetNLLoc_HypOutTypes(char* line1)
{
	int istat;

	char *pchr, hyp_type[MAXLINE];


	sprintf(MsgStr, "LOCHYPOUT:  ");

	iSaveNLLocEvent = iSaveNLLocSum = iSaveHypo71Event = iSaveHypo71Sum
		= iSaveHypoEllEvent = iSaveHypoEllSum = iSaveHypoInvSum
		= iSaveAlberto4Sum = 0;

	pchr = line1;
	do {

		/* check for blank line */
		while (*pchr == ' ')
			pchr++;
		if (isspace(*pchr))
			break;

		if ((istat = sscanf(pchr, "%s", hyp_type)) != 1)
			return(-1);

		if (strcmp(hyp_type, "SAVE_NLLOC_ALL") == 0)
			iSaveNLLocEvent = iSaveNLLocSum = 1;
		else if (strcmp(hyp_type, "SAVE_NLLOC_SUM") == 0)
			iSaveNLLocSum = 1;
		else if (strcmp(hyp_type, "SAVE_HYPO71_ALL") == 0)
			iSaveHypo71Event = iSaveHypo71Sum = 1;
		else if (strcmp(hyp_type, "SAVE_HYPO71_SUM") == 0)
			iSaveHypo71Sum = 1;
		else if (strcmp(hyp_type, "SAVE_HYPOELL_ALL") == 0)
			iSaveHypoEllEvent = iSaveHypoEllSum = 1;
		else if (strcmp(hyp_type, "SAVE_HYPOELL_SUM") == 0)
			iSaveHypoEllSum = 1;
		else if (strcmp(hyp_type, "SAVE_HYPOINV_SUM") == 0)
			iSaveHypoInvSum = 1;
		else if (strcmp(hyp_type, "SAVE_ALBERTO_3D_4") == 0)
			iSaveAlberto4Sum = 1;
		else
			return(-1);

		strcat(MsgStr, hyp_type);
		strcat(MsgStr, " ");

	} while ((pchr = strchr(pchr + 1, ' ')) != NULL);

	putmsg(2, MsgStr);

	return(0);
}




/*** function to read method ***/

int GetNLLoc_Method(char* line1)
{
	int istat, ierr;

	char loc_method[MAXLINE];


	istat = sscanf(line1, "%s %lf %d %d %d %lf %d", loc_method,
		&DistStaGridMax, &MinNumArrLoc, &MaxNumArrLoc, &MinNumSArrLoc,
		&VpVsRatio, &MaxNum3DGridMemory);

	sprintf(MsgStr,
"LOCMETHOD:  method: %s  maxDistStaGrid: %lf  minNumberPhases: %d  maxNumberPhases: %d  minNumberSphases: %d  VpVsRatio: %lf,   max3DGridMemory: %d",
		loc_method, DistStaGridMax, MinNumArrLoc, MaxNumArrLoc,
		MinNumSArrLoc, VpVsRatio, MaxNum3DGridMemory);
	putmsg(2, MsgStr);

	ierr = 0;

	if (ierr < 0 || istat != 7)
		return(-1);

	if (strcmp(loc_method, "GAU_ANALYTIC") == 0)
		LocMethod = METH_GAU_ANALYTIC;
	else if (strcmp(loc_method, "GAU_TEST") == 0)
		LocMethod = METH_GAU_TEST;
	else {
		LocMethod = METH_UNDEF;
		puterr2("ERROR: unrecognized location method:", loc_method);
	}

	if (MaxNumArrLoc < 1)
		MaxNumArrLoc = MAX_NUM_ARRIVALS;

	return(0);
}



/*** function to read fixed origin time parameters ***/

int GetNLLoc_FixOriginTime(char* line1)
{
	int istat;


	istat = sscanf(line1, "%d %d %d %d %d %lf",
		&Hypocenter.year, &Hypocenter.month, &Hypocenter.day,
		&Hypocenter.hour, &Hypocenter.min, &Hypocenter.sec);

	sprintf(MsgStr,
	   "LOCFIXOTIME:  %4.4d%2.2d%2.2d %2.2d%2.2d %5.2lf",
		Hypocenter.year, Hypocenter.month, Hypocenter.day,
		Hypocenter.hour, Hypocenter.min, Hypocenter.sec);
	putmsg(2, MsgStr);

	if (istat != 6)
		return(-1);

	FixOriginTimeFlag = 1;

	return(0);
}




/** function to read grid params */

int GetNLLoc_Grid(char* input_line)
{
	int istat;
	char str_save[20];

	istat = sscanf(input_line, "%d %d %d %lf %lf %lf %lf %lf %lf %s %s",
		&(grid_in.numx), &(grid_in.numy), &(grid_in.numz),
		&(grid_in.origx), &(grid_in.origy), &(grid_in.origz),
		&(grid_in.dx), &(grid_in.dy), &(grid_in.dz), grid_in.chr_type,
		str_save);

	convert_grid_type(&grid_in);
	if (message_flag >= 2)
		display_grid_param(&grid_in);
	sprintf(MsgStr, "LOCGRID: Save: %s", str_save);
	putmsg(2, MsgStr);

	if (istat != 11)
		return(-1);

	if (NumLocGrids < MAX_NUM_LOCATION_GRIDS) {
		LocGrid[NumLocGrids] = grid_in;
		LocGrid[NumLocGrids].autox = 0;
		LocGrid[NumLocGrids].autoy = 0;
		LocGrid[NumLocGrids].autoz = 0;
		if (LocGrid[NumLocGrids].origx < -LARGE_DOUBLE)
			LocGrid[NumLocGrids].autox = 1;
		if (LocGrid[NumLocGrids].origy < -LARGE_DOUBLE)
			LocGrid[NumLocGrids].autoy = 1;
		if (LocGrid[NumLocGrids].origz < -LARGE_DOUBLE)
			LocGrid[NumLocGrids].autoz = 1;
		if (strcmp(str_save, "SAVE") == 0)
			LocGridSave[NumLocGrids] = 1;
		else
			LocGridSave[NumLocGrids] = 0;
		NumLocGrids++;
	} else
		puterr("WARNING: maximum number of location grids exceeded.");

	return(0);
}



/*** function to read phase identification values ***/

int GetPhaseID(char* line1)
{
	int istat, ilen;
	char *substr, *cpos;


	if (NumPhaseID >= MAX_NUM_PHASE_ID) {
		puterr(
"LOCPHASEID: WARNING: maximum number of PhaseIDs reached, ignoring phase ID.");
		return(-1);
	}

	if ((istat = sscanf(line1, "%s", PhaseID[NumPhaseID].phase)) != 1)
		return(-1);

	substr = strstr(line1, PhaseID[NumPhaseID].phase);

	/* save phase id values with spaces at each end */
	if ((cpos = strchr(substr, '\n')) != NULL)
		*cpos = '\0';
	sprintf(PhaseID[NumPhaseID].id_string, " %s ", substr + 1);

	if ((ilen = strlen(PhaseID[NumPhaseID].id_string)) == 0)
		return(-1);

	sprintf(MsgStr, "LOCPHASEID:");
	putmsg(2, MsgStr);
	sprintf(MsgStr, "  Phase: %s  PhaseID: <%s>",
			PhaseID[NumPhaseID].phase,
			PhaseID[NumPhaseID].id_string);
	putmsg(2, MsgStr);

	NumPhaseID++;

	return(0);
}






/*** function to read gaussian params ***/

int GetNLLoc_Gaussian(char* line1)
{
	int istat, ierr;


	istat = sscanf(line1, "%lf %lf", &(Gauss.SigmaT), &(Gauss.CorrLen));

	sprintf(MsgStr, "LOCGAUSS:  SigmaT: %lf  CorrLen: %lf",
		Gauss.SigmaT, Gauss.CorrLen);
	putmsg(2, MsgStr);

	ierr = 0;
	if (checkRangeDouble("LOCGAU", "SigmaT",
			Gauss.SigmaT, 1, 0.0, 0, 0.0) != 0)
		ierr = -1;
	if (checkRangeDouble("LOCGAU", "CorrLen",
			Gauss.CorrLen, 1, 0.0, 0, 0.0) != 0)
		ierr = -1;

	if (ierr < 0 || istat != 2)
		return(-1);

	return(0);
}





/*** function to read magnitude calculation type ***/

int GetNLLoc_Magnitude(char* line1)
{
	int istat, ierr;

	char mag_type[MAXLINE];

	if (NumMagnitudeMethods >= MAX_NUM_MAG_METHODS) {
		puterr2("ERROR: maximum number of LOCMAG statements read: ignoring: ", line1);
		return(-1);
	}

	istat = sscanf(line1, "%s", mag_type);

	if (istat != 1)
		return(-1);

	if (strcmp(mag_type, "ML_HB") == 0) {
	
		// default values
		Magnitude[NumMagnitudeMethods].hb_Ro = 100.0;
		Magnitude[NumMagnitudeMethods].hb_Mo = 3.0;

		Magnitude[NumMagnitudeMethods].type = MAG_ML_HB;
		istat = sscanf(line1, "%s %lf %lf %lf %lf %lf",
			mag_type, &(Magnitude[NumMagnitudeMethods].amp_fact_ml_hb),
			&(Magnitude[NumMagnitudeMethods].hb_n), &(Magnitude[NumMagnitudeMethods].hb_K),
			&(Magnitude[NumMagnitudeMethods].hb_Ro), &(Magnitude[NumMagnitudeMethods].hb_Mo));
		sprintf(MsgStr, "LOCMAGNITUDE:  Type: %s  f %e  n %f  K %f  Ro %f  Mo %f",
			mag_type, Magnitude[NumMagnitudeMethods].amp_fact_ml_hb, Magnitude[NumMagnitudeMethods].hb_n,
			Magnitude[NumMagnitudeMethods].hb_K,
			Magnitude[NumMagnitudeMethods].hb_Ro, Magnitude[NumMagnitudeMethods].hb_Mo);
		putmsg(0, MsgStr);
		//putmsg(2, MsgStr);

		ierr = 0;
		if (checkRangeDouble("LOCMAG", "f", Magnitude[NumMagnitudeMethods].amp_fact_ml_hb, 1, 0.0, 0, 0.0) != 0)
			ierr = -1;

		if (istat < 4)
			return(-1);

	} else if (strcmp(mag_type, "MD_FMAG") == 0) {

		Magnitude[NumMagnitudeMethods].type = MAG_MD_FMAG;
		istat = sscanf(line1, "%s %lf %lf %lf %lf %lf",
			mag_type, &(Magnitude[NumMagnitudeMethods].fmag_c1), &(Magnitude[NumMagnitudeMethods].fmag_c2),
			&(Magnitude[NumMagnitudeMethods].fmag_c3), &(Magnitude[NumMagnitudeMethods].fmag_c4), &(Magnitude[NumMagnitudeMethods].fmag_c5));
		sprintf(MsgStr, "LOCMAGNITUDE:  Type: %s  C1 %lf  C2 %lf  C3 %lf  C4 %lf  C5 %lf",
			mag_type, Magnitude[NumMagnitudeMethods].fmag_c1, Magnitude[NumMagnitudeMethods].fmag_c2, Magnitude[NumMagnitudeMethods].fmag_c3,
			Magnitude[NumMagnitudeMethods].fmag_c4, Magnitude[NumMagnitudeMethods].fmag_c5);
		putmsg(2, MsgStr);

		if (istat != 6)
			return(-1);

	} else {
		Magnitude[NumMagnitudeMethods].type = MAG_UNDEF;
		puterr2("ERROR: unrecognized magnitude calculation type:", mag_type);
	}

	NumMagnitudeMethods++;


	return(0);
}




/*** function to read phase statistics params ***/

int GetNLLoc_PhaseStats(char* line1)
{
	int istat;


	istat = sscanf(line1, "%lf %d %lf %lf %lf",
		&RMS_Max, &NRdgs_Min, &Gap_Max, &P_ResidualMax, &S_ResidualMax);

	sprintf(MsgStr,
"LOCPHSTAT:  RMS_Max: %lf  NRdgs_Min: %d  Gap_Max: %lf  P_ResidualMax: %lf S_ResidualMax: %lf",
			RMS_Max, NRdgs_Min, Gap_Max,
			P_ResidualMax, S_ResidualMax);
	putmsg(2, MsgStr);

	if (istat != 5)
		return(-1);

	return(0);
}


/*** function to read angles mode params ***/

int GetNLLoc_Angles(char* line1)
{
	char strAngleMode[MAXLINE];


	sscanf(line1, "%s %d", strAngleMode, &iAngleQualityMin);

	sprintf(MsgStr, "LOCANGLES:  %s  %d", strAngleMode, iAngleQualityMin);
	putmsg(3, MsgStr);

	if (strcmp(strAngleMode, "ANGLES_YES") == 0)
		angleMode = ANGLE_MODE_YES;
	else if (strcmp(strAngleMode, "ANGLES_NO") == 0)
		angleMode = ANGLE_MODE_NO;
	else {
		angleMode = ANGLE_MODE_UNDEF;
		puterr("ERROR: unrecognized angle mode");
		return(-1);
	}

	return(0);

}



/*** function to read component description ***/

int GetCompDesc(char* line1)
{
	int istat, ierr;

	if (NumCompDesc >= MAX_NUM_COMP_DESC) {
		sprintf(MsgStr, "%s", line1);
		putmsg(1, MsgStr);
		sprintf(MsgStr,
"WARNING: maximum number of component descriptions reached, ignoring description.");
		putmsg(1, MsgStr);
		return(-1);
	}

	Component[NumCompDesc].sta_corr_md_fmag = 1.0;	// fmag default

	istat = sscanf(line1, "%s %s %s %lf %lf %lf",
		Component[NumCompDesc].label, Component[NumCompDesc].inst,
		Component[NumCompDesc].comp, &(Component[NumCompDesc].amp_fact_ml_hb),
		&(Component[NumCompDesc].sta_corr_ml_hb),
		&(Component[NumCompDesc].sta_corr_md_fmag));

	sprintf(MsgStr,
"LOCCMP:  Label: %s  Inst: %s  Comp: %s  Afact: %lf  StaCorr_ML_HB: %lf  StaCorr_MD_FMAG: %lf",
		Component[NumCompDesc].label, Component[NumCompDesc].inst,
		Component[NumCompDesc].comp, Component[NumCompDesc].amp_fact_ml_hb,
		Component[NumCompDesc].sta_corr_ml_hb,
		Component[NumCompDesc].sta_corr_md_fmag);
	putmsg(2, MsgStr);

	ierr = 0;
		if (checkRangeDouble("LOCCMP", "amp_fact_ml_hb",
				Component[NumCompDesc].amp_fact_ml_hb, 1, 0.0, 0, 0.0) != 0)
			ierr = -1;

	if (ierr < 0 || istat < 5)
		return(-1);

	NumCompDesc++;

	return(0);
}




/*** function to read arrival label alias ***/

int GetLocAlias(char* line1)
{

	if (NumLocAlias >= MAX_NUM_LOC_ALIAS) {
		sprintf(MsgStr, "%s", line1);
		putmsg(1, MsgStr);
		sprintf(MsgStr,
"WARNING: maximum number of aliases reached, ignoring alias.");
		putmsg(1, MsgStr);
		return(-1);
	}

	sscanf(line1, "%s %s  %d %d %d  %d %d %d",
		LocAlias[NumLocAlias].name, LocAlias[NumLocAlias].alias,
		&(LocAlias[NumLocAlias].byr), &(LocAlias[NumLocAlias].bmo),
		&(LocAlias[NumLocAlias].bday),
		&(LocAlias[NumLocAlias].eyr), &(LocAlias[NumLocAlias].emo),
		&(LocAlias[NumLocAlias].eday));

	sprintf(MsgStr,
"LOCALIAS:  Name: %s  Alias: %s  Valid: %4.4d %2.2d %2.2d -> %4.4d %2.2d %2.2d",
		LocAlias[NumLocAlias].name, LocAlias[NumLocAlias].alias,
		LocAlias[NumLocAlias].byr, LocAlias[NumLocAlias].bmo,
		LocAlias[NumLocAlias].bday,
		LocAlias[NumLocAlias].eyr, LocAlias[NumLocAlias].emo,
		LocAlias[NumLocAlias].eday);
	putmsg(2, MsgStr);

	NumLocAlias++;

	return(0);
}



/*** function to read exclude arrival label and phase ***/

int GetLocExclude(char* line1)
{

	if (NumLocExclude >= MAX_NUM_LOC_EXCLUDE) {
		sprintf(MsgStr, "%s", line1);
		putmsg(1, MsgStr);
		sprintf(MsgStr,
"WARNING: maximum number of excluded phases reached, ignoring exclude.");
		putmsg(1, MsgStr);
		return(-1);
	}

	sscanf(line1, "%s %s  %d %d %d  %d %d %d",
		LocExclude[NumLocExclude].label, LocExclude[NumLocExclude].phase);

	sprintf(MsgStr, "LOCEXCLUDE:  Name: %s  Phase: %s",
		LocExclude[NumLocExclude].label, LocExclude[NumLocExclude].phase);
	putmsg(2, MsgStr);

	NumLocExclude++;

	return(0);
}




/*** function to read station phase time delays ***/

int GetTimeDelays(char* line1)
{

	if (NumTimeDelays >= MAX_NUM_STA_DELAYS) {
		sprintf(MsgStr, "%s", line1);
		putmsg(2, MsgStr);
		sprintf(MsgStr,
"WARNING: maximum number of station delays reached, ignoring alias.");
		putmsg(2, MsgStr);
		return(-1);
	}

	sscanf(line1, "%s %s %d %lf",
		TimeDelay[NumTimeDelays].label, TimeDelay[NumTimeDelays].phase,
		&(TimeDelay[NumTimeDelays].n_residuals),
		&(TimeDelay[NumTimeDelays].delay));

	sprintf(MsgStr,
"LOCDELAY:  Label: %s  Phase: %s  NumResiduals: %d  TimeDelay: %lf",
		TimeDelay[NumTimeDelays].label, TimeDelay[NumTimeDelays].phase,
		TimeDelay[NumTimeDelays].n_residuals,
		TimeDelay[NumTimeDelays].delay);
	putmsg(2, MsgStr);

	NumTimeDelays++;

	return(0);
}






/*** function to open summary output files */

int OpenSummaryFiles(char *path_output)
{

	int ngrid;
	char fname[FILENAME_MAX];



	for (ngrid = 0; ngrid < NumLocGrids; ngrid++) {

		if (!LocGridSave[ngrid])
			continue;

		/* Grid Hyp format */

		pSumFileHypNLLoc[ngrid] = NULL;
		sprintf(fname, "%s.sum.grid%d.loc.hyp",
						path_output, ngrid);
		if ((pSumFileHypNLLoc[ngrid] = fopen(fname, "w")) == NULL) {
			puterr2("ERROR: opening summary output file", fname);
			return(-1);
		} else {
			NumFilesOpen++;
		}


		iWriteHypHeader[ngrid] = 1;

		/* Hypo71 format */
		pSumFileHypo71[ngrid] = NULL;
		if (iSaveHypo71Sum) {
			sprintf(fname, "%s.sum.grid%d.loc.hypo_71",
						path_output, ngrid);
			if ((pSumFileHypo71[ngrid] = fopen(fname, "w"))
					== NULL) {
				puterr2(
"ERROR: opening HYPO71 summary output file",
					fname);
				return(-1);
			} else {
				NumFilesOpen++;
			}
			fprintf(pSumFileHypo71[ngrid], "%s\n",
				Hypocenter.comment);
		}


		/* HypoEllipse format */
		pSumFileHypoEll[ngrid] = NULL;
		if (iSaveHypoEllSum) {
			sprintf(fname, "%s.sum.grid%d.loc.hypo_ell",
						path_output, ngrid);
			if ((pSumFileHypoEll[ngrid] = fopen(fname, "w"))
					== NULL) {
				puterr2(
"ERROR: opening HypoEllipse summary output file",
					fname);
				return(-1);
			} else {
				NumFilesOpen++;
			}
			fprintf(pSumFileHypoEll[ngrid], "%s\n",
				Hypocenter.comment);
		}


		/* HypoInverse Archive format */
		pSumFileHypoInv[ngrid] = NULL;
		if (iSaveHypoInvSum) {
			sprintf(fname, "%s.sum.grid%d.loc.hypo_inv",
						path_output, ngrid);
			if ((pSumFileHypoInv[ngrid] = fopen(fname, "w"))
					== NULL) {
				puterr2(
"ERROR: opening HypoInverse Archive (FPFIT) summary output file",
					fname);
				return(-1);
			} else {
				NumFilesOpen++;
			}
		}

		/* Alberto 3D 4 chr sta SIMULPS format */
		pSumFileAlberto4[ngrid] = NULL;
		if (iSaveAlberto4Sum) {
			sprintf(fname, "%s.sum.grid%d.loc.sim",
						path_output, ngrid);
			if ((pSumFileAlberto4[ngrid] = fopen(fname, "w"))
					== NULL) {
				puterr2(
"ERROR: opening Alberto 3D, 4 chr sta, SIMULPS output file",
					fname);
				return(-1);
			} else {
				NumFilesOpen++;
			}
		}

	}


	return(0);

}



/*** function to close summary output files */

int CloseSummaryFiles()
{
	int ngrid;


	for (ngrid = 0; ngrid < NumLocGrids; ngrid++) {

		if (!LocGridSave[ngrid])
			continue;

		/* Grid Hyp format */

		if (pSumFileHypNLLoc[ngrid] != NULL) {
			fclose(pSumFileHypNLLoc[ngrid]);
			NumFilesOpen--;
		}

		/* Hypo71 format */

		if (pSumFileHypo71[ngrid] != NULL) {
			fclose(pSumFileHypo71[ngrid]);
			NumFilesOpen--;
		}

		/* HypoEll format */

		if (pSumFileHypoEll[ngrid] != NULL) {
			fclose(pSumFileHypoEll[ngrid]);
			NumFilesOpen--;
		}

		/* HypoInv format */

		if (pSumFileHypoInv[ngrid] != NULL) {
			fclose(pSumFileHypoInv[ngrid]);
			NumFilesOpen--;
		}

	}

	return(0);

}



/*** function to write hypocenter and arrivals to file (Alberto 4 SIMULPS) */

int WriteHypoAlberto4(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals, char* filename)
{

	int ifile = 0, narr;
	char fname[FILENAME_MAX];
	double mag;
	ArrivalDesc* parr;
	int nlat, nlon;
	double tpobs, resid;


	/* set hypocenter parameters */
	if (phypo->amp_mag != MAG_NULL)
		mag = phypo->amp_mag;
	else if (phypo->dur_mag != MAG_NULL)
		mag = phypo->dur_mag;
	else
		mag = 0.0;


	/* write hypocenter to file */

	if (fpio == NULL) {
		sprintf(fname, "%s.loc.sim", filename);
		if ((fpio = fopen(fname, "w")) == NULL) {
			puterr("ERROR: opening Alberto 4 hypocenter output file.");
			return(-1);
		} else {
			NumFilesOpen++;
		}
		ifile = 1;
	}


	/* write hypocenter parameters */

//87 1 1  1 2  0.00 35N14.26 120W44.00   5.15   0.00
	nlat = (int) fabs(phypo->dlat);
	nlon = (int) fabs(phypo->dlong);
	fprintf(fpio,
"%2.2d%2.2d%2.2d %2.2d%2.2d%6.2f %2.2d%c%5.2f %3.3d%c%5.2f %6.2f %6.2f",
		phypo->year % 100, phypo->month, phypo->day,
		phypo->hour, phypo->min, phypo->sec,
		nlat, phypo->dlat > 0.0 ? 'N' : 'S',
			60.0 * (fabs(phypo->dlat) - (double) nlat),
		nlon, phypo->dlong > 0.0 ? 'E' : 'W',
			60.0 * (fabs(phypo->dlong) - (double) nlon),
		phypo->depth, mag
	);
	
	// write arrivals

	for (narr = 0; narr < phypo->nreadings; narr++) {
		if (narr % 5 == 0)
			fprintf(fpio, "\n");
		parr = parrivals + narr;
		fprintf(fpio, "%4s%1s%1s%2.2d%7.4f",
			parr->label,
			strcmp(parr->onset, "?") == 0 ? "i" : parr->onset,
			parr->phase,
			parr->min, parr->sec
		);
	}

	fprintf(fpio, "\n");

	if (ifile) {
		fclose(fpio);
		NumFilesOpen--;
	}

}



/*** function to write hypocenter summary to file (quasi HypoEllipse format) */

int WriteHypoEll(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals,
		char* filename, int write_header, int write_arrivals)
{

	int ifile = 0, narr;
	char fname[FILENAME_MAX];
	double mag;
	ArrivalDesc* parr;
	double tpobs, resid;


	/* set hypocenter parameters */
	if (phypo->amp_mag != MAG_NULL)
		mag = phypo->amp_mag;
	else if (phypo->dur_mag != MAG_NULL)
		mag = phypo->dur_mag;
	else
		mag = 0.0;


	/* write hypocenter to file */

	if (fpio == NULL) {
		sprintf(fname, "%s.loc.hypo_ell", filename);
		if ((fpio = fopen(fname, "w")) == NULL) {
			puterr("ERROR: opening hypocenter output file.");
			return(-1);
		} else {
			NumFilesOpen++;
		}
		ifile = 1;
	}


	/* write hypocenter parameters */

	if (write_header) {
		fprintf(fpio,
		    "DATE     ORIGIN     LAT         LONG         DEPTH   ");
		fprintf(fpio,
			"MAG  NO  GAP D1     RMS   ");
		fprintf(fpio,
			"AZ1  DIP1 SE1    AZ2  DIP2 SE2    SE3    \n");
/*"ERH  ERZ Q SQD  ADJ IN NR  AVR  AAR NM AVXM SDXM NF AVFM SDFM I\n");*/
	}
	fprintf(fpio,
"%4.4d%2.2d%2.2d %2.2d%2.2d %5.2lf %3d %1c %5.2lf %4d %1c %5.2lf %7.3lf ",
		phypo->year, phypo->month, phypo->day,
		phypo->hour, phypo->min, phypo->sec,
		(int) fabs(phypo->dlat), (phypo->dlat >= 0.0 ? 'N' : 'S'),
		(fabs(phypo->dlat) - (int) fabs(phypo->dlat)) * 60.0,
		(int) fabs(phypo->dlong), (phypo->dlong >= 0.0 ? 'E' : 'W'),
		(fabs(phypo->dlong) - (int) fabs(phypo->dlong)) * 60.0,
		phypo->depth);
	fprintf(fpio, "%4.2lf %3d %3d %6.2lf %5.2lf ",
		mag, phypo->nreadings, phypo->gap, phypo->dist, phypo->rms);
	fprintf(fpio, "%4d %4d %6.2lf %4d %4d %6.2lf %6.2lf ",
		(int) (0.5 + phypo->ellipsoid.az1),
			(int) (0.5 + phypo->ellipsoid.dip1),
			phypo->ellipsoid.len1,
		(int) (0.5 + phypo->ellipsoid.az2),
			(int) (0.5 + phypo->ellipsoid.dip2),
			phypo->ellipsoid.len2,
		phypo->ellipsoid.len3);
	fprintf(fpio, "\n");



    if (write_arrivals) {

	fprintf(fpio, "\n");


	fprintf(fpio,
"  STN  DIST AZM AIN PRMK HRMN P-SEC TPOBS TPCAL DLY/H1 P-RES P-WT AMX PRX CALX K XMAG RMK FMP FMAG\n");
/*"  STN  DIST AZM AIN PRMK HRMN P-SEC TPOBS TPCAL DLY/H1 P-RES P-WT AMX PRX CALX K XMAG RMK FMP FMAG SRMK S-SEC TSOBS S-RES  S-WT    DT\n"); */

	for (narr = 0; narr < phypo->nreadings; narr++) {
		parr = parrivals + narr;
		tpobs = parr->obs_travel_time > -9.99 ? parr->obs_travel_time : 0.0;
		resid = parr->residual > -99.99 ? parr->residual : -99.99;
		fprintf(fpio,
"%5s %5.1lf %3d %3d %2s%1s%1d %2.2d%2.2d %5.2lf %5.2lf %5.2lf       %-6.2lf %5.2lf\n",
			parr->label, parr->dist,
			(int) (0.5 + rect2latlonAngle(0, parr->ray_azim)),
			(int) (0.5 + parr->ray_dip),
			parr->phase, parr->first_mot, parr->quality,
			parr->hour, parr->min, parr->sec,
			tpobs, parr->pred_travel_time,
			resid, parr->weight
		);
	}

    }


	if (ifile) {
		fclose(fpio);
		NumFilesOpen--;
	}

}




/*** function to write hypocenter summary to file (HYPO71 format) */

int WriteHypo71(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals,
		char* filename, int write_header, int write_arrivals)
{

	int ifile = 0, narr;
	char fname[FILENAME_MAX];
	ArrivalDesc* parr;
	double tpobs, resid;
	double mag, erh, erz;
	char qualS, qualD, qual;
	int pha_qual;
	double xmag, fmag;

	/* set hypocenter parameters */
	if (phypo->amp_mag != MAG_NULL)
		mag = phypo->amp_mag;
	else if (phypo->dur_mag != MAG_NULL)
		mag = phypo->dur_mag;
	else
		mag = 0.0;

	/* write hypocenter to file */

	if (fpio == NULL) {
		sprintf(fname, "%s.loc.h71", filename);
		if ((fpio = fopen(fname, "w")) == NULL) {
			puterr("ERROR: opening hypocenter output file.");
			return(-1);
		} else {
			NumFilesOpen++;
		}
		ifile = 1;
	}


	/* write hypocenter parameters */

	if (write_header) {
		fprintf(fpio,
"  DATE    ORIGIN    LAT      LONG      DEPTH    ");
		fprintf(fpio,
"MAG NO DM GAP M  RMS  ERH  ERZ Q SQD  ADJ IN NR  AVR  AAR NM AVXM SDXM NF AVFM SDFM I\n");
	}

	fprintf(fpio,
	   " %2.2d%2.2d%2.2d %2.2d%2.2d %5.2lf%3d %5.2lf%4d %5.2lf %6.2lf",
		phypo->year % 100, phypo->month, phypo->day,
		phypo->hour, phypo->min, phypo->sec,
		(int) phypo->dlat, /*(phypo->dlat >= 0.0 ? 'N' : 'S'),*/
		(phypo->dlat - (int) phypo->dlat) * 60.0,
		(int) phypo->dlong, /*(phypo->dlong >= 0.0 ? 'E' : 'W'),*/
		(phypo->dlong - (int) phypo->dlong) * 60.0,
		phypo->depth);

	fprintf(fpio, " %6.2lf%3d%3d %3d 0%5.2lf",
		mag, phypo->nreadings, (int) (0.5 + phypo->dist),
				phypo->gap, phypo->rms);

	erh = sqrt(phypo->cov.xx * phypo->cov.xx
				+ phypo->cov.yy * phypo->cov.yy);
	erz = sqrt(phypo->cov.zz * phypo->cov.zz);
	fprintf(fpio, "%5.1lf%5.1lf", erz, erh);

	/* ABCD quality levels (from HYPO71 open file report 75-311, p. 27) */
	if (     phypo->rms < 0.15 && erh <= 1.0 && erz <= 2.0)
		qualS = 'A';
	else if (phypo->rms < 0.30 && erh <= 2.5 && erz <= 5.0)
		qualS = 'B';
	else if (phypo->rms < 0.50 && erh <= 5.0)
		qualS = 'C';
	else
		qualS = 'D';
	if (     phypo->nreadings >= 6 && phypo->gap <= 90
			&& (phypo->dist <= phypo->depth || phypo->dist <= 5.0))
		qualD = 'A';
	else if (phypo->nreadings >= 6 && phypo->gap <= 135
			&& (phypo->dist <= 2.0 * phypo->depth
						|| phypo->dist <= 10.0))
		qualD = 'B';
	else if (phypo->nreadings >= 6 && phypo->gap <= 180
			&& phypo->dist <= 50.0)
		qualD = 'C';
	else
		qualD = 'D';
	/*if (abs(qualS - qualD) == 1)
		qual = qualS < qualD ? qualS : qualD;
	else*/
	qual = 1 + (qualS + qualD) / 2;
	fprintf(fpio, " %1c %1c %1c", qual, qualS, qualD);

	/* dummy values for remaining fields */
	fprintf(fpio,
" %4.2lf %2d %2d-%4.2lf %4.2lf %2d %4.1lf %4.1lf %2d %4.1lf %4.1lf%2d\n",
		0.0, 0, 0, 0.0, 0.0, 0, 0.0, 0.0, 0, 0.0, 0.0, 0);



    if (write_arrivals) {

	fprintf(fpio, "\n");


	fprintf(fpio,
"  STN  DIST AZM AIN PRMK HRMN P-SEC TPOBS TPCAL DLY/H1 P-RES P-WT AMX PRX CALX K XMAG RMK FMP FMAG SRMK S-SEC TSOBS S-RES  S-WT    DT\n");

	for (narr = 0; narr < phypo->nreadings; narr++) {
		parr = parrivals + narr;
		pha_qual = (parr->quality >= 0 && parr->quality <= 4) ?
				parr->quality : Err2Qual(parr);
		if (pha_qual < 0)
			pha_qual = 4;
		tpobs = parr->obs_travel_time > -9.99 ?
					parr->obs_travel_time : 0.0;
		resid = parr->residual > -99.99 ? parr->residual : -99.99;
		fprintf(fpio,
/*" %-4s %5.1lf %3d %3d %2s%1s%1d %2.2d%2.2d %5.2lf %5.2lf %5.2lf  0.00 %6.3lf %5.2lf\n", */
" %-4s %5.1lf %3d %3d %2s%1s%1d %2.2d%2.2d %5.2lf %5.2lf %5.2lf  0.00 %-6.2lf %4.2lf",
			parr->label, parr->dist,
			(int) (0.5 + rect2latlonAngle(0, parr->ray_azim)),
				(int) (0.5 + parr->ray_dip),
			parr->phase, parr->first_mot, pha_qual,
			parr->hour, parr->min, parr->sec,
			tpobs, parr->pred_travel_time,
			resid, parr->weight
		);

		/* set magnitudes */
		xmag = (parr->amp_mag != MAG_NULL) ? parr->amp_mag : 0.0;
		fmag = (parr->dur_mag != MAG_NULL) ? parr->dur_mag : 0.0;
		fprintf(fpio,
" 0.0 0.0 0.00 0 %3.2lf 000 00.0 %3.2lf ??4 00.00 00.00 00.00   0.0      \n",
			xmag, fmag
		);
	}

    }


	if (ifile) {
		fclose(fpio);
		NumFilesOpen--;
	}

}




/*** function to write hypocenter summary to file (HYPOINVERSE archive format) */

int WriteHypoInverseArchive(FILE *fpio, HypoDesc *phypo, ArrivalDesc *parrivals, int narrivals,
		char *filename, int write_header, int write_arrivals)
{

	int ifile = 0, narr;
	char fname[FILENAME_MAX];
	ArrivalDesc *parr, *psarr;
	double resid;

	/* hypocenter dummies */
	char loc_remark[] = "XXX", aux_remark[] = "  ";
	int num_S_wt = 0, num_P_fmot = 0;
	double xmag1, fmag1;
	double err_horiz = 0.0, err_vert = 0.0;
	/* arrival dummies */
	char sta_remark[] = " ";
	int amp_mag_wt_code = 0, dur_mag_wt_code = 0;
	double dev_pp_amp = 0.0, per_amp = 0.0, coda_duration = 0.0;
	double dur_mag = 0.0, amp_mag = 0.0;
	double p_imporance = 0.0, s_importance = 0.0;



	/* set hypocenter parameters */
//printf("Ma %4.2f  Md %4.2f  MAG_NULL %4.2f\n", phypo->amp_mag, phypo->dur_mag, MAG_NULL);
	if (fabs(phypo->amp_mag - MAG_NULL) > 0.01)
		xmag1 = phypo->amp_mag;
	else
		xmag1 = 0.0;
	if (fabs(phypo->dur_mag - MAG_NULL) > 0.01)
		fmag1 = phypo->dur_mag;
	else
		fmag1 = xmag1;
//printf("-> xmag1 %4.2f  fmag1 %4.2f\n", xmag1, fmag1);



	/* write hypocenter to file */

	if (fpio == NULL) {
		sprintf(fname, "%s.loc.h71", filename);
		if ((fpio = fopen(fname, "w")) == NULL) {
			puterr("ERROR: opening hypocenter output file.");
			return(-1);
		} else {
			NumFilesOpen++;
		}
		ifile = 1;
	}


	/* write hypocenter parameters */

	fprintf(fpio,
	   "%2.2d%2.2d%2.2d%2.2d%2.2d%4.0lf%2.2d%1c%4.0lf%3.3d%1c%4.0lf%5.0lf",
		phypo->year % 100, phypo->month, phypo->day,
		phypo->hour, phypo->min, 100.0 * phypo->sec,
		(int) fabs(phypo->dlat), (phypo->dlat >= 0.0 ? ' ' : 'S'),
		(fabs(phypo->dlat) - (int) fabs(phypo->dlat)) * 6000.0,
		(int) fabs(phypo->dlong), (phypo->dlong >= 0.0 ? 'E' : ' '),
		(fabs(phypo->dlong) - (int) fabs(phypo->dlong)) * 6000.0,
		100.0 * phypo->depth);
	fprintf(fpio, "%2.0lf%3.3d%3.3d%3.0lf%4.0lf",
		10.0 * xmag1, phypo->nreadings, phypo->gap, phypo->dist,
		100.0 * phypo->rms);
/*	fprintf(fpio, "%3.3d%2.2d%4.0lf%3.3d%2.2d%4.0lf",
		(int) (0.5 + phypo->ellipsoid.az1),
		90 + (int) (0.5 + phypo->ellipsoid.dip1),
		100.0 * phypo->ellipsoid.len1,
		(int) (0.5 + phypo->ellipsoid.az2),
		90 + (int) (0.5 + phypo->ellipsoid.dip2),
		100.0 * phypo->ellipsoid.len2);
*/
	fprintf(fpio, "000000000000000000");
	fprintf(fpio, "%2.0lf%3s%4.0lf%2s%2.2d%4.0lf%4.0lf%2.2d",
		10.0 * fmag1, loc_remark, 100.0 * phypo->ellipsoid.len3,
		aux_remark, num_S_wt,
		err_horiz, err_vert, num_P_fmot);
	fprintf(fpio, "                   D   X");

	fprintf(fpio, "\n");



    if (write_arrivals) {

// AJL 21JUN2000 bug fix (not all phasese are written,
//   looses fm's from phases not used for misfit!
//	for (narr = 0; narr < phypo->nreadings; narr++) {
	for (narr = 0; narr < narrivals; narr++) {
		parr = parrivals + narr;

		/* skip non-P phases */
		if (!IsPhaseID(parr->phase, "P"))
			continue;

		/* skip phases with take-off angle quality < min qual*/
		if (parr->ray_qual < iAngleQualityMin)
			continue;

		/* check for following S arrival at same station */
		psarr = NULL;
// AJL 21JUN2000
//		if (narr + 1 < phypo->nreadings) {
		if (narr + 1 < narrivals) {
			psarr = parr + 1;
			if (!IsPhaseID(psarr->phase, "S")
				    || strcmp(parr->label, psarr->label) != 0)
				psarr = NULL;
		}


		fprintf(fpio,
			"%4s%1s%1s%1s%1d%1s",
			parr->label, " ", parr->phase, parr->first_mot,
			parr->quality, " ");
		fprintf(fpio,
			"%2.2d%2.2d%2.2d%2.2d%2.2d%5.0lf",
			parr->year % 100, parr->month, parr->day,
			parr->hour, parr->min, 100.0 * parr->sec);
		resid = parr->residual > -9.99 ? parr->residual : 0.0;
		fprintf(fpio,
			"%4.0lf%3.0lf",
			100.0 * resid, 100.0 * parr->weight);
		if (psarr != NULL) {
			resid = psarr->residual > -9.99 ? psarr->residual : 0.0;
			fprintf(fpio,
			"%5.0lf%1s%1s %1d%4.0lf%3.0lf%3.0lf%4.0lf%4.0lf",
				100.0 * psarr->sec, " ", psarr->phase,
				psarr->quality, 100.0 * resid,
				dev_pp_amp, 100.0 * psarr->weight,
				100.0 * parr->delay, 100.0 * psarr->delay);
		} else {
			fprintf(fpio,
			"                   %4.0lf%4.0lf",
				100.0 * parr->delay, 0.0);
		}

		fprintf(fpio,
			"%4.0lf%3d%1d%1d%3.0lf%1s%4.0lf%3d",
			10.0 * parr->dist, (int) (0.5 + parr->ray_dip),
			amp_mag_wt_code, dur_mag_wt_code, 100.0 * per_amp,
			sta_remark, coda_duration,
			(int) (0.5 + rect2latlonAngle(0, parr->ray_azim)));
		fprintf(fpio,
			"%2.0lf%2.0lf%4.0lf%4.0lf",
			10.0 * dur_mag, 10.0 * amp_mag,
			1000.0 * p_imporance, 1000.0 * s_importance);


		fprintf(fpio, "\n");
	}

    }


	fprintf(fpio, "\n");
	if (ifile) {
		fclose(fpio);
		NumFilesOpen--;
	}
}



/*** function to determine azimuth gap -
			largest azimuth separation between
			stations as seen from epicenter (deg) */

double CalcAzimuthGap(ArrivalDesc *arrival, int num_arrivals)
{

	int narr;
	double az_last, az, gap, gap_max = -1.0;
	static double azimuth[MAX_NUM_ARRIVALS];

	/* load azimuths to working array */
	for (narr = 0; narr < num_arrivals; narr++) {
		azimuth[narr] = (arrival + narr)->azim;
	}

	/* sort */
	SortDoubles(azimuth, num_arrivals);

	/* find largest gap */
	az_last = azimuth[0];
	for (narr = 1; narr < num_arrivals; narr++) {
		az = azimuth[narr];
		gap = az - az_last;
		if (gap > gap_max)
			gap_max = gap;
		az_last = az;
	}
	az = azimuth[0] + 360.0;
	gap = az - az_last;
	if (gap > gap_max)
		gap_max = gap;

	return(gap_max);

}




/*** function to estimate Vp/Vs ration */

/*  follows chapter 5 of Lahr, J.C., 1989. HYPOELLIPSE/Version 2.0: A computer program for
determining local earthquake hypocentral parameters, magnitude and first motion pattern, U.S.
Geological Survey Open-File Report 89-116, 92p.
*/

double CalculateVpVsEstimate(HypoDesc* phypo, ArrivalDesc* arrival, int narrivals)
{

	int narr, npair;
	int k;
	double az_last, az, gap, gap_max = -1.0;
	double p_time[MAX_NUM_ARRIVALS], p_time_wt;
	double p_time_cent, p_error[MAX_NUM_ARRIVALS];
	double s_time[MAX_NUM_ARRIVALS], s_time_wt;
	double s_time_cent, s_error[MAX_NUM_ARRIVALS];
	double weight[MAX_NUM_ARRIVALS], weight_sum;
	double B, Btest, Bmin, dB, Tmin, T, temp;
	double tsp, tsp_min_max_diff;
	double tsp_min = VERY_LARGE_DOUBLE;
	double tsp_max = -VERY_LARGE_DOUBLE;

//printf("CALCULATING VpVs\n");

	/* load P S pairs to working arrays (assumes arrivals sorted by distance) */

	npair = 0;
	for (narr = 1; narr < narrivals; narr++) {
		if (strcmp(arrival[narr - 1].label, arrival[narr].label) == 0
				&& IsPhaseID(arrival[narr - 1].phase, "P")
				&& IsPhaseID(arrival[narr].phase, "S")) {
// DELAY_CORR
			p_time[npair] = arrival[narr - 1].obs_time + arrival[narr - 1].delay;
			p_error[npair] = arrival[narr - 1].error;
// DELAY_CORR
			s_time[npair] = arrival[narr].obs_time + arrival[narr].delay;
			s_error[npair] = arrival[narr].error;
//printf("PAIR %s %s (%lf +/- %lf) + %s %s (%lf +/- %lf)\n", arrival[narr - 1].label, arrival[narr - 1].phase, p_time[npair], p_error[npair], arrival[narr].label, arrival[narr].phase, s_time[npair], s_error[npair]);

			tsp = s_time[npair] - p_time[npair];
			tsp_min = tsp < tsp_min ? tsp : tsp_min;
			tsp_max = tsp > tsp_max ? tsp : tsp_max;

			npair++;
		}
	}

	tsp_min_max_diff = tsp_max - tsp_min;
	phypo->tsp_min_max_diff = tsp_min_max_diff;

	/* not enough pairs found */
	if (npair < 2) {
		phypo->VpVs = -1.0;
		phypo->nVpVs = npair;
		return(-1.0);
	}


	/* search for optimal VpVs */

	B = 2.0;
	Tmin = VERY_LARGE_DOUBLE;
	for (dB = 0.5; dB > 0.00001; dB *= 0.4) {

		for (k = -3; k < 4; k++) {

			Btest = B + (double) k * dB;

			// form weights
			p_time_wt = 0.0;
			s_time_wt = 0.0;
			weight_sum = 0.0;
			for (narr = 1; narr < npair; narr++) {
				weight[narr] = 1.0 /
					(s_error[narr] * s_error[narr]
						+ Btest * p_error[narr] * p_error[narr]);
				p_time_wt += weight[narr] * p_time[narr];
				s_time_wt += weight[narr] * s_time[narr];
				weight_sum += weight[narr];
			}

			// form centered times and accumulate T
			T = 0.0;
			for (narr = 1; narr < npair; narr++) {
				p_time_cent = p_time[narr] - p_time_wt / weight_sum;
				s_time_cent = s_time[narr] - s_time_wt / weight_sum;
				temp = s_time_cent - Btest * p_time_cent;
				T += weight[narr] * temp * temp;
			}

			if (T < Tmin) {
				Tmin = T;
				Bmin = Btest;
//printf("  NEW MIN VpVs = %lf dB = %lf T = %le  k = %d\n", Btest, dB, T, k);
			}

		}

		B = Bmin;
	}


//printf("  OPTIMAL VpVs = %.3lf last dB = %lf npair = %d\n", B, dB / 0.4, npair);

	phypo->VpVs = B;
	phypo->nVpVs = npair;

	return(B);

}




/*** function to calculate magintude */

int CalculateMagnitude(HypoDesc* phypo, ArrivalDesc* parrivals,
		int narrivals, CompDesc* pcomp, int nCompDesc, MagDesc* pmagnitude)
{

	int nmag , narr, nComp;
	double amp_mag_sum, dur_mag_sum;
	double amp_fact_ml_hb, sta_corr;
	ArrivalDesc* parr;



	/* calculate magnitude for specified calculation type */

	if (pmagnitude == MAG_UNDEF) {
		/* no magnitude calculation type given */

		return(0);

	}
	else if (pmagnitude->type == MAG_ML_HB) {
		/* ML - Hutton & Boore, BSSA, v77, n6, Dec 1987 */

	/* write message */
	sprintf(MsgStr, "\nComponent results for: ML - Hutton & Boore, BSSA, v77, n6, Dec 1987:");
	putmsg(2, MsgStr);

		nmag = 0;
		amp_mag_sum = 0.0;
		for (narr = 0; narr < narrivals; narr++) {

			/* check for valid amplitude */
			if ((parr = parrivals + narr)->amplitude > 0.0 && parr->amplitude != AMPLITUDE_NULL) {

				nComp = findStaInstComp(parr, pcomp, nCompDesc);
				if (nComp >= 0) {
					/* sta/inst/comp found */
					amp_fact_ml_hb = (pcomp + nComp)->amp_fact_ml_hb;
					sta_corr = (pcomp + nComp)->sta_corr_ml_hb;
				} else {
					amp_fact_ml_hb = 1.0;
					sta_corr = 0.0;
				}

				/* calc ML */
				parr->amp_mag = Calc_ML_HuttonBoore(
					parr->amplitude * pmagnitude->amp_fact_ml_hb * amp_fact_ml_hb,
					parr->dist, phypo->depth, sta_corr,
					pmagnitude->hb_n, pmagnitude->hb_K,
					pmagnitude->hb_Ro, pmagnitude->hb_Mo);

	/* write message */
	sprintf(MsgStr, "%s %s %s amp %.2e f %.2e f_sta %.2e dist %.2f depth %.2f sta_corr %.4f hb_n %.2f hb_K %.5f mag %.2f",
parr->label, parr->inst, parr->comp, parr->amplitude, pmagnitude->amp_fact_ml_hb, amp_fact_ml_hb, parr->dist,
phypo->depth, sta_corr, pmagnitude->hb_n, pmagnitude->hb_K, parr->amp_mag);
	putmsg(2, MsgStr);


				/* update event magnitude */
				amp_mag_sum += parr->amp_mag;
				nmag++;
			}
		}
		if (nmag > 0) {
			phypo->amp_mag = amp_mag_sum / (double) nmag;
			phypo->num_amp_mag = nmag;
		}
	}
	else if (pmagnitude->type == MAG_MD_FMAG) {
		/* coda duration (FMAG) - HYPOELLIPSE users manual chap 4;
			Lee et al., 1972; Lahr et al., 1975; Bakun and Lindh, 1977 */

		nmag = 0;
		dur_mag_sum = 0.0;
		for (narr = 0; narr < narrivals; narr++) {

			/* check for valid amplitude */
			if ((parr = parrivals + narr)->coda_dur > 0.0 && parr->coda_dur != CODA_DUR_NULL) {

				nComp = findStaInstComp(parr, pcomp, nCompDesc);
				if (nComp >= 0) {
					/* sta/inst/comp found */
					sta_corr = (pcomp + nComp)->sta_corr_md_fmag;
				} else {
					sta_corr = 1.0;
				}

				/* calc ML */
				parr->dur_mag = Calc_MD_FMAG(
					parr->coda_dur, parr->dist, phypo->depth, sta_corr,
					pmagnitude->fmag_c1, pmagnitude->fmag_c2, pmagnitude->fmag_c3,
					pmagnitude->fmag_c4, pmagnitude->fmag_c5);

/*printf("%s %s %s coda_dur %lf dist %lf depth %lf sta_corr %lf c1 %lf c2 %lf c3 %lf c4 %lf c5 %lf mag %lf\n",
parr->label, parr->inst, parr->comp, parr->coda_dur, parr->dist,
phypo->depth, sta_corr, pmagnitude->fmag_c1, pmagnitude->fmag_c2, pmagnitude->fmag_c3,
pmagnitude->fmag_c4, pmagnitude->fmag_c5, parr->dur_mag);
*/

				/* update event magnitude */
				dur_mag_sum += parr->dur_mag;
				nmag++;
			}
		}
		if (nmag > 0) {
			phypo->dur_mag = dur_mag_sum / (double) nmag;
			phypo->num_dur_mag = nmag;
		}
	}


}


/*** function to find component parameters */

int findStaInstComp(ArrivalDesc* parr, CompDesc* pcomp, int nCompDesc) {

	int nComp;
	char *pchr, test_label[ARRIVAL_LABEL_LEN];
	
	strcpy(test_label, parr->time_grid_label);

	for (nComp = 0; nComp < nCompDesc; nComp++) {
		strcpy(test_label, parr->time_grid_label);
		if ((pchr = strrchr(test_label, '_')) != NULL)
			*pchr = '\0';
//printf("comp %s  arr_test %s %s %s\n", (pcomp + nComp)->label, parr->label, parr->time_grid_label, test_label);
		if ((pcomp + nComp)->label[0] != '?' &&
				(pcomp + nComp)->label[0] != '*' &&
				strcmp((pcomp + nComp)->label,
				test_label) != 0)
			continue;
		if ((pcomp + nComp)->inst[0] != '?' &&
				(pcomp + nComp)->inst[0] != '*' &&
				strcmp((pcomp + nComp)->inst,
				parr->inst) != 0)
			continue;
		if ((pcomp + nComp)->comp[0] != '?' &&
				(pcomp + nComp)->comp[0] != '*' &&
				strcmp((pcomp + nComp)->comp,
				parr->comp) != 0)
			continue;
		// found
		return(nComp);
	}

	// not found
	return(-1);
}


/*** function to calculate local magnitude ML */
/* ML - Hutton & Boore, BSSA, v77, n6, Dec 1987 */

double Calc_ML_HuttonBoore(double amplitude, double dist, double depth, double sta_corr,
			double hb_n, double hb_K, double hb_Ro, double hb_Mo)
{

	double hyp_dist, magnitude;

	hyp_dist = sqrt(dist * dist + depth * depth);

	if (hyp_dist < SMALL_DOUBLE)
		return(MAG_NULL);

	magnitude = log10(amplitude)
			+ hb_n * log10(hyp_dist / hb_Ro) + hb_K * (hyp_dist - hb_Ro)
			+ hb_Mo + sta_corr;

	return(magnitude);

}



/*** function to calculate coda duration magnitude MD */
/* coda duration (FMAG) - HYPOELLIPSE users manual chap 4;
	Lee et al., 1972; Lahr et al., 1975; Bakun and Lindh, 1977 */

double Calc_MD_FMAG(double coda_dur, double dist, double depth, double sta_corr,
	double fmag_c1, double fmag_c2, double fmag_c3, double fmag_c4, double fmag_c5)
{

	double hyp_dist, magnitude;

	hyp_dist = sqrt(dist * dist + depth * depth);

	if (coda_dur * sta_corr < SMALL_DOUBLE)
		return(MAG_NULL);

	magnitude = fmag_c1
			+ fmag_c2 * log10(coda_dur * sta_corr)
			+ fmag_c3 * dist
			+ fmag_c4 * depth
			+ fmag_c5 * pow(log10(coda_dur * sta_corr), 2);

	return(magnitude);

}






/*------------------------------------------------------------/ */
/** hashtable routines for accumulating station statistics */

/* from Kernigham and Ritchie, C prog lang, 2nd ed, 1988, sec 6.6 */


/*** funtion to form ordered hash value from firt char of label */

unsigned hash(char* label, char* phase)
{
	unsigned hashval;

	if (isdigit(label[0]))
		hashval = label[0] - '0';
	else if (isalpha(label[0]))
		hashval = 10 + toupper(label[0]) - 'A';
	else
		hashval = 36 + label[0] % 10;

	hashval = hashval % HASHSIZE;

	return hashval;
}


/*** function to lookup labelphase in hashtable */

StaStatNode *lookup(int ntable, char* label, char* phase)
{
	StaStatNode *np;

	for (np = hashtab[ntable][hash(label, phase)]; np != NULL;
							np = np->next)
		if (strcmp(label, np->label) == 0
				&& strcmp(phase, np->phase) == 0)
			return(np);	/* found */

	return(NULL);		/* not found */

}


/*** function to install or add (labelphase, residual, weight) in hashtable */

StaStatNode *InstallStaStatInTable(int ntable, char* label, char* phase, int flag_ignore,
	double residual, double weight,
	double pdf_residual_sum, double pdf_weight_sum, double delay)
{
	int icomp;
	StaStatNode *np, *npcheck, *nplast;
	unsigned hashval;

	if ((np = lookup(ntable, label, phase)) == NULL)
	{
		/* not found, create new StaStatNode */
		if ((np = (StaStatNode *) malloc(sizeof(StaStatNode))) == NULL)
			return(NULL);
		strcpy(np->label, label);
		strcpy(np->phase, phase);
		np->flag_ignore = flag_ignore;
		np->residual_min = residual;
		np->residual_max = residual;
		np->residual_sum = residual * weight;
		np->residual_square_sum = residual * residual * weight;
		np->weight_sum = weight;
		np->num_residuals = 1;
		if (pdf_weight_sum > VERY_SMALL_DOUBLE) {
			np->pdf_residual_sum =
				pdf_residual_sum / pdf_weight_sum;
			np->pdf_residual_square_sum =
				pdf_residual_sum * pdf_residual_sum / (pdf_weight_sum * pdf_weight_sum);
			np->num_pdf_residuals = 1;
		} else {
			np->num_pdf_residuals = 0;
		}
		np->delay = delay;

		/* put in table in alphabetical order */
		hashval = hash(label, phase);
		nplast = NULL;
		npcheck = hashtab[ntable][hashval];
		while (npcheck != NULL &&
				(icomp = strcmp(npcheck->label, label)) <= 0) {
			if (icomp == 0 && strcmp(npcheck->phase, phase) >= 0)
				break;
			nplast = npcheck;
			npcheck = npcheck->next;
		}
		np->next = npcheck;
		if (nplast != NULL)
			nplast->next = np;
		else
			hashtab[ntable][hashval] = np;
	}
	else
	{
		/* already there */
		if (residual < np->residual_min)
			np->residual_min = residual;
		if (residual > np->residual_max)
			np->residual_max = residual;
		np->residual_sum += residual * weight;
		np->residual_square_sum += residual * residual * weight;
		np->weight_sum += weight;
		np->num_residuals++;
		if (pdf_weight_sum > VERY_SMALL_DOUBLE) {
			np->pdf_residual_sum +=
				pdf_residual_sum / pdf_weight_sum;
			np->pdf_residual_square_sum +=
				pdf_residual_sum * pdf_residual_sum / (pdf_weight_sum * pdf_weight_sum);
			np->num_pdf_residuals++;
		}
	}

	return(np);

}


/*** function to output hashtable values */

int WriteStaStatTable(int ntable, FILE *fpio,
		double rms_max, int nRdgs_min, double gap_max,
		double p_residual_max, double s_residual_max, int imode)
{
	int nnodes;
	unsigned hashval;
	char frmt1[MAXLINE], frmt2[MAXLINE];
	double res_temp, res_std_temp;
	StaStatNode *np;

	sprintf(frmt1, "LOCDELAY  %%-%ds %%-%ds %%-8d %%-12lf\n",
			ARRIVAL_LABEL_LEN, ARRIVAL_LABEL_LEN);
	sprintf(frmt2, "LOCDELAY  %%-%ds %%-%ds %%-8d %%-12lf %%-12lf %%-12lf %%-12lf %%d\n",
			ARRIVAL_LABEL_LEN, ARRIVAL_LABEL_LEN);

	if (imode == WRITE_RESIDUALS) {
		fprintf(fpio,
"\nAverage Phase Residuals (CalcResidual)  RMS_Max: %lf  NRdgs_Min: %d  Gap_Max: %lf  P_Res_Max: %lf  S_Res_Max: %lf\n",
			rms_max, nRdgs_min, gap_max,
			p_residual_max, s_residual_max);
		fprintf(fpio,
"          ID      Phase   Nres      AveRes       StdDev       ResMin       ResMax     ignored\n");
	} else if (imode == WRITE_RES_DELAYS) {
		fprintf(fpio,
"\nTotal Phase Corrections (CalcResidual + InputDelay)  RMS_Max: %lf  NRdgs_Min: %d  Gap_Max: %lf  P_Res_Max: %lf  S_Res_Max: %lf\n",
			rms_max, nRdgs_min, gap_max,
			p_residual_max, s_residual_max);
		fprintf(fpio,
			"          ID      Phase   Nres      TotCorr\n");
	} else if (imode == WRITE_PDF_RESIDUALS) {
		fprintf(fpio,
"\nAverage Phase Residuals PDF (CalcPDFResidual)  RMS_Max: %lf  NRdgs_Min: %d  Gap_Max: %lf  P_Res_Max: %lf  S_Res_Max: %lf\n",
			rms_max, nRdgs_min, gap_max,
			p_residual_max, s_residual_max);
		fprintf(fpio,
"          ID      Phase   Nres      AveRes       StdDev       ResMin       ResMax     ignored\n");
	} else if (imode == WRITE_PDF_DELAYS) {
		fprintf(fpio,
"\nTotal Phase Corrections PDF (CalcPDFResidual + InputDelay)  RMS_Max: %lf  NRdgs_Min: %d  Gap_Max: %lf  P_Res_Max: %lf  S_Res_Max: %lf\n",
			rms_max, nRdgs_min, gap_max,
			p_residual_max, s_residual_max);
		fprintf(fpio,
			"          ID      Phase   Nres      TotCorr\n");
	}

	for (hashval = 0; hashval < HASHSIZE; hashval++) {
		for (np = hashtab[ntable][hashval]; np != NULL; np = np->next) {
			if (imode == WRITE_RESIDUALS) {
				res_temp = np->residual_sum / np->weight_sum;
				res_std_temp = np->residual_square_sum / np->weight_sum
							- res_temp * res_temp;
				if (np->num_residuals > 1)
					res_std_temp = sqrt(np->residual_square_sum / np->weight_sum
								- res_temp * res_temp);
				else
					res_std_temp = -1.0;
				fprintf(fpio, frmt2, np->label, np->phase,
					np->num_residuals, res_temp, res_std_temp,
					np->residual_min, np->residual_max, np->flag_ignore);
			} else if (imode == WRITE_RES_DELAYS) {
				res_temp = np->delay +
					np->residual_sum / np->weight_sum;
				fprintf(fpio, frmt1, np->label, np->phase,
					np->num_residuals, res_temp);
			} else if (imode == WRITE_PDF_RESIDUALS) {
				if (np->num_pdf_residuals > 0) {
					res_temp = np->pdf_residual_sum
						/ (double) np->num_pdf_residuals;
				} else {
					res_temp = 0.0;
				}
				if (np->num_pdf_residuals > 1)
					res_std_temp = sqrt(np->pdf_residual_square_sum
						/ (double) (np->num_pdf_residuals - 1)
							- res_temp * res_temp);
				else
					res_std_temp = -1.0;
				fprintf(fpio, frmt2, np->label, np->phase,
					np->num_pdf_residuals, res_temp, res_std_temp,
					np->residual_min, np->residual_max, np->flag_ignore);
			} else if (imode == WRITE_PDF_DELAYS) {
				if (np->num_pdf_residuals > 0) {
					res_temp = np->pdf_residual_sum
					    / (double) np->num_pdf_residuals;
				} else
					res_temp = 0.0;
				res_temp += np->delay;
				fprintf(fpio, frmt1, np->label, np->phase,
					np->num_pdf_residuals, res_temp);
			}
			nnodes++;
		}
	}


	return(nnodes);

}


/*** function to update hashtable */

void UpdateStaStat(int ntable, ArrivalDesc *arrival, int num_arrivals,
		double p_residual_max, double s_residual_max)
{
	int narr;
	StaStatNode *np;
	double weight = 1.0;

	for (narr = 0; narr < num_arrivals; narr++)
		if (IsPhaseID((arrival + narr)->phase, "P") &&
				fabs((arrival + narr)->residual)
							<= p_residual_max
			|| IsPhaseID((arrival + narr)->phase, "S") &&
				fabs((arrival + narr)->residual)
							<= s_residual_max)
			if ((np = InstallStaStatInTable(ntable,
					(arrival + narr)->label,
					(arrival + narr)->phase,
					(arrival + narr)->flag_ignore,
					(arrival + narr)->residual,
					weight,
					(arrival + narr)->pdf_residual_sum,
					(arrival + narr)->pdf_weight_sum,
					(arrival + narr)->delay)) == NULL)
				puterr(
"ERROR: cannot put arrival statistics in table");

}


/** end of hashtable routines */
/*------------------------------------------------------------/ */




/*------------------------------------------------------------/ */
/** 3D grid memory management routines to allow persistence of grids in memory */

#define USE_GRID_LIST 1
#define GRIDMEM_MESSAGE 3

/*** wrapper function to allocate buffer for 3D grid ***/

float* NLL_AllocateGrid(GridDesc* pgrid)
{
	int index, nactive, ngrid_read, n;
	float* fptr = NULL;
	GridMemStruct* pGridMemStruct;

//printf("IN: NLL_AllocateGrid\n");

	if (USE_GRID_LIST) {

		if ((index = GridMemList_IndexOfGridDesc(0, pgrid)) >= 0){
			// already in list
			pGridMemStruct = GridMemList_ElementAt(index);
			pGridMemStruct->active = 1;
			fptr = pGridMemStruct->buffer;
if (message_flag >= GRIDMEM_MESSAGE)
printf("GridMemManager: Grid exists in mem (%d/%d): %s %ld\n", index, GridMemListNumElements, pGridMemStruct->pgrid->title, fptr);
			return(fptr);
		} else {
			// check number of active grids in list
			nactive = 0;
			ngrid_read = 0;
			for (n = 0; n < GridMemList_NumElements(); n++) {
				pGridMemStruct = GridMemList_ElementAt(n);
				nactive += pGridMemStruct->active;
				ngrid_read += pGridMemStruct->grid_read;
			}
			// list already full of active grids, do normal allocation
			if (MaxNum3DGridMemory > 0 && nactive >= MaxNum3DGridMemory) {
				fptr = AllocateGrid(pgrid);
if (message_flag >= GRIDMEM_MESSAGE)
printf("GridMemManager: Memory full (%d/%d): %s\n", index, GridMemListNumElements, pGridMemStruct->pgrid->title);
				return(fptr);
			}
			// ok
			// remove an inactive grid if necessary
			if (MaxNum3DGridMemory > 0 && ngrid_read >= MaxNum3DGridMemory) {
				for (n = GridMemList_NumElements() - 1; n >= 0; n--) {
					pGridMemStruct = GridMemList_ElementAt(n);
					if (!pGridMemStruct->active && pGridMemStruct->grid_read) {
						GridMemList_RemoveElementAt(n);
						break;
					}
				}
			}
			// create new list element
			pGridMemStruct = GridMemList_AddGridDesc(pgrid);
			fptr = pGridMemStruct->buffer;
			if (fptr == NULL) {
				// error allocating grid memory or out of memory
				GridMemList_RemoveElementAt(GridMemList_NumElements() -1);
			}
			return(fptr);
		}

	} else {
		fptr = AllocateGrid(pgrid);
		return(fptr);
	}

}



/*** wrapper function to free buffer for 3D grid ***/

void NLL_FreeGrid(GridDesc* pgrid)
{
	int index;
	GridMemStruct* pGridMemStruct;

//printf("IN: NLL_FreeGrid\n");
	if (USE_GRID_LIST  && (index = GridMemList_IndexOfGridDesc(0, pgrid)) >= 0) {
		pGridMemStruct = GridMemList_ElementAt(index);
		pGridMemStruct->active = 0;
		pgrid->buffer = NULL;
		return;
	}

	FreeGrid(pgrid);
}



/*** wrapper function to create array for accessing 3D grid ***/

float*** NLL_CreateGridArray(GridDesc* pgrid)
{

	float*** fptr = NULL;

	int index;
	GridMemStruct* pGridMemStruct;

//printf("IN: NLL_CreateGridArray\n");
	if (USE_GRID_LIST  && (index = GridMemList_IndexOfGridDesc(0, pgrid)) >= 0) {
		pGridMemStruct = GridMemList_ElementAt(index);
		fptr = pGridMemStruct->array;
	}
	else
		fptr = CreateGridArray(pgrid);

	return(fptr);
}


/*** wrapper function to free array for accessing 3D grid ***/

void NLL_DestroyGridArray(GridDesc* pgrid)
{

	int index;

//printf("IN: NLL_DestroyGridArray\n");
	if (USE_GRID_LIST  && (index = GridMemList_IndexOfGridDesc(0, pgrid)) >= 0) {
		pgrid->array = NULL;
		return;
	}

	DestroyGridArray(pgrid);
}


/*** wrapper function to read entire grid buffer from disk ***/

int NLL_ReadGrid3dBuf(GridDesc* pgrid, FILE* fpio)
{

	int istat;

	int index;
	GridMemStruct* pGridMemStruct;

//printf("IN: NLL_ReadGrid3dBuf\n");
	if (USE_GRID_LIST  && (index = GridMemList_IndexOfGridDesc(0, pgrid)) >= 0) {
		pGridMemStruct = GridMemList_ElementAt(index);
		if (!pGridMemStruct->grid_read) {
			istat = ReadGrid3dBuf(pGridMemStruct->pgrid, fpio);
			pGridMemStruct->grid_read = 1;
		}
	} else {
		istat = ReadGrid3dBuf(pgrid, fpio);
	}

	return(0);
}


/*** add GridDescription to GridMemList ***/

GridMemStruct* GridMemList_AddGridDesc(GridDesc* pgrid)
{
	GridMemStruct* pnewGridMemStruct;

//printf("IN: GridMemList_AddGridDesc\n");
	pnewGridMemStruct = (GridMemStruct*) malloc(sizeof(GridMemStruct));
	pnewGridMemStruct->pgrid = (GridDesc*) malloc(sizeof(GridDesc));
	*(pnewGridMemStruct->pgrid) = *pgrid;
	strcpy(pnewGridMemStruct->pgrid->chr_type, pgrid->chr_type);
	strcpy(pnewGridMemStruct->pgrid->title, pgrid->title);
	pnewGridMemStruct->buffer = AllocateGrid(pnewGridMemStruct->pgrid);
	pnewGridMemStruct->array = CreateGridArray(pnewGridMemStruct->pgrid);
	pnewGridMemStruct->active = 1;
	pnewGridMemStruct->grid_read = 0;

	GridMemList_AddElement(pnewGridMemStruct);

	return(pnewGridMemStruct);

}



/*** add element to GridMemList ***/

#define LIST_SIZE_INCREMENT 10

void GridMemList_AddElement(GridMemStruct* pnewGridMemStruct)
{
	int n;
	int newGridMemListSize;
	GridMemStruct** newGridMemList;

//printf("IN: GridMemList_AddElement\n");
	if (GridMemListSize <= GridMemListNumElements) {
		// allocate enlarged list
		newGridMemListSize =  GridMemListSize + LIST_SIZE_INCREMENT;
		newGridMemList = (GridMemStruct**)
			malloc(newGridMemListSize * sizeof(GridMemStruct*));
		// load old list to new list
		for (n = 0; n < GridMemListSize; n++)
			newGridMemList[n] = GridMemList[n];
		for (n = GridMemListSize; n < newGridMemListSize; n++)
			newGridMemList[n] = NULL;
		GridMemListSize = newGridMemListSize;
		if (GridMemList != NULL)
			free(GridMemList);
		GridMemList = newGridMemList;
	}

	// load new element
	GridMemList[GridMemListNumElements] = pnewGridMemStruct;
	GridMemListNumElements++;

if (message_flag >= GRIDMEM_MESSAGE)
printf("GridMemManager: Add grid (%d): %s\n", GridMemListNumElements - 1, pnewGridMemStruct->pgrid->title);
}



/*** remove element from GridMemList ***/

void GridMemList_RemoveElementAt(int index)
{
	int n;
	GridMemStruct* pGridMemStruct;

//printf("IN: GridMemList_RemoveElementAt\n");
	// check ranges
	if (index < 0 || index >= GridMemListNumElements)
		return;

	// free allocated memory
	pGridMemStruct = GridMemList[index];
if (message_flag >= GRIDMEM_MESSAGE)
printf("GridMemManager: Remove grid (%d/%d): %s\n", index, GridMemListNumElements, pGridMemStruct->pgrid->title);
	DestroyGridArray(pGridMemStruct->pgrid);
	FreeGrid(pGridMemStruct->pgrid);
	free(pGridMemStruct->pgrid);
	free(pGridMemStruct);


	// shift down element references
	for (n = index; n < GridMemListNumElements - 1; n++)
		GridMemList[n] = GridMemList[n + 1];

	GridMemList[n] = NULL;
	GridMemListNumElements--;

	return;
}


/*** return element from GridMemList ***/

GridMemStruct* GridMemList_ElementAt(int index)
{

//printf("IN: GridMemList_ElementAt\n");
	// check ranges
	if (index < 0 || index >= GridMemListNumElements)
		return(NULL);

	return(GridMemList[index]);
}


/*** find index of grid desc in GridMemList ***/

int GridMemList_IndexOfGridDesc(int verbose, GridDesc* pgrid)
{

	int n;

//printf("IN: GridMemList_IndexOfGridDesc\n");
	for (n = 0; n < GridMemListNumElements; n++) {
if (verbose) printf("indexOf: %s ==? %s\n", GridMemList[n]->pgrid->title, pgrid->title);
		if (strcmp(GridMemList[n]->pgrid->title, pgrid->title) == 0)
			return(n);
	}

if (verbose) printf("indexOf: NOT FOUND\n");
	return(-1);

}


/*** return size of GridMemList ***/

int GridMemList_NumElements()
{

//printf("IN: GridMemList_NumElements\n");
	return(GridMemListNumElements);

}



/** end of 3D grid memory management routines */
/*------------------------------------------------------------/ */



/*** function to add station arrival to station list */

int addToStationList(SourceDesc *stations, int numStations, ArrivalDesc *arrival, int nArrivals) {

	int i, n, nAdded = 0;;


	for (i = 0; i < nArrivals; i++) {

		if (numStations >= MAX_NUM_ARRIVALS)
			return(0);

		for (n = 0; n < numStations; n++) {
			if (strcmp((stations + n)->label, (arrival + i)->label) == 0)
				break;	// already in list
		}

		if (n == numStations) {
			*(stations + n) = (arrival + i)->station;
			strcpy((stations + n)->label, (arrival + i)->label);
			nAdded++;
			numStations++;
printf("Added to station list: %s (%lf,%lf,%lf)\n", (stations + n)->label, (stations + n)->x, (stations + n)->y, (stations + n)->z);
		}

	}

	return(nAdded);


}

/*** function to write station list to file */

int WriteStationList(FILE* fpio, SourceDesc *stations, int numStations) {

	int n;

	for (n = 0; n < numStations; n++) {
		fprintf(fpio, "%s %lf %lf %lf\n",
			(stations + n)->label, (stations + n)->x, (stations + n)->y, (stations + n)->z);
	}

}




/** function to get travel times for all observed arrivals */

int getTravelTimes(ArrivalDesc *arrival, int num_arr_loc, double xval, double yval, double zval)
{
	int nReject;
	int narr, n_compan;
	FILE* fp_grid;
	double yval_grid;


	/* loop over observed arrivals */

	nReject = 0;
	for (narr = 0; narr < num_arr_loc; narr++) {
		/* check for companion */
		if ((n_compan = arrival[narr].n_companion) >= 0) {
			if ((arrival[narr].pred_travel_time =
					arrival[n_compan].pred_travel_time) < 0.0)
				nReject++;
		/* else check grid type */
		} else if (arrival[narr].gdesc.type == GRID_TIME) {
			/* 3D grid */
			if (arrival[narr].gdesc.buffer == NULL) {
				/* read time grid from disk */
				fp_grid = arrival[narr].fpgrid;
			} else {
				/* read time grid from memory buffer */
				fp_grid = NULL;
			}
			if ((arrival[narr].pred_travel_time =
				(double) ReadAbsInterpGrid3d(
					fp_grid,
					&(arrival[narr].gdesc),
					xval, yval, zval)) < 0.0)
				nReject++;
		} else {
			/* 2D grid (1D model) */
			yval_grid = GetEpiDist(
				&(arrival[narr].station),
				xval, yval);
			if ((arrival[narr].pred_travel_time =
				ReadAbsInterpGrid2d(
					NULL,
					&(arrival[narr].sheetdesc),
					yval_grid, zval)) < 0.0)
				nReject++;
		}
		arrival[narr].pred_travel_time *= arrival[narr].tfact;
	}

	return(nReject);

}


/*------------------------------------------------------------/ */
/** Octree search routines */



/** function to initialize Octree search */

Tree3D*  InitializeOcttree(GridDesc* ptgrid, ArrivalDesc* parrivals,
		int numArrLoc,  OcttreeParams* pParams)
{

	double dx, dy, dz;
	Tree3D* newTree;

	// set up octree x, y, x grid
	dx = ptgrid->dx * (double) (ptgrid->numx - 1)
		/ (double) (pParams->init_num_cells_x);
	dy = ptgrid->dy * (double) (ptgrid->numy - 1)
		/ (double) (pParams->init_num_cells_y);
	dz = ptgrid->dz * (double) (ptgrid->numz - 1)
		/ (double) (pParams->init_num_cells_z);

	newTree = newTree3D(pParams->init_num_cells_x,
			pParams->init_num_cells_y, pParams->init_num_cells_z,
		ptgrid->origx, ptgrid->origy, ptgrid->origz,
		dx, dy, dz, OCTREE_UNDEF_VALUE);

	return(newTree);
}



/*** function to perform Octree location */

#define MAX_NUM_MET_TRIES 1000

int LocOctree(int ngrid, int num_arr_total, int num_arr_loc,
		ArrivalDesc *arrival,
		GridDesc* ptgrid, GaussLocParams* gauss_par, HypoDesc* phypo,
		OcttreeParams* pParams, Tree3D* pOctTree, float* fdata,
		double* poct_node_value_max)
{

	int istat;
	int nSamples, narr, ipos;
	int nInitial;
	int writeMessage = 0;
	int iGridType;
	int nReject;
	int iReject = 0;
	double xval, yval, zval;

	long double value, dlike, dlike_max = (long double) -VERY_LARGE_DOUBLE;
	double misfit;
	double misfit_min = VERY_LARGE_DOUBLE, misfit_max = -VERY_LARGE_DOUBLE;
	double hypo_dx;

	int nScatterSaved;

	Vect3D expect = {0.0, 0.0, 0.0};

	int ix, iy, iz;
	double volume, log_value;
	ResultTreeNode* presult_node;
	OctNode* poct_node;
	OctNode* pparent_oct_node;

	int n_neigh;
	OctNode* neighbor_node;
	Vect3D coords;
	double ds_jitter;
	double smallest_node_size;



	iGridType = GRID_PROB_DENSITY;

	putmsg(3, "");
	putmsg(3, "Calculating solution in Octree...");


	/* first get solutions at each cell in Tree3D */

	nSamples = 0;
	resultTreeRoot = NULL;
	for (ix = 0; ix < pOctTree->numx; ix++) {
		for (iy = 0; iy < pOctTree->numy; iy++) {
			for (iz = 0; iz < pOctTree->numz; iz++) {
				poct_node = pOctTree->nodeArray[ix][iy][iz];
// $$$ NOTE: this block must be identical to block $$$ below
		// evaluate solution at node
		xval = poct_node->center.x;
		yval = poct_node->center.y;
		zval = poct_node->center.z;
		/* get travel times for observed arrivals */
		nReject = getTravelTimes(arrival, num_arr_loc, xval, yval, zval);
		if (nReject) {
			sprintf(MsgStr,
"ERROR: octree sample at (%lf,%lf,%lf) is outside of %d travel time grids.",
				xval, yval, zval, nReject);
			puterr(MsgStr);
		}
		/* calc misfit and prob density */
		value = CalcSolutionQuality(num_arr_loc, arrival, gauss_par, iGridType, &misfit);
		dlike = gauss_par->WtMtrxSum * exp(value);
		// calculate cell volume * prob density and put cell in resultTree
		volume = poct_node->ds.x * poct_node->ds.y * poct_node->ds.z;
		log_value =  log(gauss_par->WtMtrxSum) + value + log(volume);
		poct_node->value = log(gauss_par->WtMtrxSum) + value;
		resultTreeRoot = addResult(resultTreeRoot, log_value, poct_node);
		nSamples++;
// END - this block must be identical to block $$$ below
			}
		}
	}
	nInitial = nSamples;


	/* loop over octree nodes */

	nScatterSaved = 0;
	ipos = 0;
	smallest_node_size = VERY_LARGE_DOUBLE;

	while (nSamples < pParams->max_num_nodes) {

		presult_node = getHighestLeafValue(resultTreeRoot);

if (presult_node == NULL)
fprintf(stderr, "\npresult_node == NULL!!\n");

		pparent_oct_node = presult_node->pnode;

		// subdivide all HighestLeafValue neighbors

		for (n_neigh = 0; n_neigh < 7; n_neigh++) {

			if (n_neigh == 0) {
		  		neighbor_node = pparent_oct_node;
			} else {
				ds_jitter = octtreeParams.min_node_size / 2.0;
				coords.x = pparent_oct_node->center.x;
				coords.y = pparent_oct_node->center.y;
				coords.z = pparent_oct_node->center.z;
				if (n_neigh == 1) {
			  		coords.x = pparent_oct_node->center.x
						+ pparent_oct_node->ds.x / 2.0 + ds_jitter;
				} else if (n_neigh == 2) {
			  		coords.x = pparent_oct_node->center.x
						- pparent_oct_node->ds.x / 2.0 - ds_jitter;
				} else if (n_neigh == 3) {
			  		coords.y = pparent_oct_node->center.y
						+ pparent_oct_node->ds.y / 2.0 + ds_jitter;
				} else if (n_neigh == 4) {
			  		coords.y = pparent_oct_node->center.y
						- pparent_oct_node->ds.y / 2.0 - ds_jitter;
				} else if (n_neigh == 5) {
			  		coords.z = pparent_oct_node->center.z
						+ pparent_oct_node->ds.z / 2.0 + ds_jitter;
				} else if (n_neigh == 6) {
			  		coords.z = pparent_oct_node->center.z
						- pparent_oct_node->ds.z / 2.0 - ds_jitter;
				}
				neighbor_node = getLeafNodeContaining(pOctTree, coords);
				// outside of octTree volume
				if (neighbor_node == NULL)
					continue;
				// already subdivided
				if (neighbor_node->ds.x < 0.99 * pparent_oct_node->ds.x)
					continue;
			}


			// subdivide node and evaluate solution at each child
			subdivide(neighbor_node, OCTREE_UNDEF_VALUE);

			for (ix = 0; ix < 2; ix++) {
			  for (iy = 0; iy < 2; iy++) {
			    for (iz = 0; iz < 2; iz++) {

				poct_node = neighbor_node->child[ix][iy][iz];

//if (poct_node->ds.x < pParams->min_node_size || poct_node->ds.y < pParams->min_node_size || poct_node->ds.z < pParams->min_node_size)
//fprintf(stderr, "\nnode size too small!! %lf %lf %lf\n", poct_node->ds.x, poct_node->ds.y, poct_node->ds.z);

				// save node size if smallest so far
				if (poct_node->ds.x < smallest_node_size)
					smallest_node_size = poct_node->ds.x;
				if (poct_node->ds.y < smallest_node_size)
					smallest_node_size = poct_node->ds.y;
				if (poct_node->ds.z < smallest_node_size)
					smallest_node_size = poct_node->ds.z;

// $$$ NOTE: this block must be identical to block $$$ below
		// evaluate solution at node
		xval = poct_node->center.x;
		yval = poct_node->center.y;
		zval = poct_node->center.z;
		/* get travel times for observed arrivals */
		nReject = getTravelTimes(arrival, num_arr_loc, xval, yval, zval);
		if (nReject) {
			sprintf(MsgStr,
"ERROR: octree sample at (%lf,%lf,%lf) is outside of %d travel time grids.",
				xval, yval, zval, nReject);
			puterr(MsgStr);
		}
		/* calc misfit and prob density */
		value = CalcSolutionQuality(num_arr_loc, arrival, gauss_par, iGridType, &misfit);
		dlike = (long double) gauss_par->WtMtrxSum * (long double) exp(value);
		// calculate cell volume * prob density and put cell in resultTree
		volume = poct_node->ds.x * poct_node->ds.y * poct_node->ds.z;
		log_value =  log(gauss_par->WtMtrxSum) + value + log(volume);
		poct_node->value = log(gauss_par->WtMtrxSum) + value;
		resultTreeRoot = addResult(resultTreeRoot, log_value, poct_node);
		nSamples++;
// END - this block must be identical to block $$$ below

				/* check for maximum likelihood */
				if (dlike > dlike_max) {
				//if (misfit < misfit_min) {
//printf("MAX_LIKE: %le -> %le (%lf,%lf,%lf)\n", dlike_max, dlike, xval, yval, zval);
					dlike_max = dlike;
					misfit_min = misfit;
					phypo->misfit = misfit;
					phypo->x = xval;
					phypo->y = yval;
					phypo->z = zval;
					hypo_dx = poct_node->ds.x;
					for (narr = 0; narr < num_arr_loc; narr++)
						arrival[narr].pred_travel_time_best =
							arrival[narr].pred_travel_time;
					*poct_node_value_max = poct_node->value;
				}
				if (misfit > misfit_max)
					misfit_max = misfit;

/* set to TRUE to save all samples
   REMEMBER to set OCT num_scatter high enough in control file */
if (0) {

	/* save sample to scatter file */
	fdata[ipos++] = xval;
	fdata[ipos++] = yval;
	fdata[ipos++] = zval;
	fdata[ipos++] = dlike;

	/* update  probabilitic residuals */
	if (1)
	UpdateProbabilisticResiduals(num_arr_loc, arrival, 1.0);

	nScatterSaved++;
}



			    } // end triple loop over node children
			  }
			}

		}  // end loop over HighestLeafValue neighbors

	}



	/* give warning if ... */
/*
	if () {
		sprintf(MsgStr,
"WARNING: %d Metropolis samples rejected; travel times for an average of %.2lf arrival observations were not valid.",
			numGridReject, (double) numStaReject / numGridReject);
		putmsg(1, MsgStr);
	}
*/

	/* check reject location conditions */

	/* maximum like hypo on edge of grid */

	if (isOnGridBoundary(phypo->x, phypo->y, phypo->z,
			ptgrid, hypo_dx)) {
		sprintf(MsgStr,
"WARNING: max prob location on grid boundary, rejecting location.");
		putmsg(1, MsgStr);
		sprintf(phypo->locStatComm, "%s", MsgStr);
		iReject = 1;
	}


	/* construct search information string */
	sprintf(phypo->searchInfo, "OCTREE nInitial %d nEvaluated %d smallestNodeSide %lf\0",
			nInitial, nSamples, smallest_node_size);
	/* write message */
	putmsg(1, phypo->searchInfo);


	/* check for termination */
	if (iReject) {
		sprintf(Hypocenter.locStat, "REJECTED");
	}


	/* re-calculate solution and arrival statistics for best location */

	SaveBestLocation(num_arr_total, num_arr_loc, arrival,  ptgrid,
		gauss_par, phypo, misfit_max, iGridType);


	return(nScatterSaved);

}



/*** function to generate sample (scatter) of OctTree results */

int GenEventScatterOcttree(OcttreeParams* pParams, double oct_node_value_max, float* fscatterdata)
{

	int tot_npoints;
	int fdata_index;
	double integral;


	/* return if no scatter samples requested */
	if (pParams->num_scatter < 1)
		return(0);

	/* write message */
	putmsg(2, "");
	putmsg(2, "Generating event scatter file...");


	integral = integrateResultTree(resultTreeRoot, 0.0, oct_node_value_max);
	sprintf(MsgStr, "Octree integral= %le", integral);
	putmsg(1, MsgStr);


	/* generate scatter points at uniformly-randomly chosen locations in each leaf node */

	tot_npoints = 0;
	fdata_index = 0;
	tot_npoints = getScatterSampleResultTree(resultTreeRoot, pParams, integral,
		fscatterdata, tot_npoints, &fdata_index, oct_node_value_max);

	/* write message */
	sprintf(MsgStr, "  %d points generated, %d points requested.",
		tot_npoints, pParams->num_scatter);
	putmsg(1, MsgStr);

	return(tot_npoints);

}


/*** function to integrate exp(val) * volume of all leafs in results tree */

double integrateResultTree(ResultTreeNode* prtree, double sum, double oct_node_value_max)
{
        OctNode* pnode;

        if (prtree->left != NULL)
                sum = integrateResultTree(prtree->left, sum, oct_node_value_max);

        pnode = prtree->pnode;
        if (pnode->isLeaf)
                sum += exp(prtree->value - oct_node_value_max);

        if (prtree->right != NULL)
                sum = integrateResultTree(prtree->right, sum, oct_node_value_max);

        return(sum);
}



/*** function to get scatter sample for all leafs in results tree */

int getScatterSampleResultTree(ResultTreeNode* prtree, OcttreeParams* pParams,
		double integral, float* fdata, int npoints, int* pfdata_index,
		double oct_node_value_max)
{
	OctNode* pnode;
	double xnpoints, xval, yval, zval;
	double dx, dy, dz;

	if (prtree->right != NULL)
		npoints = getScatterSampleResultTree(
			prtree->right, pParams, integral,
			fdata, npoints, pfdata_index, oct_node_value_max);

	pnode = prtree->pnode;
	if (pnode->isLeaf /*&& npoints < pParams->num_scatter*/) {

		xnpoints = (double) pParams->num_scatter
			* exp(prtree->value - oct_node_value_max) / integral;

		xval = pnode->center.x;
		yval = pnode->center.y;
		zval = pnode->center.z;
		dx = pnode->ds.x / 2.0;
		dy = pnode->ds.y / 2.0;
		dz = pnode->ds.z / 2.0;

		while (xnpoints > 0.0 /*&& npoints < pParams->num_scatter*/) {

		    if (xnpoints > 1.0 ||
				xnpoints - (double) ((int) xnpoints) > get_rand_double(0.0, 1.0)) {
			fdata[*pfdata_index + 0] = xval + get_rand_double(-dx, dx);
//printf("npoints %d  *pfdata_index %d  %lf  dx %lf  exp(prtree->value) %le  integral %le\n", npoints, *pfdata_index, fdata[*pfdata_index + 0], dx, exp(prtree->value), integral);
			fdata[*pfdata_index + 1] = yval + get_rand_double(-dy, dy);
			fdata[*pfdata_index + 2] = zval + get_rand_double(-dz, dz);
			fdata[*pfdata_index + 3] = pnode->value;
			npoints++;
			*pfdata_index += 4;
		    }

		    xnpoints -= 1.0;

		}
	}

	if (prtree->left != NULL)
		npoints = getScatterSampleResultTree(
			prtree->left, pParams, integral,
			fdata, npoints, pfdata_index, oct_node_value_max);

	return(npoints);
}






/** end of Octree search routines */
/*------------------------------------------------------------/ */





/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

