/*
 * Copyright (C) 1999-2006 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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

	Program to do global search earthquake location in 3-D models

*/

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
	history:	(see also http://alomax.net/nlloc -> Updates)

	ver 01    26SEP1997  AJL  Original version
	ver 02    08JUN1998  AJL  Metropolis added
	ver  2         2000  AJL  Oct-Tree added
	ver  3      DEC2003  AJL  EDT added
	ver  4    10MAY2004  AJL  Added following changes from S. Husen:
		  27MAR2002  *SH   VELEST phase format added (changes in GetNextObservation)
		  28AUG2002  *SH   event origin time is now calculated relative to
		                  the second (and not minute as before); initial
				  OT seconds are now read in and added to arrival
				  time (changes in GetNextObservation)
		  17NOV2002  *SH   UUSS phase format added (changes in GetObservations
		  		  and GetNextObservation)
		  12JUN2003  *SH   Fixed bug with arrival times > 100s in UUSS phase format
		  01OCT2003  *SH   added SED format "SED_LOC"
		                  code was written by A. Lomax; bug fix by Danijel Schorlemmer
		  02MAR2004  *SH   modifications to make NLLoc compatible for routine earthquake
		                  location with SNAP (SED):
		                  - introduced 2nd argument snap_pid, which is the snap_pid of SNAP; needed
		                    to form filename of outputfile hyprint{snap_pid}; usage of NLLoc
		                     now is:
		                       NLLoc <control file> <snap_pid>
		                  - added new subroutine WriteSnapSum: output location results
		                    into file hyprint{snap_pid} in format readable by SNAP; format
		                    is identical to output format of program grid_search by
		                    M. Baer of SED;
		                    file hyprint{snap_pid} will be written if control parameter
		                    LOCHYPOUT is set to SAVE_SNAP_SUM
		                  - added subroutine get_region_names_nr and associated subroutines
		                    to convert lat/lon into Swiss coordinates in km and to find
		                    region name for local earthquakes in Switzerland
	for more history, see NLLocLib.c


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





#ifdef CUSTOM_ETH
#define PNAME  "NLLoc(ETH)"
#include "custom_eth/eth_functions.h"
#else
#define PNAME  "NLLoc"
#endif

#include "GridLib.h"
#include "ran1.h"
#include "velmod.h"
#include "NLLocLib.h"
#include "GridMemLib.h"


// function declarations
int Locate(int ngrid, char* fn_root_out, int numArrivalsReject);


/*** program to perform global search event locations */

#ifdef CUSTOM_ETH
#define NARGS_MIN 3
#define ARG_DESC "<control file> <snap_pid> <snap_param_file>"
#else
#define NARGS_MIN 2
#define ARG_DESC "<control file>"
#endif

int main(int argc, char *argv[])
{

	int istat, n;
	int i_end_of_input, iLocated;
	int narr, ngrid, nObsFile;
	int numArrivalsIgnore, numSArrivalsLocation;
	int numArrivalsReject;
	int maxArrExceeded = 0;
	int n_file_root_count = 0;
	char fn_root_out[FILENAME_MAX], fname[FILENAME_MAX], fn_root_out_last[FILENAME_MAX];
	char sys_command[MAXLINE_LONG];
	char *chr;
	FILE *fp_obs, *fpio;


	char *ppath;


	/* set program name */
	strcpy(prog_name, PNAME);

	/* check command line for correct usage */

	if (argc < NARGS_MIN) {
		disp_usage(prog_name, ARG_DESC);
		exit(EXIT_ERROR_USAGE);
	}



	// DD
	nll_mode = MODE_ABSOLUTE;

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
	ApplyCrustElevCorrFlag = 0;
	MinDistCrustElevCorr = 2.0;   // deg
	ApplyElevCorrFlag = 0;
	NumTimeDelaySurface = 0;
	iRejectDuplicateArrivals = 1;
	
	// station distance weighting 
	iSetStationDistributionWeights = 0;
	stationDistributionWeightCutoff = -1;

	// GridMemLib
	MaxNum3DGridMemory = -1;
	GridMemList = NULL;
	GridMemListSize = 0;
	GridMemListNumElements = 0;

	// GLOBAL
	NumSources = 0;
	NumStations = 0;


	/* open control file */

	strcpy(fn_control, argv[1]);
	if ((fp_control = fopen(fn_control, "r")) == NULL) {
		puterr("FATAL ERROR: opening control file.");
		exit(EXIT_ERROR_FILEIO);
	} else {
		NumFilesOpen++;
	}


#ifdef CUSTOM_ETH
	/* SH 02/27/2004
	    added snap_pid   */
	if (argc > 2)
		strcpy(snap_pid, argv[2]);
	else
		strcpy(snap_pid, "000");
	// AJL 20040527 added snap param file
	if (argc > 3)
		strcpy(snap_param_file, argv[3]);
	else
		strcpy(snap_param_file, "snap_param.txt");
#endif

	/* read NLLoc control statements from control file */

	if ((istat = ReadNLLoc_Input(fp_control)) < 0) {
		puterr("FATAL ERROR: reading control file.");
		exit(EXIT_ERROR_FILEIO);
	}
	fclose(fp_control);
	NumFilesOpen--;



	// get path to output files
	strcpy(f_outpath, fn_path_output);
	if ((ppath = strrchr(f_outpath, '/')) != NULL
			|| (ppath = strrchr(f_outpath, '\\')) != NULL)
		*(ppath + 1) = '\0';
	else
		strcpy(f_outpath, "");


	// copy control file to output directory
	strcpy(fname, fn_control);
	chr = strrchr(fname, '/');
	if (chr != NULL)
		strcpy(fname, chr +1);
	sprintf(sys_command, "cp -p %s %s_%s", fn_control, fn_path_output, fname);
	system(sys_command);
	sprintf(sys_command, "cp -p %s %slast.in", fn_control, f_outpath);
	system(sys_command);
//printf("sys_command: %s\n", sys_command);


	/* convert source location coordinates  */
	istat = ConvertSourceLoc(0, Source, NumSources, 1, 1);


	/* initialize random number generator */

	SRAND_FUNC(RandomNumSeed);
	if (message_flag >= 4)
		test_rand_int();


	/* open summary output file */

	if ((istat = OpenSummaryFiles(fn_path_output, "grid")) < 0) {
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
				&numArrivalsReject,
				MaxNumArrLoc, &Hypocenter,
				&maxArrExceeded, &numSArrivalsLocation, 0)) == 0)
			break;

		if (NumArrivals < 0)
			goto cleanup;


		/* set number of arrivals to be used in location */

		NumArrivalsLocation = NumArrivals - numArrivalsIgnore;

		putmsg(2, "");
		// AJL 20040720 SetOutName(Arrival + 0, fn_path_output, fn_root_out, fn_root_out_last, 1);
                SetOutName(Arrival + 0, fn_path_output, fn_root_out, fn_root_out_last, iSaveDecSec, n_file_root_count);
		//strcpy(fn_root_out_last, fn_root_out); /* save filename */
		sprintf(MsgStr,
"... %d observations read, %d will be used for location (%s).",
			NumArrivals + numArrivalsReject, NumArrivalsLocation, fn_root_out);
		putmsg(1, MsgStr);


		/* sort to get rejected arrivals at end of arrivals array */

		if ((istat = SortArrivalsIgnore(Arrival, NumArrivals + numArrivalsReject)) < 0) {
			puterr("ERROR: sorting arrivals by ignore flag.");
				goto cleanup;
		}


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

		if (iSetStationDistributionWeights || iSaveNLLocSum) {
printf(">>>>>>>>>>> NumStations %d, NumArrivals %d, numArrivalsReject %d\n", NumStations, NumArrivals, numArrivalsReject);
			NumStations = addToStationList(StationList, NumStations, Arrival, NumArrivals + numArrivalsReject);
			if (iSetStationDistributionWeights)
				setStationDistributionWeights(StationList, NumStations, Arrival, NumArrivals);

		}

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
			if ((istat = Locate(ngrid, fn_root_out, numArrivalsReject)) < 0) {
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
		n_file_root_count++;

		/* release grid buffer or sheet storage */

		for (narr = 0; narr < NumArrivalsLocation; narr++) {
			if (Arrival[narr].n_time_grid < 0) {	// check has opened time grid
				DestroyGridArray(&(Arrival[narr].sheetdesc));
				FreeGrid(&(Arrival[narr].sheetdesc));
				NLL_DestroyGridArray(&(Arrival[narr].gdesc));
				NLL_FreeGrid(&(Arrival[narr].gdesc));
			}
		}

		/* close time grid files (opened in function GetObservations) */

		for (narr = 0; narr < NumArrivalsLocation; narr++)
			CloseGrid3dFile(&(Arrival[narr].fpgrid), &(Arrival[narr].fphdr));

		if (iLocated) {
			putmsg(2, "");
	    		sprintf(MsgStr,
				"Finished event location, output files: %s.*",
				fn_root_out);
			putmsg(0, MsgStr);
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
	putmsg(0, MsgStr);
	putmsg(2, "");


	/* write cumulative arrival statistics */
	for (ngrid = 0; ngrid < NumLocGrids; ngrid++)
		if (LocGridSave[ngrid]) {
			sprintf(fname, "%s.sum.grid%d.loc.stat", fn_path_output, ngrid);
			if ((fpio = fopen(fname, "w")) == NULL) {
				puterr2(
"ERROR: opening cumulative phase statistics output file", fname);
				return(-1);
			} else {
				NumFilesOpen++;
			}
			WriteStaStatTable(ngrid, fpio,
				RMS_Max, NRdgs_Min, Gap_Max,
				P_ResidualMax, S_ResidualMax, Ell_Len3_Max, Hypo_Depth_Min, Hypo_Depth_Max, WRITE_RESIDUALS);
			WriteStaStatTable(ngrid, fpio,
				RMS_Max, NRdgs_Min, Gap_Max,
				P_ResidualMax, S_ResidualMax, Ell_Len3_Max, Hypo_Depth_Min, Hypo_Depth_Max, WRITE_RES_DELAYS);
			WriteStaStatTable(ngrid, fpio,
				RMS_Max, NRdgs_Min, Gap_Max,
				P_ResidualMax, S_ResidualMax, Ell_Len3_Max, Hypo_Depth_Min, Hypo_Depth_Max, WRITE_PDF_RESIDUALS);
			WriteStaStatTable(ngrid, fpio,
				RMS_Max, NRdgs_Min, Gap_Max,
				P_ResidualMax, S_ResidualMax, Ell_Len3_Max, Hypo_Depth_Min, Hypo_Depth_Max, WRITE_PDF_DELAYS);
			fclose(fpio);
			// save to last
			sprintf(sys_command, "cp %s %slast.stat", fname, f_outpath);
			system(sys_command);
			// write delays only
			sprintf(fname, "%s.sum.grid%d.loc.stat_totcorr", fn_path_output, ngrid);
			if ((fpio = fopen(fname, "w")) == NULL) {
				puterr2(
"ERROR: opening total phase corrections output file", fname);
				return(-1);
			} else {
				NumFilesOpen++;
			}
			WriteStaStatTable(ngrid, fpio,
				RMS_Max, NRdgs_Min, Gap_Max,
				P_ResidualMax, S_ResidualMax, Ell_Len3_Max, Hypo_Depth_Min, Hypo_Depth_Max, WRITE_RES_DELAYS);
			fclose(fpio);
			// save to last
			sprintf(sys_command, "cp %s %slast.stat_totcorr", fname, f_outpath);
			system(sys_command);
		}

	/* write station list */
	for (ngrid = 0; ngrid < NumLocGrids; ngrid++)
		if (LocGridSave[ngrid]) {
			sprintf(fname, "%s.sum.grid%d.loc.stations", fn_path_output, ngrid);
			if ((fpio = fopen(fname, "w")) == NULL) {
				puterr2(
"ERROR: opening cumulative phase statistics output file", fname);
				return(-1);
			} else {
				NumFilesOpen++;
			}
			WriteStationList(fpio, StationList, NumStations);
			fclose(fpio);
			// save to last
			sprintf(sys_command, "cp %s %slast.stations", fname, f_outpath);
			system(sys_command);
		}

	CloseSummaryFiles();

	exit(EXIT_NORMAL);

}



/*** function to perform grid search location */

int Locate(int ngrid, char* fn_root_out, int numArrivalsReject)
{

	int istat, n, narr;
	char fnout[FILENAME_MAX];

	FILE *fpio;
	char fname[FILENAME_MAX];
	int nScatterSaved = -1;
	float *fdata = NULL;
	float ftemp;
	int iSizeOfFdata;
	double alpha_2;
	double oct_node_value_max, oct_tree_integral = 0.0;



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
		if ((fdata = (float *) malloc(iSizeOfFdata)) == NULL)
			return(EXIT_ERROR_LOCATE);
		NumAllocations++;

	} else if (SearchType == SEARCH_OCTTREE) {
		
		// initialize memory/arrays for regular, initial oct-tree search grid
		// this is and x, y, z array of octtree root nodes, 
		// a true oct-tree is created at each of these roots
		octTree = InitializeOcttree(LocGrid + ngrid, Arrival, NumArrivalsLocation, &octtreeParams);
		NumAllocations++;

		// allocate scatter array for saved samples
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
					IsSameArrival(Arrival, narr, narr, "P")) < 0) {
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

		/* do Octree location (importance sampling) */
		if ((nScatterSaved =
			LocOctree(ngrid, NumArrivals, NumArrivalsLocation,
				Arrival, LocGrid + ngrid,
				&Gauss, &Hypocenter, &octtreeParams,
				octTree, fdata, &oct_node_value_max)) < 0) {
			puterr("ERROR: in Octree location.");
			return(EXIT_ERROR_LOCATE);
		}

	}



	/* clean up dates, caclulate rms */
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

		if (SearchType == SEARCH_OCTTREE) {
		
			// determine integral of all oct-tree leaf node pdf values
			oct_tree_integral = integrateResultTree(resultTreeRoot, 0.0, oct_node_value_max);
			sprintf(MsgStr, "Octree oct_node_value_max= %le oct_tree_integral= %le", oct_node_value_max, oct_tree_integral);
			putmsg(1, MsgStr);
		
			// generate scatter sample
			if (nScatterSaved == 0) // not saved during search
			nScatterSaved = GenEventScatterOcttree(&octtreeParams, oct_node_value_max, fdata, oct_tree_integral);
			
		}

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

		if (SearchType == SEARCH_OCTTREE && iSaveNLLocOctree) {
		
			if (LocGrid[ngrid].type == GRID_PROB_DENSITY)  {
				// convert oct tree values to likelihood
				convertOcttreeValuesToProb(resultTreeRoot, 0.0, oct_node_value_max);
				octTree->data_code = GRID_LIKELIHOOD;
				// create new result tree sorted by node values only, without multiplication by volume
//				resultTreeLikelihoodRoot = NULL;
//				resultTreeLikelihoodRoot = createResultTree(resultTreeRoot, resultTreeLikelihoodRoot);
				sprintf(MsgStr, "Oct tree structure converted to probability.");
				putmsg(1, MsgStr);
				// convert oct tree values to confidence
				//convertOcttreeValuesToConfidence(resultTreeRoot, 0.0);			
			}
				
			// write oct tree structure to file 
			sprintf(fname, "%s.loc.octree", fnout);
			if ((fpio = fopen(fname, "w")) != NULL) {
				istat = writeTree3D(fpio, octTree);
				fclose(fpio);
				sprintf(MsgStr, "Oct tree structure written to file : %d nodes", istat);
				putmsg(1, MsgStr);
			} else {
				puterr("ERROR: opening oct tree structure output file.");
				return(EXIT_ERROR_IO);
			}
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
			if ((istat = WriteGrid3dHdr(LocGrid + ngrid, NULL, fnout, "loc")) < 0) {
				puterr("ERROR: writing grid header to disk.");
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
		if ((istat = SaveLocation(&Hypocenter, ngrid, fnout, numArrivalsReject, "grid", 1)) < 0)
			return(istat);
		/* update station statistics table */
		if (strncmp(Hypocenter.locStat, "LOCATED", 7) == 0
			&& Hypocenter.rms <= RMS_Max
			&& Hypocenter.nreadings >= NRdgs_Min
			&& Hypocenter.gap <= Gap_Max
			&& Hypocenter.ellipsoid.len3 <= Ell_Len3_Max
			&& Hypocenter.z >= Hypo_Depth_Min
				  && Hypocenter.z <= Hypo_Depth_Max) {
				UpdateStaStat(ngrid, Arrival, NumArrivals, P_ResidualMax, S_ResidualMax);
//printf("INSTALLED in Stat Table: ");
				  } else {
//printf("NOT INSTALLED in Stat Table: ");
				  }
//printf("Hypo: %s %f %d %d %f %f\n", Hypocenter.locStat, Hypocenter.rms, Hypocenter.nreadings, Hypocenter.gap, Hypocenter.ellipsoid.len3, Hypocenter.z);
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
		freeTree3D(octTree, 1);
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

	return(0);


}






