/*
 * Copyright (C) 1999-2005 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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
	history:

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
		    NOV2004  AJL   Split off NLLocLib and created NLDiffLoc (non-linear double-difference location)


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





#include "GridLib.h"
#include "ran1.h"
#include "octtree.h"
#include "velmod.h"
#include "NLLocLib.h"
#include "GridMemLib.h"
#include "calc_crust_corr.h"



/** function to initialize Metropolis walk */

void InitializeMetropolisWalk(GridDesc* ptgrid, ArrivalDesc* parrivals, int
	numArrLoc, WalkParams* pMetrop, int numSamples, double initStep)
{
	int narr;
	double xlen, ylen, zlen, dminlen;
	double xmin, xmax, ymin, ymax;
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
	putmsg(4, MsgStr);

	/* set likelihood */
	pMetrop->likelihood = -1.0;

}



/** function to display and save minimum misfit location to file */

int SaveLocation(HypoDesc* hypo, int ngrid, char *fnout, int numArrivalsReject,
	char* typename, int isave_phases)
{
	int istat;
	char *pchr;
	char sys_command[MAXLINE_LONG];
	char fname[FILENAME_MAX], frootname[FILENAME_MAX];
	FILE *fp_tmp;

		/* set signature string */
		sprintf(hypo->signature, "%s   %s:v%s %s",
			LocSignature, prog_name, PVER, CurrTimeStr());
		while((pchr = strchr(hypo->signature, '\n')))
			*pchr = ' ';

		/* display hypocenter to std out */
		if (message_flag >= 3)
			WriteLocation(stdout, hypo, Arrival,
				NumArrivals + numArrivalsReject, fnout,
				isave_phases, 1, 0, LocGrid + ngrid, 0);

		/*  save requested hypocenter/phase formats */

#ifdef CUSTOM_ETH
/* SH 02/26/2004
    added new routine WriteSnapSum to write hypocenter summary to file (SNAP format) */
// must call WriteSnapSum before other output because it sets ETH magnitudes
		if (iSaveSnapSum) {
		    WriteSnapSum(NULL, hypo, Arrival, NumArrivals);
		}
#endif

		if (iSaveNLLocEvent) {
			/* write NLLoc hypocenter to event file */
			sprintf(frootname, "%s.loc", fnout);
			sprintf(fname, "%s.hyp", frootname);
			if ((istat = WriteLocation(NULL, hypo, Arrival,
					NumArrivals + numArrivalsReject, fname, isave_phases, 1, 0,
					LocGrid + ngrid, 0)) < 0) {
				puterr("ERROR: writing location to event file.");
				return(EXIT_ERROR_IO);
			}
			/* copy event file to last.hyp */
			sprintf(sys_command, "cp %s %slast.hyp", fname, f_outpath);
			system(sys_command);
			sprintf(fname, "%s.hdr", frootname);
			sprintf(sys_command, "cp %s %slast.hdr", fname, f_outpath);
			system(sys_command);
			sprintf(fname, "%s.scat", frootname);
			if ((fp_tmp = fopen(fname, "r")) != NULL) {
				fclose(fp_tmp);
				sprintf(sys_command, "cp %s %slast.scat", fname, f_outpath);
				system(sys_command);
			}
		}
		if (iSaveNLLocSum) {
			/* write NLLoc hypocenter to summary file */
			if ((istat = WriteLocation(pSumFileHypNLLoc[ngrid],
					hypo,
					Arrival, NumArrivals, fnout, 0, 1, 0,
					LocGrid + ngrid, 0)) < 0) {
				puterr("ERROR: writing location to summary file.");
				return(EXIT_ERROR_IO);
			}
			/* copy event grid header to .sum header */
			sprintf(sys_command,
				"cp %s.loc.hdr %s.sum.%s%d.loc.hdr",
				fnout, fn_path_output, typename, ngrid);
			system(sys_command);
		}


		if (iSaveHypo71Event) {
			/* write HYPO71 hypocenter to event file */
			WriteHypo71(NULL, hypo, Arrival, fnout, 1, 1);
		}
		if (iSaveHypo71Sum) {
			/* write HYPO71 hypocenter to summary file */
			WriteHypo71(pSumFileHypo71[ngrid], hypo,
				Arrival, fnout, iWriteHypHeader[ngrid], 0);
		}

		if (iSaveHypoEllEvent) {
			/* write pseudo-HypoEllipse hypo to event file */
			WriteHypoEll(NULL, hypo, Arrival, fnout, 1, 1);
		}
		if (iSaveHypoEllSum) {
			/* write pseudo-HypoEllipse hypo to summary file */
			WriteHypoEll(pSumFileHypoEll[ngrid], hypo,
				Arrival,
				fnout, iWriteHypHeader[ngrid], 0);
		}

		if (iSaveHypoInvSum) {
			/* write HypoInverseArchive (FPFIT) hypocenter to summary file */
			iWriteHypHeader[ngrid] = 0;
			WriteHypoInverseArchive(pSumFileHypoInv[ngrid], hypo,
				Arrival, NumArrivals, fnout, iWriteHypHeader[ngrid], 1);
			/* also write to last.hypo_inv */
			sprintf(fname, "%slast.hypo_inv", f_outpath);
//printf(">>%s\n", fname);
			if ((fp_tmp = fopen(fname, "w")) != NULL) {
				WriteHypoInverseArchive(fp_tmp, hypo,
					Arrival, NumArrivals, fnout, iWriteHypHeader[ngrid], 1);
				fclose(fp_tmp);
			}
		}

		if (iSaveAlberto4Sum) {
			/* write Alberto 4 SIMULPS format */
			WriteHypoAlberto4(pSumFileAlberto4[ngrid], hypo,
				Arrival, fnout);
		}

		iWriteHypHeader[ngrid] = 0;

	return(0);
}




/*** function to check and initialize an observation */

int checkObs(ArrivalDesc *arrival, int nobs)
{

	int istat;

	/* check for aliased arrival label */
	if ((istat = EvaluateArrivalAlias(arrival + nobs)) < 0)
		;
	/* set some fields */
	InitializeArrivalFields(arrival + nobs);
	/* check some fields */
	if (!isgraph(arrival[nobs].phase[0]))
		strcpy(arrival[nobs].phase, ARRIVAL_NULL_STR);
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
	putmsg(3, MsgStr);
	/* remove blanks/whitespace from phase string */
	removeSpace(arrival[nobs].phase);
	/* check for reject phase code */
	if (IsPhaseID(arrival[nobs].phase, "$")) {
		sprintf(MsgStr,
"WARNING: phase code is $, rejecting observation: %s %s", arrival[nobs].label, arrival[nobs].phase);
		putmsg(2, MsgStr);
		return(-1);
	}
	/* check for valid P or S phase code */
/*INGV		if (!IsPhaseID(arrival[nobs].phase, "P")
				&& !IsPhaseID(arrival[nobs].phase, "S")) {
			sprintf(MsgStr,
"WARNING: phase code not in P or S phase id list, rejecting observation: %s %s", arrival[nobs].label, arrival[nobs].phase);
			putmsg(2, MsgStr);
			return(-1);
		}
INGV*/
	/* check for duplicate arrival */
	if (nll_mode != MODE_DIFFERENTIAL) {
		// AJL 20041201 - alogoritm changed, less strict
		//if ( IsDuplicateArrival(arrival, nobs + 1, nobs, NULL) >= 0
		//		|| IsDuplicateArrival(arrival, nobs + 1, nobs, arrival[nobs].phase) >= 0 )
		//		) {
		if ( IsDuplicateArrival(arrival, nobs + 1, nobs) >= 0 ) {
			sprintf(MsgStr,
	"WARNING: duplicate arrival, rejecting observation: %s %s", arrival[nobs].label, arrival[nobs].phase);
			putmsg(2, MsgStr);
			return(-1);
		}
	}

	// all OK
	return(1);

	}



/*** function read observation file and open station grid files */

int GetObservations(FILE* fp_obs, char* ftype_obs, char* fn_grids,
	ArrivalDesc *arrival, int *pi_end_of_input,
	int* pnignore, int* pnreject, int maxNumArrivals, HypoDesc* phypo,
	int* pMaxArrExceeded, int *pnumSArrivals, int nobs_prev)
{

	int nobs, nobs_read, nobs_total;
	int istat, ntry, nLocate, n_NoAbsTime, n_compan, n_time_grid;
	char filename[FILENAME_MAX];
	char eval_phase[PHASE_LABEL_LEN];
	int isZeroWeight;
	char arrival_phase[PHASE_LABEL_LEN];
	int read_2d_sheets;
	int i_need_elev_corr, n_phs_try;

	SourceDesc* pstation;

	*pnumSArrivals = 0;


	/** read observations to arrival array */

	nobs = nobs_prev;
	*pnignore = 0;
	*pnreject = 0;
	nLocate = 0;
	ntry = 0;
	while ((istat = GetNextObs(fp_obs, arrival + nobs, ftype_obs, ntry++ == 0)) != EOF) {


		if (istat == OBS_FILE_FORMAT_ERROR) {
			putmsg(1, "ERROR: format error in observation file.");
			break;
		}
		if (istat == OBS_FILE_END_OF_EVENT) {
			if (nobs == 0)
				continue;
			else
				break;
		}
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
		/* SH if comment line found -> read next obs */
		if (istat == OBS_IS_COMMENT_LINE)
		        continue;
		if (nobs == MAX_NUM_ARRIVALS - 1) {
			if (!*pMaxArrExceeded) {
				*pMaxArrExceeded = 1;
				putmsg(1, "WARNING: maximum number of arrivals exceeded.");
			}
			continue;
		}

		// check and initialize observation
		if (checkObs(arrival, nobs) > 0)
			nobs++;
		// check and initialize second observation
		// SH if istat = OBS_FILE_TWO_ARRIVALS_READ then S-arrival has been read in -> two observations per line
		if (istat == OBS_FILE_TWO_ARRIVALS_READ)
			if (checkObs(arrival, nobs + 1) > 0)
				nobs++;

	}


	// save number of obs read
	nobs_total = nobs;
	nobs_read = nobs - nobs_prev;


	/** check for minimum number of arrivals */

//printf("*pi_end_of_input %d  nobs_read %d  MinNumArrLoc %d\n", *pi_end_of_input, nobs_read, MinNumArrLoc);

	if (!(*pi_end_of_input) && (nobs_total > 0 && nobs_total < MinNumArrLoc)) {
		sprintf(MsgStr,
"WARNING: too few observations to locate, skipping event.");
		putmsg(1, MsgStr);
		sprintf(MsgStr,
"INFO: %d observations needed (specified in control file entry LOCMETH).",
		MinNumArrLoc);
		putmsg(3, MsgStr);
		return(-1);
	}



	/** process arrivals, determine if reject or ignore */

	// check for no absolute timing (inst begins with '*')
	n_NoAbsTime = CheckAbsoluteTiming(arrival + nobs_prev, nobs_read);

	/* homogenize date/time */
	if ((istat = HomogDateTime(arrival, nobs_total, phypo)) < 0) {
		puterr("ERROR: in arrival date/times.");
		if (istat == OBS_FILE_ARRIVALS_CROSS_YEAR_BOUNDARY)
			puterr("ERROR: arrivals cross year boundary.");
		return(-1);
	}


	/* sort to get arrivals in time order */
	if (nll_mode != MODE_DIFFERENTIAL && (istat = SortArrivalsTime(arrival, nobs_total)) < 0) {
		puterr("ERROR: sorting arrivals by time.");
		return(-1);
	}


	/* check each arrival and initialize */

	Num3DGridReadToMemory = 0;
	for (nobs = nobs_prev; nobs < nobs_total; nobs++) {

		if (message_flag > 2) {
			sprintf(MsgStr, "Checking Arrival %d:  %s (%s)  %s %s %s %d",
				nobs,
				arrival[nobs].label,
				arrival[nobs].time_grid_label,
				arrival[nobs].onset,
				arrival[nobs].phase,
				arrival[nobs].first_mot,
				arrival[nobs].quality);
			putmsg(3, MsgStr);
		} else if (nobs_read > 1000 && message_flag > 0 && nobs % 100 == 0) {
			fprintf(stdout, "Checking Arrival %d/%d\r", nobs, nobs_read);
			fflush(stdout);
		}

		// make sure station location is null
		arrival[nobs].station.x = arrival[nobs].station.y = arrival[nobs].station.z = 0.0;


		// set some flags
		read_2d_sheets = 1;
		i_need_elev_corr = 0;

		// check for zero weight phase
		strcpy(arrival_phase, arrival[nobs].phase);
		isZeroWeight = 0;
		if (arrival_phase[0] == '*') {
			isZeroWeight = 1;
			strcpy(arrival_phase, arrival_phase + 1);
		}

		strcpy(arrival[nobs].fileroot, "\0");
		arrival[nobs].n_companion = -1;
		arrival[nobs].n_time_grid = -1;
		arrival[nobs].tfact = 1.0;

		/* if Vp/Vs > 0, check for companion phase, its time grid will be used for times */
		if (VpVsRatio > 0.0) {
			/* try finding previously initialized companion phase */
			if (IsPhaseID(arrival_phase, "S") &&
					(n_compan =
					IsSameArrival(arrival, nobs, nobs, "P")) >= 0 &&
					arrival[n_compan].flag_ignore == 0) {
				arrival[nobs].tfact = VpVsRatio;
				arrival[nobs].gdesc.type = arrival[n_compan].gdesc.type;
				arrival[nobs].station = arrival[n_compan].station;
				arrival[nobs].n_companion = n_compan;
				sprintf(MsgStr,
"INFO: S phase: %d %s %s using companion phase %d travel time grids.",
				nobs, arrival[nobs].label, arrival[nobs].phase, arrival[nobs].n_companion);
				putmsg(3, MsgStr);
				// save filename as grid identifier (needed for GridMemList)
				strcpy(arrival[nobs].gdesc.title, arrival[n_compan].gdesc.title);
			}
		}

		/* if MODE_DIFFERENTIAL, check for same station phase, its time grid will be used for times */
		if (nll_mode == MODE_DIFFERENTIAL && arrival[nobs].n_companion < 0) {
			/* try finding previously initialized companion phase */
			if ((n_compan = IsSameArrival(arrival, nobs, nobs, NULL)) >= 0 &&
					arrival[n_compan].flag_ignore == 0) {
				arrival[nobs].gdesc.type = arrival[n_compan].gdesc.type;
				arrival[nobs].station = arrival[n_compan].station;
				arrival[nobs].n_companion = n_compan;
				sprintf(MsgStr,
"INFO: MODE_DIFFERENTIAL: %d %s %s using companion phase %d travel time grids.",
				nobs, arrival[nobs].label, arrival[nobs].phase, arrival[nobs].n_companion);
				putmsg(3, MsgStr);
				// save filename as grid identifier (needed for GridMemList)
				strcpy(arrival[nobs].gdesc.title, arrival[n_compan].gdesc.title);
			}
		}

		/* no previously initialized companion phase */
		if (arrival[nobs].n_companion < 0) {

			// try to open time grid file using original phase ID
			sprintf(arrival[nobs].fileroot, "%s.%s.%s", fn_grids,
				arrival_phase, arrival[nobs].time_grid_label);
			sprintf(filename, "%s.time", arrival[nobs].fileroot);
			// try opening time grid file for this phase
			istat = OpenGrid3dFile(filename,
					&(arrival[nobs].fpgrid),
					&(arrival[nobs].fphdr),
					&(arrival[nobs].gdesc), "time",
					&(arrival[nobs].station),
					arrival[nobs].gdesc.iSwapBytes);

			if (istat < 0 ) {
				// try to open time grid file using LOCPHASEID mapped phase ID
				EvalPhaseID(eval_phase, arrival_phase);
				sprintf(arrival[nobs].fileroot, "%s.%s.%s", fn_grids,
					eval_phase, arrival[nobs].time_grid_label);
				sprintf(filename, "%s.time", arrival[nobs].fileroot);
				/* try opening time grid file for this phase */
				istat = OpenGrid3dFile(filename,
						&(arrival[nobs].fpgrid),
						&(arrival[nobs].fphdr),
						&(arrival[nobs].gdesc), "time",
						&(arrival[nobs].station),
						arrival[nobs].gdesc.iSwapBytes);
			}

			/* try opening P time grid file for S if no P companion phase */
			if (istat < 0 && VpVsRatio > 0.0 && IsPhaseID(arrival_phase, "S")) {
				arrival[nobs].tfact = VpVsRatio;
				sprintf(arrival[nobs].fileroot, "%s.%s.%s", fn_grids,
						"P", arrival[nobs].time_grid_label);
				sprintf(filename, "%s.time", arrival[nobs].fileroot);
				istat = OpenGrid3dFile(filename,
					&(arrival[nobs].fpgrid),
					&(arrival[nobs].fphdr),
					&(arrival[nobs].gdesc), "time",
					&(arrival[nobs].station),
					arrival[nobs].gdesc.iSwapBytes);
				sprintf(MsgStr,
"INFO: S phase: using P phase travel time grid file: %s", filename);
					putmsg(3, MsgStr);
			}

			// check if station/source in grid hdr file was DEFAULT
			if (istat >= 0 && strcmp(arrival[nobs].station.label, "DEFAULT") == 0) {
				// get station/source coordinates, etc
				pstation = FindSource(arrival[nobs].time_grid_label);
				if (pstation != NULL) {
					arrival[nobs].station = *pstation;
					i_need_elev_corr = 1;
				} else
					istat = -1;
			}

			// try opening DEFAULT time grid file
			if (istat < 0) {
				pstation = FindSource(arrival[nobs].time_grid_label);
				if (pstation != NULL) {
					// open DEFAULT time grid for this phase
					n_phs_try = 0;
					while (istat < 0 && n_phs_try < 2) {
						n_phs_try++;
						if (n_phs_try == 1)	// try to open time grid file using original phase ID
							strcpy(eval_phase, arrival_phase);
						else	// try to open time grid file using LOCPHASEID mapped phase ID
							EvalPhaseID(eval_phase, arrival_phase);
						sprintf(arrival[nobs].fileroot, "%s.%s.%s", fn_grids, eval_phase, "DEFAULT");
						sprintf(filename, "%s.time", arrival[nobs].fileroot);

						/* check if time grid already read */
						if ((n_time_grid = FindDuplicateTimeGrid(arrival, nobs, nobs)) >= 0 &&
								arrival[n_time_grid].flag_ignore == 0) {
							arrival[nobs].gdesc.type = arrival[n_time_grid].gdesc.type;
							arrival[nobs].sheetdesc = arrival[n_time_grid].sheetdesc;
							arrival[nobs].station = *pstation;
							arrival[nobs].n_time_grid = n_time_grid;
							sprintf(MsgStr,
			"INFO: DEFAULT travel time: %d %s %s using previous phase %d travel time grids.",
							nobs, arrival[nobs].label, arrival[nobs].phase, arrival[nobs].n_time_grid);
							putmsg(3, MsgStr);
							// save filename as grid identifier (needed for GridMemList)
							strcpy(arrival[nobs].gdesc.title, arrival[n_time_grid].gdesc.title);
						istat = 1;
						} else {

							istat = OpenGrid3dFile(filename,
								&(arrival[nobs].fpgrid),
								&(arrival[nobs].fphdr),
								&(arrival[nobs].gdesc), "time",
								&(arrival[nobs].station),
								arrival[nobs].gdesc.iSwapBytes);
							sprintf(MsgStr,
			"INFO: using DEFAULT travel time grid file: %s", filename);
								putmsg(3, MsgStr);
							arrival[nobs].station = *pstation;
						}
					}
					i_need_elev_corr = 1;
					//read_2d_sheets = 0; // too slow!
				}
			}

			// save filename as grid identifier (needed for GridMemList)
			strcpy(arrival[nobs].gdesc.title, filename);

			if (istat < 0) {
				sprintf(MsgStr,
"WARNING: cannot open time grid file: %s: rejecting observation: %s %s",
filename, arrival[nobs].label, arrival[nobs].phase);
				putmsg(2, MsgStr);
				CloseGrid3dFile(&(arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
				strcpy(arrival[nobs].fileroot, "\0");
				goto RejectArrival;
			}


			/* check that search grid is inside time grid (3D grids) */

			if (arrival[nobs].gdesc.type == GRID_TIME &&
				!IsGridInside(LocGrid + 0, &(arrival[nobs].gdesc), 0)) {
				sprintf(MsgStr,
"WARNING: initial location search grid not contained inside arrival time grid, rejecting observation: %s %s", arrival[nobs].label, arrival[nobs].phase);
				putmsg(1, MsgStr);
				CloseGrid3dFile(&(arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
				goto RejectArrival;
			}


			/* (2D grids) check that
				(1) greatest distance from search grid to station
					is inside time grid, and
				(2) distance from center of search grid to station
					is greater than DistStaGridMin (if DistStaGridMin > 0)
					is less than DistStaGridMax (if DistStaGridMax > 0) */

			if (arrival[nobs].gdesc.type == GRID_TIME_2D &&
					(istat = IsGrid2DBigEnough(LocGrid + 0,
					&(arrival[nobs].gdesc),
					&(arrival[nobs].station),
					DistStaGridMin, DistStaGridMax,
					LocGrid[0].origx + (LocGrid[0].dx
						* (double) (LocGrid[0].numx - 1)) / 2.0,
					LocGrid[0].origy + (LocGrid[0].dy
						* (double) (LocGrid[0].numy - 1)) / 2.0)
				) != 1)
			{
				CloseGrid3dFile(&(arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
				if (istat == -1) {
					sprintf(MsgStr,
"WARNING: greatest distance from initial 3D location search grid to station \n\texceeds 2D time grid size, rejecting observation: %s %s",
						arrival[nobs].label, arrival[nobs].phase);
					putmsg(2, MsgStr);
					goto RejectArrival;
				} else if (istat == -2) {
					sprintf(MsgStr,
"WARNING: distance from grid center to station \n\texceeds maximum station distance, ignoring observation in misfit calculation: %s %s",
						arrival[nobs].label, arrival[nobs].phase);
					putmsg(2, MsgStr);
					arrival[nobs].flag_ignore = 1;
					goto IgnoreArrival;
				}
				if (istat == -3) {
					sprintf(MsgStr,
"WARNING: depth range of initial 3D location search grid exceeds that of 2D time grid size, rejecting observation: %s %s",
						arrival[nobs].label, arrival[nobs].phase);
					putmsg(2, MsgStr);
					goto RejectArrival;
				}
			}

		} /* if (arrival[nobs].n_companion < 0) */


		/* check for time delays */
		ApplyTimeDelays(arrival + nobs);
			;
		// calculate elevation correction if needed
		if (ApplyElevCorrFlag && i_need_elev_corr) {
			arrival[nobs].elev_corr = CalcSimpleElevCorr(arrival, nobs, ElevCorrVelP, ElevCorrVelS);
		}



		/* check if arrival is excluded */

		// zero weight or excluded in control file
		if (isZeroWeight || isExcluded(arrival[nobs].label, arrival[nobs].phase)) {
			sprintf(MsgStr,
"INFO: arrival is excluded, ignoring observation in misfit calculation: %s %s",
				arrival[nobs].label, arrival[nobs].phase);
			arrival[nobs].flag_ignore = 1;
			CloseGrid3dFile(&(arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
			putmsg(2, MsgStr);
			goto IgnoreArrival;
		}
		// no absolute time and not EDT
		if (!arrival[nobs].abs_time && LocMethod != METH_EDT) {
			sprintf(MsgStr,
"INFO: arrival does not have absolute timing, ignoring observation in misfit calculation: %s %s %s",
				arrival[nobs].label, arrival[nobs].phase, arrival[nobs].inst);
			arrival[nobs].flag_ignore = 1;
			CloseGrid3dFile(&(arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
			putmsg(2, MsgStr);
			goto IgnoreArrival;
		}
		// EDT_BOX but not BOX error
		if (!strcmp(arrival[nobs].error_type, "BOX") && LocMethod != METH_EDT_BOX) {
			sprintf(MsgStr,
"INFO: method is EDT_BOX but arrival does not have error type BOX, ignoring observation in misfit calculation: %s %s %s",
				arrival[nobs].label, arrival[nobs].phase, arrival[nobs].inst);
			arrival[nobs].flag_ignore = 1;
			CloseGrid3dFile(&(arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
			putmsg(2, MsgStr);
			goto IgnoreArrival;
		}


		/* check if maximum number of arrivals for location exceeded */

		if (nLocate + nobs_prev >= maxNumArrivals) {
			putmsg(2,
"WARNING: maximum number of arrivals for location exceeded, \n\tignoring observation in misfit calculation.");
			arrival[nobs].flag_ignore = 1;
			CloseGrid3dFile(&(arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
			goto IgnoreArrival;
		}



		/** arrival to be used for location */

		/* no further processing for arrivals with companion time grids */
		if (arrival[nobs].n_companion >= 0 || arrival[nobs].n_time_grid >= 0)
			goto AcceptArrival;

		/** prepare time grids access in memory or on disk */

		/* construct dual-sheet description
			(2D grids for all search types,
				and 3D grids for grid-search) */

		if (SearchType == SEARCH_GRID
				|| (read_2d_sheets && arrival[nobs].gdesc.type == GRID_TIME_2D)) {

			arrival[nobs].sheetdesc = arrival[nobs].gdesc;
	//INGV ??
	//if (arrival[nobs].gdesc.numx > 1)
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
				CloseGrid3dFile(&(arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
				Num3DGridReadToMemory++;
			}
		}


		/* read time grid and close file (2D grids)*/

		if (read_2d_sheets && arrival[nobs].gdesc.type == GRID_TIME_2D) {
			istat = ReadArrivalSheets(1, &(arrival[nobs]), 0.0);
			CloseGrid3dFile(&(arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
			if (istat  < 0) {
				sprintf(MsgStr,
"ERROR: reading arrival travel time sheets (2D grid), rejecting observation: %s %s",
					arrival[nobs].label, arrival[nobs].phase);
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

		(*pnreject)++;
		arrival[nobs].flag_ignore = 999;
		sprintf(MsgStr, "   Rejected Arrival %d:  %s (%s)  %s %s %s %d",
			nobs,
			arrival[nobs].label,
			arrival[nobs].time_grid_label,
			arrival[nobs].onset,
			arrival[nobs].phase,
			arrival[nobs].first_mot,
			arrival[nobs].quality);
		putmsg(3, MsgStr);

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

 	arrival->abs_time = 1;
 	arrival->obs_centered = 0.0;
 	arrival->pred_travel_time = 0.0;
 	arrival->pred_travel_time_best = -LARGE_DOUBLE;
	arrival->pred_centered = 0.0;
 	arrival->cent_resid = 0.0;
 	arrival->obs_travel_time = 0.0;
// 	arrival->delay = 0.0;	// do not initialize here!
 	arrival->elev_corr = 0.0;
 	arrival->residual = 0.0;
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
	arrival->gdesc.iSwapBytes = iSwapBytesOnInput;
	arrival->sheetdesc.buffer = NULL;
	
	arrival->station_weight = 1.0;
	

//DD
	if (nll_mode != MODE_DIFFERENTIAL) {
		arrival->weight = 0.0;
	}


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



/*** function to determine type P or S of last leg of phase */

// returns 'P' if P, 'S' if S, ' ' otherwise

char lastLegType(ArrivalDesc *arrival)
{
	char *c_last_p, *c_last_P, *c_last_s, *c_last_S;
	int i_last_p, i_last_P, i_last_s, i_last_S;


	// get offsets to last char of types p/P or s/S
	c_last_p = strrchr(arrival->phase, 'p');
	i_last_p = c_last_p == NULL ? -1 : (int) (c_last_p - arrival->phase);
	c_last_P = strrchr(arrival->phase, 'P');
	i_last_P = c_last_P == NULL ? -1 : (int) (c_last_P - arrival->phase);
	c_last_s = strrchr(arrival->phase, 's');
	i_last_s = c_last_s == NULL ? -1 : (int) (c_last_s - arrival->phase);
	c_last_S = strrchr(arrival->phase, 'S');
	i_last_S = c_last_S == NULL ? -1 : (int) (c_last_S - arrival->phase);

	// determine type of last char
	i_last_p = i_last_p > i_last_P ? i_last_p : i_last_P;
	i_last_s = i_last_s > i_last_S ? i_last_s : i_last_S;
	if (i_last_p >= 0 && i_last_p > i_last_s)
		return('P');
	if (i_last_s >= 0 && i_last_s > i_last_p)
		return('S');
	return(' ');

}



/*** function to calculuate simple elevation correction based on surface velcoity */

double CalcSimpleElevCorr(ArrivalDesc *arrival, int narr, double pvel, double svel)
{
	double yval_grid, t_surface, t_elev, elev_corr;
	int n_compan, diagnostic;


	// Switch debug messages from this function on/off (1/0).
	diagnostic = message_flag >= 3;
diagnostic = 1;

	// check for companion
	if ((n_compan = arrival[narr].n_companion) >= 0) {
		if (diagnostic) {
			sprintf(MsgStr,"CalcSimpleElevCorr: n_compan=%d", n_compan);
			putmsg(1, MsgStr);
		}
		if ((elev_corr = arrival[n_compan].elev_corr) < 0.0)
			return(0.0);
	// check for specific P or S vel
	} else if (pvel > 0.0 && lastLegType(arrival + narr) == 'P') {
		elev_corr = -arrival[narr].station.depth / pvel;
	} else if (svel > 0.0 && lastLegType(arrival + narr) == 'S') {
		elev_corr = -arrival[narr].station.depth / svel;
	// else check grid type
	} else {
		if (arrival[narr].gdesc.type == GRID_TIME) {
			if (diagnostic) {
				sprintf(MsgStr,"CalcSimpleElevCorr: GRID_TIME");
				putmsg(1, MsgStr);
			}
			/* 3D grid */
			if ((t_surface = (double) ReadAbsInterpGrid3d(arrival[narr].fpgrid,
					&(arrival[narr].gdesc), 0.0, 0.0, 0.0)) < 0.0)
				return(0.0);
			// read pos along X dir because grid may not extend above surface
			if ((t_elev = (double) ReadAbsInterpGrid3d(arrival[narr].fpgrid,
					&(arrival[narr].gdesc), fabs(arrival[narr].station.depth), 0.0, 0.0)) < 0.0)
				return(0.0);
		} else {
			if (diagnostic) {
				sprintf(MsgStr,"CalcSimpleElevCorr: GRID_TIME_2D");
				putmsg(1, MsgStr);
			}
			/* 2D grid (1D model) */
			if ((t_surface = ReadAbsInterpGrid2d(arrival[narr].fpgrid,
					&(arrival[narr].gdesc), 0.0, 0.0)) < 0.0)
				return(0.0);
			// read pos along Horiz dir because grid may not extend above surface
			yval_grid = fabs(arrival[narr].station.depth);
			if (GeometryMode == MODE_GLOBAL)
				yval_grid *= KM2DEG;
			if ((t_elev = ReadAbsInterpGrid2d(arrival[narr].fpgrid,
					&(arrival[narr].gdesc), yval_grid, 0.0)) < 0.0)
				return(0.0);
		}
		if (arrival[narr].station.depth > 0.0)	// below surface
			t_elev *= -1.0;
		elev_corr = t_elev - t_surface;
	}

	elev_corr *= arrival[narr].tfact;

	if (diagnostic) {
		sprintf(MsgStr,"CalcSimpleElevCorr: lat=%.3f  lon=%.3f  depth=%.3f  elev_corr=%.3f",
			arrival[narr].station.dlat, arrival[narr].station.dlong, arrival[narr].station.depth, elev_corr);
		putmsg(1, MsgStr);
	}

	return (elev_corr);

}



/*** function to evaluate arrival label aliases */

int EvaluateArrivalAlias(ArrivalDesc *arrival)
{
	int nAlias;
	int checkAgain = 1, icount = 0, aliasApplied = 0;
	char *pchr, tmpLabel[MAXLINE];


	strcpy(tmpLabel, arrival->label);

	sprintf(MsgStr, "Checking for station name alias: %s", tmpLabel);
	putmsg(4, MsgStr);


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
			putmsg(4, MsgStr);

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
		putmsg(4, "");
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
		putmsg(4, "");
		return(0);
	}


	putmsg(4, "");


	return(0);
}



/*** function to apply time delays */

int ApplyTimeDelays(ArrivalDesc *arrival)
{
	int nDelay, i;
	int ifound = 0;
	double tmp_delay;


	sprintf(MsgStr, "Checking for time delay: %s %s",
		arrival->label, arrival->phase);
	putmsg(4, MsgStr);


	/* check time delays for match to label/phase */

	for (nDelay = 0; nDelay < NumTimeDelays; nDelay++) {

		if (strcmp(TimeDelay[nDelay].label, arrival->label) == 0
			&& strcmp(TimeDelay[nDelay].phase, arrival->phase) == 0)
		{
		    /* apply station delays */

		    tmp_delay = TimeDelay[nDelay].delay;
		    arrival->delay = 0.0;
		    if (fabs(tmp_delay) > VERY_SMALL_DOUBLE) {
// DELAY_CORR			arrival->sec -= tmp_delay;
			arrival->delay = tmp_delay;
		arrival->obs_time -= (long double) arrival->delay;   // DELAY_CORR	- incorporating delay so subtract (Tcorr = Tobs - (O-C))
			sprintf(MsgStr,
			    "   delay of %lf sec subtracted from obs time.",
				tmp_delay);
			putmsg(4, MsgStr);
			ifound = 1;
		    }
		    break;
		}
	}
	putmsg(4, "");


	// if time delay not found, check for delay surface
	if (!ifound && NumTimeDelaySurface) {
		tmp_delay = LARGE_FLOAT;
		for (i = 0; i < NumTimeDelaySurface; i++) {
			if (strcmp(arrival->phase, TimeDelaySurfacePhase[i]) == 0) {
				tmp_delay = ApplySurfaceTimeDelay(i, arrival);
				tmp_delay *= TimeDelaySurfaceMultiplier[i];
				break;
			}
		}
		if (i < NumTimeDelaySurface && tmp_delay < LARGE_FLOAT / 2.0) {
			arrival->delay = tmp_delay;
		arrival->obs_time -= (long double) arrival->delay;   // DELAY_CORR	- incorporating delay so subtract (Tcorr = Tobs - (O-C))
printf("%s %s %s, ", arrival->label, arrival->phase, TimeDelaySurfacePhase[i]);
			sprintf(MsgStr,
			"    %s surface delay of %lf sec at lat %f, long %f subtracted from obs time.",
				TimeDelaySurfacePhase[i], tmp_delay,
				arrival->station.dlat, arrival->station.dlong);
			putmsg(4, MsgStr);
putmsg(0, MsgStr);
		}
	}


	return(0);
}




/*** function to get time delay from time delay surface */

double ApplySurfaceTimeDelay(int nsurface, ArrivalDesc *arrival)
{
	double x, y;

	if (arrival->station.is_coord_latlon) {
		x = arrival->station.dlong;
		y = arrival->station.dlat;
	} else {
		return(LARGE_FLOAT);
	}

	return(get_surface_z(nsurface, x, y));
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
	int istat = 0, iloop;
	char *cstat;
	char chr, instruction[MAXSTRING];
	char chrtmp[MAXSTRING];

	double psec, ssec;

	double weight, ttime;

	int ioff;

	static char line[MAXLINE_LONG];

	static int date_saved, year_save, month_save, day_save;
	static int check_for_S_arrival;

	static int in_hypocenter_event;

	int ifound;
	
	// ETH LOC format
	char eth_line_key[MAXSTRING];
	int eth_use_loc, eth_use_mag;
	double  vpvs;

	// NEIC format
	static char last_label[10];
	char cmonth[4];
	char* pchr;

	// ISC format
	char isc_time_str[10];

	// DD
	static int hypo_cc_flag;
	static long int dd_event_id_1, dd_event_id_2;
	static double dd_otime_corr;
	double tt_sta1, tt_sta2;



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
printf("1 %s", line);
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
printf("2 %s", line);

				iloop = 1;
				continue;
			} else {		// skip this line
				iloop = 1;
				continue;
			}
		}
		if (istat != 1) {
			puterr2("WARNING: unrecognized control statement", line);
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
	strcpy(arrival->error_type, "GAU");
 	arrival->error = 0.0;
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
	

/* SH Seismograph Station University of Utah format (PING)
	first line: reftime (yymmddhhmm)
	next lines arrival times rel. to reftime incl. duration
	*/

	else if (strcmp(ftype_obs,"UUSS") == 0)
	{
		/* interpret line: reftime or phase information or comment line etc */
		chr = fgetc(fp_obs);
		/* fprintf(stderr,"GetNextObs: chr %d\n",chr); */
		if (chr == EOF) {
			return(OBS_FILE_END_OF_INPUT);
		} else if (chr == 99 || chr == 67) {
		  /* comment line starts with 'c' */
		     /* read until end of line */
		      ungetc(chr,fp_obs);
		      cstat = fgets(line, MAXLINE_LONG, fp_obs);
		      return(OBS_IS_COMMENT_LINE);
		} else if (chr == 35) {
		  /* end of event marked by '#' */
		     /* read until end of line */
		      ungetc(chr,fp_obs);
		      cstat = fgets(line, MAXLINE_LONG, fp_obs);
		  return(OBS_FILE_END_OF_EVENT);
		} else if (chr == 32) {
			/* line contains reference time starting with a space */
			ungetc(chr,fp_obs);
			cstat = fgets(line, MAXLINE_LONG, fp_obs);
			istat = ReadFortranInt(line,2,2,&EventTime.year);
			if (EventTime.year < 20)
				EventTime.year += 100;
			EventTime.year += 1900;
			istat += ReadFortranInt(line,4,2,&EventTime.month);
			istat += ReadFortranInt(line,6,2,&EventTime.day);
			istat += ReadFortranInt(line,8,2,&EventTime.hour);
			istat += ReadFortranInt(line,10,2,&EventTime.min);
			arrival->sec = 0.0;
		} else {
			/* line contains phase information */
			ungetc(chr,fp_obs);
		}

		/* read next line */
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);
		/* now read phase info */
		istat = ReadFortranString(line,1,4,arrival->label);
		/* check if comment line */
		if ( (strncmp(arrival->label,"c",1) == 0) || (strncmp(arrival->label,"C",1) == 0) ) {
			return(OBS_IS_COMMENT_LINE);
		}
		TrimString(arrival->label);
		strcpy(arrival->phase,"P");
		if ( (istat += ReadFortranString(line,6,1,arrival->first_mot)) != 2) {
		/* read next lines until arrival time is found */
		  while ((cstat = fgets(line, MAXLINE_LONG, fp_obs)) != NULL) {
		    if (strncmp(line,"#",1) == 0 )
		      /* end of event detected */
		      return (OBS_FILE_END_OF_EVENT);
		    istat = 0;
		    /* fprintf(stderr,"GetNextObs: %s\n",line); */
		    /* fprintf(stderr,"GetNextObs: %d %s\n",istat,arrival->first_mot); */
		    if ( (istat += ReadFortranString(line,6,1,arrival->first_mot)) == 1 ) {
			istat = ReadFortranString(line,1,4,arrival->label);
			if ( (strncmp(arrival->label,"c",1) == 0) || (strncmp(arrival->label,"C",1) == 0) ){
			  return(OBS_IS_COMMENT_LINE);
			  }
			TrimString(arrival->label);
			break;
		    }
		   }
		}
		
		if ( cstat == NULL )
		     return (OBS_FILE_END_OF_INPUT);

		/* continue with reading of remaining phase information */
		istat += ReadFortranInt(line,10,1,&arrival->quality);
		istat += ReadFortranReal(line,12,6,&arrival->sec);

		arrival->year = EventTime.year;
		arrival->month = EventTime.month;
		arrival->day = EventTime.day;
		arrival->hour = EventTime.hour;
		arrival->min = EventTime.min;

		/* convert quality to error */
		Qual2Err(arrival);

		/* fprintf(stderr,"GetNextObs: phase %4s %1s %1s %1d %5.2f\n",arrival->label,
			arrival->phase,arrival->first_mot,arrival->quality,arrival->sec); */

		/* check for possible S-arrival */
		istat += ReadFortranReal(line, 21, 6, &ssec);
		if ( ssec > 0.0) {
			arrival++;
			istat = ReadFortranString(line,1,4,arrival->label);
			TrimString(arrival->label);
			strcpy(arrival->phase,"S");
			istat += ReadFortranInt(line,19,1,&arrival->quality);
			arrival->year = EventTime.year;
			arrival->month = EventTime.month;
			arrival->day = EventTime.day;
			arrival->hour = EventTime.hour;
			arrival->min = EventTime.min;
			arrival->sec = ssec;
			/* convert quality to error */
			Qual2Err(arrival);
			return(OBS_FILE_TWO_ARRIVALS_READ);
		} else
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


		check_for_S_arrival = 1;
		/* check for valid phase code */
//		if (IsPhaseID(arrival->phase, "P")) {
			//strcpy(arrival->phase, "P");
			check_for_S_arrival = 1;

//		} else if (IsPhaseID(arrival->phase, "S")) {
			//strcpy(arrival->phase, "S");
//			check_for_S_arrival = 0;
//		} else
//			return(OBS_FILE_END_OF_EVENT);

		if (!IsGoodDate(arrival->year, arrival->month, arrival->day))
			return(OBS_FILE_END_OF_EVENT);


		/* check for integer sec format */
		if (arrival->sec > 99.999)
			arrival->sec /= 100.0;

		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}



	else if (strcmp(ftype_obs, "INGV_BOLL") == 0 ||
		strcmp(ftype_obs, "INGV_BOLL_LOCAL") == 0 ||
		strcmp(ftype_obs, "INGV_ARCH") == 0)
	{
/*
   1 0
DOI 1016   EZPG  03243563 EZSG  243882                             23
  11 0
SEI 1016   EZPG  04132382 EZSG  132881
VMG 1016   EZPG  04132514 EZSG  133068                             56
ZCCA1016   EZPG  04132655                                          47
*/
		if (check_for_S_arrival > 0) {
			/* check for additional phase input in last input line read */

			ioff = 0;
			if (strcmp(ftype_obs, "INGV_ARCH") == 0)
				ioff = 13 * (check_for_S_arrival - 1);
			/* check for zero or blank S phase time */
			istat = ReadFortranString(line, 33 + ioff, 6, chrtmp);
			if (istat > 0
				&& strncmp(chrtmp, "      ", 6) != 0
				&& strncmp(chrtmp, "     0", 6) != 0)
			{
				/* read S phase input in last input line read */
				istat = ReadFortranString(line, 1, 4, arrival->label);
				TrimString(arrival->label);
				//istat += ReadFortranInt(line, 10, 2, &arrival->year);
			// no year in INGV boll!
			arrival->year = 2003;
				istat += ReadFortranInt(line, 5, 2, &arrival->month);
				istat += ReadFortranInt(line, 7, 2, &arrival->day);
				istat += ReadFortranInt(line, 18, 2, &arrival->hour);
				istat += ReadFortranInt(line, 33 + ioff, 2, &arrival->min);
				istat += ReadFortranReal(line, 35 + ioff, 4, &arrival->sec);
				arrival->sec /= 100.0;
				istat += ReadFortranReal(line, 22, 4, &psec);
				psec /= 100.0;
		// check for P second >= 60.0 (correct for INGV???)
		if (psec >= 60.0 && arrival->sec < 60.0)
			arrival->sec += 60.0;
				istat += ReadFortranString(line, 27 + ioff, 1, arrival->onset);
				// set weight
				arrival->quality = arrival->onset[0] == 'i' ? 1 : 2;
				istat += ReadFortranString(line, 28 + ioff, 1, arrival->comp);
				if (arrival->comp[0] == ' ')
					arrival->comp[0] = ARRIVAL_NULL_CHR;
				arrival->phase[0] = '\0';
				istat += ReadFortranString(line, 29 + ioff, 4, arrival->phase);
				TrimString(arrival->phase);
			}
		}

		/* check for S arrival input found */
		if (check_for_S_arrival > 0 && istat == 10
//				&& IsPhaseID(arrival->phase, "S")
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

			// increment additional arrival count
			if (strcmp(ftype_obs, "INGV_ARCH") == 0 && check_for_S_arrival < 3) {
				check_for_S_arrival++;
			} else {
				check_for_S_arrival = 0;
				line[0] = '\0';
			}
			return(istat);
		} else if (check_for_S_arrival > 0) {
			check_for_S_arrival = 0;
			return(OBS_FILE_SKIP_INPUT_LINE);
		}



/*
   1 0
DOI 1016   EZPG  03243563 EZSG  243882                             23
  11 0
SEI 1016   EZPG  04132382 EZSG  132881
VMG 1016   EZPG  04132514 EZSG  133068                            156
ZCCA1016   EZPG  04132655                                          47
*/
		// read next line
		cstat = fgets(line, MAXLINE_LONG, fp_obs);

		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);

		/* read formatted P (or S) arrival input */
		istat = ReadFortranString(line, 1, 4, arrival->label);
		if (isspace(arrival->label[0])) {
			// check to skip next event(s) if is teleseism
			if (strcmp(ftype_obs, "INGV_BOLL_LOCAL") == 0) {
				while(1) {
					TrimString(line);
					if (strlen(line) > 10) {
						while(1) {		// read to next event
							cstat = fgets(line, MAXLINE_LONG, fp_obs);
							if (cstat == NULL)
								return(OBS_FILE_END_OF_INPUT);
							if(isspace(line[0]))
								break;
						}
					} else {
						break;
					}
				}
			}
			line[0] = '\0';
			return(OBS_FILE_END_OF_EVENT);
		}
		TrimString(arrival->label);
		istat += ReadFortranString(line, 12, 1, arrival->onset);
		TrimString(arrival->onset);
		istat += ReadFortranString(line, 13, 1, arrival->comp);
		if (arrival->comp[0] == ' ')
			arrival->comp[0] = ARRIVAL_NULL_CHR;
		// set weight
		arrival->quality = arrival->onset[0] == 'i' ? 0 : 1;
		istat += ReadFortranString(line, 14, 4, arrival->phase);
		TrimString(arrival->phase);
		istat += ReadFortranString(line, 10, 1, arrival->first_mot);
		TrimString(arrival->first_mot);
	// no year in INGV boll!
	arrival->year = 2003;
		istat += ReadFortranInt(line, 5, 2, &arrival->month);
		istat += ReadFortranInt(line, 7, 2, &arrival->day);
		istat += ReadFortranInt(line, 18, 2, &arrival->hour);
		istat += ReadFortranInt(line, 20, 2, &arrival->min);
		istat += ReadFortranReal(line, 22, 4, &arrival->sec);
		arrival->sec /= 100.0;


		if (istat != 10) {
			line[0] = '\0';
			return(OBS_FILE_END_OF_EVENT);
		}

		/* read optional amplitude/period fields */
		istat += ReadFortranReal(line, 66, 4, &arrival->amplitude);
//printf("AMPLITUDE: %lf (%s)\n", arrival->amplitude, line);
		//istat += ReadFortranReal(line, 48, 3, &arrival->period);
// 20030826 AJL coda_dur added
		//istat += ReadFortranReal(line, 71, 5, &arrival->coda_dur);


		check_for_S_arrival = 1;
		/* check for valid phase code */
//		if (IsPhaseID(arrival->phase, "P")) {
			//strcpy(arrival->phase, "P");
//			check_for_S_arrival = 1;

//		} else if (IsPhaseID(arrival->phase, "S")) {
			//strcpy(arrival->phase, "S");
//			check_for_S_arrival = 0;
//		} else
//			return(OBS_FILE_END_OF_EVENT);

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
		if (strcmp(line, "    ") == 0) {
			line[0] = '\0';
			return(OBS_FILE_END_OF_EVENT);
		}

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

		if (!IsGoodDate(arrival->year, arrival->month, arrival->day))
			return(OBS_FILE_END_OF_EVENT);

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
			return(OBS_FILE_SKIP_INPUT_LINE);

		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}
	
	

	else if (strcmp(ftype_obs, "PAOLO_OV") == 0 ||
		strcmp(ftype_obs, "ALBERTO_3D") == 0 ||
		strcmp(ftype_obs, "ALBERTO_3D_4") == 0 ||
		strcmp(ftype_obs, "SIMULPS") == 0 ||
		strcmp(ftype_obs, "VELEST") == 0 )
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

"VELEST"
 phase format of VELEST program; from Stephan Husen email: stephan@tomo.ig.erdw.ethz.ch
summary line containing event origin time, location, and magnitude
 format:a4,a2,1x,a2,i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,f7.2,f7.2
        iyrmo, iday, ihr, mino, seco, ltde, cns, lnde, cew, dep, mag
observed travel times (six sets per line)
 each set is in format a4,a1,a1,f6.2
                       sta, phaseID, weight, tt

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

		  } else if (strcmp(ftype_obs, "VELEST") == 0) {
			/* try to read phase arrival input */
			istat = fscanf(fp_obs, "%12c", line);
			line[12]= '\0';
			if (istat == EOF)
			return(OBS_FILE_END_OF_INPUT);
			/* printf("<%s>\n", line); */
			istat += ReadFortranString(line, 1, 4, arrival->label);
			istat += ReadFortranString(line, 5, 1, arrival->phase);
			istat += ReadFortranInt(line, 6, 1, &arrival->quality);
			istat += ReadFortranReal(line, 7, 6, &arrival->sec);

			if (istat != 5)
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



	else if (strcmp(ftype_obs, "SED_LOC") == 0)
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
	          /* SH 03/05/2003
	          read control parameters, which are specified in first line of SED_LOC format
	          so far only VpVSRatio is of use
	          the value of VpVsRatio specified in NLLoc control file will be overwritten!  */
	       		if (strcmp(eth_line_key, "INST") == 0) {
	       		     istat = ReadFortranReal(line, 23, 6, &vpvs);
	       		     if (istat != 1) {
	       		       puterr(
	       		       "WARNING: could not read VpVsRatio! Use value of control file\n");
	       		     } else {
	       		       VpVsRatio = vpvs;
	       		     }
	       		}
				if (strcmp(eth_line_key, "Loca") == 0) {
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
/* SH 03/05/2004  added event number  */
			istat += ReadFortranInt(line, 71, 8, &EventTime.ev_nr);
			if (istat == 6) {
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
//		istat += ReadFortranReal(line, 21, 6, &arrival->sec);
		istat += ReadFortranInt(line, 28, 1, &eth_use_loc);
		istat += ReadFortranReal(line, 29, 5, &arrival->period);
		istat += ReadFortranReal(line, 34, 9, &arrival->amplitude);
/* SH 02/25/2004 bug fix
     arrival->inst is only of type char[5] and not char[9]!
		istat += ReadFortranString(line, 44,9, arrival->inst);  */
		istat += ReadFortranString(line, 44,4, arrival->inst);
	/* SH 07232004 added */
		//istat += ReadFortranInt(line, 54, 1, &eth_use_mag);
		istat += ReadFortranInt(line, 54, 1, &arrival->clipped);


		if (istat != 10) {
			return(OBS_FILE_END_OF_EVENT);
		}

		/* remove blanks/whitespace */
		removeSpace(arrival->label);
		removeSpace(arrival->phase);
		removeSpace(arrival->inst);
		
/* AJL EvalPhaseID not used here anymore, used later */
/*		if (EvalPhaseID(arrival->phase) < 0) {
			puterr2("WARNING: phase ID not found", arrival->phase);
			return(OBS_FILE_INVALID_PHASE);
		}
*/

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
			arrival->amplitude = AMPLITUDE_NULL;

		return(istat);
	}



//  ETH_LOC replaced by SED_LOC
	//else if (strcmp(ftype_obs, "ETH_LOC") == 0)
	else if (0)
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
		//20040226 S Husen bug fix  istat += ReadFortranString(line, 44,9, arrival->inst);
		istat += ReadFortranString(line, 44, 4, arrival->inst);
		istat += ReadFortranInt(line, 54, 1, &eth_use_mag);


		if (istat != 10) {
			return(OBS_FILE_END_OF_EVENT);
		}

		/* remove blanks/whitespace */
		removeSpace(arrival->label);
		removeSpace(arrival->phase);
		removeSpace(arrival->inst);
		
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
			arrival->amplitude = AMPLITUDE_NULL;

		return(istat);
	}



	else if (strcmp(ftype_obs, "NEIC") == 0)
	{

/* example:
     20 JUN 2003  (171)

     ot  = 06:19:38.51   +/-   0.17              AMAZONAS, BRAZIL
     lat =      -7.513   +/-    3.5
     lon =     -71.642   +/-    4.6              MAGNITUDE 7.1 (HRV)
     dep =       555.2  (depth phases)

     110 km (70 miles) E of Cruzeiro do Sul, Brazil (pop N/A)
     335 km (205 miles) ENE of Pucallpa, Peru
     410 km (255 miles) SSW of Leticia, Colombia
     2730 km (1700 miles) WNW of BRASILIA, Brazil

     nph =  194 of 401    se = 0.83        FE=113                      A

     error ellipse = ( 65.1,155.1,  0.0;  0.0,  0.0,  0.0;  6.3,  4.0,  0.0)

 mb = 6.4 ( 93)  ML = 0.0 (  0)  mblg = 0.0 (  0)  md = 0.0 (  0)  MS = 0.0 (  0)

 sta  phase     arrival     res   dist azm    amp  per mag     amp  per mag  sta
 SAML iPc     06:21:44.12   2.0    8.5 100                                   SAML
      eS      06:23:25.14  103.X
 LPAZ eP      06:21:51.01  -0.7    9.4 159 L:3.7+0 .93 6.2X g:3.0+0 .93 5.4X LPAZ
 OTAV ePc     06:22:02.21   1.7   10.3 318 g:9.9+0 1.0 6.0X                  OTAV
      eS      06:23:47.88  107.X
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
			//     20 JUN 2003  (171)
			istat = ReadFortranInt(line, 13, 4, &EventTime.year);
			istat += ReadFortranString(line, 9, 3, cmonth);
			istat += ReadFortranInt(line, 6, 2, &EventTime.day);
			if (istat == 3) {
				in_hypocenter_event = 1;
				EventTime.month = Month2Int(cmonth);
				// find phs line
				while (strncmp(line, " sta", 4) != 0) {
					// read next line
					cstat = fgets(line, MAXLINE_LONG, fp_obs);
					if (cstat == NULL)
						return(OBS_FILE_END_OF_INPUT);
				}
				/* read next line */
				cstat = fgets(line, MAXLINE_LONG, fp_obs);
				if (cstat == NULL)
					return(OBS_FILE_END_OF_INPUT);
			} else {
				return(OBS_FILE_END_OF_EVENT);
			}
		}


		/* read phase arrival input */
/*
 SAML iPc     06:21:44.12   2.0    8.5 100                                   SAML
      eS      06:23:25.14  103.X
 LPAZ eP      06:21:51.01  -0.7    9.4 159 L:3.7+0 .93 6.2X g:3.0+0 .93 5.4X LPAZ
*/
		istat = ReadFortranString(line, 2, 4, arrival->label);
		if (strncmp(arrival->label, "    ", 4) == 0)
			strcpy(arrival->label, last_label);
		else
			strcpy(last_label, arrival->label);
		TrimString(arrival->label);
		TrimString(last_label);
		istat += ReadFortranString(line, 7, 1, arrival->onset);
		istat += ReadFortranString(line, 8, 6, arrival->phase);
		TrimString(arrival->phase);
		// check for valid onset
		pchr = arrival->onset;
		// onset 'x' added by AJL to allow zero weighting of phases
		if (!(*pchr == ' ' || *pchr == 'i' || *pchr == 'e' || *pchr == 'q' || *pchr == 'x')) {
			sprintf(chrtmp, "%c%s", *pchr, arrival->phase);
			strcpy(arrival->phase, chrtmp);
			*pchr = ARRIVAL_NULL_CHR;
		}
		// check for first motion
		pchr = &(arrival->phase[strlen(arrival->phase) - 1]);
		if (*pchr == 'c' || *pchr == 'd') {
			arrival->first_mot[0] = *pchr;
			*pchr = '\0';
		}
		// set weight, check for P - reset weight
		if (arrival->onset[0] == 'i')
			arrival->quality = 0;
		else if (arrival->onset[0] == 'e')
			arrival->quality = 1;
		else if (arrival->onset[0] == 'x')
			arrival->quality = 4;
		else 
			arrival->quality = 2;
		if (strncmp(arrival->phase, "P", 1) != 0) {
			arrival->quality += 1;
		}
		istat += ReadFortranInt(line, 15, 2, &arrival->hour);
		istat += ReadFortranInt(line, 18, 2, &arrival->min);
		istat += ReadFortranReal(line, 21, 5, &arrival->sec);

		if (istat != 6) {
			return(OBS_FILE_END_OF_EVENT);
		}

		arrival->year = EventTime.year;
		arrival->month = EventTime.month;
		arrival->day = EventTime.day;


		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}



	else if (strcmp(ftype_obs, "ISC_IMS1.0") == 0)
	{

/* example:
...
DATA_TYPE BULLETIN IMS1.0:short
ISC Comprehensive Bulletin
Event  2957613 Sicily

   Date       Time        Err   RMS Latitude Longitude  Smaj  Smin  Az Depth   Err Ndef Nsta Gap  mdist  Mdist Qual   Author      OrigID
2002/04/05 04:52:15.70         0.60  38.7500   14.2800                  14.0               6                          BJI        4180045
2002/04/05 04:52:20.23   0.05        37.9440   16.4600  38.9  17.4  -1  10.0   1.1         9 348                      PDG        4681048
...
Sta     Dist  EvAz Phase        Time      TRes  Azim AzRes   Slow   SRes Def   SNR       Amp   Per Qual Magnitude    ArrID
MSI     0.45 121.8 Pg       04:52:30.8     1.0                           T__                        _e            71903810
SCLL    0.55 110.1 Pg       04:52:32.5     1.0                           T__                        _e            71903811
MTTG    0.67 131.6 Pg       04:52:34.5     0.6                           T__                        _e            71903812
...
SGO     2.12   5.0 PN       04:52:59.0     2.6                           T__                        _e            71903828
LVI     2.19 258.8 PN       04:52:58.3     0.8                           T__                        _e            71903826
LVI     2.19 258.8                                                       ___           736.0  0.60  __            71903827
MRLC    2.33   8.0 PN       04:53:01.3     1.9                           T__                        _e            71903830
CSSN    2.41 359.2 PN       04:53:01.0     0.4                           T__                        _e            71903829
...
AVF    11.99 317.8 P        04:55:12.9   -1.37                           T__            20.8  1.26  _e            71869701
LOR    12.04 320.6 P        04:55:13.3   -1.73                           T__            22.9  1.39  _e            71869702
LOR    12.04 320.6 R                                                     ___           126.3 18.75  _e            71869703
SSF    12.08 319.1 P        04:55:15.3    -0.3                           T__            12.1  1.04  _e            71869705
BGF    12.11 315.9 P        04:55:14.6   -1.37                           T__            34.9  1.38  _e            71869704
MOX    12.44 349.8 P        04:55:31.3   10.96                           T__                        _e            75622714
MOX    12.44 349.8 LR                                                    ___           700.0 16.00  __        3.8 75622715
MOX    12.44 349.8 LR                                                    ___           300.0 16.00  __        3.8 75622717
MOX    12.44 349.8 LR                                                    ___           600.0 15.00  __            75622716
BRG    12.45 356.7 P        04:55:21.5     1.0                           T__             9.7  1.20  __            75622621
BRG    12.45 356.7          04:55:40.5                                   ___                        _i            75622622
BRG    12.45 356.7 LR                                                    ___           870.0 13.80  __        3.9 75622623
BRG    12.45 356.7 LR                                                    ___            60.0 13.80  __        3.8 75622625
BRG    12.45 356.7 LR                                                    ___           690.0 13.80  __            75622624
ETSF   12.66 295.5 P        04:55:22.4    -1.0                           T__            31.9  1.41  _e            71869706
CLL    12.94 354.2 P        04:55:26     -1.05                           T__                        __            79343809
...


STOP

*/
		/* read line */
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
//printf("1 %s", line);
		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);
		/* check for end of event (assumes Event or STOP at end of event) */
		if ((in_hypocenter_event && strncmp(line, "Event", 5) == 0)
				|| strncmp(line, "STOP", 4) == 0) {
			/* end of event */
			return(OBS_FILE_END_OF_EVENT);
		}


		/* not yet in event, find and read event hypocenter line */
		if (!in_hypocenter_event) {
			// assume may alredy be in "Event" block
			// find date
			while (strncmp(line, "   Date", 7) != 0) {
				// read next line
				cstat = fgets(line, MAXLINE_LONG, fp_obs);
//printf("2 %s", line);
				if (cstat == NULL)
					return(OBS_FILE_END_OF_INPUT);
			}
			// read next line
			cstat = fgets(line, MAXLINE_LONG, fp_obs);
//printf("3 %s", line);
			if (cstat == NULL)
				return(OBS_FILE_END_OF_INPUT);
			/* read hypocenter date */
//2002/04/05 04:52:15.70         0.60  38.7500   14.2800                  14.0               6                          BJI        4180045
			istat = ReadFortranInt(line, 1, 4, &EventTime.year);
			istat += ReadFortranInt(line, 6, 2, &EventTime.month);
			istat += ReadFortranInt(line, 9, 2, &EventTime.day);
			if (istat == 3) {
				in_hypocenter_event = 1;
				// find phases lines
				while (strncmp(line, "Sta     Dist", 12) != 0) {
					// read next line
					cstat = fgets(line, MAXLINE_LONG, fp_obs);
//printf("4 %s", line);
					if (cstat == NULL)
						return(OBS_FILE_END_OF_INPUT);
				}
				/* read next line */
				cstat = fgets(line, MAXLINE_LONG, fp_obs);
//printf("5 %s", line);
				if (cstat == NULL)
					return(OBS_FILE_END_OF_INPUT);
			} else {
				return(OBS_FILE_END_OF_EVENT);	// blank line or ignored line
			}
		}


		/* read phase arrival input */
/*
LVI     2.19 258.8 PN       04:52:58.3     0.8                           T__                        _e            71903826
LVI     2.19 258.8                                                       ___           736.0  0.60  __            71903827
*/

		// check if is a reading with a seconds time
		//while (ReadFortranString(line, 44, 3, isc_time_str) != 1
		while (ReadFortranString(line, 35, 3, isc_time_str) != 1
				|| strncmp(isc_time_str, "   ", 3) == 0) {
			// not a time, check if end
			if (cstat == NULL || strncmp(line, "Event", 5) == 0
					|| strncmp(line, "STOP", 4) == 0)
				return(OBS_FILE_END_OF_EVENT);
			// not a time, not end, read next line
			cstat = fgets(line, MAXLINE_LONG, fp_obs);
//printf("6 %s", line);
			if (cstat == NULL)
				return(OBS_FILE_END_OF_INPUT);
		}

		// read phase time reading
		istat = ReadFortranString(line, 1, 5, arrival->label);
		if (strncmp(arrival->label, "     ", 5) == 0)
			strcpy(arrival->label, last_label);
		else
			strcpy(last_label, arrival->label);
		TrimString(arrival->label);
		TrimString(last_label);
		istat += ReadFortranString(line, 20, 8, arrival->phase);
		TrimString(arrival->phase);
		istat += ReadFortranInt(line, 29, 2, &arrival->hour);
		istat += ReadFortranInt(line, 32, 2, &arrival->min);
		istat += ReadFortranReal(line, 35, 6, &arrival->sec);
		istat += ReadFortranString(line, 101, 1, arrival->first_mot);
		istat += ReadFortranString(line, 102, 1, arrival->onset);
		// set weight, check for P - reset weight
		if (arrival->onset[0] == 'i')
			arrival->quality = 0;
		else if (arrival->onset[0] == 'e')
			arrival->quality = 1;
		else if (arrival->onset[0] == 'q')
			arrival->quality = 2;
		else if (arrival->onset[0] == '_')
			arrival->quality = 2;
		else
			arrival->quality = 3;
		// check if "P" arrival
		if (strncmp(arrival->phase, "P", 1) != 0) {
			arrival->quality += 1;
		}

		if (istat != 7) {
			return(OBS_FILE_INVALID_PHASE);
		}

		arrival->year = EventTime.year;
		arrival->month = EventTime.month;
		arrival->day = EventTime.day;


		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}



	else if (strcmp(ftype_obs, "CSEM_GSE2.0") == 0)
	{

/* example:
...
BEGIN GSE2.0
MSG_TYPE DATA
MSG_ID 041203222858 EMSC
DATA_TYPE BULLETIN GSE2.0

EVENT 20041203222858
   Date       Time       Latitude Longitude    Depth    Ndef Nsta Gap    Mag1  N    Mag2  N    Mag3  N  Author          ID
       rms   OT_Error      Smajor Sminor Az        Err   mdist  Mdist     Err        Err        Err     Quality

2004/12/03 22:28:58.5     44.2742    7.5419      8.0      66   53 109  ML 3.3 11                        MIX       00000001
      1.22   +-  0.10       2.6    1.6  123               0.42   7.13   +-0.6                           a i ke

 ( Data from stations operated by :     LDG,   MDD,   ZUR )

NORTHERN ITALY             
Sta     Dist  EvAz     Phase      Date       Time     TRes  Azim  AzRes  Slow  SRes Def   SNR       Amp   Per   Mag1   Mag2       ID
SBF     0.42 190.6 a   Pn      2004/12/03 22:29:09.3   1.3                          T     1.0                                    LDG
MBDF    0.71 309.7 a   Pg      2004/12/03 22:29:09.8  -2.4                          T     1.0                                    LDG
FRF     0.96 222.5 a   Pg      2004/12/03 22:29:16.4  -0.2                          T     1.0                                    LDG
FRF     0.96 222.5 a   Sg      2004/12/03 22:29:27.9  -1.3                          T     1.0     134.3  0.16 ML 4.0             LDG
...
BAIF    6.22 339.7 a   Pn      2004/12/03 22:30:30.2   1.9                          T     1.0                                    LDG
ERTA    6.27 240.5 m   P       2004/12/03 22:30:31.9  -1.8                          T               0.8  0.16                    MAD
LDF     6.83 311.9 a   Pn      2004/12/03 22:30:38.8   2.1                          T     1.0                                    LDG
GRR     7.12 308.2 a   Pn      2004/12/03 22:30:41.8   1.1                          T     1.0                                    LDG
FLN     7.13 311.9 a   Pn      2004/12/03 22:30:41.8   1.1                          T     1.0                                    LDG

STOP

...

*/
		/* read line */
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
//printf("1 %s", line);
		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);
		/* check for end of event (assumes Event or STOP at end of event) */
		if ((in_hypocenter_event && strncmp(line, "EVENT", 5) == 0)
				|| strncmp(line, "STOP", 4) == 0) {
			/* end of event */
			return(OBS_FILE_END_OF_EVENT);
		}


		/* not yet in event, find and read event hypocenter line */
		if (!in_hypocenter_event) {
			// assume may alredy be in "Event" block
			// find phases lines
			while (strncmp(line, "Sta     Dist", 12) != 0) {
				// read next line
				cstat = fgets(line, MAXLINE_LONG, fp_obs);
//printf("4 %s", line);
				if (cstat == NULL)
					return(OBS_FILE_END_OF_INPUT);
			}
			/* read next line */
			cstat = fgets(line, MAXLINE_LONG, fp_obs);
//printf("5 %s", line);
			if (cstat == NULL)
				return(OBS_FILE_END_OF_INPUT);
			in_hypocenter_event = 1;
		}


		/* read phase arrival input */



/* this block needed for ISC_IMSL because of amplitude readings, etc.
		// check if is a reading with a seconds time
		//while (ReadFortranString(line, 44, 3, isc_time_str) != 1
		while (ReadFortranString(line, 35, 3, isc_time_str) != 1
				|| strncmp(isc_time_str, "   ", 3) == 0) {
			// not a time, check if end
			if (cstat == NULL || strncmp(line, "Event", 5) == 0
					|| strncmp(line, "STOP", 4) == 0)
				return(OBS_FILE_END_OF_EVENT);
			// not a time, not end, read next line
			cstat = fgets(line, MAXLINE_LONG, fp_obs);
//printf("6 %s", line);
			if (cstat == NULL)
				return(OBS_FILE_END_OF_INPUT);
		}
*/

		// read phase time reading
		istat = ReadFortranString(line, 1, 5, arrival->label);
		if (strncmp(arrival->label, "     ", 5) == 0)
			strcpy(arrival->label, last_label);
		else
			strcpy(last_label, arrival->label);
		TrimString(arrival->label);
		TrimString(last_label);
		istat += ReadFortranString(line, 24, 8, arrival->phase);
		TrimString(arrival->phase);
		istat += ReadFortranInt(line, 32, 4, &arrival->year);
		istat += ReadFortranInt(line, 37, 2, &arrival->month);
		istat += ReadFortranInt(line, 40, 2, &arrival->day);
		istat += ReadFortranInt(line, 43, 2, &arrival->hour);
		istat += ReadFortranInt(line, 46, 2, &arrival->min);
		istat += ReadFortranReal(line, 49, 5, &arrival->sec);
		istat += ReadFortranString(line, 21, 1, arrival->first_mot);
		istat += ReadFortranString(line, 22, 1, arrival->onset);
		// set weight, check for P - reset weight
		if (arrival->onset[0] == 'i')
			arrival->quality = 0;
		else if (arrival->onset[0] == ' ')
			arrival->quality = 0;
		else if (arrival->onset[0] == 'e')
			arrival->quality = 1;
		else if (arrival->onset[0] == 'q')
			arrival->quality = 2;
		else if (arrival->onset[0] == '_')
			arrival->quality = 2;
		else
			arrival->quality = 3;
		// check if "P" arrival
		if (strncmp(arrival->phase, "P", 1) != 0) {
			arrival->quality += 1;
		}

		if (istat != 10) {
			return(OBS_FILE_INVALID_PHASE);
		}

		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}



	else if (strcmp(ftype_obs, "CSEM_ALERT") == 0)
	{

/* example:

...
Sta       Phase     Date     Time       Res     Dist Azm    Net
---------------------------------------------------------------
PRNI  m   Pn    2004/02/11 08:15:27.7   0.4     1.41 197    GII
PRNI  m   S     2004/02/11 08:15:46.6   1.5     1.41 197    GII
KSDI  m   Pn    2004/02/11 08:15:29.9   1.6     1.49   5    GII
KSDI  a   P     2004/02/11 08:15:28.8   0.3     1.49   5    ODC
MAMC  ad  P     2004/02/11 08:16:03.9   0.4     3.95 331    CYP
UZH   m i P     2004/02/11 08:19:33.0   1.1    19.64 333    LVV
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


		if (!in_hypocenter_event) {
			// find phase line
			while (strncmp(line, "Sta       Phase", 15) != 0) {
				// read next line
				cstat = fgets(line, MAXLINE_LONG, fp_obs);
				if (cstat == NULL)
					return(OBS_FILE_END_OF_INPUT);
			}
			// skip next line */
			cstat = fgets(line, MAXLINE_LONG, fp_obs);
			if (cstat == NULL)
				return(OBS_FILE_END_OF_INPUT);
			in_hypocenter_event = 1;
		}


		/* read phase arrival input */
/*
PRNI  m   Pn    2004/02/11 08:15:27.7   0.4     1.41 197    GII
PRNI  m   S     2004/02/11 08:15:46.6   1.5     1.41 197    GII
KSDI  m   Pn    2004/02/11 08:15:29.9   1.6     1.49   5    GII
KSDI  a   P     2004/02/11 08:15:28.8   0.3     1.49   5    ODC
MAMC  ad  P     2004/02/11 08:16:03.9   0.4     3.95 331    CYP
UZH   m i P     2004/02/11 08:19:33.0   1.1    19.64 333    LVV
*/
		istat = ReadFortranString(line, 1, 4, arrival->label);
		TrimString(arrival->label);
		istat += ReadFortranString(line, 11, 5, arrival->phase);
		TrimString(arrival->phase);
		istat += ReadFortranString(line, 8, 1, arrival->first_mot);
		if (arrival->first_mot[0] == ' ')
			arrival->first_mot[0] = ARRIVAL_NULL_CHR;
		istat += ReadFortranString(line, 9, 1, arrival->onset);
		if (arrival->onset[0] == ' ')
			arrival->onset[0] = ARRIVAL_NULL_CHR;
		// set weight, check for P - reset weight
		arrival->quality = arrival->onset[0] == 'i' ? 0 : 1;
		if (strncmp(arrival->phase, "P", 1) != 0) {
			arrival->quality += 1;
		}
		istat += ReadFortranInt(line, 17, 4, &arrival->year);
		istat += ReadFortranInt(line, 22, 2, &arrival->month);
		istat += ReadFortranInt(line, 25, 2, &arrival->day);
		istat += ReadFortranInt(line, 28, 2, &arrival->hour);
		istat += ReadFortranInt(line, 31, 2, &arrival->min);
		istat += ReadFortranReal(line, 34, 5, &arrival->sec);

		if (istat != 10) {
			return(OBS_FILE_END_OF_EVENT);
		}


		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}



	else if (strcmp(ftype_obs, "NCEDC_UCB") == 0)
	{

/* example:
...
19990818 010619.0000 +37.8950 -122.7000 010.0000 0.00 000 4.87 008 4.99 040 0.00 000 0.000e+00 000 099 180 009.47 00.0300 00.2000 00.2000 00.4000 00.1700 x F
$COM Bolinas, 19 km WSW of San Rafael, CA
$COM BRK: Mo=1.23E+23 dyne-cm, Mw=4.7.  Minor damage reported in Bolinas and San Rafael; felt from Santa Rosa in the north to San Jose in the
$COM south.
$AMP CMSB CL   2 WAS  407.40 001.10 0.00e+00
$AMP CMSB CL   3 WAS  596.20 000.77 0.00e+00
$PHS CMSB CL   1 i P        d 19990818 010625.9300 0039.535 093.58 000.00 000.00 0.00e+00 Y
$PHS CMSB CL   3 e S        x 19990818 010630.8700 0039.535 093.58 000.00 000.00 0.00e+00 N
$AMP CRQB CL   2 WAS  972.90 000.76 0.00e+00
$AMP CRQB CL   3 WAS  1122.00 000.82 0.00e+00
$PHS CRQB CL   1 i P        d 19990818 010627.5600 0045.409 066.72 000.00 000.00 0.00e+00 N
$PHS CRQB CL   3 i S        x 19990818 010633.6500 0045.409 066.72 000.00 000.00 0.00e+00 N
*/
		/* read line */
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
		if (cstat == NULL )
			return(OBS_FILE_END_OF_INPUT);
		/* check for end of event (assumes empty lines between event) */
		if (LineIsBlank(line)  || strncmp(line, "$END", 4) == 0) {
			/* end of event */
			return(OBS_FILE_END_OF_EVENT);
		}


		// find next phase line
		while (strncmp(line, "$PHS", 4) != 0) {
			// read next line
			cstat = fgets(line, MAXLINE_LONG, fp_obs);
			if (cstat == NULL || strncmp(line, "$END", 4) == 0)
				return(OBS_FILE_END_OF_EVENT);
		}


		/* read phase arrival input */
/* example:
$PHS CRQB CL   1 i P        d 19990818 010627.5600 0045.409 066.72 000.00 000.00 0.00e+00 N
$PHS CRQB CL   3 i S        x 19990818 010633.6500 0045.409 066.72 000.00 000.00 0.00e+00 N
$PHS MHC  VSP  Z i P        c 19990818 010636.8500 0111.833 123.00 000.00 000.00 0.00e+00 N
$PHS MHC  VSP  E i S        x 19990818 010651.9000 0111.833 123.00 000.00 000.00 0.00e+00 N
*/
/*
$PHS ZSP  TS13 Z i P        d 19890212 195041.9000 0032.674 068.58 000.00 000.00 0.00e+00 N
$PHS ZSP  TS13 Z i S        x 19890212 195047.4000 0032.674 068.58 000.00 000.00 0.00e+00 N
$PHS NFI  x    Z 0 P        c 19890212 195042.2500 0038.265 246.22 000.00 000.00 0.00e+00 Y
$PHS JEG  x    Z 0 P        d 19890212 195042.4500 0038.037 161.07 000.00 000.00 0.00e+00 N
$PHS PCC  B14  Z i P        c 19890212 195042.8000 0232.545 146.11 000.00 000.00 0.00e+00 Y
$PHS PCC  B14  Z i S        x 19890212 195048.7000 0232.545 146.11 000.00 000.00 0.00e+00 N
$PHS NPR  x    Z 3 P        c 19890212 195044.1300 0040.492 295.89 000.00 000.00 0.00e+00 N
*/
		istat = ReadFortranString(line, 6, 4, arrival->label);
		TrimString(arrival->label);
		istat += ReadFortranString(line, 11, 4, arrival->inst);
		TrimString(arrival->inst);
		istat += ReadFortranString(line, 16, 1, arrival->comp);
		TrimString(arrival->comp);
		istat += ReadFortranString(line, 18, 1, arrival->onset);
		if (arrival->onset[0] == '?')
			arrival->onset[0] = ARRIVAL_NULL_CHR;
		istat += ReadFortranString(line, 20, 8, arrival->phase);
		TrimString(arrival->phase);
		istat += ReadFortranString(line, 29, 1, arrival->first_mot);
		if (arrival->first_mot[0] == 'x')
			arrival->first_mot[0] = ARRIVAL_NULL_CHR;
		// set weight, check for P - reset weight
		if (arrival->onset[0] == 'x' || arrival->onset[0] == 'X')
			arrival->quality = 1;
		if (arrival->onset[0] == 'i')
			arrival->quality = 0;
		if (arrival->onset[0] == 'e')
			arrival->quality = 1;
		if (arrival->onset[0] == '?')
			arrival->quality = 1;
		if (arrival->onset[0] == '0')
			arrival->quality = 0;
		if (arrival->onset[0] == '1')
			arrival->quality = 1;
		if (arrival->onset[0] == '2')
			arrival->quality = 2;
		if (arrival->onset[0] == '3')
			arrival->quality = 3;
		if (arrival->onset[0] == '4')
			arrival->quality = 4;
		if (strncmp(arrival->phase, "P", 1) != 0) {
			arrival->quality += 1;
		}
		istat += ReadFortranInt(line, 31, 4, &arrival->year);
		istat += ReadFortranInt(line, 35, 2, &arrival->month);
		istat += ReadFortranInt(line, 37, 2, &arrival->day);
		istat += ReadFortranInt(line, 40, 2, &arrival->hour);
		istat += ReadFortranInt(line, 42, 2, &arrival->min);
		istat += ReadFortranReal(line, 44, 7, &arrival->sec);

		if (istat != 12) {
			return(OBS_FILE_END_OF_EVENT);
		}


		/* convert quality to error */
		Qual2Err(arrival);

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


		arrival->year = EventTime.year;
		arrival->month = EventTime.month;
		arrival->day = EventTime.day;


		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}



	else if (strcmp(ftype_obs, "HYPODD_PHA") == 0)
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
//		if (strlen(chrtmp) > 4 && strncmp(chrtmp, "NC", 2) == 0)
//			strcpy(arrival->label, chrtmp + 2);
//		else
			strcpy(arrival->label, chrtmp);
		
		arrival->hour = EventTime.hour;
		arrival->min = EventTime.min;
		arrival->sec = EventTime.sec + ttime;

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


	// DD
	else if (strcmp(ftype_obs, "HYPODD_CC") == 0 || strcmp(ftype_obs, "HYPODD_CT") == 0
		|| strcmp(ftype_obs, "HYPODD_") == 0)
	{

/* HYPODD_CC example:
#    28136    46442     -0.174000
NCCMO    -0.038245201    0.63    P
NCCRP    -0.062200069    0.78    P
NCJBG    -0.020649910    0.60    P
NCJPS    -0.011000156    0.50    P
*/
/* OTC - Origin time correction relative to the event origin time reported
	in the catalog data. If not known for all available data,
	then cross correlation data can not be used in combination
	with catalog data. Set OTC to -999 if not known for individual
	observations; Set to 0.0 if cross-correlation and catalog origin
	times are identical, or if only cross-correlation data is used.
*/
/* HYPODD_CT example:
#     38542     38520
NCCBW     3.430   3.430 1.0000 P
NCCSP     2.850   2.840 1.0000 P
NCNLN     8.950   8.950 1.0000 P
NCCAI     3.440   3.420 1.0000 P
*/


		/* read line */
		cstat = fgets(line, MAXLINE_LONG, fp_obs);
		if (cstat == NULL)
			return(OBS_FILE_END_OF_INPUT);


		/* read events line */
		if (line[0] == '#') {
			hypo_cc_flag = 1;
			istat = sscanf(line, "# %ld %ld %lf",
				&dd_event_id_1, &dd_event_id_2, &dd_otime_corr);
			if (dd_event_id_1 == dd_event_id_2)
				printf("ERROR: dd_event_id_1 == &dd_event_id_2  %ld %ld\n", dd_event_id_1, dd_event_id_2);
			if (istat == 2)
				hypo_cc_flag = 0;

			if (istat >= 2) {
				in_hypocenter_event = 1;
				/* read next line */
				cstat = fgets(line, MAXLINE_LONG, fp_obs);
				if (cstat == NULL)
					return(OBS_FILE_END_OF_INPUT);
			} else {
				printf("ERROR: format error: %s\n", line);
				return(OBS_FILE_FORMAT_ERROR);
			}
		}


		// read phase arrival input

		arrival->xcorr_flag = hypo_cc_flag;

		if (hypo_cc_flag) {	// HYPODD_CC

			istat = sscanf(line, "%s %lf %lf %s",
				arrival->label, &arrival->dd_dtime, &arrival->weight, arrival->phase);
	//printf("%s %lf %lf %s\n", arrival->label, arrival->dd_dtime, arrival->weight, arrival->phase);
			if (istat != 4) {
				printf("ERROR: HYPODD_CC format error: %s\n", line);
				return(OBS_FILE_FORMAT_ERROR);
			}

			arrival->dd_event_id_1 = dd_event_id_1;
			arrival->dd_event_id_2 = dd_event_id_2;
			// incorporate OTC into dd_dtime when reading dt file
			arrival->dd_dtime -= dd_otime_corr;
			strcpy(arrival->error_type, "XCC");

		} else {	// HYPODD_CT

			istat = sscanf(line, "%s %lf %lf %lf %s",
				arrival->label, &tt_sta1, &tt_sta2, &arrival->weight, arrival->phase);
	//printf("%s %lf %lf %lf %s\n", arrival->label, tt_sta1, tt_sta2, arrival->weight, arrival->phase);
			if (istat != 5) {
				printf("ERROR: HYPODD_CT format error: %s\n", line);
				return(OBS_FILE_FORMAT_ERROR);
			}

			arrival->dd_event_id_1 = dd_event_id_1;
			arrival->dd_event_id_2 = dd_event_id_2;
			arrival->dd_dtime = tt_sta1 - tt_sta2;
			strcpy(arrival->error_type, "CAT");

		}

// add noise
//arrival->dd_dtime += 0.02 * get_rand_double(-1.0, 1.0);
//

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


		arrival->quality = 0;
		arrival->year = EventTime.year;
		arrival->month = EventTime.month;
		arrival->day = EventTime.day;
		if (arrival->hour != EventTime.hour)
			puterr("WARNING: filename and arrival hours do not match.");
		if (arrival->min != EventTime.min)
			puterr("WARNING: filename and arrival minutes do not match.");


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

		/* convert quality to error */
		Qual2Err(arrival);

		return(istat);
	}


	else
	{
		puterr2("ERROR: unrecognized observation file type", ftype_obs);
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


/*** function to check if phase_in matches phase_check
       true if phase_in == phase_check
       true if phase_in is in phase ID list for phase_check */

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

	// check if phase_in == phase_check (AJL 20041201 added)
	if (strcmp(phase_in, phase_check) == 0)
		return(1);

	// check if phase_in is in phase ID list for phase_check
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

int EvalPhaseID(char *phase_out, char *phase_in)
{
	int npha_id;

	for (npha_id = 0; npha_id < NumPhaseID; npha_id++) {
		if (IsPhaseID(phase_in, PhaseID[npha_id].phase)) {
			if (phase_out != NULL)
				strcpy(phase_out, PhaseID[npha_id].phase);
			return(1);
		}
	}

	if (phase_out != NULL)
		strcpy(phase_out, phase_in);

	return(0);	// returned unchanged
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
//			- (long double) arrival[narr].delay	// DELAY_CORR	- incorporating delay so subtract (Tcorr = Tobs - (O-C))
			+ 60.0L * ((long double) arrival[narr].min
			+ 60.0L * (long double) arrival[narr].hour);

	if (!FixOriginTimeFlag) {
		/* initialize hypocenter year/month/day if origin time not fixed */
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




/*** function to check for arrivals with no absolute timing */

int CheckAbsoluteTiming(ArrivalDesc *arrival, int num_arrivals)
{
	int narr;
	int nNoAbs = 0;

	for (narr = 0; narr < num_arrivals; narr++) {
		if (arrival[narr].inst[0] == '*') {
			arrival[narr].abs_time = 0;
			nNoAbs++;
		} else {
			arrival[narr].abs_time = 1;
		}

	}

	return(nNoAbs);


}




/*** function to standardize date / time of arrivals */

int StdDateTime(ArrivalDesc *arrival, int num_arrivals, HypoDesc* phypo)
{
	int narr;
	double rms_resid = 0.0, weight_sum = 0.0;
	long double sec_tmp, hyp_time_tmp;


	for (narr = 0; narr < num_arrivals; narr++) {
		/* calc obs travel time and residual */
		if (arrival[narr].abs_time) {
			arrival[narr].obs_travel_time =
				arrival[narr].obs_time - phypo->time;
			arrival[narr].residual = arrival[narr].obs_travel_time -
				arrival[narr].pred_travel_time;
			rms_resid += arrival[narr].weight *
				arrival[narr].residual * arrival[narr].residual;
			weight_sum += arrival[narr].weight;
		} else {
			arrival[narr].obs_travel_time = 0.0;
			arrival[narr].residual = 0.0;
		}
		/* convert time to year/month/day/hour/min */
// DELAY_CORR		sec_tmp = arrival[narr].obs_time;
		// removing delay so add (Tobs = Tcorr + (O-C))
		sec_tmp = arrival[narr].obs_time + arrival[narr].delay;
		arrival[narr].hour = (int) (sec_tmp / 3600.0L);
		sec_tmp -= (long double) arrival[narr].hour * 3600.0L;
		arrival[narr].min = (int) (sec_tmp / 60.0L);
		sec_tmp -= (long double) arrival[narr].min * 60.0L;
		arrival[narr].sec = (double) sec_tmp;
		MonthDay(arrival[narr].year, arrival[narr].day_of_year,
			&(arrival[narr].month), &(arrival[narr].day));
	}

	phypo->rms = 999.99;
	if (weight_sum > 0.0)
		phypo->rms = sqrt(rms_resid / weight_sum);

	hyp_time_tmp = phypo->time;
	phypo->hour = (int) (hyp_time_tmp / 3600.0L);
	hyp_time_tmp -= (long double) phypo->hour * 3600.0L;
	phypo->min = (int) (hyp_time_tmp / 60.0L);
	hyp_time_tmp -= (long double) phypo->min * 60.0L;
	phypo->sec = (double) hyp_time_tmp;

	return(0);

}




/*** function to set output file root name using arrival time */

int SetOutName(ArrivalDesc *arrival, char* out_file_root, char* out_file,
	char lastfile[FILENAME_MAX], int isec)
{

  /*	if (isec)
		sprintf(out_file, "%s.%4.4d%2.2d%2.2d.%2.2d%2.2d%2.2d",
		out_file_root, arrival->year, arrival->month, arrival->day,
		arrival->hour, arrival->min, (int) arrival->sec);
	else
		sprintf(out_file, "%s.%4.4d%2.2d%2.2d.%2.2d%2.2d",
		out_file_root, arrival->year, arrival->month, arrival->day,
		arrival->hour, arrival->min);  */
  /* SH 03/28/02 change to include digits of sec to construct filename */
        if (isec)
		sprintf(out_file, "%s.%4.4d%2.2d%2.2d.%2.2d%2.2d%05.2f",
		out_file_root, arrival->year, arrival->month, arrival->day,
		arrival->hour, arrival->min, arrival->sec);
	else
		sprintf(out_file, "%s.%4.4d%2.2d%2.2d.%2.2d%2.2d%2.2d",
		out_file_root, arrival->year, arrival->month, arrival->day,
		arrival->hour, arrival->min, (int) arrival->sec);
  /* SH 04/08/02 check if same filename as previous event;
        if so append 'b' to filename;
        identical filenames can happen with swarm data */
        if (strcmp(out_file, lastfile) == 0)
		strcat(out_file, "b");

	return(0);

}


/*** function to check for duplicate label and phase in arrival */

int IsDuplicateArrival(ArrivalDesc *arrival, int num_arrivals, int ntest)
{
	int narr;

	for (narr = 0; narr < num_arrivals; narr++) {
		if (narr != ntest
				&& !strcmp(arrival[narr].time_grid_label, arrival[ntest].time_grid_label)
				&& !strcmp(arrival[narr].phase, arrival[ntest].phase))
			return(narr);
	}

	return(-1);

}



/*** function to check for duplicate label and phase in arrival */

int IsSameArrival(ArrivalDesc *arrival, int num_arrivals, int ntest, char *phase_test)
{
	int narr;

	if (phase_test == NULL) {
		for (narr = 0; narr < num_arrivals; narr++) {
			if (narr != ntest
					&& ((IsPhaseID(arrival[narr].phase, "P") &&
							IsPhaseID(arrival[ntest].phase, "P"))
						|| (IsPhaseID(arrival[narr].phase, "S") &&
							IsPhaseID(arrival[ntest].phase, "S")))
					&& !strcmp(arrival[narr].time_grid_label, arrival[ntest].time_grid_label))
				return(narr);
		}
	} else {
		for (narr = 0; narr < num_arrivals; narr++) {
			if (narr != ntest
					&& !strcmp(arrival[narr].time_grid_label, arrival[ntest].time_grid_label)
					&& IsPhaseID(arrival[narr].phase, phase_test))
				return(narr);
		}
	}

	return(-1);

}



/*** function to check for duplicate label and phase in arrival */

int FindDuplicateTimeGrid(ArrivalDesc *arrival, int num_arrivals, int ntest)
{
	int narr;

	for (narr = 0; narr < num_arrivals; narr++) {
		if (narr != ntest
				&& !strcmp(arrival[narr].fileroot, arrival[ntest].fileroot)
				&& arrival[narr].flag_ignore == 0
				)
			return(narr);
	}

	return(-1);

}



/*** function to perform grid search location */

int LocGridSearch(int ngrid, int num_arr_total, int num_arr_loc,
		ArrivalDesc *arrival,
		GridDesc* ptgrid, GaussLocParams* gauss_par, HypoDesc* phypo)
{

	int istat;
	int ix = -1, iy = -1, iz = -1, narr, n_compan;
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

	putmsg(4, "");
	putmsg(4, "Calculating solution over grid...");

	iGridType = ptgrid->type;

	xval = ptgrid->origx;

	/* loop over grid points */

	for (ix = 0; ix <  ptgrid->numx; ix++) {

		/* read y-z sheets for arrival travel-times (3D grids) */
		if ((istat = ReadArrivalSheets(num_arr_loc, arrival, xval)) < 0)
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
				} else if (arrival[narr].sheetdesc.type == GRID_TIME) {
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
						iGridType, &misfit, NULL);
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
	sprintf(phypo->searchInfo, "GRID nPts %d%c", ix * iy *iz, '\0');
	/* write message */
	/*putmsg(2, phypo->searchInfo);*/


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
	long int ngenerated;
	int maxNumTries;
	int writeMessage = 0;
	int iGridType;
	int nReject, numClipped = 0, numGridReject = 0, numStaReject = 0;
	int iAbort = 0, iReject = 0;
	int iBoundary = 0;
	int iAccept, numAcceptDeepMinima = 0;
	double xval, yval, zval;
	double currentMetStepFact;

	double value, dlike, dlike_max = -VERY_LARGE_DOUBLE;

	double misfit;
	double misfit_min = VERY_LARGE_DOUBLE, misfit_max = -VERY_LARGE_DOUBLE;

	double xmin, xmax, ymin, ymax, zmin, zmax;
	double dx_init, dx_test;

	int nScatterSaved;

	double xmean_sum = 0.0, ymean_sum = 0.0, zmean_sum = 0.0;
	double xvar_sum = 0.0, yvar_sum = 0.0, zvar_sum = 0.0;
	double xvar = 0.0, yvar = 0.0, zvar = 0.0;
	double dsamp = 0.0, dsamp2;



	/* get solution quality at each sample on random walk */

	putmsg(4, "");
	putmsg(4, "Calculating solution along Metropolis walk...");

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
					iGridType, &misfit, NULL);
			dlike = gauss_par->WtMtrxSum * exp(value);

			/* apply Metropolis test */
			iAccept = MetropolisTest(pMetrop->likelihood, dlike);

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
		putmsg(4, MsgStr);
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
		sprintf(MsgStr, "WARNING: %d Metropolis samples clipped at search grid boundary.",
			numClipped);
		putmsg(1, MsgStr);
	}


	/* give warning if grid points rejected */

	if (numGridReject > 0) {
		sprintf(MsgStr, "WARNING: %d Metropolis samples rejected; travel times for an average of %.2lf arrival observations were not valid.",
			numGridReject, (double) numStaReject / numGridReject);
		putmsg(1, MsgStr);
	}


	/* check reject location conditions */

	/* maximum like hypo on edge of grid */
	if ((iBoundary = isOnGridBoundary(phypo->x, phypo->y, phypo->z,
			ptgrid, pMetrop->dx, pMetrop->dx, 0))) {
		sprintf(MsgStr, "WARNING: max prob location on grid boundary %d, rejecting location.", iBoundary);
		putmsg(1, MsgStr);
		sprintf(phypo->locStatComm, "%s", MsgStr);
		iReject = 1;
	}

	/* construct search information string */
	sprintf(phypo->searchInfo,
"METROPOLIS nSamp %ld nAcc %d nSave %d nClip %d Dstep0 %lf Dstep %lf%c",
		ngenerated, nSamples, nScatterSaved, numClipped, dx_init, pMetrop->dx, '\0');
	/* write message */
	putmsg(2, phypo->searchInfo);


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

INLINE int GetNextMetropolisSample(WalkParams* pMetrop, double xmin, double xmax,
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

INLINE int MetropolisTest(double likelihood_last, double likelihood_new)
{

	double prob;

	/* compare with last sample using Mosegaard & Tarantola eq (17) */

	if (likelihood_new >= likelihood_last)
		return(1);
	else if ((prob = get_rand_double(0.0, 1.0)) < likelihood_new / likelihood_last)
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
	double value, misfit, otime;

	SourceDesc station;



//printf("SaveBestLocation num_arr_total %d num_arr_loc %d\n", num_arr_total, num_arr_loc);

	/* loop over observed arrivals */
	for (narr = 0; narr < num_arr_total; narr++) {
		arrival[narr].dist = GetEpiDist(
				&(arrival[narr].station), phypo->x, phypo->y);
		if (GeometryMode == MODE_GLOBAL)
			arrival[narr].dist *= KM2DEG;
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
			// save staion information (will be overwritten in OpenGrid3dFile()
			station = arrival[narr].station;
			sprintf(filename, "%s.time", arrival[narr].fileroot);
			if ((istat = OpenGrid3dFile(filename,
					&(arrival[narr].fpgrid),
					&(arrival[narr].fphdr),
					&(arrival[narr].gdesc), "time",
					&(arrival[narr].station),
					iSwapBytesOnInput)) < 0)
				continue;
			arrival[narr].station = station;
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
			arrival[narr].pred_travel_time *= arrival[narr].tfact;
			// apply crustal correction
			if (ApplyCrustElevCorrFlag
					&& GeometryMode == MODE_GLOBAL
					&& arrival[narr].pred_travel_time > 0.0) {
				if (arrival[narr].dist > MinDistCrustElevCorr)
					arrival[narr].pred_travel_time +=
						applyCrustElevCorrection(arrival + narr, phypo->x, phypo->y, phypo->z);
			}
			else if (ApplyElevCorrFlag) {
				arrival[narr].pred_travel_time += arrival[narr].elev_corr;
			}
			CloseGrid3dFile(&(arrival[narr].fpgrid), &(arrival[narr].fphdr));

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
					&(arrival[narr].ray_qual), -1.0, iSwapBytesOnInput);
			} else {
				/* 2D grid (1D model) */
				ReadTakeOffAnglesFile(filename,
					0.0, arrival[narr].dist, phypo->z,
					&(arrival[narr].ray_azim),
					&(arrival[narr].ray_dip),
					&(arrival[narr].ray_qual), arrival[narr].azim, iSwapBytesOnInput);
			}
		}

//		/* close time grid file for ignored arrivals */
//		if (iopened)
//			CloseGrid3dFile(&(arrival[narr].fpgrid),
//				&(arrival[narr].fphdr));
	}

	/* calc misfit or prob density */
	value = CalcSolutionQuality(num_arr_loc, arrival, gauss_par, iGridType, &misfit, &otime);

	/* set origin time */
	if (!FixOriginTimeFlag)
		phypo->time = otime;
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
		if (arrival[narr].gdesc.numx > 1)
			ixsheet = (int) ((xsheet - arrival[narr].gdesc.origx)
						/ arrival[narr].gdesc.dx);
		else
			ixsheet = 0;
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

	return(0);

}




/*** function to consruct weight matrix (inverse of covariance matrix) */

int ConstWeightMatrix(int num_arrivals, ArrivalDesc *arrival,
		GaussLocParams* gauss_par)
{
	int istat, nrow, ncol;
	double sigmaT2, corr_len2;
	double dx, dy, dz, dist2;
	double weight_sum;
	double sta_wt;
	DMatrix wt_matrix, edt_matrix, null_mtrx = NULL;
	SourceDesc *sta1, *sta2;


	/* allocate square matrix  */

	edt_matrix = dmatrix(0, num_arrivals - 1, 0, num_arrivals - 1);
	wt_matrix = dmatrix(0, num_arrivals - 1, 0, num_arrivals - 1);


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
		// same phase types, include spatial correlation  
			dx = sta1->x - sta2->x;
			dy = sta1->y - sta2->y;
			dz = sta1->z - sta2->z;
			dist2 = dx * dx + dy * dy + dz * dz;
			if (GeometryMode == MODE_GLOBAL)
				dist2 *= DEG2KM * DEG2KM;
			// EDT
			if (ncol == nrow)	// diagonal of EDT gets gaussian model time error
				edt_matrix[nrow][ncol] = sigmaT2 * exp(-0.5 * dist2 / corr_len2);
			else			// off-diagonal of EDT gets gaussian model weight
				edt_matrix[nrow][ncol] = edt_matrix[ncol][nrow] = exp(-0.5 * dist2 / corr_len2);
			// LS/L2 gaussian model time error
			wt_matrix[nrow][ncol] = wt_matrix[ncol][nrow] = sigmaT2 * exp(-0.5 * dist2 / corr_len2);
		} else { 
		// different phase types, assumed no spatial correlation  
			edt_matrix[nrow][ncol] = edt_matrix[ncol][nrow] = 0.0;
			wt_matrix[nrow][ncol] = wt_matrix[ncol][nrow] = 0.0;
		}

		/* obs time error */
		if (ncol == nrow) {
			edt_matrix[nrow][ncol] += arrival[nrow].error * arrival[nrow].error;
			wt_matrix[nrow][ncol] += arrival[nrow].error * arrival[nrow].error;
		}

		}
	}

	if (message_flag >= 5)
		DisplayDMatrix("Covariance", wt_matrix, num_arrivals, num_arrivals);


	/* invert covariance matrix to obtain weight matrix */

	if ((istat = dgaussj(wt_matrix, num_arrivals, null_mtrx, 0)) < 0) {
		puterr("ERROR: inverting covariance matrix.");
		return(-1);
	}

	if (message_flag >= 5)
		DisplayDMatrix("Weight", wt_matrix, num_arrivals, num_arrivals);
			

	// station distance weighting
	
	for (nrow = 0; nrow < num_arrivals; nrow++) {
		for (ncol = 0; ncol <= nrow; ncol++) {
			sta_wt = (arrival[nrow].station_weight + arrival[ncol].station_weight) / 2.0;
			wt_matrix[nrow][ncol] *= sta_wt;
			wt_matrix[ncol][nrow] *= sta_wt;
		}
	}
	
	/* get row weights & sum of weights */

	weight_sum = 0.0;
	for (nrow = 0; nrow < num_arrivals; nrow++) {
		arrival[nrow].weight = 0.0;
		for (ncol = 0; ncol < num_arrivals; ncol++) {
			arrival[nrow].weight += wt_matrix[nrow][ncol];
			weight_sum += wt_matrix[nrow][ncol];
		}
	}
	for (nrow = 0; nrow < num_arrivals; nrow++) {
		arrival[nrow].weight = (double) num_arrivals * arrival[nrow].weight / weight_sum;
		if (arrival[nrow].weight < 0.0) {
			sprintf(MsgStr,
"ERROR: negative observation weight: %s %s %s weight: %lf",
				arrival[nrow].label, arrival[nrow].inst,
				arrival[nrow].comp, arrival[nrow].weight);
			puterr(MsgStr);
			puterr("   Gaussian model error (see LOCGAU) may be too large relative to obs uncertainty (see LOCQUAL2ERR, or NLL-Phase format ErrMag).");
		}
	}
	sprintf(MsgStr, "Weight Matrix sum: %lf", weight_sum);
	putmsg(4, MsgStr);

	/* set global variables */

	gauss_par->EDTMtrx = edt_matrix;
	gauss_par->WtMtrx = wt_matrix;
	gauss_par->WtMtrxSum = weight_sum;

	return(0);

}



/*** function to calculate weighted mean of observed arrival times */
/*		(TV82, eq. A-38) */

void CalcCenteredTimesObs(int num_arrivals, ArrivalDesc *arrival,
				GaussLocParams* gauss_par, HypoDesc* phypo)
{

	int nrow, ncol, narr;
	long double sum, weighted_mean;
	DMatrix wtmtx;
	double *wtmtxrow, wt_sum;


	if (!FixOriginTimeFlag) {

		/* calculate weighted mean of observed times */

		wtmtx = gauss_par->WtMtrx;
		sum = 0.0L;
		wt_sum = 0.0;

		for (nrow = 0; nrow < num_arrivals; nrow++) {
			if (!arrival[nrow].abs_time)
				continue;	// ignore obs without absolute timing
			wtmtxrow =  wtmtx[nrow];
			for (ncol = 0; ncol < num_arrivals; ncol++) {
				if (!arrival[ncol].abs_time)
					continue;	// ignore obs without absolute timing
				sum += (long double) *(wtmtxrow + ncol) * arrival[ncol].obs_time;
				wt_sum += (double) *(wtmtxrow + ncol);
			}
		}
		if (wt_sum > 0.0)
			weighted_mean = sum / (long double) wt_sum;
		else
			weighted_mean = (long double) arrival[0].obs_time;

	} else {

		/* use fixed origin time as reference */

		weighted_mean = phypo->time;
	}


	/* set centered observed times */

	putmsg(3, "");
	putmsg(3, "Delayed, Sorted, Centered Observations:");
	for (narr = 0; narr < num_arrivals; narr++) {
		arrival[narr].obs_centered =
			(double) (arrival[narr].obs_time - weighted_mean);
		sprintf(MsgStr,
"  %3d  %-6s %-6s %2.2d:%2.2d:%7.4lf - %7.4lfs -> %8.4lf (%10.4lf)",
				narr, arrival[narr].label, arrival[narr].phase,
				arrival[narr].hour, arrival[narr].min,
				arrival[narr].sec, arrival[narr].delay, arrival[narr].obs_centered,
				((double) arrival[narr].obs_time));
		putmsg(3, MsgStr);
	}

	gauss_par->meanObs = weighted_mean;

}


/*** function to calculate weighted mean of predicted travel times */
/*		(TV82, eq. A-38) */

INLINE void CalcCenteredTimesPred(int num_arrivals, ArrivalDesc *arrival,
				GaussLocParams* gauss_par)
{

	int nrow, ncol, narr;
	double sum, weighted_mean, pred_time_row;
	DMatrix wtmtx;
	double *wtmtxrow, wt_sum;


	if (!FixOriginTimeFlag) {

		wtmtx = gauss_par->WtMtrx;
		sum = 0.0;
		wt_sum = 0.0;

		for (nrow = 0; nrow < num_arrivals; nrow++) {
			// AJL 20041115 bug fix!
			if (arrival[nrow].pred_travel_time <= 0.0)
				continue;	// ignore obs without predicted times
			// END
			if (!arrival[nrow].abs_time)
				continue;	// ignore obs without absolute timing
			wtmtxrow =  wtmtx[nrow];
			pred_time_row = arrival[nrow].pred_travel_time;
			for (ncol = 0; ncol < num_arrivals; ncol++) {
				// AJL 20041115 bug fix!
				if (arrival[ncol].pred_travel_time <= 0.0)
					continue;	// ignore obs without predicted times
				// END
				if (!arrival[ncol].abs_time)
					continue;	// ignore obs without absolute timing
				sum += (double) *(wtmtxrow + ncol) * pred_time_row;
				wt_sum += (double) *(wtmtxrow + ncol);
			}
		}

		if (wt_sum > 0.0)
			weighted_mean = sum / wt_sum;
		else
			weighted_mean = (long double) arrival[0].pred_travel_time;

	} else {

		/* for fixed origin time use travel time directly */
		weighted_mean = 0.0;

	}



	/* set centered predicted times */

	for (narr = 0; narr < num_arrivals; narr++) {
		// AJL 20041115 bug fix!
		if (arrival[narr].pred_travel_time <= 0.0)
			continue;	// ignore obs without predicted times
		// END
		arrival[narr].pred_centered =
			arrival[narr].pred_travel_time - weighted_mean;
	}


	gauss_par->meanPred = (double) weighted_mean;

}



/*** function to calculate probability density */
/*		(MEN92, eq. 14) */

INLINE double CalcSolutionQuality(int num_arrivals, ArrivalDesc *arrival,
	GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime)
{

	if (LocMethod == METH_GAU_ANALYTIC) {
		return (CalcSolutionQuality_GAU_ANALYTIC(num_arrivals, arrival,
			gauss_par, itype, pmisfit, potime));
	} else if (LocMethod == METH_GAU_TEST) {
		return (CalcSolutionQuality_GAU_TEST(num_arrivals, arrival,
			gauss_par, itype, pmisfit, potime));
	} else if (LocMethod == METH_EDT) {
		return (CalcSolutionQuality_EDT(num_arrivals, arrival,
			gauss_par, itype, pmisfit, potime, 0));
	} else if (LocMethod == METH_EDT_BOX) {
		return (CalcSolutionQuality_EDT(num_arrivals, arrival,
			gauss_par, itype, pmisfit, potime, 1));
	} else {
		return(-1.0);
	}

}



/*** function to calculate probability density */
/*	EDT - sum of probabilities of difference of obs - difference of travel times
		for all pairs of obs
*/

INLINE double CalcSolutionQuality_EDT(int num_arrivals, ArrivalDesc *arrival,
	GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime, int method_box)
{

	int nrow, ncol;

	long double edt_sum, prob;

	double edt_misfit, edt_weight;
	double weight, weight2;
	double edtSumMisfit;
	double ln_prob_density, rms_misfit;

	DMatrix edtmtx;
	double sigma2_row;
	double obs_minus_pred;

	int no_abs_time_row;

	long double ot_sum, ot_weight;
	int icalc_otime = 0;
	
	double error_row, amp_row, unc_limit;



	if (potime != NULL)
		icalc_otime = 1;

	edtmtx = gauss_par->EDTMtrx;


	/* calculate weighted mean of predicted travel times  */
	/*		(TV82, eq. A-38) */

	CalcCenteredTimesPred(num_arrivals, arrival, gauss_par);


	/* calculate EDT prop sum */

	edt_sum = 0.0L;
	edt_weight = 0.0;
	ot_sum = 0.0L;
	ot_weight = 0.0L;
	for (nrow = 1; nrow < num_arrivals; nrow++) {
		// AJL 20041115 bug fix!
		if (arrival[nrow].pred_travel_time <= 0.0)
			continue;	// ignore obs without predicted times
		// END
		sigma2_row = edtmtx[nrow][nrow];
		error_row = arrival[nrow].error;
		amp_row = arrival[nrow].amplitude;
		obs_minus_pred = arrival[nrow].obs_centered - arrival[nrow].pred_centered;
		no_abs_time_row = !arrival[nrow].abs_time;
		for (ncol = 0; ncol < nrow; ncol++) {
			// AJL 20041115 bug fix!
			if (arrival[ncol].pred_travel_time <= 0.0)
				continue;	// ignore obs without predicted times
			// END
			// check abslolute timing
			if (no_abs_time_row) {
				if (arrival[ncol].abs_time)	// cannot be same station/inst
					continue;
				if (strcmp(arrival[nrow].label, arrival[ncol].label) != 0
						|| strcmp(arrival[nrow].inst, arrival[ncol].inst) != 0)
					continue;	// not same sta/inst
			}
			// calculate EDT misfit:  (obs1 - obs2) - (pred1 - pred2)
			edt_misfit = (double) (obs_minus_pred + arrival[ncol].pred_centered - arrival[ncol].obs_centered);
			// calculate probability
			if (method_box) {
				unc_limit = amp_row + arrival[ncol].amplitude;	// sum of mean pick unc for each box
				//unc_limit = error_row + arrival[ncol].error;	// sum of box widths
				prob = fabs(edt_misfit) <= unc_limit ? 1.0 : 0.0;
				weight = amp_row * arrival[ncol].amplitude;	// product of mean pick unc for each box
				weight *= (1.0 - edtmtx[nrow][ncol]);		// correlation coeff
			} else {
				weight2 = 1.0 / (sigma2_row + edtmtx[ncol][ncol]);	// sum of errors**2
				prob = exp(-0.5 * edt_misfit * edt_misfit * weight2);
				weight = sqrt(weight2);		// errors factor
				weight *= (1.0 - edtmtx[nrow][ncol]);		// correlation coeff

			}
			if (iSetStationDistributionWeights)
				weight *= (arrival[nrow].station_weight + arrival[ncol].station_weight) / 2.0;
			prob *= weight;
			edt_sum += prob;
			edt_weight += weight;
			// otime
			if (icalc_otime) {
				ot_sum += prob * (
					arrival[nrow].obs_time - (long double) arrival[nrow].pred_travel_time
						+ arrival[ncol].obs_time - (long double) arrival[ncol].pred_travel_time);
				ot_weight += 2.0L * prob;
			}
		}
	}

	if (edt_sum < SMALL_DOUBLE)
		edt_sum = SMALL_DOUBLE;

	if (edt_weight > 0.0 && num_arrivals > 0) {
		edtSumMisfit = -log(edt_sum / edt_weight);
		//edtSumMisfit *= 2.0;
		edtSumMisfit *= num_arrivals;
		//edtSumMisfit *= sqrt(num_arrivals);	// works much better for INGV Italy region scale
		if (icalc_otime)
			*potime = ot_sum / ot_weight;
	} else {
		edtSumMisfit = VERY_LARGE_DOUBLE;
		if (potime != NULL)
			*potime = VERY_LARGE_DOUBLE;
	}
//printf("edt_sum %lf  edtSumMisfit %lf\n", edt_sum, edtSumMisfit);

	/* return misfit or ln(prob density) */

	if (itype == GRID_MISFIT) {
		/* convert misfit to rms misfit */
		rms_misfit = sqrt(edtSumMisfit / num_arrivals);
		*pmisfit = rms_misfit;
		return(rms_misfit);
	} else if (itype == GRID_PROB_DENSITY) {
		ln_prob_density = -edtSumMisfit;
		//ln_prob_density = -0.5 * edtSumMisfit;
		rms_misfit = sqrt(edtSumMisfit / num_arrivals);
		*pmisfit = rms_misfit;
		return(ln_prob_density);
	} else {
		return(-1.0);
	}



}


/*** function to calculate probability density */
/*	sum of individual L2 residual probablities */

INLINE double CalcSolutionQuality_GAU_TEST(int num_arrivals, ArrivalDesc *arrival,
			GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime)
{

	int nrow, ncol, narr;

	double misfit;
	double ln_prob_density, rms_misfit;

	double arr_row_res;
	DMatrix wtmtx;
	double *wtmtxrow;

	double pred_min = VERY_LARGE_DOUBLE, pred_max = -VERY_LARGE_DOUBLE;
	double prob_max = -VERY_LARGE_DOUBLE;
	double tshift, tshift_at_prob_max = 0.0;
	double tstep, tstart, tstop;
	double misfit_min = VERY_LARGE_DOUBLE;
	double misfit_tmp, prob;


	wtmtx = gauss_par->WtMtrx;


	/* calculate weighted mean of predicted travel times  */
	/*		(TV82, eq. A-38) */

	CalcCenteredTimesPred(num_arrivals, arrival, gauss_par);


	/* calculate residuals (TV82, eqs. 10-12, 10-13; MEN92, eq. 15) */

	for (narr = 0; narr < num_arrivals; narr++) {
		// AJL 20041115 bug fix!
		if (arrival[narr].pred_travel_time <= 0.0)
			arrival[narr].cent_resid = 0.0;		// ignore obs without predicted times
		else {
			arrival[narr].cent_resid = arrival[narr].obs_centered - arrival[narr].pred_centered;
			if (arrival[narr].pred_centered < pred_min)
				pred_min = arrival[narr].pred_centered;
			if (arrival[narr].pred_centered > pred_max)
				pred_max = arrival[narr].pred_centered;
		}
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
		// AJL 20041115 bug fix!
		if (arrival[nrow].pred_travel_time <= 0.0)
			continue;	// ignore obs without predicted times
		// END
				wtmtxrow = wtmtx[nrow];
				arr_row_res = arrival[nrow].cent_resid + tshift;
				for (ncol = 0; ncol <= nrow; ncol++) {
		// AJL 20041115 bug fix!
		if (arrival[ncol].pred_travel_time <= 0.0)
			continue;	// ignore obs without predicted times
		// END
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

	// get origin time
	if (potime != NULL)
		*potime = CalcMaxLikeOriginTime(num_arrivals, arrival, gauss_par);

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

INLINE double CalcSolutionQuality_GAU_ANALYTIC(int num_arrivals, ArrivalDesc *arrival,
			GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime)
{

	int nrow, ncol, narr;
	int narr_misfit;

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

	for (narr = 0; narr < num_arrivals; narr++) {
		// AJL 20041115 bug fix!
		if (arrival[narr].pred_travel_time <= 0.0)
			arrival[narr].cent_resid = 0.0;		// ignore obs without predicted times
		else
			arrival[narr].cent_resid = arrival[narr].obs_centered - arrival[narr].pred_centered;
	}

	narr_misfit = 0;
	for (nrow = 0; nrow < num_arrivals; nrow++) {
		// AJL 20041115 bug fix!
		if (arrival[nrow].pred_travel_time <= 0.0) {
//printf("IGNORE: %s %s\n", arrival[nrow].label, arrival[nrow].phase);
			continue;	// ignore obs without predicted times
		}
		// END
		if (!arrival[nrow].abs_time)
			continue;
		narr_misfit++;
		wtmtxrow = wtmtx[nrow];
		arr_row_res = arrival[nrow].cent_resid;
		for (ncol = 0; ncol <= nrow; ncol++) {
			// AJL 20041115 bug fix!
			if (arrival[ncol].pred_travel_time <= 0.0)
				continue;	// ignore obs without predicted times
			// END
			if (!arrival[ncol].abs_time)
				continue;
			if (ncol != nrow)
			    misfit += 2.0 * (double) *(wtmtxrow + ncol) *
				arr_row_res * arrival[ncol].cent_resid;
			else
			    misfit += (double) *(wtmtxrow + ncol) *
				arr_row_res * arrival[ncol].cent_resid;
		}
	}


	if (potime != NULL)
		*potime = CalcMaxLikeOriginTime(num_arrivals, arrival, gauss_par);

	/* return misfit or ln(prob density) */

	if (itype == GRID_MISFIT) {
		/* convert misfit to rms misfit */
		rms_misfit = sqrt(misfit / narr_misfit);
		*pmisfit = rms_misfit;
		return(rms_misfit);
	} else if (itype == GRID_PROB_DENSITY) {
		ln_prob_density = -0.5 * misfit;
		rms_misfit = sqrt(misfit / narr_misfit);
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

	return(gauss_par->meanObs - (long double) gauss_par->meanPred);

}



/*** function to update arrival probabilistic residuals */

INLINE void UpdateProbabilisticResiduals(int num_arrivals,
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
	putmsg(3, "");
	putmsg(3, "Calculating confidence intervals over grid...");


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
	double xnpoints;
	double probmax, prob_den;



	/* return if no scatter samples requested */
	if (pscat->npts < 1)
		return(0);

	/* write message */
	putmsg(3, "");
	putmsg(3, "Generating event scatter file...");

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
	putmsg(3, MsgStr);
	sprintf(MsgStr, "  (%d points requested, dvol= %lf, probmax=%lf)",
			pscat->npts, dvol, probmax);
	putmsg(3, MsgStr);

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
		flag_phstat = 0, flag_phase_id = 0, flag_sta_wt = 0, flag_qual2err = 0,
		flag_mag = 0, flag_alias = 0, flag_exclude = 0, flag_time_delay = 0,
		flag_time_delay_surface = 0, flag_elev_corr = 0,
		flag_otime = 0, flag_angles = 0, flag_source = 0;
	int flag_include = 1;


	// param lines not needed by NLDiffLoc
	if (nll_mode == MODE_DIFFERENTIAL) {
		flag_search = 1;
	}


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

		if (strcmp(param, "CONTROL") == 0) {
			if ((istat = get_control(strchr(line, ' '))) < 0)
				puterr("ERROR: readingng control params.");
			else
				flag_control = 1;
		}


		/* read grid params */

		if (strcmp(param, "LOCGRID") == 0) {
    			if ((istat = GetNLLoc_Grid(strchr(line, ' '))) < 0)
				puterr("ERROR: reading grid parameters.");
			else
				flag_grid = 1;
		}


		/* read file names */

		if (strcmp(param, "LOCFILES") == 0) {
			if ((istat = GetNLLoc_Files(strchr(line, ' '))) < 0)
			  puterr("ERROR: reading NLLoc output file name.");
			else
				flag_outfile = 1;
		}


		/* read output file types names */

		if (strcmp(param, "LOCHYPOUT") == 0) {
			if ((istat = GetNLLoc_HypOutTypes(strchr(line, ' '))) < 0)
				puterr("ERROR: reading NLLoc hyp output file types.");
			else
				flag_hyptype = 1;
		}


		/* read search type */

		if (strcmp(param, "LOCSEARCH") == 0) {
			if ((istat = GetNLLoc_SearchType(strchr(line, ' '))) < 0)
				puterr("ERROR: reading NLLoc search type.");
			else
				flag_search = 1;
		}


		/* read method */

		if (strcmp(param, "LOCMETH") == 0) {
			if ((istat = GetNLLoc_Method(strchr(line, ' '))) < 0)
				puterr("ERROR: reading NLLoc method.");
			else
				flag_method = 1;
		}


		/* read fixed origin time parameters */

		if (strcmp(param, "LOCFIXOTIME") == 0) {
			if ((istat = GetNLLoc_FixOriginTime(
						strchr(line, ' '))) < 0)
				puterr("ERROR: reading NLLoc fixed origin time params.");
			else
				flag_otime = 1;
		}


		/* read phase identifier values */

		if (strcmp(param, "LOCPHASEID") == 0) {
			if ((istat = GetPhaseID(strchr(line, ' '))) < 0)
				puterr("ERROR: reading phase identifier values.");
			else
				flag_phase_id = 1;
		}


		/* read station distance weighting values */

		if (strcmp(param, "LOCSTAWT") == 0) {
			if ((istat = GetStaWeight(strchr(line, ' '))) < 0)
				puterr("ERROR: reading station distance weighting values.");
			else
				flag_sta_wt = 1;
		}


		/* read quality2error values */

		if (strcmp(param, "LOCQUAL2ERR") == 0) {
			if ((istat = GetQuality2Err(strchr(line, ' '))) < 0)
				puterr("ERROR: reading quality2error values.");
			else
				flag_qual2err = 1;
		}


		/* read comment */

		if (strcmp(param, "LOCCOM") == 0) {
			strcpy(Hypocenter.comment, strchr(line, ' ') + 1);
			*(strchr(Hypocenter.comment, '\n')) = '\0';
			sprintf(MsgStr, "LOCCOMMENT:  %s\n",
				Hypocenter.comment);
			putmsg(3, MsgStr);
			flag_comment = 1;
		}


		/* read signature */

		if (strcmp(param, "LOCSIG") == 0) {
			strcpy(LocSignature, strchr(line, ' ') + 1);
			*(strchr(LocSignature, '\n')) = '\0';
			sprintf(MsgStr, "LOCSIGNATURE:  %s\n",
				LocSignature);
			putmsg(3, MsgStr);
			flag_signature = 1;
		}


		/* read gauss params */

		if (strcmp(param, "LOCGAU") == 0) {
			if ((istat = GetNLLoc_Gaussian(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Gaussian parameters.");
			else
				flag_gauss = 1;
		}



		/* read phase statistics params */

		if (strcmp(param, "LOCPHSTAT") == 0) {
			if ((istat = GetNLLoc_PhaseStats(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Phase Statistics parameters.");
			else
				flag_phstat = 1;
		}


		/* read take-off angles params */

		if (strcmp(param, "LOCANGLES") == 0) {
			if ((istat = GetNLLoc_Angles(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Take-off Angles parameters.");
			else
				flag_angles = 1;
		}


		/* read magnitude calculation params */

		if (strcmp(param, "LOCMAG") == 0) {
			if ((istat = GetNLLoc_Magnitude(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Magnitude Calculation parameters.");
			else
				flag_mag = 1;
		}


		/* read component params */

		if (strcmp(param, "LOCCMP") == 0) {
			if ((istat = GetCompDesc(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Component description parameters.");
			else
				flag_comp = 1;
		}


		/* read alias params */

		if (strcmp(param, "LOCALIAS") == 0) {
			if ((istat = GetLocAlias(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Alias parameters.");
			else
				flag_alias = 1;
		}


		/* read exclude params */

		if (strcmp(param, "LOCEXCLUDE") == 0) {
			if ((istat = GetLocExclude(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Exclude parameters.");
			else
				flag_exclude = 1;
		}


		/* read station time delay params */

		if (strcmp(param, "LOCDELAY") == 0) {
			if ((istat = GetTimeDelays(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Time Delay parameters.");
			else
				flag_time_delay = 1;
		}


		/* read surface time delay params */

		if (strcmp(param, "LOCDELAY_SURFACE") == 0) {
			if ((istat = GetTimeDelaySurface(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Time Delay Surface parameters.");
			else
				flag_time_delay_surface = 1;
		}


		/* read station elevation correction params */

		if (strcmp(param, "LOCELEVCORR") == 0) {
			if ((istat = GetElevCorr(strchr(line, ' '))) < 0)
				puterr("ERROR: reading Elevation Correction parameters.");
			else
				flag_elev_corr = 1;
		}


		/* read source params */

		if (strcmp(param, "LOCSRCE") == 0 || strcmp(param, "GTSRCE") == 0) {
			if ((istat = GetNextSource(strchr(line, ' '))) < 0) {
				puterr("ERROR: reading source params:");
				puterr(line);
			} else
				flag_source = 1;
		}



		/*read transform params */

		if (strcmp(param, "TRANS") == 0) {
    			if ((istat = get_transform(0, strchr(line, ' '))) < 0)
			    puterr("ERROR: reading transformation parameters.");
			else
				flag_trans = 1;
		}



		/* unrecognized input */

		if (istat < 0) {
			if ((pchr = strchr(line, '\n')) != NULL)
				*pchr = '\0';
			sprintf(MsgStr, "Skipping input: %s", line);
			putmsg(5, MsgStr);
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
		putmsg(2, MsgStr);
		Hypocenter.comment[0] = '\0';
	}
	if (!flag_signature) {
		sprintf(MsgStr, "INFO: no signature (LOCSIG) params read.");
		putmsg(2, MsgStr);
		LocSignature[0] = '\0';
	}
	if (!flag_hyptype) {
		sprintf(MsgStr,
"INFO: no hypocenter output file type (LOCHYPOUT) params read.");
		putmsg(2, MsgStr);
		sprintf(MsgStr,
"INFO: DEFAULT: \"LOCHYPOUT SAVE_NLLOC_ALL SAVE_HYPOINV_SUM\"");
		putmsg(2, MsgStr);
		iSaveNLLocEvent = iSaveNLLocSum = 1;
		iSaveHypo71Sum = iSaveHypoEllSum = 0;
		iSaveHypo71Event = iSaveHypoEllEvent = 0;
		iSaveHypoInvSum = 1;
		iSaveAlberto4Sum = 1;
		iSaveSnapSum = 0;
		iSaveDecSec = 0;
	}
	if (!flag_phase_id) {
		sprintf(MsgStr, "INFO: no phase identifier (LOCPHASEID) values read.");
		putmsg(2, MsgStr);
	}
	if (!flag_sta_wt) {
		sprintf(MsgStr, "INFO: no station distance weighting (LOCSTAWT) values read.");
		putmsg(2, MsgStr);
	}
	if (!flag_mag) {
		sprintf(MsgStr, "INFO: no Magnitude Calculation (LOCMAG) params read.");
		putmsg(2, MsgStr);
	}
	if (!flag_phstat) {
		sprintf(MsgStr,
"INFO: no PhaseStatistics (LOCPHSTAT) params read.");
		putmsg(2, MsgStr);
		RMS_Max = VERY_LARGE_DOUBLE;
		NRdgs_Min = -1;
		Gap_Max = VERY_LARGE_DOUBLE;
		P_ResidualMax = VERY_LARGE_DOUBLE;
		S_ResidualMax = VERY_LARGE_DOUBLE;
	}
	if (!flag_comp) {
		sprintf(MsgStr, "INFO: no Component Descirption (LOCCMP) params read.");
		putmsg(2, MsgStr);
	}
	if (!flag_alias) {
		sprintf(MsgStr, "INFO: no Alias (LOCALIAS) params read.");
		putmsg(2, MsgStr);
	}
	if (!flag_exclude) {
		sprintf(MsgStr, "INFO: no Exclude (LOCEXCLUDE) params read.");
		putmsg(2, MsgStr);
	}
	if (!flag_time_delay) {
		sprintf(MsgStr, "INFO: no Time Delay (LOCDELAY) params read.");
		putmsg(2, MsgStr);
	}
	if (!flag_time_delay_surface) {
		sprintf(MsgStr, "INFO: no Time Delay Surface (LOCDELAY_SURFACE) params read.");
		putmsg(2, MsgStr);
	}
	if (!flag_elev_corr) {
		sprintf(MsgStr, "INFO: no Elevation Correction (LOCELEVCORR) params read.");
		putmsg(2, MsgStr);
	}
	if (!flag_otime) {
		sprintf(MsgStr, "INFO: no Fixed Origin Time (LOCFIXOTIME) params read.");
		putmsg(2, MsgStr);
	}
	if (!flag_angles) {
		sprintf(MsgStr, "INFO: no Take-off Angles (LOCANGLES) params read,");
		putmsg(2, MsgStr);
		sprintf(MsgStr, "      default is angleMode=ANGLES_NO, qualtiyMin=5 .");
		putmsg(2, MsgStr);
		angleMode = ANGLE_MODE_NO;
		iAngleQualityMin = 5;
	}
	if (!flag_source) {
		sprintf(MsgStr, "INFO: no Station (LOCSRCE or GTSRCE) params read.");
		putmsg(2, MsgStr);
	}



	return (flag_include * flag_control * flag_outfile * flag_grid *
		flag_search *
		flag_method * flag_gauss * flag_qual2err * flag_trans - 1);
}



/*** function to read output file name ***/

int GetNLLoc_Files(char* line1)
{
	int istat, nObsFile;
	char fnobs[FILENAME_MAX];

	istat = sscanf(line1, "%s %s %s %s %d", fnobs, ftype_obs, fn_loc_grids,
			fn_path_output, &iSwapBytesOnInput);
	if (istat < 5)
	        iSwapBytesOnInput = 0;

	/* check for wildcards in observation file name */
	NumObsFiles = ExpandWildCards(fnobs, fn_loc_obs, MAX_NUM_OBS_FILES);

	sprintf(MsgStr,
		"LOCFILES:  ObsType: %s  InGrids: %s.*  OutPut: %s.* iSwapBytesOnInput: %d",
		 ftype_obs, fn_loc_grids, fn_path_output, iSwapBytesOnInput);
	putmsg(3, MsgStr);
	for (nObsFile = 0; nObsFile < NumObsFiles; nObsFile++) {
		sprintf(MsgStr,
			"   Obs File: %3d  %s", nObsFile, fn_loc_obs[nObsFile]);
			putmsg(3, MsgStr);
		}
	if (NumObsFiles == MAX_NUM_OBS_FILES)
		putmsg(1, "LOCFILES: WARNING: maximum number of files/events reached");

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
		putmsg(3, MsgStr);

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
		putmsg(3, MsgStr);

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
"LOCSEARCH:  Type: %s  init_num_cells_x %d  init_num_cells_y %d  init_num_cells_z %d  min_node_size %lf  max_num_nodes %d  num_scatter %d",
			search_type, octtreeParams.init_num_cells_x, octtreeParams.init_num_cells_y,
			octtreeParams.init_num_cells_z,
			octtreeParams.min_node_size, octtreeParams.max_num_nodes,
			octtreeParams.num_scatter);
		putmsg(3, MsgStr);

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
		/* SH 02/26/2004 added SNAP summary file */
		else if (strcmp(hyp_type, "SAVE_SNAP_SUM") == 0)
			iSaveSnapSum = 1;
		/* filename format, int sec or 5.2 decimal seconds */
		else if (strcmp(hyp_type, "FILENAME_DEC_SEC") == 0)
			iSaveDecSec = 1;
		else
			return(-1);

		strcat(MsgStr, hyp_type);
		strcat(MsgStr, " ");

	} while ((pchr = strchr(pchr + 1, ' ')) != NULL);

	putmsg(3, MsgStr);

	return(0);
}




/*** function to read method ***/

int GetNLLoc_Method(char* line1)
{
	int istat, ierr;

	char loc_method[MAXLINE];


	istat = sscanf(line1, "%s %lf %d %d %d %lf %d %lf", loc_method,
		&DistStaGridMax, &MinNumArrLoc, &MaxNumArrLoc, &MinNumSArrLoc,
		&VpVsRatio, &MaxNum3DGridMemory, &DistStaGridMin);
	if (istat < 8)
		DistStaGridMin = -1.0;

	sprintf(MsgStr,
"LOCMETHOD:  method: %s  minDistStaGrid: %lf  maxDistStaGrid: %lf  minNumberPhases: %d  maxNumberPhases: %d  minNumberSphases: %d  VpVsRatio: %lf  max3DGridMemory: %d",
		loc_method, DistStaGridMin, DistStaGridMax, MinNumArrLoc, MaxNumArrLoc,
		MinNumSArrLoc, VpVsRatio, MaxNum3DGridMemory);
	putmsg(3, MsgStr);

	ierr = 0;

	if (ierr < 0 || istat < 7)
		return(-1);

	if (strcmp(loc_method, "GAU_ANALYTIC") == 0)
		LocMethod = METH_GAU_ANALYTIC;
	else if (strcmp(loc_method, "GAU_TEST") == 0)
		LocMethod = METH_GAU_TEST;
	else if (strcmp(loc_method, "EDT") == 0 || strcmp(loc_method, "EDT_TEST") == 0)
		LocMethod = METH_EDT;
	else if (strcmp(loc_method, "EDT_BOX") == 0)
		LocMethod = METH_EDT_BOX;
	else {
		LocMethod = METH_UNDEF;
		puterr2("ERROR: unrecognized location method:", loc_method);
	}

	if (MaxNumArrLoc < 1)
		MaxNumArrLoc = MAX_NUM_ARRIVALS;

	if (VpVsRatio > 0.0 && GeometryMode == MODE_GLOBAL) {
		puterr("ERROR: cannot use VpVsRatio>0 with TRANSFORM GLOBAL.");
		return(EXIT_ERROR_LOCATE);
	}

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
	putmsg(3, MsgStr);

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
	putmsg(3, MsgStr);

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
	putmsg(3, MsgStr);
	sprintf(MsgStr, "  Phase: %s  PhaseID: <%s>",
			PhaseID[NumPhaseID].phase,
			PhaseID[NumPhaseID].id_string);
	putmsg(3, MsgStr);

	NumPhaseID++;

	return(0);
}




/*** function to read station distance weighting params ***/

int GetStaWeight(char* line1)
{
	int istat, ierr;


	istat = sscanf(line1, "%d %lf", &iSetStationDistributionWeights, &stationDistributionWeightCutoff);

	sprintf(MsgStr, "LOCSTAWT:  flag: %d  CutoffDist: %f",
		iSetStationDistributionWeights, stationDistributionWeightCutoff);
	putmsg(3, MsgStr);

	ierr = 0;
	//if (checkRangeDouble("LOCSTAWT", "Station distribution weight cutoff distance",
	//		stationDistributionWeightCutoff, 1, 0.0, 0, 0.0) != 0)
	//	ierr = -1;

	if (ierr < 0 || istat != 2)
		return(-1);

	return(0);
}



/*** function to read gaussian params ***/

int GetNLLoc_Gaussian(char* line1)
{
	int istat, ierr;


	istat = sscanf(line1, "%lf %lf", &(Gauss.SigmaT), &(Gauss.CorrLen));

	sprintf(MsgStr, "LOCGAUSS:  SigmaT: %lf  CorrLen: %lf",
		Gauss.SigmaT, Gauss.CorrLen);
	putmsg(3, MsgStr);

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
		putmsg(3, MsgStr);

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
		putmsg(3, MsgStr);

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
	putmsg(3, MsgStr);

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
	putmsg(4, MsgStr);

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
	putmsg(3, MsgStr);

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
	putmsg(3, MsgStr);

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

	sscanf(line1, "%s %s",
		LocExclude[NumLocExclude].label, LocExclude[NumLocExclude].phase);

	sprintf(MsgStr, "LOCEXCLUDE:  Name: %s  Phase: %s",
		LocExclude[NumLocExclude].label, LocExclude[NumLocExclude].phase);
	putmsg(3, MsgStr);

	NumLocExclude++;

	return(0);
}




/*** function to read station phase time delays ***/

int GetTimeDelays(char* line1)
{

	if (NumTimeDelays >= MAX_NUM_STA_DELAYS) {
		sprintf(MsgStr, "%s", line1);
		putmsg(3, MsgStr);
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
	putmsg(3, MsgStr);

	NumTimeDelays++;

	return(0);
}





/*** function to read time delay surface (GMT GRD file with x=long, y=lat, z=delay in sec ***/

int GetTimeDelaySurface(char* line1)
{

	sscanf(line1, "%s %lf %s",
		TimeDelaySurfacePhase[NumTimeDelaySurface],
		&TimeDelaySurfaceMultiplier[NumTimeDelaySurface],
		model_surface[NumTimeDelaySurface].grd_file);

	sprintf(MsgStr, "LOCDELAY_SURFACE:  Phase: %s  Mult: %f  GMT GRD File: %s",
		TimeDelaySurfacePhase[NumTimeDelaySurface],
		TimeDelaySurfaceMultiplier[NumTimeDelaySurface],
		model_surface[NumTimeDelaySurface].grd_file);
	putmsg(3, MsgStr);
putmsg(0, MsgStr);

	if (read_grd(&model_surface[NumTimeDelaySurface], message_flag > 2) < 0) {
		puterr2("ERROR: reading Surface Delay GMT GRD File: ",
			model_surface[NumTimeDelaySurface].grd_file);
		return(-1);
	}

	NumTimeDelaySurface++;

	return(0);
}





/*** function to read elevation correction params ***/

int GetElevCorr(char* line1)
{

	int istat;

	istat = sscanf(line1, "%d %lf %lf",
		&ApplyElevCorrFlag, &ElevCorrVelP, &ElevCorrVelS);

	sprintf(MsgStr, "LOCELEVCORR:  Flag: %d  VelP: %lf  VelS: %lf",
		ApplyElevCorrFlag, ElevCorrVelP, ElevCorrVelS);
putmsg(1, MsgStr);
	putmsg(3, MsgStr);

	if (istat != 3)
		return(-1);

	return(0);
}






/*** function to open summary output files */

int OpenSummaryFiles(char *path_output, char* typename)
{

	int ngrid;
	char fname[FILENAME_MAX];



	for (ngrid = 0; ngrid < NumLocGrids; ngrid++) {

		if (!LocGridSave[ngrid])
			continue;

		/* Grid Hyp format */

		pSumFileHypNLLoc[ngrid] = NULL;
		sprintf(fname, "%s.sum.%s%d.loc.hyp", path_output, typename, ngrid);
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
			sprintf(fname, "%s.sum.%s%d.loc.hypo_71", path_output, typename, ngrid);
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
			sprintf(fname, "%s.sum.%s%d.loc.hypo_ell", path_output, typename, ngrid);
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
			sprintf(fname, "%s.sum.%s%d.loc.hypo_inv", path_output, typename, ngrid);
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
			sprintf(fname, "%s.sum.%s%d.loc.sim", path_output, typename, ngrid);
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
			pSumFileHypNLLoc[ngrid] = NULL;
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

		/* Alberto 3D 4 chr sta SIMULPS format */
		if (pSumFileAlberto4[ngrid] != NULL) {
			fclose(pSumFileAlberto4[ngrid]);
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

	return(0);

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

	return(0);

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
	fprintf(fpio, "%5.1lf%5.1lf", erh, erz);

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

	return(0);

}




/*** function to write hypocenter summary to file (HYPOINVERSE archive format) */

int WriteHypoInverseArchive(FILE *fpio, HypoDesc *phypo, ArrivalDesc *parrivals, int narrivals,
		char *filename, int write_header, int write_arrivals)
{

	int ifile = 0, narr;
	char fname[FILENAME_MAX];
	ArrivalDesc *parr, *psarr;
	double resid;

	char first_mot;

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

// AJL 21JUN2000 bug fix (not all phases are written,
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

		if (strpbrk(parr->first_mot, "cCuU+"))
			first_mot = 'U';
		else if (strpbrk(parr->first_mot, "dD-"))
			first_mot = 'D';
		else
			first_mot = ' ';

		fprintf(fpio,
			"%4s%1s%1s%1c%1d%1s",
			parr->label, " ", parr->phase, first_mot,
			parr->quality > 9 ? 9 : parr->quality, " ");
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

	return(0);

}



/*** function to determine azimuth gap -
			largest azimuth separation between
			stations as seen from epicenter (deg) */

double CalcAzimuthGap(ArrivalDesc *arrival, int num_arrivals)
{

	int narr, naz;
	double az_last, az, gap, gap_max = -1.0;
	double azimuth[MAX_NUM_ARRIVALS];

	/* load azimuths to working array */
	naz = 0;
	for (narr = 0; narr < num_arrivals; narr++) {
		if ((arrival + narr)->weight > 0.001)
			azimuth[naz++] = (arrival + narr)->azim;
	}

	/* sort */
	SortDoubles(azimuth, naz);

	/* find largest gap */
	az_last = azimuth[0];
	for (narr = 1; narr < naz; narr++) {
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
			// removing delay so add (Tobs = Tcorr + (O-C))
			p_time[npair] = arrival[narr - 1].obs_time + arrival[narr - 1].delay;
			p_error[npair] = arrival[narr - 1].error;
// DELAY_CORR
			// removing delay so add (Tobs = Tcorr + (O-C))
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

	B = Bmin = 2.0;
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
	putmsg(3, MsgStr);

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
	putmsg(3, MsgStr);


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

	return(0);
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
		if ((pcomp + nComp)->label[0] != ARRIVAL_NULL_CHR &&
				(pcomp + nComp)->label[0] != '*' &&
				strcmp((pcomp + nComp)->label,
				test_label) != 0)
			continue;
		if ((pcomp + nComp)->inst[0] != ARRIVAL_NULL_CHR &&
				(pcomp + nComp)->inst[0] != '*' &&
				strcmp((pcomp + nComp)->inst,
				parr->inst) != 0)
			continue;
		if ((pcomp + nComp)->comp[0] != ARRIVAL_NULL_CHR &&
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
"\n#Average Phase Residuals (CalcResidual)  RMS_Max: %lf  NRdgs_Min: %d  Gap_Max: %lf  P_Res_Max: %lf  S_Res_Max: %lf\n",
			rms_max, nRdgs_min, gap_max,
			p_residual_max, s_residual_max);
		fprintf(fpio,
"#         ID      Phase   Nres      AveRes       StdDev       ResMin       ResMax     ignored\n");
	} else if (imode == WRITE_RES_DELAYS) {
		fprintf(fpio,
"\n#Total Phase Corrections (CalcResidual + InputDelay)  RMS_Max: %lf  NRdgs_Min: %d  Gap_Max: %lf  P_Res_Max: %lf  S_Res_Max: %lf\n",
			rms_max, nRdgs_min, gap_max,
			p_residual_max, s_residual_max);
		fprintf(fpio,
			"#         ID      Phase   Nres      TotCorr\n");
	} else if (imode == WRITE_PDF_RESIDUALS) {
		fprintf(fpio,
"\n#Average Phase Residuals PDF (CalcPDFResidual)  RMS_Max: %lf  NRdgs_Min: %d  Gap_Max: %lf  P_Res_Max: %lf  S_Res_Max: %lf\n",
			rms_max, nRdgs_min, gap_max,
			p_residual_max, s_residual_max);
		fprintf(fpio,
"#         ID      Phase   Nres      AveRes       StdDev       ResMin       ResMax     ignored\n");
	} else if (imode == WRITE_PDF_DELAYS) {
		fprintf(fpio,
"\n#Total Phase Corrections PDF (CalcPDFResidual + InputDelay)  RMS_Max: %lf  NRdgs_Min: %d  Gap_Max: %lf  P_Res_Max: %lf  S_Res_Max: %lf\n",
			rms_max, nRdgs_min, gap_max,
			p_residual_max, s_residual_max);
		fprintf(fpio,
			"#         ID      Phase   Nres      TotCorr\n");
	}

	nnodes = 0;
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
//printf("LOCDELAY  %s %s %d %f = %f + %f/%f\n", np->label, np->phase, np->num_residuals, res_temp,
//np->delay, np->residual_sum, np->weight_sum);
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
		if ((IsPhaseID((arrival + narr)->phase, "P") &&
				fabs((arrival + narr)->residual)
							<= p_residual_max)
			|| (IsPhaseID((arrival + narr)->phase, "S") &&
				fabs((arrival + narr)->residual)
							<= s_residual_max))
			if ((np = InstallStaStatInTable(ntable,
					(arrival + narr)->label,
					(arrival + narr)->phase,
					(arrival + narr)->flag_ignore,
					(arrival + narr)->residual,
					weight,
					(arrival + narr)->pdf_residual_sum,
					(arrival + narr)->pdf_weight_sum,
					(arrival + narr)->delay)) == NULL)
				puterr("ERROR: cannot put arrival statistics in table");

}


/** end of hashtable routines */
/*------------------------------------------------------------/ */






/*** function to add arrival station to station list */

int addToStationList(SourceDesc *stations, int numStations, ArrivalDesc *arrival, int nArrivals) {

	int i, n, nAdded = 0;;


	// set default ignored flag
	for (n = 0; n < MAX_NUM_ARRIVALS; n++) {
		(stations + n)->ignored = 1;
		(stations + n)->station_weight = 1.0;
	}

			
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
			if (!((arrival + i)->flag_ignore))
				(stations + n)->ignored = 0;
			nAdded++;
			numStations++;
			sprintf(MsgStr, "Added to station list: %s (%lf,%lf,%lf)",
				(stations + n)->label, (stations + n)->x, (stations + n)->y,
				(stations + n)->z);
			putmsg(4, MsgStr);
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

	return(0);

}



/*** function to set station distribution weights */

int setStationDistributionWeights(SourceDesc *stations, int numStations, ArrivalDesc *arrival, int nArrivals) {

	int i, nsta, ndist, n, m, nError = 0;
	double x, y, dist, station_weight, weight;
	double dist_sum, dist_ave, cutoff2;
	double station_weight_sum;

	
	// get cutoff distance
	if (stationDistributionWeightCutoff > 0.0) {
		cutoff2 = stationDistributionWeightCutoff * stationDistributionWeightCutoff;
	} else {
		// use average distance
		dist_sum = 0.0;
		ndist= 0;
		for (n = 0; n < numStations; n++) {
			x = (stations + n)->x;
			y = (stations + n)->y;
			if (x == 0.0 && y == 0.0)	// station location not known
				continue;
			for (m = 0; m < numStations; m++) {
				if (m == n || (stations + m)->ignored)
					continue;
				dist_sum += GetEpiDist((stations + m), x, y);
				ndist++;
			}
		}
		if (ndist <= 0)
			return(-1);
		dist_ave = dist_sum / (double) ndist;
		sprintf(MsgStr, "Station Dist Weight:  Ave Station Distance: %lf (N=%d)", dist_ave, ndist);
		putmsg(2, MsgStr);
		cutoff2 = dist_ave * dist_ave;
	}
	
	
		
		// calculate station weights
		
	station_weight_sum = 0.0;
	nsta= 0;
	for (n = 0; n < numStations; n++) {
		station_weight = 0.0;
		x = (stations + n)->x;
		y = (stations + n)->y;
		if (x == 0.0 && y == 0.0)	// station location not known
			continue;
		for (m = 0; m < numStations; m++) {
			if ((stations + m)->ignored)
				continue;
			dist = GetEpiDist((stations + m), x, y);
			weight = exp(-(dist * dist) / cutoff2);	// 0.0 for farthest -> 1.0 for same postion
			station_weight += weight;	// count of number close stations
		}
		nsta++;
		station_weight = 1.0 / station_weight;	// weight inversely prop to number of close stations
		(stations + n)->station_weight = station_weight;
//printf("Station Weight: %s %lf (%lf,%lf,%lf)\n", (stations + n)->label, (stations + n)->station_weight, (stations + n)->x, (stations + n)->y, (stations + n)->z);
		station_weight_sum += station_weight;
	}
	if (nsta > 0) {
		// normalize
		station_weight_sum /= (double) nsta;
		for (n = 0; n < numStations; n++) {
			(stations + n)->station_weight = (stations + n)->station_weight / station_weight_sum;
			sprintf(MsgStr, "Station Dist Weight: %s %lf (%lf,%lf,%lf)", (stations + n)->label, (stations + n)->station_weight, (stations + n)->x, (stations + n)->y, (stations + n)->z);
			putmsg(2, MsgStr);
		}
	}
	

	// set arrival weight
	
	for (i = 0; i < nArrivals; i++) {

		// find station for this arrival in station list
		for (n = 0; n < numStations; n++) {
			if (strcmp((stations + n)->label, (arrival + i)->label) == 0) {
				(arrival + i)->station_weight = (stations + n)->station_weight;
//printf("Arrival Sta Dist Weight: %s %lf (%lf,%lf,%lf)\n", (arrival + i)->label, (arrival + i)->station_weight, (stations + n)->x, (stations + n)->y, (stations + n)->z);
				break;
			}
		}

		// not found (should not occur)
		if (n == numStations) {
			sprintf(MsgStr, "ERROR: setting Station Distance Weights: station not found in station list: %s", 
				(arrival + i)->label);
			putmsg(1, MsgStr);
			nError++;
			continue;
		}
		

	}

	return(nError);


}





/** function to get travel times for all observed arrivals */

int getTravelTimes(ArrivalDesc *arrival, int num_arr_loc, double xval, double yval, double zval)
{
	int nReject;
	int narr, n_compan;
	FILE* fp_grid;
	double yval_grid = 0.0;
	GridDesc* ptgrid;


	/* loop over observed arrivals */

	nReject = 0;
	for (narr = 0; narr < num_arr_loc; narr++) {
		/* check for companion */
		if ((n_compan = arrival[narr].n_companion) >= 0) {
			if ((arrival[narr].pred_travel_time =
					arrival[n_compan].pred_travel_time) < 0.0)
				nReject++;
			arrival[narr].pred_travel_time *= arrival[narr].tfact;
		/* else check grid type */
		} else {
			if (arrival[narr].gdesc.type == GRID_TIME) {
				/* 3D grid */
				if (arrival[narr].gdesc.buffer == NULL) {
					/* read time grid from disk */
					fp_grid = arrival[narr].fpgrid;
				} else {
					/* read time grid from memory buffer */
					fp_grid = NULL;
				}
				if ((arrival[narr].pred_travel_time =
					(double) ReadAbsInterpGrid3d(fp_grid, &(arrival[narr].gdesc),
						xval, yval, zval)) < 0.0)
					nReject++;
			} else {
				/* 2D grid (1D model) */
				yval_grid = GetEpiDist(&(arrival[narr].station), xval, yval);
				if (GeometryMode == MODE_GLOBAL)
					yval_grid *= KM2DEG;
				if (arrival[narr].sheetdesc.buffer == NULL) {
					/* read time grid from disk */
					fp_grid = arrival[narr].fpgrid;
					ptgrid = &(arrival[narr].gdesc);
				} else {
					/* read time grid from memory buffer */
					fp_grid = NULL;
					ptgrid = &(arrival[narr].sheetdesc);
				}
				if ((arrival[narr].pred_travel_time =
						ReadAbsInterpGrid2d(fp_grid, ptgrid, yval_grid, zval)) < 0.0)
					nReject++;
	//printf("getTT:  xval %lf yval %lf yval_grid %lf zval %lf t %lf \n", xval, yval, yval_grid, zval, arrival[narr].pred_travel_time);
	//display_grid_param(&(arrival[narr].sheetdesc));
			}
			arrival[narr].pred_travel_time *= arrival[narr].tfact;
			// apply crustal correction
			if (ApplyCrustElevCorrFlag
					&& GeometryMode == MODE_GLOBAL
					&& arrival[narr].pred_travel_time > 0.0) {
//printf("arrival[narr].pred_travel_time before: %f", arrival[narr].pred_travel_time);
				if (yval_grid > MinDistCrustElevCorr)
					arrival[narr].pred_travel_time +=
						applyCrustElevCorrection(arrival + narr, xval, yval, zval);
//printf(" -> after: %f\n", arrival[narr].pred_travel_time);
			}
			else if (ApplyElevCorrFlag) {
				arrival[narr].pred_travel_time += arrival[narr].elev_corr;
			}
		}
	}

	return(nReject);

}



/** function to apply crustal correction and elevation correction */

// assumes vertical ray (dtdd = 0.0) !!!

double applyCrustElevCorrection(ArrivalDesc* parrival, double xval, double yval, double zval) {

	double correction = 0.0;
	double dtdd = 0.0;
	char cphase;

	if (IsPhaseID(parrival->phase, "P"))
		cphase = 'P';
	else if (IsPhaseID(parrival->phase, "S"))
		cphase = 'S';
	else
		return(0.0);


	// source
	correction = calc_crust_corr(cphase, yval, xval, zval, VERY_LARGE_DOUBLE, dtdd);
	// receiver
	correction +=
		calc_crust_corr(cphase, parrival->station.dlat,
			parrival->station.dlong, 0.0, -1000.0 * parrival->station.depth, dtdd);

	return(correction);

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

	int nSamples, narr, ipos;
	int nInitial;
	int iGridType;
	int nReject;
	int iReject = 0;
	int iBoundary = 0;
	double xval, yval, zval;

	long double value, dlike, value_max = (long double) -VERY_LARGE_DOUBLE;
	double misfit;
	double misfit_min = VERY_LARGE_DOUBLE, misfit_max = -VERY_LARGE_DOUBLE;
	double hypo_dx = -1.0, hypo_dz = -1.0;

	int nScatterSaved;

	int ix, iy, iz;
	double logWtMtrxSum;
	double volume, log_value;
	ResultTreeNode* presult_node;
	OctNode* poct_node;
	OctNode* pparent_oct_node;

	int n_neigh;
	OctNode* neighbor_node;
	Vect3D coords;
	double min_node_size_x;
	double min_node_size_y;
	double min_node_size_z;
	double smallest_node_size_x = -1.0;
	double smallest_node_size_y = -1.0;
	double smallest_node_size_z = -1.0;



	iGridType = GRID_PROB_DENSITY;

	putmsg(4, "");
	putmsg(4, "Calculating solution in Octree...");

	logWtMtrxSum = log(gauss_par->WtMtrxSum);

	// set min node size
	min_node_size_x = min_node_size_y = min_node_size_z = pParams->min_node_size;
	if (GeometryMode == MODE_GLOBAL)
		min_node_size_x = min_node_size_y = pParams->min_node_size * KM2DEG;

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
		if (nReject && GeometryMode != MODE_GLOBAL) {
			sprintf(MsgStr,
"ERROR: octree sample at (%lf,%lf,%lf) is outside of %d travel time grids.",
				xval, yval, zval, nReject);
			puterr(MsgStr);
		}
		/* calc misfit and prob density */
		value = CalcSolutionQuality(num_arr_loc, arrival, gauss_par, iGridType, &misfit, NULL);
		// calculate cell volume * prob density and put cell in resultTree
		volume = poct_node->ds.x * poct_node->ds.y * poct_node->ds.z;
		if (GeometryMode == MODE_GLOBAL) {
			volume *= DEG2KM * DEG2KM;
			volume *= (ERAD - poct_node->center.z) / ERAD;
			volume *= cos(DE2RA * poct_node->center.y);
		}
		log_value =  logWtMtrxSum + value + log(volume);
		poct_node->value = logWtMtrxSum + value;
		resultTreeRoot = addResult(resultTreeRoot, log_value, poct_node);
		nSamples++;
// END - this block must be identical to block $$$ below

		// save node size
		smallest_node_size_x = poct_node->ds.x;
		smallest_node_size_y = poct_node->ds.y;
		smallest_node_size_z = poct_node->ds.z;

		if (message_flag >= 1 && nSamples % 1000 == 0) {
			fprintf(stdout,
				"OctTree num samples = %d / %d\r", nSamples, pParams->max_num_nodes);
			fflush(stdout);
		}

			}
		}
	}
	nInitial = nSamples;


	/* loop over octree nodes */

	nScatterSaved = 0;
	ipos = 0;

	while (nSamples < pParams->max_num_nodes) {

		presult_node = getHighestLeafValue(resultTreeRoot);
//if (nSamples % 100 == 0)
//fprintf(stderr, "%d getHighestLeafValue %lf\n", nSamples, presult_node->value);

if (presult_node == NULL)
fprintf(stderr, "\npresult_node == NULL!!\n");

		pparent_oct_node = presult_node->pnode;

		// subdivide all HighestLeafValue neighbors

		for (n_neigh = 0; n_neigh < 7; n_neigh++) {

			if (n_neigh == 0) {
		  		neighbor_node = pparent_oct_node;
			} else {
				coords.x = pparent_oct_node->center.x;
				coords.y = pparent_oct_node->center.y;
				coords.z = pparent_oct_node->center.z;
				if (n_neigh == 1) {
			  		coords.x = pparent_oct_node->center.x
						+ (pparent_oct_node->ds.x + smallest_node_size_x) / 2.0;
				} else if (n_neigh == 2) {
			  		coords.x = pparent_oct_node->center.x
						- (pparent_oct_node->ds.x + smallest_node_size_x) / 2.0;
				} else if (n_neigh == 3) {
			  		coords.y = pparent_oct_node->center.y
						+ (pparent_oct_node->ds.y + smallest_node_size_y) / 2.0;
				} else if (n_neigh == 4) {
			  		coords.y = pparent_oct_node->center.y
						- (pparent_oct_node->ds.y + smallest_node_size_y) / 2.0;
				} else if (n_neigh == 5) {
			  		coords.z = pparent_oct_node->center.z
						+ (pparent_oct_node->ds.z + smallest_node_size_z) / 2.0;
				} else if (n_neigh == 6) {
			  		coords.z = pparent_oct_node->center.z
						- (pparent_oct_node->ds.z + smallest_node_size_z) / 2.0;
				}
				// check for longitude wrap-around
				if (GeometryMode == MODE_GLOBAL) {
					if (coords.x > 180.0)
						coords.x -= 360.0;
					if (coords.x < -180.0)
						coords.x += 360.0;
				}

				// find neighbor node
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
				if (poct_node->ds.x < smallest_node_size_x)
					smallest_node_size_x = poct_node->ds.x;
				if (poct_node->ds.y < smallest_node_size_y)
					smallest_node_size_y = poct_node->ds.y;
				if (poct_node->ds.z < smallest_node_size_z)
					smallest_node_size_z = poct_node->ds.z;

// $$$ NOTE: this block must be identical to block $$$ above
		// evaluate solution at node
		xval = poct_node->center.x;
		yval = poct_node->center.y;
		zval = poct_node->center.z;
		/* get travel times for observed arrivals */
		nReject = getTravelTimes(arrival, num_arr_loc, xval, yval, zval);
		if (nReject && GeometryMode != MODE_GLOBAL) {
			sprintf(MsgStr,
"ERROR: octree sample at (%lf,%lf,%lf) is outside of %d travel time grids.",
				xval, yval, zval, nReject);
			puterr(MsgStr);
		}
		/* calc misfit and prob density */
		value = CalcSolutionQuality(num_arr_loc, arrival, gauss_par, iGridType, &misfit, NULL);

		// calculate cell volume * prob density and put cell in resultTree
		volume = poct_node->ds.x * poct_node->ds.y * poct_node->ds.z;
		if (GeometryMode == MODE_GLOBAL) {
			volume *= DEG2KM * DEG2KM;
			volume *= (ERAD - poct_node->center.z) / ERAD;
			volume *= cos(DE2RA * poct_node->center.y);
		}
		log_value =  logWtMtrxSum + value + log(volume);
		poct_node->value = logWtMtrxSum + value;
		resultTreeRoot = addResult(resultTreeRoot, log_value, poct_node);
		nSamples++;
// END - this block must be identical to block $$$ above

		if (message_flag >= 1 && nSamples % 1000 == 0) {
			fprintf(stdout,
				"OctTree num samples = %d / %d\r", nSamples, pParams->max_num_nodes);
			fflush(stdout);
		}

				/* check for maximum likelihood */
				if (value > value_max) {
					value_max = value;
					misfit_min = misfit;
					phypo->misfit = misfit;
					phypo->x = xval;
					phypo->y = yval;
					phypo->z = zval;
					hypo_dx = poct_node->ds.x;
					hypo_dz = poct_node->ds.z;
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
	dlike = (long double) gauss_par->WtMtrxSum * (long double) exp(value);
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

		// check if minimum node size reached
		if (smallest_node_size_x < min_node_size_x
				|| smallest_node_size_y < min_node_size_y
				|| smallest_node_size_z < min_node_size_z) {
			if (message_flag >= 1)
				fprintf(stdout, "\nINFO: Min node size reached, terminating Octree search.");
			break;
		}

	}  // end while (nSamples < pParams->max_num_nodes)

	if (message_flag >= 1)
		fprintf(stdout, "\n");


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

	if ((iBoundary = isOnGridBoundary(phypo->x, phypo->y, phypo->z,
			ptgrid, hypo_dx, hypo_dz, 0))) {
		sprintf(MsgStr,
"WARNING: max prob location on grid boundary %d, rejecting location.", iBoundary);
		putmsg(1, MsgStr);
		sprintf(phypo->locStatComm, "%s", MsgStr);
		iReject = 1;
	}


	// construct search information string
	if (GeometryMode == MODE_GLOBAL) {
		smallest_node_size_x *= DEG2KM;
		smallest_node_size_y *= DEG2KM;
	}
	sprintf(phypo->searchInfo, "OCTREE nInitial %d nEvaluated %d smallestNodeSide %lf/%lf/%lf%c",
			nInitial, nSamples, smallest_node_size_x, smallest_node_size_y, smallest_node_size_z, '\0');
	/* write message */
	putmsg(2, phypo->searchInfo);


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
	putmsg(3, "");
	putmsg(3, "Generating event scatter file...");


	integral = integrateResultTree(resultTreeRoot, 0.0, oct_node_value_max);
	sprintf(MsgStr, "Octree oct_node_value_max= %le integral= %le", oct_node_value_max, integral);
	putmsg(2, MsgStr);


	/* generate scatter points at uniformly-randomly chosen locations in each leaf node */

	tot_npoints = 0;
	fdata_index = 0;
	tot_npoints = getScatterSampleResultTree(resultTreeRoot, pParams, integral,
		fscatterdata, tot_npoints, &fdata_index, oct_node_value_max);

	/* write message */
	sprintf(MsgStr, "  %d points generated, %d points requested.",
		tot_npoints, pParams->num_scatter);
	putmsg(2, MsgStr);

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



