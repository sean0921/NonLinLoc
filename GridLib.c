/*
 * Copyright (C) 1999 Anthony Lomax <lomax@faille.unice.fr>
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


/*   GridLib.c

	grid library functions

*/

/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */


/*
	history:

	ver 01    22SEP1997  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#define EXTERN_MODE 1

#include "GridLib.h"


/*** function to set constants */

void SetConstants(void)
{

	Ellipsoid3D EllNULL = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};

	strcpy(package_name, PACKAGE);
	strcpy(prog_ver, PVER);
	strcpy(prog_date, PDATE);
	strcpy(prog_copyright, PCOPYRIGHT);
	message_flag = 0;

	pi = 4. * atan(1.);		/* PI */
	rpd = pi /180.;			/* radians per degree */
	c111 = 10000.0 / 90.0;		/* kilometers per degree */

	fp_include = NULL;		/* set include file ptr */

	/* bookkeeping */
	NumFilesOpen = 0;
	NumGridBufFilesOpen = 0;
	NumGridHdrFilesOpen = 0;
	NumAllocations = 0;

	/* program variables */
	NumQuality2ErrorLevels = 0;

	/* set null angles indicator */
	AnglesNULL = SetTakeOffAngles(400.0, 200.0, 0);

	/* set null elliposid */
	EllipsoidNULL = EllNULL;

}



/*** function to read control params ***/

int get_control(char* line1)
{
	int istat;

	istat = sscanf(line1, "%d", &message_flag, &RandomNumSeed);

	if (istat == 1)
		RandomNumSeed = 837465;

	/* display program information */
	DispProgInfo();

	sprintf(MsgStr, "CONTROL:  MessageFlag: %d  RandomNumSeed: %d",
		message_flag, RandomNumSeed);
	putmsg(3, MsgStr);

	if (checkRangeInt("CONTROL", "MessageFlag", message_flag, 1, 0, 0, 0) != 0)
		return(-1);

	if (istat != 1 && istat != 2)
		return(-1);

	return(0);

}



/*** function to read output file name ***/

int get_outfile(char* line1)
{

	int istat;

	istat = sscanf(line1, "%s", fn_output);

	sprintf(MsgStr, "OUTPUT FILES: %s.*", fn_output);
	putmsg(3, MsgStr);

	if (istat != 1)
		return(-1);

	return(0);
}



/*** function to read include file name and reset input to include file ***/

int GetIncludeFile(char* line1, FILE **fp_io)
{

	/* read include file name */

	sscanf(line1, "%s", fn_include);

	sprintf(MsgStr, "Reading from INCLUDE FILE: %s", fn_include);
	putmsg(3, MsgStr);


	/* open include file */

	if ((fp_include = fopen(fn_include, "r")) == NULL) {
		puterr2("ERROR: opening INCLUDE file", fn_include);
		return(-1);
	}
	NumFilesOpen++;


	/* swap file pointers */

	fp_input_save = *fp_io;
	*fp_io = fp_include;


	return(0);
}


/*** function to reset input from include file to inlput file ***/

void SwapBackIncludeFP(FILE **fp_io)
{


	/* swap file pointers */

	*fp_io = fp_input_save;
	fp_include = NULL;

	sprintf(MsgStr, "Returning from INCLUDE FILE: %s.*", fn_include);
	putmsg(3, MsgStr);

}


/** function to read source params fom input line */

int GetNextSource(char* in_line)
{
	int istat;
	SourceDesc *srce_in;


	/* check number of sources */
	if (NumSources >= MAX_NUM_SOURCES) {
		puterr("ERROR: to many sources, ignoring source.");
		return(0);
	}

	srce_in = Source + NumSources;
	istat = GetSource(in_line, srce_in, NumSources);
	NumSources++;

	return(istat);

}


/** function to read source params fom input line */

int GetSource(char* in_line, SourceDesc *srce_in, int num_sources)
{
	int istat, ierr;
	char chr1, chr2, coord_type[MAXLINE];
	double val1, val1a, val1b, val2, val2a, val2b, val3, val4;
	double sign;
	char label[10 * ARRIVAL_LABEL_LEN];


	/* initialize some source fields */
	srce_in->is_coord_xyz = 0;
	srce_in->is_coord_latlon = 0;
	srce_in->otime = 0.0;


	/* read coordinate type */

	istat = sscanf(in_line, "%*s %s", coord_type);

	/* read coordinate type and coordinates */

	if (strncmp(coord_type, "XYZ", 3) == 0)
	{
		istat = sscanf(in_line, "%s %s %lf %lf %lf %lf",
			label,
			coord_type, &val1, &val2, &val3, &val4);
		if (istat != 6)
			return(-1);
		strncpy(srce_in->label, label, ARRIVAL_LABEL_LEN - 1);
		srce_in->is_coord_xyz = 1;
		srce_in->x = val1;
		srce_in->y = val2;
		srce_in->z = val3 - val4;
		srce_in->is_coord_xyz = 1;
		sprintf(MsgStr,
"SOURCE: %3d  Name: %s  Loc:  type: %s  X(east) %.2lf  Y(north) %.2lf  Z(pos DOWN) %.2lf",
			num_sources, srce_in->label,
			coord_type, srce_in->x, srce_in->y, srce_in->z);
		putmsg(3, MsgStr);
	}
	else if (strcmp(coord_type, "LATLON") == 0)
	{
		istat = sscanf(in_line, "%s %s %lf %lf %lf %lf", label,
			coord_type, &val1, &val2, &val3, &val4);
		strncpy(srce_in->label, label, ARRIVAL_LABEL_LEN - 1);
		srce_in->dlat = val1;
		srce_in->dlong = val2;
		srce_in->depth = val3 - val4;
		srce_in->is_coord_latlon = 1;
		ierr = 0;
		if (checkRangeDouble("SRCE",
				"Lat", srce_in->dlat, 1, -90.0, 1, 90.0) != 0)
			ierr = -1;
		if (checkRangeDouble("SRCE",
				"Long", srce_in->dlong, 1, -180.0, 1, 180.0) != 0)
			ierr = -1;
		sprintf(MsgStr,
"SOURCE:  %d  Name: %s  Loc:  type: %s  Lat %.2lf  Long %.2lf  Depth %.2lf",
			num_sources, srce_in->label,
			coord_type, srce_in->dlat, srce_in->dlong,
			srce_in->depth);
		putmsg(3, MsgStr);
		if (ierr < 0 || istat != 6)
			return(-1);
	}
	else if (strcmp(coord_type, "LATLONDM") == 0)
	{
		istat = sscanf(in_line,
				"%s %s %lf %lf %c %lf %lf %c %lf %lf",
			label, coord_type, &val1, &val1a, &chr1,
			&val2, &val2a, &chr2, &val3, &val4);
		strncpy(srce_in->label, label, ARRIVAL_LABEL_LEN - 1);
		if ((toupper(chr1) != 'N' && toupper(chr1) != 'S')
				|| (toupper(chr2) != 'E' && toupper(chr2) != 'W'))
			return(-1);
		sign = toupper(chr1) == 'N' ? 1.0 : -1.0;
		srce_in->dlat = sign * (val1 + val1a / 60.0);
		sign = toupper(chr2) == 'E' ? 1.0 : -1.0;
		srce_in->dlong = sign * (val2 + val2a / 60.0);
		srce_in->depth = val3 - val4;
		srce_in->is_coord_latlon = 1;
		ierr = 0;
		if (checkRangeDouble("SRCE",
				"Lat", srce_in->dlat, 1, -90.0, 1, 90.0) != 0)
			ierr = -1;
		if (checkRangeDouble("SRCE",
				"Long", srce_in->dlong, 1, -180.0, 1, 180.0) != 0)
			ierr = -1;
		sprintf(MsgStr,
"SOURCE:  %d  Name: %s  Loc:  type: %s  Lat %.2lf  Long %.2lf  Depth %.2lf",
			num_sources, srce_in->label,
			coord_type, srce_in->dlat, srce_in->dlong,
			srce_in->depth);
		putmsg(3, MsgStr);
		if (ierr < 0 || istat != 10)
			return(-1);
	}
	else if (strcmp(coord_type, "LATLONDS") == 0)
	{
		istat = sscanf(in_line,
				"%s %s %lf %lf %lf %c %lf %lf %lf %c %lf %lf",
			label, coord_type,
			&val1, &val1a, &val1b, &chr1,
			&val2, &val2a, &val2b, &chr2, &val3, &val4);
		strncpy(srce_in->label, label, ARRIVAL_LABEL_LEN - 1);
		if ((toupper(chr1) != 'N' && toupper(chr1) != 'S')
				|| (toupper(chr2) != 'E' && toupper(chr2) != 'W'))
			return(-1);
		sign = toupper(chr1) == 'N' ? 1.0 : -1.0;
		srce_in->dlat = sign * (val1 + (val1a + val1b / 60.0) / 60.0);
		sign = toupper(chr2) == 'E' ? 1.0 : -1.0;
		srce_in->dlong = sign * (val2 + (val2a + val2b / 60.0) / 60.0);
		srce_in->depth = val3 - val4;
		srce_in->is_coord_latlon = 1;
		ierr = 0;
		if (checkRangeDouble("SRCE",
				"Lat", srce_in->dlat, 1, -90.0, 1, 90.0) != 0)
			ierr = -1;
		if (checkRangeDouble("SRCE",
				"Long", srce_in->dlong, 1, -180.0, 1, 180.0) != 0)
			ierr = -1;
		sprintf(MsgStr,
"SOURCE:  %d  Name: %s  Loc:  type: %s  Lat %.2lf  Long %.2lf  Depth %.2lf",
			num_sources, srce_in->label,
			coord_type, srce_in->dlat, srce_in->dlong,
			srce_in->depth);
		putmsg(3, MsgStr);
		if (ierr < 0 || istat != 12)
			return(-1);
	}
	else
	{
		puterr("ERROR: unrecognized coordinate type.");
		return(-1);
	}




	return(0);
}



/** function to find source from label */

SourceDesc* FindSource(char* label)
{
	int nsrce;

	for (nsrce = 0; nsrce < NumSources; nsrce++) {
		if(strncmp((Source + nsrce)->label, label, ARRIVAL_LABEL_LEN) == 0)
			return(Source + nsrce);
	}
	return(NULL);
}



/*** function to read map transformation parameters from input line ***/

int get_transform(int n_proj, char* in_line)
{
	int istat, ierr;
	double angle;


	map_itype[n_proj] = MAP_TRANS_UNDEF;
	GeometryMode = MODE_RECT;

	/* read transform input line */

	sscanf(in_line, "%s", map_trans_type[n_proj]);

	if (strcmp(map_trans_type[n_proj], "GLOBAL") == 0) {

		// mode
		GeometryMode = MODE_GLOBAL;

		map_itype[n_proj] = MAP_TRANS_GLOBAL;
		istat = sscanf(in_line, "%s",
			map_trans_type[n_proj]);

		angle = 0.0;
		map_cosang[n_proj] = cos(angle);
		map_sinang[n_proj] = sin(angle);

		sprintf(MapProjStr[n_proj],
			"TRANSFORM  %s",
			map_trans_type[n_proj]);
		putmsg(3, MapProjStr[n_proj]);

		ierr = 0;
		if (ierr < 0 || istat != 1) {
			puterr("ERROR: reading GLOBAL transformation parameters");
			return(-1);
		}


	} else if (strcmp(map_trans_type[n_proj], "SIMPLE") == 0) {

		map_itype[n_proj] = MAP_TRANS_SIMPLE;
		istat = sscanf(in_line, "%s %lf %lf %lf",
			map_trans_type[n_proj], &map_orig_lat[n_proj], &map_orig_long[n_proj],
			&map_rot[n_proj]);

		angle = -rpd * map_rot[n_proj];
		map_cosang[n_proj] = cos(angle);
		map_sinang[n_proj] = sin(angle);

		sprintf(MapProjStr[n_proj],
			"TRANSFORM  %s LatOrig %lf  LongOrig %lf  RotCW %lf",
			map_trans_type[n_proj], map_orig_lat[n_proj], map_orig_long[n_proj], map_rot[n_proj]);
		putmsg(3, MapProjStr[n_proj]);

		ierr = 0;
		if (checkRangeDouble("TRANS",
				"LatOrig", map_orig_lat[n_proj], 1, -90.0, 1, 90.0) != 0)
			ierr = -1;
		if (checkRangeDouble("TRANS",
				"LongOrig", map_orig_long[n_proj], 1, -180.0, 1, 180.0) != 0)
			ierr = -1;
		if (checkRangeDouble("TRANS",
				"RotCW", map_rot[n_proj], 1, -360.0, 1, 360.0) != 0)
			ierr = -1;
		if (ierr < 0 || istat != 4) {
			puterr("ERROR: reading SIMPLE transformation parameters");
			return(-1);
		}



	} else if (strcmp(map_trans_type[n_proj], "LAMBERT") == 0) {

		map_itype[n_proj] = MAP_TRANS_LAMBERT;
		istat = sscanf(in_line, "%s %s %lf %lf %lf %lf %lf",
			map_trans_type[n_proj], map_ref_ellipsoid[n_proj],
			&map_orig_lat[n_proj], &map_orig_long[n_proj],
			&map_lambert_1st_std_paral[n_proj], &map_lambert_2nd_std_paral[n_proj],
			&map_rot[n_proj]);

		ierr = 0;
		if (checkRangeDouble("TRANS",
				"LatOrig", map_orig_lat[n_proj], 1, -90.0, 1, 90.0) != 0)
			ierr = -1;
		if (checkRangeDouble("TRANS",
				"LongOrig", map_orig_long[n_proj], 1, -180.0, 1, 180.0) != 0)
			ierr = -1;
		if (checkRangeDouble("TRANS",
				"FirstStdParal", map_lambert_1st_std_paral[n_proj], 1, -90.0, 1, 90.0) != 0)
			ierr = -1;
		if (checkRangeDouble("TRANS",
				"SecondStdParal", map_lambert_2nd_std_paral[n_proj], 1, -90.0, 1, 90.0) != 0)
			ierr = -1;
		if (checkRangeDouble("TRANS",
				"RotCW", map_rot[n_proj], 1, -360.0, 1, 360.0) != 0)
			ierr = -1;

		angle = -rpd * map_rot[n_proj];
		map_cosang[n_proj] = cos(angle);
		map_sinang[n_proj] = sin(angle);

		/* initialize GMT projection values */
		if (map_setup_proxy(n_proj, map_ref_ellipsoid[n_proj]) < 0) {
			puterr(
"ERROR: initializing general transformation parameters, RefEllipsoid may be invalid");
			return(-1);
		}

		/* initialize lambert projection */
		vlamb(n_proj, map_orig_long[n_proj], map_orig_lat[n_proj],
			map_lambert_1st_std_paral[n_proj], map_lambert_2nd_std_paral[n_proj]);

		sprintf(MapProjStr[n_proj],
"TRANSFORM  %s RefEllipsoid %s  LatOrig %lf  LongOrig %lf  FirstStdParal %lf  SecondStdParal %lf  RotCW %lf",
			map_trans_type[n_proj], map_ref_ellipsoid[n_proj],
			map_orig_lat[n_proj], map_orig_long[n_proj],
			map_lambert_1st_std_paral[n_proj], map_lambert_2nd_std_paral[n_proj],
			map_rot[n_proj]);
		putmsg(3, MapProjStr[n_proj]);

		if (ierr < 0 || istat != 7) {
			puterr("ERROR: reading LAMBERT transformation parameters");
			return(-1);
		}

	} else {

		puterr("ERROR: unrecognized map transformation type");
		return(-1);
	}

    return (0);
}



/** function to read grid params */

int get_grid(char* input_line)
{
	int istat, ierr = 0;

	istat = sscanf(input_line, "%d %d %d %lf %lf %lf %lf %lf %lf %s",
		&(grid_in.numx), &(grid_in.numy), &(grid_in.numz),
		&(grid_in.origx), &(grid_in.origy), &(grid_in.origz),
		&(grid_in.dx), &(grid_in.dy), &(grid_in.dz), grid_in.chr_type);

	grid_in.iSwapBytes = 0;

	convert_grid_type(&grid_in);
	if (message_flag >= 2)
		display_grid_param(&grid_in);

	if (checkRangeInt("GRID", "xNum", grid_in.numx, 1, 2, 0, 0) != 0)
		ierr = -1;
	if (checkRangeInt("GRID", "yNum", grid_in.numy, 1, 2, 0, 0) != 0)
		ierr = -1;
	if (checkRangeInt("GRID", "zNum", grid_in.numz, 1, 2, 0, 0) != 0)
		ierr = -1;

	if (ierr < 0)
		return(-1);

	if (istat != 10)
		return(-1);

	return(0);
}


/** function to convert grid type from char to numeric */

int convert_grid_type(GridDesc* pgrid)
{
	int istat;

	/* check grid type */

	if (strcmp(pgrid->chr_type, "VELOCITY") == 0)
		pgrid->type = GRID_VELOCITY;
	else if (strcmp(pgrid->chr_type, "VELOCITY_METERS") == 0)
		pgrid->type = GRID_VELOCITY_METERS;
	else if (strcmp(pgrid->chr_type, "SLOWNESS") == 0)
		pgrid->type = GRID_SLOWNESS;
	else if (strcmp(pgrid->chr_type, "SLOW_LEN") == 0)
		pgrid->type = GRID_SLOW_LEN;
	else if (strcmp(pgrid->chr_type, "VEL2") == 0)
		pgrid->type = GRID_VEL2;
	else if (strcmp(pgrid->chr_type, "SLOW2") == 0)
		pgrid->type = GRID_SLOW2;
	else if (strcmp(pgrid->chr_type, "SLOW2_METERS") == 0)
		pgrid->type = GRID_SLOW2_METERS;

	else if (strcmp(pgrid->chr_type, "TIME") == 0)
		pgrid->type = GRID_TIME;
	else if (strcmp(pgrid->chr_type, "TIME2D") == 0)
		pgrid->type = GRID_TIME_2D;

	else if (strcmp(pgrid->chr_type, "ANGLE") == 0)
		pgrid->type = GRID_ANGLE;
	else if (strcmp(pgrid->chr_type, "ANGLE2D") == 0)
		pgrid->type = GRID_ANGLE_2D;

	else if (strcmp(pgrid->chr_type, "PROB_DENSITY") == 0)
		pgrid->type = GRID_PROB_DENSITY;
	else if (strcmp(pgrid->chr_type, "MISFIT") == 0)
		pgrid->type = GRID_MISFIT;

	else if (strcmp(pgrid->chr_type, "DEPTH") == 0)
		pgrid->type = GRID_DEPTH;

	else {
		pgrid->type = GRID_UNDEF;
		puterr2("WARNING: unrecognized grid type", pgrid->chr_type);
	}



	return(pgrid->type);
}




/** function to display grid params */

int display_grid_param(GridDesc* pgrid)
{
	int istat;

	fprintf(stdout,
"GRID: {x, y, z}\n  Num: {%d, %d, %d}\n  Orig: {%.3le, %.3le, %.3le}\n  LenSide: {%.4lf, %.4lf, %.4lf}\n",
		pgrid->numx, pgrid->numy, pgrid->numz,
		pgrid->origx, pgrid->origy, pgrid->origz,
		pgrid->dx, pgrid->dy, pgrid->dz);
	fprintf(stdout,"  Type: %s\n", pgrid->chr_type);

	return(0);
}






/*** function to allocate buffer for 3D grid ***/

float* AllocateGrid(GridDesc* pgrid)
{

	pgrid->buffer = (float *) malloc((size_t)
		(pgrid->numx * pgrid->numy * pgrid->numz * sizeof(float)));
	if (pgrid->buffer != NULL)
		NumAllocations++;

	return(pgrid->buffer);
}



/*** function to free buffer for 3D grid ***/

void FreeGrid(GridDesc* pgrid)
{
	if (pgrid->buffer != NULL) {
		free(pgrid->buffer);
		NumAllocations--;
	}
	pgrid->buffer = NULL;
}



/*** function to initialize buffer for 3D grid ***/

int InitializeGrid(GridDesc* pgrid, float init_value)
{

	float *gbuf;

	gbuf = pgrid->buffer + pgrid->numx * pgrid->numy * pgrid->numz;

	while(gbuf-- > pgrid->buffer)
		*gbuf = init_value;

	return(0);
}



/*** function to create array for accessing 3D grid ***/

float*** CreateGridArray(GridDesc* pgrid)
{

	int ix, iy, numyz;
	float ***garray;

	if ((garray = (float ***)
			malloc((size_t) pgrid->numx * sizeof(float **)))
			== NULL)
		return(garray);
	NumAllocations++;

	numyz = pgrid->numy * pgrid->numz;
	for (ix = 0; ix < pgrid->numx; ix++) {
        	if ((garray[ix] = (float **) malloc((size_t) pgrid->numy *
				sizeof(float *))) == NULL)
			return(NULL);
		NumAllocations++;
		for (iy = 0; iy < pgrid->numy; iy++) {
        		garray[ix][iy] = pgrid->buffer +
			ix * numyz + iy * pgrid->numz;
		}
	}

	pgrid->array = garray;


	return(garray);
}


/*** function to free array for accessing 3D grid ***/

void DestroyGridArray(GridDesc* pgrid)
{

	int ix;

	if (pgrid->array != NULL) {

		for (ix = 0; ix < pgrid->numx; ix++) {
        		free(pgrid->array[ix]);
			NumAllocations--;
		}

		free(pgrid->array);
		NumAllocations--;

		pgrid->array = NULL;

	}

}




/*** function to duplicate a 3D grid description and allocate grid memory */

void DuplicateGrid(GridDesc* pnew_grid, GridDesc* pold_grid, char *new_chr_type)
{

	/* copy grid description */
	*pnew_grid = *pold_grid;

	/* set grid type */
	strcpy(pnew_grid->chr_type, new_chr_type);
	convert_grid_type(pnew_grid);


	/* allocate grid */
	pnew_grid->buffer = AllocateGrid(pnew_grid);
	if (pnew_grid->buffer == NULL) {
		puterr(
"ERROR: allocating memory for duplicate 3D grid buffer.");
		exit(EXIT_ERROR_MEMORY);
	}

	/* create grid array access pointers */
	pnew_grid->array = CreateGridArray(pnew_grid);
	if (pnew_grid->array == NULL) {
		puterr(
"ERROR: creating array for accessing duplicate 3D grid buffer.");
		exit(EXIT_ERROR_MEMORY);
	}

}




/*** function to perform several checks on grid data ***/

int CheckGridArray(GridDesc* pgrid, double gridMax, double gridMaxReplace,
		double gridMin,double gridMinReplace)
{

	int ix, iy, iz;
	int ierror = 0, inegative = 0, imax = 0, imin = 0;
	float val;

	for (ix = 0; ix < pgrid->numx; ix++) {
		for (iy = 0; iy < pgrid->numy; iy++) {
			for (iz = 0; iz < pgrid->numz; iz++) {

				/* check for negative values */
        			if ((val = pgrid->array[ix][iy][iz]) < 0.0)
					inegative++;

				/* check for out of range values */
        			if (val > gridMax) {
					val = gridMaxReplace;
					imax++;
				} else if (val < gridMin) {
					val = gridMinReplace;
					imin++;
				}

			}
		}
	}

	if (inegative) {
		sprintf(MsgStr,
			"WARNING: %d negative values in grid.", inegative);
		putmsg(1, MsgStr);
			ierror = -1;
	}
	if (imax) {
		sprintf(MsgStr,
"WARNING: %d values > %e in grid replaced with %e",
				imax, gridMax, gridMaxReplace);
		putmsg(1, MsgStr);
			ierror = -1;
	}
	if (imin) {
		sprintf(MsgStr,
"WARNING: %d values < %e in grid replaced with %e",
				imin, gridMin, gridMinReplace);
		putmsg(1, MsgStr);
			ierror = -1;
	}

	return(ierror);
}




/*** function to sum 2 grids ***/

/* reads from  array if fp_grid_new == NULL */

int SumGrids(GridDesc* pgrid_sum, GridDesc* pgrid_new, FILE* fp_grid_new)
{

	int ix, iy, iz;
	float xval, yval, zval, newval;

	xval = pgrid_sum->origx;
	for (ix = 0; ix < pgrid_sum->numx; ix++) {

		yval = pgrid_sum->origy;
		for (iy = 0; iy < pgrid_sum->numy; iy++) {

			zval = pgrid_sum->origz;
			for (iz = 0; iz < pgrid_sum->numz; iz++) {

				if ((newval =
					ReadAbsInterpGrid3d(fp_grid_new,
						pgrid_new, xval, yval, zval))
							> -LARGE_FLOAT)
					pgrid_sum->array[ix][iy][iz] += newval;

				zval += pgrid_sum->dz;
			}

			yval += pgrid_sum->dy;
		}

		xval += pgrid_sum->dx;
	}


	return(0);
}




/** function to check if a point is on grid boundary */

/* return 1 if within tolerance of boundary, return 0 otherwise */

int isOnGridBoundary(double xloc, double yloc, double zloc, GridDesc* pgrid,
	double tolerance_xy, double tolerance_z, int i_check_top)
{

	if (fabs(xloc - pgrid->origx) <= tolerance_xy)
		return(1);
	if (fabs(xloc - (pgrid->origx + (double) (pgrid->numx - 1) * pgrid->dx))
			<= tolerance_xy)
		return(1);
	if (fabs(yloc - pgrid->origy) <= tolerance_xy)
		return(1);
	if (fabs(yloc - (pgrid->origy + (double) (pgrid->numy - 1) * pgrid->dy))
			<= tolerance_xy)
		return(1);
	if (i_check_top && fabs(zloc - pgrid->origz) <= tolerance_z)
		return(1);
	if (fabs(zloc - (pgrid->origz + (double) (pgrid->numz - 1) * pgrid->dz))
			<= tolerance_z)
		return(1);

	return(0);

}


/** function to check if a point is contained within a grid */

int IsPointInsideGrid(GridDesc* pgrid, double xloc, double yloc, double zloc)
{

	if (xloc < pgrid->origx ||
		xloc > pgrid->origx + (double) (pgrid->numx - 1) * pgrid->dx)
		return(0);

	if (yloc < pgrid->origy ||
		yloc > pgrid->origy + (double) (pgrid->numy - 1) * pgrid->dy)
		return(0);

	if (zloc < pgrid->origz ||
		zloc > pgrid->origz + (double) (pgrid->numz - 1) * pgrid->dz)
		return(0);

	return(1);

}




/** function to check if a grid is entirely contained within another grid */

/* return 1 if pgrid_inside is inside pgrid, return 0 otherwise */
/* if iShiftFlag==1 shifts grid to attempt to get it inside */

int IsGridInside(GridDesc* pgrid_inside, GridDesc* pgrid, int iShiftFlag)
{

	double xmin_in, xmax_in, ymin_in, ymax_in, zmin_in, zmax_in;
	double xmin, xmax, ymin, ymax, zmin, zmax;

	if (pgrid_inside == pgrid)
		return(1);

	xmin_in = pgrid_inside->origx;
	xmax_in = pgrid_inside->origx +
		(double) (pgrid_inside->numx - 1) * pgrid_inside->dx;
	ymin_in = pgrid_inside->origy;
	ymax_in = pgrid_inside->origy +
		(double) (pgrid_inside->numy - 1) * pgrid_inside->dy;
	zmin_in = pgrid_inside->origz;
	zmax_in = pgrid_inside->origz +
		(double) (pgrid_inside->numz - 1) * pgrid_inside->dz;

	xmin = pgrid->origx;
	xmax = pgrid->origx + (double) (pgrid->numx - 1) * pgrid->dx;
	ymin = pgrid->origy;
	ymax = pgrid->origy + (double) (pgrid->numy - 1) * pgrid->dy;
	zmin = pgrid->origz;
	zmax = pgrid->origz + (double) (pgrid->numz - 1) * pgrid->dz;

	if (!iShiftFlag) {
		if (xmin_in < xmin || xmax_in > xmax ||
				ymin_in < ymin || ymax_in > ymax ||
				zmin_in < zmin || zmax_in > zmax)
			return(0);
		else
			return(1);

	} else {
		if (xmin_in < xmin)
			pgrid_inside->origx += xmin - xmin_in;
		else if (xmax_in > xmax)
			pgrid_inside->origx -= xmax_in - xmax;
		if (ymin_in < ymin)
			pgrid_inside->origy += ymin - ymin_in;
		else if (ymax_in > ymax)
			pgrid_inside->origy -= ymax_in - ymax;
		if (zmin_in < zmin)
			pgrid_inside->origz += zmin - zmin_in;
		else if (zmax_in > zmax)
			pgrid_inside->origz -= zmax_in - zmax;
		return(IsGridInside(pgrid_inside, pgrid, 0));
	}

}




/** function to check if greatest distance from a station to a 3D grid
		is less than the horizontal extent of a 2D grid and
		if depth range of the 3D grid is less than that of the 2D grid
		and if station is within dist_horiz_max of grid location
		xcent, ycent */

/* return 1 if yes, return 0 otherwise */

int IsGrid2DBigEnough(GridDesc* pgrid_3D, GridDesc* pgrid_2D,
		SourceDesc* station,
		double dist_horiz_max, double xcent, double ycent)
{

	double extent_horiz, extent_vert, distance;

	double xmin_in, xmax_in, ymin_in, ymax_in, zmin_in, zmax_in;
	double ymin, ymax, zmin, zmax;


	// GLOBAL
	if (GeometryMode == MODE_GLOBAL) {
		return(1);	//INGV do not worry about grid size
	}


	/* check distance from grid location xcent, ycent */
	if (dist_horiz_max > SMALL_DOUBLE) {
		if ((distance = GetEpiDist(station, xcent, ycent))
				> dist_horiz_max)
			return(-2);
	}

	/* get extent of 2D grid */
	ymin = pgrid_2D->origy;
	ymax = pgrid_2D->origy + (double) (pgrid_2D->numy - 1) * pgrid_2D->dy;
	extent_horiz = ymax - ymin;
	zmin = pgrid_2D->origz;
	zmax = pgrid_2D->origz + (double) (pgrid_2D->numz - 1) * pgrid_2D->dz;

	/* get bounds of 3D grid */
	xmin_in = pgrid_3D->origx;
	xmax_in = pgrid_3D->origx +
		(double) (pgrid_3D->numx - 1) * pgrid_3D->dx;
	ymin_in = pgrid_3D->origy;
	ymax_in = pgrid_3D->origy +
		(double) (pgrid_3D->numy - 1) * pgrid_3D->dy;
	zmin_in = pgrid_3D->origz;
	zmax_in = pgrid_3D->origz +
		(double) (pgrid_3D->numz - 1) * pgrid_3D->dz;

	// adjust for GLOBAL (grid units = deg, GetEpiDist units = km)
	if (GeometryMode == MODE_GLOBAL)
		extent_horiz /= KM2DEG;

	/* check each horizontal corner */
	if ((distance = GetEpiDist(station, xmin_in, ymin_in)) > extent_horiz)
		return(-1);
	if ((distance = GetEpiDist(station, xmin_in, ymax_in)) > extent_horiz)
		return(-1);
	if ((distance = GetEpiDist(station, xmax_in, ymax_in)) > extent_horiz)
		return(-1);
	if ((distance = GetEpiDist(station, xmax_in, ymin_in)) > extent_horiz)
		return(-1);
	/* check depth range */
	if (zmin_in < zmin || zmax_in > zmax)
		return(-3);


	return(1);
}



/** function to calculate epicentral (horizontal) distance from a
	source to an x-y poistion */

INLINE double GetEpiDist(SourceDesc* psrce, double xval, double yval)
{
	double xtmp, ytmp;

	if (GeometryMode == MODE_GLOBAL) {
		return(GCDistance(yval, xval, psrce->y, psrce->x));
	} else {
		xtmp = xval - psrce->x;
		ytmp = yval - psrce->y;
		return(sqrt(xtmp * xtmp + ytmp * ytmp));
	}
}



/** function to calculate epicentral (horizontal) distance from a
	station to an x-y poistion */

INLINE double GetEpiDistSta(StationDesc* psta, double xval, double yval)
{
	double xtmp, ytmp;

	if (GeometryMode == MODE_GLOBAL) {
		return(GCDistance(yval, xval, psta->y, psta->x));
	} else {
		xtmp = xval - psta->x;
		ytmp = yval - psta->y;
		return(sqrt(xtmp * xtmp + ytmp * ytmp));
	}
}



/** function to calculate epicentral (horizontal) azimuth from a
	x-y position to a source */

double GetEpiAzim(SourceDesc* psrce, double xval, double yval)
{
	double xtmp, ytmp, azim;

	if (GeometryMode == MODE_GLOBAL) {
		return(GCAzimuth(yval, xval, psrce->y, psrce->x));
	} else {
		xtmp = psrce->x - xval;
		ytmp = psrce->y - yval;
		azim = atan2(xtmp, ytmp) / rpd;
		if (azim < 0.0)
			azim += 360.0;
		return(azim);
	}
}


/** function to calculate epicentral (horizontal) azimuth from a
	x-y position to a station */

double GetEpiAzimSta(StationDesc* psta, double xval, double yval)
{
	double xtmp, ytmp, azim;

	if (GeometryMode == MODE_GLOBAL) {
		return(GCAzimuth(yval, xval, psta->y, psta->x));
	} else {
		xtmp = psta->x - xval;
		ytmp = psta->y - yval;
		azim = atan2(xtmp, ytmp) / rpd;
		if (azim < 0.0)
			azim += 360.0;
		return(azim);
	}
}



/** function to calculate distance between 2 XYZ points */

INLINE double Dist3D(double x1, double x2, double y1, double y2,
							double z1, double z2)
{
	double dx, dy, dz;

	dx = x1 - x2;
	dy = y1 - y2;
	dz = z1 - z2;
	return(sqrt(dx * dx + dy * dy + dz * dz));
}




/*** function to write grid buffer and header to disk ***/

int WriteGrid3dBuf(GridDesc* pgrid, SourceDesc* psrce,
			char* filename, char* file_type)
{

	int istat;
	FILE *fpio;
	char fname[FILENAME_MAX];


	/* write buffer file */

	sprintf(fname, "%s.%s.buf", filename, file_type);
	if ((fpio = fopen(fname, "w")) == NULL) {
		puterr("ERROR: opening buffer output file.");
		return(-1);
	}
	NumFilesOpen++;

	if (fwrite((char *) pgrid->buffer,
		(pgrid->numx * pgrid->numy * pgrid->numz * sizeof(float)),
		1, fpio) != 1) {
		puterr("ERROR: writing grid output file.");
		return(-1);
	}

	fclose(fpio);
	NumFilesOpen--;


	/* write header file */

	istat = WriteGrid3dHdr(pgrid, psrce, filename, file_type);

	return(istat);
}




/*** function to write grid header to disk ***/

int WriteGrid3dHdr(GridDesc* pgrid, SourceDesc* psrce,
			char* filename, char* file_type)
{

	FILE *fpio;
	char fname[FILENAME_MAX];


	/* write header file */

	if (file_type != NULL)
		sprintf(fname, "%s.%s.hdr", filename, file_type);
	else
		sprintf(fname, "%s.hdr", filename);

	if ((fpio = fopen(fname, "w")) == NULL) {
		puterr("ERROR: opening grid output header file.");
		return(-1);
	}
	NumFilesOpen++;

	fprintf(fpio, "%d %d %d  %lf %lf %lf  %lf %lf %lf %s\n",
		pgrid->numx, pgrid->numy, pgrid->numz,
		pgrid->origx, pgrid->origy, pgrid->origz,
		pgrid->dx, pgrid->dy, pgrid->dz, pgrid->chr_type);

	if (pgrid->type == GRID_TIME || pgrid->type == GRID_TIME_2D
		|| pgrid->type == GRID_ANGLE || pgrid->type == GRID_ANGLE_2D)
		fprintf(fpio, "%s %lf %lf %lf\n",
			psrce->label, psrce->x, psrce->y, psrce->z);

	fprintf(fpio, "%s\n", MapProjStr[0]);


	fclose(fpio);
	NumFilesOpen--;

	return(0);
}




/*** function to swap bytes in grid buffer ***/

int swapBytes(float *buffer, long bufsize)
{

	float *bufpos;
	char ctmp;
	unsigned short itmp;

	union
	{
		char cval[4];
		unsigned short ival[2];
		float fval;
	}
	bufval;


	bufpos = buffer;
	while (bufpos < buffer + bufsize) {
		bufval.fval = *bufpos;
		ctmp = bufval.cval[0];
		bufval.cval[0] = bufval.cval[1];
		bufval.cval[1] = ctmp;
		ctmp = bufval.cval[2];
		bufval.cval[2] = bufval.cval[3];
		bufval.cval[3] = ctmp;
		itmp = bufval.ival[0];
		bufval.ival[0] = bufval.ival[1];
		bufval.ival[1] = itmp;
		*bufpos = bufval.fval;
		bufpos++;
	}

	return(0);
}



/*** function to read entire grid buffer from disk ***/

int ReadGrid3dBuf(GridDesc* pgrid, FILE* fpio)
{

	long readsize;



	readsize = pgrid->numx * pgrid->numy * pgrid->numz * sizeof(float);

	/* read from grid file to buffer */

	if (fread((char *) pgrid->buffer, readsize, 1, fpio) != 1) {
		puterr2("ERROR: reading grid file", pgrid->title);
		return(-1);
	}

	if (pgrid->iSwapBytes)
	        swapBytes(pgrid->buffer, readsize / sizeof(float));


	return(0);
}



/*** function to read y-z sheet of grid buffer from disk ***/

int ReadGrid3dBufSheet(float* sheetbuf, GridDesc* pgrid_disk,
				FILE* fpio, int ix)
{

	long offset;
	long readsize;


	/* check indexes in range */

	if (ix < 0 || ix >= pgrid_disk->numx) {
		sprintf(MsgStr,
			"WARNING: grid file x-sheet index %d out of range (%d,%d)",
			ix, 0, pgrid_disk->numx - 1);
		return(-1);
	}


	/* calculate offset in bytes */

	offset = sizeof(float) * (ix * (pgrid_disk->numy * pgrid_disk->numz));
	fseek(fpio, offset, SEEK_SET);


	/* calculate number of bytes to read */

	readsize = pgrid_disk->numy * pgrid_disk->numz * sizeof(float);


	/* read sheet from grid file to buffer */

	if (fread((char *) sheetbuf, readsize, 1, fpio) != 1) {
		puterr("ERROR: reading x-sheet grid file.");
		return(-1);
	}

	if (pgrid_disk->iSwapBytes)
	        swapBytes(sheetbuf, readsize / sizeof(float));


	return(0);
}



/*** function to read grid header file ***/

int ReadGrid3dHdr(GridDesc* pgrid, SourceDesc* psrce,
			char* filename, char* file_type)
{

	FILE *fpio;
	char fname[FILENAME_MAX];



	/* read header file */

	sprintf(fname, "%s.%s.hdr", filename, file_type);
	if ((fpio = fopen(fname, "r")) == NULL) {
		if (message_flag >= 1)
			puterr("ERROR: opening grid header file.");
		return(-1);
	}
	NumFilesOpen++;

	fscanf(fpio, "%d %d %d  %lf %lf %lf  %lf %lf %lf %s\n",
		&(pgrid->numx), &(pgrid->numy), &(pgrid->numz),
		&(pgrid->origx), &(pgrid->origy), &(pgrid->origz),
		&(pgrid->dx), &(pgrid->dy), &(pgrid->dz), &(pgrid->chr_type));

	if (strncmp(file_type, "time", 4) == 0
			|| strncmp(file_type, "angle", 4) == 0)
		fscanf(fpio, "%s %lf %lf %lf\n",
			psrce->label, &(psrce->x), &(psrce->y), &(psrce->z));


	fclose(fpio);
	NumFilesOpen--;

	return(0);
}




/*** function to open grid file and read header ***/

int OpenGrid3dFile(char *fname, FILE **fp_grid, FILE **fp_hdr,
		GridDesc* pgrid, char* file_type, SourceDesc* psrce, int iSwapBytes)
{

	char fn_grid[FILENAME_MAX], fn_hdr[FILENAME_MAX];

	/* open grid file and header file */

	sprintf(fn_grid, "%s.buf", fname);
	sprintf(MsgStr, "Opening Grid File: %s", fn_grid);
	putmsg(3, MsgStr);
	if ((*fp_grid = fopen(fn_grid, "r")) == NULL) {
		sprintf(MsgStr,
			    "WARNING: cannot open grid buffer file", fn_grid);
		putmsg(3, MsgStr);
		//return(-1);	// sometimes only header is wanted
	} else {
		NumGridBufFilesOpen++;
		NumFilesOpen++;
	}
	sprintf(fn_hdr, "%s.hdr", fname);
	if ((*fp_hdr = fopen(fn_hdr, "r")) == NULL) {
		sprintf(MsgStr,
			    "WARNING: cannot open grid header file", fn_grid);
		putmsg(3, MsgStr);
		if (*fp_grid != NULL) {
			fclose(*fp_grid);
			NumGridBufFilesOpen--;
			NumFilesOpen--;
		}
		return(-1);
	}
	NumGridHdrFilesOpen++;
	NumFilesOpen++;


	/* read header file */

	pgrid->iSwapBytes = iSwapBytes;
	fscanf(*fp_hdr, "%d %d %d  %lf %lf %lf  %lf %lf %lf  %s\n",
		&(pgrid->numx), &(pgrid->numy), &(pgrid->numz),
		&(pgrid->origx), &(pgrid->origy), &(pgrid->origz),
		&(pgrid->dx), &(pgrid->dy), &(pgrid->dz), pgrid->chr_type);

	// make sure that dx for 2D grids is non-zero
	if (pgrid->numx == 1)
		pgrid->dx = 1.0;


	convert_grid_type(pgrid);
	if (message_flag >= 4)
		display_grid_param(pgrid);

	if (psrce != NULL && (strncmp(file_type, "time", 4) == 0
			|| strncmp(file_type, "angle", 4) == 0))
		fscanf(*fp_hdr, "%s %lf %lf %lf\n",
			psrce->label, &(psrce->x), &(psrce->y), &(psrce->z));

	// save filename as grid identifier
	strcpy(pgrid->title, fname);

	return(0);

}


/*** function to close grid file and header ***/

void CloseGrid3dFile(FILE **fp_grid, FILE **fp_hdr)
{

	if (*fp_grid != NULL) {
		fclose(*fp_grid);
		*fp_grid = NULL;
		NumGridBufFilesOpen--;
		NumFilesOpen--;
	}

	if (*fp_hdr != NULL) {
		fclose(*fp_hdr);
		*fp_hdr = NULL;
		NumGridHdrFilesOpen--;
		NumFilesOpen--;
	}

}


/*** function to read grid data from disk or array at index location ***/

INLINE float ReadGrid3dValue(FILE *fpgrid, int ix, int iy, int iz,
				GridDesc* pgrid)
{
	int numyz;
	long offset;
	float fvalue;

	/* check indexes in range */

	if (ix < 0 || ix >= pgrid->numx || iy < 0 || iy >= pgrid->numy
			|| iz < 0 || iz >= pgrid->numz) {
		//puterr("WARNING: grid file index out of range.");
		return(-VERY_LARGE_FLOAT);
	}

	/* get fvalue */

	if (fpgrid != NULL) {
		/* calculate offset in bytes */
		numyz = pgrid->numy * pgrid->numz;
		offset = sizeof(float) * (ix * numyz + iy * pgrid->numz + iz);
		fseek(fpgrid, offset, SEEK_SET);
		/* read fvalue */
		if (fread(&fvalue, sizeof(float), 1, fpgrid) != 1) {
			puterr2("ERROR: reading grid file", pgrid->title);
			return(-VERY_LARGE_FLOAT);
		}
		if (pgrid->iSwapBytes)
	                swapBytes(&fvalue, 1);
	} else
		fvalue = pgrid->array[ix][iy][iz];


	return(fvalue);
}


/*** function to read grid data from disk or array at absolute location ***/

INLINE float ReadAbsGrid3dValue(FILE *fpgrid, GridDesc* pgrid,
		double xloc, double yloc, double zloc, int ifloor)
{
	int ix, iy, iz;
	float fvalue;
	double shift;


	if (ifloor)
		shift = 0.0;
	else
		shift = 0.5;


	/* calculate nearest grid location */

	ix = (int) (shift + (xloc - pgrid->origx) / pgrid->dx);
	iy = (int) (shift + (yloc - pgrid->origy) / pgrid->dy);
	iz = (int) (shift + (zloc - pgrid->origz) / pgrid->dz);


	/* read grid data */

	fvalue = ReadGrid3dValue(fpgrid, ix, iy, iz, pgrid);


	return(fvalue);
}


/*** function to read grid data from disk or buffer at absolute location
		with interpolation ***/

/* read from file if fpgrid != NULL, otherwise read from grid buffer */

INLINE float ReadAbsInterpGrid3d(FILE *fpgrid, GridDesc* pgrid,
		double xloc, double yloc, double zloc)
{
	int nx, ny, nz;
	int ix0, ix1, iy0, iy1, iz0, iz1;
	float value;
	DOUBLE 	vval000, vval001, vval010, vval011,
		vval100, vval101, vval110, vval111;

	DOUBLE xoff, yoff, zoff;
	DOUBLE xdiff, ydiff, zdiff;


	/* calculate grid locations on edge of solid containing point */

	xoff = (xloc - pgrid->origx) / pgrid->dx;
	ix0 = (int) (xoff - VERY_SMALL_DOUBLE);
	yoff = (yloc - pgrid->origy) / pgrid->dy;
	iy0 = (int) (yoff - VERY_SMALL_DOUBLE);
	zoff = (zloc - pgrid->origz) / pgrid->dz;
	iz0 = (int) (zoff - VERY_SMALL_DOUBLE);

	ix1 = (ix0 < pgrid->numx - 1) ? ix0 + 1 : ix0;
	iy1 = (iy0 < pgrid->numy - 1) ? iy0 + 1 : iy0;
	iz1 = (iz0 < pgrid->numz - 1) ? iz0 + 1 : iz0;

/*	if (ix1 < 0 || ix1 >= pgrid->numx) fprintf(stderr, "GRID INDEX ERROR: ix1 %d (%d-%d) xloc %lf ix0 %d ix1 %d\n", ix1, 0, pgrid->numx, xloc, ix0, ix1);
	if (iy1 < 0 || iy1 >= pgrid->numy) fprintf(stderr, "GRID INDEX ERROR: iy1 %d (%d-%d) yloc %lf iy0 %d iy1 %d\n", iy1, 0, pgrid->numy, yloc, iy0, iy1);
	if (iz1 < 0 || iz1 >= pgrid->numz) fprintf(stderr, "GRID INDEX ERROR: iz1 %d (%d-%d) zloc %lf iz0 %d iz1 %d\n", iz1, 0, pgrid->numz, zloc, iz0, iz1);
	if (ix0 < 0 || ix0 >= pgrid->numx) fprintf(stderr, "GRID INDEX ERROR: ix0 %d (%d-%d) xloc %lf ix0 %d ix1 %d\n", ix0, 0, pgrid->numx, xloc, ix0, ix1);
	if (iy0 < 0 || iy0 >= pgrid->numy) fprintf(stderr, "GRID INDEX ERROR: iy0 %d (%d-%d) yloc %lf iy0 %d iy1 %d\n", iy0, 0, pgrid->numy, yloc, iy0, iy1);
	if (iz0 < 0 || iz0 >= pgrid->numz) fprintf(stderr, "GRID INDEX ERROR: iz0 %d (%d-%d) zloc %lf iz0 %d iz1 %d\n", iz0, 0, pgrid->numz, zloc, iz0, iz1);
*/

/*	if (ix1 < 0 || ix1 >= pgrid->numx) ix1 = ix0;
	if (iy1 < 0 || iy1 >= pgrid->numy) iy1 = iy0;
	if (iz1 < 0 || iz1 >= pgrid->numz) iz1 = iz0;

	if (ix0 < 0 || ix0 >= pgrid->numx) ix0 = ix1;
	if (iy0 < 0 || iy0 >= pgrid->numy) iy0 = iy1;
	if (iz0 < 0 || iz0 >= pgrid->numz) iz0 = iz1;
*/

	xdiff = xoff - (DOUBLE) ix0;
	ydiff = yoff - (DOUBLE) iy0;
	zdiff = zoff - (DOUBLE) iz0;

if (xdiff < 0.0 || xdiff > 1.0)
	return(-VERY_LARGE_DOUBLE);
if (ydiff < 0.0 || ydiff > 1.0)
	return(-VERY_LARGE_DOUBLE);
if (zdiff < 0.0 || zdiff > 1.0)
	return(-VERY_LARGE_DOUBLE);

	/* location at grid node */

	if (xdiff + ydiff + zdiff < SMALL_FLOAT) {
		if (fpgrid != NULL)
			value = ReadGrid3dValue(fpgrid,
				ix0, iy0, iz0, pgrid);
		else
			value = pgrid->array[ix0][iy0][iz0];
		return(value);
	}


	/* read vertex values from grid file or array */

	if (fpgrid != NULL) {
		vval000 = ReadGrid3dValue(fpgrid, ix0, iy0, iz0, pgrid);
		vval001 = ReadGrid3dValue(fpgrid, ix0, iy0, iz1, pgrid);
		vval010 = ReadGrid3dValue(fpgrid, ix0, iy1, iz0, pgrid);
		vval011 = ReadGrid3dValue(fpgrid, ix0, iy1, iz1, pgrid);
		vval100 = ReadGrid3dValue(fpgrid, ix1, iy0, iz0, pgrid);
		vval101 = ReadGrid3dValue(fpgrid, ix1, iy0, iz1, pgrid);
		vval110 = ReadGrid3dValue(fpgrid, ix1, iy1, iz0, pgrid);
		vval111 = ReadGrid3dValue(fpgrid, ix1, iy1, iz1, pgrid);
	} else {
		vval000 = pgrid->array[ix0][iy0][iz0];
		vval001 = pgrid->array[ix0][iy0][iz1];
		vval010 = pgrid->array[ix0][iy1][iz0];
		vval011 = pgrid->array[ix0][iy1][iz1];
		vval100 = pgrid->array[ix1][iy0][iz0];
		vval101 = pgrid->array[ix1][iy0][iz1];
		vval110 = pgrid->array[ix1][iy1][iz0];
		vval111 = pgrid->array[ix1][iy1][iz1];
	}

	/* interpolate values */

	if (pgrid->type == GRID_ANGLE || pgrid->type == GRID_ANGLE_2D) {
		value = InterpCubeAngles(xdiff, ydiff, zdiff,
			vval000, vval001, vval010, vval011,
			vval100, vval101, vval110, vval111);
	} else {
// INGV
// check for invalid / mask nodes
if (vval000 < 0.0 || vval010 < 0.0 || vval100 < 0.0 || vval110 < 0.0
 || vval001 < 0.0 || vval011 < 0.0 || vval101 < 0.0 || vval111 < 0.0)
	return(-VERY_LARGE_DOUBLE);
		value = InterpCubeLagrange(xdiff, ydiff, zdiff,
			vval000, vval001, vval010, vval011,
			vval100, vval101, vval110, vval111);
	}


	return(value);
}


/*** function to find value inside a cube using Lagrange interpolation***/

/* 	0.0 <= vvalKLM <= 1.0 */
/*	returns interp value at (xdiff, ydiff, zdiff) */

INLINE DOUBLE InterpCubeLagrange(DOUBLE xdiff, DOUBLE ydiff, DOUBLE zdiff,
		DOUBLE vval000, DOUBLE vval001, DOUBLE vval010, DOUBLE vval011,
		DOUBLE vval100, DOUBLE vval101, DOUBLE vval110, DOUBLE vval111)
{

	DOUBLE value;

	value =   vval000 * (1.0 - xdiff) * (1.0 - ydiff)  * (1.0 - zdiff)
		+ vval001 * (1.0 - xdiff) * (1.0 - ydiff)  * zdiff
		+ vval010 * (1.0 - xdiff) * ydiff          * (1.0 - zdiff)
		+ vval011 * (1.0 - xdiff) * ydiff          * zdiff
		+ vval100 * xdiff         * (1.0 - ydiff)  * (1.0 - zdiff)
		+ vval101 * xdiff         * (1.0 - ydiff)  * zdiff
		+ vval110 * xdiff         * ydiff          * (1.0 - zdiff)
		+ vval111 * xdiff         * ydiff          * zdiff;

	return (value);

}



/*** function to find angles inside a cube */

INLINE float InterpCubeAngles(DOUBLE xdiff, DOUBLE ydiff, DOUBLE zdiff,
		DOUBLE vval000, DOUBLE vval001, DOUBLE vval010, DOUBLE vval011,
		DOUBLE vval100, DOUBLE vval101, DOUBLE vval110, DOUBLE vval111)
{

	int nx, ny, nz;
	double azim[2][2][2], dip[2][2][2], azim_interp, dip_interp;
	int iqual[2][2][2], iqual_interp, iqual_low;
	float value;
	TakeOffAngles angles;


	/* decode angles on vertices of cube */

	SetAnglesFloat(&angles, vval000);
	iqual[0][0][0] = GetTakeOffAngles(&angles,
		&(azim[0][0][0]), &(dip[0][0][0]), &(iqual[0][0][0]));
	SetAnglesFloat(&angles, vval001);
	iqual[0][0][1] = GetTakeOffAngles(&angles,
		&(azim[0][0][1]), &(dip[0][0][1]), &(iqual[0][0][1]));
	SetAnglesFloat(&angles, vval010);
	iqual[0][1][0] = GetTakeOffAngles(&angles,
		&(azim[0][1][0]), &(dip[0][1][0]), &(iqual[0][1][0]));
	SetAnglesFloat(&angles, vval011);
	iqual[0][1][1] = GetTakeOffAngles(&angles,
		&(azim[0][1][1]), &(dip[0][1][1]), &(iqual[0][1][1]));
	SetAnglesFloat(&angles, vval100);
	iqual[1][0][0] = GetTakeOffAngles(&angles,
		&(azim[1][0][0]), &(dip[1][0][0]), &(iqual[1][0][0]));
	SetAnglesFloat(&angles, vval101);
	iqual[1][0][1] = GetTakeOffAngles(&angles,
		&(azim[1][0][1]), &(dip[1][0][1]), &(iqual[1][0][1]));
	SetAnglesFloat(&angles, vval110);
	iqual[1][1][0] = GetTakeOffAngles(&angles,
		&(azim[1][1][0]), &(dip[1][1][0]), &(iqual[1][1][0]));
	SetAnglesFloat(&angles, vval111);
	iqual[1][1][1] = GetTakeOffAngles(&angles,
		&(azim[1][1][1]), &(dip[1][1][1]), &(iqual[1][1][1]));


	/* check for lowest quality angles */

	iqual_low = 999;
	for (nx = 0; nx < 2 ; nx++)
		for (ny = 0; ny < 2 ; ny++)
			for (nz = 0; nz < 2 ; nz++)
				if (iqual[nx][ny][nz] < iqual_low)
					iqual_low = iqual[nx][ny][nz];


	/* determine angles to return */

	if (iqual_low < ANGLE_QUALITY_CUTOFF) {
		/* if lowest quality is too low, use nearest node */
		value = vval000;
	} else {
		/* otherwise interpolate */
		azim_interp = InterpCubeLagrange(xdiff, ydiff, zdiff,
			azim[0][0][0], azim[0][0][1], azim[0][1][0],
			azim[0][1][1], azim[1][0][0], azim[1][0][1],
			azim[1][1][0], azim[1][1][1]);
		dip_interp = InterpCubeLagrange(xdiff, ydiff, zdiff,
			dip[0][0][0], dip[0][0][1], dip[0][1][0],
			dip[0][1][1], dip[1][0][0], dip[1][0][1],
			dip[1][1][0], dip[1][1][1]);
		angles = SetTakeOffAngles(azim_interp, dip_interp, iqual_low);
		value = angles.fval;
	}


	return (value);

}




/*** function to read grid data from disk or buffer at absolute location
		with interpolation ***/

/* 2D version - ix assumed = 0 */

/* read from file if fpgrid != NULL, otherwise read from grid buffer */

INLINE DOUBLE ReadAbsInterpGrid2d(FILE *fpgrid, GridDesc* pgrid,
		double yloc, double zloc)
{
	int nx, ny, nz;
	int ix0, ix1, iy0, iy1, iz0, iz1;
	DOUBLE value;
	DOUBLE vval00, vval01, vval10, vval11;

	DOUBLE yoff, zoff;
	DOUBLE ydiff, zdiff;

	/* calculate grid locations on edge of solid containing point */

	ix0 = 0;
	yoff = (yloc - pgrid->origy) / pgrid->dy;
	iy0 = (int) (yoff - VERY_SMALL_DOUBLE);
	zoff = (zloc - pgrid->origz) / pgrid->dz;
	iz0 = (int) (zoff - VERY_SMALL_DOUBLE);

	ix1 = 0;
	iy1 = (iy0 < pgrid->numy - 1) ? iy0 + 1 : iy0;
	iz1 = (iz0 < pgrid->numz - 1) ? iz0 + 1 : iz0;


	ydiff = yoff - (DOUBLE) iy0;
	zdiff = zoff - (DOUBLE) iz0;

//INGV
if (iy0 < 0 || iy1 >= pgrid->numy )
	return(-VERY_LARGE_DOUBLE);
if (iz0 < 0 || iz1 >= pgrid->numz )
	return(-VERY_LARGE_DOUBLE);

if (ydiff < 0.0 || ydiff > 1.0)
	return(-VERY_LARGE_DOUBLE);
if (zdiff < 0.0 || zdiff > 1.0)
	return(-VERY_LARGE_DOUBLE);

	/* location at grid node */

	if (ydiff + zdiff < SMALL_FLOAT) {
		if (fpgrid != NULL)
			value = ReadGrid3dValue(fpgrid,
				ix0, iy0, iz0, pgrid);
		else
			value = pgrid->array[ix0][iy0][iz0];
		return(value);
	}


	/* read vertex values from grid file or array */

	if (fpgrid != NULL) {
		vval00 = ReadGrid3dValue(fpgrid, ix0, iy0, iz0, pgrid);
		vval01 = ReadGrid3dValue(fpgrid, ix0, iy0, iz1, pgrid);
		vval10 = ReadGrid3dValue(fpgrid, ix0, iy1, iz0, pgrid);
		vval11 = ReadGrid3dValue(fpgrid, ix0, iy1, iz1, pgrid);
	} else {
		vval00 = pgrid->array[ix0][iy0][iz0];
		vval01 = pgrid->array[ix0][iy0][iz1];
		vval10 = pgrid->array[ix0][iy1][iz0];
		vval11 = pgrid->array[ix0][iy1][iz1];
	}

// INGV
// check for invalid / mask nodes
if (vval00 < 0.0 || vval01 < 0.0 || vval10 < 0.0 || vval11 < 0.0)
	return(-VERY_LARGE_DOUBLE);

	/* interpolate values */

	value = InterpSquareLagrange(ydiff, zdiff,
			vval00, vval01, vval10, vval11);


	return(value);
}


/*** function to find value inside a square using Lagrange interpolation***/

/* 	0.0 <= vvalKLM <= 1.0 */
/*	returns interp value at (xdiff, zdiff) */

INLINE DOUBLE InterpSquareLagrange(DOUBLE xdiff, DOUBLE zdiff,
		DOUBLE vval00, DOUBLE vval01, DOUBLE vval10, DOUBLE vval11)
{

	DOUBLE value;

	value =   vval00 * (1.0 - xdiff)  * (1.0 - zdiff)
		+ vval01 * (1.0 - xdiff)  * zdiff
		+ vval10 * xdiff          * (1.0 - zdiff)
		+ vval11 * xdiff          * zdiff;
	return (value);

}




/*** function to write hypocenter/arrivals to output */

int WriteLocation(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals,
		int narrivals, char* filename,
		int iWriteArrivals, int iWriteEndLoc, int iWriteMinimal,
		GridDesc* pgrid, int n_proj)
{
	return(_WriteLocation(fpio, phypo, parrivals,
		narrivals, filename,
		iWriteArrivals, iWriteEndLoc,iWriteMinimal,
		pgrid, n_proj, IO_ARRIVAL_ALL));
}


/*** function to write arrivals to output */

int WritePhases(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals,
		int narrivals, char* filename,
		int iWriteArrivals, int iWriteEndLoc, int iWriteMinimal,
		GridDesc* pgrid, int n_proj, int io_arrival_mode)
{
	return(_WriteLocation(fpio, phypo, parrivals,
		narrivals, filename,
		iWriteArrivals, iWriteEndLoc,iWriteMinimal,
		pgrid, n_proj, io_arrival_mode));
}


/*** function to write hypocenter/arrivals to output */

int _WriteLocation(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals,
		int narrivals, char* filename,
		int iWriteArrivals, int iWriteEndLoc, int iWriteMinimal,
		GridDesc* pgrid, int n_proj, int io_arrival_mode)
{

	int istat, ifile = 0, narr;
	ArrivalDesc* parr;
	double stat_dlat, stat_dlong;


	/* write hypocenter to file */

	if (fpio == NULL) {
		if ((fpio = fopen(filename, "w")) == NULL) {
			puterr("ERROR: opening hypocenter output file:");
			puterr(filename);
			return(-1);
		}
		NumFilesOpen++;
		ifile = 1;
	}

	/* write hypocenter parameters */

    if (iWriteMinimal) {
	fprintf(fpio, "NLLOC \"%s\" \"%s\" \" \"\n",
		phypo->fileroot, phypo->locStat);
    } else {
	fprintf(fpio, "NLLOC \"%s\" \"%s\" \"%s\"\n",
		phypo->fileroot, phypo->locStat, phypo->locStatComm);
	fprintf(fpio, "SIGNATURE \"%s\"\n", phypo->signature);
	fprintf(fpio, "COMMENT \"%s\"\n", phypo->comment);
	if (pgrid != NULL )
		fprintf(fpio, "GRID  %d %d %d  %lf %lf %lf  %lf %lf %lf %s\n",
			pgrid->numx, pgrid->numy, pgrid->numz,
			pgrid->origx, pgrid->origy, pgrid->origz,
			pgrid->dx, pgrid->dy, pgrid->dz, pgrid->chr_type);
	else
		fprintf(fpio, "GRID  %d %d %d  %lf %lf %lf  %lf %lf %lf %s\n",
			-1, -1, -1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "NULLGRID");
	fprintf(fpio, "SEARCH %s\n", phypo->searchInfo);
	fprintf(fpio,
"HYPOCENTER  x %lf y %lf z %lf  OT %lf  ix %d iy %d iz %d\n",
		phypo->x, phypo->y, phypo->z, (double) phypo->sec,
		phypo->ix, phypo->iy, phypo->iz);
    }
	fprintf(fpio,
"GEOGRAPHIC  OT %4.4d %2.2d %2.2d  %2.2d %2.2d %9.6lf  Lat %lf Long %lf Depth %lf\n",
//"GEOGRAPHIC  OT %4.4d %2.2d %2.2d  %2.2d %2.2d %lf  Lat %lf Long %lf Depth %lf\n",
		phypo->year, phypo->month, phypo->day,
		phypo->hour, phypo->min, (double) phypo->sec,
		phypo->dlat, phypo->dlong, phypo->depth);
	fprintf(fpio,
	    "QUALITY  Pmax %.3le MFmin %lf MFmax %lf RMS %lf Nphs %d Gap %d Dist %lf Mamp %5.2lf %d Mdur %5.2lf %d\n",
		phypo->probmax, phypo->misfit, phypo->grid_misfit_max,
			phypo->rms, phypo->nreadings, phypo->gap, phypo->dist,
			phypo->amp_mag, phypo->num_amp_mag,
			phypo->dur_mag, phypo->num_dur_mag
			);

	fprintf(fpio,
	    "VPVSRATIO  VpVsRatio %.3le  Npair %d  Diff %.3le\n",
		phypo->VpVs, phypo->nVpVs, phypo->tsp_min_max_diff
			);

    if (!iWriteMinimal) {

	fprintf(fpio,
		"STATISTICS  ExpectX %.4le Y %.4le Z %.4le",
		phypo->expect.x, phypo->expect.y, phypo->expect.z);
	fprintf(fpio, "  CovXX %.2le XY %.2le XZ %.2le",
		phypo->cov.xx, phypo->cov.xy, phypo->cov.xz);
	fprintf(fpio, " YY %.2le YZ %.2le",
		phypo->cov.yy, phypo->cov.yz);
	fprintf(fpio, " ZZ %.2le",
		phypo->cov.zz);
	fprintf(fpio, " EllAz1  %.1lf Dip1  %.1lf Len1  %.2le",
		phypo->ellipsoid.az1, phypo->ellipsoid.dip1,
				phypo->ellipsoid.len1);
	fprintf(fpio, " Az2  %.1lf Dip2  %.1lf Len2  %.2le",
		phypo->ellipsoid.az2, phypo->ellipsoid.dip2,
				phypo->ellipsoid.len2);
	fprintf(fpio, " Len3  %.2le\n", phypo->ellipsoid.len3);

	fprintf(fpio,
		"STAT_GEOG  ExpectLat %lf Long %lf Depth %lf\n",
		phypo->expect_dlat, phypo->expect_dlong, phypo->expect.z);

	fprintf(fpio, "%s\n", MapProjStr[n_proj]);
    }


	/* write arrival parameters */
	/* !!! Remember to modify the sscanf in function GetHypLoc()
						when modifying this fprintf */

	if (iWriteArrivals) {

		fprintf(fpio,
"PHASE ID Ins Cmp On Pha  FM Date     HrMn   Sec     Err  ErrMag    Coda      Amp       Per");
		if (io_arrival_mode == IO_ARRIVAL_ALL)
			fprintf(fpio,
"      >   TTpred    Res       Weight    StaLoc(X  Y         Z)        SDist    SAzim  RAz   RDip  RQual\n");
		else
			fprintf(fpio, "\n");
		for (narr = 0; narr < narrivals; narr++) {
			parr = parrivals + narr;
			WriteArrival(fpio, parr, io_arrival_mode);
		}
		fprintf(fpio, "END_PHASE\n");
	}


	/* write end line and blank line */
	if (iWriteEndLoc)
		fprintf(fpio, "END_NLLOC\n\n");

	if (ifile) {
		fclose(fpio);
		NumFilesOpen--;
	}
}



/*** function to read hypocenter/arrival parameters */

int GetHypLoc(FILE *fpio, char* filein, HypoDesc* phypo,
		ArrivalDesc* parrivals, int *pnarrivals, int iReadArrivals,
		GridDesc* pgrid, int n_proj)
{

	int istat, ifile = 0;
	int lineLength;
	char fn_in[FILENAME_MAX];
	char line[MAXLINE_LONG], *pstr, *pstr2;
	double hypo_sec, templat, templong;
	ArrivalDesc* parr;


	/* open hypocenter file */

	if (fpio == NULL) {
		sprintf(fn_in, "%s.hyp", filein);
		if ((fpio = fopen(fn_in, "r")) == NULL)
		{
			puterr("ERROR: opening hypocenter file.");
			return(-1);
		}
		NumFilesOpen++;
		ifile = 1;
	}


	/* read hypocenter paramters */

	/* search for first line of hypocenter description */
	do {
		if (fgets(line, MAXLINE_LONG, fpio) == NULL)
			goto eof_exit;
	} while (strncmp("NLLOC", line, 5));
	phypo->fileroot[0], '\0';
	if ((pstr = strchr(line, '"')) != NULL) {
		if ((pstr2 = strchr(pstr + 1, '"')) != NULL) {
			strncpy( phypo->fileroot, pstr + 1, pstr2 - pstr - 1);
			phypo->fileroot[pstr2 - pstr - 1] = '\0';
		}
	}
	phypo->locStat[0], '\0';
	if ((pstr = strchr(pstr2 + 1, '"')) != NULL) {
		if ((pstr2 = strchr(pstr + 1, '"')) != NULL) {
			strncpy( phypo->locStat, pstr + 1, pstr2 - pstr - 1);
			phypo->locStat[pstr2 - pstr - 1] = '\0';
		}
	}
	phypo->locStatComm[0], '\0';
	if ((pstr = strchr(pstr2 + 1, '"')) != NULL) {
		if ((pstr2 = strchr(pstr + 1, '"')) != NULL) {
			strncpy( phypo->locStatComm, pstr + 1,
							pstr2 - pstr - 1);
			phypo->locStatComm[pstr2 - pstr - 1] = '\0';
		}
	}


	/* read hypocenter description lines until END_NLLOC reached */

	while(1) {

		/* read next line */
		if (fgets(line, MAXLINE_LONG, fpio) == NULL)
			goto eof_exit;

		/* end of NLLoc hypocenter descrition */
		if (strncmp(line, "END_NLLOC", 9) == 0)
			break;

		lineLength = TrimString(line);

		/* blank line */
		if (lineLength == 0)
			continue;

		/* identify and read line */

		if (strncmp(line, "SIGNATURE", 9) == 0) {
			/* SIGNATURE */
			strcpy(phypo->signature, strchr(line, '"') + 1);
			*(strchr(phypo->signature, '"')) = '\0';
		}

		else if (strncmp(line, "COMMENT", 7) == 0) {
			/* COMMENT */
			strcpy(phypo->comment, strchr(line, '"') + 1);
			*(strchr(phypo->comment, '"')) = '\0';
		}

		else if (strncmp(line, "GRID", 4) == 0) {
			/* GRID */
			if (pgrid != NULL ) {
			if (sscanf(line,
					"%*s  %d %d %d  %lf %lf %lf  %lf %lf %lf %s",
					&pgrid->numx, &pgrid->numy, &pgrid->numz,
					&pgrid->origx, &pgrid->origy, &pgrid->origz,
					&pgrid->dx, &pgrid->dy, &pgrid->dz,
					pgrid->chr_type)
						== EOF)
				goto eof_exit;
			}
		}

		else if (strncmp(line, "SEARCH", 6) == 0) {
			/* SEARCH */
			line[MAXLINE_LONG - 1] = '\0';
			strcpy(phypo->searchInfo, strchr(line, ' ') + 1);
		}

		else if (strncmp(line, "HYPOCENTER", 10) == 0) {
			/* HYPOCENTER */
			if (sscanf(line,
"%*s %*s %lf %*s %lf %*s %lf %*s %lf %*s %d %*s %d %*s %d",
				&(phypo->x), &(phypo->y), &(phypo->z),
				&hypo_sec,
				&(phypo->ix), &(phypo->iy), &(phypo->iz))
					== EOF)
				goto eof_exit;
			phypo->sec = (long double) hypo_sec;
		}

		else if (strncmp(line, "GEOGRAPHIC", 10) == 0) {
			/* GEOGRAPHIC */
			if (sscanf(line,
"%*s %*s %d %d %d   %d %d %lf %*s %lf %*s %lf %*s %lf",
				&(phypo->year), &(phypo->month), &(phypo->day),
				&(phypo->hour), &(phypo->min), &hypo_sec,
				&(phypo->dlat), &(phypo->dlong),
				&(phypo->depth)) == EOF)
			goto eof_exit;
			phypo->sec = (long double) hypo_sec;
		}

		else if (strncmp(line, "QUALITY", 7) == 0) {
			/* QUALITY */
			if (sscanf(line,
"%*s %*s %lf %*s %lf %*s %lf %*s %lf %*s %d %*s %d %*s %lf %*s %lf %d %*s %lf %d",
				&phypo->probmax, &phypo->misfit,
				&phypo->grid_misfit_max,
				&phypo->rms, &phypo->nreadings, &phypo->gap,
				&phypo->dist,
				&phypo->amp_mag, &phypo->num_amp_mag,
				&phypo->dur_mag, &phypo->num_dur_mag
					) == EOF)
			goto eof_exit;
		}

		else if (strncmp(line, "VPVSRATIO", 9) == 0) {
			/* VPVSRATIO */
			if (sscanf(line,
"%*s %*s %lf %*s %d %*s %lf",
				&phypo->VpVs, &phypo->nVpVs, &phypo->tsp_min_max_diff
					) == EOF)
			goto eof_exit;
		}

		else if (strncmp(line, "STATISTICS", 10) == 0) {
			/* STATISTICS */
			if (sscanf(line,
				"%*s %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf",
				&phypo->expect.x, &phypo->expect.y,
				&phypo->expect.z,
				&phypo->cov.xx, &phypo->cov.xy, &phypo->cov.xz,
				&phypo->cov.yy, &phypo->cov.yz, &phypo->cov.zz,
				&phypo->ellipsoid.az1, &phypo->ellipsoid.dip1,
					&phypo->ellipsoid.len1,
				&phypo->ellipsoid.az2, &phypo->ellipsoid.dip2,
					&phypo->ellipsoid.len2,
				&phypo->ellipsoid.len3) == EOF)
					goto eof_exit;
			phypo->cov.yx = phypo->cov.xy;
			phypo->cov.zx = phypo->cov.xz;
			phypo->cov.zy = phypo->cov.yz;
		}

		else if (strncmp(line, "STAT_GEOG", 9) == 0) {
			/* STATISTICS */
			if (sscanf(line,
				"%*s %*s %lf %*s %lf %*s %lf",
				&phypo->expect_dlat, &phypo->expect_dlong, &phypo->expect.z) == EOF)
					goto eof_exit;
		}

		else if (strncmp(line, "ELLIPSOID", 9) == 0) {
			/* ELLIPSOID */
			if (sscanf(line,
				"%*s %*s %lf %lf %lf  %*s %lf %lf %lf  %*s %lf %lf %lf  %*s %lf",
				&templat, &templong,
				&phypo->expect.z,
				&phypo->ellipsoid.az1, &phypo->ellipsoid.dip1,
					&phypo->ellipsoid.len1,
				&phypo->ellipsoid.az2, &phypo->ellipsoid.dip2,
					&phypo->ellipsoid.len2,
				&phypo->ellipsoid.len3
				) == EOF)
					goto eof_exit;
				latlon2rect(0, templat, templong,
					&phypo->expect.x, &phypo->expect.y);
		}

		else if (strncmp(line, "FOCALMECH", 9) == 0) {
			/* FOCALMECH */
			/* Hyp dlat dlong depth Mech dipDir dipAng rake mf misfit nObs nObs */
			if (sscanf(line,
	"%*s %*s %lf %lf %lf %*s %lf %lf %lf %*s %lf %*s %d",
				&phypo->focMech.dlat, &phypo->focMech.dlong,
					&phypo->focMech.depth,
				&phypo->focMech.dipDir, &phypo->focMech.dipAng,
					&phypo->focMech.rake,
				&phypo->focMech.misfit, &phypo->focMech.nObs
				) == EOF)
					goto eof_exit;
		}

		else if (strncmp(line, "TRANSFORM", 9) == 0) {
			/* TRANSFORM */
			line[MAXLINE_LONG - 1] = '\0';
			strcpy(MapProjStr[n_proj], line);
		}

		else if (strncmp(line, "PHASE", 5) == 0) {

			/* read arrival/phase parameters */
			/* !!! Remember to modify this sscanf when modifying
				fprintf format  in function WriteLocation() */

			*pnarrivals = 0;

			/* if requested to read arrivals */
			if (iReadArrivals) {

				/* read arrivals */
				 while ((pstr = fgets(line, MAXLINE_LONG, fpio))
						!= NULL
					&& strncmp("END_PHASE", line, 9))
				{
					if (*pnarrivals >= MAX_NUM_ARRIVALS) {
						puterr(
"WARNING: maximum number of arrivals exceeded.");
						(*pnarrivals)--;
						break;
					}
					parr = parrivals + *pnarrivals;
					istat = ReadArrival(line, parr,
								IO_ARRIVAL_ALL);
					(*pnarrivals)++;
				};
			}

			/* not requested to read arrivals */
			else {
				 while ((pstr = fgets(line, MAXLINE_LONG, fpio))
						!= NULL
					&& strncmp("END_PHASE", line, 9))
						;
			}

			if (pstr == NULL)
				goto eof_exit;
		}

		else if (strncmp(line, "SCATTER", 7) == 0) {

			/* skip scatter points */
			 while ((pstr = fgets(line, MAXLINE_LONG, fpio))
					!= NULL
					&& strncmp("END_SCATTER", line, 11))
				;
			if (pstr == NULL)
				goto eof_exit;
		}

		else {

			putmsg(1,
"WARNING: unrecognized line in NLLOC hypocenter description:");
			sprintf(MsgStr, "   <%s>", line);
			putmsg(1, MsgStr);
		}

	}	/* while(1) */



	/* normal return */
	if (ifile) {
		fclose(fpio);
		NumFilesOpen--;
	}
	return(0);


	/* return at end of file */
	eof_exit:
	if (ifile) {
		fclose(fpio);
		NumFilesOpen--;
	}
	return(EOF);

}



/*** function to read arrival */

/* returns:	1    if only observation part of phase read
		2    if observation and calculated parts of phase read
		EOF  if EOF or error occurs before any values read
		-1   otherwise
*/

int ReadArrival(char* line, ArrivalDesc* parr, int iReadType)
{

	int istat;
	long int idate, ihrmin;
	char *line_calc;
	char label[10 * ARRIVAL_LABEL_LEN];


	/* read observation part of phase line */

	istat = sscanf(line, "%s %s %s %s %s %s %ld %ld %lf %s %lf %lf %lf %lf",
				label,
				parr->inst,
				parr->comp,
				parr->onset,
				parr->phase,
				parr->first_mot,
				/*&parr->quality, */
				&idate, &ihrmin,
				&(parr->sec),
				parr->error_type, &parr->error,
				&(parr->coda_dur),
				&(parr->amplitude),
				&(parr->period)
		);

	strncpy(parr->label, label, ARRIVAL_LABEL_LEN - 1);
/*printf("%s %s %s %s %s %s %ld %ld %lf %s %lf %lf %lf %lf\n",
	parr->label,
	parr->inst,
	parr->comp,
	parr->onset,
	parr->phase,
	parr->first_mot,
	//&parr->quality,
	idate, ihrmin,
	(parr->sec),
	parr->error_type, &parr->error,
	(parr->coda_dur),
	(parr->amplitude),
	(parr->period)
);
*/

	if (istat == EOF)
		return(istat);
	if (istat != 14)
		return(-1);

	/* decode data and time integers */
	parr->year = idate / 10000;
	idate = idate % 10000;
	parr->month = idate / 100;
	parr->day = idate % 100;
	parr->hour = ihrmin / 100;
	parr->min = ihrmin % 100;

	/* set null quality value */
	parr->quality = -1;

	if (iReadType == IO_ARRIVAL_ALL) {

		/* read calculated part of phase line */

		if ((line_calc = strchr(line, '>')) == NULL)
			return(1);

		istat = sscanf(line_calc + 1,
"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
				&parr->pred_travel_time, &parr->residual,
				&parr->weight,
				&parr->station.x, &parr->station.y,
				&parr->station.z,
				&parr->dist, &parr->azim,
				&parr->ray_azim, &parr->ray_dip, &parr->ray_qual
				);
		if (istat == EOF)
			return(istat);
		if (istat != 11)
			return(-1);

		/* convert ray aziumuth to grid coords direction */

		parr->ray_azim = latlon2rectAngle(0, parr->ray_azim);
	}

	return(2);

}






/*** function to write arrival */

int WriteArrival(FILE* fpio, ArrivalDesc* parr, int iWriteType)
{

	int istat;
	long int idate, ihrmin;
	char *line_calc;
	double sta_azim, ray_azim;


	/* code data and time integers */
	idate = parr->year * 10000 + parr->month * 100
		+ parr->day;
	ihrmin = parr->hour * 100 + parr->min;

	/* write observation part of phase line */

	istat = fprintf(fpio,
"%-6.6s %-4.4s %-4.4s %-1.1s %-6.6s %-1.1s %8.8ld %4.4ld %9.4lf %-3.3s %9.2le %9.2le %9.2le %9.2le",
				parr->label,
				parr->inst,
				parr->comp,
				parr->onset,
				parr->phase,
				parr->first_mot,
				/*parr->quality, */
				idate, ihrmin,
				parr->sec,
				parr->error_type, parr->error,
				parr->coda_dur,
				parr->amplitude,
				parr->period
		);
	if (istat < 0)
		return(-1);


	if (iWriteType == IO_ARRIVAL_ALL) {

		/* convert ray aziumuth to geographic direction */

		sta_azim = rect2latlonAngle(0, parr->azim);
		ray_azim = rect2latlonAngle(0, parr->ray_azim);

		/* write calculated part of phase line */

		istat = fprintf(fpio,
" > %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %6.2lf %5.1lf %5.1lf %2d",
				parr->pred_travel_time, parr->residual,
				parr->weight,
				parr->station.x, parr->station.y,
				parr->station.z,
				parr->dist, sta_azim,
				ray_azim, parr->ray_dip, parr->ray_qual
				);

		if (istat < 0)
			return(-1);

	}

	istat = fprintf(fpio, "\n");
	if (istat < 0)
		return(-1);

	return(0);

}



/*** function to write arrival in hypo71/hypoellipse format */

int WriteArrivalHypo(FILE* fpio, ArrivalDesc* arrival, int iwriteEOL)
{

	int istat;
	int pha_qual;


	/* write phase */

	pha_qual = (arrival->quality >= 0 && arrival->quality <= 4) ?
			arrival->quality : Err2Qual(arrival);
	if (pha_qual < 0)
		pha_qual = 0; /* !! not necessarily a good choice */


	if (iwriteEOL)
		istat = fprintf(fpio, "\n");

	/* P phase */
	if (strcmp(arrival->phase, "P") == 0) {

		/* write P arrival output */
		istat = fprintf(fpio, "%4.4s", arrival->label);
		istat = fprintf(fpio, "%1s", arrival->onset);
		istat = fprintf(fpio, "%1s", arrival->phase);
		istat = fprintf(fpio, "%1s", arrival->first_mot);
		istat = fprintf(fpio, "%1.1d", pha_qual);
		istat = fprintf(fpio, " %2.2d", arrival->year % 100);
		istat = fprintf(fpio, "%2.2d", arrival->month);
		istat = fprintf(fpio, "%2.2d", arrival->day);
		istat = fprintf(fpio, "%2.2d", arrival->hour);
		istat = fprintf(fpio, "%2.2d", arrival->min);
		istat = fprintf(fpio, "%5.2f", arrival->sec);

	/* S phase */
	} else if (strcmp(arrival->phase, "S") == 0) {

		/* write S phase output (assumes written directly after P) */
		istat = fprintf(fpio, "       %5.2f", arrival->sec);
		istat = fprintf(fpio, " %1s ", arrival->phase);
		istat = fprintf(fpio, "%1.1d", pha_qual);

	}

	if (istat < 0)
		return(-1);

	return(0);

}



/** function to generate GMT JVAL from map transformation parameters */

double getGMTJVAL(int n_proj, char* jval_string, double xlen, double vxmax, double vxmin,
	double ylen, double vymax, double vymin)
{

	double gmt_scale;

	jval_string[0] = '\0';

	if (map_itype[n_proj] == MAP_TRANS_GLOBAL) {

		/* -Jmscale or -JMwidth (Mercator [C])
                   Give scale along equator (1:xxxx or inch/degree).
		*/

		gmt_scale = xlen / (vxmax - vxmin);

		sprintf(jval_string, "-Jm%lf", gmt_scale);

		return(gmt_scale);

	} else if (map_itype[n_proj] == MAP_TRANS_SIMPLE) {

		/* -Jmscale or -JMwidth (Mercator [C])
                   Give scale along equator (1:xxxx or inch/degree).
		*/

		gmt_scale = xlen / (vxmax - vxmin);

		sprintf(jval_string, "-Jm%lf", gmt_scale);

		return(gmt_scale);

	} else if (map_itype[n_proj] == MAP_TRANS_LAMBERT) {

		/* -JLlon0/lat0/lat1/lat2/width (Lambert [C])
                   Give origin, 2 standard parallels, and scale along
                   these (1:xxxx or inch/degree)
		*/

		gmt_scale = (ylen / (vymax - vymin)); // * (xmaxrect0 - xminrect0) / (xmaxrect - xminrect);
		sprintf(jval_string, "-JL%lf/%lf/%lf/%lf/%lf",
			map_orig_long[n_proj], map_orig_lat[n_proj],
			map_lambert_1st_std_paral[n_proj], map_lambert_2nd_std_paral[n_proj],
			xlen);

		return(gmt_scale);

	}

	return(-1.0);

}


/** function to convert between coord systems */

int convertCoordsRect(int proj_index_from, int proj_index_to, double x, double y, double *pxnew, double *pynew)
{

	double dlat, dlong;

	if (proj_index_from < 0 || proj_index_to < 0)
		return(-1);

	if (proj_index_from == proj_index_to) {
		*pxnew = x;
		*pynew = y;
		return(0);
	}

	rect2latlon(proj_index_from, x, y, &dlat, &dlong);
	latlon2rect(proj_index_to, dlat, dlong, pxnew, pynew);

}



/** function to convert lat/long to rectangular km coord */

/* rotation about rect coord origin */

INLINE int latlon2rect(int n_proj, double dlat, double dlong, double* pxrect, double* pyrect)
{

	double xtemp, ytemp;


	if (map_itype[n_proj] == MAP_TRANS_GLOBAL) {
		*pxrect = dlong;
		*pyrect = dlat;
		return(0);

	} else if (map_itype[n_proj] == MAP_TRANS_SIMPLE) {
		xtemp = dlong - map_orig_long[n_proj];
		if (xtemp > 180.0)
			xtemp -= 360.0;
		if (xtemp < -180.0)
			xtemp += 360.0;
		xtemp = xtemp * c111 * cos(rpd * dlat);
		ytemp = (dlat - map_orig_lat[n_proj]) * c111;
		*pxrect = xtemp * map_cosang[n_proj] - ytemp * map_sinang[n_proj];
		*pyrect = ytemp * map_cosang[n_proj] + xtemp * map_sinang[n_proj];
		return(0);

	} else if (map_itype[n_proj] == MAP_TRANS_LAMBERT) {

		lamb(n_proj, dlong, dlat, &xtemp, &ytemp);
		xtemp /= 1000.0;	/* m -> km */
		ytemp /= 1000.0;	/* m -> km */
		*pxrect = xtemp * map_cosang[n_proj] - ytemp * map_sinang[n_proj];
		*pyrect = ytemp * map_cosang[n_proj] + xtemp * map_sinang[n_proj];

		return(0);
	}

	return(-1);

}



/** function to convert rectangular km coord to lat/long */

INLINE int rect2latlon(int n_proj, double xrect, double yrect, double* pdlat, double* pdlong)
{

	double xtemp, ytemp;

	if (map_itype[n_proj] == MAP_TRANS_GLOBAL) {
		xtemp = xrect * map_cosang[n_proj] + yrect * map_sinang[n_proj];
		ytemp = yrect * map_cosang[n_proj] - xrect * map_sinang[n_proj];
		*pdlat = yrect;
		*pdlong = xrect;

		return(0);

	} else if (map_itype[n_proj] == MAP_TRANS_SIMPLE) {
		xtemp = xrect * map_cosang[n_proj] + yrect * map_sinang[n_proj];
		ytemp = yrect * map_cosang[n_proj] - xrect * map_sinang[n_proj];
		*pdlat = map_orig_lat[n_proj] + ytemp / c111;
		*pdlong = map_orig_long[n_proj] + xtemp / (c111 * cos(rpd * *pdlat));

		return(0);

	} else if (map_itype[n_proj] == MAP_TRANS_LAMBERT) {

		xtemp = xrect * map_cosang[n_proj] + yrect * map_sinang[n_proj];
		ytemp = yrect * map_cosang[n_proj] - xrect * map_sinang[n_proj];
		ilamb(n_proj, pdlong, pdlat, xtemp * 1000.0, ytemp * 1000.0);

		return(0);
	}

	return(-1);
}



/** function to convert rectangular km angle to lat/long angle */

INLINE double rect2latlonAngle(int n_proj, double rectAngle)
{
	double angle;

	if (map_itype[n_proj] == MAP_TRANS_SIMPLE || map_itype[n_proj] == MAP_TRANS_LAMBERT) {
		angle = rectAngle - map_rot[n_proj];
		if (angle < 0.0)
			angle += 360.0;
		else if (angle > 360.0)
			angle -= 360.0;
		return(angle);
	} else
		return(rectAngle);
}



/** function to convert lat/long km angle to rectangular angle */

INLINE double latlon2rectAngle(int n_proj, double latlonAngle)
{
	double angle;

	if (map_itype[n_proj] == MAP_TRANS_SIMPLE || map_itype[n_proj] == MAP_TRANS_LAMBERT) {
		angle = latlonAngle + map_rot[n_proj];
		if (angle < 0.0)
			angle += 360.0;
		else if (angle > 360.0)
			angle -= 360.0;
		return(angle);
	} else
		return(latlonAngle);
}



/** function to convert source location parameters between coord systems */

int ConvertSourceLoc(int n_proj, SourceDesc *source, int numSources,
					int toXY, int toLatLon)
{

	int istat, nsource;
	SourceDesc *srce_in;


	for (nsource = 0; nsource < numSources; nsource++) {

		srce_in = source + nsource;

		if (toXY && srce_in->is_coord_latlon &&
						!srce_in->is_coord_xyz) {
			istat = latlon2rect(n_proj, srce_in->dlat, srce_in->dlong,
				&(srce_in->x), &(srce_in->y));
			srce_in->z = srce_in->depth;
		}
		if (toLatLon && srce_in->is_coord_xyz &&
					!srce_in->is_coord_latlon) {
			istat = rect2latlon(n_proj, srce_in->x, srce_in->y,
				&(srce_in->dlat), &(srce_in->dlong));
			srce_in->depth = srce_in->z;
		}
	}

	return(istat);

}



/*** function to set model coordinates mode to rect or latlon */

int SetModelCoordsMode(int num_surfaces)
{
	double xloc, yloc;

	/* if surfaces read, assume lat/long */
	if (num_surfaces > 0) {
		ModelCoordsMode = COORDS_LATLON;
		/* check geographic transformation */
		if (rect2latlon(0, 0.0, 0.0, &yloc, &xloc) < 0) {
			puterr(
"FATAL ERROR: geographic transformation required with SURFACE options,\n\tbut transformation (TRANS) not initialized.");
			exit(-1);
		}
	} else
		ModelCoordsMode = COORDS_RECT;

}




/*** date functions */


/*** function to convert character month to integer month */

int Month2Int(char* cmonth)
{
	int i;

	for (i = 0; i < strlen(cmonth); i++)
		cmonth[i] = toupper(cmonth[i]);

	if (strcmp(cmonth, "JAN") == 0)
		return(1);
	if (strcmp(cmonth, "FEB") == 0)
		return(2);
	if (strcmp(cmonth, "MAR") == 0)
		return(3);
	if (strcmp(cmonth, "APR") == 0)
		return(4);
	if (strcmp(cmonth, "MAY") == 0)
		return(5);
	if (strcmp(cmonth, "JUN") == 0)
		return(6);
	if (strcmp(cmonth, "JUL") == 0)
		return(7);
	if (strcmp(cmonth, "AUG") == 0)
		return(8);
	if (strcmp(cmonth, "SEP") == 0)
		return(9);
	if (strcmp(cmonth, "OCT") == 0)
		return(10);
	if (strcmp(cmonth, "NOV") == 0)
		return(11);
	if (strcmp(cmonth, "DEC") == 0)
		return(12);

	puterr2("ERROR: unrecognized charcter month", cmonth);

	return(0);

}


/*** date functions */

static char  daytab[2][13] = {
	{0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
	{0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
};

/*** function to set day of year */

int DayOfYear(int year, int month, int day)
{
	int i, leap;

	leap = year%4 == 0 && year%100 != 0 || year%400 == 0;
	for (i = 1; i < month; i++)
		day += daytab[leap][i];

	return day;

}

/*** function to set month / day from day of year */

void MonthDay(int year, int yearday, int* pmonth, int* pday)
{
	int i, leap;

	leap = year%4 == 0 && year%100 != 0 || year%400 == 0;
	for (i = 1; yearday > daytab[leap][i]; i++)
		yearday -= daytab[leap][i];
	*pmonth = i;
	*pday = yearday;

}


/** function to construct current date/time string */

char* CurrTimeStr(void)
{
	static char timestr[MAXLINE];
	time_t curr_time;

	curr_time = time(NULL);

	strftime(timestr, (size_t) MAXLINE,
				"%d%b%Y %Hh%Mm%S", localtime(&curr_time));
	return(timestr);

}



/** function to check for and expand wild card characters in filenames
	and to return a list of equivalent files */

int ExpandWildCards(char* fileName, char fileList[][FILENAME_MAX_SMALL],
	int maxNumFiles)
{
	int istat;
	int nfiles = 0;
	char system_str[MAXLINE];
	char list_file[] = "filelist.tmp";
	char *pchr;
	FILE* fpio;


	/* check for no '*' or '?' character */

	if ((pchr = strchr(fileName, '*')) == NULL
			&& (pchr = strchr(fileName, '?')) == NULL) {
		strcpy(fileList[0], fileName);
		nfiles = 1;
		return(nfiles);
	}


	/* expand wildcard file names into list of files */

	sprintf(system_str, "ls %s > %s", fileName, list_file);
	system(system_str);

	if ((fpio = fopen(list_file, "r")) == NULL) {
		puterr("ERROR: opening fileList:  temporary file.");
		return(-1);
	}
	NumFilesOpen++;

	nfiles = 0;
	while(nfiles < maxNumFiles &&
			(istat = fscanf(fpio, "%s", fileList[nfiles])) != EOF
			&& istat == 1) {
		nfiles++;
	}


	fclose(fpio);
	NumFilesOpen--;

	return(nfiles);


}



/*** function to sort arrivals by obs_time field */

int SortArrivalsTime(ArrivalDesc* arrival, int num_arrivals)
{

	qsort((void *) arrival, (size_t) num_arrivals, sizeof(ArrivalDesc),
		(int (*)(const void *, const void * )) CmpArrivalsTime);


	return(0);

}

/*** function to compare arrivals by obs_time field */

int CmpArrivalsTime(const ArrivalDesc *keyval, const ArrivalDesc *datum )
{
	if (keyval->obs_time < datum->obs_time)
		return(-1);
	if (keyval->obs_time > datum->obs_time)
		return(1);

	return(0);
}


/*** function to sort arrivals by flag_ignore field */

int SortArrivalsIgnore(ArrivalDesc* arrival, int num_arrivals)
{

	qsort((void *) arrival, (size_t) num_arrivals, sizeof(ArrivalDesc),
		(int (*)(const void *, const void * )) CmpArrivalsIgnore);


	return(0);

}

/*** function to compare arrivals by flag_ignore field */

int CmpArrivalsIgnore(const ArrivalDesc *keyval, const ArrivalDesc *datum )
{
	if (keyval->flag_ignore < datum->flag_ignore)
		return(-1);
	if (keyval->flag_ignore > datum->flag_ignore)
		return(1);

	return(0);
}


/*** function to sort arrivals by dist field */

int SortArrivalsDist(ArrivalDesc* arrival, int num_arrivals)
{

	qsort((void *) arrival, (size_t) num_arrivals, sizeof(ArrivalDesc),
		(int (*)(const void *, const void * )) CmpArrivalsDist);


	return(0);

}

/*** function to compare arrivals by dist field */

int CmpArrivalsDist(const ArrivalDesc *keyval, const ArrivalDesc *datum )
{
	if (keyval->dist < datum->dist)
		return(-1);
	if (keyval->dist > datum->dist)
		return(1);

	return(0);
}


/*** function to sort array of doubles */

int SortDoubles(double* array, int num_elements)
{

	qsort((void *) array, (size_t) num_elements, sizeof(double),
		(int (*)(const void *, const void * )) CmpDoubles);


	return(0);

}

/*** function to compare doubles */

int CmpDoubles(const double *keyval, const double *datum )
{
	if (*keyval < *datum)
		return(-1);
	if (*keyval > *datum)
		return(1);

	return(0);
}




/*** function to set angle values in take-off angles union */

TakeOffAngles SetTakeOffAngles(double azim, double dip, int iqual)
{
	TakeOffAngles angles;

	angles.ival[1] = (unsigned short) (0.5 + 10.0 * azim);
	angles.ival[0] = (unsigned short) iqual
		+ (unsigned short) ANGLES_OFFSET
		* (unsigned short) (0.5 + 10.0 * dip);

	return(angles);
}


/*** function to set float values in take-off angles union */

void SetAnglesFloat(TakeOffAngles* pangles, float fvalue)
{
	pangles->fval = fvalue;
}


/*** function to get values in take-off angles union */

int GetTakeOffAngles(TakeOffAngles *pangles,
			double *pazim, double *pdip, int *piqual)
{
	*pazim = ((double) pangles->ival[1]) / 10.0;
	*pdip = ((double) (pangles->ival[0] / (int) ANGLES_OFFSET)) / 10.0;
	*piqual = (int) pangles->ival[0] % (int) ANGLES_OFFSET;

	return(*piqual);
}


/*** function to read take-off angles from file */

int ReadTakeOffAnglesFile(char *fname, double xloc, double yloc, double zloc,
		double *pazim, double *pdip, int *piqual, double sta_azim, int iSwapBytes)
{
	int istat;
	FILE *fp_grid, *fp_hdr;
	float fvalue;
	GridDesc gdesc;
	TakeOffAngles angles;

	/* open angle grid file */
	if ((istat =  OpenGrid3dFile(fname, &fp_grid, &fp_hdr,
			&gdesc, "angle", NULL, iSwapBytes)) < 0) {
		putmsg(3,
"WARNING: cannot open angle grid file, ignoring angles.");
		angles = SetTakeOffAngles(0.0, 0.0, 0);
		GetTakeOffAngles(&angles, pazim, pdip, piqual);
		return(-1);
	}

	/* get angles float value on grid */
	fvalue = ReadAbsInterpGrid3d(fp_grid, &gdesc, xloc, yloc, zloc);

	/* get angles */
	SetAnglesFloat(&angles, fvalue);
	GetTakeOffAngles(&angles, pazim, pdip, piqual);

	/* determine azimuth (2D grids) */
	if (gdesc.type == GRID_ANGLE_2D) {
		if (*pazim > 0.0)
			*pazim = sta_azim;
		else {
			*pazim = sta_azim - 180.0;
			if (*pazim < 0.0)
				*pazim += 360.0;
		}
	}

	/* close angle grid file */
	CloseGrid3dFile(&fp_grid, &fp_hdr);

	return(0);


}



/** function to generate and display "traditional" statistics from grid file */

int GenTraditionStats(GridDesc *pgrid, Vect3D *pexpect, Mtrx3D *pcov,
	FILE* fpgrid)
{
	int istat;


	/* allocate grid buffer */

	pgrid->buffer = AllocateGrid(pgrid);
	if (pgrid->buffer == NULL) {
		puterr(
"ERROR: allocating memory for 3D PDF grid buffer.");
		exit(EXIT_ERROR_MEMORY);
	}

	/* create grid array access pointers */

	pgrid->array = CreateGridArray(pgrid);
	if (pgrid->array == NULL) {
		puterr(
"ERROR: creating array for accessing 3D PDF grid buffer.");
		exit(EXIT_ERROR_MEMORY);
	}


	/* read PDF grid */

	if (istat =
		        ReadGrid3dBuf(pgrid, fpgrid)< 0) {
		puterr("ERROR: reading PDF grid from disk.");
		exit(EXIT_ERROR_IO);
	}


	/* calculate expectation */

	*pexpect = CalcExpectation(pgrid, NULL);
	sprintf(MsgStr, "EXPECTATION { x %lf  y %lf  z %lf }",
		pexpect->x, pexpect->y, pexpect->z);
	putmsg(3, MsgStr);

	/* calculate covariance matrix */

	*pcov = CalcCovariance(pgrid, pexpect, NULL);
	sprintf(MsgStr, "COVARIANCE: {");
	putmsg(3, MsgStr);
	sprintf(MsgStr, "   xx: %lf  xy: %lf  xz: %lf",
		pcov->xx, pcov->xy, pcov->xz);
	putmsg(3, MsgStr);
	sprintf(MsgStr, "   yx: %lf  yy: %lf  yz: %lf",
		pcov->yx, pcov->yy, pcov->yz);
	putmsg(3, MsgStr);
	sprintf(MsgStr, "   zx: %lf  zy: %lf  zz: %lf",
		pcov->zx, pcov->zy, pcov->zz);
	putmsg(3, MsgStr);
	sprintf(MsgStr, "}");
	putmsg(3, MsgStr);
	/* clean up */

	FreeGrid(pgrid);
	DestroyGridArray(pgrid);

	return(0);

}



/** function to calculate the expectation (mean) of a PDF grid */

Vect3D CalcExpectation(GridDesc* pgrid, FILE* fpgrid)
{

	int ix, iy, iz;

	float val;
	double volume;
	Vect3D expect = {0.0, 0.0, 0.0};


	/* cannot calculate for misfit grid */
	if (pgrid->type == GRID_MISFIT) {
		expect.x = expect.y = expect.z = -LARGE_DOUBLE;
		return(expect);
	}


	for (ix = 0; ix < pgrid->numx; ix++) {
		for (iy = 0; iy < pgrid->numy; iy++) {
			for (iz = 0; iz < pgrid->numz; iz++) {

				if (fpgrid != NULL)
					val = ReadGrid3dValue(fpgrid,
						ix, iy, iz, pgrid);
				else
					val = pgrid->array[ix][iy][iz];

				expect.x += (double) val * (double) ix;
				expect.y += (double) val * (double) iy;
				expect.z += (double) val * (double) iz;

			}
		}
	}

	volume = pgrid->dx * pgrid->dy * pgrid->dz;
	expect.x = pgrid->origx + expect.x * pgrid->dx * volume;
	expect.y = pgrid->origy + expect.y * pgrid->dy * volume;
	expect.z = pgrid->origz + expect.z * pgrid->dz * volume;


	return(expect);
}




/** function to calculate the covariance a PDF grid */

Mtrx3D CalcCovariance(GridDesc* pgrid, Vect3D* pexpect, FILE* fpgrid)
{

	int ix, iy, iz;

	float val;
	double x, y, z, xx, xy, xz, yy, yz, zz;
	double volume;

	Mtrx3D cov = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};



	/* cannot calculate for misfit grid */
	if (pgrid->type == GRID_MISFIT) {
		cov.xx =  cov.xy =  cov.xz =
			cov.yx = cov.yy =  cov.yz =
			cov.zx =  cov.zy = cov.zz
				= -LARGE_DOUBLE;
		return(cov);
	}



	/* calculate covariance following eq. (6-12), T & V, 1982 */


	for (ix = 0; ix < pgrid->numx; ix++) {
		x = pgrid->origx + (double) ix * pgrid->dx;
		xx = x * x;

		for (iy = 0; iy < pgrid->numy; iy++) {
			y = pgrid->origy + (double) iy * pgrid->dy;
			yy = y * y;
			xy = x * y;

			for (iz = 0; iz < pgrid->numz; iz++) {
				z = pgrid->origz + (double) iz * pgrid->dz;
				xz = x * z;
				yz = y * z;
				zz = z * z;

				if (fpgrid != NULL)
					val = ReadGrid3dValue(fpgrid,
						ix, iy, iz, pgrid);
				else
					val = pgrid->array[ix][iy][iz];

				cov.xx += (double) val * xx;
				cov.xy += (double) val * xy;
				cov.xz += (double) val * xz;

				cov.yy += (double) val * yy;
				cov.yz += (double) val * yz;

				cov.zz += (double) val * zz;

			}
		}
	}

	volume = pgrid->dx * pgrid->dy * pgrid->dz;

	cov.xx = cov.xx * volume - pexpect->x * pexpect->x;
	cov.xy = cov.xy * volume - pexpect->x * pexpect->y;
	cov.xz = cov.xz * volume - pexpect->x * pexpect->z;

	cov.yx = cov.xy;
	cov.yy = cov.yy * volume - pexpect->y * pexpect->y;
	cov.yz = cov.yz * volume - pexpect->y * pexpect->z;

	cov.zx = cov.xz;
	cov.zy = cov.yz;
	cov.zz = cov.zz * volume - pexpect->z * pexpect->z;


	return(cov);
}



/** function to calculate the expectation (mean)  of a set of samples */

Vect3D CalcExpectationSamples(float* fdata, int nSamples)
{

	int nsamp, ipos;

	float x, y, z, prob;
	Vect3D expect = {0.0, 0.0, 0.0};


	ipos = 0;
	for (nsamp = 0; nsamp < nSamples; nsamp++) {
		x = fdata[ipos++];
		y = fdata[ipos++];
		z = fdata[ipos++];
		prob = fdata[ipos++];
		expect.x += (double) x;
		expect.y += (double) y;
		expect.z += (double) z;
	}

	expect.x /= (double) nSamples;
	expect.y /= (double) nSamples;
	expect.z /= (double) nSamples;

	return(expect);
}




/** function to calculate the covariance of a set of samples */

Mtrx3D CalcCovarianceSamples(float* fdata, int nSamples, Vect3D* pexpect)
{

	if (GeometryMode == MODE_GLOBAL)
		return(CalcCovarianceSamplesGlobal(fdata, nSamples, pexpect));
	else
		return(CalcCovarianceSamplesRect(fdata, nSamples, pexpect));

}



/** function to calculate the covariance of a set of samples */

Mtrx3D CalcCovarianceSamplesRect(float* fdata, int nSamples, Vect3D* pexpect)
{

	int nsamp, ipos;

	float x, y, z, prob;
	Mtrx3D cov = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


	/* calculate covariance following eq. (6-12), T & V, 1982 */

	ipos = 0;
	for (nsamp = 0; nsamp < nSamples; nsamp++) {
		x = fdata[ipos++];
		y = fdata[ipos++];
		z = fdata[ipos++];
		prob = fdata[ipos++];

		cov.xx += (double) (x * x);
		cov.xy += (double) (x * y);
		cov.xz += (double) (x * z);

		cov.yy += (double) (y * y);
		cov.yz += (double) (y * z);

		cov.zz += (double) (z * z);

	}

	cov.xx = cov.xx / (double) nSamples - pexpect->x * pexpect->x;
	cov.xy = cov.xy / (double) nSamples - pexpect->x * pexpect->y;
	cov.xz = cov.xz / (double) nSamples - pexpect->x * pexpect->z;

	cov.yx = cov.xy;
	cov.yy = cov.yy / (double) nSamples - pexpect->y * pexpect->y;
	cov.yz = cov.yz / (double) nSamples - pexpect->y * pexpect->z;

	cov.zx = cov.xz;
	cov.zy = cov.yz;
	cov.zz = cov.zz / (double) nSamples - pexpect->z * pexpect->z;


	return(cov);
}



/** function to calculate the covariance of a set of samples in long(deg)/lat(deg)/depth(km) coordinates */

Mtrx3D CalcCovarianceSamplesGlobal(float* fdata, int nSamples, Vect3D* pexpect)
{

	int nsamp, ipos;

	float x, y, z, prob;
	Mtrx3D cov = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


	/* calculate covariance following eq. (6-12), T & V, 1982 */

	ipos = 0;
	for (nsamp = 0; nsamp < nSamples; nsamp++) {
		x = fdata[ipos++] * DEG2KM;
		y = fdata[ipos++] * DEG2KM;
		z = fdata[ipos++];
		prob = fdata[ipos++];

		cov.xx += (double) (x * x);
		cov.xy += (double) (x * y);
		cov.xz += (double) (x * z);

		cov.yy += (double) (y * y);
		cov.yz += (double) (y * z);

		cov.zz += (double) (z * z);

	}

	cov.xx = cov.xx / (double) nSamples - pexpect->x * pexpect->x * DEG2KM * DEG2KM;
	cov.xy = cov.xy / (double) nSamples - pexpect->x * pexpect->y * DEG2KM * DEG2KM;
	cov.xz = cov.xz / (double) nSamples - pexpect->x * pexpect->z * DEG2KM;

	cov.yx = cov.xy;
	cov.yy = cov.yy / (double) nSamples - pexpect->y * pexpect->y * DEG2KM * DEG2KM;
	cov.yz = cov.yz / (double) nSamples - pexpect->y * pexpect->z * DEG2KM;

	cov.zx = cov.xz;
	cov.zy = cov.yz;
	cov.zz = cov.zz / (double) nSamples - pexpect->z * pexpect->z;


	return(cov);
}



/*** function to calculate confidence ellipsoid from covariance matrix */

/* 	finds confidence ellipsoid from SVD of Cov mtrx.  See
		Num Rec, 2nd ed, secs 2.6 & 15.6

		del_chi_2 is delta Chi-square (see Num Rec, 2nd ed, fig 15.6.5)
*/


Ellipsoid3D CalcErrorEllipsoid(Mtrx3D *pcov, double del_chi_2)
{
	int istat, ndx, iSwitched;
	Matrix A_matrix, V_matrix;
	Vector W_vector;
	float wtemp, vtemp;
	Ellipsoid3D ell;



	/* allocate A mtrx */
	A_matrix = matrix(0, 3, 0, 3);

	/* load A matrix in NumRec format */
	A_matrix[0][0] = pcov->xx;
	A_matrix[0][1] = A_matrix[1][0] = pcov->xy;
	A_matrix[0][2] = A_matrix[2][0] = pcov->xz;
	A_matrix[1][1] = pcov->yy;
	A_matrix[1][2] = A_matrix[2][1] = pcov->yz;
	A_matrix[2][2] = pcov->zz;


	/* allocate V mtrx and W vector */
	V_matrix = matrix(0, 2, 0, 2);
	W_vector = vector(0, 2);

	/* do SVD */
	if ((istat = svdcmp0(A_matrix, 3, 3, W_vector, V_matrix)) < 0) {
		puterr("ERROR: doing SVD for confidence ellipsoids.");
		return(EllipsoidNULL);
	}

	/* check
	if (W_vector[0] < SMALL_DOUBLE || W_vector[1] < SMALL_DOUBLE
			|| W_vector[2] < SMALL_DOUBLE) {
		puterr(
"ERROR: invalid SVD singular value for confidence ellipsoids.");
		return(EllipsoidNULL);
	}


	/* sort by singular values W */
	iSwitched = 1;
	while (iSwitched) {
		iSwitched = 0;
		for (ndx = 0; ndx < 2; ndx++) {
			if (W_vector[ndx] > W_vector[ndx + 1]) {
				wtemp = W_vector[ndx];
				W_vector[ndx] = W_vector[ndx + 1];
				W_vector[ndx + 1] = wtemp;
				vtemp = V_matrix[0][ndx];
				V_matrix[0][ndx] = V_matrix[0][ndx + 1];
				V_matrix[0][ndx + 1] = vtemp;
				vtemp = V_matrix[1][ndx];
				V_matrix[1][ndx] = V_matrix[1][ndx + 1];
				V_matrix[1][ndx + 1] = vtemp;
				vtemp = V_matrix[2][ndx];
				V_matrix[2][ndx] = V_matrix[2][ndx + 1];
				V_matrix[2][ndx + 1] = vtemp;
				iSwitched = 1;
			}
		}
	}


	/* calculate ellipsoid axes */
	/* length: w in Num Rec, 2nd ed, fig 15.6.5 must be replaced
		by 1/sqrt(w) since we are using SVD of Cov mtrx and not
		SVD of A mtrx (compare eqns 2.6.1  & 15.6.10) */

	ell.az1 = atan2(V_matrix[0][0], V_matrix[1][0]) / rpd;
	if (ell.az1 < 0.0)
		ell.az1 += 360.0;
	ell.dip1 = asin(V_matrix[2][0]) / rpd;
	ell.len1 = sqrt(del_chi_2) / sqrt(1.0 / W_vector[0]);
	ell.az2 = atan2(V_matrix[0][1], V_matrix[1][1]) / rpd;
	if (ell.az2 < 0.0)
		ell.az2 += 360.0;
	ell.dip2 = asin(V_matrix[2][1]) / rpd;
	ell.len2 = sqrt(del_chi_2) / sqrt(1.0 / W_vector[1]);
	ell.len3 = sqrt(del_chi_2) / sqrt(1.0 / W_vector[2]);

	return(ell);

}




/*** function to read expectation and covariance from a hyp file */

int ReadHypStatistics(FILE **pfpio, char* fnroot_in,
			Vect3D* pmax_like, Vect3D* pexpect,
			Mtrx3D* pcov, Ellipsoid3D *pellipsoid,
			ArrivalDesc* parrivals, int *pnarrivals)
{
	int idummy;
	char fn_in[FILENAME_MAX];
	static HypoDesc hypo;


	/* open hypocenter file if necessary */

	if (*pfpio == NULL) {
		sprintf(fn_in, "%s.hyp", fnroot_in);
		if ((*pfpio = fopen(fn_in, "r")) == NULL)
		{
			puterr("ERROR: opening hypocenter file.");
			return(EOF);
		}
		NumFilesOpen++;
	}


	/* read next hypocenter */

	if (GetHypLoc(*pfpio, fnroot_in, &hypo, parrivals, pnarrivals, 1, NULL, 0) != EOF) {
		pmax_like->x = hypo.x;
		pmax_like->y = hypo.y;
		pmax_like->z = hypo.z;
		*pexpect = hypo.expect;
		*pcov = hypo.cov;
		*pellipsoid = hypo.ellipsoid;
		return(0);
	}


	/* end of file */

	fclose(*pfpio);
	NumFilesOpen--;

	return(EOF);

}




/*** function to read mechanism from a hyp file */

int ReadFocalMech(FILE **pfpio, char* fnroot_in,
			FocalMech* pfocalMech,
			ArrivalDesc* parrivals, int *pnarrivals)
{
	int idummy;
	char fn_in[FILENAME_MAX];
	static HypoDesc hypo;


	/* open hypocenter file if necessary */

	if (*pfpio == NULL) {
		sprintf(fn_in, "%s.hyp", fnroot_in);
		if ((*pfpio = fopen(fn_in, "r")) == NULL)
		{
			puterr("ERROR: opening hypocenter file.");
			return(EOF);
		}
		NumFilesOpen++;
	}


	/* read next hypocenter */

	if (GetHypLoc(*pfpio, fnroot_in, &hypo, parrivals, pnarrivals, 1, NULL, 0) != EOF) {
		*pfocalMech = hypo.focMech;
		return(0);
	}


	/* end of file */

	fclose(*pfpio);
	NumFilesOpen--;

	return(EOF);

}




/*** function to trim leading and trailing blanks from character strings */

int TrimString(char* line)
{
	char *line0, *line1, *pendchr;

	/* find end of string */
	pendchr = strchr(line, '\0');
	if (pendchr == NULL)
		return(-1);

	/* remove leading spaces */
	while(isspace(line[0])) {
		line0 = line + 1;
		line1 = line;
		do
			*(line1++) = *line0;
		while (*(line0++) != '\0');
	}

	/* re-find end of string */
	pendchr = strchr(line, '\0');
	if (pendchr == NULL)
		return(-1);

	/* remove trailing spaces */
	while (--pendchr > line && isspace(*pendchr))
		*pendchr = '\0';

	return(pendchr - line);

}


/*** function to check for blank line */

int LineIsBlank(char *line)
{
	char* pstr;

	pstr = line;
	while (*pstr) {
		if (!isspace(*pstr++))
			return(0);
	}

	return(1);
}



/*** function to read FORTRAN character format fields */

int ReadFortranString(char* line, int istart, int ilen, char* string_in)
{
	char chrtmp[MAXLINE], *chrpos;
	int istat, n;


	/* check for null char before end of read */
	chrpos = strchr(line, '\0');
	if (chrpos - line < istart + ilen - 1)
		return(-1);

	/* copy field to temp string */
	strncpy(chrtmp, line + (istart - 1), ilen);
	chrtmp[ilen] = '\0';

	/* check for blank field */
/*	for (n = 0; n < ilen; n++) {
		if (chrtmp[n] != ' ')
			break;
	}
	if (n == ilen) {
		for (n = 0; n < ilen; n++)
			string_in[n] = ' ';
		return(1);
	}
*/

	strncpy(string_in, chrtmp, ilen);
	string_in[ilen] = '\0';

	return(1);

}


/*** function to read FORTRAN integer format fields */

int ReadFortranInt(char* line, int istart, int ilen, int* pintval)
{
	char chrtmp[MAXLINE], *chrpos;
	int istat, n;


	/* check for null char before end of read */
	chrpos = strchr(line, '\0');
	if (chrpos - line < istart + ilen - 1)
		return(-1);

	/* copy field to temp string */
	strncpy(chrtmp, line + (istart - 1), ilen);
	chrtmp[ilen] = '\0';

	/* check for blank field */
	for (n = 0; n < ilen; n++) {
		if (chrtmp[n] != ' ')
			break;
	}
	if (n == ilen) {
		*pintval = 0;
		return(1);
	}


	istat = sscanf(chrtmp, "%d", pintval);

	return(istat);

}


/*** function to read FORTRAN real format fields */

int ReadFortranReal(char* line, int istart, int ilen, double* pdblval)
{
	char chrtmp[MAXLINE], *chrpos;
	int istat, n;


	/* check for null char before end of read */
	chrpos = strchr(line, '\0');
	if (chrpos - line < istart + ilen - 1)
		return(-1);

	/* copy field to temp string */
	strncpy(chrtmp, line + (istart - 1), ilen);
	chrtmp[ilen] = '\0';

	/* check for blank field */
	for (n = 0; n < ilen; n++) {
		if (chrtmp[n] != ' ')
			break;
	}
	if (n == ilen) {
		*pdblval = 0.0;
		return(1);
	}


	istat = sscanf(chrtmp, "%lf", pdblval);

	return(istat);

}



/*** function to convert phase error to hypo71 style quality */

int Err2Qual(ArrivalDesc *arrival)
{
	int errLevel;

	/* find error level with error >= arrival error */
	for (errLevel = 0; errLevel < NumQuality2ErrorLevels; errLevel++) {
		if (arrival->error <= Quality2Error[errLevel]) {
			arrival->quality = errLevel;
			return(errLevel);
		}
	}

	/* not found or no quality/error values available */
	return(-1);
}

/*** function to convert quality to error */

void Qual2Err(ArrivalDesc *arrival)
{

	/* set error fields */
	strcpy(arrival->error_type, "GAU");
	if (arrival->quality >= 0 &&
		arrival->quality < NumQuality2ErrorLevels) {
			arrival->error =
				Quality2Error[arrival->quality];
	} else {
		arrival->error =
			Quality2Error[NumQuality2ErrorLevels - 1];
		puterr("WARNING: invalid arrival quality.");
	}

}



/*** function to read Quality2Error values ***/

int GetQuality2Err(char* line1)
{

	int istat, ierr, nlev;
	double qual2err;
	char frmt[MAXLINE] = "%lf";
	char frmttmp[MAXLINE];

	NumQuality2ErrorLevels = 0;

	while ((istat = sscanf(line1, frmt, &qual2err)) == 1) {
		Quality2Error[NumQuality2ErrorLevels++] =  qual2err;
		sprintf(frmttmp, "%%*lf %s", frmt);
		strcpy(frmt, frmttmp);
	}


	sprintf(MsgStr, "NLLoc LOCQUAL2ERR:");
	putmsg(2, MsgStr);
	ierr = 0;
	for (nlev = 0; nlev < NumQuality2ErrorLevels; nlev++) {
		sprintf(MsgStr, " %d ->  %lf", nlev, Quality2Error[nlev]);
		putmsg(2, MsgStr);
		if ((istat = checkRangeDouble("QUAL2ERR", "Quality2Error",
				Quality2Error[nlev], 1, 0.0, 0, 0.0)) != 0)
			ierr = -1;
	}

	if (ierr < 0)
		return(-1);

	return(0);
}







/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

