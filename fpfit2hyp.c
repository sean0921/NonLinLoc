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


/*   fpfit2hyp.c

	Program to convert fpfit summary files to NLLoc .hyp format

*/

/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */


/*
	history:

	ver 01    24Nov1998  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "GridLib.h"


/* defines */



/* globals */



/* functions */

int Fpfit2Hyp(int , char **);
int ReadFpfitSum(FILE *fp_in, HypoDesc *phypo);



/*** program to sum event scatter files */

#define PNAME  "fpfit2hyp"


main(int argc, char *argv[])
{

	int istat, narg;


	/* set program name */

	strcpy(prog_name, PNAME);


	/* check command line for correct usage */

	fprintf(stdout, "\n%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc < 3) {
		puterr("ERROR wrong number of command line arguments.");
		disp_usage(PNAME, 
"<fpfit_file> <out_hyp_file> [MechMisfitMax [RMSMax [NRdgsMin [GapMax]]]]");
		exit(-1);
	}

	if ((istat = Fpfit2Hyp(argc, argv)) < 0) {
		puterr("ERROR converting fpfitise summary file.");
		exit(-1);
	}



	exit(0);

}



int Fpfit2Hyp(int argc, char *argv[])
{
	int istat, narg;
	int nLocWritten, nLocRead;
	char fn_hyp_out[FILENAME_MAX];
	char fn_fpfit_in[FILENAME_MAX];
	FILE *fp_hyp_out, *fp_fpfit_in;


	double MechMisfitMax, RMSMax;
	int NRdgsMin, GapMax;

	GridDesc Grid, locgrid;
	SourceDesc* Srce;
	HypoDesc Hypo, *phypo;



	/* get command line parameters */
	strcpy(fn_fpfit_in, argv[1]);
	strcpy(fn_hyp_out, argv[2]);

	MechMisfitMax = 1.0e6;
	if (argc > 3) {
		sscanf(argv[3], "%lf", &MechMisfitMax);
	}
	fprintf(stdout, "  Mechanism Misfit Maximum: %lf\n", MechMisfitMax);
	RMSMax = 1.0e6;
	if (argc > 4) {
		sscanf(argv[4], "%lf", &RMSMax);
	}
	fprintf(stdout, "  RMS Maximum: %lf\n", RMSMax);
	NRdgsMin = 0;
	if (argc > 5) {
		sscanf(argv[5], "%d", &NRdgsMin);
	}
	fprintf(stdout, "  Num Readings Minimum: %d\n", NRdgsMin);
	GapMax = 360;
	if (argc > 6) {
		sscanf(argv[6], "%d", &GapMax);
	}
	fprintf(stdout, "  Gap Maximum: %d\n", GapMax);



	/* open fpfit summary file */

	if ((fp_fpfit_in = fopen(fn_fpfit_in, "r")) == NULL) {
		puterr("ERROR: opening scatter output file.");
		return(-1);
	}

	/* open ascii hypocenter output file */
	if ((fp_hyp_out = fopen(fn_hyp_out, "w")) == NULL) {
		puterr("ERROR: opening scatter ascii output file.");
		return(-1);
	}



	nLocWritten = 0;
	nLocRead = 0;
	while (1) {


		strcpy(Hypo.fileroot, fn_fpfit_in);
		strcpy(Hypo.locStat, "Converted from FPFIT");
		Hypo.grid_misfit_max = -1.0;

		/* read next hypocenter */
	    	if ((istat = ReadFpfitSum(fp_fpfit_in, &Hypo)) == EOF)
			break;
		else if (istat < 0) {
			puterr2(
"ERROR: reading fpfit summary file", fn_fpfit_in);
			break;
		}
		nLocRead++;

		if (Hypo.focMech.misfit > MechMisfitMax) {
/*			puterr(
"WARNING: solution misfit is greater than MFMax, ignoring event");*/
			continue;
		} else if (Hypo.rms > RMSMax) {
/*			puterr(
"WARNING: location RMS is Greater than RMSMax, ignoring event");*/
			continue;
		} else if (Hypo.nreadings < NRdgsMin) {
/*			puterr(
"WARNING: location num readings is less than NRdgsMin, ignoring event");*/
			continue;
		} else if (Hypo.gap > GapMax) {
/*			puterr(
"WARNING: location gap is greater than GapMax, ignoring event");*/
			continue;
		} else {
			NumArrivals = 0;
			WriteLocation(fp_hyp_out, &Hypo,
				Arrival, NumArrivals, fn_hyp_out, 0, 0, 1,
				&locgrid, 0);
		}


		/* write ellipsoid */
		phypo = &Hypo;
		fprintf(fp_hyp_out, 
			"ELLIPSOID  Hyp  %lf %lf %lf",
			phypo->dlat, phypo->dlong, phypo->depth);
		fprintf(fp_hyp_out, 
			" Ell1  %.1lf %.1lf %.2le",
			phypo->ellipsoid.az1, phypo->ellipsoid.dip1, 
				phypo->ellipsoid.len1);
		fprintf(fp_hyp_out, " Ell2  %.1lf %.1lf %.2le",
			phypo->ellipsoid.az2, phypo->ellipsoid.dip2, 
				phypo->ellipsoid.len2);
		fprintf(fp_hyp_out, " Ell3  %.2le", phypo->ellipsoid.len3);
		fprintf(fp_hyp_out, "\n");


		/* write mechanism */
		phypo = &Hypo;
		fprintf(fp_hyp_out, 
			"FOCALMECH  Hyp  %lf %lf %lf",
			phypo->dlat, phypo->dlong, phypo->depth);
		fprintf(fp_hyp_out, 
			" Mech  %.1lf %.1lf %.1lf",
			phypo->focMech.dipDir, phypo->focMech.dipAng, 
			phypo->focMech.rake);
		fprintf(fp_hyp_out, " mf  %.2lf nObs %d",
			phypo->focMech.misfit, phypo->focMech.nObs);
		fprintf(fp_hyp_out, "\n");


		/* write end line and blank line */
		fprintf(fp_hyp_out, "END_NLLOC\n\n");

		nLocWritten++;

	}


	fclose(fp_hyp_out);

	/* write message */
	fprintf(stdout, 
"%d locations read, %d written to ascii hyp file <%s>\n", 
		nLocRead, nLocWritten, fn_hyp_out);


	return(0);

}



/*** function to read fpfit summary record to HypoDesc structure */

int ReadFpfitSum(FILE *fp_in, HypoDesc *phypo)
{

	int istat;
	char *cstat;
	double mag, dtemp;

	double deg, dmin;
	char strNS[2], strMagType[2];

	static char line[MAXLINE_LONG];


	/* read next line */
	cstat = fgets(line, MAXLINE_LONG, fp_in);
 	if (cstat == NULL)
		return(EOF);

	/* read hypocenter parameters */

	istat = 0;
	istat += ReadFortranInt(line, 1, 2, &phypo->year);
	if (phypo->year < 100)
		phypo->year += 1900;
	istat += ReadFortranInt(line, 3, 2, &phypo->month);
	istat += ReadFortranInt(line, 5, 2, &phypo->day);    
	istat += ReadFortranInt(line, 8, 2, &phypo->hour);
	istat += ReadFortranInt(line, 10, 2, &phypo->min);
	istat += ReadFortranReal(line, 12, 6, &phypo->sec);

	istat += ReadFortranReal(line, 18, 3, &deg);
	istat += ReadFortranString(line, 21, 1, strNS);
	istat += ReadFortranReal(line, 22, 5, &dmin);
	phypo->dlat = deg + dmin / 60.0;
	if (strncmp(strNS, "S", 1) == 0)
		phypo->dlat = -phypo->dlat;
	istat += ReadFortranReal(line, 27, 4, &deg);
	istat += ReadFortranString(line, 31, 1, strNS);
	istat += ReadFortranReal(line, 32, 5, &dmin);
	phypo->dlong = deg + dmin / 60.0;
	if (strncmp(strNS, "W", 1) == 0)
		phypo->dlong = -phypo->dlong;

	istat += ReadFortranReal(line, 37, 7, &phypo->depth);

	istat += ReadFortranReal(line, 46, 5, &mag);

	istat += ReadFortranInt(line, 51, 3, &phypo->nreadings);
	istat += ReadFortranReal(line, 54, 4, &dtemp);
	phypo->gap = (int) 0.5 + dtemp;
	istat += ReadFortranReal(line, 58, 5, &phypo->dist);
	istat += ReadFortranReal(line, 63, 5, &phypo->rms);

	/* hypoinverse horiz and vertical error are converted to ellipsoid,
			this is not correct statistically  */
	istat += ReadFortranReal(line, 68, 5, &(phypo->ellipsoid.len1));
	phypo->ellipsoid.az1 = 0.0;
	phypo->ellipsoid.dip1 = 0.0;
	phypo->ellipsoid.len2 = phypo->ellipsoid.len1;
	phypo->ellipsoid.az2 = 90.0;
	phypo->ellipsoid.dip2 = 0.0;
	istat += ReadFortranReal(line, 73, 5, &(phypo->ellipsoid.len3));

	istat += ReadFortranString(line, 80, 1, strMagType);

	/* focal mechanism parameters */

	istat += ReadFortranReal(line, 82, 3, &phypo->focMech.dipDir);
	istat += ReadFortranReal(line, 86, 2, &phypo->focMech.dipAng);
	istat += ReadFortranReal(line, 88, 4, &phypo->focMech.rake);
	istat += ReadFortranReal(line, 94, 4, &phypo->focMech.misfit);
	istat += ReadFortranInt(line, 99, 3, &phypo->focMech.nObs);
	istat += ReadFortranReal(line, 103, 5, &phypo->focMech.misfit90);
	istat += ReadFortranReal(line, 109, 4, &phypo->focMech.staDist);
	istat += ReadFortranReal(line, 114, 4, &phypo->focMech.ratioMH);
	istat += ReadFortranReal(line, 120, 2, &phypo->focMech.conf90strike);
	istat += ReadFortranReal(line, 123, 2, &phypo->focMech.conf90dip);
	istat += ReadFortranReal(line, 126, 2, &phypo->focMech.conf90rake);
	istat += ReadFortranString(line, 128, 1, phypo->focMech.convFlag);
	istat += ReadFortranString(line, 129, 1, phypo->focMech.multSolFlag);

	return(istat);



}


/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

