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


/* util.h 

   AJL utility functions.

*/


/*
	by Anthony Lomax
	Geosciences Azur, Valbonne, France
*/



#include <stdlib.h>
#include <stdio.h>


#ifndef MAXLINE
#define MAXLINE 101
#endif

#ifndef VERY_SMALL_DOUBLE
#define VERY_SMALL_DOUBLE 1.0e-30
#endif

EXTERN_TXT char package_name[MAXLINE];
EXTERN_TXT char prog_name[MAXLINE];
EXTERN_TXT char prog_ver[MAXLINE];
EXTERN_TXT char prog_date[MAXLINE];
EXTERN_TXT char prog_copyright[MAXLINE];
EXTERN_TXT int message_flag;
EXTERN_TXT char MsgStr[2 * MAXLINE];


/*** function to display command usage */
void disp_usage(char* , char*);

/*** function to display error message */
void puterr(char *);

/*** function to display error message */
void puterr2(char *, char *);

/*** function to display message */
void putmsg(int , char *);

/*** function to display program information */
void DispProgInfo();

/*** function to check that int val is in range */
int checkRangeInt(char* name, char* param, int val, 
	int checkMin, int min, int checkMax, int max);

/*** function to check that double val is in range */
int checkRangeDouble(char* name, char* param, double val, 
	int checkMin, double min, int checkMax, double max);

/* misc structures */

/* 3D vector or point */
struct Vect3D {
	double x;
	double y;
	double z;
};




