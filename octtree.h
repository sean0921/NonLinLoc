/* 
 * Copyright (C) 1999-2000 Anthony Lomax <lomax@geoazur.unice.fr>
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


/*  octree.h

	include file for octree search

*/



/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@geoazur.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#ifdef EXTERN_MODE
#define	EXTERN_TXT extern
#else
#define EXTERN_TXT
#endif

	/* misc defines */

#ifndef SMALL_DOUBLE
#define SMALL_DOUBLE 1.0e-20
#endif
#ifndef LARGE_DOUBLE
#define LARGE_DOUBLE 1.0e20
#endif
#ifndef VERY_SMALL_DOUBLE
#define VERY_SMALL_DOUBLE 1.0e-30
#endif
#ifndef VERY_LARGE_DOUBLE
#define VERY_LARGE_DOUBLE 1.0e30
#endif


/*------------------------------------------------------------/ */
/* structures */
/*------------------------------------------------------------/ */


/* octree node */

typedef struct octnode* OctNodePtr;
typedef struct octnode
{
	OctNodePtr parent;		/* parent node */
	Vect3D center;			/* absolute coordinates of center */
	Vect3D ds;			/* length of sides */
	double value;			/* node value */
	OctNodePtr child[2][2][2];	/* child nodes */
	int isLeaf;			/* leaf flag, 1=leaf */
} OctNode;



/* 3D tree with Nx, Ny, Nz arbitrary */

typedef struct
{
	OctNode**** nodeArray;		/* parent nodes */
	int numx, numy, numz;		/* grid size */
 	Vect3D orig; 	/* orig (km) */
	Vect3D ds;		/* len side (km) */
}
Tree3D;


/* structure for storing results */

typedef struct resultTreeNode* ResultTreeNodePtr;
typedef struct resultTreeNode {
	ResultTreeNodePtr left;		/* address of left node */
	ResultTreeNodePtr right;	/* address of right node */
	double value;			/* prob * volume */
	OctNode* pnode;			/* correspnding octree node */
} ResultTreeNode;



/* */
/*------------------------------------------------------------/ */



/*------------------------------------------------------------/ */
/* globals  */
/*------------------------------------------------------------/ */

//EXTERN_TXT char fn_control[MAXLINE];	/* control file name */

/* */
/*------------------------------------------------------------/ */




/*------------------------------------------------------------/ */
/* function declarations */
/*------------------------------------------------------------/ */

Tree3D* newTree3D(int numx, int numy, int numz, 
	double origx, double origy, double origz,
	double dx,  double dy,  double dz, double value);
OctNode* newOctNode(OctNode* parent, Vect3D center, Vect3D ds, double value);
void subdivide(OctNode* parent, double value);
void freeTree3D(Tree3D* tree);
void freeNode(OctNode* parent);
OctNode* getLeafNodeContaining(Tree3D* tree, Vect3D coords);
OctNode* getLeafContaining(OctNode* node, double x, double y, double z);

ResultTreeNode* addResult(ResultTreeNode* prtn, double value, OctNode* pnode);
void freeResultTree(ResultTreeNode* prtn);
ResultTreeNode*  getHighestValue(ResultTreeNode* prtn);
ResultTreeNode* getHighestLeafValue(ResultTreeNode* prtree);


/* */
/*------------------------------------------------------------/ */


