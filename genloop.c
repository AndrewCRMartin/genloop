/*************************************************************************

   Program:    genloop
   File:       genloop.c
   
   Version:    V2.4
   Date:       03.11.17
   Function:   Build loop backbones from phi/psi/omega data
   
   Copyright:  (c) SciTech Software 1993-2017
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be freely copied
   and distributed for no charge providing this header is included.
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! The code may not be sold commercially without prior permission 
   from the author, although it may be given away free with commercial 
   products, providing it is made clear that this program is free and 
   that the source code is provided with the program.

   Special Provisions
   ==================
   Oxford Molecular Limited (OML) are granted rights to distribute this
   code as part of the AbM package providing all copyright information
   is retained in all source code files pertaining to this program.

**************************************************************************

   Description:
   ============
   genloop reads an AbM searchdb hit file and torsion file and creates
   a backbone (N,CA,C) or CA trace for the loops. The output file is in
   PDB format or in a smaller format if the -s flag is used.

   This version assumes that the torsion file also contains omega angles
   (at present, 10.11.93, it does not). Using the -o flag causes omega
   angles all to be set at 180.0 and allows reading of the current
   searchdb torsion files.

   19.11.93 A modified version of OML's AbM searchdb now exists which will
   generate the omega angles.

   07.06.94 A slight re-write to use output files from Martin's neural
   networks.

**************************************************************************

   Usage:
   ======
   genloop [-c] [-o] [-s] [-n] <hitfile> <torfile> <outfile>
   -c causes only C-alpha coordinates to be written
   -o causes no omega angles to be used; all 180.0
   -s causes small format output
   -n use neural network data

Normal input format is:

nres nentries            (Free format integers)

title line
---------------------------------------------------------
   X    phi      psi      omega
   X    phi      psi      omega
   X    phi      psi      omega
   ...

title line
---------------------------------------------------------
   X    phi      psi      omega
   X    phi      psi      omega
   X    phi      psi      omega
   ...

Rigid columns are required: (3X,A1,4X,3(F8,1X))


The neural network file format is the same but has pairs of
actual and predicted loop conformations. Each pair is separated 
by a > sign


The HITS FILE is simply a set of title lines which are copied
into the output file, but which must contain the sequence in
parentheses.

EXAMPLE HITS FILE
=================
First hit (RASEN)
Second hit (RASQS)


EXAMPLE TORSIONS FILE
=====================
     5    2   (Loop Length and number of entries)

Residue    PHI      PSI     OMEGA
----------------------------------
   R        -     127.593  180.000
   A     -82.379  146.780  178.243
   S     -65.588  -17.103  182.203
   E    -160.615 -164.887  179.372
   N     -97.135     -        -

Residue    PHI      PSI     OMEGA
----------------------------------
   R        -     138.362  177.265
   A    -122.742  108.730  182.157
   S     -49.318   -3.193  180.275
   Q    -142.030  152.142  179.286
   S     -51.111     -        -

**************************************************************************

   Notes:
   ======
   
**************************************************************************

   Revision History:
   =================
   V1.0  09.11.93 Original
   V1.1  10.11.93 Added omega handling and -o and -s options
   V1.2  19.11.93 Now creates the correct sequence in PDB style output
   V2.0  07.06.94 Rewritten to use Martin's network output files
   V2.1  10.06.94 Bug fix in fixed omegas (=180.0)
   V2.2  28.07.94 Bug fix in call to GetTorsions()
   V2.3  25.10.94 Another Bug fix in call to GetTorsions()
   V2.4  03.11.17 Updated for new Bioplib and various tidying for clean
                  compile for ISO C90

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/seq.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines
*/
#define MAXTOR 20   /* Max number of torsions (and therefore residues)  */

/* Parameters from CHARMm                                               */
#define ANGLE_C   117.5*PI/180.0
#define ANGLE_N   120.0*PI/180.0
#define ANGLE_CA  111.6*PI/180.0
#define DIST_CA_C 1.52
#define DIST_C_N  1.33
#define DIST_N_CA 1.49

/************************************************************************/
/* Macros
*/

/************************************************************************/
/* Type definitions
*/
struct _entry
{
   struct _entry *next,
                 *prev;
   REAL          x,
                 y,
                 z;
   char          atnam[8];
}  ;
typedef struct _entry ENTRY;

/************************************************************************/
/* Globals
*/
FILE *gHitfp = NULL,
     *gTorfp = NULL,
     *gOutfp = NULL;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ProcessCmdLine(int argc, char **argv, BOOL *calphas, char *hitfile,
                    char *torfile, char *outfile, BOOL *DoOmega,
                    BOOL *SmallOut, BOOL *DoNNFile);
BOOL OpenFiles(char *hitfile, char *torfile, char *outfile);
void DoProcessing(BOOL calphas, BOOL DoOmega, BOOL SmallOut, 
                  BOOL DoNNFile);
int GetTorsions(FILE *fp, int nres, REAL *phi, REAL *psi, BOOL DoOmega,
                REAL *omega);
ENTRY *BuildFragment(int ntor, REAL *phi, REAL *psi, REAL *omega);
BOOL BuildAtom(ENTRY *p, ENTRY *q, ENTRY *r, REAL theta, REAL bond, 
               REAL phi, REAL *coords);
void WriteFragment(FILE *gOutfp, ENTRY *entry, BOOL calphas, 
                   BOOL SmallOut, char *seq);
char *GetSequence(char *buffer);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for building loops from torsion data produced by searchdb

   08.11.93 Original   By: ACRM
   10.11.93 Modified to support omega angle
   07.06.94 Modified to support neural nets files
*/
int main(int argc, char **argv)
{
   char hitfile[160],
        torfile[160],
        outfile[160];
   BOOL calphas  = FALSE,
        DoOmega  = TRUE,
        SmallOut = FALSE,
        DoNNFile = FALSE;

   if(ProcessCmdLine(argc, argv, &calphas, hitfile, torfile, outfile, 
                     &DoOmega, &SmallOut, &DoNNFile))
   {
      /* Open the files                                                 */
      if(OpenFiles(hitfile, torfile, outfile))
      {
         DoProcessing(calphas, DoOmega, SmallOut, DoNNFile);
      }
   }
   else
   {
      printf("\ngenloop V2.4 (c) 1993-2017, Dr. Andrew C.R. Martin\n");
      printf("\nUsage: genloop [-c] [-o] [-s] [-n] <hitfile> <torfile> \
<outfile>\n");
      printf("       -c causes only C-alpha coordinates to be written\n");
      printf("       -o causes no omega angles to be used; all 180.0\n");
      printf("       -s causes a small format output to be created\n");
      printf("       -n read Martin Reczko's neural net output\n");
      
      printf("\nBuilds a file containing coordinates from the hit \
file and\n");
      printf("torsion file created by AbM or Martin Reczko's neural \
nets\n");

      printf("\n\n<hitfile> is simply a set of title lines which are \
copied\n");
      printf("into the output file.\n");

      printf("\n<torfile> format is:\n");

      printf("\nnres nentries            (Free format integers)\n");

      printf("\n   title line\n");
      printf("-----------------------------------------------------\
----\n");
      printf("   X    phi      psi      omega\n");
      printf("   X    phi      psi      omega\n");
      printf("   X    phi      psi      omega\n");
      printf("   ...\n");

      printf("\n   title line\n");
      printf("-----------------------------------------------------\
----\n");
      printf("   X    phi      psi      omega\n");
      printf("   X    phi      psi      omega\n");
      printf("   X    phi      psi      omega\n");
      printf("   ...\n");
      printf("\n");
      printf("Rigid columns are required: (3X,A1,4X,3(F8,1X))\n");

      printf("\nThe neural network file format is the same but has \
pairs of\n");
      printf("actual and predicted loop conformations. Each pair is \
separated \n");
      printf("by a > sign\n");

      printf("\n\nEXAMPLE TORSIONS FILE\n");
      printf("=====================\n");
      printf("     5    2\n");

      printf("\nResidue    PHI      PSI     OMEGA\n");
      printf("----------------------------------\n");
      printf("   R        -     127.593  180.000\n");
      printf("   A     -82.379  146.780  178.243\n");
      printf("   S     -65.588  -17.103  182.203\n");
      printf("   E    -160.615 -164.887  179.372\n");
      printf("   N     -97.135     -        -\n");

      printf("\nResidue    PHI      PSI     OMEGA\n");
      printf("----------------------------------\n");
      printf("   R        -     138.362  177.265\n");
      printf("   A    -122.742  108.730  182.157\n");
      printf("   S     -49.318   -3.193  180.275\n");
      printf("   Q    -142.030  152.142  179.286\n");
      printf("   S     -51.111     -        -\n");
      printf("\n");

      printf("NOTE: The OMEGA angle is placed with the amino acid \
BEFORE the peptide\n");
      printf("      bond instead of the one after as is more normal\n");

      printf("\nThis program is freely distributable providing no \
profit is made in so doing.\n\n");
   }

   return(0);
}

/************************************************************************/
/*>BOOL ProcessCmdLine(int argc, char **argv, BOOL *calphas, 
                       char *hitfile, char *torfile, char *outfile, 
                       BOOL *DoOmega, BOOL *SmallOut, BOOL *DoNNFile)
   ------------------------------------------------------------------
   Process the command line

   08.11.93 Original   By: ACRM
   10.11.93 Added -o and -s flags
   07.06.94 Added -n flag
*/
BOOL ProcessCmdLine(int argc, char **argv, BOOL *calphas, char *hitfile,
                    char *torfile, char *outfile, BOOL *DoOmega,
                    BOOL *SmallOut, BOOL *DoNNFile)
{
   argc--; argv++;

   while(argc > 3)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
	 case 'c': case 'C':
            *calphas = TRUE;
            break;
	 case 'o': case 'O':
            *DoOmega = FALSE;
            break;
	 case 's': case 'S':
            *SmallOut = TRUE;
            break;
         case 'n': case 'N':
            *DoNNFile = TRUE;
            break;
	 default:
            fprintf(stderr,"Unknown switch %s, ignored\n",argv[0]);
            break;
         }
      }
      else
      {
         fprintf(stderr,"Command line confused at: %s\n",argv[0]);
      }

      argc--; argv++;
   }

   if(argc!=3)
      return(FALSE);

   /* Copy in the filename arguments                                    */
   strcpy(hitfile,argv[0]);
   strcpy(torfile,argv[1]);
   strcpy(outfile,argv[2]);

   return(TRUE);
}

/************************************************************************/
/*>BOOL OpenFiles(char *hitfile, char *torfile, char *outfile)
   -----------------------------------------------------------
   Open files reporting any errors

   08.11.93 Original   By: ACRM
*/
BOOL OpenFiles(char *hitfile, char *torfile, char *outfile)
{
   if((gHitfp = fopen(hitfile,"r")) == NULL)
   {
      printf("Unable to open input file: %s\n", hitfile);
      return(FALSE);
   }
   if((gTorfp = fopen(torfile,"r")) == NULL)
   {
      printf("Unable to open input file: %s\n", torfile);
      return(FALSE);
   }
   if((gOutfp = fopen(outfile,"w")) == NULL)
   {
      printf("Unable to open output file: %s\n", outfile);
      return(FALSE);
   }

   return(TRUE);
}

/************************************************************************/
/*>void DoProcessing(BOOL calphas, BOOL DoOmega, BOOL SmallOut, 
                     BOOL DoNNFile)
   ------------------------------------------------------------
   Read headers then read info line from hits file and torsion data
   from torsion file then call routine to build conformation.

   08.11.93 Original   By: ACRM
   10.11.93 Added omega processing and small output flag
   19.11.93 Modified to create correct sequence in output PDB type
            listings
   07.06.94 Added DoNNFile parameter handling
   08.06.94 Added reading of loop length from torsion file
   28.07.94 Fixed bug in call to GetTorsions() By: MR
   25.10.94 The change to ntor rather than nres meant that it only
            worked for Neural net files as the fix was in the wrong
            place. Now done correctly. 
*/
void DoProcessing(BOOL calphas, BOOL DoOmega, BOOL SmallOut, 
                  BOOL DoNNFile)
{
   REAL  phi[MAXTOR],
         psi[MAXTOR],
         omega[MAXTOR];
   ENTRY *entry;
   int   ntor,
         nres,
         nentries;
   char  HitBuffer[160],
         buffer[160],
         *seq;

   /* Read header from torsions file                                    */
   fgets(HitBuffer,160,gTorfp);
   sscanf(HitBuffer,"%d %d",&nres,&nentries);

   while(fgets(HitBuffer,160,gHitfp))
   {
      fprintf(gOutfp,"REMARK %s",HitBuffer);
      seq = GetSequence(HitBuffer);

      /* If it's a neural net file, we get the number of residues next  */
      if(DoNNFile)
      {
         fgets(buffer,160,gTorfp);
         TERMINATE(buffer);
         sscanf(buffer,"%d",&nres);   /* ACRM: nres not ntor            */
      }

      /* Get the torsion data (if a neural net file, this is the 
         target data)
         MR: ntor not nres in paramters
         25.10.94 ACRM: nres again, fix should have been in the sscanf() 
                        above
      */
      if((ntor = GetTorsions(gTorfp,nres,phi,psi,DoOmega,omega)) != 0)
      {
         entry = BuildFragment(ntor,phi,psi,omega);
         WriteFragment(gOutfp,entry,calphas,SmallOut,seq);
         FREELIST(entry,ENTRY);
      }

      /* If it's a neural net file we now get the predicted data        */
      if(DoNNFile)
      {
         /* MR: ntor not nres in paramters 
            25.10.94 ACRM: fixed back to nres
         */
         if((ntor = GetTorsions(gTorfp,nres,phi,psi,DoOmega,omega)) != 0)
         {
            entry = BuildFragment(ntor,phi,psi,omega);
            fprintf(gOutfp,"TER\n");
            
            WriteFragment(gOutfp,entry,calphas,SmallOut,seq);
            fprintf(gOutfp,"TER\n");
            FREELIST(entry,ENTRY);
         }
      }

      /* If it's a neural net file, skip along until we read a > sign   */
      if(DoNNFile)
      {
         while(fgets(buffer,160,gTorfp))
         {
            if(buffer[0] == '>') break;
         }
      }
   }
}

/************************************************************************/
/*>int GetTorsions(FILE *fp, int nres, REAL *phi, REAL *psi, 
                   BOOL DoOmega, REAL *omega)
   ---------------------------------------------------------
   Get torsion angle data for a loop from torsion file

   09.11.93 Original   By: ACRM
   10.11.93 Added omega input
   07.06.94 Modified to skip all lines until the line of ----s
   10.06.94 Corrected fixed omegas to PI not 180.0
*/
int GetTorsions(FILE *fp, int nres, REAL *phi, REAL *psi, BOOL DoOmega,
                REAL *omega)
{
   char TorBuffer[160],
        *pphi   = TorBuffer+5,
        *ppsi   = TorBuffer+17,
        *pomega = TorBuffer+26;
   int  i;

   /* Skip any lines up to and including the line of ------s
      This includes a blank line and the title line; also a line starting
      with > to separate pairs
    */
   for(TorBuffer[0] = '\0'; TorBuffer[0] != '-';)
      if(!fgets(TorBuffer,160,fp)) return(0);

   for(i=0; i<nres; i++)
   {
      if(!fgets(TorBuffer,160,fp)) return(0);
      TERMINATE(TorBuffer);

      /* Split the buffer with NULL characters                          */
      TorBuffer[4]  = '\0';
      TorBuffer[16] = '\0';
      TorBuffer[25] = '\0';

      /* Read the phi, psi & omega values treating the first and last 
         cases specially as only some torsions will be defined.
      */
      if(i==0)
      {
         phi[i] = (REAL)0.00;
         psi[i] = (REAL)(PI*atof(ppsi)/180.0);

         if(DoOmega)
            omega[i] = (REAL)(PI*atof(pomega)/180.0);
         else
            omega[i] = (REAL)PI;
      }
      else if(i==(nres-1))
      {
         phi[i]   = (REAL)(PI*atof(pphi)/180.0);
         psi[i]   = (REAL)0.0;
         omega[i] = (REAL)0.0;
      }
      else
      {
         phi[i]   = (REAL)(PI*atof(pphi)/180.0);
         psi[i]   = (REAL)(PI*atof(ppsi)/180.0);

         if(DoOmega)
            omega[i] = (REAL)(PI*atof(pomega)/180.0);
         else
            omega[i] = (REAL)PI;
      }
   }

   return(nres);
}


/************************************************************************/
/*>ENTRY *BuildFragment(int ntor, REAL *phi, REAL *psi, REAL *omega)
   -----------------------------------------------------------------
   Builds a fragment from a set of torsion angles. Assuming a fixed
   trans omega angle for the peptide bond:
   each phi angle defines the coordinates of this residues's C (and O)
   each psi angle defines the coordinates of the next residue's N
   each omega angle defines the coordinates of the next residue's CA
   The data are built into a ENTRY linked list.

   09.11.93 Original   By: ACRM
   10.11.93 Added omega parameter
*/
ENTRY *BuildFragment(int ntor, REAL *phi, REAL *psi, REAL *omega)
{
   ENTRY *entry, *p;
   REAL  coords[3];
   int   i;

   /* Check input                                                       */
   if(ntor == 0) return(NULL);

   /* Initialise an ENTRY structure for the start N                     */
   INITPREV(entry,ENTRY);
   if(entry == NULL) return(NULL);

   /* Fill in the data for the start N                                  */
   p = entry;
   strcpy(p->atnam,"N   ");
   p->x = (REAL)(-DIST_C_N);
   p->y = (REAL)0.0;
   p->z = (REAL)0.0;

   /* And for the C-alpha                                               */
   ALLOCNEXTPREV(p,ENTRY);
   if(p == NULL)
   {
      FREELIST(entry,ENTRY);
      return(NULL);
   }
   strcpy(p->atnam,"CA  ");
   p->x = (REAL)0.0;
   p->y = (REAL)0.0;
   p->z = (REAL)0.0;
   
   /* And for the C                                                     */
   ALLOCNEXTPREV(p,ENTRY);
   if(p == NULL)
   {
      FREELIST(entry,ENTRY);
      return(NULL);
   }
   strcpy(p->atnam,"C   ");
   p->x = (REAL)(DIST_CA_C * sin((double)(PI*ANGLE_CA/180.0)));
   p->y = (REAL)(DIST_CA_C * cos((double)(PI*ANGLE_CA/180.0)));
   p->z = (REAL)0.0;
   
   for(i=0; i<ntor; i++)
   {
      if(i!=0)
      {
         D("Built phi C atom\n");
         /* Build the C position using the phi data                     */
         if(BuildAtom(p->prev->prev, p->prev, p,   /* C, N, CA coords   */
                      ANGLE_CA, DIST_CA_C,         /* Constants         */
                      phi[i],                      /* Torsion angle     */
                      coords))                     /* Results           */
	 {
            ALLOCNEXTPREV(p,ENTRY);
            if(p == NULL)
            {
               FREELIST(entry,ENTRY);
               return(NULL);
            }
            strcpy(p->atnam,"C   ");
            p->x = coords[0];
            p->y = coords[1];
            p->z = coords[2];
	 }
      }

      if(i != (ntor-1))
      {
         /* Build the N position using the psi data                     */
         if(BuildAtom(p->prev->prev, p->prev, p,   /* N, CA, C coords   */
                      ANGLE_C, DIST_C_N,           /* Constants         */
                      psi[i],                      /* Torsion angle     */
                      coords))                     /* Results           */
	 {
            D("Built psi N atom\n");
            ALLOCNEXTPREV(p,ENTRY);
            if(p == NULL)
            {
               FREELIST(entry,ENTRY);
               return(NULL);
            }
            strcpy(p->atnam,"N   ");
            p->x = coords[0];
            p->y = coords[1];
            p->z = coords[2];
	 }

         /* Build the CA position using a 180 degree omega angle        */
         if(BuildAtom(p->prev->prev, p->prev, p,   /* CA, C, N coords   */
                      ANGLE_N, DIST_N_CA,          /* Constants         */
                      omega[i],                    /* Torsion angle     */
                      coords))                     /* Results           */
	 {
            D("Built omega CA atom\n");
            ALLOCNEXTPREV(p,ENTRY);
            if(p == NULL)
            {
               FREELIST(entry,ENTRY);
               return(NULL);
            }
            strcpy(p->atnam,"CA  ");
            p->x = coords[0];
            p->y = coords[1];
            p->z = coords[2];
	 }
      }
   }
   return(entry);
}

/************************************************************************/
/*>BOOL BuildAtom(ENTRY *p, ENTRY *q, ENTRY *r, REAL ang, REAL dist, 
                  REAL tor, REAL coords[3])
   -----------------------------------------------------------------
   Build coords for an atom given the coords of 3 anticedants, an angle,
   distance and torsion. Based on code from CARTX2 in CONGEN

   09.11.93 Original   By: ACRM
*/
#define ETA 1.0e-7
#define ETA2 ETA*ETA


BOOL BuildAtom(ENTRY *p, ENTRY *q, ENTRY *r, REAL theta, REAL bond, 
               REAL phi, REAL *coords)
{
   REAL stht, ctht,
        sphi, cphi,
        bsin,
        x1, y1, z1,
        x2, y2, z2,
        x3, y3, z3,
        x4, y4, z4,
        lyz1, ovlyz1,
        yy4, zz4,
        y1o, z1o,
        lxz22, l2, lxz2,
        ovl2, ovlxz2,
        x2o, z2o, xz2o, y2o, xx1, xx4;

   if(p==NULL || q==NULL || r==NULL) return(FALSE);

   stht  = (REAL)sin((double)(PI-theta));
   ctht  = (REAL)cos((double)(PI-theta));
   sphi  = (REAL)sin((double)phi);
   cphi  = (REAL)cos((double)phi);
   bsin  = bond * stht;

   x4    = bond * ctht;
   y4    = bsin * cphi;
   z4    = bsin * sphi;

   x3    = r->x;
   y3    = r->y;
   z3    = r->z;

   x1    = p->x - x3;
   y1    = p->y - y3;
   z1    = p->z - z3;

   x2    = q->x - x3;
   y2    = q->y - y3;
   z2    = q->z - z3;

   lxz22 = x2*x2 + z2*z2;
   l2    = (REAL)sqrt((double)(lxz22+y2*y2));
   lxz2  = (REAL)sqrt((double)lxz22);

   if(l2 < ETA)
   {
      fprintf(stderr,"Atoms 2 & 3 too close!\n");
      ovl2 = (REAL)1.0/ETA;
   }
   else
   {
      ovl2 = (REAL)1.0/l2;
   }

   if(lxz2 < ETA)
   {
      xx1 = x1;
      x2o = (REAL)1.0;
      z2o = (REAL)0.0;
   }
   else
   {
      ovlxz2 = (REAL)1.0/lxz2;
      x2o    = x2 * ovlxz2;
      z2o    = z2 * ovlxz2;
      xx1    = x1*x2o + z1*z2o;
      z1     = z1*x2o - x1*z2o;
   }

   xz2o   = lxz2 * ovl2;
   y2o    = y2   * ovl2;

   x1     = -xx1*xz2o - y1*y2o;
   y1     = xx1*y2o   - y1*xz2o;

   lyz1   = (REAL)sqrt((double)(y1*y1 + z1*z1));
   ovlyz1 = (REAL)1.0 / lyz1;

   y1o    = y1 * ovlyz1;
   z1o    = z1 * ovlyz1;

   yy4    = y1o*y4 - z1o*z4;
   zz4    = y1o*z4 + z1o*y4;
   xx4    = y2o*yy4 - xz2o*x4;

   y4     = -xz2o*yy4 - y2o*x4;
   x4     = x2o*xx4 - z2o*zz4;
   z4     = z2o*xx4 + x2o*zz4;

   coords[0] = x4 + x3;
   coords[1] = y4 + y3;
   coords[2] = z4 + z3;

   return(TRUE);
}

/************************************************************************/
/*>void WriteFragment(FILE *gOutfp, ENTRY *entry, BOOL calphas, 
                      BOOL SmallOut, char *seq)
   ------------------------------------------------------------
   Writes the atom names and coordinates from the ENTRY linked list to 
   the output file. If the calphas flag is set, only CA's will be written.
   
   09.11.93 Original   By: ACRM
   10.11.93 Added SmallOut option
   19.11.93 Added seq parameter; now write correct sequence in PDB type
            output.
   03.11.17 Initialize resnam and change %lf to %f
*/
void WriteFragment(FILE *gOutfp, ENTRY *entry, BOOL calphas, 
                   BOOL SmallOut, char *seq)
{
   ENTRY *p;
   int    atnum  = 1,
          resnum = 0;
   char  *resnam = NULL;

   for(p=entry; p!=NULL; NEXT(p))
   {
      if(SmallOut)
      {
         if(!calphas)
         {
            if(!strcmp(p->atnam,"N   "))
               fprintf(gOutfp,"N    %8.3f%8.3f%8.3f\n",p->x,p->y,p->z);
            if(!strcmp(p->atnam,"C   "))
               fprintf(gOutfp,"C    %8.3f%8.3f%8.3f\n",p->x,p->y,p->z);
         }
         if(!strcmp(p->atnam,"CA  "))
            fprintf(gOutfp,"CA   %8.3f%8.3f%8.3f\n",p->x,p->y,p->z);
      }
      else   /* PDB output */
      {
         if(!strcmp(p->atnam,"N   ")) 
         {
            resnam = blOnethr(seq[resnum]);
            resnum++;
	 }

         if(!calphas)
         {
            if(!strcmp(p->atnam,"N   "))
               fprintf(gOutfp,"ATOM  %5d  %-4s%-4s %4d    \
%8.3f%8.3f%8.3f\n",
                       atnum++,p->atnam,resnam,resnum,p->x,p->y,p->z);

            if(!strcmp(p->atnam,"C   "))
               fprintf(gOutfp,"ATOM  %5d  %-4s%-4s %4d    \
%8.3f%8.3f%8.3f\n",
                       atnum++,p->atnam,resnam,resnum,p->x,p->y,p->z);
         }
         if(!strcmp(p->atnam,"CA  "))
            fprintf(gOutfp,"ATOM  %5d  %-4s%-4s %4d    \
%8.3f%8.3f%8.3f\n",
                    atnum++,p->atnam,resnam,resnum,p->x,p->y,p->z);
      }
   }
}

/************************************************************************/
/*>char *GetSequence(char *buffer)
   -------------------------------
   Pulls the 1-letter code sequence out of a record from the hits file

   19.11.93 Original   By: ACRM
*/
char *GetSequence(char *buffer)
{
   char *ptr;

   /* Terminate the buffer at the )                                     */
   if((ptr = strchr(buffer,')')) != NULL) *ptr = '\0';

   /* Now find the (                                                    */
   ptr = strchr(buffer,'(');

   if(ptr == NULL) return(NULL);

   return(ptr+1);
}
