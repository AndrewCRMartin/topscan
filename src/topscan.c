/*************************************************************************

   Program:    
   File:       
   
   Version:    
   Date:       
   Function:   
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1995
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 419 3890
   EMail:      INTERNET: martin@biochem.ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "bioplib/general.h"
#include "bioplib/seq.h"
#include "bioplib/macros.h"
#include "bioplib/fsscanf.h"

/************************************************************************/
/* Defines and macros
*/
#define DSSP    "dssp"
#define MAXBUFF 160
#define GAPPEN  8
#define MATFILE "topmat.mat"
#define ELEN    4
#define HLEN    4

/************************************************************************/
/* Globals
*/
BOOL   gVerbose = FALSE;
char   gBest1[MAXBUFF],
       gBest2[MAXBUFF];

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
int RunAlignment(char *top1, char *top2);
void TurnAboutX(char *top);
void TurnAboutY(char *top);
void TurnAboutZ(char *top);
BOOL ParseCmdLine(int argc, char **argv, char *infile1, char *infile2, 
                  char *matfile, int *ELen, int *HLen, BOOL *RunDSSP);
void Usage(void);
char *ReadDSSP(FILE *fp, int ELen, int HLen);
char CalcElement(char struc, REAL x1, REAL y1, REAL z1, 
                 REAL x2, REAL y2, REAL z2);
int CalcIDScore(char *seq1, char *seq2);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for the topscan program

   13.01.97 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE  *fdssp1,
         *fdssp2;
   char  infile1[MAXBUFF],
         infile2[MAXBUFF],
         matfile[MAXBUFF],
         *top1,
         *top2;
   int   score, 
         IDScore, 
         ELen = ELEN, 
         HLen = HLEN;
   BOOL  RunDSSP = FALSE;
#ifdef __linux__
   __pid_t pid;
#else
   pid_t pid;
#endif
   
   if(ParseCmdLine(argc, argv, infile1, infile2, matfile, &ELen, &HLen,
                   &RunDSSP))
   {
      if(RunDSSP)
      {
         char ofile1[MAXBUFF],
              ofile2[MAXBUFF],
              cmd[MAXBUFF];
         
         pid = getpid();
         sprintf(ofile1,"/tmp/file1.%d",pid);
         sprintf(ofile2,"/tmp/file2.%d",pid);
         sprintf(cmd,"%s %s %s", DSSP, infile1, ofile1);
         system(cmd);
         sprintf(cmd,"%s %s %s", DSSP, infile2, ofile2);
         system(cmd);
         strcpy(infile1, ofile1);
         strcpy(infile2, ofile2);
      }
      
      /* Open the DSSP files                                            */
      if((fdssp1=fopen(infile1,"r"))==NULL)
      {
         fprintf(stderr,"Can't read %s\n",infile1);
         return(1);
      }
      if((fdssp2=fopen(infile2,"r"))==NULL)
      {
         fprintf(stderr,"Can't read %s\n",infile2);
         return(1);
      }

      /* Read the DSSP files                                            */
      if((top1 = ReadDSSP(fdssp1, ELen, HLen))==NULL)
      {
         fprintf(stderr,"Unable to read topology from %s\n",infile1);
         return(1);
      }
      if((top2 = ReadDSSP(fdssp2, ELen, HLen))==NULL)
      {
         fprintf(stderr,"Unable to read topology from %s\n",infile2);
         return(1);
      }

      /* Read the Matrix file                                           */
      if(!ReadMDM(matfile))
      {
         fprintf(stderr,"Unable to read matrix file %s\n",matfile);
         return(1);
      }

      IDScore = CalcIDScore(top1, top2);

      if((score = RunAlignment(top1, top2))==(-1))
         return(1);
      
      /* Print the result                                               */
      if(gVerbose)
      {
         printf("%s\n%s\n",gBest1,gBest2);
      }
      printf("%f\n",(REAL)100.0 * (REAL)score / (REAL)IDScore);

      if(RunDSSP)
      {
         unlink(infile1);
         unlink(infile2);
      }
   }
   else
   {
      Usage();
   }
   

   return(0);
}


/************************************************************************/
/*>int RunAlignment(char *top1, char *top2)
   ----------------------------------------
   Input:   char     *top1   First topology string
            char     *top2   Second topology string
   Returns: int              Alignment score
   Globals: char     *gBest1 Best alignment of top1 (maybe rotated)
                     *gBest2 Best alignment of top2

   Does the alignment in all 24 rotations

   13.01.97 Original   By: ACRM
*/
int RunAlignment(char *top1, char *top2)
{
   int  length1,
        length2,
        score,
        maxscore,
        i, j,
        align_len;
   char *align1,
        *align2;
   
   
   /* Allocate memory to store the alignment and run the N&W            */
   length1 = strlen(top1);
   length2 = strlen(top2);
   if((align1 = (char *)malloc((length1+length2)*sizeof(char)))==NULL)
   {
      fprintf(stderr,"No memory for alignment1\n");
      return(-1);
   }
   if((align2 = (char *)malloc((length1+length2)*sizeof(char)))==NULL)
   {
      fprintf(stderr,"No memory for alignment2\n");
      return(-1);
   }

   /* Native position                                                   */
   maxscore = align(top1,length1,top2,length2,FALSE,
                    FALSE,GAPPEN,align1,align2,&align_len);
   strcpy(gBest1,align1);
   strcpy(gBest2,align2);
   gBest1[align_len] = '\0';
   gBest2[align_len] = '\0';
   
   for(i=0; i<4; i++)
   {
      TurnAboutX(top1);
      for(j=0; j<4; j++)
      {
         TurnAboutZ(top1);
         
         score = align(top1,length1,top2,length2,FALSE,
                       FALSE,GAPPEN,align1,align2,&align_len);
         if(score > maxscore)
         {
            maxscore = score;
            strcpy(gBest1,align1);
            strcpy(gBest2,align2);
            gBest1[align_len] = '\0';
            gBest2[align_len] = '\0';
         }
      }
   }
   for(i=0; i<2; i++)
   {
      TurnAboutY(top1);
      for(j=0; j<4; j++)
      {
         TurnAboutZ(top1);
         
         score = align(top1,length1,top2,length2,FALSE,
                       FALSE,GAPPEN,align1,align2,&align_len);
         if(score > maxscore)
         {
            maxscore = score;
            strcpy(gBest1,align1);
            strcpy(gBest2,align2);
            gBest1[align_len] = '\0';
            gBest2[align_len] = '\0';
         }
      }
      if(i==0) TurnAboutY(top1);
   }

   return(maxscore);
}


/************************************************************************/
/*>void TurnAboutX(char *top)
   --------------------------
   I/O:     char     *top    Topology string to rotate

   Modifies the topology string by rotating about X

   13.01.97 Original   By: ACRM
*/
void TurnAboutX(char *top)
{
   static char *old = "ABCDEFGHIJKL";
   static char *new = "EBFDCAKHLJIG";
   
   while(*top)
   {
      *top = new[chindex(old,*top)];
      top++;
   }
}


/************************************************************************/
/*>void TurnAboutY(char *top)
   --------------------------
   I/O:     char     *top    Topology string to rotate

   Modifies the topology string by rotating about Y

   13.01.97 Original   By: ACRM
*/
void TurnAboutY(char *top)
{
   static char *old = "ABCDEFGHIJKL";
   static char *new = "AECFDBGKILJH";
   
   while(*top)
   {
      *top = new[chindex(old,*top)];
      top++;
   }
}


/************************************************************************/
/*>void TurnAboutZ(char *top)
   --------------------------
   I/O:     char     *top    Topology string to rotate

   Modifies the topology string by rotating about Z

   13.01.97 Original   By: ACRM
*/
void TurnAboutZ(char *top)
{
   static char *old = "ABCDEFGHIJKL";
   static char *new = "BCDAEFHIJGKL";
   
   while(*top)
   {
      *top = new[chindex(old,*top)];
      top++;
   }
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile1, char *infile2, 
                     char *matfile, int *ELen, int *HLen, BOOL *RunDSSP)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile1     DSSP Input file 1
            char   *infile2     DSSP Input file 2
            char   *matfile     Matrix file      
            int    *ELen        Minimum strand length
            int    *HLen        Minimum helix length
            BOOL   *RunDSSP     Run the DSSP program on input files
   Returns: BOOL                Success?

   Parse the command line
   
   13.01.97 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile1, char *infile2, 
                  char *matfile, int *ELen, int *HLen, BOOL *RunDSSP)
{
   argc--;
   argv++;

   infile1[0] = infile2[0] = '\0';
   strcpy(matfile,MATFILE);

   if(!argc)
      return(FALSE);
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'm':
            argc--;
            argv++;
            if(argc>0)
               strcpy(matfile,argv[0]);
            break;
         case 'h':
            argc--;
            argv++;
            if(argc>0)
               sscanf(argv[0],"%d",HLen);
            break;
         case 'e':
            argc--;
            argv++;
            if(argc>0)
               sscanf(argv[0],"%d",ELen);
            break;
         case 'v':
            gVerbose = TRUE;
            break;
         case 'p':
            *RunDSSP = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 2 arguments left                       */
         if(argc != 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile1, argv[0]);
         strcpy(infile2, argv[1]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   13.01.97 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\ntopscan V1.0 (c) Dr. Andrew C.R. Martin, UCL\n");

   fprintf(stderr,"Usage: topscan [-m matrix] [-p] [-v] [-h hlen] \
[-e elen] file1.dssp file2.dssp\n");
   fprintf(stderr,"       -m Specify the matrix file [Default: %s]\n",
           MATFILE);
   fprintf(stderr,"       -v Verbose mode\n");
   fprintf(stderr,"       -p Input files are PDB and the DSSP program \
will be run first\n");
   fprintf(stderr,"       -h Specify minimum helix length [Default: \
%d]\n", HLEN);
   fprintf(stderr,"       -e Specify minimum strand length [Default: \
%d]\n", ELEN);

   fprintf(stderr,"\nCondenses a protein structure into a topology \
string by reading from\n");
   fprintf(stderr,"a DSSP file. The string has a 12-letter alphabet for \
the 6 orientations\n");
   fprintf(stderr,"of sheet and helix. The comparison is performed using \
24 orientations\n");
   fprintf(stderr,"for one of the structures compared with the other in \
a fixed position.\n");
   fprintf(stderr,"Output is the best score obtained - the alignment is \
also given if the\n");
   fprintf(stderr,"verbose option is selected\n\n");
}


/************************************************************************/
/*>char *ReadDSSP(FILE *fp, int ELen, int HLen)
   --------------------------------------------
   Input:   FILE     *fp     DSSP file pointer
            int      ELen    Minimum length of strand
            int      HLen    Minimum length of helix
   Returns: char *           Topology string

   Reads a DSSP file returning a string representing the topology.

   13.01.97 Original   By: ACRM
*/
char *ReadDSSP(FILE *fp, int ELen, int HLen)
{
   char *top,
        buffer[MAXBUFF*2],
        struc,
        LastStruc    = ' ';
   REAL x,  y,  z,
        x1, y1, z1,
        xp, yp, zp;
   BOOL InBody       = FALSE,
        InElement    = FALSE;
   int  i            = 0,
        EleLength    = 0;

   
   if((top = (char *)malloc(MAXBUFF * sizeof(char)))==NULL)
      return(NULL);
   
   while(fgets(buffer,MAXBUFF*2,fp))
   {
      TERMINATE(buffer);
      if(InBody)
      {
         fsscanf(buffer,"%16x%c%90x%7lf%7lf%7lf",&struc,&x,&y,&z);
         if((struc=='E' || struc=='H') && (struc==LastStruc))
         {
            EleLength++;
         }
         
         if((struc=='E' || struc=='H') && (struc!=LastStruc))
         {
            /* Start of new element                                     */
            if(InElement)
            {
               if((LastStruc=='E' && EleLength>=ELen) ||
                  (LastStruc=='H' && EleLength>=HLen))
                  top[i++] = CalcElement(LastStruc,x1,y1,z1,xp,yp,zp);
            }
            
            InElement = TRUE;
            EleLength = 1;
            x1 = x;
            y1 = y;
            z1 = z;
         }
         else if(InElement && struc!='E' && struc!='H')
         {
            /* Just come out of an element                              */
            if((LastStruc=='E' && EleLength>=ELen) ||
               (LastStruc=='H' && EleLength>=HLen))
               top[i++] = CalcElement(LastStruc,x1,y1,z1,xp,yp,zp);
            InElement = FALSE;
         }
            
         LastStruc = struc;
         xp = x;
         yp = y;
         zp = z;
      }
      else
      {
         if(!strncmp(buffer, "  #",3))
         {
            InBody = TRUE;
         }
      }
   }

   if(InElement)
   {
      /* File ended with an element                                     */
      top[i++] = CalcElement(LastStruc,x1,y1,z1,xp,yp,zp);
   }
   top[i] = '\0';
   
   return(top);
}


/************************************************************************/
/*>char CalcElement(char struc, REAL x1, REAL y1, REAL z1, 
                    REAL x2, REAL y2, REAL z2)
   -------------------------------------------------------
   Input:   char     struc   DSSP structure assignment
            REAL     x1      Coordinates of SS element start
            REAL     y1      
            REAL     z1      
            REAL     x2      Coordinates of SS element end
            REAL     y2      
            REAL     z2      
   Returns: char             Code for element (structure & direction)

   Given a structure code from DSSP (E or H) and the coordinates of the
   ends of the element, returns a combined code of structure plus
   direction.

   13.01.97 Original   By: ACRM
*/
char CalcElement(char struc, REAL x1, REAL y1, REAL z1, 
                 REAL x2, REAL y2, REAL z2)
{
   char code;
   int  dirn;
   REAL dx,    dy,    dz,
        absdx, absdy, absdz;
   
   /* Base coding is A for sheet, G for helix                           */
   if(struc == 'E')
      code = 'A';
   else if(struc == 'H')
      code = 'G';
   else
      return(' ');

   /* Calculate deltas along x, y and z                                 */
   dx = x2 - x1;
   dy = y2 - y1;
   dz = z2 - z1;

   /* Find absolute values                                              */
   absdx = ABS(dx);
   absdy = ABS(dy);
   absdz = ABS(dz);

   /* Find the direction modifier                                       */
   if((dx >= absdy) && (dx >= absdz))
      dirn = 1;
   else if((dx <= (-absdy)) && (dx <= (-absdz)))
      dirn = 3;
   else if((dy >= absdx) && (dy >= absdz))
      dirn = 0;
   else if((dy <= (-absdx)) && (dy <= (-absdz)))
      dirn = 2;
   else if((dz >= absdx) && (dz >= absdy))
      dirn = 4;
   else if((dz <= (-absdx)) && (dz <= (-absdy)))
      dirn = 5;
   else
      return(' ');
   
   return(code + (char)dirn);
}

/************************************************************************/
/*>int CalcIDScore(char *seq1, char *seq2)
   ---------------------------------------
   Input:   char     *seq1   First sequence
            char     *seq2   Second sequence
   Returns: int              Max score of each sequence vs itself

   Calculates the maximum possible score resulting from the identical
   sequence

   (Lifted from nw.c)

   11.07.96 Original   By: ACRM
*/
int CalcIDScore(char *seq1, char *seq2)
{
   int score1 = 0,
       score2 = 0,
       i,
       seqlen1 = strlen(seq1),
       seqlen2 = strlen(seq2);

   for(i=0; i<seqlen1; i++)
   {
      if(isalpha(seq1[i]))
         score1 += CalcMDMScore(seq1[i], seq1[i]);
   }
   for(i=0; i<seqlen2; i++)
   {
      if(isalpha(seq2[i]))
         score2 += CalcMDMScore(seq2[i], seq2[i]);
   }
   
   return(MIN(score1, score2));
}


