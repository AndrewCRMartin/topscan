/*************************************************************************

   Program:    topscan
   File:       topscan.c
   
   Version:    V1.2
   Date:       12.03.98
   Function:   Compare protein topologies
   
   Copyright:  (c) UCL, Dr. Andrew C. R. Martin 1998
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
   V1.0  13.01.98 Original
   V1.1  15.01.98 Fixed for where there is no secondary structure found
   V1.2  12.03.98 Added option to use STRIDE rather than DSSP

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
#define DSSP        "dssp"
#define STRIDE      "stride"
#define MERGESTRIDE "mergestride"
#define MAXBUFF     160
#define GAPPEN      8
#define MATFILE     "topmat.mat"
#define ELEN        4
#define HLEN        4

/************************************************************************/
/* Globals
*/
BOOL   gVerbose = FALSE;        /* Should we display alignments?        */
char   gBest1[MAXBUFF],         /* Used to store the best alignment     */
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
                  char *matfile, int *ELen, int *HLen, BOOL *RunDSSP,
                  BOOL *BuildOnly, BOOL *ScanMode, BOOL *UseBoth,
                  BOOL *UseStride);
void Usage(void);
char *ReadTopology(FILE *fp, int ELen, int HLen, BOOL UseStride);
char CalcElement(char struc, REAL x1, REAL y1, REAL z1, 
                 REAL x2, REAL y2, REAL z2);
int CalcIDScore(char *seq1, char *seq2, BOOL UseBoth);
char *ReadDSSP(FILE *fp, int ELen, int HLen);
char *ReadStride(FILE *fp, int ELen, int HLen);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for the topscan program

   13.01.98 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE  *fdssp1,
         *fdssp2;
   char  infile1[MAXBUFF],
         infile2[MAXBUFF],
         sourcefile[MAXBUFF],
         matfile[MAXBUFF],
         *top1,
         *top2;
   int   score, 
         IDScore, 
         ELen      = ELEN, 
         HLen      = HLEN;
   BOOL  RunDSSP   = FALSE,
         BuildOnly = FALSE,
         ScanMode  = FALSE,
         UseBoth   = FALSE,
         UseStride = FALSE;
#ifdef __linux__
   __pid_t pid;
#else
   pid_t pid;
#endif

   gBest1[0] = '\0';
   gBest2[0] = '\0';
   
   if(ParseCmdLine(argc, argv, infile1, infile2, matfile, &ELen, &HLen,
                   &RunDSSP, &BuildOnly, &ScanMode, &UseBoth, &UseStride))
   {
      strcpy(sourcefile,infile1);
      
      if(RunDSSP)
      {
         char ofile1[MAXBUFF],
              ofile2[MAXBUFF],
              cmd[MAXBUFF];
         
         pid = getpid();

         sprintf(ofile1,"/tmp/file1.%d",pid);
         if(UseStride)
         {
            sprintf(cmd,"%s %s | %s %s > %s", STRIDE, infile1, 
                    MERGESTRIDE, infile1, ofile1);
            system(cmd);
         }
         else
         {
            sprintf(cmd,"%s %s %s >/dev/null", DSSP, infile1, ofile1);
            system(cmd);
         }
         strcpy(infile1, ofile1);
         

         if(!BuildOnly && !ScanMode)
         {
            sprintf(ofile2,"/tmp/file2.%d",pid);
            if(UseStride)
            {
               sprintf(cmd,"%s %s | %s %s > %s", STRIDE, infile2, 
                       MERGESTRIDE, infile2, ofile2);
               system(cmd);
            }
            else
            {
               sprintf(cmd,"%s %s %s >/dev/null", DSSP, infile2, ofile2);
               system(cmd);
            }
            strcpy(infile2, ofile2);
         }
      }
      
      /* Open the DSSP files                                            */
      if((fdssp1=fopen(infile1,"r"))==NULL)
      {
         fprintf(stderr,"Can't read %s\n",infile1);
         return(1);
      }
      if(!BuildOnly)
      {
         if((fdssp2=fopen(infile2,"r"))==NULL)
         {
            fprintf(stderr,"Can't read %s\n",infile2);
            return(1);
         }
      }

      /* Read the DSSP files                                            */
      if((top1 = ReadTopology(fdssp1, ELen, HLen, UseStride))==NULL)
      {
         fprintf(stderr,"Unable to read topology from %s\n",infile1);
         return(1);
      }

      if(BuildOnly)
      {
         printf("%s %s\n",sourcefile,top1);
      }
      else
      {
         /* Read the Matrix file                                        */
         if(!ReadMDM(matfile))
         {
            fprintf(stderr,"Unable to read matrix file %s\n",matfile);
            return(1);
         }
            
         if(ScanMode)
         {
            char buffer[MAXBUFF], 
                 name[MAXBUFF],
                 *ptr;

            if((top2 = (char *)malloc(MAXBUFF * sizeof(char)))==NULL)
            {
               fprintf(stderr,"No memory for topology buffer\n");
               return(1);
            }
            
            while(fgets(buffer,MAXBUFF,fdssp2))
            {
               TERMINATE(buffer);
               
               ptr = buffer;
               while(*ptr == ' ' || *ptr == '\t')
                  ptr++;
               if(strlen(ptr) && (*ptr != '!') && (*ptr != '#'))
               {
                  name[0]   = '\0';
                  top2[0]   = '\0';
                  strcpy(gBest1, top1);
                  gBest2[0] = '\0';
                  
                  sscanf(ptr,"%s %s",name,top2);

                  IDScore = CalcIDScore(top1, top2, UseBoth);
            
                  if((score = RunAlignment(top1, top2))==(-1))
                     return(1);
            
                  /* Print the result                                   */
                  if(gVerbose)
                  {
                     printf("! %s\n! %s\n",gBest1,gBest2);
                  }
                  printf("%s %f\n", name,
                         (REAL)100.0 * (REAL)score / (REAL)IDScore);
                  
               }
            }
         }
         else
         {
            if((top2 = ReadTopology(fdssp2, ELen, HLen, UseStride))==NULL)
            {
               fprintf(stderr,"Unable to read topology from %s\n",
                       infile2);
               return(1);
            }
            
            IDScore = CalcIDScore(top1, top2, UseBoth);
            
            if((score = RunAlignment(top1, top2))==(-1))
               return(1);
            
            /* Print the result                                         */
            if(gVerbose)
            {
               printf("%s\n%s\n",gBest1,gBest2);
            }
            printf("%f\n",(REAL)100.0 * (REAL)score / (REAL)IDScore);
         }
      }
      
      if(RunDSSP)
      {
         unlink(infile1);
         if(!BuildOnly && !ScanMode)
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
   If both are of length 0, returns a score of 100. If only one is
   of length zero, returns 0

   13.01.98 Original   By: ACRM
   15.01.98 Added check for 0-length topology strings
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

   /* 15.01.98 Added this check on 0-length topology strings            */
   if((length1 == 0) && (length2 == 0))
      return(100);
   if((length1 == 0) || (length2 == 0))
      return(0);
   
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

   13.01.98 Original   By: ACRM
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

   13.01.98 Original   By: ACRM
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

   13.01.98 Original   By: ACRM
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
                     char *matfile, int *ELen, int *HLen, BOOL *RunDSSP,
                     BOOL *BuildOnly, BOOL *ScanMode, BOOL *UseBoth,
                     BOOL *UseStride)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile1     DSSP Input file 1
            char   *infile2     DSSP Input file 2
            char   *matfile     Matrix file      
            int    *ELen        Minimum strand length
            int    *HLen        Minimum helix length
            BOOL   *RunDSSP     Run the DSSP program on input files?
            BOOL   *BuildOnly   Just create the topology string for the
                                structure?
            BOOL   *ScanMode    Run in scan mode?
            BOOL   *UseBoth     Calculate the percentages wrt max score
                                from either struc rather than just the 
                                first
            BOOL   *UseStride   Use Stride rather than DSSP
   Returns: BOOL                Success?

   Parse the command line
   
   13.01.98 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile1, char *infile2, 
                  char *matfile, int *ELen, int *HLen, BOOL *RunDSSP,
                  BOOL *BuildOnly, BOOL *ScanMode, BOOL *UseBoth,
                  BOOL *UseStride)
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
            if(argv[0][2] == 's')
               *UseStride = TRUE;
            break;
         case 'b':
            *BuildOnly = TRUE;
            break;
         case 's':
            *ScanMode = TRUE;
            break;
         case 'w':
            *UseBoth = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 2 arguments left                       */
         if(*BuildOnly && argc != 1)
            return(FALSE);
         if(!(*BuildOnly) && argc != 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile1, argv[0]);
         if(!(*BuildOnly))
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

   13.01.98 Original   By: ACRM
   15.01.98 V1.1
   13.03.98 V1.2
*/
void Usage(void)
{
   fprintf(stderr,"\ntopscan V1.2 (c) 1998, Dr. Andrew C.R. Martin, \
UCL\n");

   fprintf(stderr,"\nUsage: topscan [-v] [-p[s]] [-w] [-h hlen] \
[-e elen] [-m matrix] file1.{dssp|pdb} file2.{dssp|pdb}\n");
   fprintf(stderr,"       topscan -b [-p[s]] [-h hlen] [-e elen] \
file1.{dssp|pdb}\n");
   fprintf(stderr,"       topscan -s [-v] [-p[s]] [-w] [-h hlen] \
[-e elen] [-m matrix] file1.{dssp|pdb} file2.top\n");
   fprintf(stderr,"\n       -b Build the topology string for a file\n");
   fprintf(stderr,"          (Don't actually run a comparison)\n");
   fprintf(stderr,"       -s Scan a DSSP or PDB file against a library \
of topology strings\n");
   fprintf(stderr,"          stored in the second file\n");
   fprintf(stderr,"       -m Specify the matrix file [Default: %s]\n",
           MATFILE);
   fprintf(stderr,"       -v Verbose mode\n");
   fprintf(stderr,"       -p Input files are PDB and the DSSP program \
will be run first\n");
   fprintf(stderr,"          If -ps is specified, then STRIDE will be \
run rather than DSSP\n");
   fprintf(stderr,"       -w Calculate score as percentage from both \
topology strings\n");
   fprintf(stderr,"          rather than just the first\n");
   fprintf(stderr,"       -h Specify minimum helix length [Default: \
%d]\n", HLEN);
   fprintf(stderr,"       -e Specify minimum strand length [Default: \
%d]\n", ELEN);

   fprintf(stderr,"\nCondenses a protein structure into a topology \
string by reading from\n");
   fprintf(stderr,"a DSSP or STRIDE file. The string has a 12-letter \
alphabet for the 6\n");
   fprintf(stderr,"orientations of sheet and helix. The comparison is \
performed using 24\n");
   fprintf(stderr,"orientations for one of the structures compared with \
the other in a fixed\n");
   fprintf(stderr,"position. Output is the best score obtained - the \
alignment is also\n");
   fprintf(stderr,"given if the verbose option is selected.\n");

   fprintf(stderr,"\nOutput is the best score obtained. The score is \
presented as a percentage\n");
   fprintf(stderr,"of the score obtained by aligning the first structure \
with itself; thus\n");
   fprintf(stderr,"the first structure is treated as a probe being \
tested against the \n");
   fprintf(stderr,"second structure. If the -w option is given, the \
maximum possible\n");
   fprintf(stderr,"score is calculated from both topology strings. \
If the verbose option\n");
   fprintf(stderr,"is selected, the alignment is also given.\n");

   fprintf(stderr,"\nIn scan mode (-s), the topology from the first \
file is scanned against\n");
   fprintf(stderr,"a library of topology strings stored in the second \
file. Entries for \n");
   fprintf(stderr,"the topology library file may be generated by using \
the program in\n");
   fprintf(stderr,"build mode (-b).\n\n");
}




/************************************************************************/
/*>char *ReadTopology(FILE *fp, int ELen, int HLen, BOOL UseStride)
   ----------------------------------------------------------------
   Input:   FILE   *fp        DSSP file pointer
            int    ELen       Minimum length of strand
            int    HLen       Minimum length of helix
            BOOL   UseStride  Read from Stride rather than DSSP file
   Returns: char *            Topology string

   Reads a the topology from a DSSP or Stride file, returning a string 
   representing the topology.

   13.03.98 Original   By: ACRM
*/
char *ReadTopology(FILE *fp, int ELen, int HLen, BOOL UseStride)
{
   if(UseStride)
   {
      return(ReadStride(fp, ELen, HLen));
   }
   else
   {
      return(ReadDSSP(fp, ELen, HLen));
   }
}


/************************************************************************/
/*>char *ReadDSSP(FILE *fp, int ELen, int HLen)
   --------------------------------------------
   Input:   FILE     *fp     DSSP file pointer
            int      ELen    Minimum length of strand
            int      HLen    Minimum length of helix
   Returns: char *           Topology string

   Reads a DSSP file returning a string representing the topology.

   13.01.98 Original   By: ACRM
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

   13.01.98 Original   By: ACRM
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
/*>int CalcIDScore(char *seq1, char *seq2, BOOL UseBoth)
   -----------------------------------------------------
   Input:   char     *seq1   Sequence 1
            char     *seq2   Sequence 2
            BOOL     UseBoth Calculate score as Max of both sequences
   Returns: int              Max score of each sequence vs itself

   Calculates the maximum possible score resulting from the identical
   sequence

   14.01.98 Original   By: ACRM
   
*/
int CalcIDScore(char *seq1, char *seq2, BOOL UseBoth)
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
   if(UseBoth)
   {
      for(i=0; i<seqlen2; i++)
      {
         if(isalpha(seq2[i]))
            score2 += CalcMDMScore(seq2[i], seq2[i]);
      }
      if(score2 > score1)
         score1 = score2;
   }

   /* 15.01.98 Added check for zero length strings                      */
   if(UseBoth)
   {
      if((seqlen1==0) && (seqlen2==0))
         score1 = 100;
   }
   else
   {
      if(seqlen1==0)
         score1 = 100;
   }
   
   return(score1);
}

/************************************************************************/
/*>char *ReadStride(FILE *fp, int ELen, int HLen)
   ----------------------------------------------
   Input:   FILE     *fp     Stride file pointer
            int      ELen    Minimum length of strand
            int      HLen    Minimum length of helix
   Returns: char *           Topology string

   Reads a Stride file returning a string representing the topology.

   13.03.98 Original   By: ACRM
*/
char *ReadStride(FILE *fp, int ELen, int HLen)
{
   char *top,
        buffer[MAXBUFF*2],
        struc,
        LastStruc    = ' ';
   REAL x,  y,  z,
        x1, y1, z1,
        xp, yp, zp;
   BOOL InElement    = FALSE;
   int  i            = 0,
        EleLength    = 0;

   
   if((top = (char *)malloc(MAXBUFF * sizeof(char)))==NULL)
      return(NULL);
   
   while(fgets(buffer,MAXBUFF*2,fp))
   {
      TERMINATE(buffer);
      fsscanf(buffer,"%15x%8lf%1x%8lf%1x%8lf%1x%c",&x,&y,&z,&struc);
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
   
   if(InElement)
   {
      /* File ended with an element                                     */
      top[i++] = CalcElement(LastStruc,x1,y1,z1,xp,yp,zp);
   }
   top[i] = '\0';
   
   return(top);
}

