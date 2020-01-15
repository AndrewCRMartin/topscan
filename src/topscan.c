/*************************************************************************

   Program:    topscan
   File:       topscan.c
   
   Version:    V2.0
   Date:       13.03.2000
   Function:   Compare protein topologies
   
   Copyright:  (c) UCL, Reading, Dr. Andrew C. R. Martin 1998-2000
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 419 3890
   EMail:      INTERNET: A.C.R.Martin@reading.ac.uk
               
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
   V1.3  26.10.99 Added -g option (treat 3_10 helix as alpha helix)
   V1.4  10.11.99 Added -1 option (primary topology only)
   V1.5  17.11.99 Added -n option (include neighbour information)
   V1.6  23.11.99 Added -a option (include accessibility)
   V1.7  19.01.00 Added -t option (give topology strings on command line)
   V1.8  26.01.00 Added -l option (include length)
   V2.0  13.03.00 Modified to work with numeric topology strings so we can
                  have a much larger alphabet

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
#include "bioplib/MathUtil.h"

/************************************************************************/
/* Defines and macros
*/
#define DSSP                  "dssp"
#define STRIDE                "stride"
#define MERGESTRIDE           "mergestride"
#define MAXBUFF               320
#define GAPPEN                8
#define MATFILE               "numtopmat.mat"
#define ELEN                  4
#define HLEN                  4
#define MARKER                -9999.0
#define ADJACENT_DIST         12.0
#define BUFFCHUNK             24

#define STRAND_MEAN_ACCESS_3  33.836
#define STRAND_MEAN_ACCESS_4  32.394
#define HELIX_MEAN_ACCESS_3   52.807
#define HELIX_MEAN_ACCESS_3G  56.547
#define HELIX_MEAN_ACCESS_4   52.795
#define HELIX_MEAN_ACCESS_4G  53.304

#define HELIX_MEAN_LENGTH     12.5
#define STRAND_MEAN_LENGTH    5.4
/************************************************************************/
/* Globals
*/
BOOL   gVerbose = FALSE;        /* Should we display alignments?        */
char   gBest1[MAXBUFF],         /* Used to store the best alignment     */
       gBest2[MAXBUFF];
REAL   gHelixMeanAccess = 0.0,  /* Mean accessibilities                 */
       gStrandMeanAccess = 0.0;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
int RunAlignment(int *top1, int *top2, BOOL PrimaryTopology);
void TurnAboutX(int *top);
void TurnAboutY(int *top);
void TurnAboutZ(int *top);
BOOL ParseCmdLine(int argc, char **argv, char *infile1, char *infile2, 
                  char *matfile, int *ELen, int *HLen, BOOL *RunDSSP,
                  BOOL *BuildOnly, BOOL *ScanMode, BOOL *UseBoth,
                  BOOL *UseStride, BOOL *Do3_10, BOOL *PrimaryTopology,
                  BOOL *DoNeighbour, BOOL *DoAccess, BOOL *GivenTopString,
                  BOOL *DoLength);
void Usage(void);
int *ReadTopology(FILE *fp, int ELen, int HLen, BOOL UseStride,
                   BOOL Do3_10, BOOL PrimaryTopology, BOOL DoNeighbour,
                   BOOL DoAccess, BOOL DoLength);
int CalcElement(char struc, REAL x1, REAL y1, REAL z1, 
                REAL x2, REAL y2, REAL z2, BOOL PrimaryTopology,
                BOOL DoNeighbour, BOOL DoAccess, REAL meanAccess,
                int  EleLength);
int CalcIDScore(int *seq1, int *seq2, BOOL UseBoth);
int *ReadDSSP(FILE *fp, int ELen, int HLen, BOOL Do3_10, 
              BOOL PrimaryTopology, BOOL DoNeighbour, BOOL DoAccess,
              BOOL DoLength);
int *ReadStride(FILE *fp, int ELen, int HLen, BOOL Do3_10, 
                BOOL PrimaryTopology, BOOL DoNeighbour, BOOL DoAccess,
                BOOL DoLength);
BOOL IsNeighbour(REAL x1, REAL y1, REAL z1,
                 REAL x2, REAL y2, REAL z2,
                 REAL prevx1, REAL prevy1, REAL prevz1,
                 REAL prevx2, REAL prevy2, REAL prevz2);
int FindArrayLength(int *array);
char *NumArrayToString(int *numarr);
int MakeIntArray(int *array1, char *inarray);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for the topscan program

   13.01.98 Original   By: ACRM
   26.10.99 Added 3_10 helix handling
   19.01.00 Added -t handling
   10.03.00 Changed to use integer coded topology array
   13.03.00 Initialise fdssp1, fdssp2 only to silence warnings with -O2
*/
int main(int argc, char **argv)
{
   FILE  *fdssp1 = NULL,
         *fdssp2 = NULL;
   char  infile1[MAXBUFF],
         infile2[MAXBUFF],
         sourcefile[MAXBUFF],
         matfile[MAXBUFF],
         top2str[MAXBUFF];
   int   *top1 = NULL,
         *top2 = NULL;
   int   score, 
         IDScore, 
         ELen            = ELEN, 
         HLen            = HLEN;
   BOOL  RunDSSP         = FALSE,
         BuildOnly       = FALSE,
         ScanMode        = FALSE,
         UseBoth         = FALSE,
         UseStride       = FALSE,
         Do3_10          = FALSE,
         PrimaryTopology = FALSE,
         DoNeighbour     = FALSE,
         DoAccess        = FALSE,
         DoLength        = FALSE,
         GivenTopString  = FALSE;
#ifdef __linux__
   __pid_t pid;
#else
   pid_t pid;
#endif

   gBest1[0] = '\0';
   gBest2[0] = '\0';
   
   if(ParseCmdLine(argc, argv, infile1, infile2, matfile, &ELen, &HLen,
                   &RunDSSP, &BuildOnly, &ScanMode, &UseBoth, &UseStride,
                   &Do3_10, &PrimaryTopology, &DoNeighbour, &DoAccess,
                   &GivenTopString, &DoLength))
   {
      if(GivenTopString)
      {
         if((top1 = (int *)malloc((1+strlen(infile1)) * sizeof(int)))
            ==NULL)
         {
            fprintf(stderr,"No memory for copying topology string\n");
            return(1);
         }
         MakeIntArray(top1, infile1);
         
         if(!ScanMode)
         {
            if((top2 = (int *)malloc((1+strlen(infile2)) * sizeof(int)))
               ==NULL)
            {
               fprintf(stderr,"No memory for copying topology string\n");
               return(1);
            }
            MakeIntArray(top2, infile2);
         }
      }
      else
      {
         strcpy(sourcefile,infile1);
         
         /* Choose which mean accessibilities to use if we are doing
            accessibilities
         */
         if(DoAccess)
         {
            BOOL Warn = FALSE;
            
            if(ELen == 3)
            {
               gStrandMeanAccess = STRAND_MEAN_ACCESS_3;
            }
            else
            {
               gStrandMeanAccess = STRAND_MEAN_ACCESS_4;
               if(ELen != 4) Warn = TRUE;
            }
            
            if(HLen == 3)
            {
               gHelixMeanAccess = ((Do3_10)?HELIX_MEAN_ACCESS_3G:
                                   HELIX_MEAN_ACCESS_3);
            }
            else
            {
               gHelixMeanAccess = ((Do3_10)?HELIX_MEAN_ACCESS_4G:
                                   HELIX_MEAN_ACCESS_4);
               if(ELen != 4) Warn = TRUE;
            }
            
            if(Warn)
            {
               fprintf(stderr,"Mean accessibilities not known for \
specified secondary structure length\nUsing values for length 4\n");
            }
         }
         
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
                  sprintf(cmd,"%s %s %s >/dev/null", 
                          DSSP, infile2, ofile2);
                  system(cmd);
               }
               strcpy(infile2, ofile2);
            }
         }
         
         /* Open the DSSP files                                         */
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
         
         /* Read the DSSP files                                         */
         if((top1 = ReadTopology(fdssp1, ELen, HLen, UseStride, Do3_10,
                                 PrimaryTopology, DoNeighbour,
                                 DoAccess, DoLength))==NULL)
         {
            fprintf(stderr,"Unable to read topology from %s\n",infile1);
            return(1);
         }
         
         if(BuildOnly)
         {
            char *ts = NumArrayToString(top1);
            if(ts==NULL)
            {
               fprintf(stderr,"No memory to create topology string from \
numeric array\n");
               return(1);
            }
            
            printf("%s %s\n",sourcefile,ts);
            free(ts);
         }
      }
      
      if(!BuildOnly)
      {
         /* Read the Matrix file                                     */
         if(!NumericReadMDM(matfile))
         {
            fprintf(stderr,"Unable to read matrix file %s\n",matfile);
            return(1);
         }
      
         if(ScanMode)
         {
            char buffer[MAXBUFF], 
                 name[MAXBUFF],
                 *ptr;

            if((top2 = (int *)malloc(MAXBUFF * sizeof(int)))==NULL)
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
                  char *ts;
                  
                  name[0]   = '\0';
                  top2[0]   = -1;
                  ts = NumArrayToString(top1);
                  strcpy(gBest1, ts);
                  free(ts);
                  
                  gBest2[0]  = '\0';
                  name[0]    = '\0';
                  top2str[0] = '\0';
                  
                  sscanf(ptr,"%s %s",name,top2str);
                  MakeIntArray(top2, top2str);

                  IDScore = CalcIDScore(top1, top2, UseBoth);
            
                  if((score = RunAlignment(top1, top2, 
                                           PrimaryTopology))==(-1))
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
            if(top2==NULL)
            {
               if((top2 = ReadTopology(fdssp2, ELen, HLen, UseStride, 
                                       Do3_10, PrimaryTopology, 
                                       DoNeighbour, DoAccess,
                                       DoLength))==NULL)
               {
                  fprintf(stderr,"Unable to read topology from %s\n",
                          infile2);
                  return(1);
               }
            }
            
            IDScore = CalcIDScore(top1, top2, UseBoth);
            
            if((score = RunAlignment(top1, top2, PrimaryTopology))==(-1))
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
/*>int RunAlignment(int *top1, int *top2, BOOL PrimaryTopology)
   ------------------------------------------------------------
   Input:   int      *top1           First topology string
            int      *top2           Second topology string
            BOOL     PrimaryTopology Primary topology only
   Returns: int                      Alignment score
   Globals: char     *gBest1         Best alignment of top1 (maybe 
                                     rotated)
                     *gBest2         Best alignment of top2

   Does the alignment in all 24 rotations
   If both are of length 0, returns a score of 100. If only one is
   of length zero, returns 0
   If PrimaryTopology is set, then only does the raw strings since
   no directions are encoded

   13.01.98 Original   By: ACRM
   15.01.98 Added check for 0-length topology strings
   10.11.99 Added PrimaryTopology
   08.03.00 Changed calls to align() to NumericAffineAlign()
            top1 and top2 now integer arrays
*/
int RunAlignment(int *top1, int *top2, BOOL PrimaryTopology)
{
   int  length1,
        length2,
        score,
        maxscore,
        i, j,
        align_len;
   int  *align1,
        *align2;
   char *ts;
   
   
   /* Allocate memory to store the alignment and run the N&W            */
   length1 = FindArrayLength(top1);
   length2 = FindArrayLength(top2);

   /* 15.01.98 Added this check on 0-length topology strings            */
   if((length1 == 0) && (length2 == 0))
      return(100);
   if((length1 == 0) || (length2 == 0))
      return(0);
   
   if((align1 = (int *)malloc((length1+length2)*sizeof(int)))==NULL)
   {
      fprintf(stderr,"No memory for alignment1\n");
      return(-1);
   }
   if((align2 = (int *)malloc((length1+length2)*sizeof(int)))==NULL)
   {
      fprintf(stderr,"No memory for alignment2\n");
      return(-1);
   }

   /* Native position                                                   */
   maxscore = NumericAffineAlign(top1,length1,top2,length2,FALSE,
                    FALSE,GAPPEN,0,align1,align2,&align_len);
   align1[align_len] = (-1);
   align2[align_len] = (-1);

   ts = NumArrayToString(align1);
   strcpy(gBest1, ts);
   free(ts);
   
   ts = NumArrayToString(align2);
   strcpy(gBest2, ts);
   free(ts);
   
   /* If we aren't doing direction information then we don't need to do
      the permutations of the string for different orientations
   */
   if(PrimaryTopology)
      return(maxscore);
   
   for(i=0; i<4; i++)
   {
      TurnAboutX(top1);
      for(j=0; j<4; j++)
      {
         TurnAboutZ(top1);
         
         score = NumericAffineAlign(top1,length1,top2,length2,FALSE,
                       FALSE,GAPPEN,0,align1,align2,&align_len);
         if(score > maxscore)
         {
            maxscore = score;

            align1[align_len] = (-1);
            align2[align_len] = (-1);
            
            ts = NumArrayToString(align1);
            strcpy(gBest1, ts);
            free(ts);
            
            ts = NumArrayToString(align2);
            strcpy(gBest2, ts);
            free(ts);
         }
      }
   }
   for(i=0; i<2; i++)
   {
      TurnAboutY(top1);
      for(j=0; j<4; j++)
      {
         TurnAboutZ(top1);
         
         score = NumericAffineAlign(top1,length1,top2,length2,FALSE,
                       FALSE,GAPPEN,0,align1,align2,&align_len);
         if(score > maxscore)
         {
            maxscore = score;

            align1[align_len] = (-1);
            align2[align_len] = (-1);
            
            ts = NumArrayToString(align1);
            strcpy(gBest1, ts);
            free(ts);
            
            ts = NumArrayToString(align2);
            strcpy(gBest2, ts);
            free(ts);
         }
      }
      if(i==0) TurnAboutY(top1);
   }

   return(maxscore);
}


/************************************************************************/
/*>void TurnAboutX(int *top)
   -------------------------
   I/O:     int      *top    Topology string to rotate

   Modifies the topology string by rotating about X

   13.01.98 Original   By: ACRM
   17.11.99 Added characters M-X for adjacent secondary structures
            Modified so array is indexed by the alphabet so it no
            longer needs to use chindex() which is painfully slow
   23.11.99 Added lower case letters
   10.03.00 Changed to use integer coded topology array
*/
void TurnAboutX(int *top)
{
   static int new[] = {0,5,2,6,4,3,1};
   
   while(*top >= 0)
   {
      *top = (6*(int)((*top - 1)/6)) + new[1+((*top - 1)%6)];
      top++;
   }
}


/************************************************************************/
/*>void TurnAboutY(int *top)
   --------------------------
   I/O:     int     *top    Topology string to rotate

   Modifies the topology string by rotating about Y

   13.01.98 Original   By: ACRM
   17.11.99 Added characters M-X for adjacent secondary structures
            Modified so array is indexed by the alphabet so it no
            longer needs to use chindex() which is painfully slow
   23.11.99 Added lower case letters
   10.03.00 Changed to use integer coded topology array
*/
void TurnAboutY(int *top)
{
   static int new[] = {0,1,5,3,6,4,2};
   
   while(*top >= 0)
   {
      *top = (6*(int)((*top - 1)/6)) + new[1+((*top - 1)%6)];
      top++;
   }
}


/************************************************************************/
/*>void TurnAboutZ(int *top)
   --------------------------
   I/O:     int     *top    Topology string to rotate

   Modifies the topology string by rotating about Z

   13.01.98 Original   By: ACRM
   17.11.99 Added characters M-X for adjacent secondary structures
            Modified so array is indexed by the alphabet so it no
            longer needs to use chindex() which is painfully slow
   23.11.99 Added lower case letters
   10.03.00 Changed to use integer coded topology array
*/
void TurnAboutZ(int *top)
{
   static int new[] = {0,2,3,4,1,5,6};
   
   while(*top >= 0)
   {
      *top = (6*(int)((*top - 1)/6)) + new[1+((*top - 1)%6)];
      top++;
   }
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile1, char *infile2, 
                     char *matfile, int *ELen, int *HLen, BOOL *RunDSSP,
                     BOOL *BuildOnly, BOOL *ScanMode, BOOL *UseBoth,
                     BOOL *UseStride, BOOL *Do3_10, BOOL *PrimaryTopology,
                     BOOL *DoNeighbour, BOOL *DoAccess, 
                     BOOL *GivenTopString, BOOL *DoLength)
   -----------------------------------------------------------------------
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
            BOOL   *Do3_10      Merge 3_10 helix with alpha helix
            BOOL   *PrimaryTopology Do only primary topology
            BOOL   *DoNeighbour Add neighbour information
            BOOL   *DoAccess    Add accessibility information
            BOOL   *GivenTopString  Given the topology string on the 
                                command line instead of a file
            BOOL   *DoLength    Add length information
   Returns: BOOL                Success?

   Parse the command line
   
   13.01.98 Original    By: ACRM
   26.10.99 Added -g / Do3_10
   10.11.99 Added -1 / PrimaryTopology
   16.11.99 Added DoNeighbour
   23.11.99 Added DoAccess
   19.01.00 Added -t (GivenTopString)
   26.01.00 Added -l (DoLength)
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile1, char *infile2, 
                  char *matfile, int *ELen, int *HLen, BOOL *RunDSSP,
                  BOOL *BuildOnly, BOOL *ScanMode, BOOL *UseBoth,
                  BOOL *UseStride, BOOL *Do3_10, BOOL *PrimaryTopology,
                  BOOL *DoNeighbour, BOOL *DoAccess,
                  BOOL *GivenTopString, BOOL *DoLength)
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
         case 'g':
            *Do3_10 = TRUE;
            break;
         case '1':
            *PrimaryTopology = TRUE;
            break;
         case 'n':
            *DoNeighbour = TRUE;
            break;
         case 'a':
            *DoAccess = TRUE;
            break;
         case 't':
            *GivenTopString = TRUE;
            break;
         case 'l':
            *DoLength = TRUE;
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
   26.10.99 V1.3
   10.11.99 V1.4
   16.11.99 V1.5
   19.11.99 V1.7
   13.03.00 V2.0
*/
void Usage(void)
{
   fprintf(stderr,"\ntopscan V2.0 (c) 1998-2000, Dr. Andrew C.R. Martin, \
UCL & Reading\n");

   fprintf(stderr,"\nUsage: topscan [-t] [-v] [-1] [-n] [-a] [-l] [-p[s]] \
[-w] [-h hlen] [-e elen] [-m matrix] [-g] file1.{dssp|pdb} \
file2.{dssp|pdb}\n");
   fprintf(stderr,"       topscan -b [-1] [-n] [-a] [-l] [-p[s]] [-h hlen] \
[-e elen] [-g] file1.{dssp|pdb}\n");
   fprintf(stderr,"       topscan -s [-t] [-1] [-n] [-a] [-l] [-v] [-p[s]] \
[-w] [-h hlen] [-e elen] [-m matrix] [-g] file1.{dssp|pdb} file2.top\n");
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
   fprintf(stderr,"       -g Treat 3_10 helix as alpha helix\n");
   fprintf(stderr,"       -1 Only do primary topology (ignore \
direction)\n");
   fprintf(stderr,"       -n Add neighbour information\n");
   fprintf(stderr,"       -a Add accessibility information\n");
   fprintf(stderr,"       -l Add element length information\n");
   fprintf(stderr,"       -t Command line has a topology string instead \
of a filename\n");

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
/*>int *ReadTopology(FILE *fp, int ELen, int HLen, BOOL UseStride,
                     BOOL Do3_10, BOOL PrimaryTopology,
                     BOOL DoNeighbour, BOOL DoAccess, BOOL DoLength)
   -----------------------------------------------------------------
   Input:   FILE   *fp             DSSP file pointer
            int    ELen            Minimum length of strand
            int    HLen            Minimum length of helix
            BOOL   UseStride       Read from Stride rather than DSSP file
            BOOL   Do3_10          Merge 3_10 helix with alpha helix
            BOOL   PrimaryTopology Only do primary topology
            BOOL   DoNeighbour     Add neighbour information
            BOOL   DoAccess        Add accessibility information
            BOOL   DoLength        Add length information
   Returns: int *                  Topology string

   Reads the topology from a DSSP or Stride file, returning a string 
   representing the topology.

   13.03.98 Original   By: ACRM
   26.10.99 Added Do3_10
   10.11.99 Added PrimaryTopology
   16.11.99 Added DoNeighbour
   23.11.99 Addes DoAccess
   26.01.00 Added DoLength
   10.03.00 Changed to use integer coded topology array
*/
int *ReadTopology(FILE *fp, int ELen, int HLen, BOOL UseStride,
                  BOOL Do3_10, BOOL PrimaryTopology, BOOL DoNeighbour,
                  BOOL DoAccess, BOOL DoLength)
{
   /* If we are doing neighbours, then reset the internal neighbour
      information
   */
   if(DoNeighbour)
   {
      CalcElement('\0', 0,0,0, 0,0,0, 0, 1, 0, 0.0, 0);
   }
   
   if(UseStride)
   {
      return(ReadStride(fp, ELen, HLen, Do3_10, PrimaryTopology,
                        DoNeighbour, DoAccess, DoLength));
   }
   else
   {
      return(ReadDSSP(fp, ELen, HLen, Do3_10, PrimaryTopology,
                      DoNeighbour, DoAccess, DoLength));
   }
}


/************************************************************************/
/*>int *ReadDSSP(FILE *fp, int ELen, int HLen, BOOL Do3_10,
                 BOOL PrimaryTopology, BOOL DoNeighbour, BOOL DoAccess,
                 BOOL DoLength)
   --------------------------------------------------------------------
   Input:   FILE     *fp             DSSP file pointer
            int      ELen            Minimum length of strand
            int      HLen            Minimum length of helix
            BOOL     Do3_10          Merge 3_10 helix with alpha
            BOOL     PrimaryTopology Only do primary topology
            BOOL     DoNeighbour     Add neighbour information
            BOOL     DoAccess        Add accessibility information
            BOOL     DoLength        Add length information
   Returns: int *                    Topology string

   Reads a DSSP file returning a string representing the topology.

   13.01.98 Original   By: ACRM
   26.10.99 Added Do3_10 handling
   10.11.99 Added PrimaryTopology
   16.11.99 Added DoNeighbour
   23.11.99 Added DoAccess
            Fixed bug - file ending with an SS element wasn't checking
            element length
   26.01.00 Added DoLength
   10.03.00 Changed to use integer coded topology array
   13.03.00 Initialise x1,y1,z1,xp,yp,zp only to silence warnings with -O2
*/
int *ReadDSSP(FILE *fp, int ELen, int HLen, BOOL Do3_10, 
              BOOL PrimaryTopology, BOOL DoNeighbour, BOOL DoAccess,
              BOOL DoLength)
{
   int  *top;
   char buffer[MAXBUFF*2],
        struc,
        LastStruc    = ' ';
   REAL x,  y,  z,
        x1 = MARKER, 
        y1 = MARKER, 
        z1 = MARKER,
        xp = MARKER, 
        yp = MARKER, 
        zp = MARKER,
        access,
        sumaccess    = 0.0;
   BOOL InBody       = FALSE,
        InElement    = FALSE;
   int  i            = 0,
        EleLength    = 0;

#ifdef UCL
   fprintf(stderr,"Code needs to be modified to support reading \
accessibility from\nUCL DSSP files\n");
#endif

   
   if((top = (int *)malloc(MAXBUFF * sizeof(int)))==NULL)
      return(NULL);
   
   while(fgets(buffer,MAXBUFF*2,fp))
   {
      TERMINATE(buffer);
      if(InBody)
      {
#ifdef UCL
         fsscanf(buffer,"%16x%c%90x%7lf%7lf%7lf",&struc,&x,&y,&z);
         access=0.0;
#else
         fsscanf(buffer,"%16x%c%17x%4lf%77x%7lf%7lf%7lf",
                 &struc,&access,&x,&y,&z);
#endif
         if(Do3_10 && struc=='G')                  /* 26.10.99          */
            struc = 'H';
         if((struc=='E' || struc=='H') && (struc==LastStruc))
         {
            EleLength++;
            sumaccess += access;
         }
         
         if((struc=='E' || struc=='H') && (struc!=LastStruc))
         {
            /* Start of new element                                     */
            if(InElement)
            {
               if((LastStruc=='E' && EleLength>=ELen) ||
                  (LastStruc=='H' && EleLength>=HLen))
                  top[i++] = CalcElement(LastStruc,x1,y1,z1,xp,yp,zp,
                                         PrimaryTopology, DoNeighbour,
                                         DoAccess, sumaccess/EleLength,
                                         (DoLength?EleLength:0));
            }
            
            InElement = TRUE;
            EleLength = 1;
            sumaccess = access;
            x1 = x;
            y1 = y;
            z1 = z;
         }
         else if(InElement && struc!='E' && struc!='H')
         {
            /* Just come out of an element                              */
            if((LastStruc=='E' && EleLength>=ELen) ||
               (LastStruc=='H' && EleLength>=HLen))
               top[i++] = CalcElement(LastStruc,x1,y1,z1,xp,yp,zp,
                                      PrimaryTopology, DoNeighbour,
                                      DoAccess, sumaccess/EleLength,
                                      (DoLength?EleLength:0));
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
      if((LastStruc=='E' && EleLength>=ELen) ||
         (LastStruc=='H' && EleLength>=HLen))
         top[i++] = CalcElement(LastStruc,x1,y1,z1,xp,yp,zp,
                                PrimaryTopology, DoNeighbour,
                                DoAccess, sumaccess/EleLength,
                                (DoLength?EleLength:0));
   }
   top[i] = (-1);
   
   return(top);
}


/************************************************************************/
/*>int CalcElement(char struc, REAL x1, REAL y1, REAL z1, 
                   REAL x2, REAL y2, REAL z2, BOOL PrimaryTopology,
                   BOOL DoNeighbour, BOOL DoAccess, REAL meanAccess,
                   int EleLength)
   -----------------------------------------------------------------
   Input:   char  struc           DSSP structure assignment
            REAL  x1              Coordinates of SS element start
            REAL  y1      
            REAL  z1      
            REAL  x2              Coordinates of SS element end
            REAL  y2      
            REAL  z2      
            BOOL  PrimaryTopology Do Primary topology only (no direction)
            BOOL  DoNeighbour     Add neighbour information
            BOOL  DoAccess        Add accessibility information
            REAL  meanAccess      Mean access for element
            int   EleLength       Element length (0 if we are ignoring
                                  these)
   Returns: char                  Code for element (structure & direction)

   Given a structure code from DSSP (E or H) and the coordinates of the
   ends of the element, returns a combined code of structure plus
   direction. Also generates combined codes with neighbour, accessibility
   and element length.

   13.01.98 Original   By: ACRM
   10.11.99 Added PrimaryTopology
   16.11.99 Added DoNeighbour
   23.11.99 Added DoAccess/meanAccess
   26.01.00 Added EleLength
   10.03.00 Changed to use integer coded topology array. i.e. returns
            a number rather than a character
*/
int CalcElement(char struc, REAL x1, REAL y1, REAL z1, 
                REAL x2, REAL y2, REAL z2, BOOL PrimaryTopology,
                BOOL DoNeighbour, BOOL DoAccess, REAL meanAccess,
                int EleLength)
{
   int         code;
   int         dirn,
               buried    = 0,
               lengthmod = 0;
   REAL        dx,    dy,    dz,
               absdx, absdy, absdz;
   static REAL prevx1 = MARKER,
               prevy1 = MARKER,
               prevz1 = MARKER,
               prevx2 = MARKER,
               prevy2 = MARKER,
               prevz2 = MARKER;
   
   /* Special condition to reset the statics                            */
   if(DoNeighbour && struc == '\0')
   {
      prevx1 = MARKER;
      prevy1 = MARKER;
      prevz1 = MARKER;
      prevx2 = MARKER;
      prevy2 = MARKER;
      prevz2 = MARKER;
      return(0);
   }

   /* Base coding is A for sheet, G for helix                           */
   if(PrimaryTopology)   /* We ignore the direction information         */
   {
      if(struc == 'E')
         return(1);
      else if(struc == 'H')
         return(7);
      else
         return(0);
   }

   if(struc == 'E')
      code = 1;
   else if(struc == 'H')
      code = 7;
   else
      return(0);

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
      return(0);
   
   if(DoNeighbour)
   {
      if((prevx1 != MARKER) && (prevy1 != MARKER) && (prevz1 != MARKER) &&
         (prevx2 != MARKER) && (prevy2 != MARKER) && (prevz2 != MARKER))
      {
         if(IsNeighbour(x1, y1, z1,
                        x2, y2, z2,
                        prevx1, prevy1, prevz1,
                        prevx2, prevy2, prevz2))
            dirn += 12;
      }
      prevx1 = x1;
      prevy1 = y1;
      prevz1 = z1;

      prevx2 = x2;
      prevy2 = y2;
      prevz2 = z2;
   }

   if(DoAccess)
   {
      if(((struc == 'E') && (meanAccess < gStrandMeanAccess)) ||
         ((struc == 'H') && (meanAccess < gHelixMeanAccess)))
         buried = 24;
   }

   if(EleLength)
   {
      if(((struc == 'E') && (EleLength > HELIX_MEAN_LENGTH)) ||
         ((struc == 'H') && (EleLength > STRAND_MEAN_LENGTH)))
         lengthmod = 48;
   }
   
   return(code + dirn + buried + lengthmod);
}

/************************************************************************/
/*>int CalcIDScore(int *seq1, int *seq2, BOOL UseBoth)
   ---------------------------------------------------
   Input:   int      *seq1   Sequence 1
            int      *seq2   Sequence 2
            BOOL     UseBoth Calculate score as Max of both sequences
   Returns: int              Max score of each sequence vs itself

   Calculates the maximum possible score resulting from the identical
   sequence

   14.01.98 Original   By: ACRM
   10.03.00 Changed to use integer coded topology array
*/
int CalcIDScore(int *seq1, int *seq2, BOOL UseBoth)
{
   int score1 = 0,
       score2 = 0,
       i,
       seqlen1 = FindArrayLength(seq1),
       seqlen2 = FindArrayLength(seq2);

   for(i=0; i<seqlen1; i++)
   {
      if(seq1[i])
         score1 += NumericCalcMDMScore(seq1[i], seq1[i]);
   }
   if(UseBoth)
   {
      for(i=0; i<seqlen2; i++)
      {
         if(seq2[i])
            score2 += NumericCalcMDMScore(seq2[i], seq2[i]);
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
/*>int *ReadStride(FILE *fp, int ELen, int HLen, BOOL Do3_10, 
                   BOOL PrimaryTopology, BOOL DoNeighbour, BOOL DoAccess,
                   BOOL DoLength)
   -----------------------------------------------------------------------
   Input:   FILE     *fp             Stride file pointer
            int      ELen            Minimum length of strand
            int      HLen            Minimum length of helix
            BOOL     Do3_10          Merge 3_10 helix with alpha
            BOOL     PrimaryTopology Only do primary topology
            BOOL     DoNeighbour     Add neighbour information
            BOOL     DoAccess        Add accessibility information
            BOOL     DoLength        Add length information
   Returns: int *                    Topology string

   Reads a Stride file returning a string representing the topology.

   13.03.98 Original   By: ACRM
   26.10.99 Added Do3_10 handling
   10.11.99 Added PrimaryTopology
   16.11.99 Added DoNeighbour
   23.11.99 Added DoAccess
            Fixed bug - file ending with an SS element wasn't checking
            element length
   26.01.00 Added DoLength
   10.03.00 Changed to use integer coded topology array
   13.03.00 Initialise x1,y1,z1,xp,yp,zp only to silence warnings with -O2
*/
int *ReadStride(FILE *fp, int ELen, int HLen, BOOL Do3_10, 
                BOOL PrimaryTopology, BOOL DoNeighbour, BOOL DoAccess,
                BOOL DoLength)
{
   int  *top;
   char buffer[MAXBUFF*2],
        struc,
        LastStruc    = ' ';
   REAL x,  y,  z,
        x1 = MARKER, 
        y1 = MARKER, 
        z1 = MARKER,
        xp = MARKER, 
        yp = MARKER, 
        zp = MARKER,
        access, 
        sumaccess    = 0.0;
   BOOL InElement    = FALSE;
   int  i            = 0,
        EleLength    = 0;

   
   if((top = (int *)malloc(MAXBUFF * sizeof(int)))==NULL)
      return(NULL);
   
   while(fgets(buffer,MAXBUFF*2,fp))
   {
      TERMINATE(buffer);
      fsscanf(buffer,"%15x%8lf%1x%8lf%1x%8lf%1x%c%1x%8lf",
              &x,&y,&z,&struc,&access);
      
      if(Do3_10 && struc=='G')                     /* 26.10.99          */
         struc = 'H';
      
      if((struc=='E' || struc=='H') && (struc==LastStruc))
      {
         EleLength++;
         sumaccess += access;
      }
      
      if((struc=='E' || struc=='H') && (struc!=LastStruc))
      {
         /* Start of new element                                        */
         if(InElement)
         {
            if((LastStruc=='E' && EleLength>=ELen) ||
               (LastStruc=='H' && EleLength>=HLen))
               top[i++] = CalcElement(LastStruc,x1,y1,z1,xp,yp,zp,
                                      PrimaryTopology, DoNeighbour,
                                      DoAccess, sumaccess/EleLength,
                                      (DoLength?EleLength:0));
         }
         
         InElement = TRUE;
         EleLength = 1;
         sumaccess = access;
         x1 = x;
         y1 = y;
         z1 = z;
      }
      else if(InElement && struc!='E' && struc!='H')
      {
         /* Just come out of an element                                 */
         if((LastStruc=='E' && EleLength>=ELen) ||
            (LastStruc=='H' && EleLength>=HLen))
            top[i++] = CalcElement(LastStruc,x1,y1,z1,xp,yp,zp,
                                   PrimaryTopology, DoNeighbour,
                                   DoAccess, sumaccess/EleLength,
                                   (DoLength?EleLength:0));
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
      if((LastStruc=='E' && EleLength>=ELen) ||
         (LastStruc=='H' && EleLength>=HLen))
         top[i++] = CalcElement(LastStruc,x1,y1,z1,xp,yp,zp,
                                PrimaryTopology, DoNeighbour,
                                DoAccess, sumaccess/EleLength,
                                (DoLength?EleLength:0));
   }
   top[i] = (-1);
   
   return(top);
}


/************************************************************************/
/*>BOOL IsNeighbour(REAL x1, REAL y1, REAL z1,
                    REAL x2, REAL y2, REAL z2,
                    REAL prevx1, REAL prevy1, REAL prevz1,
                    REAL prevx2, REAL prevy2, REAL prevz2)
   -------------------------------------------------------
   Input:   REAL    x1       Nter of this SS element
            REAL    y1 
            REAL    z1
            REAL    x2       Cter of this SS element
            REAL    y2
            REAL    z2
            REAL    prevx1   Nter of previous SS element
            REAL    prevy1
            REAL    prevz1
            REAL    prevx2   Cter of previous SS element
            REAL    prevy2
            REAL    prevz2
   Returns: BOOL             Are they neighbours?

   Determines whether the current secondary structure element is a direct
   neighbour of the preceeding element.

   17.11.99 Original   By: ACRM
   13.03.00 Fixed parameter variable types to REAL (were BOOL!)
*/
BOOL IsNeighbour(REAL x1, REAL y1, REAL z1,
                 REAL x2, REAL y2, REAL z2,
                 REAL prevx1, REAL prevy1, REAL prevz1,
                 REAL prevx2, REAL prevy2, REAL prevz2)
{
   REAL d[4], f[4],
        mindist = 9999.0;
   int  i;
   
   
   /* Calculate distance and position of Nter end against previous
      element
   */
   d[0] = PointLineDistance(x1, y1, z1, 
                            prevx1, prevy1, prevz1,
                            prevx2, prevy2, prevz2,
                            NULL, NULL, NULL,
                            &f[0]);

   /* Calculate distance and position of Cter end against previous
      element
   */
   d[1] = PointLineDistance(x2, y2, z2, 
                            prevx1, prevy1, prevz1,
                            prevx2, prevy2, prevz2,
                            NULL, NULL, NULL,
                            &f[1]);

   /* Calculate distance and position of Nter end of previous element
      against current
   */
   d[2] = PointLineDistance(prevx1, prevy1, prevz1,
                            x1, y1, z1, 
                            x2, y2, z2,
                            NULL, NULL, NULL,
                            &f[2]);

   /* Calculate distance and position of Cter end of previous element
      against current
   */
   d[3] = PointLineDistance(prevx2, prevy2, prevz2,
                            x1, y1, z1, 
                            x2, y2, z2,
                            NULL, NULL, NULL,
                            &f[3]);

   /* See which are in line with the other element                      */
   for(i=0; i<4; i++)
   {
      if((f[i] >= 0.0) && (f[i] <= 1.0))
      {
         if(d[i] < mindist)
            mindist = d[i];
      }
   }
#ifdef DEBUG
   fprintf(stderr,"Minimum distance = %f\n", mindist);
#endif

   /* If the minimum distance is less than our cutoff, then we return
      TRUE as they are adjacent
   */
   if(mindist < ADJACENT_DIST)
      return(TRUE);
   
   return(FALSE);
}


/************************************************************************/
/*>int FindArrayLength(int *array)
   -------------------------------
   Input:   int   *array     Array of integers terminated with -1
   Returns: int              Number of elements

   Takes an integer array and counts how many elements there are up to but
   excluding the first negative number

   08.03.00 Original   By: ACRM
*/
int FindArrayLength(int *array)
{
   int len = 0;
   
   while(*array++ >= 0) len++;
   
   return(len);
}

/************************************************************************/
/*>int MakeIntArray(int *array1, char *inarray)
   --------------------------------------------
   Input:   int    *array1    Integer array to be filled
            char   *inarray   Dash-delimited list of positive integers
   Returns: int               Number of integers copied into array

   Builds an integer array terminated with a (-1) from a dash separated 
   list of numbers
   (e.g. 1-5-7-23-7-31 would go into a 7 element array containing
   1,5,7,23,7,31,-1)

   08.03.00 Original   By: ACRM
*/
int MakeIntArray(int *array1, char *inarray)
{
   int pos = 0;
   char tempbuff[16],
        *chp,
        *buffp;
   
   chp = inarray;
   
   while(*chp)
   {
      buffp = tempbuff;
      while(*chp != '-' && *chp)
      {
         *buffp++ = *chp++;
      }
      *buffp = '\0';
      sscanf(tempbuff,"%d", &(array1[pos++]));
      if(*chp) chp++;
   }
   array1[pos] = (-1);
   
   return(pos);
}


/************************************************************************/
/*>char *NumArrayToString(int *numarr)
   -----------------------------------
   Input:   int    *numarray   Integer array terminated with -1
   Returns: char   *           Pointer to allocated character string
                               containing dash deliminated numbers

   Creates a character representation of an integer array. Each element
   is printed as three characters zero-padded and separated by a dash.
   The first negative number in the array represents the end of the
   array.

   08.03.00 Original   By: ACRM
*/
char *NumArrayToString(int *numarr)
{
   char *buff = NULL,
        tmpbuff[16],
        *tbp;
   int  buffsize = BUFFCHUNK,
        i,
        buffused = 0,
        strsize;
   
   if((buff = (char *)malloc(BUFFCHUNK * sizeof(char)))==NULL)
      return(NULL);
   buff[0] = '\0';

   for(i=0; numarr[i] >= 0; i++)
   {
      sprintf(tmpbuff,"%03d",numarr[i]);
      KILLLEADSPACES(tbp,tmpbuff);

      /* Check whether the number string will fit, if not increase the
         array size
      */
      strsize = strlen(tbp)+1;
      if((buffused + strsize + 1) > buffsize)
      {
         buffsize += BUFFCHUNK;
         buff = (char *)realloc(buff, buffsize);
      }

      /* Add a dash if this isn't the first number added                */
      if(i) strcat(buff, "-");
      /* Now add the number itself                                      */
      strcat(buff,tbp);
      buffused = strlen(buff);
   }
   
   return(buff);
}

