/*************************************************************************

   Program:    mergestride
   File:       mergestride.c
   
   Version:    V2.0
   Date:       06.08.18
   Function:   Merge original PDB file with STRIDE secondary structure
               assignments
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 1998-2018
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      +44 (0)207 679 7034
   EMail:      andrew@bioinf.org.uk
               
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
   Generally used as:

   stride <pdbfile> | mergestride <pdbfile> > output

**************************************************************************

   Revision History:
   =================
   V1.0  13.03.98 Original
   V2.0  06.08.18 Updated for new Bioplib

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160

typedef struct _stride
{
   struct _stride *next,
                  *prev;
   int            resnum;
   char           resnam[8],
                  insert,
                  chain,
                  ss;
} STRIDE;

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL DoMerge(FILE *pdbfp, FILE *stridefp, FILE *out);
BOOL ParseCmdLine(int argc, char **argv, char *pdbfile, char *infile, 
                  char *outfile);
STRIDE *ReadStride(FILE *fp);
void Usage(void);

/************************************************************************/
int main(int argc, char **argv)
{
   FILE *in    = stdin,
        *out   = stdout,
        *pdbfp = NULL;
   char infile[MAXBUFF],
        pdbfile[MAXBUFF],
        outfile[MAXBUFF];

   if(ParseCmdLine(argc, argv, pdbfile, infile, outfile))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdbfp = fopen(pdbfile,"r"))!=NULL)
         {
            DoMerge(pdbfp, in, out);
         }
         else
         {
            fprintf(stderr,"mergestride: Error opening PDB file\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"mergestride: Error opening stride input or \
output file\n");
         return(1);
      }
   }
   else
   {
      Usage();
   }

   return(0);
}


/************************************************************************/
/*>BOOL DoMerge(FILE *pdbfp, FILE *stridefp, FILE *out)
   ----------------------------------------------------
   Input:   FILE    *pdbfp       PDB file pointer
            FILE    *stridefp    STRIDE output file pointer
            FILE    *out         File pointer for output
   Returns: BOOL                 Success?

   Merges STRIDE output data with coordinate information from a PDB file

   13.03.98 Orginal   By: ACRM
*/
BOOL DoMerge(FILE *pdbfp, FILE *stridefp, FILE *out)
{
   PDB    *pdb, 
          *p;
   STRIDE *stride, 
          *s;
   int    natoms;

   /* Read the PDB file                                                 */
   if((pdb = blReadPDB(pdbfp, &natoms))==NULL)
   {
      fprintf(stderr,"mergestride: No atoms read from PDB file\n");
      return(FALSE);
   }

   /* Reduce to CA atoms                                                */
   pdb = blSelectCaPDB(pdb);

   /* Read the STRIDE file                                              */
   if((stride = ReadStride(stridefp))==NULL)
   {
      fprintf(stderr,"mergestride: No records read from STRIDE file\n");
      return(FALSE);
   }

   /* Replace a chain name of ' ' in the PDB file with '-' as used by
      stride
   */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->chain[0] == ' ')
         p->chain[0] = '-';
   }
   

   /* Run through the PDB linked list finding the matching record from
      the stride link list and outputting the composite record. Delete
      the record from the Stride list to speed things up
   */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      for(s=stride; s!=NULL; NEXT(s))
      {
         /* If they match                                               */
         if((p->resnum    == s->resnum) &&
            (p->insert[0] == s->insert) &&
            (p->chain[0]  == s->chain))
         {
            fprintf(out, "%4s %c %5d %c %8.3f %8.3f %8.3f %c\n",
                    p->resnam,
                    p->chain[0],
                    p->resnum,
                    p->insert[0],
                    p->x,
                    p->y,
                    p->z,
                    s->ss);

            /* Now unlink this from the list                            */
            if((s->next != NULL) && (s->prev != NULL))
            {
               /* Middle of list                                        */
               s->prev->next = s->next;
               s->next->prev = s->prev;
               free(s);
            }               
            else if(s->next == NULL)
            {
               /* End of list                                           */
               if(s->prev == NULL)
               {
                  stride = NULL;
               }
               else
               {
                  s->prev->next = NULL;
               }
               free(s);
            }
            else
            {
               /* Beginning of list                                     */
               if(s->next != NULL)
                  s->next->prev = NULL;
               stride = s->next;
               free(s);
            }

            /* Break out of the search through stride                   */
            break;
         }
      }
   }
   return(TRUE);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *pdbfile, char *infile, 
                     char *outfile)
   ----------------------------------------------------------------------
   Input:   int    argc        Argument count
            char   **argv      Argument array
   Output:  char   *pdbfile    Input PDB file
            char   *infile     Input filename (or blank string)
            char   *outfile    Output filename (or blank string)
   Returns: BOOL               Success

   Parse the command line

   16.08.94 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *pdbfile, char *infile, 
                  char *outfile)
{
   argc--;
   argv++;
   
   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            return(FALSE);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 1, 2 or 3 arguments left               */
         if((argc < 1) || (argc > 3))
            return(FALSE);
         
         /* Copy the first to pdbfile                                   */
         strcpy(pdbfile, argv[0]);
         argc--;
         argv++;
         
         if(argc)
         {
            /* Copy the next to infile                                  */
            strcpy(infile, argv[0]);
            
            /* If there's another, copy it to outfile                   */
            argc--;
            argv++;
            if(argc)
               strcpy(outfile, argv[0]);
         }
         
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>STRIDE *ReadStride(FILE *fp)
   ----------------------------
   Input:     FILE     *fp         STRIDE output file to be read
   Returns:   STRIDE   *           Linked list of STRIDE assignment data

   Reads a STRIDE file and creates a doubly linked list of items from the
   file

   13.03.98 Original   By: ACRM
*/
STRIDE *ReadStride(FILE *fp)
{
   char   buffer[MAXBUFF],
          junk[8],
          ResnumBuff[8];
   STRIDE *stride = NULL,
          *s;
   int    iResnum,
          lastchar;
   
   
   while(fgets(buffer,MAXBUFF,fp))
   {
      TERMINATE(buffer);

      if(!strncmp(buffer,"ASG ",4))
      {
         if(stride==NULL)
         {
            INITPREV(stride,STRIDE);
            s=stride;
         }
         else
         {
            ALLOCNEXTPREV(s,STRIDE);
         }
         if(s==NULL)
         {
            return(NULL);
         }

         sscanf(buffer,"%s %s %c %s %d %c",
                junk,
                s->resnam,
                &(s->chain),
                ResnumBuff,
                &iResnum,
                &(s->ss));

         lastchar = strlen(ResnumBuff) - 1;
         if(isalpha(ResnumBuff[lastchar]))
         {
            s->insert = ResnumBuff[lastchar];
            ResnumBuff[lastchar] = '\0';
         }
         else
         {
            s->insert = ' ';
         }
         s->resnum = atoi(ResnumBuff);
      }
   }

   return(stride);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   13.03.98 Original   By: ACRM
   06.08.18 V2.0
*/
void Usage(void)
{
   fprintf(stderr,"\nmergestride V2.0 (c) 1998 UCL, Dr. Andrew C.R. \
Martin\n");

   fprintf(stderr,"\nUsage: mergestride pdbfile [stridefile \
[outputfile]]\n");

   fprintf(stderr,"\nMerges coordinate data from the PDB file with \
STRIDE secondary structure\n");
   fprintf(stderr,"assignments.\n");

   fprintf(stderr,"\nIf the stridefile and outputfile are not specified, \
stdin and stdout\n");
   fprintf(stderr,"are used.\n\n");
}
