/*************************************************************************

   Program:    mergepdbsecstr
   File:       mergepdbsecstr.c
   
   Version:    V1.0
   Date:       06.08.18
   Function:   Merge original PDB file with PDBSECSTR secondary structure
               assignments
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2018
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

   pdbsecstr <pdbfile> | mergepdbsecstr <pdbfile> > output

**************************************************************************

   Revision History:
   =================
   V1.0  06.08.18 Original

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
#define MAXBUFF     160
#define MAXNAMEBUFF 8

typedef struct _pdbsecstr
{
   struct _pdbsecstr *next,
                     *prev;
   int               resnum;
   char              resnam[MAXNAMEBUFF],
                     insert[MAXNAMEBUFF],
                     chain[MAXNAMEBUFF],
                     ss;
} SECSTR;

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL DoMerge(FILE *pdbfp, FILE *pdbsecstrfp, FILE *out);
BOOL ParseCmdLine(int argc, char **argv, char *pdbfile, char *infile, 
                  char *outfile);
SECSTR *ReadPdbsecstr(FILE *fp);
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
            fprintf(stderr,"mergepdbsecstr: Error opening PDB file\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"mergepdbsecstr: Error opening PDBSECSTR \
input or output file\n");
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
/*>BOOL DoMerge(FILE *pdbfp, FILE *pdbsecstrfp, FILE *out)
   -------------------------------------------------------
   Input:   FILE    *pdbfp       PDB file pointer
            FILE    *pdbsecstrfp PDBSECSTR output file pointer
            FILE    *out         File pointer for output
   Returns: BOOL                 Success?

   Merges PDBSECSTR output data with coordinate information from a 
   PDB file

   13.03.98 Orginal   By: ACRM
*/
BOOL DoMerge(FILE *pdbfp, FILE *pdbsecstrfp, FILE *out)
{
   PDB    *pdb, 
          *p;
   SECSTR *pdbsecstr, 
          *s;
   int    natoms;

   /* Read the PDB file                                                 */
   if((pdb = blReadPDB(pdbfp, &natoms))==NULL)
   {
      fprintf(stderr,"mergepdbsecstr: No atoms read from PDB file\n");
      return(FALSE);
   }

   /* Reduce to CA atoms                                                */
   pdb = blSelectCaPDB(pdb);

   /* Read the PDBSECSTR file                                           */
   if((pdbsecstr = ReadPdbsecstr(pdbsecstrfp))==NULL)
   {
      fprintf(stderr,"mergepdbsecstr: No records read from PDBSECSTR \
file\n");
      return(FALSE);
   }

   /* Replace a chain name of ' ' in the PDB file with '-' as used by
      pdbsecstr
   */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->chain[0] == ' ')
         p->chain[0] = '-';
   }

   /* Run through the PDB linked list finding the matching record from
      the pdbsecstr link list and outputting the composite record. Delete
      the record from the Pdbsecstr list to speed things up
   */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      for(s=pdbsecstr; s!=NULL; NEXT(s))
      {
         /* If they match                                               */
         if((p->resnum    == s->resnum)       &&
            INSERTMATCH(p->insert, s->insert) &&
            CHAINMATCH(p->chain, s->chain))
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
                  pdbsecstr = NULL;
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
               pdbsecstr = s->next;
               free(s);
            }

            /* Break out of the search through SECSTR                   */
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
/*>SECSTR *ReadPdbsecstr(FILE *fp)
   -------------------------------
   Input:     FILE     *fp      PDBSECSTR output file to be read
   Returns:   SECSTR   *        Linked list of SECSTR assignment data

   Reads a PDBSECSTR file and creates a doubly linked list of items \
   from the file

   13.03.98 Original   By: ACRM
*/
SECSTR *ReadPdbsecstr(FILE *fp)
{
   char   buffer[MAXBUFF],
          ResnumBuff[8];
   SECSTR *pdbsecstr = NULL,
          *s;
   
   
   while(fgets(buffer,MAXBUFF,fp))
   {
      TERMINATE(buffer);

      if(pdbsecstr==NULL)
      {
         INITPREV(pdbsecstr,SECSTR);
         s=pdbsecstr;
      }
      else
      {
         ALLOCNEXTPREV(s,SECSTR);
      }
      if(s==NULL)
      {
         return(NULL);
      }

      sscanf(buffer,"%s %s %c",
             ResnumBuff,
             s->resnam,
             &(s->ss));

      blParseResSpec(ResnumBuff, s->chain, &(s->resnum), s->insert);
   }

   return(pdbsecstr);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   13.03.98 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nmergepdbsecstr V1.0 (c) 1998 UCL, Dr. Andrew C.R. \
Martin\n");

   fprintf(stderr,"\nUsage: mergepdbsecstr pdbfile [pdbsecstrfile \
[outputfile]]\n");

   fprintf(stderr,"\nMerges coordinate data from the PDB file with \
PDBSECSTR secondary structure\n");
   fprintf(stderr,"assignments.\n");

   fprintf(stderr,"\nIf the pdbsecstrfile and outputfile are not specified, \
stdin and stdout\n");
   fprintf(stderr,"are used.\n\n");
}
