#!/usr/local/bin/perl
#*************************************************************************
#
#   Program:    analyse
#   File:       analyse.perl
#   
#   Version:    V1.0
#   Date:       23.02.98
#   Function:   
#   
#   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 1998
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   Phone:      (Home) +44 (0)1372 275775
#               (Work) +44 (0)171 419 3890
#   EMail:      martin@biochem.ucl.ac.uk
#               andrew@stagleys.demon.co.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#   Analyses the output from topscan by reprinting as the domain ID,
#   CAT code and score.
#   Normally one would feed the output of topscan through something like
#      | sort +1 -n | tail -50 | analyse.perl
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#
#*************************************************************************
use Pg;

#*************************************************************************
if(!defined $ENV{'dbname'})
{
    die "You must define the dbname environment variable to the database name";
}
$dbname = $ENV{'dbname'};

#*************************************************************************
$conn   = Pg::connectdb("dbname=$dbname");

while(<>)
{
    chomp;
    ($name, $score) = split;
    if($name =~ /pdb.*\.ent/)
    {
        $name =~ s/.*\/pdb//;
        $name =~ s/\.ent//;
        $name .= "00";
    }
    else
    {
        $name =~ s/.*\///g;
    }
    $query  = "SELECT cat FROM domains WHERE domid = '$name'";
    $result = $conn->exec($query);
    $ntups  = $result->ntuples;
    if($ntups != 1)
    {
        $cat = "???";
    }
    else
    {
        $cat   = $result->getvalue(0, $result->fnumber("cat"));
    }

    print "$name $cat $score\n";
}



