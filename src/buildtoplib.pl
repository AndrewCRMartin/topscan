#!/usr/local/bin/perl
#*************************************************************************
#
#   Program:    
#   File:       
#   
#   Version:    
#   Date:       
#   Function:   
#   
#   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 1997
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
$minhelix  = 3;
$minstrand = 4;

#*************************************************************************
$domdir   = "/nfs/cathdata/dompdb";
$pdbprep  = "/data/pdb_release/all/pdb";
$pdbext   = ".ent";
$topscan  = "/home/bsm/martin/bin/topscan";
$getchain = "/home/bsm/martin/bin/getchain";

#*************************************************************************
use Pg;

#*************************************************************************
if(!defined $ENV{'dbname'})
{
    die "You must define the dbname environment variable to the database name";
}
$dbname = $ENV{'dbname'};

$cleanup = 0;

$conn   = Pg::connectdb("dbname=$dbname");
$query  = "SELECT domid, pdbcode, chainid, domain FROM domains WHERE nrep = 't'";
$result = $conn->exec($query);
$ntups  = $result->ntuples;

for($i=0; $i<$ntups; $i++)
{
    $domid   = $result->getvalue($i, $result->fnumber("domid"));
    $pdbcode = $result->getvalue($i, $result->fnumber("pdbcode"));
    $chainid = $result->getvalue($i, $result->fnumber("chainid"));
    $domain  = $result->getvalue($i, $result->fnumber("domain"));

    print STDERR "Processing $domid...";

    if($domain eq "0")
    {
        if($chainid eq "0")
        {
            $file = $pdbprep . $pdbcode . $pdbext;
        }
        else
        {
            $file = $pdbprep . $pdbcode . $pdbext;
            `$getchain $chainid $file /tmp/$domid`;
            $cleanup = 1;
            $file = "/tmp/$domid";
        }
    }
    else
    {
        $file = $domdir . "/" . $domid;
    }

    $topstring = `$topscan -b -p -h $minhelix -e $minstrand $file`;
    print $topstring;
    if($cleanup)
    {
        unlink $file;
        $cleanup = 0;
    }

    print STDERR "done\n";
}


