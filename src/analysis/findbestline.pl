#!/usr/bin/perl

$ignore = shift;
$term = shift;

($ignore,$junk) = split(/[\._]/,$ignore);

$line = 0;
$theline = 0;
$thescore = -1;

while(<>)
{
    chomp;
    $line++;
    if(/$term/)
    {
        if(!/$ignore/)
        {
            $theline = $line;
            ($domid,$cat,$score) = split;
            last;
        }
    }
}

print "$theline $score\n";
