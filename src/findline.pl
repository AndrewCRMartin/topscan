#!/usr/local/bin/perl

$term = shift;
$line = 0;
$theline = 0;
$thescore = -1;

while(<>)
{
    chomp;
    $line++;
    if(/$term/)
    {
        $theline = $line;
        ($domid,$cat,$score) = split;
    }
}

print "$theline $score\n";
