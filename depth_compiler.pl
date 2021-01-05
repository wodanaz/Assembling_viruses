#!/usr/bin/env perl
#use Bio::Perl;
my %bigh;

$usage="

data_compiler.pl : 

assembles data into a single table

Arguments:
arg1: [pattern to counts files] 

Example:
data_compiler.pl *.tab > alldata.tab

";

if (!$ARGV[0]) { die $usage;} 

print "\t";
foreach $argi (0..$#ARGV)
        {
        print $ARGV[$argi], "\t";
        open(TAB, $ARGV[$argi]);
        while(<TAB>)
                {
                chomp;
                @cols = split("\t", $_);
                $bigh{$cols[0]}{$argi} = $cols[1];
                }
        }
print "\n";

foreach $r (sort(keys(%bigh)))
        {
        print $r, "\t";
        foreach $argi (0..$#ARGV)
                {
                if(exists($bigh{$r}{$argi}))
                        {print $bigh{$r}{$argi}, "\t";
                        }
                else    {print 0, "\t";}
                }
        print "\n";
        }

