#!/usr/bin/perl

##########################################################################
#
#     BDIST v. 1.0
#
#     Author:  Todd Wood
#
#     Purpose:  An implementation of Robinson and Cavanaugh's
#       baraminic distance method.  Please refer to CRSQ 34(4):196-208
#       for more information on the method.
#
#     Usage: bdist.pl [datamatrix]
#
#     Restrictions:
#       1.  "datamatrix" is the file that contains the Cladistic dataset
#           that you want to examine.
#       2.  Matrix must be in proper format.  First column should be taxon
#           name, followed by characters in separate columns.
#       3.  Output is always written to bdist.txt by default.
#       4.  Please see the BSG website (http://www.bryancore.org/bsg for
#           more information.
#
##########################################################################




$infile=$ARGV[0];
$outfile=">bdist.txt";
$version="1.0";
$date=scalar localtime;

open OUT, $outfile;
print OUT "BDIST performs Baraminic Distance analyses on cladistic datasets.\n";
print OUT "  version $version October 3, 2001\n";
print OUT "\n  Please cite:\n\n";
print OUT "  Robinson, D.A. and D.P. Cavanaugh.  1998.  A quantitative approach to\n";
print OUT "\tbaraminology with examples from the primates.  CRSQ 34:196-208.\n";
print OUT "  Wood, T.C.  2001.  BDIST software, v. $version.  Center for Origins Research\n";
print OUT "\tand Education, Bryan College.  Distributed by the author.\n\n\n";

open IN, $infile;
$ntaxa=0;
while (<IN>) {
    $ntaxa++;
    chomp;@m=split;$name=shift(@m);
    @{ $dat{$name} }=@m; $nc{$name}=@m;$nchar=@m;
}
close IN;

@tname=sort keys %dat;

foreach (sort keys %nc) {
    if($nc{$_} != $nchar) {
	print "Taxon $_ has wrong number of characters.  Please check your input file.\n";
	exit(1);
    }
}

print OUT " Input file: $infile\n Date: $date\n Number of taxa:  $ntaxa\n Number of characters:  $nchar\n\n";

$nchar--;
$novern1=$ntaxa/($ntaxa-1);

print OUT "Character\tRelevance\n";
print OUT "---------\t---------\n";

$badchar=0;
for$j(0..$nchar) {
    foreach $t(@tname) {
	if($dat{$t}[$j] ne "?") {
	    $a[$j]++;
	    $csfreq[$j]{$dat{$t}[$j]}++;
	}
    }
    $a[$j]/=$ntaxa;if($a[$j]<0.95) {$badchar++;}
    $x=$j+1;printf OUT "%5d\t\t%5.3f\n", $x, $a[$j];
}
print OUT "\nCharacters with a<0.95:  $badchar\n";
$goodchar=$nchar-$badchar+1;
print OUT "$goodchar characters will be used in the BDIST analysis\n\n";

print OUT "Character\tDiversity\n";
print OUT "---------\t---------\n";

for$j(0..$nchar) {
    if($a[$j]<0.95) {next;}
    @char=sort keys %{ $csfreq[$j] };
    foreach$k(@char) {
	$x=$csfreq[$j]{$k}/$ntaxa;
	$x=($x**2)*$novern1;
	$c[$j]+=$x;
    }
    $c[$j]=1-$c[$j];
    $x=$j+1;printf OUT "%5d\t\t%5.3f\n", $x, $c[$j];
}

print OUT "\n\n--------------------Distance Matrix--------------------\n\n";

$ntaxa--;
$totchar=$goodchar;
for$j(0..$ntaxa) {
    $next=$j+1;$d[$j][$j]=0;
    for$k($j..$ntaxa) {
	$questions=0;
	for$x(0..$nchar) {
	    if($a[$x]<0.95) {next;}
	    if($dat{$tname[$j]}[$x]=~/\?/ || $dat{$tname[$k]}[$x]=~/\?/) {
		$questions++;
		next;
	    }
	    if($dat{$tname[$j]}[$x]=~/\{/ && $dat{$tname[$k]}[$x]!~/\{/) {
		if($dat{$tname[$j]}[$x]!~/$dat{$tname[$k]}[$x]/) {
		    $mismatch[$j][$k]++;$mismatch[$k][$j]++;
		}
	    } elsif($dat{$tname[$k]}[$x]=~/\{/ && $dat{$tname[$j]}[$x]!~/\{/) {
		if($dat{$tname[$k]}[$x]!~/$dat{$tname[$j]}[$x]/) {
		    $mismatch[$j][$k]++;$mismatch[$k][$j]++;
		}
	    } elsif($dat{$tname[$k]}[$x]=~/\{/ && $dat{$tname[$j]}[$x]=~/\{/) {
		$match=0;
		$charlen=length($dat{$tname[$k]}[$x])-1;
		for$q(0..$charlen) {
		    $thisstate=substr($dat{$tname[$k]}[$x],$q,1);
		    if($thisstate ne "{" && $thisstate ne "}" && $dat{$tname[$j]}[$x]=~/$thisstate/) {
			$match=1;
		    }
		}
		if($match==0) {$mismatch[$j][$k]++;$mismatch[$k][$j]++;}
	    } elsif($dat{$tname[$j]}[$x] ne $dat{$tname[$k]}[$x]) {
		$mismatch[$j][$k]++;$mismatch[$k][$j]++;
	    }
	}
	$d[$j][$k]=$mismatch[$j][$k]/($totchar-$questions);
	$d[$k][$j]=$mismatch[$k][$j]/($totchar-$questions);
    }
}

for$j(0..$ntaxa) {
    print OUT "\t$tname[$j]";
}
print OUT "\n";
for$j(0..$ntaxa) {
    print OUT "$tname[$j]";
    for$k(0..$ntaxa) {
	printf OUT "\t%5.3f", $d[$j][$k];
    }
    print OUT "\n";
}
