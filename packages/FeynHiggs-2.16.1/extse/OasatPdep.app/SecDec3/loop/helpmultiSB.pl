#! /usr/bin/perl -w
#
my $graph=$ARGV[0];
my $pointname=$ARGV[1];
my $ms1=$ARGV[2]**2;
my $ms2=$ARGV[3]**2;
my $ms3=$ARGV[4]**2;
my $ms4=$ARGV[5]**2;
my $mhig1s=$ARGV[6];
my $mhig2s=$ARGV[7];
my $mas=$ARGV[8]**2;


$dirbase=`pwd`;
chomp $dirbase;

#Specify special points:
#
@specp =($mhig1s,$mhig2s,$mas);
$specp = @specp;
#print "$specp and @specp\n";
#
#
print "\ngraph= $graph \n";
print "pointname= $pointname\n";
print "masses=$ms1 $ms2 $ms3 $ms4\n";
print "m_H^2=$mhig1s m_h^2=$mhig2s m_A^2=$mas\n";

#$0 = 0;
#$ms1="$$ms1";
#$ms2="$$ms2";
#$ms3="$$ms3";
#$ms4="$$ms4";
#
$destination=$dirbase;
$paramfile="${graph}.input";
$filename="multi${pointname}.input";
$lines=3; #+2 because of 1 special additional point
#
open(EWRITE, ">","$destination/$filename") || die "cannot open $destination/$filename\n";
print EWRITE "paramfile=$paramfile\n";
print EWRITE "pointname=$pointname\n";
print EWRITE "# only lines 1 to [lines] will be calculated\n";
print EWRITE "lines=$lines\n";
print EWRITE "#---------------------------------------------\n";
print EWRITE "# each line below contains numsij sij's, then numpi2 pi^2, then numms2 mi^2\n";
print EWRITE "numsij=0\n";
print EWRITE "numpi2=1\n";
print EWRITE "numms2=4\n";
print EWRITE "xplot=1\n";
print EWRITE "#---------------------------------------------\n";
foreach $pnt (@specp) {
    print EWRITE "$pnt,$ms1,$ms2,$ms3,$ms4\n";
}
close(EWRITE);
print "$filename written.\n";







 

