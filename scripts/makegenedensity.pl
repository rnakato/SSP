#!/usr/bin/perl -w

if($#ARGV != 2){
    print " makegenedensity.pl <genometable> <refFlat> <windowsize>\n";
    exit;
}

$gtfile = $ARGV[0];
$genefile = $ARGV[1];
$width = $ARGV[2];

open(InputFile, $gtfile) ||die "error: can't open $gtfile.\n";
while(<InputFile>){
    next if($_ eq "\n");
    chomp;
    my @clm= split(/\t/, $_);
    $name{$clm[0]} = $clm[0];
    $len{$clm[0]} = $clm[1];
    $nwin{$clm[0]}=$len{$clm[0]}/$width;
    for($i=0; $i<$nwin{$clm[0]}; $i++){ $array{$clm[0]}[$i]=0;}
}
close (InputFile);

open(ListFile, $genefile) ||die "error: can't open $genefile.\n";
while(<ListFile>){
    next if($_ eq "\n");
    chomp;
    my @clm = split(/\t/, $_);
    $Hash_chr{$clm[0]} = $clm[2];
    $Hash_start{$clm[0]} = $clm[4];
#    $Hash_end{$clm[0]} = $clm[5];
}
close (ListFile);

foreach $name (keys(%Hash_chr)){
    my $chr = $Hash_chr{$name};
    my $s = $Hash_start{$name};
 #   my $e = $Hash_end{$name};
    $array{$chr}[int($s/$width)]++;
}

foreach $chr (keys(%name)){
    open(OUT, ">$chr-bs$width") ||die "error: can't open $chr-bs$width.\n";
    for($i=0;$i<$nwin{$chr};$i++){
	printf OUT "%d\t%d\n", $i*$width, $array{$chr}[$i];
    }
    close(OUT);
}
