#! /usr/bin/perl -w

open(IN, $ARGV[0]) || die;
while(<IN>) {
    next if($_ eq "\n");

    my @clm = split(/\t/, $_);
    my $head = "";
    my $chr = $clm[0];
    if($chr =~ /chr(.+)/) {
      $head="chr";
      $chr = $1;
    }

    if($chr eq "I"){ $chr = "1"; }
    elsif($chr eq "II"){ $chr = "2"; }
    elsif($chr eq "III"){ $chr = "3"; }
    elsif($chr eq "IV"){ $chr = "4"; }
    elsif($chr eq "V"){ $chr = "5"; }
    elsif($chr eq "VI"){ $chr = "6"; }
    elsif($chr eq "VII"){ $chr = "7"; }
    elsif($chr eq "VIII"){ $chr = "8"; }
    elsif($chr eq "IX"){ $chr = "9"; }
    elsif($chr eq "X"){ $chr = "10"; }
    elsif($chr eq "XI"){ $chr = "11"; }
    elsif($chr eq "XII"){ $chr = "12"; }
    elsif($chr eq "XIII"){ $chr = "13"; }
    elsif($chr eq "XIV"){ $chr = "14"; }
    elsif($chr eq "XV"){ $chr = "15"; }
    elsif($chr eq "XVI"){ $chr = "16"; }

    print "$head$chr";
    for($i=1;$i<=$#clm;$i++) {
      print "\t$clm[$i]";
    }
  }
close IN;
