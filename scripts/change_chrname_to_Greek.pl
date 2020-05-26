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

    if($chr eq "1"){ $chr = "I"; }
    elsif($chr eq "2"){ $chr = "II"; }
    elsif($chr eq "3"){ $chr = "III"; }
    elsif($chr eq "4"){ $chr = "IV"; }
    elsif($chr eq "5"){ $chr = "V"; }
    elsif($chr eq "6"){ $chr = "VI"; }
    elsif($chr eq "7"){ $chr = "VII"; }
    elsif($chr eq "8"){ $chr = "VIII"; }
    elsif($chr eq "9"){ $chr = "IX"; }
    elsif($chr eq "10"){ $chr = "X"; }
    elsif($chr eq "11"){ $chr = "XI"; }
    elsif($chr eq "12"){ $chr = "XII"; }
    elsif($chr eq "13"){ $chr = "XIII"; }
    elsif($chr eq "14"){ $chr = "XIV"; }
    elsif($chr eq "15"){ $chr = "XV"; }
    elsif($chr eq "16"){ $chr = "XVI"; }

    print "$head$chr";
    for($i=1;$i<=$#clm;$i++) {
      print "\t$clm[$i]";
    }
  }
close IN;
