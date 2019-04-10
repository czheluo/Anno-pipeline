#!/usr/bin/perl -w
use strict;
use warnings;
use Carp;

if(@ARGV < 1 ) {
    print STDERR "mergeGO.pl <1.go.list> <2.go.list> .. <n.go.list>\n";
    print STDERR "mergeGO.pl <*.go.list>\n";
    exit;
}

my %gos;

foreach my $f (@ARGV){
	#print "$f\n";
	open(GO, "< $f") || die "Can't open $f\n";
	while(<GO>){
		chomp;
		next if(/^\s*#/);
		/^(\S+)\s*(.*)$/;
		my $a=$1;
		my $b=$2;
		if(defined $b){
			#print $b."\n";
			while($b=~/(GO:\d+)/gi){
				#print "$1\n";
				$gos{$a}{$1}=1;
			}			
		}
	}
	close GO;
}

foreach my $x (keys %gos){
		print "$x\t".join(";",keys(%{$gos{$x}}))."\n";
}