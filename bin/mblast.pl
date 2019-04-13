#!/usr/bin/perl -w
use strict;
use warnings;
use Carp;
use FindBin qw($RealBin);
use File::Basename;
use Getopt::Long;

my %opts;
my $VERSION="2.0";
GetOptions( \%opts,"i=s", "p=s","o=s",'d=s','e=s','v=i',"h!");

my $usage = <<"USAGE";
       Program : $0
       Version : $VERSION
       Contact : meng.luo\@majorbio.com
       Lastest modify:2019-04-13
       Discription: run mblast
       Usage :perl $0 [options]
                -i*	fasta		input fasta file     
				-p*	blastx/blastp blasttype
				-o*	output file	
				-d*	database		find under system Variable \$MBLASTX_DATADIR or use absolute path
				-e	evalue			default 1.0E-05
				-v 	10				Maximum number of MSPs to print.  Default value is 10
                -h					Display this usage information
                * 					must be given Argument
                exmaple:perl $0 -p blastx -i test.fa -o out.xml -d nt
USAGE


die $usage if ((!$opts{i} || !$opts{p} || !$opts{d}|| ! $opts{o})||$opts{h});

$opts{e}=$opts{e}?$opts{e}:'1.0E-05';
$opts{v}=$opts{v}?$opts{v}:'10';

my @startmem=(0,20,64,96,140,196,256);
my $mem=140;#shift(@startmem);
my $setp=24;
my $up=0;
my $down=0;

die("Input file $opts{i} not exists!\n") unless(-f $opts{i});
my $blasttype;
if($opts{p}=~/blastx/i  ){
	$blasttype="mblastx";
}elsif($opts{p}=~/blastp/i){
	$blasttype="mblastp";
}else{
	die("Blast type mast be blastp or blastx!\n") 
}

&runblast();
warn("Convert to XML fomart ... \n");
system("$RealBin/mblast2blastxml.pl -i $opts{o}.mblast -db ".basename($opts{d})." -p $blasttype -o $opts{o} ");
system("rm $opts{o}.mblast");


sub runblast{
	my $mblastpath="$RealBin/mblastx";
	warn("try with $mem G  memory paramas !\n");
	my $cmd;
	#="rm -rf ./key.lic;ln -s $mblastpath/key.lic .;";	
	my $db=$opts{d};

	if($opts{d}=~/^\//){
		$cmd.="source /etc/profile;export MBLASTX_DATADIR=".dirname($opts{d}).";";
		$db=basename($opts{d});
	}else{
		#die("Please set system variable MBLASTX_DATADIR first!\n") unless($ENV{'MBLASTX_DATADIR'});
		unless( $ENV{'MBLASTX_DATADIR'}  && $ENV{'MBLASTX_DATADIR'} =~ /^\//){
			$cmd.="source /etc/profile;export MBLASTX_DATADIR=/mnt/ilustre/app/rna/database/mblast;";
		}
	}
	
	system("cp -rf $mblastpath/key.lic .");
	if($mem==0){
		$cmd.="$mblastpath/$blasttype -c $db -q $opts{i} -f F -T 16 -e $opts{e} -m $opts{v} -o $opts{o}.mblast  2>&1";
	}else{
		$cmd.="$mblastpath/$blasttype -c $db -q $opts{i} -f F -T 16 -M $mem -e $opts{e} -m $opts{v} -o $opts{o}.mblast  2>&1";
	}
	
	#print($cmd);
	warn("Run mblast ... \n");
	print "$cmd\n";
	my $out=`$cmd`;
	unless($out=~/completed successfully/i){
		if($out=~/ERROR|Trail version expired|Query length is bigger/){
			warn($out);
		}elsif($out=~/failed/){			
			warn($out);			
			if($mem>1000 || $mem<0||$setp<2){
				warn("Run mblast failed\n!");
			}else{
				if($out=~/hHits_UberQueryOffset|malloc\(\) failed|uProt allocUProt/){
					$setp=$setp/2   if($down);
					$up=1;
					$down=0;
				 	$mem+=$setp;
				}
				if($out=~/hMspLog|hitsInHitBuffer/){
					$setp=$setp/2   if($up);
					$up=0;
					$down=1;
					$mem-=$setp;
				}


				&runblast();
				if($setp<2){
					if(scalar(@startmem)>0){
						$mem=140;#shift(@startmem);
						$setp=24;
						$up=0;
						$down=0;
					}					
				}
			}
		}else{
			warn("Mblast not correctly completed, rerun!\n");
			&runblast();
		}
		
	}
	system("rm -rf ./key.lic");
}
