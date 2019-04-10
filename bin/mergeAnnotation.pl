#!/usr/bin/perl -w
use strict;
use warnings;
use Carp;
use DBI;
use Getopt::Long;
use File::Basename;
my %opts;
my $VERSION="1.0";
GetOptions( \%opts,"fa=s", "nt=s","nr=s","string=s","cog=s","swissport=s","uniprot=s","uniprottype=s","go=s","goslim=s","kobas=s","interpro=s","target=s","diff=s","other=s","h!");

my $usage = <<"USAGE";
       Program : $0
       Version : $VERSION
       Contact : quan.guo\@majorbio.com
       Lastest modify:2012-04-13
       Discription:merge all anntations into a big table
       Usage :perl $0 [options]
			-fa*	fastafile	fasta file for run blast
			-nt	fasta_vs_nt.table.xls	blast table export by Blast2table -format 10
			-nr	fasta_vs_nr.table.xls	blast table export by Blast2table -format 10
			-swissport fasta_vs_swissport.table.xls blast table export by Blast2table -format 10,blast with  uniprot_sprot
			-uniprot	fasta_vs_uniprot.table.xls	blast table export by Blast2table -format 10
			-uniprottype	UniProtKB	type of uniprot database which you run blast with,support UniProtKB,UniRef100,UniRef90,UniRef50,Default:UniProtKB
			-interpro	test.gene.interpro.txt	Interproscan result export by Interpro_extract.pl
			-go	go.list	GO list file 
			-goslim	goslim.list	goslim list file
			-string	fasta_vs_string.table.xls	table create by /data/users/guoquan/workspace/perl/String2Cog.pl
			-cog	cog.list	COG list file
			-kobas	kobas.out		output file of kobas2
			-target	target.txt	miRNA target prediction file created by target_check.pl
			-diff	A_vs_B.txt,B_vs_C.txt	diff genelist,file name must like sample1_vs_sample2.txt,split by ",",can use wildcard character in bash,like '*.txt',but must use '' when using  wildcard character
			-other	annotation.xls	other annotation table ,first column must be fastaname, title must be in the first row
			-h			display help message
			
      
USAGE

die $usage if(!$opts{fa} || $opts{h});
$opts{uniprottype}=$opts{uniprottype}?$opts{uniprottype}:"UniProtKB";
my %fasta;
my %genes;
my %orf;
my ($fastaname,$fastaseq);
open(FA, "< $opts{fa}") || die "Can't open $opts{fa}\n";
while(<FA>){
	chomp;
	if(/>(\S*)/){
		$fasta{$fastaname}=length($fastaseq) if($fastaname);
		$fastaname=$1;
		$fastaseq="";
		if(/gene=(\S+)/){
			$genes{$fastaname}=$1;	
		}
		if(/orf=(\S+)/){
			$orf{$fastaname}=$1;	
		}
	}else{
		$fastaseq.=$_;
	}
}
$fasta{$fastaname}=length($fastaseq);		
close FA;

my (%output,$title);
$title="#Qeury_name\tQeury_length";
foreach my $f (keys(%fasta)){
	$output{$f}="$f\t$fasta{$f}";
}
my $isgene=0;
$isgene=1 if scalar(keys %genes)>0;
if($isgene){
	$title.="\tgene_name";
	foreach my $f (keys(%fasta)){
		$output{$f}.="\t$genes{$f}";
	}
}

if(scalar(keys %orf)>0){
	$title.="\tORF";
	foreach my $f (keys(%orf)){
		$output{$f}.="\t$orf{$f}";
	}
}

if($opts{target}){	
my %miRNA;
	open(TA, "< $opts{target}") || die "Can't open $opts{target}\n";
	while(<TA>){
		chomp;
		my @a=split /\t/;
		my @b=split(";",$a[2]);
		foreach my $x (@b){
			$miRNA{$x}{$a[0]}=1;
		}
	}
	close TA;
	$title.="\tmiRNA_number\tmiRNA" if scalar(keys %miRNA)>0;
	foreach my $f (keys(%fasta)){
		if(exists($miRNA{$f})){
			$output{$f}.="\t".scalar(keys(%{$miRNA{$f}}))."\t".join(";",keys(%{$miRNA{$f}}));
		}else{
			$output{$f}.="\t_\t_";
		}
		
	}
}

if($opts{diff}){
	my %difflist;
	my @diff=split(",",$opts{diff});
	foreach my $d (@diff){
		my @file= glob $d;
		foreach my $f (@file){
			my $fname=basename($f);
			my $s;
			if($fname=~/(\w+)_vs_(\w+)/i){
				$s=$1.":".$2;
			}else{
				$fname=~/^([^\.]+)/;
				$s=$1;
			}
			open(DIFF, "< $f") || die "Can't open $f}\n";
			while(<DIFF>){
				chomp;
				/^(\S+)/;
				my $a=$1;
				$difflist{$1}{$s}=1;
			}
			close DIFF;
		}
	}
		$title.="\tdiff_exp" if scalar(keys %difflist)>0;
	foreach my $f (keys(%fasta)){
		if(exists($difflist{$f})){
			$output{$f}.="\t".join(";",keys(%{$difflist{$f}}));
		}else{
			$output{$f}.="\t_"
		}
	}
}

if($opts{nt}){
	my %nt;
	open(NT, "< $opts{nt}") || die "Can't open $opts{nt}\n";
	while(<NT>){
		chomp;
		my @a=split /\t/;
		next if exists($nt{$a[5]});
		$nt{$a[5]}="\t$a[11]\t$a[16]\t$a[4]";
	}
	close NT;
	$title.="\tNT_tophit_name\tNT_tophit_description\tNT_topHSP_%-Simil" if scalar(keys %nt)>0;
	foreach my $f (keys(%fasta)){
		if(exists($nt{$f})){
			$output{$f}.=$nt{$f};
		}else{
			$output{$f}.="\t_\t_\t_";
		}
		
	}
}

if($opts{nr}){
	my %nr;
	open(NR, "< $opts{nr}") || die "Can't open $opts{nr}\n";
	while(<NR>){
		chomp;
		my @a=split /\t/;
		next if exists($nr{$a[5]});
		$nr{$a[5]}="\t$a[11]\t$a[16]\t$a[4]";
	}
	close NR;
	$title.="\tNR_tophit_name\tNR_tophit_description\tNR_topHSP_%-Simil" if scalar(keys %nr)>0;
	foreach my $f (keys(%fasta)){
		if(exists($nr{$f})){
			$output{$f}.=$nr{$f};
		}else{
			$output{$f}.="\t_\t_\t_";
		}
		
	}
}

if($opts{swissport}){
	my %swissport;
	open(SW, "< $opts{swissport}") || die "Can't open $opts{swissport}\n";
	while(<SW>){
		chomp;
		my @a=split /\t/;
		next if exists($swissport{$a[5]});
		$swissport{$a[5]}="\t$a[11]\t$a[16]\t$a[4]";
	}
	close SW;
	$title.="\tSwissprot_tophit_name\tSwissprot_tophit_description\tSwissprot_tophsp_%-Simil" if scalar(keys %swissport)>0;
	foreach my $f (keys(%fasta)){
		if(exists($swissport{$f})){
			$output{$f}.=$swissport{$f};
		}else{
			$output{$f}.="\t_\t_\t_";
		}
	}
}

if($opts{uniprot}){
	my %uniprot;
	open(UN, "< $opts{uniprot}") || die "Can't open $opts{uniprot}\n";
	while(<UN>){
		chomp;
		my @a=split /\t/;
		next if exists($uniprot{$a[5]});
		$uniprot{$a[5]}="\t$a[11]\t$a[16]\t$a[4]";
	}
	close UN;
	$opts{uniprottype}=$opts{uniprottype}=~/UniProtKB/?"uniprot_trembl":$opts{uniprottype};
	$title.="\tUniprot_tophit_name($opts{uniprottype})\tUniprot_tophit_description\tUniprot_topHSP_%-Simil" if scalar(keys %uniprot)>0;
	foreach my $f (keys(%fasta)){
		if(exists($uniprot{$f})){
			$output{$f}.=$uniprot{$f};
		}else{
			$output{$f}.="\t_\t_\t_";
		}
		
	}
}
if($opts{interpro}){
	my %interpro;
	open(IPR, "< $opts{interpro}") || die "Can't open $opts{interpro}\n";
	while(<IPR>){
		chomp;
		my @a=split /\t/;
		while(/\%(IPR\d+)/g ){
   			$interpro{$a[0]}{$1}=1;
   		}
	}
	close IPR;
	$title.="\tInterpro_IDs" if scalar(keys %interpro)>0;
	foreach my $f (keys(%fasta)){
		if(exists($interpro{$f})){
			$output{$f}.="\t".join(";",keys(%{$interpro{$f}}));
		}else{
			$output{$f}.="\t_";
		}		
	}
}
if($opts{go}){
	my %go;
	open(GO, "< $opts{go}") || die "Can't open $opts{go}\n";
	while(<GO>){
		chomp;
		my @a=split /\t/;
		next if(scalar(@a)<2);
		my @b=split(";",$a[1]);
		$go{$a[0]}=\@b;
	}
	close GO;
	$title.="\tGOs" if scalar(keys %go)>0;
	foreach my $f (keys(%fasta)){
		if(exists($go{$f})){
			$output{$f}.="\t".join(";",@{$go{$f}});
		}else{
			$output{$f}.="\t_";
		}
		
	}
}

if($opts{goslim}){
	my %goslim;
	open(GOSLIM, "< $opts{goslim}") || die "Can't open $opts{goslim}\n";
	while(<GOSLIM>){
		chomp;
		my @a=split /\t/;
		next if(scalar(@a)<2);
		my @b=split(";",$a[1]);
		$goslim{$a[0]}=\@b;
	}
	close GOSLIM;
	$title.="\tGOslim" if scalar(keys %goslim)>0;
	foreach my $f (keys(%fasta)){
		if(exists($goslim{$f})){
			$output{$f}.="\t".join(";",@{$goslim{$f}});
		}else{
			$output{$f}.="\t_";
		}
		
	}
}

if($opts{string}){
	my %string;
	open(STR, "< $opts{string}") || die "Can't open $opts{string}\n";
	while(<STR>){
		chomp;
		my @a=split /\t/;
		next if exists($string{$a[0]});
		$string{$a[0]}="\t$a[5]\t$a[6]\t$a[17]";
	}
	close STR;
	$title.="\tString_tophit_name\tString_tophit_description\tString_topHSP_%-Simil" if scalar(keys %string)>0;
	foreach my $f (keys(%fasta)){
		if(exists($string{$f})){
			$output{$f}.=$string{$f};
		}else{
			$output{$f}.="\t_\t_\t_";
		}
		
	}
}

if($opts{cog}){
	my %cog;
	open(COG, "< $opts{cog}") || die "Can't open $opts{cog}\n";
	my $x=<COG>;
	chomp($x);
	my @c=split(/\t/,$x);
	while(<COG>){
		chomp;		
		my @a=split /\t/;
		$cog{$a[0]}{$c[1]}=$a[1]?$a[1]:'_';
		$cog{$a[0]}{$c[2]}=$a[2]?$a[2]:'_';
		$cog{$a[0]}{$c[3]}=$a[3]?$a[3]:'_';
	}
	close COG;
	$title.="\t$c[1]\t$c[2]\t$c[3]" if scalar(keys %cog)>0;
	foreach my $f (keys(%fasta)){
		if(exists($cog{$f})){
			$output{$f}.="\t$cog{$f}{$c[1]}\t$cog{$f}{$c[2]}\t$cog{$f}{$c[3]}";
		}else{
			$output{$f}.="\t_\t_\t_";
		}
		
	}
}

if($opts{kobas}){
	my %kobas;
	open(KO, "< $opts{kobas}") || die "Can't open $opts{kobas}\n";
	while(<KO>){
		chomp;
		next if(/^\s*#/);
		last if(/^\/\/\/\//);
		if(/^(\S*)\s+(.*)$/){
			my $a=$1;
			my $b=$2;
			#print "$a\t$b\n";
			my @c=split(/\|/,$b);
			$c[1]="_" unless(defined $c[1]);
			$kobas{$a}="\t$c[0]\t$c[1]";
		}
	}
	close KO;
	$title.="\tKO/Gene_ID\tKEGG_GENE_NAME" if scalar(keys %kobas)>0;
	foreach my $f (keys(%fasta)){
		if(exists($kobas{$f})){
			$output{$f}.=$kobas{$f};
		}else{
			$output{$f}.="\t_\t_";
		}
		
	}
}
if($opts{other}){
	my %other;
	open(OT, "< $opts{other}") || die "Can't open $opts{other}\n";
	my $x=<COG>;
	chomp($x);
	my @c=split(/\t/,$x);
	shift(@c);
	while(<OT>){
		chomp;
		my @b=split /\t/;
		my $n=shift(@b);
		$other{$n}=join("\t",@b);
	}
	close OT;
	$title.="\t".join("\t",@c) if scalar(keys %other)>0;
	foreach my $f (keys(%fasta)){
		$output{$f}.="\t".$other{$f};
	}
}
print "$title\n";
foreach my $f (keys(%fasta)){
	print $output{$f}."\n";
}
