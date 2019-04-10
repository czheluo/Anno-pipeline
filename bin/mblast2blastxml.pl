#! /usr/bin/perl -w

use strict;
use warnings;


use Getopt::Long;
my %opts;

GetOptions( \%opts,"i=s", "o=s","db=s","p=s","h!");

my $usage = <<"USAGE";
       Program : $0
       Version : 1.0
       Contact:
       Lastest modify: 20121018
       Discription: translate mblast result into xml format
       Usage :perl $0 [options]
                -i		Gallus_gallus.pep_vs_nr.blast.bsp    	result of mblast
                -o		Gallus_gallus.pep_vs_nr.blast.bsp.xml	output.file
				-db		nr|genes|string|kobas_ko|kobas_d.rerio|kobas_g.gallus|kobas_h.sapiens|kobas_o.sativa|kobas_s.cerevisiae				blast database
				-p		mblastp|mblastx				program type		
                -h      	Display this usage information
                
                
USAGE




my $xml_version="1.0";
my $par_filter="F";
my $blast_ref;




die $usage if ( !( $opts{i} &&  $opts{db} && $opts{p} ) || $opts{h} );

$opts{o}=$opts{o}?$opts{o}:$opts{i}.".xml";
$opts{db}=$opts{db}?$opts{db}:"unknown";



if($opts{p} eq "mblastx"){
	$par_filter='L;';
	$blast_ref= 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402.';
	
}else{
	$par_filter="F";
	$blast_ref='~Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, ~Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), ~&quot;,  Nucleic Acids Res. 25:3389-3402.';
}


my $param = <<"PARAM";
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>1e-05</Parameters_expect>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>$par_filter</Parameters_filter>
    </Parameters> 
  </BlastOutput_param>
PARAM

my $program=$opts{p};
my $program_version=uc($program)." 2.2.24+ [Feb-01-2011]";

open (MBLAST, "< $opts{i}") or die "Error: Couldn't open $opts{i}\n";
open (XML, ">$opts{o}") or die "Error: Couldn't open $opts{o}\n";

print XML '<?xml version="'.$xml_version.'"?>'."\n".'<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">'."\n";
print XML '<BlastOutput>'."\n";
print XML '  <BlastOutput_program>'.$program.'</BlastOutput_program>'."\n";
print XML '  <BlastOutput_version>'.$program_version.'</BlastOutput_version>'."\n";
print XML '  <BlastOutput_reference>'.$blast_ref.'</BlastOutput_reference>'."\n"; 
print XML '  <BlastOutput_db>'.$opts{db}.'</BlastOutput_db>'."\n";


my $label=0;
my $alignlabel=0;
my $hitnamelabel=0;
my $hitid=200000;

my $db_num=1000;
my $db_len=100000;
my $lcl = 0;
my $queryName = "";
my $queryLength = 0;
my $hit_num = 0;
my %hit=();
if ($opts{db} eq "genes"){
	$db_num=5502463;
	$db_len=1973349375;	
}elsif ($opts{db} eq "string"){
	$db_num=5214234;
	$db_len=1903618658;
}elsif ($opts{db} eq "string"){
	$db_num=21062489;
	$db_len=7218481314;
}elsif ($opts{db} eq "kobas_ko"){
	$db_num=6973788;
	$db_len=2496150079;
}elsif ($opts{db} eq "kobas_d.rerio"){
	$db_num=26566;
	$db_len=14075278
}elsif ($opts{db} eq "kobas_g.gallus"){
	$db_num=18135;
	$db_len=8605407;
}elsif ($opts{db} eq "kobas_h.sapiens"){
	$db_num=19761;
	$db_len=10876086;
}elsif ($opts{db} eq "kobas_o.sativa"){
	$db_num=28453;
	$db_len=10265149
}elsif ($opts{db} eq "kobas_s.cerevisiae"){
	$db_num=5882;
	$db_len=2914779;
}

my $stat = <<"STAT";
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>$db_num</Statistics_db-num>
          <Statistics_db-len>$db_len</Statistics_db-len>
          <Statistics_hsp-len>0</Statistics_hsp-len>
          <Statistics_eff-space>0</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>0.267</Statistics_lambda>
          <Statistics_entropy>0.14</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
STAT

while(<MBLAST>){
	chomp;
	if ($_ =~ /^Query\t:\s(.*)$/){
		if ($lcl == 1){
			&PrintHead();
		}
		if ($lcl != 0){
			&PrintQuery();
		}
		
		$queryName = $1;
		$lcl++;
		$label=0;
		$alignlabel = 1;		
		$hit_num=0;
		%hit=();
		next;
	}
	if ($_ =~ /^Length\t:\s([\d]*)/){
		$queryLength = $1;
		next;
	}
	if ($_=~/Sequences producing significant alignments:/){
		$label=1;
		next;
	}
	if ($_=~/^>/){
		$label=2;
		$alignlabel = 0;
		
	}

	if($label == 2){

		if ($_=~/^>\s([\S]*)/){
			$hit_num++;
			$hitnamelabel=1;
				if($_=~/^>\s([\s\S]*)/){
					$hit{$hit_num}{name}=$1;
					$hitid += 1;
					$hit{$hit_num}{acess}=$hitid;
					$hit{$hit_num}{id}="gnl|BL_ORD_ID|".$hitid;					
				}				
				

			$hit{$hit_num}{query}="";
			$hit{$hit_num}{sbjct}="";
			$hit{$hit_num}{mid}="";
			$hit{$hit_num}{query_from}=9999999;
			$hit{$hit_num}{query_to}=0;
			$hit{$hit_num}{sbjct_from}=9999999;
			$hit{$hit_num}{sbjct_to}=0;
			next;
		}		
		
		if ($_=~/Length=([\d]*)/){
			$hit{$hit_num}{len}=$1;
			$hitnamelabel=0;
			next;
			
		}
		if ($hitnamelabel==1){
			$hit{$hit_num}{name} .= $_;
		}
		if ($_=~/Score = ([\S]*) bits \(([\d]*)\), Expect = (.*)/){
			$hit{$hit_num}{bit_Score}=$1;
			$hit{$hit_num}{Score}=$2;
			$hit{$hit_num}{e_value}=$3;	
			next;
		}
		if ($_=~/Identities = (\d*)\/\d* \([\s\d]*%\), Positives = (\d*)\/\d* \([\s\d]*%\), Gaps = (\d*)\/\d*/){
			$1=~ s/\s//g;
			$2=~ s/\s//g;
			$hit{$hit_num}{identity}=$1;
			$hit{$hit_num}{positive}=$2;
			$hit{$hit_num}{Hsp_gaps}=$3;
			next;
		}

		# ��ʼλ�á���ֹλ����Ҫ�Ľ�
		if ($_=~/Query\s*([\d-]*)\s*(\S*)\s*([\d-]*)/){
			$alignlabel = 1;
				if ( $hit{$hit_num}{query_from}>$1 ){
				$hit{$hit_num}{query_from}=$1;
			}
			if ( $hit{$hit_num}{query_to}<$3 ){
				$hit{$hit_num}{query_to}=$3;
			}
			$hit{$hit_num}{query}.=$2;
			next;
			
		}
		if ($alignlabel==1){
			if ($_=~/^                     ([\s\S]*)$/){
				$hit{$hit_num}{mid} .= $1;
			}
		}
		if ($_=~/Sbjct\s*([\d-]*)\s*(\S*)\s*([\d-]*)/){
			$alignlabel = 0;
			if ( $hit{$hit_num}{sbjct_from}>$1 ){
				$hit{$hit_num}{sbjct_from}=$1;
			}
			if ( $hit{$hit_num}{sbjct_to}<$3 ){
				$hit{$hit_num}{sbjct_to}=$3;
			}
			$hit{$hit_num}{sbjct}.=$2;
		}		
	}	
}

&PrintQuery();
print XML '  </BlastOutput_iterations>'."\n".'</BlastOutput>'."\n";

close XML;
close MBLAST;

sub PrintQuery(){
	print XML '    <Iteration>'."\n";
	print XML '      <Iteration_iter-num>'.$lcl.'</Iteration_iter-num>'."\n";
	if ($opts{p} eq "mblastx" || $opts{db} eq "string" ){
		print XML '      <Iteration_query-ID>Query_'.$lcl.'</Iteration_query-ID>'."\n";
	}else{
		print XML '      <Iteration_query-ID>lcl|'.$lcl.'_0</Iteration_query-ID>'."\n";
	}
	
	print XML '      <Iteration_query-def>'.$queryName.'</Iteration_query-def>'."\n";
	print XML '      <Iteration_query-len>'.$queryLength.'</Iteration_query-len>'."\n";
	print XML '      <Iteration_hits>'."\n";
	
	foreach(sort by_number keys %hit){
		print XML '        <Hit>'."\n";
		print XML '          <Hit_num>'.$_.'</Hit_num>'."\n";
		print XML '          <Hit_id>'.$hit{$_}{id}.'</Hit_id>'."\n";
		
		if( $hit{$_}{name} =~ /[\S\s]*\]/ ){
			my @names =split(/\]/,$hit{$_}{name});
			$hit{$_}{name}=$names[0].']';
		}elsif( $hit{$_}{name} =~ /([\S\s]*)\001([\S\s]*)/ ){
			$hit{$_}{name} = $1;
		}
		
		if ($hit{$_}{name} =~/([^\&]*)[\&]([^\&]*)/){
			$hit{$_}{name} = $1." and ".$2;
		}
		
		$hit{$_}{name}=~ s/\>/\)/g;
		$hit{$_}{name}=~ s/\</\(/g;
		#
		print XML '          <Hit_def>'.$hit{$_}{name}.'</Hit_def>'."\n";
		print XML '          <Hit_accession>'.$hit{$_}{acess}.'</Hit_accession>'."\n";
		print XML '          <Hit_len>'.$hit{$_}{len}.'</Hit_len>'."\n";
		print XML '          <Hit_hsps>'."\n".'            <Hsp>'."\n";
		print XML '              <Hsp_num>1</Hsp_num>'."\n";
		print XML '              <Hsp_bit-score>'.$hit{$_}{bit_Score}.'</Hsp_bit-score>'."\n";
		print XML '              <Hsp_score>'.$hit{$_}{Score}.'</Hsp_score>'."\n";
		print XML '              <Hsp_evalue>'.$hit{$_}{e_value}.'</Hsp_evalue>'."\n";
		print XML '              <Hsp_query-from>'.$hit{$_}{query_from}.'</Hsp_query-from>'."\n";
		print XML '              <Hsp_query-to>'.$hit{$_}{query_to}.'</Hsp_query-to>'."\n";
		print XML '              <Hsp_hit-from>'.$hit{$_}{sbjct_from}.'</Hsp_hit-from>'."\n";
		print XML '              <Hsp_hit-to>'.$hit{$_}{sbjct_to}.'</Hsp_hit-to>'."\n";
		print XML '              <Hsp_query-frame>1</Hsp_query-frame>'."\n";
		print XML '              <Hsp_hit-frame>1</Hsp_hit-frame>'."\n";
		print XML '              <Hsp_identity>'.$hit{$_}{identity}.'</Hsp_identity>'."\n";
		print XML '              <Hsp_positive>'.$hit{$_}{positive}.'</Hsp_positive>'."\n";
		if ($opts{p} eq "mblastx"){
			print XML '              <Hsp_gaps>0</Hsp_gaps>'."\n";
		}
		
		my $len = rindex($hit{$_}{mid}."\$", "\$");
		print XML '              <Hsp_align-len>'.$len.'</Hsp_align-len>'."\n";
		print XML '              <Hsp_qseq>'.$hit{$_}{query}.'</Hsp_qseq>'."\n";
		print XML '              <Hsp_hseq>'.$hit{$_}{sbjct}.'</Hsp_hseq>'."\n";
		print XML '              <Hsp_midline>'.$hit{$_}{mid}.'</Hsp_midline>'."\n";
		print XML '           </Hsp>'."\n".'          </Hit_hsps>'."\n".'        </Hit>'."\n";
	}	
	print XML $stat;
}

sub PrintHead(){
	my $query_ID='lcl|1_0';
	if($opts{p} eq "mblastx" || $opts{db} eq "string"){
		$query_ID='Query_1';
	}
	print XML '  <BlastOutput_query-ID>'.$query_ID.'</BlastOutput_query-ID>'."\n";
	print XML '  <BlastOutput_query-def>'.$queryName.'</BlastOutput_query-def>'."\n";
	print XML '  <BlastOutput_query-len>'.$queryLength.'</BlastOutput_query-len>'."\n";
	print XML $param;
	print XML '  <BlastOutput_iterations>'."\n";
	
}

sub by_number { $a<=>$b }




