#!/usr/bin/perl -w
use strict;
use warnings;
use Carp;
use Bio::SearchIO;
use Getopt::Long;
use FindBin qw($RealBin);
use DBI qw(:sql_types);
use Try::Tiny;
use LWP::Simple;
use LWP::UserAgent;
use HTML::TreeBuilder;
use Math::Round qw(:all);
use lib "$RealBin/PerlLib/";
use PBS::Host;
use GD;
my %opts;
my $VERSION="3.0";
GetOptions( \%opts,"i=s", "format=s","o=s","maxEvalue=f","remote=s","minIdentity=i","org=s","fresh!","rank=i","database=s","QminCoverage=i","HminCoverage=i","parse_id!","use_proxy!","proxy_server=s","h!");

my $usage = <<"USAGE";
       Program : $0
       Version : $VERSION
       Contact : meng.luo\@majorbio.com
       Lastest modify:2019-04-13
       Discription:parse blast to genes databse result and get kegg pathway info and map
       Usage :perl $0 [options]
                -i*		blastn.out		blast to genes database output,can use wildcard character in bash,like '*_blast.out',but must use '' when using  wildcard character     
                -format		blastformat		the format of blast output
												 kobas	kobas2 anntation file
                								 blast      BLAST (WUBLAST, NCBIBLAST,bl2seq)   defualt
 												 fasta      FASTA -m9 and -m0
  												 blasttable BLAST -m9 or -m8 output (both NCBI and WUBLAST tabular)
  												 megablast  MEGABLAST
  												 psl        UCSC PSL format
												 waba       WABA output
												 axt        AXT format
  												 sim4       Sim4
  												 hmmer      HMMER2 hmmpfam and hmmsearch or HMMER3 hmmscan and hmmsearch
												 exonerate  Exonerate CIGAR and VULGAR format
												 blastxml   NCBI BLAST XML
												 wise       Genewise -genesf format
				-maxEvalue	1e-5			Max E value,default:1e-6
				-minIdentity	75			Min percent(Positives) identity over the alignment,default:75
				-QminCoverage	0			Min Converage percent of Query ,Suggest set to >70 when Query seq is protien
				-HminCoverage	30			Min Converage percent of Hit,Defualt:30
				-rank	10		rank cutoff for valid hits from BLAST result, Defualt:10
				-parse_id			parse_id from query description not the Query ID
                -o		dir			output dir,defualt kegg_out under current dir                 
                -org		organism		organism name of three letters ,list in http://www.genome.jp/kegg/catalog/org_list.html ,like hsa
                								defuat:ko
                								also can use:map
                -fresh					fresh database from network
                -database	database path		defaut:/state/partition1/kegg/kegg.db or /mnt/ilustre/app/rna/database/kegg/kegg.db
                -use_proxy	whether use http proxy
                -proxy_server		http proxy server address default:http://100.168.10.80:8888/
				-remote				yes|no   yes:get picture from kegg web  no:get picture from local
                -h					Display this usage information
                * 					must be given Argument
                exmaple:perl $0 -i 'unfinish_*.out' -format blastxml -minIdentity 70
USAGE

######################
#default paramas
die $usage if ((!$opts{i})||$opts{h});
$opts{format}=$opts{format}?$opts{format}:"blast";
$opts{o}=$opts{o}?$opts{o}:"./kegg_out";
$opts{org}=$opts{org}?$opts{org}:"ko";
# if(-f "/state/partition1/kegg/kegg.db"){
	# $opts{database}="/state/partition1/kegg/kegg.db";
# }else{
	# $opts{database}="/mnt/ilustre/app/rna/database/kegg/kegg.db";
# }
unless($opts{database}){
	 $opts{database}="/mnt/ilustre/app/rna/database/kegg/kegg.db";
	 $opts{database}="/state/partition1/kegg/kegg.db" unless(-e $opts{database});
}

$opts{proxy_server}=$opts{proxy_server}?$opts{proxy_server}:"http://101.168.10.80:8888/";
$opts{maxEvalue}=$opts{maxEvalue}?$opts{maxEvalue}:"1e-6";
$opts{minIdentity}=$opts{minIdentity}?$opts{minIdentity}:75;
$opts{HminCoverage}=$opts{HminCoverage}?$opts{HminCoverage}:30;
$opts{rank}=$opts{rank}?$opts{rank}:"10";

my $remote = "no";
if($opts{org} ne "ko" && $opts{org} ne "map"){
	$remote = "yes";
}
my $local_img_dir = "$RealBin/kegg";

unless(-f $opts{database}){
	warn("Database not exists,Create new ...\n");
	$opts{fresh}=1;
}

my $ua =LWP::UserAgent->new();
$ua->timeout(20);
if($opts{use_proxy}){
	my $host=PBS::Host->new({hostname=>'101.168.10.80'});
	if($host->checkstat() eq 'alive'){
		$ua->proxy('http',$opts{proxy_server});
	}else{
		$ua->proxy('http',"http://101.168.10.81:8888/");
	}	
}


my $dbh = DBI->connect("dbi:SQLite:dbname=$opts{database}","","",{AutoCommit => 1});


#####################################################
# check databases

my $check=$dbh->prepare("select count(*) from sqlite_master where type='table' and name='pathway_".$opts{org}."'");

$check->execute();
my @row_ary  = $check->fetchrow_array;
if ($row_ary[0]<=0){
	$opts{fresh}=1;
	warn("Local database has no info of this organism,getting from kegg network ...\n");
}

if($opts{org} =~ /^ko$|^map$/i){
	$check=$dbh->prepare("select count(*) from sqlite_master where type='table' and name='ko_pathway_".$opts{org}."'");
}else{
	$check=$dbh->prepare("select count(*) from sqlite_master where type='table' and name='gene_pathway_".$opts{org}."'");
	$check->execute();
	my @row_ary  = $check->fetchrow_array;
	if ($row_ary[0]<=0){
		$opts{fresh}=1;
		warn("Local database has no info of this organism,getting from kegg network ...\n");
	}
}

##################################################
#fresh database

&freshdatabase($opts{org}) if($opts{fresh});

###################################################
#read blast input
my @file= glob $opts{i};
warn("Input blast result files:\n");
warn(join("\n",@file)."\n");

my $pathway = &getpathways($opts{org});
my $kos =&getpathwaykos($opts{org});

my %seqkos;
mkdir("$opts{o}","493") or die "Can't create dir at $opts{o}\n" unless( -e $opts{o});

#############################################################
# read kobas result
open(KEGG, "> $opts{o}/kegg_table.xls") || die "Can't open $opts{o}/kegg_table.xls\n";
if($opts{format} eq 'kobas'){
	print KEGG "Query\tKo id(Gene id)\tKo name(Gene name)\thyperlink\tPaths\n";
	foreach my $f (@file){
		open ANNOT, "< $f" or die "Error:Cannot open file  $f : $! \n";
		my $line=<ANNOT>;
		my $kobas_specie;
		if($line=~/##(\w+)\s+KEGG Orthology|##Species:\s+(\w{2,4})\s+/){
				$kobas_specie=$1;
		}
		if($opts{org} =~ /^ko$|^map$/i){
			die "The  species in kobas file  $f is $kobas_specie but your input is  $opts{org} \n! "  unless($kobas_specie =~ /^ko$|^map$/i)
		}else{
			die "The  species in kobas file  $f is $kobas_specie but your input is  $opts{org} \n! "  if($kobas_specie ne $opts{org} );
		}
		
		while(<ANNOT>){
			chomp;
			next if(/^\s*#/);
			last if(/^\/\/\/\//);
			if(/^(\S*)\s+(.*)$/){
				my $query_name=$1;
				my $annot=$2;
				next if $annot=~/^None/;	
				#print $annot."\n";			
				my @a=split(/\|/,$annot);				
				my $ko=$a[0];
				#print $ko."\n";
				my @paths;
				my @kos;
				push(@kos,$ko);	
				$seqkos{$query_name}=\@kos;
				my $koid;
				if($opts{org} =~ /^ko$|^map$/i){
					$koid="ko:".$ko;
				}else{
					$koid=$ko;
				}				
				if(exists($kos->{$koid})){
					foreach my $p (keys(%{$kos->{$koid}})){
						#print "$p\n";	
						push(@paths,$p);
						push(@{$pathway->{$p}{'kos'}},$koid);
						push(@{$pathway->{$p}{'seqs'}},$query_name);					
					}
				}else{
					warn("warn:$ko is not in all pathways in this organism $opts{org},if your database is newest, this is ok ... \n");
				}
				my $paths_ref=&uniq(\@paths);
				print KEGG "$query_name\t$a[0]\t$a[1]\t$a[2]\t".join(";",@$paths_ref)."\n";
			}
		}
		close ANNOT;
	}
}else{
###########################################################
# pharse blast result
	print KEGG "Queryname\tHitname\tDBLink\tHit_discription\tevalue\tScore\ttopHSP_strand\tMax_identity\tQuery_length\ttopHSP_Query_converage\tHit_length\tmaxHSP_Hit_coverage\tkos/genes\tecs\tpathway\n";
	foreach my $f (@file){
		warn("Parsing blast result file $f ...\n");
		my $searchio= Bio::SearchIO->new(-format => $opts{format},
									 -file => $f,
									 -best => 1,
									);
		while(my $result = $searchio->next_result){
				my $algorithm=$result->algorithm();
				die "Only support blastp and blastx result!\n" unless($algorithm=~/blastx|blastp/i);
				my $query_name=$result->query_name;
				if($opts{parse_id}||$query_name=~/^Query_\d+$/){
					$query_name=$result->query_description;
					$query_name=~/^\s*(\S+)/;
					$query_name=$1;				
				}else{
					$query_name=$result->query_name;
					$query_name=~/^\s*(\S+)/;
					$query_name=$1;
				}
				my $query_length=$result->query_length;
	
				my @quiery_ko;
			while(my $hit = $result->next_hit){
				last if $hit->rank() > $opts{rank};
				my $hit_length=$hit->length();
				my $score=$hit->score();
				my @paths;
				my @kos;
				my @ecs;
				my $hsp= $hit->hsp; #Bio::Search::HSP::HSPI
				my ($query_hsp_length,$hit_hsp_length);
				
					$query_hsp_length=$hsp->length('qeury');
					$hit_hsp_length=$hsp->length('hit');
					#print "$query_name\t$b\n";
	
				my ($query_coverage,$hit_coverage);
				$query_coverage=$query_hsp_length/$query_length;
				$hit_coverage=$hit_hsp_length/$hit_length;
				if($opts{'QminCoverage'}){
					next if $query_coverage <$opts{'QminCoverage'}/100;
				}
				if($opts{'HminCoverage'}){
					next if $hit_coverage <$opts{'HminCoverage'}/100;
				}
				if($opts{'maxEvalue'}){
					last if $hsp->evalue > $opts{'maxEvalue'};
				}
				my $identity=nearest(.01, $hsp->frac_conserved('total')*100);
	
				
				if($opts{'minIdentity'}){
					last if $identity < $opts{'minIdentity'};
				}
				$identity=$identity."%";
				
				#$hash{$result->query_name}{des}=$hit->description;
				my $des=$hit->description;
				#$hash{$result->query_name}{evalue}=$hsp->evalue;
				my $evalue=$hsp->evalue;
				#$hash{$result->query_name}{hitname}=$hit->name;
				my $hitname=$hit->name;
				$hitname=~s/^\s*//g;
				$hitname=~s/\s*$//g;
				#$hash{$result->query_name}{strand}=$hit->strand("query")==1?"+":"-";
				my $strand;
				if($algorithm=~/blastx/i){
					$strand=$hit->strand("query")==1?"+":"-";
				}else{
					$strand=" ";
				}
				
				if($opts{org} =~ /^ko$|^map$/i){
					while($des =~ /\s+(K\d{4,6})\s+/g ){
						#$hash{$result->query_name}{ko}=$1;
						my $ko=$1;
						push(@kos,$ko);	
						#print "$ko\n";			
						if(exists($kos->{"ko:".$ko})){
							foreach my $p (keys(%{$kos->{"ko:".$ko}})){
								#print "$p\n";	
								push(@paths,$p);
								push(@{$pathway->{$p}{'kos'}},"ko:".$ko);
								push(@{$pathway->{$p}{'seqs'}},$query_name);					
							}
						}else{
							warn("warn:$ko is not in all pathways in this organism $opts{org},if your database is newest, this is ok ... \n");
						}
					}
				}else{
					if(exists($kos->{$hitname})){
						foreach my $p (keys(%{$kos->{$hitname}})){
							push(@paths,$p);
							push(@{$pathway->{$p}{'kos'}},$hitname);
							push(@{$pathway->{$p}{'seqs'}},$query_name);
						}
					}else{
						warn("warn:$hitname is not in all pathways in this organism $opts{org},if your database is newest, this is ok ... \n");
					}
				}			
				
				
				while($des =~ /[\[\(](EC:[\d\.\-\s\,]+)[\]\)]/g){
					#$hash{$result->query_name}{ec}=$1;
					my $ec=$1;
					push(@ecs,$ec);
				}
				my $paths_ref=&uniq(\@paths);
				my $kos_ref=&uniq(\@kos);
				my $ecs_ref=&uniq(\@ecs);
				push(@quiery_ko,@$kos_ref);
				#print $result->query_name."\t".$hit->name."\t".$hash{$result->query_name}{ko}."\t".$hash{$result->query_name}{ec}."\t".$hit->description."\t".$hash{$result->query_name}{strand}."\t".$hsp->strand('hit')."\t".$hsp->evalue."\n";
				print KEGG "$query_name\t$hitname\thttp:\/\/www.genome.jp\/dbget-bin\/www_bget\?$hitname\t$des\t$evalue\t$score\t$strand\t$identity\t$query_length\t".sprintf("%.2f",$query_coverage*100)."%"."\t$hit_length\t".sprintf("%.2f",$hit_coverage*100)."%"."\t".join(";",@$kos_ref)."\t".join(";",@$ecs_ref)."\t".join(";",@$paths_ref)."\n";
				#print "$query_name\t$hitname\t$des\t$evalue\t$strand\t$identity\t".join(";",@$kos_ref)."\t".join(";",@$ecs_ref)."\t".join(";",@$paths_ref)."\n";
						
			}
			$seqkos{$query_name}=&uniq(\@quiery_ko);
		}
	}
}
close KEGG;

##########################################################
# output pathways
open(PATHWATY,"> $opts{o}/pathway_table.xls" ) || die "Can't open $opts{o}/pathway_table.xls\n";
print PATHWATY "PathWay\tPathway_definition\tnumber_of_seqs\tseqs_kos/genes_list\tpathway_imagename\n";
warn("outputing Pathway table ...\n");
#print (keys %$pathway);

foreach my $p (keys(%$pathway)){
	# if($p eq "path:ko00920"){
	#	  print "test";
	 #}

	next unless exists($pathway->{$p}{'kos'});
	#print $p."****\n";
	my $kolist=&uniq(\@{$pathway->{$p}{'kos'}});
	my $pathway_id = $p;
	$pathway_id=~s/^path://ig;
	my %skippath;
	$skippath{'ko01100'}=1;
	my $pathfile;
	if(exists($skippath{$pathway_id})){
	    $pathfile="http://www.kegg.jp/pathway/$pathway_id";   
	}else{
	    $pathfile="http://www.kegg.jp/pathway/$pathway_id\+".join("+",@$kolist);
	}
	

	my $imgname=&getimgname($p);
	my $htmlfile=&getimgname1($p);	
	
	if ($remote eq "yes"){	
		warn("Geting pathway image from  $pathfile  ...\n");
		&savekegg($pathfile,$imgname,$htmlfile);
	}else{
		&Mark_Pathway_local($imgname, $htmlfile, $kolist);
	}	

	
	
	my $seqlist=&uniq(\@{$pathway->{$p}{'seqs'}});
	
	my $seq_ko_list;
	
	foreach my $n (@$seqlist){
		$seq_ko_list.=$n."(".join(",",@{$seqkos{$n}}).");";
	}
	
	print PATHWATY "$p\t".$pathway->{$p}{'definition'}."\t".scalar(@$seqlist)."\t".$seq_ko_list."\t".$imgname."\n"; 
}
close PATHWATY;
warn("All done!\n");


######################################################
sub Mark_Pathway_local(){
	my $img_name = shift;
	my $html_name= shift;	
	my $ko_list = shift;
	
	
	
	my $img_file = "";
	my $html_file = "";
	
	my %ko_annot;
	foreach(@$ko_list){
		my $ko= $_;
		$ko =~ s/ko://;
		$ko_annot{$ko}=1;
	}
	
	if ($opts{org} eq "ko"){
		$img_file = $local_img_dir."/ko/".$img_name;
		$html_file = $local_img_dir."/ko/".$html_name;
	}elsif($opts{org} eq "map"){
		$img_file = $local_img_dir."/map/".$img_name;
		$html_file = $local_img_dir."/map/".$html_name;
	}else{
		die "can not support $opts{org}";
	}
	warn("annot $img_name \n");
	
	if(-e $img_file){
	       
	}else{
	 return 1;
	}
	
	#print $html_file."\n";

	my $datapage=HTML::TreeBuilder->new_from_file($html_file);
	my @data=$datapage->find_by_attribute("shape","rect");		

	
	my $image = GD::Image->newFromPng($img_file);
	my $red = $image->colorAllocate(255,1,1); 
	my $black = $image->colorAllocate(0,0,0);
	my $width = 2;
	$image->setThickness($width);
	foreach(@data){
		my $position = $_ ->attr("coords");
		my @coords = split(/,/,$position);
		my $href = $_ ->attr("href");
		
		my $isanno = 0;
		
		if($href =~ /http.*\?(.*)$/){
			my @ko_rect = split(/\+/,$1);
			foreach(@ko_rect){
				if(exists $ko_annot{$_}){
					$isanno = 1;
				}
			}
			
		}else{
			next;			
		}
		
		
		if($isanno == 1){
		        #$image->rectangle($coords[0],$coords[1],$coords[2],$coords[3],$black);
			$image->rectangle($coords[0],$coords[1],$coords[2],$coords[3],$red);
		}
		# $image->rectangle($coords[0],$coords[1],$coords[2]+2,$coords[3],$red);
		# $image->rectangle($coords[2]-2,$coords[1],$coords[2],$coords[3],$red);
		# $image->rectangle($coords[0],$coords[3]-2,$coords[2],$coords[3],$red);
	}
	
		 
	open (PNG,"> $opts{o}/$img_name");
	binmode PNG;
	print PNG $image->png;
	close PNG;
	system ("cp $html_file $opts{o}/$html_name");	
	
}

######################################################
#
sub savekegg(){
	my $url=shift;
	my $imgfile=shift;
	my $htmlfile=shift;
	try{
		my $response =$ua->get($url);
		my $filepath=$opts{o}."/".$imgfile;
		my $html;
		if($response->is_success){
		  $html = $response->decoded_content;
		  $html=&formathtml($html);
	 	   
		   my $imgurl;
		   #try{
		   		my $datapage=HTML::TreeBuilder->new_from_content($html);		
		   		my @data=$datapage->find_by_attribute("usemap","#mapdata"); 
		   		
		   			$imgurl=$data[0]->attr("src") or return &savekegg($url,$imgfile,$htmlfile);
		   		#}catch{
		   			#warn("Parsing html error! retrying ... \n");
		  			
		   		#}
		   		$datapage->delete();
		   #}catch{
		   		 #warn("Parsing html error! retrying ... \n");
		  		 #&savekegg($url,$imgfile,$htmlfile);
		   #}
		   
		    $html=&formathtml1($html,$imgfile);
		   	 open (HTML,"> $opts{o}/$htmlfile") or die "Can't create file $opts{o}/$htmlfile\n";
	 		 print HTML  $html;
	 		 close HTML;
	 		 
		 #  getstore($imgurl,$filepath) or &savekegg($url,$imgfile,$htmlfile);
		 	my $r=$ua->mirror($imgurl,$filepath);
		 	unless($r->is_success||$r->code eq '304'){
		 		 warn("Saveing image file error!:".$r->status_line." retrying get from $imgurl ... \n");
		 		 return &savekegg($url,$imgfile,$htmlfile);
		 	}
		 	if($r->code eq '304'){
		 		warn("It seems to image haven't being coloring,download uncolored image form $imgurl.This problem is form kegg server.\n");
		 	}
	
		}else{
		   if($response->code eq '414' || $response->code eq '500'){
		   		warn($response->status_line."  ,server not accept or error,skiping ... \n");
		   }else{
		   		warn($response->status_line." retrying ... \n");		   
		   		return  &savekegg($url,$imgfile,$htmlfile);
		   }
			
		   
		}
	}catch{
		warn("Server connection serious error:$_,Geting pathway image from  $url  ...\n");
		return &savekegg($url,$imgfile,$htmlfile);
	}
}

sub formathtml1(){
	my $htm =shift;
	my $imgname=shift;	
	$htm =~ s/<img src\=\".*\" usemap\=\"#mapdata\" border\=\"0\" \/>/<img src\=\"$imgname\" usemap\=\"#mapdata\" border\=\"0\" \/>/g;
	return $htm;
}

sub formathtml(){
	my($htm)=@_;
	$htm =~ s/\"\//\"http\:\/\/www.kegg.jp\//g;
	$htm =~ s/\'\//\'http\:\/\/www.kegg.jp\//g;
	return $htm;
}

sub getimgname(){
	my $path=shift;
	my @a=split(":",$path);
	return $a[1].".png";
}

sub getimgname1(){
	my $path=shift;
	my @a=split(":",$path);
	return $a[1].".html";
}

####################################
# get all pathways from database
sub getpathways(){
	my $org=shift;
	my $pw=$dbh->prepare(<<SQL
select class,definition from pathway_$org;
SQL
			  );
	$pw->execute();
	my $ref = $pw->fetchall_hashref('class');
	$pw->finish;
	return $ref;
}


########################################
# get all kos/genes from database
sub getpathwaykos(){
	my $org=shift;
	my %kolist;
	my $mm;
	if($opts{org} =~ /^ko$|^map$/i){
		$mm=$dbh->prepare(<<SQL
	select * from ko_pathway_$org;
SQL
			  );
	}else{		
		$mm=$dbh->prepare(<<SQL
	select * from gene_pathway_$org;
SQL
			  );
	}
	
	$mm->execute();
	my $ref = $mm->fetchall_hashref('id');
	foreach my $ids ( keys(%$ref) ) {
			#push(@{$list->{$ref->{$id}->{'ko'}}},$ref->{$id}->{'pathway'});
			#print $ref->{$ids}->{'ko'}."\t".$ref->{$ids}->{'pathway'}." $ids\n";
			my $k=$ref->{$ids}->{'ko'};
			my $p=$ref->{$ids}->{'pathway'};
			$kolist{$k}{$p}=1;
	}
	$mm->finish;
	return \%kolist;
}

###########################################
#fresh database from remote kegg server
sub freshdatabase(){
	my $org=shift;
	warn("Freshing database from kegg netwok,please wating ...\n");
	$dbh->do(<<SQL
	drop table if exists pathway_$org;
SQL
 		 );
	$dbh->do(<<SQL
		 CREATE TABLE  pathway_$org( 
 			 id  INTEGER PRIMARY KEY ASC, 
			 class  varchar(50) NOT NULL,
			 definition	 varchar(10) NOT NULL
		 );
SQL
 );
 	$dbh->do(<<SQL
	CREATE INDEX IF NOT EXISTS i_pathway_$org\_class ON pathway_$org(class);
SQL
 		 );
 	my $insert;
    if($opts{org} =~ /^ko$|^map$/i){
    	  	$dbh->do(<<SQL
	drop table if exists ko_pathway_$org;
SQL
 		 );
 	$dbh->do(<<SQL
	 CREATE TABLE ko_pathway_$org( 
 			 id  INTEGER PRIMARY KEY ASC, 
			 pathway  varchar(50) NOT NULL,
			 ko	varchar(10) NOT NULL
	);
SQL
 		 );	 
	 $dbh->do(<<SQL
	 CREATE INDEX IF NOT EXISTS i_ko_pathway_$org ON ko_pathway_$org(pathway);
SQL
 		 );	
    }else{
		   $dbh->do(<<SQL
	drop table if exists gene_pathway_$org;
SQL
 		 );
 			$dbh->do(<<SQL
	 CREATE TABLE gene_pathway_$org( 
 			 id  INTEGER PRIMARY KEY ASC, 
			 pathway  varchar(50) NOT NULL,
			 ko	varchar(10) NOT NULL
	);
SQL
 		 );	 
	 		$dbh->do(<<SQL
	 CREATE INDEX IF NOT EXISTS i_gene_pathway_$org ON gene_pathway_$org(pathway);
SQL
 		 );	
 	#$dbh->commit;
 
   }
 	$insert = $dbh->prepare(<<SQL
INSERT INTO pathway_$org(class,definition) VALUES (?,?);
SQL
);
	warn("getting pathway list ....\n");
	
 	my $pathway=&listpathways($org);
 	foreach my $p (@$pathway){
 		my $n=$p->{'entry_id'};
 		my $m=$p->{'definition'};
 		#print "$n $m \n";
 		$insert->execute($n,$m);
 	}
	#$dbh->commit;
	if($opts{org} =~ /^ko$|^map$/i){
		$insert = $dbh->prepare(<<SQL
INSERT INTO ko_pathway_$org(pathway,ko) VALUES (?,?);
SQL
); 
	}else{
		$insert = $dbh->prepare(<<SQL
INSERT INTO gene_pathway_$org(pathway,ko) VALUES (?,?);
SQL
); 
	}
	
	warn("getting kos/genes list for each pathway ....\n");
	
	foreach my $p (@$pathway){
		warn("getting kos/genes list for $p->{'entry_id'} ....\n");		
 		my $kos=&listkos($p->{'entry_id'});
 		foreach my $x (@$kos){
 			$insert->execute($p->{'entry_id'},$x);	
 		}
 	}
 	#undef $service;
 	#$dbh->commit;
}

#########################################
#list all pathways for organism
sub listpathways(){
    my $organism = shift;
	try{
		my $response=$ua->get("http://rest.kegg.jp/list/pathway/$organism");
		if($response->is_success){
			my $result=$response->decoded_content;
			my @lines=split(/\n+/,$result);
			my @a;
			foreach(@lines){
				if(/^(.+)\t+(.*)$/){
					my %h=('entry_id'=>$1,'definition'=>$2);
					push(@a,\%h);
				}
			}
			return \@a;
		}else{
			warn "Server response error:".$response->status_line."\n";
			warn("Server return error,retrying getting $organism ...\n");
			return &listpathways($organism);
		}
	}catch{		
		warn("Server connection serious error:$_,retrying getting $organism ...\n");
		return &listpathways($organism);
	}	
}

#####
#get all kos for a pathway
sub listkos(){
	my $pathway_id=shift;

	my $response;
        try{
		if($opts{org} =~ /^ko$|^map$/i){
			$response=$ua->get("http://rest.kegg.jp/link/ko/$pathway_id");
		
		}else{
			$response=$ua->get("http://rest.kegg.jp/link/genes/$pathway_id");
		}
		
		if($response->is_success){
				my $result=$response->decoded_content;
				my @lines=split(/\n+/,$result);
				my @a;
				foreach(@lines){
					if(/^(.+)\t+(.*)$/){
						push(@a,$2);
					}
				}
				return \@a;
		}else{
				warn "Server response error:".$response->status_line."\n";
				warn("Server return error,retrying  getting $pathway_id...\n");
				return &listkos($pathway_id);
		}		
	}catch{
		warn("Server connection serious error:$_,retrying  getting $pathway_id...\n");
		return &listkos($pathway_id);
    }
}

sub uniq {
	my $array = shift;
	my %hash = map { $_ => 1 } @$array;
	my @uniq_array = sort( keys %hash );
	return \@uniq_array;
}
