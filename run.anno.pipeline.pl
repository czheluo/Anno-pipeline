#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fout,$fa,$wsh,$step,$stop);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"out:s"=>\$fout,
	"fa:s"=>\$fa,
	"wsh:s"=>\$wsh,
	"step:s"=>\$step,
	"stop:s"=>\$stop,
			) or &USAGE;
&USAGE unless ($fout);
$step||=1;
$stop||=-1;
$fout=ABSOLUTE_DIR($fout);
mkdir $fout if (!-d $fout);
mkdir $wsh if (!-d $wsh);
my $ann||="$fout/run_anno";
mkdir $ann if (!-d $ann);
my $tmp||="$ann/tmp";
mkdir $tmp if (!-d $tmp);
my $blast||="$tmp/blast";
mkdir $blast if (!-d $blast);
my $qsub="perl /mnt/ilustre/centos7users/dna/.env/bin/qsub-slurm.pl";
my $faa=ABSOLUTE_DIR($fa);
open LOG,">$fout/$wsh/anno.$BEGIN_TIME.log";
if ($step == 1) {
	print LOG "########################################\n";
	print LOG "split fa\n"; my $time=time();
	print LOG "########################################\n";
	open SH1,">$wsh/splitfa.sh";
	print SH1 "cd $blast && perl $Bin/bin/splitfasta.pl $faa 100 && ";
	close SH1;
	my $job="$qsub $wsh/splitfa.sh ";
	`$job`;
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}

if ($step == 2) {
	print LOG "########################################\n";
	print LOG "blast to nr kegg go database\n"; my $time=time();
	print LOG "########################################\n";
	##blast to nr
	my @split=glob("$blast/*.split");
	my $nr=1;
	foreach my $split (@split) {
		open SH2,">$wsh/blastnr.$nr.sh";
		my $fsplit=basename($split);
		my $name=basename($split);
		$name=~s/\.split//g;
		print SH2 "perl $Bin/bin/mblast.pl -p blastx -i $blast/$fsplit -d nr -o $blast/$name\_vs_nr.xml -e 1E-5 -v  5 && ";
		print SH2 "perl $Bin/bin/Blast2table -format 10 -xml $blast/$name\_vs_nr.xml >$blast/$name\_vs_nr.xml.table.xls ";
		close SH2;
		my $job="$qsub --Resource mem=80G --CPU 16 $wsh/blastnr.$nr.sh ";
		`$job`;
		$nr++;
	}
	##blast to string
	my $ns=1;
	foreach my $split (@split) {
		my $fsplit=basename($split);
		my $name=basename($split);
		$name=~s/\.split//g;
		open SH3,">$wsh/blaststring.$ns.sh";
		print SH3 "perl $Bin/bin/mblast.pl -p blastx -i $blast/$fsplit -d string -o $blast/$name\__vs_string.xml -e 1E-5 -v  5 \n";
		close SH3;
		my $job="$qsub --Resource mem=80G --CPU 16 $wsh/blaststring.$ns.sh ";
		`$job`;
		$ns++;
	}
	##blast to Swissprot
	my $nsw=1;
	foreach my $split (@split) {
		my $fsplit=basename($split);
		my $name=basename($split);
		$name=~s/\.split//g;
		open SH4,">$wsh/blastswissprot.$nsw.sh";
		print SH4 "perl $Bin/bin/mblast.pl -p blastx -i $blast/$fsplit -d swissprot -o $blast/$name\_vs_swissprot.xml -e 1E-5 -v 5 && ";
		print SH4 "perl $Bin/bin/Blast2table -format 10 -xml $blast/$name\_vs_swissprot.xml >$blast/$name\_vs_swissprot.xml.table.xls ";
		close SH4;
		my $job="$qsub --Resource mem=80G --CPU 16 $wsh/blastnr.$nsw.sh ";
		`$job`;
		$nsw++;
	}
	##blast to kegg
	my $nk=1;
	foreach my $split (@split) {
		my $fsplit=basename($split);
		my $name=basename($split);
		$name=~s/\.split//g;
		open SH5,">$wsh/blastswissprot.$nk.sh";
		print SH5 "perl $Bin/bin/mblast.pl -p blastx -i $blast/$fsplit -d ko.pep.fasta -o $blast/$name\_vs_kegg.xml -e 1E-5 -v  5 && ";
		print SH5 "python $Bin/bin//kobas/kobas2/scripts/annotate.py -i $blast/$name\_vs_kegg.xml -t blastout:xml -s ko -o $blast/$name\_vs_kegg.xml.kobas ";
		close SH5;
		my $job="$qsub --Resource mem=80G --CPU 16 $wsh/blastnr.$nk.sh ";
		`$job`;
		$nk++;
	}
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
my ($annnr,$annstr,$annswi,$annke,$annotation);
my $db="/mnt/ilustre/users/bingxu.liu/workspace/annotation/db/cog.db";
if ($step == 3) {
	print LOG "########################################\n";
	print LOG "split fa\n"; my $time=time();
	print LOG "########################################\n";
	## merge annotation/NR/
	$annotation="$ann/annotation/";
	mkdir $annotation if (!-d $annotation);
	$annnr ="$ann/annotation/NR";
	mkdir $annnr if (!-d $annnr);
	open SH6,">$wsh/mergexmlnr.sh";
	print SH6 "perl $Bin/bin/mergeBlastXml.pl $blast/$fa\*_vs_nr.xml >$annnr/$fa\_vs_nr.xml \n";
	open Out,">$annnr/$fa\_vs_nr.blasttable.xls";
	print Out "Score\tE-Value\tHSP-Len\t%-Ident\t%-Simil\tQuery-Name\tNum-Rds\tQ-Len\tQ-Begin\tQ-End\tQ-Frame\tHit-Name\tH-Len\tH-Begin\tH-End\tH-Frame\tDescription\n";
	close Out ;
	print SH6 "cat $blast/$fa*_vs_nr.xml >> $annnr/$fa\_vs_nr.blasttable.xls ";
	close SH6;
	my $job="$qsub $wsh/mergexmlnr.sh ";
	`$job`;
	## merge annotation/string/
	$annstr ="$ann/annotation/COG_KOG";
	mkdir $annstr if (!-d $annstr);
	open SH7,">$wsh/mergexmlstr.sh";
	print SH7 "perl $Bin/bin/mergeBlastXml.pl $blast/$fa\*_vs_string.xml >$annstr/$fa\_vs_string.xml \n";
	close SH7;
	my $job="$qsub $wsh/mergexmlstr.sh ";
	`$job`;
	## merge annotation/Swissprot
	$annswi ="$ann/annotation/Swissprot";
	mkdir $annswi if (!-d $annswi);
	open SH8,">$wsh/mergexmlSwi.sh";
	print SH8 "perl $Bin/bin/mergeBlastXml.pl $blast/$fa\*_vs_swissprot.xml >$annswi/$fa\_vs_swissprot.xml \n";
	open Out ," > $annswi/$fa\_vs_swissprot.blasttable.xls ";
	print Out  "Score\tE-Value\tHSP-Len\t%-Ident\t%-Simil\tQuery-Name\tNum-Rds\tQ-Len\tQ-Begin\tQ-End\tQ-Frame\tHit-Name\tH-Len\tH-Begin\tH-End\tH-Frame\tDescription \n";
	close Out ;
	print SH8 "cat $blast/$fa\*_vs_swissprot.xml >> $annswi/$fa\_vs_swissprot.blasttable.xls ";
	close SH8;
	my $job="$qsub $wsh/mergexmlSwi.sh ";
	`$job`;
	## merge annotation/kegg
	$annke ="$ann/annotation/KEGG";
	mkdir $annke if (!-d $annke);
	open SH8,">$wsh/mergexmlkegg.sh";
	print SH8 "perl $Bin/bin/mergeBlastXml.pl $blast/$fa\*_vs_kegg.xml >$annke/$fa\_vs_kegg.xml \n";
	print SH8 "perl $Bin/bin/mergeKobas.pl $blast/$fa\*_vs_kegg.xml > $annke/pathway.txt ";
	close SH8;
	my $job="$qsub $wsh/mergexmlkegg.sh ";
	`$job`;
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 4) {
	print LOG "########################################\n";
	print LOG "split fa\n"; my $time=time();
	print LOG "########################################\n";
	##String2Cog
	open SH9,">$wsh/String2Cog.sh";
	print SH9 "perl $Bin/bin/String2Cog.pl -i $annstr/$fa\_vs_string.xml -parse_id --format blastxml -o $annstr -img COG,KOG -db $db -mblast";
	my $job="$qsub $wsh/String2Cog.sh ";
	`$job`;

	##b2g4pipe
	my $anngo="$ann/annotation/GO";
	mkdir $anngo if (!-d $anngo);
	my $GO="$tmp/GO";
	mkdir $GO if (!-d $GO);
	open SH10,">$wsh/b2g4pipe.sh";
	print SH10 "cd $Bin/bin/b2g4pipe && perl $Bin/bin/mergeBlastXml.pl -f $annnr/$fa\_vs_nr.xml > $GO/$fa\_vs_nr.xml \n";
	print SH10 "java -Xms128m -Xmx80000m -cp *:ext/*: es.blast2go.prog.B2GAnnotPipe -in $GO/$fa\_vs_nr.xml -fas $faa -out $anngo/blast2go -prop $Bin/bin/b2g4pipe/b2gPipe.properties -v -annot -dat -img -annex -wiki $Bin/bin/b2g4pipe/html_template.html \n";
	print SH10 "perl $Bin/bin/mergeGO.pl $anngo/blast2go.annot > $anngo/GO.list ";
	close SH10;
	my $job="$qsub --Resource mem=50G $wsh/b2g4pipe.sh ";
	`$job`;

	#getKEGG
	my $dbkegg="/mnt/ilustre/users/bingxu.liu/workspace/annotation/db/kegg.db";
	open SH11,">$wsh/getKEGG.sh";
	print SH11 "perl $Bin/bin/getKEGGfromBlast.pl -i $annke/pathway.txt -format kobas -o $annke/pathways -org ko -use_proxy -database $dbkegg ";
	my $job="$qsub $wsh/getKEGG.sh ";
	`$job`;

	#gene-onto
	open SH12,">$wsh/gene-onto.sh";
	print SH12 "perl $Bin/bin/gene-ontology.pl -i $anngo/GO.list -l 2 -list $anngo/GO.level2.list > $anngo/level2.go.txt && ";
	print SH12 "cd $anngo && $Bin/bin/go-bar.pl -i $anngo/level2.go.txt ";
	my  $job="$qsub $wsh/gene-onto.sh ";
	`$job`;

	##mergeAnnotation.pl
	open SH13,">$wsh/megerann.sh";
	print SH13 "perl $Bin/bin/mergeAnnotation.pl -fa $faa -nr $annnr/$fa\_vs_nr.blasttable.xls -string $annstr/cog_table.xls -cog $annstr/cog.list.xls -swissport $annswi/$fa\_vs_swissprot.blasttable.xls -kobas $annke/pathway.txt -go $anngo/GO.list >$annotation/annotation.table.xls";
	my $job="$qsub $wsh/megerann.sh ";
	`$job`;
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
close LOG;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        meng.luo\@majorbio.com;
Script:			$Script
Description:

	eg: perl -int filename -out filename 

Usage:
  Options:
	-fa ref_genome.gtf.exon.fa
	-out ouput file name 
	-wsh work_sh
	-step 1
	-stop -1
	-h         Help
USAGE
        print $usage;
        exit;
}
