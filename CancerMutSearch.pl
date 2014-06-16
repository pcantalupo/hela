#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use BamMutSearch;
use Getopt::Long;

my $folder;
my $bam_file = "";
my $genome;
my $out;
my $mutations;
my $aid;
my $r = GetOptions('if|id_folder=s' => \$folder,
	'b|bam=s' => \$bam_file,
	'a|aid=s' => \$aid,
	'g|genome=s' => \$genome,
	'o|out=s' => \$out,
	'm|mutation=s' => \$mutations);
die "Need genome file in fasta format\n" if !-e $genome;
die "Need folder -if|id_folder or bam file -b|bam and -a|aid\n" if !(-e $bam_file && $aid) && !-d $folder;
die "Need mutation file -m|mutation\n" if ! -e $mutations;
die "Need out -o|out for file name\n" if !$out;
my $seqI = new Bio::SeqIO(-file => $genome,-format => 'fasta');
$genome = $seqI->next_seq;

my @mutations;
open MUT, $mutations;
while (<MUT>) {
	next if /^#/;
	chomp;
	s/#.*$//;
	s/\s+$//;
	push (@mutations, $_);
	#my @mutations = map {chomp; $_} <MUT>;
}
close MUT;
my %bam;
if(defined($folder)){
	opendir DIR, $folder;
	my @files = readdir(DIR);
	closedir DIR;
	foreach my $f (@files){
		next if $f !~ /\w{8}\-\w{4}\-\w{4}\-\w{4}\-\w{12}/;
		my $bamSearch = new BamMutSearch($genome,@mutations);
		$bamSearch->BamSearch("$folder/$f/vrs_pipeline/ppr.qc.fq.vrs.hg19.bam");
		$bam{$f} = $bamSearch->get('mut');
	}
}
else{
	my $bamSearch = new BamMutSearch($genome,@mutations);
	$bamSearch->BamSearch($bam_file);
	$bam{$aid} = $bamSearch->get('mut');
}
my @stdout_keys = qw(total_align total_poly perc_poly_w_mutation range);
my %stdout_print;
open OUT, ">$out";
my @header = qw(mutation index bp);
foreach my $analysis (sort keys %bam){
	push(@header,$analysis."_delta",$analysis."_total");
}
print OUT join("\t",@header),$/;
foreach my $m (@mutations){
	next if ($m =~ /(\d+)-(\d+)/);   # skip Range
	my @toPrint = ($m);
	push(@toPrint, ($m =~ /(\d+)(\w)/));
	foreach my $a (sort keys %bam){
		$stdout_print{$a}{'total_poly'} += $bam{$a}{$toPrint[1]}{'total'};
		$stdout_print{$a}{'perc_poly_w_mutation'} += $bam{$a}{$toPrint[1]}{'delta'};
		$stdout_print{$a}{'total_align'} = $bam{$a}{range}{'total'};
		$stdout_print{$a}{'range'} = $bam{$a}{range}{'delta'};
		push(@toPrint,$bam{$a}{$toPrint[1]}{'delta'},
			$bam{$a}{$toPrint[1]}{'total'});
	}
	print OUT join("\t",@toPrint),$/;
}

print join("\t",'aid',@stdout_keys),$/;
foreach my $a (keys %stdout_print){
	$stdout_print{$a}{'perc_poly_w_mutation'} = 0 if !defined($stdout_print{$a}{'perc_poly_w_mutation'});
	$stdout_print{$a}{'range'} = 0 if !defined($stdout_print{$a}{'range'});
	my $ppwm = $stdout_print{$a}{'perc_poly_w_mutation'}/$stdout_print{$a}{'total_poly'}
		if defined($stdout_print{$a}{'total_poly'}) && $stdout_print{$a}{'total_poly'} > 0;
	$ppwm = 0 if !$ppwm;
	my $percrange = $stdout_print{$a}{'range'}/$stdout_print{$a}{'total_align'}
		if defined($stdout_print{$a}{'total_align'}) && $stdout_print{$a}{'total_align'} > 0;
	$percrange = 0 if !$percrange;
	$stdout_print{$a}{'total_align'} = 0 if !$stdout_print{$a}{'total_align'};
	print join("\t",$a,$stdout_print{$a}{'total_align'},
		$stdout_print{$a}{'total_poly'},$ppwm,$percrange),$/;
}
