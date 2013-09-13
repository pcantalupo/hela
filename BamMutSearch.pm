package BamMutSearch;
use strict;
use warnings;
use Bio::SeqIO;

sub new{
	my $class = shift;
	my $self = {};
	$self->{'genome'} = shift;
	foreach my $m (@_){
		my ($index,$bp) = ($m =~ /(\d+)(\w)/);
		$self->{'mut'}{$index}{'bp'} = $bp;
		$self->{'mut'}{$index}{'delta'} = 0;
		$self->{'mut'}{$index}{'total'} = 0;
	}
	bless $self, $class;
	return $self;
}
sub BamSearch{
	my $self = shift;
	my $bam = shift;
	my $id = $self->{'genome'}->id;
	my @lines = `samtools view $bam | grep '$id'`;

	foreach my $l (@lines){
		my @cols = split "\t", $l;
		my $r5p = $cols[3];
		my $r3p = $cols[3]+length($cols[9])-1;
		$self->{'mut'}{'gt2500'}{'delta'}++ if $r5p > 2500;
		$self->{'mut'}{'gt2500'}{'total'}++;
		foreach my $m (keys %{$self->{'mut'}}){
			next if $m eq 'gt2500';
			if($m >= $r5p && $m <= $r3p){
				my @read = split '',$cols[9];
				my $m5p = $m-$r5p;
				if($self->{'mut'}{$m}{'bp'} eq $read[$m5p]){
					my $subseq = $self->{'genome'}->subseq($r5p,$r3p);
					$self->{'mut'}{$m}{'delta'}++;
				}
				$self->{'mut'}{$m}{'total'}++;
			}
		}
	}
}
sub get{
	my $self = shift;
	my $key = shift;
	return $self->{$key};
}
1;

=cut
my $genome;
my $bam;

my $r = GetOptions('g|genome=s' => \$genome, 'b|bam=s' => \$bam);

die "Need to define g|genome with fasta file\nNeed to defined b|bam with one bam file" if !$genome || !$bam;
die "Define mutations psotions without a flag ex 281C \n" if !scalar(@ARGV);
my $seqI = new Bio::SeqIO(-file => $genome,-format => 'fasta');
my $seq = $seqI->next_seq;
$seqI = '';
my @muts = @ARGV;
my $id = $seq->id;
my %mutation;
my %out;

foreach my $m (@muts){
	my ($index,$delta) = ($m =~ /(\d+)(\w)/);
	$mutation{$index} = $delta;
	$out{$index}{'total'} = 0;
	$out{$index}{'delta'} = 0;
}

my @lines = `samtools view $bam | grep '$id'`;

foreach my $l (@lines){
	my @cols = split "\t", $l;
	my $r5p = $cols[3];
	my $r3p = $cols[3]+length($cols[9])-1;
	foreach my $m (keys %mutation){
		if($m >= $r5p && $m <= $r3p){
			my @read = split '',$cols[9];
			my $m5p = $m-$r5p;
			if($mutation{$m} eq $read[$m5p]){
				my $subseq = $seq->subseq($r5p,$r3p);
				$out{$m}{'delta'}++;
				print join("\t",$m5p,$cols[9]),$/;
				print join("\t",$r5p,$subseq),$/;
			}
			$out{$m}{'total'}++;
		}
	}
}
foreach my $k (keys %out){
	print join("\t",$k,$out{$k}{'delta'},$out{$k}{'total'}),$/;
}
1;
