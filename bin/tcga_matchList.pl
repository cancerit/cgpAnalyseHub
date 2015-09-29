#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Cwd qw(getcwd abs_path);
use Const::Fast qw(const);
use File::Find;
use Data::Dumper;
use Capture::Tiny qw(capture);


const my $RG_LINE => q{samtools view -H %s | grep '^@RG' | head -n 1};

# ref, output, threads
const my $BWA_MEM => q{bwa_mem.pl -r %s -s %s -t %d -o %s %s};

# have to use fixed vagrent files /nfs/cancer_ref02/human/GRCh37d5/vagrent/e75/
const my $PINDEL => q{pindel.pl -st WXS -as GRCh37 -sp Human -e NC_007605,hs37d5,GL%% -c %d -r %s -sf %s/softRules.lst -u %s/pindel_np.gff3.gz -s %s/simpleRepeats.bed.gz -f %s/pulldownRules.lst -g %s/codingexon_regions.indel.bed.gz -n %s -t %s -o %s};
const my $CAVEMAN => q{caveman.pl -tc cn_tum.bed -nc cn_norm.bed -td 5 -nd 2 -t %d -r %s/genome.fa.fai -ig %s/caveman/ucscHiDepth_0.01_merge1000_no_exon.bed -b %s/caveman/flagging/ -u %s/caveman/unmatched -ab %s -s HUMAN -sa hs37d5 -st pulldown -np WXS -tp WXS -in empty.bed -nb %s -tb %s -o %s};
const my $VAGRENT => q{AnnotateVcf.pl -t -c %s -i %s -o %s};

my %find_dispatch = ( mapping => \&bwa_find,
                      pindel => \&compare_find,
                      caveman => \&compare_find,
                      annot => \&annot_find,
                    );

my %alg_dispatch = (mapping => \&bwa_mem,
                    pindel => \&pindel,
                    caveman => \&caveman,
                    annot => \&vagrent,
                    );

# when true means that you need both tumour and normal to run
my %paired_algs = ( mapping => 1, # no point mapping data if unmatched
                    pindel => 1,
                    caveman => 1,
                    annot => 0,   # annotation only needs a tumour result (as it is a comparison)
                    );


my ($base_path, $alg, @alg_args) = @ARGV;

my %donors; # this is global, be careful
file_list($alg, $base_path);

my $comparisons = 0;
for my $donor(keys %donors) {
  if($paired_algs{$alg} == 1) {
    next unless(exists $donors{$donor}{'normal'} && exists $donors{$donor}{'tumour'});
  }
  my @inputs = ();
  push @inputs, @{$donors{$donor}{'tumour'}};
  if($paired_algs{$alg} == 1) {
    unshift @inputs, $donors{$donor}{'normal'};
    $comparisons += (scalar @inputs) - 1;
  }
  else {
    $comparisons += scalar @inputs;
  }

  $alg_dispatch{$alg}->(@alg_args, @inputs);
}

my $donor_count = scalar (keys %donors);

warn "$donor_count Donors\n";
warn "$comparisons Jobs\n";

sub file_list {
  my ($alg, $path) = @_;
  my $cwd = getcwd; # we don't want the absolute path as we need the info from the link area
  $path = $cwd.'/'.$path if($path !~ m/^\//);
  find( { wanted => $find_dispatch{$alg}, preprocess => \&dir_exclude }, $path);
}

sub bwa_find {
  return unless($_ =~ m/[.]bam$/);
  my $file = $File::Find::name;
  return if(-s $file == 0);
  my @bits = split '/', $file;
  my $donor = $bits[-3];
  my $type = $bits[-2];
  if($type eq 'normal') {
    $donors{$donor}{$type} = $file;
  }
  else {
    push @{$donors{$donor}{$type}}, $file;
  }
}

sub compare_find {
  return unless($_ =~ m/[.]bas$/);
  my $file = $File::Find::name;
  return if(-s $file == 0);
  $file =~ s/[.]bas$//;
  my @bits = split '/', $file;
  my $donor = $bits[-4];
  my $type = $bits[-3];
  if($type eq 'normal') {
    if(exists $donors{$donor}{$type}) {
      warn "Multiple normals available for donor '$donor'\n";
      return;
    }
    $donors{$donor}{$type} = $file;
  }
  else {
    push @{$donors{$donor}{$type}}, $file;
  }
}

sub annot_find {
  return unless($_ =~ m/[.]flagged([.]muts)?[.]vcf[.]gz[.]tbi$/);
  my $file = $File::Find::name;
  return if(-s $file == 0);
  $file =~ s/[.]tbi$//;
  my @bits = split '/', $file;
  my $donor = $bits[-5];
  my $type = $bits[-4];
  my $bc = $bits[-2];
  if($type eq 'normal') {
    warn "Analysis result found for normal BC '$bc' of donor '$donor'\n";
  }
  else {
    push @{$donors{$donor}{$type}}, $file;
  }
}

sub dir_exclude {
  return () if($File::Find::dir =~ m/tmp[^\/]+$/ || $File::Find::dir =~ m/logs([^\/]+)?$/);
  return @_;
}

sub vagrent {
  my ($vagrent_cache, @vcfs) = @_;
  for my $vcf_in(@vcfs) {
    my $out = $vcf_in;
    $out =~ s/[.]vcf[.]gz$/.annot.vcf/;
    my $exist_res = $out.'.gz.tbi';
    next if(-e $exist_res && -s $exist_res > 0);
    my $command = sprintf $VAGRENT, $vagrent_cache, $vcf_in, $out;
    print $command,"\n";
  }
}

sub pindel {
  my ($ref, $pindel_base, $vagrent_base, $threads, @bams) = @_;
  my $normal_bam = shift @bams;
  for my $tumour_bam(@bams) {
    my $sample = sample_from_bam($tumour_bam);
    my $outbase = $tumour_bam;
    $outbase =~ s|/mapped/.+.bam$||;
    $outbase .= '/pindel/'.$sample;
    my ($stdout, $stderr, $exit) = capture {system("ls -1 $outbase/*.flagged.vcf.gz.tbi");};
    chomp $stdout;
    next if($stdout =~ m/[.]flagged[.]vcf[.]gz[.]tbi$/ && !-e "$outbase/tmpPindel");

    my $command = sprintf $PINDEL, $threads, $ref, $pindel_base, $pindel_base, $pindel_base, $pindel_base, $vagrent_base, $normal_bam, $tumour_bam, $outbase;
    print $command,"\n";
  }
}

sub caveman {
  my ($ref_base, $vagrent_base, $threads, @bams) = @_;
  my $normal_bam = shift @bams;
  for my $tumour_bam(@bams) {
    my $sample = sample_from_bam($tumour_bam);
    my $outbase = $tumour_bam;
    $outbase =~ s|/mapped/.+.bam$||;
    $outbase .= '/caveman/'.$sample;

    my ($stdout, $stderr, $exit) = capture {system("ls -1 $outbase/*.flagged.muts.vcf.gz.tbi");};
    chomp $stdout;
    next if($stdout =~ m/[.]flagged[.]muts[.]vcf[.]gz[.]tbi$/ && !-e "$outbase/tmpCaveman");

    my $command = sprintf $CAVEMAN, $threads, $ref_base, $ref_base, $ref_base, $ref_base, $vagrent_base, $normal_bam, $tumour_bam, $outbase;
    print $command,"\n";
  }
}

sub bwa_mem {
  my ($ref, $threads, @bams) = @_;
  for my $bam(@bams) {
    my $sample = sample_from_bam($bam);
    my ($type) = $bam =~ m/\/(tumour|normal)\//;
    my $outbase = $bam;
    $outbase =~ s/\/(tumour|normal)\/.+/\/$1\/mapped/;
    my $existing_bas = $outbase.'/'.$sample.'.bam.bas';
    next if(-e $existing_bas && -s $existing_bas != 0);
    my $command = sprintf $BWA_MEM, $ref, $sample, $threads, $outbase, $bam;
    print $command,"\n";
  }
}


sub sample_from_bam {
  my ($bam) = @_;
  my $cmd = sprintf $RG_LINE, $bam;
  my ($sample) = `$cmd` =~ m/\tSM:([^\t\n]+)/;
  die "Unable to find sample name in $bam\n" if(!defined $sample || (length $sample) == 0);
  return $sample;
}
