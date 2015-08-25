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

my ($base_path, $alg, @alg_args) = @ARGV;


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


__END__


## GETTING DATA
cd /lustre/scratch116/casm/cancer_external/admin_root
source /software/CGP/pancan/final.csh.setup
perl-5.16.3 ~kr2/GitHub/cgpAnalyseHub/bin/cghub_pull.pl -k /lustre/scratch116/casm/cgp/pancancer/downloads/pem_keys/cghub.pem -g /software/genetorrent-3.8.7/bin -o /lustre/scratch116/casm/cancer_external/downloads -b bad_analysis.txt -s cghub_COAD.tsv -t 4

perl-5.16.3 ~kr2/GitHub/cgpAnalyseHub/bin/cghub_pull.pl -k /lustre/scratch116/casm/cgp/pancancer/downloads/pem_keys/cghub.pem -g /software/genetorrent-3.8.7/bin -o /lustre/scratch116/casm/cancer_external/downloads -b bad_analysis.txt -s cghub_READ.tsv -t 4

# downloads ongoing:
ps -fu cgppipe | grep 'sh -c' | grep -v grep | cut -d '&' -f 2 | xargs tail -n 2

# lists the number of downloaded BAMs for each project/tissue type
ls -1 /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/*/*/*/*.bam | cut -f 8-9 -d '/' | sort | uniq -c

# lists the number of remapped BAMs for each project/tissue type
ls -1 /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/*/*/*/mapped/*.bas | cut -f 8-9 -d '/' | sort | uniq -c

#################

# example of this scripts use
perl-5.16.3 ~kr2/noddy/tcga_matchList.pl donors/TCGA/ mapping /lustre/scratch116/casm/cgp/pancancer/reference/genome.fa 12

ssh cgppipe@cgp2-farm-login
cd /lustre/scratch116/casm/cancer_external/admin_root
source /software/CGP/pancan/final.csh.setup
perl-5.16.3 ~kr2/noddy/tcga_matchList.pl /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/READ mapping /lustre/scratch116/casm/cgp/pancancer/reference/genome.fa 12 > map_commands_READ.txt
wc -l map_commands_READ.txt
rm -f logs/mapping/READ.*.log
bsub -oo logs/mapping/READ.%I.log -sp 100 -q long -J 'tcga_mapREAD[1-30]%50' -P analysis-cgp -n 12 -M18000 -R"span[hosts=1] select[mem>18000] rusage[mem=18000]" '/nfs/users/nfs_k/kr2/noddy/farm_idx_exec.pl map_commands_READ.txt $LSB_JOBINDEX'
bjobs DONE
rm -rf /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/READ/*/*/mapped/logs

# pindel generations
cd /lustre/scratch116/casm/cancer_external/admin_root
source /software/CGP/pancan/final.csh.setup
mkdir -p logs/pindel
perl-5.16.3 ~kr2/noddy/tcga_matchList.pl /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/READ pindel /lustre/scratch116/casm/cgp/pancancer/reference/genome.fa $SCRATCH112/PanCancerFinal/ref/pindel /nfs/cancer_ref02/human/GRCh37d5/vagrent/e75 4 > pindel_commands_READ.txt
# warnings give number of comparisons, but not accurate if some complete
wc -l pindel_commands_READ.txt
rm -f logs/pindel/READ.*.log
bsub -oo logs/pindel/READ.%I.log -sp 100 -q long -J 'tcga_pindelREAD[1-12]%100' -P analysis-cgp -n 4 -M24000 -R"span[hosts=1] select[mem>24000] rusage[mem=24000]" '/nfs/users/nfs_k/kr2/noddy/farm_idx_exec.pl pindel_commands_READ.txt $LSB_JOBINDEX'
bjobs -A DONE
rm -rf /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/READ/*/*/pindel/*/logs


# caveman generations
cd /lustre/scratch116/casm/cancer_external/admin_root
source /software/CGP/pancan/final.csh.setup
mkdir -p logs/caveman
perl -e 'print qq{1\t0\t1\n};' > empty.bed
perl -ane 'print "$F[0]\t0\t$F[1]\t5\n";' < $SCRATCH112/PanCancerFinal/ref/genome.fa.fai > cn_tum.bed
perl -ane 'print "$F[0]\t0\t$F[1]\t2\n";' < $SCRATCH112/PanCancerFinal/ref/genome.fa.fai > cn_norm.bed
perl-5.16.3 ~kr2/noddy/tcga_matchList.pl /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/READ caveman /lustre/scratch116/casm/cgp/pancancer/reference /lustre/scratch116/casm/cgp/pancancer/reference/vagrent/e75 12 > caveman_commands_READ.txt
# warnings give number of comparisons, but not accurate if some complete
wc -l caveman_commands_READ.txt
rm -f logs/caveman/READ.*.log
bsub -oo logs/caveman/READ.%I.log -sp 100 -q long -J 'tcga_cavemanREAD[1]%100' -P analysis-cgp -n 12 -M48000 -R"span[hosts=1] select[mem>48000] rusage[mem=48000]" '/nfs/users/nfs_k/kr2/noddy/farm_idx_exec.pl caveman_commands_READ.txt $LSB_JOBINDEX'
bjobs -A DONE # to best of ability one sample fails to generate any results
rm -rf /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/READ/*/tumour/caveman/*/logs


## annotation
#ls -1 /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/READ/*/tumour/*/*/*flagged*.vcf.gz.tbi
cd /lustre/scratch116/casm/cancer_external/admin_root
source /software/CGP/pancan/final.csh.setup
mkdir -p logs/vagrent
perl-5.16.3 ~kr2/noddy/tcga_matchList.pl /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/READ annot /lustre/scratch116/casm/cgp/pancancer/reference/vagrent/e75/Human.GRCh37.vagrent.cache.gz > vagrent_commands_READ.txt
wc -l vagrent_commands_READ.txt
rm -f logs/vagrent/READ.*.log
bsub -oo logs/vagrent/READ.%I.log -sp 100 -q normal -J 'tcga_vagrentREAD[1-14]' -P analysis-cgp -n 1 -M2000 -R"span[hosts=1] select[mem>2000] rusage[mem=2000]" '/nfs/users/nfs_k/kr2/noddy/farm_idx_exec.pl vagrent_commands_READ.txt $LSB_JOBINDEX'
bjobs -A DONE

## COAD MAPPING

ssh cgppipe@cgp2-farm-login
cd /lustre/scratch116/casm/cancer_external/admin_root
source /software/CGP/pancan/final.csh.setup
perl-5.16.3 ~kr2/noddy/tcga_matchList.pl /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/COAD mapping /lustre/scratch116/casm/cgp/pancancer/reference/genome.fa 12 > map_commands_COAD.txt
wc -l map_commands_COAD.txt
rm -f logs/mapping/COAD.*.log
bsub -oo logs/mapping/COAD.%I.log -sp 100 -q long -J 'tcga_mapCOAD[1-12]%300' -P analysis-cgp -n 12 -M20000 -R"span[hosts=1] select[mem>20000] rusage[mem=20000]" '/nfs/users/nfs_k/kr2/noddy/farm_idx_exec.pl map_commands_COAD.txt $LSB_JOBINDEX'
bjobs -A DONE

# pindel generations
cd /lustre/scratch116/casm/cancer_external/admin_root
source /software/CGP/pancan/final.csh.setup
mkdir -p logs/pindel
perl-5.16.3 ~kr2/noddy/tcga_matchList.pl /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/COAD pindel /lustre/scratch116/casm/cgp/pancancer/reference/genome.fa $SCRATCH112/PanCancerFinal/ref/pindel /nfs/cancer_ref02/human/GRCh37d5/vagrent/e75 4 > pindel_commands_COAD.txt
# warnings give number of comparisons, but not accurate if some complete
wc -l pindel_commands_COAD.txt
rm -f logs/pindel/COAD.*.log
bsub -oo logs/pindel/COAD.%I.log -sp 100 -q long -J 'tcga_pindelCOAD[1-6]%100' -P analysis-cgp -n 4 -M24000 -R"span[hosts=1] select[mem>24000] rusage[mem=24000]" '/nfs/users/nfs_k/kr2/noddy/farm_idx_exec.pl pindel_commands_COAD.txt $LSB_JOBINDEX'
bjobs -A DONE
rm -rf /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/COAD/*/tumour/pindel/*/logs

# caveman generations
cd /lustre/scratch116/casm/cancer_external/admin_root
source /software/CGP/pancan/final.csh.setup
mkdir -p logs/caveman
perl -e 'print qq{1\t0\t1\n};' > empty.bed
perl -ane 'print "$F[0]\t0\t$F[1]\t5\n";' < $SCRATCH112/PanCancerFinal/ref/genome.fa.fai > cn_tum.bed
perl -ane 'print "$F[0]\t0\t$F[1]\t2\n";' < $SCRATCH112/PanCancerFinal/ref/genome.fa.fai > cn_norm.bed
perl-5.16.3 ~kr2/noddy/tcga_matchList.pl /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/COAD caveman /lustre/scratch116/casm/cgp/pancancer/reference /lustre/scratch116/casm/cgp/pancancer/reference/vagrent/e75 12 > caveman_commands_COAD.txt
# warnings give number of comparisons, but not accurate if some complete
wc -l caveman_commands_COAD.txt
rm -f logs/caveman/COAD.*.log
bsub -oo logs/caveman/COAD.%I.log -sp 100 -q long -J 'tcga_cavemanCOAD[1-410]%100' -P analysis-cgp -n 12 -M48000 -R"span[hosts=1] select[mem>48000] rusage[mem=48000]" '/nfs/users/nfs_k/kr2/noddy/farm_idx_exec.pl caveman_commands_COAD.txt $LSB_JOBINDEX'
bjobs -A DONE
rm -rf /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/COAD/*/tumour/caveman/*/logs

# bad analysis, cant proceed
grep -l Exited logs/caveman/COAD.*.log | xargs -I {} bash -c 'head -n 1 {} | cut -d " " -f 40' | xargs -I {} bash -c 'grep -Fl "got 0. at /software/CGP/pancan/bin/mergeCavemanResults line 54." {}/tmpCaveman/logs/Sanger::CGP::Caveman::Implement::caveman_merge_results.0.err'

## annotation
cd /lustre/scratch116/casm/cancer_external/admin_root
source /software/CGP/pancan/final.csh.setup
mkdir -p logs/vagrent
perl-5.16.3 ~kr2/noddy/tcga_matchList.pl /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/COAD annot /lustre/scratch116/casm/cgp/pancancer/reference/vagrent/e75/Human.GRCh37.vagrent.cache.gz > vagrent_commands_COAD.txt
wc -l vagrent_commands_COAD.txt
rm -f logs/vagrent/COAD.*.log
bsub -oo logs/vagrent/COAD.%I.log -sp 100 -q normal -J 'tcga_vagrentCOAD[1-789]' -P analysis-cgp -n 1 -M2000 -R"span[hosts=1] select[mem>2000] rusage[mem=2000]" '/nfs/users/nfs_k/kr2/noddy/farm_idx_exec.pl vagrent_commands_COAD.txt $LSB_JOBINDEX'
bjobs -A DONE



## BRCA DOWNLOAD
cd /lustre/scratch116/casm/cancer_external/admin_root
source /software/CGP/pancan/final.csh.setup
perl-5.16.3 ~kr2/GitHub/cgpAnalyseHub/bin/cghub_pull.pl -k /lustre/scratch116/casm/cgp/pancancer/downloads/pem_keys/cghub.pem -g /software/genetorrent-3.8.7/bin -o /lustre/scratch116/casm/cancer_external/downloads -b bad_analysis.txt -s cghub_ly2.tsv -t 4
DONE

# start mapping
perl-5.16.3 ~kr2/noddy/tcga_matchList.pl /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/BRCA mapping /lustre/scratch116/casm/cgp/pancancer/reference/genome.fa 12 > map_commands_BRCA.txt
mkdir -p logs/mapping
wc -l map_commands_BRCA.txt
rm -f logs/mapping/BRCA.*.log
bsub -oo logs/mapping/BRCA.%I.log -sp 100 -q long -J 'tcga_mapBRCA[1-46]%300' -P analysis-cgp -n 12 -M20000 -R"span[hosts=1] select[mem>20000] rusage[mem=20000]" '/nfs/users/nfs_k/kr2/noddy/farm_idx_exec.pl map_commands_BRCA.txt $LSB_JOBINDEX'
bjobs -A 944143
rm -rf /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/BRCA/*/*/mapped/logs_bwamem_*


# pindel generations
cd /lustre/scratch116/casm/cancer_external/admin_root
source /software/CGP/pancan/final.csh.setup
mkdir -p logs/pindel
perl-5.16.3 ~kr2/noddy/tcga_matchList.pl /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/BRCA pindel /lustre/scratch116/casm/cgp/pancancer/reference/genome.fa $SCRATCH112/PanCancerFinal/ref/pindel /nfs/cancer_ref02/human/GRCh37d5/vagrent/e75 4 > pindel_commands_BRCA.txt
# warnings give number of comparisons, but not accurate if some complete
wc -l pindel_commands_BRCA.txt
rm -f logs/pindel/BRCA.*.log
bsub -oo logs/pindel/BRCA.%I.log -sp 100 -q long -J 'tcga_pindelBRCA[1-29]%100' -P analysis-cgp -n 4 -M24000 -R"span[hosts=1] select[mem>24000] rusage[mem=24000]" '/nfs/users/nfs_k/kr2/noddy/farm_idx_exec.pl pindel_commands_BRCA.txt $LSB_JOBINDEX'
bjobs -A 948384
rm -rf /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/BRCA/*/tumour/pindel/*/logs

# caveman generations
cd /lustre/scratch116/casm/cancer_external/admin_root
source /software/CGP/pancan/final.csh.setup
mkdir -p logs/caveman
perl -e 'print qq{1\t0\t1\n};' > empty.bed
perl -ane 'print "$F[0]\t0\t$F[1]\t5\n";' < $SCRATCH112/PanCancerFinal/ref/genome.fa.fai > cn_tum.bed
perl -ane 'print "$F[0]\t0\t$F[1]\t2\n";' < $SCRATCH112/PanCancerFinal/ref/genome.fa.fai > cn_norm.bed
perl-5.16.3 ~kr2/noddy/tcga_matchList.pl /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/BRCA caveman /lustre/scratch116/casm/cgp/pancancer/reference /lustre/scratch116/casm/cgp/pancancer/reference/vagrent/e75 12 > caveman_commands_BRCA.txt
# warnings give number of comparisons, but not accurate if some complete
wc -l caveman_commands_BRCA.txt
rm -f logs/caveman/BRCA.*.log
bsub -oo logs/caveman/BRCA.%I.log -sp 100 -q long -J 'tcga_cavemanBRCA[1-426]%100' -P analysis-cgp -n 12 -M48000 -R"span[hosts=1] select[mem>48000] rusage[mem=48000]" '/nfs/users/nfs_k/kr2/noddy/farm_idx_exec.pl caveman_commands_BRCA.txt $LSB_JOBINDEX'
bjobs -A 943305
rm -rf /lustre/scratch116/casm/cancer_external/downloads/donors/TCGA/BRCA/*/tumour/caveman/*/logs
