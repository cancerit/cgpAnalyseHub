#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Carp qw(croak);
use Getopt::Long;
use Pod::Usage;
use Const::Fast qw(const);
use File::Which qw(which);

use threads;

use Data::Dumper;

use PCAP::Cli;
use Sanger::CGP::AnalyseHub;
use Sanger::CGP::AnalyseHub::Download;
use Sanger::CGP::AnalyseHub::Layout;
use Sanger::CGP::AnalyseHub::Parsers::Cart;

const my $CGHUB_BASE => 'https://cghub.ucsc.edu/cghub/metadata/analysisFull';

{
  my $options = option_builder();
  gtdl_setup($options);
  my $cart = Sanger::CGP::AnalyseHub::Parsers::Cart->new( $options->{'summary'},
                                                          $options->{'normtiss'},
                                                          $options->{'assembly'},
                                                          $options->{'bad'});


#die Dumper($cart);

  my $donors = $cart->participants;
  my $gb = $cart->total_size_gb;
  my $gb_p_donor = 0;
  $gb_p_donor = $gb / $donors if($donors);

  warn sprintf "%d donors\n", $donors;
  warn sprintf "%d tumours\n", $cart->tumours;
  warn sprintf "%.2f GB per donor\n", $gb_p_donor;
  warn sprintf "%.2f TB data total\n", $gb / 1024;

  my $layout = Sanger::CGP::AnalyseHub::Layout->new($cart); # as different cart types may be necessary

  if($options->{'info'}) {
    my $total_local = 0;
    while (my $set = $layout->next_part_set) {
      $total_local += local_count($set, $options);
    }
    warn sprintf "%d files already local\n", $total_local;
    exit 0;
  }

  my $thread_count = $options->{'threads'} || 1;

  if($thread_count == 1) {
    while (my $set = $layout->next_part_set) {
      download($set, $options);
    }
  }
  else {
    while (my $set = $layout->next_part_set) {
      if(threads->list(threads::all) < $thread_count) {
        threads->create(\&download, $set, $options);
        next if(threads->list(threads::all) < $thread_count);
      }
      sleep 10 while(threads->list(threads::joinable) == 0);
      for my $thr(threads->list(threads::joinable)) {
        $thr->join;
        if(my $err = $thr->error) { die "Thread error: $err\n"; }
      }
    }
    sleep 10 while(threads->list(threads::running) > 0);
    for my $thr(threads->list(threads::joinable)) {
      $thr->join;
      if(my $err = $thr->error) { die "Thread error: $err\n"; }
    }
  }
}

sub local_count {
  my ($set, $options) = @_;
  my $local = 0;
  for my $ds(@{$set}) {
    my $dl = Sanger::CGP::AnalyseHub::Download->new( dataset => $ds,
                                                     keyfile => $options->{'key'},
                                                     gtdownload => $options->{'gtdownload'},
                                                     proxy => $options->{'proxyon'});
    $dl->debug($options->{'debug'});
    $dl->make_links($options->{'outdir'}, 1) if($options->{'symlinks'});
    $local++ if($dl->is_local($options->{'outdir'}));
  }
  return $local;
}

sub download {
  my ($set, $options) = @_;
  for my $ds(@{$set}) {
    my $dl = Sanger::CGP::AnalyseHub::Download->new( dataset => $ds,
                                                     keyfile => $options->{'key'},
                                                     gtdownload => $options->{'gtdownload'},
                                                     proxy => $options->{'proxyon'});
    $dl->get($options->{'outdir'});
  }
}

sub gtdl_setup {
  my $options = shift;
  if(exists $options->{'gtbin'} && defined $options->{'gtbin'}) {
    my $gtbin = $options->{'gtbin'};
    $options->{'gtdownload'} = $gtbin.'/gtdownload';
    my $gtshare = $gtbin;
    $gtshare =~ s|bin/?$|share/GeneTorrent|;
    $options->{'gtdownload'} .= ' -R '.$gtshare;
  }
  else {
    $options->{'gtdownload'} = which('gtdownload');
  }
  return 1;
}


sub option_builder {
	my ($factory) = @_;

	my %opts = ();

	my $result = &GetOptions (
		'h|help|?' => \$opts{'h'},
		'm|man' => \$opts{'m'},
		'i|info' => \$opts{'info'},
		'y|symlinks' => \$opts{'symlinks'},
  	'n|normtiss' => \$opts{'normtiss'},
		'b|bad=s' => \$opts{'bad'},
		'a|assembly=s' => \$opts{'assembly'},
		's|summary=s' => \$opts{'summary'},
		'g|gtbin=s' => \$opts{'gtbin'},
		'k|key=s' => \$opts{'key'},
		't|threads=i' => \$opts{'threads'},
		'o|outdir=s' => \$opts{'outdir'},
		'd|debug' => \$opts{'debug'},
		'p|proxyon' => \$opts{'proxyon'},
    'v|version' => \$opts{'version'},
	) or pod2usage(2);

	pod2usage(1) if(defined $opts{'h'});
  pod2usage(-exitval => 0, -verbose => 2) if(defined $opts{'m'});

  if(exists $opts{'version'} && defined $opts{'version'}) {
    print sprintf "VERSION: %s\n", Sanger::CGP::AnalyseHub->VERSION;
    exit 0;
  }

  $opts{'threads'} = 1 unless(defined $opts{'threads'});

  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});
  PCAP::Cli::file_for_reading('summary', $opts{'summary'});

  $opts{'url'} = $CGHUB_BASE unless(defined $opts{'url'});

	return \%opts;
}

__END__

=head1 NAME

gnos_pull.pl - retrieve/update analysis flow results on local systems.

=head1 SYNOPSIS

./gnos_pull.pl [-h] -u http://pancancer.info/gnos_metadata/latest/ -c gnos_pull.ini -o local_mirror/

  Required input:

    --summary   (-s)  Summary TSV file from CGHub cart download.

    --outdir    (-o)  Where to save downloads.

    --key       (-k)  Required for actual downloads.

  Other options:

    --assembly  (-a)  Which assembly data to use.

    --normtiss  (-n)  Presence indicates normal tissue in preference to blood.

    --bad       (-b)  File listing analysis_ids for known bad data.

    --gtbin     (-g)  Specify gtdownload bin directory (default environment)

    --proxyon   (-p)  Don't clear https_proxy variable for download

    --symlinks  (-y)  Rebuild symlinks only

    --threads   (-t)  Number of parallel GNOS retrievals.

    --info      (-i)  Just prints how many donor's will be included in pull and some stats.

    --debug     (-d)  prints extra debug information

    --help      (-h)  Brief documentation

    --man       (-m)  More verbose usage info

=head1 OPTIONS

=over 2

=item s|summary

A TSV file downloaded from CGHubs' Cart tool.

=item o|outdir

The base output directory, this will need to be very large if pulling BAM files.

=item u|url

Provide if standard CGHub URL is not functioning

=item i|info

Find out how much data you will potentially retrieve.  Can also be used to see how data is distributed.

This doesn't take into account data already downloaded.

=item d|debug

Turns on debug information which can aid in diagnosing why a particular donor is not retrieved.

=back
