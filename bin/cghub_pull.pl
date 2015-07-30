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
  my $cart = Sanger::CGP::AnalyseHub::Parsers::Cart->new($options->{'summary'}, $options->{'assembly'});

  my $donors = $cart->participants;
  my $gb = $cart->total_size_gb;
  my $gb_p_donor = $gb / $donors;

  warn sprintf "%d donors\n", $donors;
  warn sprintf "%d tumours\n", $cart->tumours;
  warn sprintf "%.2f GB per donor\n", $gb_p_donor;
  warn sprintf "%.2f TB data total\n", $gb / 1024;

  exit 0 if($options->{'info'});

  my $layout = Sanger::CGP::AnalyseHub::Layout->new($cart); # as different cart types may be necessary

  while (my $set = $layout->next_part_set) {
    for my $ds(@{$set}) {
      my $dl = Sanger::CGP::AnalyseHub::Download->new( dataset => $ds,
                                                       keyfile => $options->{'key'},
                                                       gtdownload => $options->{'gtdownload'},
                                                       proxy => $options->{'proxyon'});
      $dl->get($options->{'outdir'});
    }
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
		'a|assembly=s' => \$opts{'assembly'},
		's|summary=s' => \$opts{'summary'},
		'g|gtbin=s' => \$opts{'gtbin'},
		'k|key=s' => \$opts{'key'},
		'u|url=s' => \$opts{'url'},
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

    --assembly  (-a)  Which assembly data to use.

    --outdir    (-o)  Where to save downloads.

  Other options:

    --key       (-k)  Required for actual downloads.

    --gtbin     (-g)  Specify gtdownload bin directory (default environment)

    --proxyon   (-p)  Don't clear https_proxy variable for download

    --symlinks  (-y)  Rebuild symlinks only

    --threads   (-t)  Number of parallel GNOS retrievals.

    --url       (-u)  The base URL to retrieve jsonl file from
                        [http://pancancer.info/gnos_metadata/latest/]

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
