package Sanger::CGP::AnalyseHub::Parsers::Cart;

########## LICENCE ##########
# Copyright (c) 2015 Genome Research Ltd.
#
# Author: Keiran Raine <cgpit@sanger.ac.uk>
#
# This file is part of cgpAnalyseHub.
#
# cgpAnalyseHub is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
########## LICENCE ##########

use strict;
use autodie qw(:all);
use warnings FATAL => 'all';
use Carp qw(croak confess);
use Const::Fast qw(const);

use Data::Dumper;

use Sanger::CGP::AnalyseHub; # import version info
use Sanger::CGP::AnalyseHub::Parsers::Metadata;

# constants for required fields etc.
const my @REQ_HEADERS => qw (study assembly analysis_id barcode disease library_type platform participant_id state barcode filename files_size);
const my @OPT_HEADERS => qw ();
const my $BC_REGEX => qr/[^\-]+\-[^\-]+\-[^\-]+\-([[:digit:]]{2})[A-Z]\-([[:digit:]]{2})[A-Z](?:\-[^\-]+\-[^\-]+)?/;

sub new {
  my ($class, $summary_file, $assembly, $verbose) = @_;
  my $self = {};
  bless $self, $class;
  $self->set_verbose($verbose);
  $self->_init($summary_file, $assembly);
  return $self;
}

sub next_participant {
  my $self = shift;
  confess "next_participant can only be used after data has been parsed." unless(defined $self->{'_part_list'});
  return undef if($self->{'_iter_val'} == -1);
  return undef if($self->{'_iter_val'} == $self->{'participants'});
  return $self->{'_part_list'}->[ $self->{'_iter_val'}++ ];
}

sub set_verbose {
  my ($self, $verbose) = @_;
  $self->{'verbose'} = $verbose if(defined $verbose);
  return 1;
}

sub verbose {
  return shift->{'verbose'};
}

sub _init {
  my ($self, $summary_file, $assembly) = @_;
  croak "ERROR: File $summary_file does not exist\n" if(!-e $summary_file);
  croak "ERROR: File $summary_file is empty\n" if(-s _ == 0);
  $self->{'summary_file'} = $summary_file;
  $self->{'assembly'} = $assembly;
  $self->_parse_file();
}

sub _parse_file {
  my ($self) = @_;
  my $file = $self->{'summary_file'};
  my ($IN, $pid);
  if($file =~ m/[.]gz$/) {
    my $command = sprintf 'gunzip -c %s', $file;
    $pid = open $IN, q{-|}, $command or croak 'Could not fork: '.$!;
  }
  else {
    open $IN, '<', $file;
  }
  $self->_parse_header($IN);
  $self->_check_header();
  $self->_load($IN);
  close $IN;
}

sub _load {
  my ($self, $FH) = @_;
  my $total = 0;
  my %map = %{$self->{'_header_map'}};
  my %studies;
  while(my $line = <$FH>) {
    chomp $line;
    ### NOTE ###
    # Empty columns in the files downloaded from CGHub have a space in them
    # this is the optimal way to handle this so that it doesn't have to be done
    # on a per field basis
    $line =~ s/\t \t/\t\t/g;
    $line =~ s/^ //g; # this ensures first column isn't blank too
    # don't clear the end of the line as split won't work properly
    # in fact add it to the end if necessary
    $line .= q{ } if($line =~ m/\t$/);
    my @elements = split /\t/, $line;
    $elements[-1] = q{} if($elements[-1] eq q{ });
    die sprintf("ERROR: File %s appears to be corrupt, header has %d columns, line %d has %d columns",
                $self->{summary_file},
                $self->{'_columns'},
                $.,
                scalar @elements)
          if(scalar @elements != $self->{'_columns'});

    $total++;

    next if($elements[ $map{'filename'} ] =~ m/_gapfillers_*/);
    next if($elements[ $map{'state'} ] ne 'Live');
    if($elements[ $map{'library_type'} ] ne 'WXS') {
      warn sprintf "Skipping %s as not Exome (WXS)\n", $elements[ $map{'barcode'} ];
      next;
    }

    my $record = $self->_convert_record($self->{'assembly'}, \@elements, \%map);
    next unless(defined $record);

    push @{$studies{$record->{'study'}}
            {$record->{'disease'}}
              {$record->{'participant_id'}}
                {$record->{'barcode'}}
          }, $record;
  }

  $self->clean_set(\%studies);
  $self->{'total'} = $total;
  return 1;
}

sub clean_set {
  my ($self, $set) = @_;
  my $retained = 0;
  my $tumour_retained = 0;
  my @iterator;
  for my $study(keys %{$set}) {
    my $study_retained = 0;
    for my $disease(keys %{$set->{$study}}) {
      my $disease_retained = 0;
      for my $partid(keys %{$set->{$study}->{$disease}}) {
        my $remove = 0;
        my $participant = $set->{$study}->{$disease}->{$partid};
        # participant must have multiple types

        my %types;
        # inefficient to loop here, but can quickly discount the rest of this block if no normal or tumour
        for my $barcode(keys %{$participant}) {
          my $barcode_warned = 0;
          my @records = @{ $participant->{$barcode} };
          for my $record(@records) {
            my $type = $self->type_from_barcode($record->{'barcode'});
            $record->{'simple_type'} = $type;
            $types{ $type }++;
          }
        }

        if(!exists $types{'tumour'} || !exists $types{'normal'}) {
          delete $set->{$study}->{$disease}->{$partid};
          next;
        }


        for my $barcode(keys %{$participant}) {
          my $barcode_warned = 0;
          my @records = @{ $participant->{$barcode} };

          if(@records > 1) {
            my %assemblies;
            for my $analysis(@records) {
              push @{$assemblies{$analysis->{'assembly'}}}, $analysis;
            }
            if(exists $assemblies{'GRCh37-lite'} && (scalar @{$assemblies{'GRCh37-lite'}}) == 1) {
              $participant->{$barcode} = $assemblies{'GRCh37-lite'}->[0];
            }
            elsif(exists $assemblies{'HG19_Broad_variant'} && (scalar @{$assemblies{'HG19_Broad_variant'}}) == 1) {
              $participant->{$barcode} = $assemblies{'HG19_Broad_variant'}->[0];
            }
            else {
              my $largest;
              my @rg_sets;
              for my $analysis(@records) {
                push @{$assemblies{$analysis->{'assembly'}}}, $analysis;
                push @rg_sets, Sanger::CGP::AnalyseHub::Parsers::Metadata::platform_library($analysis);
                if(defined $largest) {
                  $largest = $analysis if($analysis->{'files_size'} > $largest->{'files_size'});
                }
                else {
                  $largest = $analysis
                }
              }

              ######
              # This cleans out ambiguous data sets
              my %query_same;
              for my $rgs(@rg_sets) {
                for my $rg(@{$rgs}) {
                  my $chk_str = $rg->{'assembly'}.':'.$rg->{'LB'}.':'.$rg->{'PU'};
                  push @{$query_same{$chk_str}}, $rg->{'analysis_id'};
                }
              }
              my $analysis_count = scalar @rg_sets;
              for my $chks(sort keys %query_same) {
                if((scalar @{$query_same{$chks}}) == $analysis_count) {
                  print "select largest for $barcode\n";
                  $participant->{$barcode} = $largest;
                }
                else {
                  unless($barcode_warned) {
                    my $error = "Duplicate data for barcode: $records[0]->{barcode}, ambiguous data for\n\tanalysis_ids:";
                    for my $r(@records) {
                      $error .= "\n\t\t".$r->{'analysis_id'};
                    }
                    $error .= "\n\tRespecitive filenames:";
                    for my $r(@records) {
                      $error .= sprintf "\n\t\t%s (%s)",$r->{'filename'},$r->{'assembly'};
                    }
                    $error .= "\n\n";

                    warn $error;
                    $barcode_warned++;

                  }
                  $remove++;
                }
              }
            }
            #
            ######

          }
          else {
            $participant->{$barcode} = $records[0];
          }
        }

        if($remove != 0) {
          delete $set->{$study}->{$disease}->{$partid};
          next;
        }
        $retained++;
        $study_retained++;
        $disease_retained++;
        push @iterator, $participant;
        $tumour_retained += (scalar (keys $participant)) - 1;
      }
      delete $set->{$study}->{$disease} if($disease_retained == 0);
    }
    delete $set->{$study} if($study_retained == 0);
  }
  $self->{'participants'} = $retained;
  $self->{'tumours'} = $tumour_retained;
  $self->{'_part_list'} = \@iterator;
  $self->{'_iter_val'} = (scalar @iterator) == 0 ? -1 : 0;
  return 1;
}

sub _largest_data {
  my ($self, $set) = @_;
  my @new_set;
  # Need to group by bc first then keep largest of each
  my %by_bc;
  for(@{$set}) {
    push @{$by_bc{$_->{'barcode'}}}, $_;
  }

  for my $bc(keys %by_bc) {
    my $bc_set = $by_bc{$bc};
    my $largest = shift @{$bc_set};
    for my $t(@{$bc_set}) {
      $largest = $t if($t->{'files_size'} > $largest->{'files_size'});
    }
    push @new_set, $largest;
  }
  return \@new_set;
}

sub total_size_gb {
  my $self = shift;
  my $raw = $self->total_size;
  return $raw / 1024 / 1024 / 1024;
}

sub participants {
  my $self = shift;
  confess "participants can only be used after data has been parsed." unless(defined $self->{'participants'});
  return $self->{'participants'};
}

sub tumours {
  my $self = shift;
  confess "tumours can only be used after data has been parsed." unless(defined $self->{'tumours'});
  return $self->{'tumours'};
}

sub total_size {
  my $self = shift;
  unless(defined $self->{'total_size'}) {
    my $set = $self->{'_part_list'};
    confess "total_size can only be used after data has been parsed." unless(defined $set);
    my $size = 0;
    for my $participant(@{$set}) {
      for my $bc(keys $participant) {
if(ref $participant->{$bc} ne 'HASH') {
  warn Dumper($participant->{$bc});
  exit;
}
        $size += $participant->{$bc}->{'files_size'};
      }
    }
    $self->{'total_size'} = $size;
  }
  return $self->{'total_size'};
}

sub _convert_record {
  my ($self, $assembly, $elements, $map) = @_;
  my $valid;

  return $valid if(defined $assembly && $assembly ne $elements->[ $map->{'assembly'} ]);

  my %by_col;
  for my $col_name(@REQ_HEADERS) {
    my $val = $elements->[ $map->{$col_name} ];
    return $valid if($val eq q{});
    $by_col{$col_name} = $val;
  }

  # we do want them all
  for my $col_name(keys %{$map}) {
    my $val = $elements->[ $map->{$col_name} ];
    $by_col{$col_name} = $val;
  }

  return $valid if( $by_col{'filename'} =~ m/_IlluminaGA-DNASeq_exome.bam$/);
  return $valid if( $by_col{'filename'} =~ m/_capture.bam$/);

#  if($by_col{'filename'} ne $by_col{'checksum'}.'.bam' && $by_col{'filename'} ne $by_col{'barcode'}.'.bam') {
#    warn $by_col{'filename'};
#  }

  $valid = \%by_col;
  return $valid;
}

sub type_from_barcode {
  my ($self, $bc) = @_;
  my ($sample_type_code) = $bc =~ m/^$BC_REGEX$/;
  croak "ERROR: barcode '$bc' is not in TCGA format\n" unless(defined $sample_type_code);
  my $type = 'other';
  $sample_type_code += 0; # force to number
  if($sample_type_code != 0) {
    if($sample_type_code < 10 || $sample_type_code == 40) {
      $type = 'tumour';
    }
    elsif($sample_type_code < 21) {
      # control considered 'normal'
      $type = 'normal';
    }
    elsif($sample_type_code == 50) {
      $type = 'cell';
    }
    elsif($sample_type_code == 60 || $sample_type_code == 61) {
      $type = 'xeno';
    }
  }
  if($type eq 'other') {
    warn "WARNING: Unknown sample type decode ($sample_type_code) in barcode: $bc, setting to 'other'\n";
  }
  return $type;
}

sub _check_header {
  my $self = shift;
  my @absent;
  for my $required(@REQ_HEADERS) {
    push @absent, $required unless(exists $self->{'_header_map'}->{$required});
  }
  @absent = sort @absent;
  if(scalar @absent != 0) {
    croak "ERROR: File $self->{summary_file}, header line is missing the following items:\n".join("\n", @absent)."\n";
  }
  return 1;
}

sub _parse_header {
  my ($self, $FH) = @_;
  my $line = <$FH>;
  chomp $line;
  my %header_map;
  my @elements = split /\t/, $line;
  for my $idx(0..(scalar @elements)-1) {
    $header_map{ $elements[$idx] } = $idx;
  }
  $self->{'_header_map'} = \%header_map;
  $self->{'_columns'} = scalar @elements;
  return 1;
}


1;

=item new

Instantiate a Cart object for parsing of CGHub cart summary.tsv file.

Inputs:

  summary_file : file path to summary.tsv
  verbose      : set 1 for additional details

Returns:

  Sanger::CGP::AnalyseHub::Parsers::Cart

=item set_verbose

Turn on additional warnings for debugging purposes.

=item type_from_barcode

Uses the TCGA barcode to determine the simple sample type.
For full details of the barcode and the granular form of sample type see:

  https://wiki.nci.nih.gov/display/TCGA/TCGA+Barcode

Inputs:

  bc : TCGA barcode

Returns

  sample_type : tumour|normal|cell|xeno

=cut
