package Sanger::CGP::AnalyseHub::Layout;

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

use File::Spec::Functions qw(catdir catfile);

use Data::Dumper;

use Sanger::CGP::AnalyseHub; # import version info
use Sanger::CGP::AnalyseHub::Dataset;

sub new {
  my ($class, $cart) = @_;
  confess "Creation requires a Sanger::CGP::AnalyseHub::Parsers::Cart object\n" unless(defined $cart);
  my $self = {'cart' => $cart};
  bless $self, $class;
  return $self;
}

sub next_part_set {
  my $self = shift;
  my $part = $self->{'cart'}->next_participant;
  return undef unless(defined $part);
  my @set;
  my @bcs = keys %{$part};
  for my $bc(@bcs) {
    for my $item($part->{$bc}) {
      my $ds = Sanger::CGP::AnalyseHub::Dataset->new($item);

      my $structured = catdir($ds->{'study'}, $ds->{'disease'}, $ds->{'participant_id'});

      my $gt_dir = catdir('original', $structured);
      my $orig_dir = catdir($gt_dir, $ds->{'analysis_id'});
      my $orig_file = catfile($orig_dir, $ds->{'filename'});
      my $link_dir = catdir('donors', $structured, $ds->{'simple_type'});
      my $link_file = catfile($link_dir, $ds->{'filename'});

      $ds->gt_dir($gt_dir);       # this value passed to GT download as output location (create in advance)
      $ds->orig_dir($orig_dir);       # this value passed to GT download as output location (create in advance)
      $ds->orig_file($orig_file); # this is the full path to the BAM once GT download is complete
      $ds->link_dir($link_dir);   # this is the folder the bam will be linked into (create in advance)
      $ds->link_file($link_file); # prebuilt name for link.bam (need to link *.bai also)

      push @set, $ds;
    }
  }
  return \@set;
}



1;

__END__

=Name Sanger::CGP::AnalyseHub::Layout

Controls the arrangement of data on the filesystem

