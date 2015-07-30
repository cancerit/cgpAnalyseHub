package Sanger::CGP::AnalyseHub::Dataset;

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

const my @REQ_HEADERS => qw (study analysis_id barcode disease library_type platform participant_id state barcode filename simple_type);

sub new {
  my ($class, $hash) = @_;
  my $self = {};
  if(defined $hash) {
    _validate($hash);
    $self = $hash; # warning, changes input into an object
  }
  bless $self, $class;
  return $self;
}

sub _validate {
  my $unblessed = shift;
  for my $req(@REQ_HEADERS) {
    confess "Input doesn't have required field of '$req'\n" unless(exists $unblessed->{$req});
  }
  return 1;
}

sub get {
  my ($self, $query) = @_;
  return $self->{$query};
}

sub gt_dir {
  my ($self, $val) = @_;
  $self->{'_gt_dir'} = $val if(defined $val);
  return $self->{'_gt_dir'};
}

sub orig_dir {
  my ($self, $val) = @_;
  $self->{'_orig_dir'} = $val if(defined $val);
  return $self->{'_orig_dir'};
}

sub orig_file {
  my ($self, $val) = @_;
  $self->{'_orig_file'} = $val if(defined $val);
  return $self->{'_orig_file'};
}

sub link_dir {
  my ($self, $val) = @_;
  $self->{'_link_dir'} = $val if(defined $val);
  return $self->{'_link_dir'};
}

sub link_file {
  my ($self, $val) = @_;
  $self->{'_link_file'} = $val if(defined $val);
  return $self->{'_link_file'};
}



1;


__END__
