package Sanger::CGP::AnalyseHub::Parsers::Metadata;

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

use LWP::UserAgent;
use XML::Simple;

use Data::Dumper;

use Sanger::CGP::AnalyseHub; # import version info

const my $GT_META_BASE => 'https://cghub.ucsc.edu/cghub/metadata/';

sub get {
  my ($url) = @_;
  my $ua = LWP::UserAgent->new;
  $ua->timeout(10);
  $ua->env_proxy;
  return $ua->get($url);
}

sub is_live {
  my $analysis_id = shift;
  my $live = 0;
  # here we want to confirm that the dataset is still live at CGHub
  my $ua = LWP::UserAgent->new();
  $ua->timeout(10);
  $ua->env_proxy;

  my $response = get($GT_META_BASE.'analysisDetail/'.$analysis_id);
  if($response->is_success) {
    my $ref = XMLin($response->decoded_content);
    $live = 1 if($ref->{'Result'}->{'state'} eq 'live');
  }
  return $live;
}

sub platform_library {
  my ($analysis) = @_;
  my @rgs;
  my $response = get($GT_META_BASE.'analysisFull/'.$analysis->{'analysis_id'});
  if($response->is_success) {
    my $ref = XMLin($response->decoded_content, ForceArray => ['RUN']);
    for my $run(@{$ref->{'Result'}->{'run_xml'}->{'RUN_SET'}->{'RUN'}}) {
      push @rgs, {  PU => $run->{'alias'},
                    LB => $run->{'EXPERIMENT_REF'}->{'refname'},
                    assembly => $analysis->{'assembly'},
                    analysis_id => $analysis->{'analysis_id'},
                 };
    }
  }
  return \@rgs;
}



1;

