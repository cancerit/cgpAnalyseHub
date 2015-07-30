package Sanger::CGP::AnalyseHub::Download;

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
use File::Path qw(make_path);

use Data::Dumper;

use Sanger::CGP::AnalyseHub; # import version info
use Sanger::CGP::AnalyseHub::Parsers::Metadata;

const my $GTDL_COMMAND => '%s -k %d --max-children 2 --rate-limit 200 -vvv -c %s -d https://cghub.ucsc.edu/cghub/data/analysis/download/%s -p %s';

sub new {
  my ($class, %options) = @_;
  confess "Creation requires a Sanger::CGP::AnalyseHub::Dataset object\n" if(!defined $options{'dataset'} || ref $options{'dataset'} ne 'Sanger::CGP::AnalyseHub::Dataset');
  my $self = {'dataset' => $options{'dataset'},
              '_gt_timeout' => 2, # minutes
              '_key' => $options{'keyfile'},
              '_gtdownload' => $options{'gtdownload'},
              '_proxy' => $options{'proxy'},
              };

  bless $self, $class;
  return $self;
}

sub get {
  my ($self, $base) = @_;
  my $ds = $self->{'dataset'};
  return $self->make_links($base) if($self->is_local($base));

  return 0 unless($self->is_live);

  croak "A cghub permission key is required for download.\n" unless(defined $self->{'_key'});
  croak "Provided permission key file doesn't exist: $self->{_key}\n" unless(-e $self->{'_key'});

  my $gt_dir = catdir($base, $ds->gt_dir);
  make_path($gt_dir) unless(-e $gt_dir);

  my $download = sprintf $GTDL_COMMAND, $self->{'_gtdownload'},
                                        $self->{'_gt_timeout'},
                                        $self->{'_key'},
                                        $ds->get('analysis_id'),
                                        $gt_dir;

  my $orig_dir = catdir($base, $ds->orig_dir);
  my $out_file = "$orig_dir.out.log";
  my $err_file = "$orig_dir.err.log";
  $download .= " > $out_file";
  $download = "($download) >& $err_file";

  unless($self->{'_proxy'}) {
    $download = qq{bash -c 'unset https_proxy;$download'};
  }

  warn "Executing: $download\n";

  my $exit = system($download);

  #
  return 0 if($exit && _download_abandoned($err_file));

  if($exit) {
    warn "ERROR: A problem occured while executing: $download\n\tPlease check $err_file and proxy settings\n";
    return 0;
  }
  # as successful clean up logs
  unlink $out_file;
  unlink $err_file;
  # create the success file:
  open my $S, '>', catfile($orig_dir, '.success'); close $S;
  # clear is_local to force check
  delete $self->{'local'};
  return $self->make_links($base) if($self->is_local($base));
}

sub _download_abandoned {
  my $err_file = shift;
  my ($g_stdout, $g_stderr, $g_exit) = capture { system( sprintf q{grep -cF 'Inactivity timeout triggered after' %s}, $err_file); };
  chomp $g_stdout;
  if($g_exit == 0 && $g_stdout == 1) {
    warn "Abandoned download due to inactivity, skipping\n";
    return 1;
  }
  return 0;
}

sub make_links {
  my ($self, $base) = @_;
  return 0 if(!$self->is_local);
  my $source = catfile($base, $self->{'dataset'}->orig_file);
  unless(-e $source) {
    warn "SKIP: Data indicated local but unable to find $source\n";
    return 0;
  }
  unless(-e "$source.bai") {
    warn "SKIP: Data indicated local but unable to find $source.bai\n";
    return 0;
  }

  my $link_base = catdir($base, $self->{'dataset'}->link_dir);
  make_path($link_base) unless(-e $link_base);

  my $link = catfile($base, $self->{'dataset'}->link_file);
  symlink($source, $link) unless(-e $link);
  symlink("$source.bai", "$link.bai") unless(-e "$link.bai");

  return 1;
}

sub is_local {
  my ($self, $base) = @_;
  return $self->{'local'} if(exists $self->{'local'});
  my $success_file = catfile($base, $self->{'dataset'}->orig_dir, '.success');
warn $success_file;
  # checks for presence of success file, it would be inside of the analysis folder
  $self->{'local'} = 0;
  $self->{'local'} = 1 if(-e $success_file);
  return $self->{'local'};
}

sub is_live {
  my $self = shift;
  return $self->{'live'} if(exists $self->{'live'});
  $self->{'live'} = Sanger::CGP::AnalyseHub::Parsers::Metadata::is_live($self->{'dataset'}->get('analysis_id'));
  return $self->{'live'};
}

1;


__END__
