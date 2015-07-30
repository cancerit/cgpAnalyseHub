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
use warnings;
use Test::More;
use Test::Fatal;
use Test::Warn;
use Const::Fast qw(const);
use FindBin qw($Bin);

my ($script) = $0 =~ m|/(.+)[.]t$|;
my $data_root = "$Bin/data/$script";

const my $MODULE => 'Sanger::CGP::AnalyseHub::Parsers::Cart';

const my $BC_TMPL => 'TCGA-KK-A7B2-%02sA-11D-A329-08';
const my @TYPE_T => (1,2,3,4,5,6,7,8,9,40); # tumour
const my @TYPE_N => (10,11,12,13,14,15,16,19,18,19,20); # normal, 10-19 are reserved, only 10-14 in use when written
const my @TYPE_C => (50); # cell (line)
const my @TYPE_X => (60,61); # xeno (grafts)

my $obj;
subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
  $obj = new_ok( $MODULE => ["$data_root/header_good.tsv"] );
  isa_ok($obj, $MODULE);
};

subtest 'File and incorrect header checks' => sub {
  like( exception {$MODULE->new("$data_root/header_bad.tsv")},
        qr/ERROR: File .+, header line is missing the following items/,
        'Check bad header fails');
  like( exception {$MODULE->new("$data_root/wibble.tsv")},
        qr/ERROR: File .+ does not exist/,
        'Check absent file fails');
  like( exception {$MODULE->new("$data_root/empty.tsv")},
        qr/ERROR: File .+ is empty/,
        'Check empty file fails');
};

subtest 'Record checks' => sub {
  $obj = new_ok( $MODULE => ["$data_root/records_good.tsv"] );
  like( exception {$MODULE->new("$data_root/records_col_mismatch.tsv")},
        qr/ERROR: File .+ appears to be corrupt, header has [[:digit:]]+ columns, line [[:digit:]]+ has [[:digit:]]+ columns/,
        'Check column mismatch fails');
};

subtest 'Barcode checks' => sub {
  $obj = new_ok( $MODULE => ["$data_root/header_good.tsv"] );
  for(@TYPE_T) {
    my $bc = sprintf $BC_TMPL, $_;
    is($obj->type_from_barcode($bc), 'tumour', "Check tumour bc conversion: $_");
  }
  for(@TYPE_N) {
    my $bc = sprintf $BC_TMPL, $_;
    is($obj->type_from_barcode($bc), 'normal', "Check normal bc conversion: $_");
  }
  for(@TYPE_C) {
    my $bc = sprintf $BC_TMPL, $_;
    is($obj->type_from_barcode($bc), 'cell', "Check cell(line) bc conversion: $_");
  }
  for(@TYPE_X) {
    my $bc = sprintf $BC_TMPL, $_;
    is($obj->type_from_barcode($bc), 'xeno', "Check xeno(graft) bc conversion: $_");
  }

  like( exception{ $obj->type_from_barcode('XXXX-XX-XXXX-123BA-11D-XXXX-XX') },
        qr/ERROR: barcode '.+' is not in TCGA format/,
        'Check bad barcode format caught');

  my $other_warn_bc = 'XXXX-XX-XXXX-99A-11D-XXXX-XX';
  {
    local $SIG{__WARN__} = sub {}; # temporarily hide warnings
    is($obj->type_from_barcode($other_warn_bc), 'other', "Check other bc conversion when new values encountered");
  }

  warning_like {$obj->type_from_barcode($other_warn_bc);} qr/WARNING: Unknown sample type decode (.+) in barcode: .+, setting to 'other'/,
              'Check warning for unknown type';
};

done_testing();
