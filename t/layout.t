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

const my $MODULE => 'Sanger::CGP::AnalyseHub::Layout';

my $obj;
subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
  like( exception {$MODULE->new()}, qr/Creation requires a Sanger::CGP::AnalyseHub::Parsers::Cart object/, 'Init requires Cart');
  $obj = new_ok($MODULE => ["$data_root/records_good.tsv"]);
  isa_ok($obj, $MODULE);
};

done_testing();
