#!/usr/bin/perl

use strict;

my ($file, $idx) = @ARGV;

my $line = `head -n $idx $file | tail -n 1`;
chomp $line;
print "Executing: $line\n";
exec $line;
