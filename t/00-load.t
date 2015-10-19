#!perl -T
use 5.006;
use strict;
use warnings;
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'VolSurface::Calibration::Equities' ) || print "Bail out!\n";
}

diag( "Testing VolSurface::Calibration::Equities $VolSurface::Calibration::Equities::VERSION, Perl $], $^X" );
