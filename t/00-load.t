#!perl -T
use 5.006;
use strict;
use warnings;
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'VolSurface::Calibration::SABR' ) || print "Bail out!\n";
}

diag( "Testing VolSurface::Calibration::SABR $VolSurface::Calibration::SABR::VERSION, Perl $], $^X" );
