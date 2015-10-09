use strict;
use warnings;

use Test::Exception;
use Test::More tests => 1;
use Test::NoWarnings;
use VolSurface::Calibration::SABR;

my $sabr = VolSurface::Calibration::SABR->new(
    surface => {}, 
    symbol => '', 
    term_by_day => [], 
    parameterization => {}, 
    smile_points => []);

is($base_date->days_between($base_date),     0,  'base to base days_between');
