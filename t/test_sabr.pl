use strict;
use warnings;
use Data::Dumper;

use VolSurface::Calibration::SABR;

my $sabr = VolSurface::Calibration::SABR->new(
    surface => {}, 
    term_by_day => [7, 31, 61, 94, 185, 276, 367],
    smile_points => [80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120],
);

my $parameterization = $sabr->compute_parameterization;

