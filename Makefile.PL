use 5.014;
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME             => 'VolSurface::Calibration::Equities',
    AUTHOR           => q{Binary.com <support@binary.com>},
    VERSION_FROM     => 'lib/VolSurface/Calibration/Equities.pm',
    ABSTRACT         => 'lib/VolSurface/Calibration/Equities.pm',
    LICENSE          => 'artistic_2',
    PL_FILES         => {},
    MIN_PERL_VERSION => 5.006,
    CONFIGURE_REQUIRES => {
        'ExtUtils::MakeMaker' => 0,
    },
    BUILD_REQUIRES => {
        'Test::More' => 0,
    },
    PREREQ_PM => {
        'Date::Utility' => 0,
        'Format::Util::Numbers' => 0,
        'Try::Tiny' => 0,
        'Math::Trig' => 0,
        'List::Util' => 0,
        'List::MoreUtils' => 0
    },
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean => { FILES => 'VolSurface-Calibration-Equities-*' },
);
