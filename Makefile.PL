use 5.014;
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME             => 'VolSurface::Calibration::SABR',
    AUTHOR           => q{Binary.com <support@binary.com>},
    VERSION_FROM     => 'lib/VolSurface/Calibration/SABR.pm',
    ABSTRACT         => 'lib/VolSurface/Calibration/SABR.pm',
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
        #'ABC'              => 1.6,
        #'Foo::Bar::Module' => 5.0401,
    },
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean => { FILES => 'VolSurface-Calibration-SABR-*' },
);
