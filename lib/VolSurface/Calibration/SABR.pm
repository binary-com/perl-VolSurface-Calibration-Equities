package VolSurface::Calibration::SABR;

=head1 NAME

VolSurface::Calibration::SABR

=head1 DESCRIPTION

SABR is a calibration model that is used in volatility surface calibration.

=cut

use Moose::Role;

requires qw(symbol surface term_by_day smile_points parameterization);

use Date::Utility;
use List::MoreUtils qw(pairwise first_index indexes);
use List::Util qw(min max);
use Math::Trig qw(tanh);
use Try::Tiny;
use Format::Util::Numbers qw(roundnear);
use DataDog::DogStatsd::Helper qw(stats_gauge);

our $VERSION = '0.01';

has _max_iteration => (
    is      => 'ro',
    isa     => 'Int',
    default => 500,
);

has _tolerance_level => (
    is      => 'ro',
    isa     => 'Num',
    default => 1e-05,
);

has _scale => (
    is      => 'ro',
    isa     => 'Num',
    default => 10,
);

=head2 function_to_optimize
The function that we want to optimize.
=cut

sub function_to_optimize {
    my ($self, $params) = @_;

    # this number is returned when the optimization is heading the wrong path
    my $error = 1000;

    #if (not $params) {
    #    $self->_logger->logcroak('initial params are not present during optimization');
    #}

    my %params =
        pairwise { our $a => our $b } @{$self->calibration_param_names},
        @$params;

    my @tenors = @{$self->term_by_day};
    my @kurtosis =
        map { $self->_calculate_kurtosis($_ / 365, \%params) } @tenors;
    my @variance =
        map { $self->_calculate_variance($_ / 365, \%params) } @tenors;
    my @skew = map { $self->_calculate_skew($_ / 365, \%params) } @tenors;

    foreach my $param_name (keys %params) {
        my $param_value = $params{$param_name};
        return $error
            if ($param_name =~ /wing(R|L)/i && max(abs($param_value)) > 20.0);
        return $error
            if ($param_name =~ /growth/i && max(abs($param_value)) > 5.0);
        return $error
            if ($param_name !~ /wing(R|L)/i
            && $param_name !~ /growth/i
            && max(abs($param_value)) > 2.0);
    }

    return $error if (min(@variance) < 0);
    my @atmvols            = map { sqrt($_) } @variance;
    my $calibrated_surface = $self->get_calibrated_surface(\%params);
    my $actual_surface     = $self->surface;

    my $total = 0;
    foreach my $day (@{$self->term_by_day}) {
        my $sum       = 0;
        my $tenorvega = 0;
        my $atm_vol   = shift @atmvols;
        foreach my $point (@{$self->smile_points}) {
            my $atm_check      = log($point / 100) / $atm_vol / sqrt($day / 365);
            my $calibrated_vol = $calibrated_surface->{$day}->{$point};
            my $actual_vol     = $actual_surface->{$day}->{smile}->{$point};
            my $d1             = (log(100 / $point) + (0.5 * $actual_vol**2) * ($day / 365)) / ($actual_vol * ($day / 365)**0.5);
            my $nd1            = (1 / (2 * 3.1416)**0.5) * exp(-$d1 * $d1 / 2);
            my $vega           = $nd1 * ($day / 365)**0.5;
            $sum += $vega * abs($calibrated_vol - $actual_vol);
            $tenorvega += $vega;
        }
        $total += $sum / $tenorvega;
    }

    return $total;
}

=head2 calibration_param_names
Calibration parameter names.
It is hard-coded here because it needs to be in this sequence
=cut

has calibration_param_names => (
    is      => 'ro',
    isa     => 'ArrayRef',
    default => sub {
        [qw(
                atmvolshort
                atmvol1year
                atmvolLong
                atmWingL
                atmWingR
                skewshort
                skew1year
                skewlong
                skewwingL
                skewwingR
                kurtosisshort
                kurtosislong
                kurtosisgrowth
                )];
    },
);

=head2 default_initial_guess
Initial guess for parameters. We need to start with something.
=cut

has default_initial_guess => (
    is      => 'ro',
    isa     => 'HashRef',
    default => sub {
        {
            atmvolshort    => 0.13,
            atmvol1year    => 0.19,
            atmvolLong     => 0.12,
            atmWingL       => 1.27,
            atmWingR       => 1.56,
            skewshort      => -0.27,
            skew1year      => -0.08,
            skewlong       => -0.03,
            skewwingL      => 2.01,
            skewwingR      => 2.03,
            kurtosisshort  => 0.001,
            kurtosislong   => 0.001,
            kurtosisgrowth => 0.001
        };
    },
);

=head2 strike_ratio
=cut

has strike_ratio => (
    is      => 'ro',
    isa     => 'HashRef',
    default => sub {
        {
            lower  => 0.60,
            higher => 1.30,
        };
    },
);

=head2 compute_parameterization

    Returns a hash reference with new parameterization and calibration_error.
    These new parameters are calculated using the current parameterization values saved on couch if present,
    else it would use default parameterization.

    my $new_values = $self->compute_parameterization;
    my $new_params = $new_values->{values};
    my $new_calibration_error = $new_values->{calibration_error};

=cut

sub compute_parameterization {
    my $self = shift;

    #$self->_logger->debug('Starting calibration for [' . $self->symbol . ']');

    my $opt_using_params_from_couch;
    if ($self->parameterization) {
        my $params_from_couch = $self->_get_params_from($self->parameterization->{values});
        $opt_using_params_from_couch = try { $self->_optimize($params_from_couch) };
    }

    my $initial_params           = $self->_get_params_from($self->default_initial_guess);
    my $opt_using_initial_params = $self->_optimize($initial_params);

    my $new_values;
    if ($opt_using_params_from_couch) {
        $new_values =
            ($opt_using_params_from_couch->{calibration_error} < $opt_using_initial_params->{calibration_error})
            ? $opt_using_params_from_couch
            : $opt_using_initial_params;
    } else {
        $new_values = $opt_using_initial_params;
    }

    if ($self->price_with_parameterized_surface) {
        stats_gauge('index_vol_calibration_error ', $new_values->{calibration_error}, {tags => ['tag:' . $self->symbol]});
    }
    # we expect the params to be passed in the correct order
    my %calib_params =
        pairwise { $a => $b } @{$self->calibration_param_names},
        @{$new_values->{params}};
    my %computed_params = (
        values            => \%calib_params,
        calibration_error => $new_values->{calibration_error},
        date              => Date::Utility->new->datetime_iso8601,
    );

    return \%computed_params;
}

=head2 get_calibrated_surface

    # compute the calibrated surface with the parameters being passed.
    my $calibrated = $calibrator->get_calibrated_surface($parameters);

=cut

sub get_calibrated_surface {
    my $self   = shift;
    my $params = shift;
    my $surface;
    foreach my $tenor (@{$self->term_by_day}) {
        foreach my $point (@{$self->smile_points}) {
            $surface->{$tenor}->{$point} = $self->calculate_calibrated_vol($tenor, $point, $params);
        }
    }

    return $surface;
}

=head2 calculate_calibrated_vol

The method that calculates calibrated vol using given parameters

=cut

sub calculate_calibrated_vol {
    my ($self, $tenor, $point, $params) = @_;

    return unless defined $point;
    if (not $params) {
        $params = $self->parameterization->{values}
            || $self->default_initial_guess;
    }

    my $tiy      = $tenor / 365;
    my $kurtosis = $self->_calculate_kurtosis($tiy, $params);
    my $variance = $self->_calculate_variance($tiy, $params);
    my $skew     = $self->_calculate_skew($tiy, $params);

    my $atm_vol  = sqrt($variance);
    my $tmp_corr = 6 * (0.25 + $kurtosis * $atm_vol / 4 / $skew**2);
    my $corr     = ($tmp_corr < 0) ? -0.99 : (-0.99 * min(1, 1 / sqrt($tmp_corr)));
    my $volofvol = 2 * $skew / $corr;

    my $x     = log($point / 100) / $atm_vol / sqrt($tiy);
    my $x_min = log($self->strike_ratio->{lower}) / $atm_vol / sqrt($tiy);
    my $x_max = log($self->strike_ratio->{higher}) / $atm_vol / sqrt($tiy);

    my $vol;
    if (abs($x) < 0.00000001) {
        $vol = $atm_vol;
    } else {
        $x =
            ($x > 0)
            ? $x_max * tanh($x / $x_max)
            : $x_min * tanh($x / $x_min);
        my $z = (-1 * $volofvol * $x);
        return $atm_vol if (abs($z) < 0.00000001);
        my $d = log((sqrt(1 - 2 * $corr * $z + $z**2) - $corr + $z) / (1 - $corr));
        $vol = $atm_vol * $z / $d;
    }
    return $vol;
}

sub _get_params_from {
    my ($self, $param_hash) = @_;

    my @guess =
        map { roundnear(0.0001, $param_hash->{$_}) } @{$self->calibration_param_names};

    return \@guess;
}

sub _calculate_skew {
    my ($self, $tiy, $params) = @_;

    my $sig_small  = $params->{skewshort};
    my $sig_1y     = $params->{skew1year};
    my $sig_inf    = $params->{skewlong};
    my $t_small    = 1.0;
    my $atm_alpha1 = $params->{skewwingL};
    my $atm_alpha2 = $params->{skewwingR};

    my $skew =
        $sig_inf +
        ($sig_1y - $sig_inf) *
        ((exp(-$atm_alpha2 * $tiy) - exp(-$atm_alpha1 * $tiy)) / (exp(-$atm_alpha2 * $t_small) - exp(-$atm_alpha1 * $t_small))) +
        ($sig_small - $sig_inf) *
        ((exp(-$atm_alpha2 * $t_small - $atm_alpha1 * $tiy) - exp(-$atm_alpha1 * $t_small - $atm_alpha2 * $tiy)) /
            (exp(-$atm_alpha2 * $t_small) - exp(-$atm_alpha1 * $t_small)));

    return $skew;
}

sub _calculate_kurtosis {
    my ($self, $tiy, $params) = @_;

    my $kurt_small = $params->{kurtosisshort};
    my $kurt_alpha = $params->{kurtosisgrowth};
    my $kurt_inf   = $params->{kurtosislong};

    return $kurt_small + ($kurt_inf - $kurt_small) * (1 - exp(-$kurt_alpha * $tiy));
}

sub _calculate_variance {
    my ($self, $tiy, $params) = @_;

    my $sig_small  = $params->{atmvolshort}**2;
    my $sig_1y     = $params->{atmvol1year}**2;
    my $sig_inf    = $params->{atmvolLong}**2;
    my $t_small    = 1.0;
    my $atm_alpha1 = $params->{atmWingL};
    my $atm_alpha2 = $params->{atmWingR};

    my $atm_vol =
        $sig_inf +
        ($sig_1y - $sig_inf) *
        ((exp(-$atm_alpha2 * $tiy) - exp(-$atm_alpha1 * $tiy)) / (exp(-$atm_alpha2 * $t_small) - exp(-$atm_alpha1 * $t_small))) +
        ($sig_small - $sig_inf) *
        ((exp(-$atm_alpha2 * $t_small - $atm_alpha1 * $tiy) - exp(-$atm_alpha1 * $t_small - $atm_alpha2 * $tiy)) /
            (exp(-$atm_alpha2 * $t_small) - exp(-$atm_alpha1 * $t_small)));

    return $atm_vol;
}

# Algorithm change - now based on centroid calculations
# A function that optimizes a set of parameters against a function.
# This optimization method is based on Amoeba optimization

sub _optimize {
    my ($self, $params) = @_;

    my $intol         = $self->_tolerance_level;
    my $num_of_var    = scalar(@$params);
    my $num_of_points = $num_of_var;

    my $highest_value;

    my $lambda1 = ((($num_of_points + 1)**0.5) - 1 + $num_of_points) / ($num_of_points * (2**0.5)) * $self->_scale;
    my $lambda2 = ((($num_of_points + 1)**0.5) - 1) / ($num_of_points * (2**0.5)) * $self->_scale;

    my $i;
    my $j;

    #initialize the simplex with the initial guess
    my @simplex = ($params);

    my $step = 0;
    for ($i = 0; $i < $num_of_var; $i++) {
        $step = 0.1 * abs($params->[$i])
            if (0.1 * abs($params->[$i]) > $step);
    }

    for (my $i = 1; $i <= $num_of_points; $i++) {
        my @new = map { ($i - 1 == $_) ? $simplex[0][$_] + $step : 0 + $simplex[0][$_] } (0 .. $num_of_var - 1);
        push @simplex, \@new;
    }

    my @function_eval = map { $self->function_to_optimize($_) } @simplex;
    my $starting_value = $function_eval[0];

    my @current_reflection;
    my @current_expansion;
    my @current_contraction;

    my $expansion_function_eval;
    my $reflected_function_eval;
    my $contraction_function_eval;
    my $lowest_point_index = 0;

    my $calcvert = 1;
    my $counter  = 0;

    while (1) {
        if ($calcvert == 1) {
            $calcvert = 0;
            for ($j = 0; $j <= $num_of_points; $j++) {
                if ($j != $lowest_point_index) {
                    for ($i = 0; $i < $num_of_points; $i++) {
                        $current_reflection[$i] = $simplex[$j][$i];
                    }
                    $function_eval[$j] = $self->function_to_optimize(\@current_reflection);
                }
            }

        }
        my @sorted_eval = sort { $a <=> $b } @function_eval;
        $lowest_point_index = first_index { $_ == $sorted_eval[0] } @function_eval;
        my $highest_point_index = first_index { $_ == $sorted_eval[-1] } @function_eval;
        my @indexes = indexes { $_ == $sorted_eval[-2] } @function_eval;

        my @simplex_centroid = _calculate_simplex_centroid(\@simplex, $num_of_points, $num_of_var, $highest_point_index);

        $highest_value = $function_eval[$highest_point_index];

        my $convtol = $intol * ($starting_value + $intol);

        last
            if ($function_eval[$lowest_point_index] < $intol
            or ($function_eval[$highest_point_index] < ($function_eval[$lowest_point_index] + $convtol))
            or $counter > $self->_max_iteration);

        #imple reflection of the highest point
        @current_reflection =
            map { $simplex_centroid[$_] + 1 * ($simplex_centroid[$_] - $simplex[$highest_point_index][$_]) } (0 .. $num_of_var - 1);
        my $reflected_function_eval = $self->function_to_optimize(\@current_reflection);
        if ($reflected_function_eval < $function_eval[$lowest_point_index]) {

#Do simple expansion or in other words look up a little further in this direction
            @current_expansion = map { $simplex_centroid[$_] + 2 * ($current_reflection[$_] - $simplex_centroid[$_]) } (0 .. $num_of_var - 1);

            @simplex_centroid =
                map { $current_reflection[$_] } (0 .. $num_of_var - 1);
            $expansion_function_eval = $self->function_to_optimize(\@current_expansion);
            if ($expansion_function_eval < $reflected_function_eval) {

                #replace highest point with expansion point
                #assing function value to highest point
                @{$simplex[$highest_point_index]} =
                    map { $current_expansion[$_] } (0 .. $num_of_var - 1);
                $function_eval[$highest_point_index] = $expansion_function_eval;
            } else {

                #replace highest point with reflected point
                #assign highest value to reflection point
                @{$simplex[$highest_point_index]} =
                    map { $simplex_centroid[$_] } (0 .. $num_of_var - 1);
                $function_eval[$highest_point_index] = $reflected_function_eval;
            }

        } else {
            if ($reflected_function_eval < $highest_value) {

                #replace the simplex highest point with reflected point;
                @{$simplex[$highest_point_index]} =
                    map { $current_reflection[$_] } (0 .. $num_of_var - 1);
                $function_eval[$highest_point_index] = $reflected_function_eval;
            }

            @current_contraction =
                map { $simplex_centroid[$_] + 0.5 * ($simplex[$highest_point_index][$_] - $simplex_centroid[$_]) } (0 .. $num_of_var - 1);
            $contraction_function_eval = $self->function_to_optimize(\@current_contraction);
            if ($contraction_function_eval < $function_eval[$highest_point_index]) {
                @{$simplex[$highest_point_index]} =
                    map { $current_contraction[$_] } (0 .. $num_of_var - 1);
                $function_eval[$highest_point_index] = $contraction_function_eval;
            } elsif ($reflected_function_eval >= $highest_value) {
                $calcvert = 1;
                for ($i = 0; $i <= $num_of_points; $i++) {
                    if ($i != $lowest_point_index) {
                        @{$simplex[$i]} =
                            map { $simplex[$lowest_point_index][$_] + ($simplex[$i][$_] - $simplex[$lowest_point_index][$_]) / 2.0 }
                            (0 .. $num_of_var - 1);
                    }
                    $function_eval[$highest_point_index] = $self->function_to_optimize($simplex[$highest_point_index]);
                }
            }
        }
        $counter++;
    }

    my @sorted_eval = sort { $a <=> $b } @function_eval;
    $lowest_point_index = first_index { $_ == $sorted_eval[0] } @function_eval;

    my $new_params            = $simplex[$lowest_point_index];
    my $new_calibration_error = $function_eval[$lowest_point_index] * 1000;

    return {
        params            => $new_params,
        calibration_error => $new_calibration_error,
    };
}

sub _calculate_simplex_centroid {
    my ($simplex, $num_of_points, $num_of_var, $highest_point_index) = @_;

    my @simplex_sum;

    for (my $j = 0; $j < $num_of_var; $j++) {
        my $sum = 0;
        for (my $i = 0; $i <= $num_of_points; $i++) {
            $sum += $simplex->[$i]->[$j] if ($i != $highest_point_index);
        }

        #Centroied instead of sum for easy of calculation
        $simplex_sum[$j] = $sum / $num_of_var;
    }
    return @simplex_sum;
}

1;
