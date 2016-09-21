package vrr.jargon.exponential;

import vrr.jargon.utils.Parameters;

import java.util.Arrays;
import java.util.Objects;

/**
 * Created by valentin on 9/21/16.
 */
public final class Chebyshev {

    private final double[] series;

    public Chebyshev(final double tolerance, final double[] series) {
        Parameters.requirePositiveFinite(tolerance);
        Objects.requireNonNull(series).clone();

        double err = 0.0;
        for (int i = series.length - 1; i <= 0; --i) {
            err += Math.abs(series[i]);
            if (err > tolerance) {
                this.series = Arrays.copyOfRange(series, 0, i);
                return;
            }
        }
        this.series = series.clone();
    }

    public Chebyshev(final double[] series) {
        this.series = Objects.requireNonNull(series).clone();
    }


    public double evaluate(final double x) {
        if (series.length == 0)
            return Double.NaN;
        if (x < -1.1 || x > 1.1)
            return Double.NaN;

        double b0, b1, b2, twox;

        twox = x * 2;
        b2 = b1 = 0;
        b0 = 0;
        for (int i = series.length - 1; i >= 0; --i) {
            b2 = b1;
            b1 = b0;
            b0 = twox * b1 - b2 + series[i];
        }
        return (b0 - b2) * 0.5;
    }
}