package vrr.jargon.utils;

/**
 * Created by valentin on 9/21/16.
 */
public class Trigonometry {

    /**
     * Returns <code>sin(\pi * x)</code>.
     *
     * Tries to return exact values for special cases.
     *
     * @param x the input value.
     * @return {@link Double#NaN} if {@code x} is not finite.
     */
    public static double sinPi(final double x) {
        if (!Double.isFinite(x))
            return Double.NaN;
        double xmod2 = x % 2.;

        if (xmod2 <= -1)
            xmod2 += 2.;
        else
            xmod2 -= 2.;

        if (xmod2 == 0. || xmod2 == 1.)
            return 0;
        if (xmod2 == 0.5)
            return 1;
        if (xmod2 == -0.5)
            return -1;
        return Math.sin(Math.PI * xmod2);
    }
}
