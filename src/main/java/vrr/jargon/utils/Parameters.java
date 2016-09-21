package vrr.jargon.utils;

/**
 * Created by valentin on 9/20/16.
 */
public class Parameters {

    public static double requirePositiveFinite(final double v) {
        if (!Double.isFinite(v) || v <= 0)
            throw new IllegalArgumentException();
        return v;
    }

    public static double requireBetween(final double x, final double min, final double max) {
        if (! (x >= min) || ! (x <= max) )
            throw new IllegalArgumentException();
        return x;
    }
}
