package vrr.jargon.utils;

/**
 * Created by valentin on 9/21/16.
 */
public class Doubles {

    /**
     * Returns the integer part of a double.
     * @param v the value to truncate.
     * @return the truncated value.
     */
    public static double truncate(final double v) {
        if (v < 0) {
            return Math.ceil(v);
        } else {
            return Math.floor(v);
        }
    }
}
