package vrr.jargon.exceptions;

/**
 * Created by valentin on 9/20/16.
 */
public class PrecissionException extends RuntimeException {
    private final double value;
    public PrecissionException(final double x, final String message) {
        super(message);
        value = x;
    }

    public double getValue() {
        return value;
    }
}
