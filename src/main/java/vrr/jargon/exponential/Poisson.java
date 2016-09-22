package vrr.jargon.exponential;

import vrr.jargon.utils.MathConstants;

/**
 * Created by valentin on 9/20/16.
 */
public class Poisson {

    public static double density(final long x, final double lambda, final boolean log) {
        return density((double) x, lambda, log);
    }

    public static double density(final long x, final double lambda) {
        return density((double)x, lambda, false);
    }

    public static double logDensity(final long x, final double lambda) {
        return density((double) x, lambda, true);
    }

    static double density(final double x, final double lambda, final boolean log) {
        if (Double.isNaN(x) || Double.isNaN(lambda))
            return Double.NaN;
        if (lambda < 0)
            return Double.NaN;

        if (!Double.isFinite(lambda))
            return log ? Double.NEGATIVE_INFINITY : 0.0;
        else if (lambda == 0) {
            if (x == 0)
                return log ? 0.0 : 1.0;
            else
                return log ? Double.NEGATIVE_INFINITY : 0.0;
        } else if (x < 0)
            return log ? Double.NEGATIVE_INFINITY : 0.0;
        else if (x <= lambda * Double.MIN_NORMAL)
            return(log ? -lambda : Math.exp(-lambda) );
        else if (lambda < x * Double.MIN_NORMAL) {
            if (Double.POSITIVE_INFINITY == x) // lambda < x = +Inf
                return log ? Double.NEGATIVE_INFINITY : 0.0;
            else
                return(log ? (-lambda + x*Math.log(lambda) -Gamma.logGamma(x+1))
                    : Math.exp(-lambda + x*Math.log(lambda) -Gamma.logGamma(x+1)) );
        }
        return (log ? -0.5*Math.log(MathConstants.TWO_PI*x)+( -Stirling.logError(x)-Binomial.deviance(x,lambda) ) :
                Math.exp( -Stirling.logError(x)-Binomial.deviance(x,lambda) )/Math.sqrt( MathConstants.TWO_PI*x));
    }

    static double density(final double x, final double lambda) {
        return density(x, lambda, false);
    }

    static double logDensity(final double x, final double lambda) {
        return density(x, lambda, true);
    }
}
