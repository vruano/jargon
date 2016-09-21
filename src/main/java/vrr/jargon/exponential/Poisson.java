package vrr.jargon.exponential;

import vrr.jargon.utils.MathConstants;

/**
 * Created by valentin on 9/20/16.
 */
public class Poisson {

    public static double density(final double x, final double lambda, final boolean log) {
        if (Double.isNaN(lambda) || lambda < 0)
            return Double.NaN;
        if (Double.isNaN(x))
            return Double.NaN;
        if (lambda == 0)
            return x == 0 ? (log ? 0.0 : 1.0) : (log ? Double.NEGATIVE_INFINITY : 0.0);
        if (lambda == Double.POSITIVE_INFINITY)
            return log ? Double.NEGATIVE_INFINITY : 0.0;
        if (x < 0)
            return log ? Double.NEGATIVE_INFINITY : 0.0;
        if (x <= lambda * Double.MIN_NORMAL) return(log ? -lambda : Math.exp(-lambda) );
        if (lambda < x * Double.MIN_NORMAL) {
            if (Double.POSITIVE_INFINITY == x) // lambda < x = +Inf
                return log ? Double.NEGATIVE_INFINITY : 0.0;
            else
                return(log ? (-lambda + x*Math.log(lambda) -Gamma.logGamma(x+1))
                    : Math.exp(-lambda + x*Math.log(lambda) -Gamma.logGamma(x+1)) );
        }
        return (log ? -0.5*Math.log(MathConstants.TWO_PI*x)+( -Stirling.logError(x)-Binomial.deviance(x,lambda) ) :
                Math.exp( -Stirling.logError(x)-Binomial.deviance(x,lambda) )/Math.sqrt( MathConstants.TWO_PI*x));
    }
}
