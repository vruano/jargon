package vrr.jargon.exponential;

import vrr.jargon.utils.DoubleConstants;
import vrr.jargon.utils.MathConstants;

/**
 * Created by valentin on 9/20/16.
 */
public class Binomial {

    static double deviance(final double x, final double np) {
        double ej, s, s1, v;
        int j;

        if(!Double.isFinite(x) || !Double.isFinite(np) || np == 0.0)
            return Double.NaN;

        if (Math.abs(x-np) < 0.1*(x+np)) {
            v = (x-np)/(x+np);  // might underflow to 0
            s = (x-np)*v;/* s using v -- change by MM */
            if(Math.abs(s) < Double.MIN_NORMAL) return s;
            ej = 2*x*v;
            v = v*v;
            for (j = 1; j < 1000; j++) { /* Taylor series; 1000: no infinite loop
					as |v| < .1,  v^2000 is "zero" */
                ej *= v;// = v^(2j+1)
                s1 = s+ej/((j<<1)+1);
                if (s1 == s) /* last term was effectively 0 */
                    return s1 ;
                s = s1;
            }
        }
    /* else:  | x - np |  is not too small */
        return(x*Math.log(x/np)+np-x);
    }

    private static double density(final double x, final double n, final double p, final double q, final boolean log)
    {
        double lf, lc;

        if (p == 0)
            return (x == 0) ? (log ? 0.0 : 1.0) : (log ? Double.NEGATIVE_INFINITY : 0.0);
        else if (q == 0)
            return (x == n) ? (log ? 0.0 : 1.0) : (log ? Double.NEGATIVE_INFINITY : 0.0);
        else if (x == 0) {
            if(n == 0)
                return log ? 0.0 : 1.0;
            lc = (p < 0.1) ? -deviance(n,n*q) - n*p : n* Math.log(q);
            return log ? lc : Math.exp(lc);
        } else if (x == n) {
            lc = (q < 0.1) ? -deviance(n,n*p) - n*q : n* Math.log(p);
            return log ? lc : Math.exp(lc);
        } else if (x < 0 || x > n)
            return log ? Double.NEGATIVE_INFINITY : 0.0;
        else {
            lc = Stirling.logError(n) - Stirling.logError(x) - Stirling.logError(n-x) - deviance(x,n*p) - deviance(n-x,n*q);
            lf = MathConstants.LN_2PI + Math.log(x) + Math.log1p(- x/n);
            return log ? lc - 0.5*lf : Math.exp(lc - 0.5 *lf);
        }
    }

    public static double density(final int x, final int n, final double p, final boolean log) {
        if (Double.isNaN(p))
            return Double.NaN;
        else if (p < 0 || p > 1 || n < 0)
            return Double.NaN;
        else if (x < 0 || !Double.isFinite(x))
            return log ? Double.NEGATIVE_INFINITY : 0.0;
        else
            return density(x, n, p, 1 - p, log);
    }

    public static double density(final int x, final int n, final double p) {
        return density(x, n, p, false);
    }

    public static double logDensity(final int x, final int n, final double p) {
        return density(x, n, p, true);
    }

    public static int inverseCDF(double p, final int n, double pr, boolean log_p) {
        double q, mu, sigma, gamma, z;
        int y;

        if (Double.isNaN(p) || Double.isNaN(n) || Double.isNaN(pr))
            return -1;
        else if (!Double.isFinite(n) || !Double.isFinite(pr))
            return -1;
        else if (!Double.isFinite(p) && !log_p)
            return -1;
        else if (pr < 0 || pr > 1 || n < 0)
            return -1;
        else if (log_p) {
            if (p == 0.0)
                return n;
            else if (p == Double.NEGATIVE_INFINITY)
                return 0;
            else if (p > 0.0)
                return -1;
        } else {
            if (p == 1.0)
                return n;
            else if (p == 0.0)
                return 0;
            else if (p < 0 || p > 1)
                return -1;
        }

        if (pr == 1.0)
            return n;
        else if (pr == 0 || n == 0)
            return 0;
        else {
            q = 1 - pr;
            mu = n * pr;
            sigma = Math.sqrt(n * pr * q);
            gamma = (q - pr) / sigma;
            double unloggedP = log_p ? Math.exp(p) : p;
            if (unloggedP + 1.01 * DoubleConstants.DBL_EPSILON >= 1.0)
                return n;
            z = Gaussian.inverseCDF(unloggedP, 0., 1.);
            y = Math.min((int) Math.floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5), n);
            z = CDF(y, n, pr, false);
            unloggedP *= 1 - 64 * DoubleConstants.DBL_EPSILON;
            if (n < 1e5)
                return invertedCDFSearch(y, z, unloggedP, n, pr, 1);
            else {
                int oldIncr;
                int incr = (int) Math.floor(n * 0.001);
                do {
                    oldIncr = incr;
                    y = invertedCDFSearch(y, z, unloggedP, n, pr, incr);
                } while (oldIncr > 1 && incr > n * 1e-15);
                return y;
            }
        }
    }

    static int invertedCDFSearch(int y, double z, double p, final int n, double pr, final int incr)
    {
        if(z >= p) {
            while (true) {
                if (y == 0 || (z = CDF(y - incr, n, pr, false)) < p)
                    return y;
                y = Math.max(0, y - incr);
            }
        } else {
            while (true) {
                y = Math.min(y + incr, n);
                if (y == n || (z = CDF(y, n, pr, false)) >= p)
                    return y;
            }
        }
    }

    public static double CDF(final int x, final int n, final double p, final boolean log) {
        if (!Double.isFinite(p))
            return Double.NaN;
        else if(n < 0 || p < 0 || p > 1)
            return Double.NaN;
        else if (x < 0)
            return log ? Double.NEGATIVE_INFINITY : 0.0;
        else if (n <= x)
            return log ? 0.0 : 1.0;
        else {
            return Beta.CDF(p, x + 1, n - x, false, log);
        }
    }

    public static double CDF(final int x, final int n, final double pr) {
        return CDF(x, n, pr, false);
    }

    public static double logCDF(final int x, final int n, final double pr) {
        return CDF(x, n, pr, true);
    }

    public static int inverseCDF(final double q, final int n, final double pr) {
        return inverseCDF(q, n, pr, false);
    }

    public static int logInverseCDF(final double q, final int n, final double pr) {
        return inverseCDF(q, n, pr, true);
    }
}
