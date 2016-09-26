package vrr.jargon.exponential;

import vrr.jargon.utils.DoubleConstants;
import vrr.jargon.utils.Doubles;
import vrr.jargon.utils.MathConstants;

/**
 * Created by valentin on 9/21/16.
 */
public class Gaussian {

    public static double inverseCDF(final double p, final double mu, final double sigma) {
        return inverseCDF(p, mu, sigma, false);
    }


    public static double inverseCDF(final double p, final double mu, final double sigma, final boolean log)
    {
        if (Double.isNaN(p) || Double.isNaN(mu) || Double.isNaN(sigma))
            return Double.NaN;
        if (log) {
            if (p > 0.0)
                return Double.NaN;
            else if (p == 0.0)
                return Double.POSITIVE_INFINITY;
            else if (p == Double.NEGATIVE_INFINITY)
                return Double.NEGATIVE_INFINITY;
        } else {
            if (p < 0 || p > 1)
                return Double.NaN;
            else if (p == 0.0)
                return Double.NEGATIVE_INFINITY;
            else if (p == 1.0)
                return Double.POSITIVE_INFINITY;
        }

        if (sigma < 0)
            return Double.NaN;
        else if (sigma == 0)
            return mu;

        final double p_ = log ? Math.exp(p) : p;
        final double q = p_ - 0.5;
        double val;

        if (Math.abs(q) <= .425) {/* 0.075 <= p <= 0.925 */
            final double r = .180625 - q * q;
            val =
                    q * (((((((r * 2509.0809287301226727 +
                            33430.575583588128105) * r + 67265.770927008700853) * r +
                            45921.953931549871457) * r + 13731.693765509461125) * r +
                            1971.5909503065514427) * r + 133.14166789178437745) * r +
                            3.387132872796366608)
                            / (((((((r * 5226.495278852854561 +
                            28729.085735721942674) * r + 39307.89580009271061) * r +
                            21213.794301586595867) * r + 5394.1960214247511077) * r +
                            687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
        }
        else { /* closer than 0.075 from {0,1} boundary */
            double r;
	/* r = min(p, 1-p) < 0.075 */
            if (q > 0)
                r = log ?  -Math.expm1(p) : (0.5 - (p) + 0.5);/* 1-p */
            else
                r = p_;/* = R_DT_Iv(p) ^=  p */

            r = Math.sqrt(- ((log && (q <= 0)) ? p : Math.log(r)));

            if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
                r += -1.6;
                val = (((((((r * 7.7454501427834140764e-4 +
                        .0227238449892691845833) * r + .24178072517745061177) *
                        r + 1.27045825245236838258) * r +
                        3.64784832476320460504) * r + 5.7694972214606914055) *
                        r + 4.6303378461565452959) * r +
                        1.42343711074968357734)
                        / (((((((r *
                        1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + .0151986665636164571966) * r +
                        .14810397642748007459) * r + .68976733498510000455) *
                        r + 1.6763848301838038494) * r +
                        2.05319162663775882187) * r + 1.);
            }
            else { /* very close to  0 or 1 */
                r += -5.;
                val = (((((((r * 2.01033439929228813265e-7 +
                        2.71155556874348757815e-5) * r +
                        .0012426609473880784386) * r + .026532189526576123093) *
                        r + .29656057182850489123) * r +
                        1.7848265399172913358) * r + 5.4637849111641143699) *
                        r + 6.6579046435011037772)
                        / (((((((r *
                        2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                        7.868691311456132591e-4) * r + .0148753612908506148525)
                        * r + .13692988092273580531) * r +
                        .59983220655588793769) * r + 1.);
            }

            if(q < 0.0)
                val = -val;
        /* return (q >= 0.)? r : -r ;*/
        }
        return mu + sigma * val;
    }

    public static double density(final double x, final double mu, final double sigma, final boolean log) {
        if (Double.isNaN(x) || Double.isNaN(mu) || Double.isNaN(sigma))
            return Double.NaN;
        else if (!Double.isFinite(sigma))
            return log ? Double.NEGATIVE_INFINITY : 0.0;
        else if (!Double.isFinite(x) && mu == x)
            return Double.NaN;
        else if (sigma < 0)
            return Double.NaN;
        else if (sigma == 0)
            return (x == mu) ? Double.POSITIVE_INFINITY : (log ? Double.NEGATIVE_INFINITY : 0.0);
        else {
            final double standarizedX = Math.abs((x - mu) / sigma);
            if (!Double.isFinite(standarizedX))
                return log ? Double.NEGATIVE_INFINITY : 0.0;
            else if (standarizedX >= 2 * Math.sqrt(Double.MAX_VALUE))
                return log ? Double.NEGATIVE_INFINITY : 0.0;
            else if (log)
                return -(MathConstants.LN_SQRT_2PI + 0.5 * standarizedX * standarizedX + Math.log(sigma));
            else if (standarizedX < 5)
                return INV_SQRT_2PI * Math.exp(-.5 * standarizedX * standarizedX) / sigma;
            else if (standarizedX > Math.sqrt(-2 * MathConstants.LN_2 * (Double.MIN_EXPONENT +
                    1 - DoubleConstants.DBL_MANT_DIG))) {
                return 0;
            } else {
                final double x1 = Math.round(standarizedX * 65536) / 65536;
                final double x2 = standarizedX - x1;
                return INV_SQRT_2PI / sigma *
                        (Math.exp(-0.5 * x1 * x1) * Math.exp((-0.5 * x2 - x1) * x2));
            }
        }
    }

    public static double CDF(final double x, final double mu, final double sigma, final boolean log) {
        if (Double.isNaN(x) || Double.isNaN(mu) || Double.isNaN(sigma))
            return Double.NaN;
        else if (!Double.isFinite(x) && mu == x)
            return Double.NaN;
        else if (sigma < 0)
            return Double.NaN;
        else if (sigma == 0)
            return (x < mu) ? (log ? Double.NEGATIVE_INFINITY : 0.0) : (log ? 0.0 : 1.0);
        else {
            double p = (x - mu) / sigma;
            if (!Double.isFinite(p))
                return (x < mu) ? (log ? Double.NEGATIVE_INFINITY : 0.0) : (log ? 0.0 : 1.0);
            else {
                return rawCDF(p, log)[0];
            }
        }

    }
    private static final double SIXTEN = 16;

    private static final double[] RAWCDF_A = {
        2.2352520354606839287,
                161.02823106855587881,
                1067.6894854603709582,
                18154.981253343561249,
                0.065682337918207449113
    };

    private static final double[] RAWCDF_B = {
        47.20258190468824187,
                976.09855173777669322,
                10260.932208618978205,
                45507.789335026729956
    };

    private static final double[] RAWCDF_C = {
        0.39894151208813466764,
                8.8831497943883759412,
                93.506656132177855979,
                597.27027639480026226,
                2494.5375852903726711,
                6848.1904505362823326,
                11602.651437647350124,
                9842.7148383839780218,
                1.0765576773720192317e-8
    };

    private static final double[] RAWCDF_D = {
        22.266688044328115691,
                235.38790178262499861,
                1519.377599407554805,
                6485.558298266760755,
                18615.571640885098091,
                34900.952721145977266,
                38912.003286093271411,
                19685.429676859990727
    };

    private static final double[] RAWCDF_P = {
        0.21589853405795699,
                0.1274011611602473639,
                0.022235277870649807,
                0.001421619193227893466,
                2.9112874951168792e-5,
                0.02307344176494017303
    };

    private static final double[] RAWCDF_Q = {
        1.28426009614491121,
                0.468238212480865118,
                0.0659881378689285515,
                0.00378239633202758244,
                7.29751555083966205e-5
    };

    private static final double SQRT_32 = Math.sqrt(32);

    private static final double INV_SQRT_2PI = 1.0 /Math.sqrt(2 * Math.PI);

    private static double[] rawCDF(final double x, boolean log) {
        if  (Double.isNaN(x))
            return new double[] {Double.NaN, Double.NaN};
        final double eps = DoubleConstants.DBL_EPSILON * .5;
        final double min = Double.MIN_NORMAL;
        final double[] result = new double[2];


                double xden, xnum, temp, del, xsq, y;
                int i;


                y = Math.abs(x);
                if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
                    if (y > eps) {
                        xsq = x * x;
                        xnum = RAWCDF_A[4] * xsq;
                        xden = xsq;
                        for (i = 0; i < 3; ++i) {
                            xnum = (xnum + RAWCDF_A[i]) * xsq;
                            xden = (xden + RAWCDF_B[i]) * xsq;
                        }
                    } else xnum = xden = 0.0;

                    temp = x * (xnum + RAWCDF_A[3]) / (xden + RAWCDF_B[3]);
                    result[0] = log ? Math.log(0.5 + temp) : 0.5 + temp;
                    result[1] = log ? Math.log(0.5 - temp) : 0.5 - temp;
                }
                else if (y <= SQRT_32) {

	/* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

                    xnum = RAWCDF_C[8] * y;
                    xden = y;
                    for (i = 0; i < 7; ++i) {
                        xnum = (xnum + RAWCDF_C[i]) * y;
                        xden = (xden + RAWCDF_D[i]) * y;
                    }
                    temp = (xnum + RAWCDF_C[7]) / (xden + RAWCDF_D[7]);


                    xsq = Doubles.truncate(y * SIXTEN) / SIXTEN;
                    del = (y - xsq) * (y + xsq);
                    if(log) {
                        result[0] = (-xsq * xsq * 0.5) + (-del * 0.5) + Math.log(temp);
                        result[1] = Math.log1p(-Math.exp(-xsq * xsq * 0.5) *
                                Math.exp(-del * 0.5) * temp);
                    } else {
                        result[0] = Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp;
                        result[1] = 1.0 - result[0];
                    }

                    if (x > 0.) {/* swap  ccum <--> cum */
                        temp = result[0];
                        result[0] = result[1];
                        result[1] = temp;
                    }

                }

/* else	  |x| > sqrt(32) = 5.657 :
 * the next two case differentiations were really for lower=T, log=F
 * Particularly	 *not*	for  log_p !
 * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
 *
 * Note that we do want symmetry(0), lower/upper -> hence use y
 */
                else if((log && y < 1e170) /* avoid underflow below */
	/*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
	 * Then, make use of  Abramowitz & Stegun, 26.2.13, something like
	 xsq = x*x;
	 if(xsq * DBL_EPSILON < 1.)
	    del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
	 else
	    del = 0.;
	 *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
	 *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./
 	 swap_tail;
	 [Yes, but xsq might be infinite.]
	*/
                        || (-37.5193 < x  &&  x < 8.2924)
                        || (-8.2924  < x  &&  x < 37.5193)
                        ) {

	/* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
                    xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
                    xnum = RAWCDF_P[5] * xsq;
                    xden = xsq;
                    for (i = 0; i < 4; ++i) {
                        xnum = (xnum + RAWCDF_P[i]) * xsq;
                        xden = (xden + RAWCDF_Q[i]) * xsq;
                    }
                    temp = xsq * (xnum + RAWCDF_P[4]) / (xden + RAWCDF_Q[4]);
                    temp = (INV_SQRT_2PI - temp) / y;

                    del = (x - xsq) * (x + xsq);
                    if(log) {
                        result[0] = (-xsq * xsq * 0.5) + (-del * 0.5) + Math.log(temp);
                        result[1] = Math.log1p(-Math.exp(-xsq * xsq * 0.5) *
                                Math.exp(-del * 0.5) * temp);
                    } else {
                        result[0] = Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp;
                        result[1] = 1.0 - result[0];
                    }

                    if (x > 0.) {/* swap  ccum <--> cum */
                        temp = result[0];
                        result[0] = result[1];
                        result[1] = temp;
                    }
                } else { /* large x such that probs are 0 or 1 */
                    if(x > 0) {
                        result[0] = log ? 0.0 : 1.0;
                        result[1] = log ? Double.NEGATIVE_INFINITY : 0.0;
                    }
                    else {
                        result[0] = log ? Double.NEGATIVE_INFINITY : 0.0;
                        result[1] = log ? 0.0 : 1.0;
                    }
                }


    //            #ifdef NO_DENORMS
    /* do not return "denormalized" -- we do in R */
    //            if(log_p) {
    //                if(*cum > -min)	 *cum = -0.;
    //                if(*ccum > -min)*ccum = -0.;
    //            }
    //            else {
    //                if(*cum < min)	 *cum = 0.;
    //                if(*ccum < min)	*ccum = 0.;
    //            }
    //            #endif
                return result;
            }
}
