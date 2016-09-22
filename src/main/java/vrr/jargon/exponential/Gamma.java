package vrr.jargon.exponential;

import vrr.jargon.exceptions.UnderflowException;
import vrr.jargon.utils.DoubleConstants;
import vrr.jargon.utils.Doubles;
import vrr.jargon.utils.MathConstants;
import vrr.jargon.utils.Trigonometry;

import java.util.Arrays;

/**
 * Created by valentin on 9/20/16.
 */
public class Gamma {

    private static final double LN2 = Math.log(2);
    private static final double M_cutoff = LN2 * Double.MAX_EXPONENT / DoubleConstants.DBL_EPSILON;

    public static double CDF(final double unscaledX, final double shape, final double scale, final boolean log) {
        if (Double.isNaN(unscaledX))
            return Double.NaN;
        else if (!Double.isFinite(shape) || shape < 0)
            return Double.NaN;
        else if (!Double.isFinite(scale) || scale <= 0)
            return Double.NaN;
        else if (Double.isNaN(unscaledX))
            return Double.NaN;
        else if (shape == 0)
            return unscaledX <= 0 ? (log ? Double.NEGATIVE_INFINITY : 0) : (log ? 0.0 : 1.0);
        else {
            final double x = unscaledX / scale; // after this transformation we can forget about the shape entirely.

            if (x <= 0)
                return log ? Double.NEGATIVE_INFINITY : 0;
            else if (x == Double.NEGATIVE_INFINITY)
                return log ? 0.0 : 1.0;
            else {
                final double result;
                if (x < 1)
                    result = CDF_smallx(x, shape, log);
                else if (x <= shape - 1 && x < 0.8 * (shape + 50)) {
                    double sum = pd_upper_series(x, shape, log);/* = x/alph + o(x/alph) */
                    double d = dpois_wrap(shape, x, log);
                    result = log ? sum + d : sum * d;
                } else if (shape - 1 < x && shape < 0.8 * (x + 50)) {
                    double sum;
                    double d = dpois_wrap(shape, x, log);
                    if (shape < 1) {
                        if (x * DoubleConstants.DBL_EPSILON > 1 - shape)
                            sum = log ? 0.0 : 1.0;
                        else {
                            double f = pd_lower_cf(shape, x - (shape - 1)) * x / shape;
                            sum = log ? Math.log(f) : f;
                        }
                    } else {
                        sum = pd_lower_series(x, shape - 1);
                        sum = log ? Math.log1p(sum) : 1 + sum;
                    }
                    result = log
                            ? d + sum > -LN2 ? Math.log(-Math.expm1(d + sum)) : Math.log1p(-Math.exp(d + sum))
                            : 1 - d * sum;
                } else {
                    result = ppois_asymp_upper_tail(shape - 1, x, log);
                }

                if (!log && result < Double.MIN_NORMAL / DoubleConstants.DBL_EPSILON) {
                    return Math.exp(CDF(x, shape, 1.0, true));
                } else
                    return result;
            }
        }
    }

    private static double CDF_smallx(final double x, final double shape, final boolean log) {
        double sum = 0, c = shape, n = 0, term;
        do {
            n++;
            c *= -x / n;
            term = c / (shape + n);
            sum += term;
        } while (Math.abs (term) > DoubleConstants.DBL_EPSILON * Math.abs(sum));


        double f1 = log ? Math.log1p (sum) : 1 + sum;
        double f2;
        if (shape > 1) {
            f2 = Poisson.density(shape, x, log);
            f2 = log ? f2 + x : f2 * Math.exp(x);
        } else if (log)
            f2 = shape * Math.log(x) - logGammaPlus1(shape);
        else
            f2 = Math.pow (x, shape) / Math.exp (logGammaPlus1(shape));

        return log ? f1 + f2 : f1 * f2;
    }

    /**
     *
     * Based on https://github.com/wch/r-source/blob/e5b21d0397c607883ff25cca379687b86933d730/src/nmath/qgamma.c
     * where lower_tail is always true.
     *
     * @param shape the shape parameter of the gamma distribution.
     * @param scale the scale parameter of the gamma distribution.
     * @param log whether the return should be logged.
     */
    public static double inverseCDF(final double p, final double shape, final double scale, final boolean log) {
        if (Double.isNaN(p) || Double.isNaN(shape) || Double.isNaN(scale))
            return Double.NaN;
        if (!log) {
            if (p < 0.0 || p > 1.0)
                return Double.NaN;
            else if (p == 0.0)
                return 0.0;
            else if (p == 1.0)
                return Double.POSITIVE_INFINITY;
        } else {
            if (p > 0)
                return Double.NaN;
            else if (p == 0)
                return Double.POSITIVE_INFINITY;
            else if (p == Double.NEGATIVE_INFINITY)
                return 0.0;
        }

        if (shape < 0)
            return Double.NaN;
        else if (scale <= 0)
            return Double.NaN;



        final double EPS1 = 1e-2;
        final double EPS2 = 5e-7;
        final double EPS_N = 1e-15;
//Never used
//        final double LN_EPS = -36.043653389117156;
        final double MAXIT = 1000;

        final double pMIN = 1e-100;
        final double pMAX = 1-1e-14;

        final double i420  = 1./ 420.;
        final double i2520 = 1./ 2520.;
        final double i5040 = 1./ 5040.;

        double p_, a, b, c, g, ch, ch0, p1;
        double p2, q, s1, s2, s3, s4, s5, s6, t, x;
        int i, max_it_Newton = 1;

        if (shape < 1e-10) // Warning value of shape (%g) is extremely small: results may be unreliable
            max_it_Newton = 7;

        p_ = log ? Math.exp(p) : p;

        g = logGamma(shape);

        ch = chiSqApproximateInverseCDF(p, 2 * shape, g, log, EPS1);

        if (!Double.isFinite(ch)) {
            max_it_Newton = 0;
        } else if (p_ > pMAX || p_ < pMIN) {
            max_it_Newton = 20;
        } else {
            c = shape - 1;
            s6 = (120+c*(346+127*c)) * i5040;
            ch0 = ch;
            for(i=1; i <= MAXIT; i++ ) {
                q = ch;
                p1 = 0.5 * ch;
                p2 = p_ - CDF(p1, shape, 1, false);
                if (!Double.isFinite(p2) || ch <= 0) {
                    ch = ch0;
                    max_it_Newton = 27;
                    break;
                }
                t = p2* Math.exp(shape* MathConstants.LN_2 +g+p1-c*Math.log(ch));
                b = t/ch;
                a = 0.5*t - b*c;
                s1 = (210+ a*(140+a*(105+a*(84+a*(70+60*a))))) * i420;
                s2 = (420+ a*(735+a*(966+a*(1141+1278*a)))) * i2520;
                s3 = (210+ a*(462+a*(707+932*a))) * i2520;
                s4 = (252+ a*(672+1182*a) + c*(294+a*(889+1740*a))) * i5040;
                s5 = (84+2264*a + c*(1175+606*a)) * i2520;
                ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
                if (Math.abs(q - ch) < EPS2*ch) {
                    break;
                }
                if(Math.abs(q - ch) > 0.1*ch) {
                    if(ch < q) ch = 0.9 * q; else ch = 1.1 * q;
                }
            }
        }

        x = 0.5*scale*ch;
        if (max_it_Newton > 0) {
            final double loggedP = log ? p : Math.log(p);
            if (x == 0) {
                final double _1_p = 1. + 1e-7;
                x = Double.MIN_NORMAL;
                p_ = CDF(x, shape, scale, true);
                if((p_ > loggedP * _1_p))
                    return(0.);
            } else {
                p_ = CDF(x, shape, scale, true);
                if(p_ == Double.NEGATIVE_INFINITY) return 0;
                for(i = 1; i <= max_it_Newton; i++) {
                    p1 = p_ - loggedP;
                    if (Math.abs(p1) < Math.abs(EPS_N * loggedP))
                        break;
                    if ((g = density(x, shape, scale, true)) == Double.NEGATIVE_INFINITY)
                        break;
                    t = p1 * Math.exp(p_ - g);
                    t = x - t;
                    p_ = CDF(t, shape, scale, true);
                    if (Math.abs(p_ - loggedP) > Math.abs(p1) ||
                            (i > 1 && Math.abs(p_ - loggedP) == Math.abs(p1)) /* <- against flip-flop */)
                        break;
//#ifdef Harmful_notably_if_max_it_Newton_is_1
//                    if(t > 1.1*x) t = 1.1*x;
//                    else if(t < 0.9*x) t = 0.9*x;
//#endif
                    x = t;
                }
		/* no improvement */
            }
        }

        return x;
    }

    public static double density(double x, double shape, double scale, boolean log) {
        if (Double.isNaN(x) || Double.isNaN(shape) || Double.isNaN(scale))
            return Double.NaN;
        else if (shape < 0 || scale <= 0)
            return Double.NaN;
        else if (x < 0)
            return log ? Double.NEGATIVE_INFINITY : 0.0;
        else if (x == 0) {
            if (shape < 1)
                return Double.POSITIVE_INFINITY;
            else if (shape > 1)
                return log ? Double.NEGATIVE_INFINITY : 0.0;
            else
                return log ? -Math.log(scale) : 1 / scale;
        } else if (shape < 1) {
            final double pr = Poisson.density(shape, x / scale, log);
            return log ? pr + Math.log(shape / x) : pr * shape / x;
        } else {
            final double pr = Poisson.density(shape - 1, x / scale, log);
            return log ? pr - Math.log(scale) : pr / scale;
        }
    }

    private static double chiSqApproximateInverseCDF(final double p, final double nu, final double g, final boolean log,
                                              final double tol) {
        if (Double.isNaN(p) || Double.isNaN(nu))
            return Double.NaN;
        if (log && p > 0 || (!log && (p < 0 || p > 1)))
            return Double.NaN;
        if (nu <= 0)
            return Double.NaN;


        final double C7	= 4.67;
        final double C8	= 6.66;
        final double C9 = 6.73;
        final double C10 = 13.32;

        double alpha, a, c, ch, p1;
        double p2, q, t, x;

        alpha = 0.5 * nu;/* = [pq]gamma() shape */
        c = alpha-1;

        if (nu < (-1.24)* (p1 = log ? p : Math.log(p))) {	/* for small chi-squared */
            double lgam1pa = (alpha < 0.5) ? logGammaPlus1(alpha) : (Math.log(alpha) + g);
            ch = Math.exp((lgam1pa + p1)/alpha + MathConstants.LN_2);

        } else if(nu > 0.32) {	/*  using Wilson and Hilferty estimate */

            x = Gaussian.inverseCDF(p, 0, 1, log);
            p1 = 2./(9*nu);
            ch = nu*Math.pow(x* Math.sqrt(p1) + 1-p1, 3);

	/* approximation for p tending to 1: */
            if( ch > 2.2*nu + 6 )
                ch = -2*((log ? (p > - MathConstants.LN_2 ? Math.log(-Math.expm1(p)) : Math.log1p(-Math.exp(p))) : Math.log1p(-p)) - c* Math.log(0.5*ch) + g);

        } else { /* "small nu" : 1.24*(-log(p)) <= nu <= 0.32 */

            ch = 0.4;
            a =(log ? (p > - MathConstants.LN_2 ? Math.log(-Math.expm1(p)) : Math.log1p(-Math.exp(p))) : Math.log1p(-p)) + g + c* MathConstants.LN_2;

            do {
                q = ch;
                p1 = 1. / (1+ch*(C7+ch));
                p2 = ch*(C9+ch*(C8+ch));
                t = -0.5 +(C7+2*ch)*p1 - (C9+ch*(C10+3*ch))/p2;
                ch -= (1- Math.exp(a+0.5*ch)*p2*p1)/t;
            } while(Math.abs(q - ch) > tol * Math.abs(ch));
        }

        return ch;
    }

    private static double dpois_wrap(final double x_plus_1, final double lambda, final boolean give_log) {
        if (!Double.isFinite(lambda))
            return give_log ? Double.NEGATIVE_INFINITY : 0.0;
        if (x_plus_1 > 1)
            return Poisson.density(x_plus_1 - 1, lambda, give_log);
        if (lambda > Math.abs(x_plus_1 - 1) * M_cutoff)
            return give_log ? (-lambda - logGamma(x_plus_1)) : Math.exp(-lambda - logGamma(x_plus_1));
        else {
            double d = Poisson.density(x_plus_1, lambda, give_log);

            return give_log
                    ? d + Math.log(x_plus_1 / lambda)
                    : d * (x_plus_1 / lambda);
        }
    }


    private static final double[] LOG_GAMMA_P1_COEFFS = {
            0.3224670334241132182362075833230126e-0,/* = (zeta(2)-1)/2 */
            0.6735230105319809513324605383715000e-1,/* = (zeta(3)-1)/3 */
            0.2058080842778454787900092413529198e-1,
            0.7385551028673985266273097291406834e-2,
            0.2890510330741523285752988298486755e-2,
            0.1192753911703260977113935692828109e-2,
            0.5096695247430424223356548135815582e-3,
            0.2231547584535793797614188036013401e-3,
            0.9945751278180853371459589003190170e-4,
            0.4492623673813314170020750240635786e-4,
            0.2050721277567069155316650397830591e-4,
            0.9439488275268395903987425104415055e-5,
            0.4374866789907487804181793223952411e-5,
            0.2039215753801366236781900709670839e-5,
            0.9551412130407419832857179772951265e-6,
            0.4492469198764566043294290331193655e-6,
            0.2120718480555466586923135901077628e-6,
            0.1004322482396809960872083050053344e-6,
            0.4769810169363980565760193417246730e-7,
            0.2271109460894316491031998116062124e-7,
            0.1083865921489695409107491757968159e-7,
            0.5183475041970046655121248647057669e-8,
            0.2483674543802478317185008663991718e-8,
            0.1192140140586091207442548202774640e-8,
            0.5731367241678862013330194857961011e-9,
            0.2759522885124233145178149692816341e-9,
            0.1330476437424448948149715720858008e-9,
            0.6422964563838100022082448087644648e-10,
            0.3104424774732227276239215783404066e-10,
            0.1502138408075414217093301048780668e-10,
            0.7275974480239079662504549924814047e-11,
            0.3527742476575915083615072228655483e-11,
            0.1711991790559617908601084114443031e-11,
            0.8315385841420284819798357793954418e-12,
            0.4042200525289440065536008957032895e-12,
            0.1966475631096616490411045679010286e-12,
            0.9573630387838555763782200936508615e-13,
            0.4664076026428374224576492565974577e-13,
            0.2273736960065972320633279596737272e-13,
            0.1109139947083452201658320007192334e-13/* = (zeta(40+1)-1)/(40+1) */
    };

    private static double logGammaPlus1(final double x) {

    /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
        final int N = 40;


        final double c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
        final double tol_logcf = 1e-14;
        double lgam;
        int i;

        if (Math.abs (x) >= 0.5)
            return logGamma(x + 1);

    /* Abramowitz & Stegun 6.1.33 : for |x| < 2,
     * <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
     * where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
     *
     * Here, another convergence acceleration trick is used to compute
     * lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
     */
        lgam = c * logcf(-x / 2, N + 2, 1, tol_logcf);
        for (i = N - 1; i >= 0; i--)
            lgam = LOG_GAMMA_P1_COEFFS[i] - x * lgam;

        return (x * lgam - MathConstants.EULERS_CONSTANT) * x - log1pmx (x);
    }

    /* Accurate calculation of log(1+x)-x, particularly for small x.  */
    private static double log1pmx (double x)
    {
        final double minLog1Value = -0.79149064;

        if (x > 1 || x < minLog1Value)
            return Math.log1p(x) - x;
        else {
            /* -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
	    * log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
	    * ---------------------------------------------
	    * S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
	   */
            double r = x / (2 + x);
            double y = r * r;
            if (Math.abs(x) < 1e-2) {
                final double two = 2;
                return r * ((((two / 9 * y + two / 7) * y + two / 5) * y +
                        two / 3) * y - x);
            } else {
                final double tol_logcf = 1e-14;
                return r * (2 * y * logcf (y, 3, 2, tol_logcf) - x);
            }
        }
    }


    public static double logGamma(final double x) {
        double ans, y, sinpiy;

        double xmax = Double.MAX_VALUE / Math.log(Double.MAX_VALUE);
        double dxrel = Math.sqrt(DoubleConstants.DBL_EPSILON);
        if (Double.isNaN(x)) return Double.NaN;

        if (x <= 0 && x == Math.ceil(x))  /* Negative integer argument */
            return Double.POSITIVE_INFINITY;

        y = Math.abs(x);

        if (y < 1e-306)
            return - Math.log(y); // denormalized range, R change
        else if (y <= 10)
            return Math.log(Math.abs(gamma(x)));
        else if (y > xmax)
            return Double.POSITIVE_INFINITY;
        else if (x > 0) {
            if(x > 1e17)
                return(x*(Math.log(x) - 1.));
            else if(x > 4934720.)
                return(MathConstants.LN_SQRT_2PI + (x - 0.5) * Math.log(x) - x);
            else
                return MathConstants.LN_SQRT_2PI + (x - 0.5) * Math.log(x) - x + logGammaCorrection(x);
        }
    /* else: x < -10; y = -x */
        sinpiy = Math.abs(Trigonometry.sinPi(y));

        if (sinpiy == 0) { /* Negative integer argument ===
			  Now UNNECESSARY: caught above */
            return Double.NaN;
        }

        ans = MathConstants.LN_SQRT_PId2 + (x - 0.5) * Math.log(y) - x - Math.log(sinpiy) - logGammaCorrection(y);

        if (Math.abs((x - Math.ceil(x - 0.5)) * ans / x) < dxrel) {
            //WARNING/ERROR:
            //throw new PrecissionException(ans, "logGamma imprecise value as x = " + x + " is to close to a negative integer: " + ans);
        }
        return ans;

    }


    //Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
    //#define SQR(x) ((x)*(x))
    //static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
    //#undef SQR
    private static final double scalefactor;

    static {
        double x = 4294967296.0;
        for (int i = 0; i < 3; i++) {
            x *= x;
        }
        scalefactor = x;
    }

    /* Continued fraction for calculation of
     *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
     *
     * auxilary in log1pmx() and lgamma1p()
     */
    private static double  logcf (double x, double i, double d,
           double eps /* ~ relative tolerance */)
    {
        double c1 = 2 * d;
        double c2 = i + d;
        double c4 = c2 + d;
        double a1 = c2;
        double b1 = i * (c2 - i * x);
        double b2 = d * d * x;
        double a2 = c4 * c2 - b2;



            b2 = c4 * b1 - i * b2;

        while (Math.abs(a2 * b1 - a1 * b2) > Math.abs(eps * b1 * b2)) {
            double c3 = c2*c2*x;
            c2 += d;
            c4 += d;
            a1 = c4 * a2 - c3 * a1;
            b1 = c4 * b2 - c3 * b1;

            c3 = c1 * c1 * x;
            c1 += d;
            c4 += d;
            a2 = c4 * a1 - c3 * a2;
            b2 = c4 * b1 - c3 * b2;

            if (Math.abs(b2) > scalefactor) {
                a1 /= scalefactor;
                b1 /= scalefactor;
                a2 /= scalefactor;
                b2 /= scalefactor;
            } else if (Math.abs(b2) < 1 / scalefactor) {
                a1 *= scalefactor;
                b1 *= scalefactor;
                a2 *= scalefactor;
                b2 *= scalefactor;
            }
        }

        return a2 / b2;
    }

    /**
     * Original sequence from the c code; want to keep it around although we actually only use the first 5 value.
     */
    private static final double[] LOG_GAMMA_CORRECTION_SERIES = {
            +.1666389480451863247205729650822e+0,
            -.1384948176067563840732986059135e-4,
            +.9810825646924729426157171547487e-8,
            -.1809129475572494194263306266719e-10,
            +.6221098041892605227126015543416e-13,
            -.3399615005417721944303330599666e-15,
            +.2683181998482698748957538846666e-17,
            -.2868042435334643284144622399999e-19,
            +.3962837061046434803679306666666e-21,
            -.6831888753985766870111999999999e-23,
            +.1429227355942498147573333333333e-24,
            -.3547598158101070547199999999999e-26,
            +.1025680058010470912000000000000e-27,
            -.3401102254316748799999999999999e-29,
            +.1276642195630062933333333333333e-30
    };

    private static final Chebyshev LOG_GAMMA_CORRECTION_CSERIES =
            new Chebyshev(Arrays.copyOfRange(LOG_GAMMA_CORRECTION_SERIES, 0, 5));

    private static double logGammaCorrection(final double x)
    {
        final double xbig = 94906265.62425156;
        final double xmax = 3.745194030963158e306;
        double tmp;

        if (x < 10)
            return Double.NaN;
        else if (x >= xmax) {
            throw new UnderflowException("lgammacor underflow");
        }
        else if (x < xbig) {
            tmp = 10 / x;
            return LOG_GAMMA_CORRECTION_CSERIES.evaluate(tmp * tmp * 2 - 1) / x;
        }
        return 1 / (x * 12);
    }

  private static final Chebyshev GAMMA_CSERIES = new Chebyshev(DoubleConstants.DBL_EPSILON / 20, new double[] {
                +.8571195590989331421920062399942e-2,
                +.4415381324841006757191315771652e-2,
                +.5685043681599363378632664588789e-1,
                -.4219835396418560501012500186624e-2,
                +.1326808181212460220584006796352e-2,
                -.1893024529798880432523947023886e-3,
                +.3606925327441245256578082217225e-4,
                -.6056761904460864218485548290365e-5,
                +.1055829546302283344731823509093e-5,
                -.1811967365542384048291855891166e-6,
                +.3117724964715322277790254593169e-7,
                -.5354219639019687140874081024347e-8,
                +.9193275519859588946887786825940e-9,
                -.1577941280288339761767423273953e-9,
                +.2707980622934954543266540433089e-10,
                -.4646818653825730144081661058933e-11,
                +.7973350192007419656460767175359e-12,
                -.1368078209830916025799499172309e-12,
                +.2347319486563800657233471771688e-13,
                -.4027432614949066932766570534699e-14,
                +.6910051747372100912138336975257e-15,
                -.1185584500221992907052387126192e-15,
                +.2034148542496373955201026051932e-16,
                -.3490054341717405849274012949108e-17,
                +.5987993856485305567135051066026e-18,
                -.1027378057872228074490069778431e-18,
                +.1762702816060529824942759660748e-19,
                -.3024320653735306260958772112042e-20,
                +.5188914660218397839717833550506e-21,
                -.8902770842456576692449251601066e-22,
                +.1527474068493342602274596891306e-22,
                -.2620731256187362900257328332799e-23,
                +.4496464047830538670331046570666e-24,
                -.7714712731336877911703901525333e-25,
                +.1323635453126044036486572714666e-25,
                -.2270999412942928816702313813333e-26,
                +.3896418998003991449320816639999e-27,
                -.6685198115125953327792127999999e-28,
                +.1146998663140024384347613866666e-28,
                -.1967938586345134677295103999999e-29,
                +.3376448816585338090334890666666e-30,
                -.5793070335782135784625493333333e-31});

    /**
     * Minimum bound of the input of {@link #gamma}, precalulated
     * for IEEE_754.
     * <p>
     * Notice that original code includes a routine to calculate for other architectures. Might
     * bet worth to port it although Java tends to be IEEE 754.
     * </p>
     */
    private final static double IEEE_754_GAMMA_XMIN = -170.5674972726612;

    /**
     * Maximum bound of the input of {@lonk #gamma}, precalculated
     * for IEEE_754.
     * <p>
     * Notice that original code includes a routine to calculate for other architectures. Might
     * bet worth to port it although Java tends to be IEEE 754.
     * </p>
     */
    private final static double IEEE_754_GAMMA_XMAX =  171.61447887182298;

    private final static double XSML = Math.exp(Math.max(Math.log(Double.MIN_NORMAL), -Math.log(Double.MAX_VALUE)) + 0.01);

    private final static double DXREL = Math.sqrt(DoubleConstants.DBL_EPSILON);

    public static double gamma(final double x) {
        if (Double.isNaN(x))
            return Double.NaN;
        int i, n;
        double y;
        double sinpiy, value;


    /* If the argument is exactly zero or a negative integer
     * then return NaN. */
        if ((x == 0) || ((x < 0) && (x == Math.ceil(x))))
            return Double.NaN;

        y = Math.abs(x);

        if (y <= 10) {

	/* Compute gamma(x) for -10 <= x <= 10
	 * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
	 * first of all. */

            n = (int) x;
            if(x < 0) --n;
            y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
            --n;
            value = GAMMA_CSERIES.evaluate(y * 2 - 1) + .9375;
            if (n == 0)
                return value;/* x = 1.dddd = 1+y */

            if (n < 0) {
	    /* compute gamma(x) for -10 <= x < 1 */

	    /* exact 0 or "-n" checked already above */

	    /* The answer is less than half precision */
	    /* because x too near a negative integer. */
//                if (x < -0.5 && Math.abs(x - (int)(x - 0.5) / x) < DXREL) {
//                    ML_ERROR(ME_PRECISION, "gammafn");
//                }

                if (y < XSML) {
	    /* The argument is so close to 0 that the result would overflow. */
        //          ML_ERROR(ME_RANGE, "gammafn");
                    return x > 0 ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
                }

                n = -n;

                for (i = 0; i < n; i++) {
                    value /= (x + i);
                }
                return value;
            }
            else {
	    /* gamma(x) for 2 <= x <= 10 */

                for (i = 1; i <= n; i++) {
                    value *= (y + i);
                }
                return value;
            }
        }
        else {
	/* gamma(x) for	 y = |x| > 10. */

            if (x > IEEE_754_GAMMA_XMAX) {			/* Overflow */
                return Double.POSITIVE_INFINITY;
            }

            if (x < IEEE_754_GAMMA_XMIN) {			/* Underflow */
                return 0.;
            }

            if(y <= 50 && y == Doubles.truncate(y)) { /* compute (n - 1)! */
                value = 1.;
                for (i = 2; i < y; i++) value *= i;
            }
            else { /* normal case */
                value = Math.exp((y - 0.5) * Math.log(y) - y + MathConstants.LN_SQRT_2PI +
                        ((2*y == Doubles.truncate(2*y)? Stirling.logError(y) : logGammaCorrection(y))));
            }
            if (x > 0)
                return value;
        // WARN:
        //    if (Math.abs((x - (int)(x - 0.5))/x) < DXREL){
        //
	    // /* The answer is less than half precision because */
	    // /* the argument is too near a negative integer. */
        //
        //        ML_ERROR(ME_PRECISION, "gammafn");
        //    }

            sinpiy = Trigonometry.sinPi(y);
            if (sinpiy == 0) {		/* Negative integer arg - overflow */
                //WARN:
                //        ML_ERROR(ME_RANGE, "gammafn");
                return Double.POSITIVE_INFINITY;
            }

            return - Math.PI / (y * sinpiy * value);
        }
    }

    private static double pd_upper_series (final double x, final double y, final boolean log)
    {
        double denominator = y;
        double term = x / denominator;
        double sum = term;


        do {
            denominator++;
            term *= x / denominator;
            sum += term;
        } while (term > sum * DoubleConstants.DBL_EPSILON);

        return log ? Math.log (sum) : sum;
    }

    /* Continued fraction for calculation of
 *    scaled upper-tail F_{gamma}
 *  ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]
 */
    private static double  pd_lower_cf (final double y, final double d)
    {
        double f= 0.0, of, f0;
        double i, c2, c3, c4,  a1, b1,  a2, b2;

        final double max_it = 200000;

        if (y == 0) return 0;

        f0 = y/d;
    /* Needed, e.g. for  pgamma(10^c(100,295), shape= 1.1, log=TRUE): */
        if(Math.abs(y - 1) < Math.abs(d) * DoubleConstants.DBL_EPSILON) { /* includes y < d = Inf */
            return (f0);
        }

        if(f0 > 1.) f0 = 1.;
        c2 = y;
        c4 = d; /* original (y,d), *not* potentially scaled ones!*/

        a1 = 0; b1 = 1;
        a2 = y; b2 = d;

        while (b2 > scalefactor) {
            a1 /= scalefactor;
            b1 /= scalefactor;
            a2 /= scalefactor;
            b2 /= scalefactor;
        }

            i = 0; of = -1.; /* far away */
        while (i < max_it) {

            i++;	c2--;	c3 = i * c2;	c4 += 2;
	/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd */
            a1 = c4 * a2 + c3 * a1;
            b1 = c4 * b2 + c3 * b1;

            i++;	c2--;	c3 = i * c2;	c4 += 2;
	/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even */
            a2 = c4 * a1 + c3 * a2;
            b2 = c4 * b1 + c3 * b2;

            if (b2 > scalefactor) {
                a1 /= scalefactor;
                b1 /= scalefactor;
                a2 /= scalefactor;
                b2 /= scalefactor;
            }

            if (b2 != 0) {
                f = a2 / b2;
	    /* convergence check: relative; "absolute" for very small f : */
                if (Math.abs (f - of) <= DoubleConstants.DBL_EPSILON * Math.max(f0, Math.abs(f))) {
                    return f;
                }
                of = f;
            }
        }

       // MATHLIB_WARNING(" ** NON-convergence in pgamma()'s pd_lower_cf() f= %g.\n",
        //        f);
        return f;/* should not happen ... */
    } /* pd_lower_cf() */


    private static double pd_lower_series (final double lambda, double y)
    {
        double term = 1, sum = 0;

        while (y >= 1 && term > sum * DoubleConstants.DBL_EPSILON) {
            term *= y / lambda;
            sum += term;
            y--;
        }

        if (y != Math.floor(y)) {
	/*
	 * The series does not converge as the terms start getting
	 * bigger (besides flipping sign) for y < -lambda.
	 */
            final double f = pd_lower_cf (y, lambda + 1 - y);
	/* FIXME: in quite few cases, adding  term*f  has no effect (f too small)
	 *	  and is unnecessary e.g. for pgamma(4e12, 121.1) */
            sum += term * f;
        }

        return sum;
    }

    /*
 * Asymptotic expansion to calculate the probability that Poisson variate
 * has value <= x.
 * Various assertions about this are made (without proof) at
 * http://members.aol.com/iandjmsmith/PoissonApprox.htm
 */
    private static final double[] PPOIS_ASYMP_UPPER_TAIL_A = new double[] {
        -1e99, /* placeholder used for 1-indexing */
                2/3.,
                -4/135.,
                8/2835.,
                16/8505.,
                -8992/12629925.,
                -334144/492567075.,
                698752/1477701225};

    private static final double[] PPOIS_ASYMP_UPPER_TAIL_B = new double[] {
            -1e99, /* placeholder */
            1/12.,
            1/288.,
            -139/51840.,
            -571/2488320.,
            163879/209018880.,
            5246819/75246796800.,
            -534703531/902961561600.};

    private static double ppois_asymp_upper_tail (double x, double lambda, final boolean log)
    {
        double elfb, elfb_term;
        double res12, res1_term, res1_ig, res2_term, res2_ig;
        double dfm, pt_, s2pt, f, np;
        int i;

        dfm = lambda - x;
    /* If lambda is large, the distribution is highly concentrated
       about lambda.  So representation error in x or lambda can lead
       to arbitrarily large values of pt_ and hence divergence of the
       coefficients of this approximation.
    */
        pt_ = - log1pmx (dfm / x);
        s2pt = Math.sqrt (2 * x * pt_);
        if (dfm < 0) s2pt = -s2pt;

        res12 = 0;
        res1_ig = res1_term = Math.sqrt (x);
        res2_ig = res2_term = s2pt;
        for (i = 1; i < 8; i++) {
            res12 += res1_ig * PPOIS_ASYMP_UPPER_TAIL_A[i];
            res12 += res2_ig * PPOIS_ASYMP_UPPER_TAIL_B[i];
            res1_term *= pt_ / i ;
            res2_term *= 2 * pt_ / (2 * i + 1);
            res1_ig = res1_ig / x + res1_term;
            res2_ig = res2_ig / x + res2_term;
        }

        elfb = x;
        elfb_term = 1;
        for (i = 1; i < 8; i++) {
            elfb += elfb_term * PPOIS_ASYMP_UPPER_TAIL_B[i];
            elfb_term /= x;
        }

        elfb = -elfb;
        f = res12 / elfb;

        np = Gaussian.CDF(s2pt, 0.0, 1.0, log);

        if (log) {
            double n_d_over_p = dpnorm (s2pt, true, np);

            return np + Math.log1p (f * n_d_over_p);
        } else {
            double nd = Gaussian.density(s2pt, 0., 1., log);

            return np + f * nd;
        }
    }

    /*
 * Compute the following ratio with higher accuracy that would be had
 * from doing it directly.
 *
 *		 dnorm (x, 0, 1, FALSE)
 *	   ----------------------------------
 *	   pnorm (x, 0, 1, lower_tail, FALSE)
 *
 * Abramowitz & Stegun 26.2.12
 */
    private static double  dpnorm (double x, boolean lower_tail, double lp)
    {
    /*
     * So as not to repeat a pnorm call, we expect
     *
     *	 lp == pnorm (x, 0, 1, lower_tail, TRUE)
     *
     * but use it only in the non-critical case where either x is small
     * or p==exp(lp) is close to 1.
     */

        if (x < 0) {
            x = -x;
            lower_tail = !lower_tail;
        }

        if (x > 10 && !lower_tail) {
            double term = 1 / x;
            double sum = term;
            double x2 = x * x;
            double i = 1;

            do {
                term *= -i / x2;
                sum += term;
                i += 2;
            } while (Math.abs(term) > DoubleConstants.DBL_EPSILON * sum);

            return 1 / sum;
        } else {
            double d = Gaussian.density(x, 0., 1., false);
            return d / Math.exp (lp);
        }
    }


    public static double CDF(final double x, final double shape, final double scale) {
        return CDF(x, shape, scale, false);
    }

    public static double logCDF(final double x, final double shape, final double scale) {
        return CDF(x, shape, scale, true);
    }

    public static double inverseCDF(final double p, final double shape, final double scale) {
        return inverseCDF(p, shape, scale, false);
    }

    public static double logInverseCDF(final double p, final double shape, final double scale) {
        return inverseCDF(p, shape, scale, true);
    }

    public static double density(final double x, final double shape, final double scale) {
        return density(x, shape, scale, false);
    }

    public static double logDensity(final double x, final double shape, final double scale) {
        return density(x, shape, scale, true);
    }

}

