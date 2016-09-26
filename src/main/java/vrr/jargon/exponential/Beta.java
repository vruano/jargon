package vrr.jargon.exponential;

import vrr.jargon.utils.MathConstants;

import static java.lang.Math.*;

/**
 * Created by valentin on 9/22/16.
 */
public class Beta {

    private static final double LN_SQRT_2PI = log(sqrt(2 * PI));
    private static final double SQRT_PI = Math.sqrt(PI);
    private static final double acu_min = 1e-300;
    private static final double p_hi = 1 - 2.22e-16;
    private static final double fpu = 3e-308;
    private static final double p_lo = fpu;
    private static final double const1 = 2.30753;
    private static final double const2 = 0.27061;
    private static final double const3 = 0.99229;
    private static final double const4 = 0.04481;
    private static final byte DBL_MANT_DIG = 53;
    private static final short DBL_MIN_EXP = -1022;
    private static final short DBL_MAX_EXP = +1024;
    private static final double DBL_MIN = 2.2250738585072014e-308;
    private static final double DBL_very_MIN = DBL_MIN / 4.;
    private static final double DBL_log_v_MIN = MathConstants.LN_2 * (DBL_MIN_EXP - 2);
    private static final double DBL_EPSILON = 2.2204460492503131e-16;

    private static double[] return_q_0(final boolean give_log_q) {
        if (give_log_q) {
            return new double[]{Double.NEGATIVE_INFINITY, 0};
        } else {
            return new double[]{0, 1};
        }
    }

    private static double[] return_q_1(final boolean give_log_q) {
        if (give_log_q) {
            return new double[]{0, Double.NEGATIVE_INFINITY};
        } else {
            return new double[]{1, 0};
        }
    }

    private static double[] return_q_half(final boolean give_log_q) {
        if (give_log_q) {
            return new double[]{-MathConstants.LN_2, -MathConstants.LN_2};
        } else {
            return new double[]{0.5, 0.5};
        }
    }

    private static double R_DT_qIv(final boolean log_p, final boolean lower_tail, final double p) {
        return (log_p ? (lower_tail ? exp(p) : -expm1(p)) : (lower_tail ? (p) : (0.5 - (p) + 0.5)));
    }

    private static double R_DT_CIv(final boolean log_p, final boolean lower_tail, final double p) {
        return (log_p ? (lower_tail ? -expm1(p) : exp(p)) : (lower_tail ? (0.5 - (p) + 0.5) : (p)));
    }

    private static double R_D_LExp(final boolean log_p, final double p) {
        return (log_p ? ((p) > -MathConstants.LN_2 ? log(-expm1(p)) : log1p(-exp(p))) : log1p(-p));
    }

    private static double R_D_log(final boolean log_p, final double p) {
        return log_p ? (p) : log(p);
    }

    private static double R_Log1_Exp(final double p) {
        return ((p) > -MathConstants.LN_2 ? log(-expm1(p)) : log1p(-exp(p)));
    }

    private static double R_pow_di(double x, int n) {
        double pow = 1.0;

        if (Double.isNaN(x)) return Double.NaN;
        if (n != 0) {
            if (!Double.isFinite(x)) return pow(x, (double) n);
            if (n < 0) {
                n = -n;
                x = 1 / x;
            }
            while (true) {
                if ((n & 1) != 0) pow *= x;
                if ((n >>= 1) != 0) {
                    x *= x;
                } else {
                    break;
                }
            }
        }
        return pow;
    }

    public static double inverseCDF(final double q, final double alpha, final double beta, final boolean log) {
        return qbeta(q, alpha, beta, true, log);
    }

    public static double logInverseCDF(final double q, final double alpha, final double beta) {
        return inverseCDF(q, alpha, beta, true);
    }

    public static double inverseCDF(final double q, final double alpha, final double beta) {
        return inverseCDF(q, alpha, beta, false);
    }

    public static double beta(final double alpha, final double beta)
    {
        if (Double.isNaN(alpha) || Double.isNaN(beta))
            return Double.NaN;
        else if (alpha < 0 || beta < 0)
            return Double.NaN;
        else if (alpha == 0 || beta == 0)
            return Double.POSITIVE_INFINITY;
        else if (!Double.isFinite(alpha) || !Double.isFinite(beta))
            return 0.0;
        else if (alpha + beta < Gamma.IEEE_754_GAMMA_XMAX)
            return (1 / Gamma.gamma(alpha + beta)) * Gamma.gamma(alpha) * Gamma.gamma(beta);
        else {
            return Math.exp(logBeta(alpha, beta));
        }
    }

    private static double qbeta(double alpha, double p, double q, boolean lower_tail, boolean log_p) {

    /* test for admissibility of parameters */
        if (Double.isNaN(p) || Double.isNaN(q) || Double.isNaN(alpha))
            return p + q + alpha;
        if (p < 0. || q < 0.) return Double.NaN;
        // allowing p==0 and q==0  <==> treat as one- or two-point mass

        double[] qbet =
                qbeta_raw(alpha, p, q, lower_tail, log_p, -5, 4);
        return qbet[0];
    }

    private static double[] qbeta_raw(double alpha, double p, double q, boolean lower_tail, boolean log_p,
                                      double log_q_cut, /* if == Inf: return log(qbeta(..));
                   otherwise, if finite: the bound for
			       switching to log(x)-scale; see use_log_x */
                                      int n_N)  // number of "unconstrained" Newton steps before switching to constrained
    // = qb[0:1] = { qbeta(), 1 - qbeta() }
    {
        boolean
                swap_choose = true,
                swap_tail,
                log_, give_log_q = (log_q_cut == Double.POSITIVE_INFINITY),
                use_log_x = give_log_q, // or u < log_q_cut  below
                warned = false, add_N_step = true;
        int i_pb, i_inn;
        double a, la, logbeta, g, h, pp, p_, qq, r, s, t, w, y = -1.;
        double u, xinbta, tx;
        double u_n = 1.; // -Wall

        // Assuming p >= 0, q >= 0  here ...

        // Deal with boundary cases here:
        if (alpha == (log_p ? Double.NEGATIVE_INFINITY : 0.0)) {
            return return_q_0(give_log_q);
        }
        if (alpha == (log_p ? 0.0 : 1.0)) {
            return return_q_1(give_log_q);
        }

        // check alpha {*before* transformation which may all accuracy}:
        if ((log_p && alpha > 0) ||
                (!log_p && (alpha < 0 || alpha > 1))) { // alpha is outside
            return new double[] { Double.NaN, Double.NaN };
        }

        //  p==0, q==0, p = Inf, q = Inf  <==> treat as one- or two-point mass
        if (p == 0 || q == 0 || !Double.isFinite(p) || !Double.isFinite(q)) {
            if (p == 0 && q == 0) { // point mass 1/2 at each of {0,1} :
                if (alpha < (log_p ? -MathConstants.LN_2 : 0.5)) {
                    return return_q_0(give_log_q);
                } else if (alpha > (log_p ? -MathConstants.LN_2 : 0.5)) {
                    return return_q_1(give_log_q);
                } else {
                    return return_q_half(give_log_q);
                }
            } else if (p == 0 || p / q == 0) { // point mass 1 at 0 - "flipped around"
                return return_q_0(give_log_q);
            } else if (q == 0 || q / p == 0) { // point mass 1 at 0 - "flipped around"
                return return_q_1(give_log_q);
            } else {
                return return_q_half(give_log_q);
            }
        }


    /* initialize */
        p_ = R_DT_qIv(log_p, lower_tail, alpha);/* lower_tail prob (in any case) */
        // Conceptually,  0 < p_ < 1  (but can be 0 or 1 because of cancellation!)
        logbeta = logBeta(p, q);

        swap_tail = (p_ > 0.5);
        // change tail; default (swap_01 = NA): afterwards 0 < a <= 1/2
        if (swap_tail) { /* change tail, swap  p <-> q :*/
            a = R_DT_CIv(log_p, lower_tail, alpha); // = 1 - p_ < 1/2
    /* la := log(a), but without numerical cancellation: */
            la = lower_tail ? R_D_LExp(log_p, alpha) : R_D_log(log_p, alpha);
            pp = q;
            qq = p;
        } else {
            a = p_;
            la = lower_tail ? R_D_log(log_p, alpha) : R_D_LExp(log_p, alpha);
            pp = p;
            qq = q;
        }

    /* calculate the initial approximation */

    /* Desired accuracy for Newton iterations (below) should depend on  (a,p)
     * This is from Remark .. on AS 109, adapted.
     * However, it's not clear if this is "optimal" for IEEE double prec.
     * acu = fmax2(acu_min, pow(10., -25. - 5./(pp * pp) - 1./(a * a)));
     * NEW: 'acu' accuracy NOT for squared adjustment, but simple;
     * ---- i.e.,  "new acu" = sqrt(old acu)
    */
        double acu = max(acu_min, pow(10., -13. - 2.5 / (pp * pp) - 0.5 / (a * a)));
        // try to catch  "extreme left tail" early
        double u0 = (la + log(pp) + logbeta) / pp; // = log(x_0)
        final double
                log_eps_c = MathConstants.LN_2 * (1. - DBL_MANT_DIG);// = log(DBL_EPSILON) = -36.04..
        r = pp * (1. - qq) / (pp + 1.);

        t = 0.2;
        // FIXME: Factor 0.2 is a bit arbitrary;  '1' is clearly much too much.

        boolean skipToNewton = false;
        tx = u = xinbta = Double.NaN; // need to explicitly initialize these variables to avoid compilation problems.
        // however this value is actually guarateed to be overriden either in this if
        // coming or the next.

        if (MathConstants.LN_2 * DBL_MIN_EXP < u0 && // cannot allow exp(u0) = 0 ==> exp(u1) = exp(u0) = 0
                u0 < -0.01 && // (must: u0 < 0, but too close to 0 <==> x = exp(u0) = 0.99..)
                // qq <= 2 && // <--- "arbitrary"
                // u0 <  t*log_eps_c - log(fabs(r)) &&
                u0 < (t * log_eps_c - log(abs(pp * (1. - qq) * (2. - qq) / (2. * (pp + 2.))))) / 2.) {
// TODO: maybe jump here from below, when initial u "fails" ?
// L_tail_u:
            // MM's one-step correction (cheaper than 1 Newton!)
            r = r * exp(u0);// = r*x0
            if (r > -1.) {
                u = u0 - log1p(r) / pp;
            } else {
                u = u0;
            }
            tx = xinbta = exp(u);
            use_log_x = true; // or (u < log_q_cut)  ??
            skipToNewton = true;
        }
        if (!skipToNewton) {


            // y := y_\alpha in AS 64 := Hastings(1955) approximation of qnorm(1 - a) :
            r = sqrt(-2 * la);
            y = r - (const1 + const2 * r) / (1. + (const3 + const4 * r) * r);

            if (pp > 1 && qq > 1) { // use  Carter(1947), see AS 109, remark '5.'
                r = (y * y - 3.) / 6.;
                s = 1. / (pp + pp - 1.);
                t = 1. / (qq + qq - 1.);
                h = 2. / (s + t);
                w = y * sqrt(h + r) / h - (t - s) * (r + 5. / 6. - 2. / (3. * h));
                if (w > 300) { // exp(w+w) is huge or overflows
                    t = w + w + log(qq) - log(pp); // = argument of log1pexp(.)
                    u = // log(xinbta) = - log1p(qq/pp * exp(w+w)) = -log(1 + exp(t))
                            (t <= 18) ? -log1p(exp(t)) : -t - exp(-t);
                    xinbta = exp(u);
                } else {
                    xinbta = pp / (pp + qq * exp(w + w));
                    u = // log(xinbta)
                            -log1p(qq / pp * exp(w + w));
                }
            } else { // use the original AS 64 proposal, Scheffé-Tukey (1944) and Wilson-Hilferty
                r = qq + qq;
    /* A slightly more stable version of  t := \chi^2_{alpha} of AS 64
	 * t = 1. / (9. * qq); t = r * R_pow_di(1. - t + y * sqrt(t), 3);  */
                t = 1. / (3. * sqrt(qq));
                t = r * R_pow_di(1. + t * (-t + y), 3);// = \chi^2_{alpha} of AS 64
                s = 4. * pp + r - 2.;// 4p + 2q - 2 = numerator of new t = (...) / chi^2
                if (t == 0 || (t < 0. && s >= t)) { // cannot use chisq approx
                    // x0 = 1 - { (1-a)*q*B(p,q) } ^{1/q}    {AS 65}
                    // xinbta = 1. - exp((log(1-a)+ log(qq) + logbeta) / qq);
                    double l1ma;/* := log(1-a), directly from alpha (as 'la' above):
			 * FIXME: not worth it? log1p(-a) always the same ?? */
                    if (swap_tail)
                        l1ma = lower_tail ? R_D_log(log_p, alpha) : R_D_LExp(log_p, alpha);
                    else
                        l1ma = lower_tail ? R_D_LExp(log_p, alpha) : R_D_log(log_p, alpha);
                    double xx = (l1ma + log(qq) + logbeta) / qq;
                    if (xx <= 0.) {
                        xinbta = -expm1(xx);
                        u = R_Log1_Exp(xx);// =  log(xinbta) = log(1 - exp(...A...))
                    } else { // xx > 0 ==> 1 - e^xx < 0 .. is nonsense
                        xinbta = 0;
                        u = Double.NEGATIVE_INFINITY; /// FIXME can do better?
                    }
                } else {
                    t = s / t;
                    if (t <= 1.) { // cannot use chisq, either
                        u = (la + log(pp) + logbeta) / pp;
                        xinbta = exp(u);
                    } else { // (1+x0)/(1-x0) = t,  solved for x0 :
                        xinbta = 1. - 2. / (t + 1.);
                        u = log1p(-2. / (t + 1.));
                    }
                }
            }

            // Problem: If initial u is completely wrong, we make a wrong decision here
            if (swap_choose &&
                    ((swap_tail && u >= -exp(log_q_cut)) || // ==> "swap back"
                            (!swap_tail && u >= -exp(4 * log_q_cut) && pp / qq < 1000.))) { // ==> "swap now" (much less easily)
                // "revert swap" -- and use_log_x
                swap_tail = !swap_tail;
                if (swap_tail) {
                    a = R_DT_CIv(log_p, lower_tail, alpha); // needed ?
                    la = lower_tail ? R_D_LExp(log_p, alpha) : R_D_log(log_p, alpha);
                    pp = q;
                    qq = p;
                } else {
                    a = p_;
                    la = lower_tail ? R_D_log(log_p, alpha) : R_D_LExp(log_p, alpha);
                    pp = p;
                    qq = q;
                }
                // we could redo computations above, but this should be stable
                u = R_Log1_Exp(u);
                xinbta = exp(u);

/* Careful: "swap now"  should not fail if
   1) the above initial xinbta is "completely wrong"
   2) The correction step can go outside (u_n > 0 ==>  e^u > 1 is illegal)
   e.g., for
	qbeta(0.2066, 0.143891, 0.05)
*/
            }

            if (!use_log_x)
                use_log_x = (u < log_q_cut);//(per default) <==> xinbta = e^u < 4.54e-5
            boolean
                    bad_u = !Double.isFinite(u),
                    bad_init = bad_u || xinbta > p_hi;

            tx = xinbta; // keeping "original initial x" (for now)

            if (bad_u || u < log_q_cut) { /* e.g.
		    qbeta(0.21, .001, 0.05)
		    try "left border" quickly, i.e.,
		    try at smallest positive number: */
                w = CDF(DBL_very_MIN, pp, qq, true, log_p);
                if (w > (log_p ? la : a)) {
                    if (log_p || abs(w - a) < abs(0 - a)) { // DBL_very_MIN is better than 0
                        tx = DBL_very_MIN;
                        u_n = DBL_log_v_MIN;// = log(DBL_very_MIN)
                    } else {
                        tx = 0.;
                        u_n = Double.NEGATIVE_INFINITY;
                    }
                    use_log_x = log_p;
                    add_N_step = false;
                    return qbeta_raw_return(log_p, give_log_q, use_log_x, swap_tail, add_N_step, r, u_n, t,
                            a, logbeta, la, tx, pp, qq);
                } else {
                    if (u < DBL_log_v_MIN) {
                        u = DBL_log_v_MIN;// = log(DBL_very_MIN)
                        xinbta = DBL_very_MIN;
                    }
                }
            }


    /* Sometimes the approximation is negative (and == 0 is also not "ok") */
            if (bad_init && !(use_log_x && tx > 0)) {
                if (u == Double.NEGATIVE_INFINITY) {
                    u = MathConstants.LN_2 * DBL_MIN_EXP;
                    xinbta = DBL_MIN;
                } else {
                    xinbta = (xinbta > 1.1) // i.e. "way off"
                            ? 0.5 // otherwise, keep the respective boundary:
                            : ((xinbta < p_lo) ? exp(u) : p_hi);
                    if (bad_u)
                        u = log(xinbta);
                    // otherwise: not changing "potentially better" u than the above
                }
            }
        }
        boolean converged = false;
        //L_Newton:
    /* --------------------------------------------------------------------
     * Solve for x by a modified Newton-Raphson method, using CDF()
     */
        r = 1 - pp;
        t = 1 - qq;
        double wprev = 0., prev = 1., adj = 1.; // -Wall

        if (use_log_x) { // find  log(xinbta) -- work in  u := log(x) scale
            // if(bad_init && tx > 0) xinbta = tx;// may have been better

            loop0:
            for (i_pb = 0; i_pb < 1000; i_pb++) {
                // using log_p == TRUE  unconditionally here
                // FIXME: if exp(u) = xinbta underflows to 0, like different formula pbeta_log(u, *)
                y = CDF(xinbta, pp, qq, /*lower_tail = */ true, true);

	    /* w := Newton step size for   L(u) = log F(e^u)  =!= 0;   u := log(x)
	     *   =  (L(.) - la) / L'(.);  L'(u)= (F'(e^u) * e^u ) / F(e^u)
	     *   =  (L(.) - la)*F(.) / {F'(e^u) * e^u } =
	     *   =  (L(.) - la) * e^L(.) * e^{-log F'(e^u) - u}
	     *   =  ( y   - la) * e^{ y - u -log F'(e^u)}
	        and  -log F'(x)= -log f(x) =  + logbeta + (1-p) log(x) + (1-q) log(1-x)
		               = logbeta + (1-p) u + (1-q) log(1-e^u)
	     */
                w = (y == Double.NEGATIVE_INFINITY) // y = -Inf  well possible: we are on log scale!
                        ? 0. : (y - la) * exp(y - u + logbeta + r * u + t * R_Log1_Exp(u));
                if (!Double.isFinite(w))
                    break;
                if (i_pb >= n_N && w * wprev <= 0.)
                    prev = max(abs(adj), fpu);
                g = 1;
                for (i_inn = 0; i_inn < 1000; i_inn++) {
                    adj = g * w;
                    // take full Newton steps at the beginning; only then safe guard:
                    if (i_pb < n_N || abs(adj) < prev) {
                        u_n = u - adj; // u_{n+1} = u_n - g*w
                        if (u_n <= 0.) { // <==> 0 <  xinbta := e^u  <= 1
                            if (prev <= acu || abs(w) <= acu) {
			    /* R_ifDEBUG_printf(" -adj=%g, %s <= acu  ==> convergence\n", */
			    /* 	 -adj, (prev <= acu) ? "prev" : "|w|"); */
                                converged = true;
                                break loop0;
                            }
                            // if (u_n != ML_NEGINF && u_n != 1)
                            break;
                        }
                    }
                    g /= 3;
                }
                // (cancellation in (u_n -u) => may differ from adj:
                double D = min(abs(adj), abs(u_n - u));
	    /* R_ifDEBUG_printf(" delta(u)=%g\n", u_n - u); */
                if (D <= 4e-16 * abs(u_n + u)) {
                    converged = true;
                    break loop0;
                }
                u = u_n;
                xinbta = exp(u);
                wprev = w;
            } // for(i )

        } else {

            loop1:
            for (i_pb = 0; i_pb < 1000; i_pb++) {
                y = CDF(xinbta, pp, qq, /*lower_tail = */ true, log_p);
                // delta{y} :   d_y = y - (log_p ? la : a);
                if (!Double.isFinite(y) && !(log_p && y == Double.NEGATIVE_INFINITY))// y = -Inf  is ok if(log_p)
                { // ML_ERR_return_NAN :
                    throw new IllegalArgumentException("");
                }


	/* w := Newton step size  (F(.) - a) / F'(.)  or,
	 * --   log: (lF - la) / (F' / F) = exp(lF) * (lF - la) / F'
	 */
                w = log_p
                        ? (y - la) * exp(y + logbeta + r * log(xinbta) + t * log1p(-xinbta))
                        : (y - a) * exp(logbeta + r * log(xinbta) + t * log1p(-xinbta));
                if (i_pb >= n_N && w * wprev <= 0.)
                    prev = max(abs(adj), fpu);
                g = 1;
                for (i_inn = 0; i_inn < 1000; i_inn++) {
                    adj = g * w;
                    // take full Newton steps at the beginning; only then safe guard:
                    if (i_pb < n_N || abs(adj) < prev) {
                        tx = xinbta - adj; // x_{n+1} = x_n - g*w
                        if (0. <= tx && tx <= 1.) {
                            if (prev <= acu || abs(w) <= acu) {
                                converged = true;
                                break loop1;
                            }
                            if (tx != 0. && tx != 1)
                                break;
                        }
                    }
                    g /= 3;
                }
                if (abs(tx - xinbta) <= 4e-16 * (tx + xinbta)) { // "<=" : (.) == 0 {
                    converged = true;
                    break;
                }

                xinbta = tx;
                if (tx == 0) // "we have lost"
                    break;
                wprev = w;
            }
        }

    /*-- NOT converged: Iteration count --*/
        if (!converged) {
            throw new IllegalArgumentException("didn't converged");
        }

        log_ = log_p || use_log_x; // only for printing
        if ((log_ && y == Double.NEGATIVE_INFINITY) || (!log_ && y == 0)) {
            // stuck at left, try if smallest positive number is "better"
            w = CDF(DBL_very_MIN, pp, qq, true, log_);
            if (log_ || abs(w - a) <= abs(y - a)) {
                tx = DBL_very_MIN;
                u_n = DBL_log_v_MIN;// = log(DBL_very_MIN)
            }
            add_N_step = false; // not trying to do better anymore
        }
        return qbeta_raw_return(log_p, give_log_q, use_log_x, swap_tail, add_N_step, r, u_n, t,
                a, logbeta, la, tx, pp, qq);
    }

    private static double[] qbeta_raw_return(final boolean log_p, final boolean give_log_q,
                                             final boolean use_log_x,
                                             final boolean swap_tail, final boolean add_N_step,
                                             double r, final double u_n, final double t,
                                             final double a, final double logbeta, final double la,
                                             double tx, final double pp, final double qq) {
        if (give_log_q) { // ==> use_log_x , too
            if (!use_log_x) // (see if claim above is true)
                r = R_Log1_Exp(u_n);
            if (swap_tail) {
                return new double[]{r, u_n};
            } else {
                return new double[]{u_n, r};
            }
        } else {
            if (use_log_x) {
                if (add_N_step) {

		/* add one last Newton step on original x scale, e.g., for
		   qbeta(2^-98, 0.125, 2^-96) */
                    final double xinbta = exp(u_n);
                    final double y = CDF(xinbta, pp, qq, /*lower_tail = */ true, log_p);
                    final double w = log_p
                            ? (y - la) * exp(y + logbeta + r * log(xinbta) + t * log1p(-xinbta))
                            : (y - a) * exp(logbeta + r * log(xinbta) + t * log1p(-xinbta));
                    tx = xinbta - w;
                } else {
                    if (swap_tail) {
                        return new double[]{-expm1(u_n), exp(u_n)};
                    } else {
                        return new double[]{exp(u_n), -expm1(u_n)};
                    }
                }
            }
            if (swap_tail) {
                return new double[]{1 - tx, tx};
            } else {
                return new double[]{tx, 1 - tx};
            }
        }
    }

    public static double CDF(final double x, final double alpha, final double beta, final boolean log) {
        if (Double.isNaN(x) || Double.isNaN(alpha) || Double.isNaN(beta))
            return Double.NaN;
        else if (alpha < 0 || beta < 0)
            return Double.NaN;
        else if (x <= 0)
            return log ? Double.NEGATIVE_INFINITY : 0.0;
        else if (x >= 1)
            return log ? 0.0 : Double.NEGATIVE_INFINITY;
        return CDF(x, alpha, beta, true, log);
    }

    public static double logCDF(final double x, final double alpha, final double beta) {
        return CDF(x, alpha, beta, true);
    }

    public static double CDF(final double x, final double alpha, final double beta) {
        return CDF(x, alpha, beta, false);
    }

    static double CDF(double x, double a, double b, boolean lower_tail, boolean log_p) {
            // treat limit cases correctly here:
            if (a == 0 || b == 0 || !Double.isFinite(a) || !Double.isFinite(b)) {
                // NB:  0 < x < 1 :
                if (a == 0 && b == 0) // point mass 1/2 at each of {0,1} :
                    return (log_p ? -MathConstants.LN_2 : 0.5);
                if (a == 0 || a / b == 0) // point mass 1 at 0 ==> P(X <= x) = 1, all x > 0
                    return (log_p ? 0 : 1);
                if (b == 0 || b / a == 0) // point mass 1 at 1 ==> P(X <= x) = 0, all x < 1
                    return (log_p ? Double.NEGATIVE_INFINITY : 0);
                // else, remaining case:  a = b = Inf : point mass 1 at 1/2
                if (x < 0.5) return (log_p ? Double.NEGATIVE_INFINITY : 0);
                else return (log_p ? 0 : 1);
            }
            // Now:  0 < a < Inf;  0 < b < Inf

            double x1 = 0.5 - x + 0.5;
            //====
            final double[] w1 = bratio(a, b, x, x1, log_p); /* -> ./toms708.c */
            //====
            // ierr in {10,14} <==> bgrat() error code ierr-10 in 1:4; for 1 and 4, warned *there*
            return lower_tail ? w1[0] : w1[1];
    }


    static double[] bratio(double a, double b, double x, double y, boolean log_p) {
/* -----------------------------------------------------------------------
 *	      Evaluation of the Incomplete ExponentialUtils function I_x(a,b)
 *		       --------------------
 *     It is assumed that a and b are nonnegative, and that x <= 1
 *     and y = 1 - x.  Bratio assigns w and w1 the values
 *			w  = I_x(a,b)
 *			w1 = 1 - I_x(a,b)
 *     ierr is a variable that reports the status of the results.
 *     If no input errors are detected then ierr is set to 0 and
 *     w and w1 are computed. otherwise, if an error is detected,
 *     then w and w1 are assigned the value 0 and ierr is set to
 *     one of the following values ...
 *	  ierr = 1  if a or b is negative
 *	  ierr = 2  if a = b = 0
 *	  ierr = 3  if x < 0 or x > 1
 *	  ierr = 4  if y < 0 or y > 1
 *	  ierr = 5  if x + y != 1
 *	  ierr = 6  if x = a = 0
 *	  ierr = 7  if y = b = 0
 *	  ierr = 8	(not used currently)
 *	  ierr = 9  NaN in a, b, x, or y
 *	  ierr = 10     (not used currently)
 *	  ierr = 11  bgrat() error code 1 [+ warning in bgrat()]
 *	  ierr = 12  bgrat() error code 2   (no warning here)
 *	  ierr = 13  bgrat() error code 3   (no warning here)
 *	  ierr = 14  bgrat() error code 4 [+ WARNING in bgrat()]
 * --------------------
 *     Written by Alfred H. Morris, Jr.
 *	  Naval Surface Warfare Center
 *	  Dahlgren, Virginia
 *     Revised ... Nov 1991
* ----------------------------------------------------------------------- */

        boolean do_swap;
        int n, ierr1 = 0;
        double z, a0, b0, x0, y0, lambda;

/*  eps is a machine dependent constant: the smallest
 *      floating point number for which   1. + eps > 1.
 * NOTE: for almost all purposes it is replaced by 1e-15 (~= 4.5 times larger) below */
        double eps = DBL_EPSILON; /* == DBL_EPSILON (in R, Rmath) */

/* ----------------------------------------------------------------------- */
        double w = log_p ? Double.NEGATIVE_INFINITY : 0;
        double w1 = log_p ? Double.NEGATIVE_INFINITY : 0;

        // safeguard, preventing infinite loops further down
        if (Double.isNaN(x) || Double.isNaN(y) ||
                Double.isNaN(a) || Double.isNaN(b)) {
            throw new IllegalArgumentException();
        }
        if (a < 0. || b < 0.) {
            throw new IllegalArgumentException();
        }
        if (a == 0. && b == 0.) {
            throw new IllegalArgumentException();
        }
        if (x < 0. || x > 1.) {
            throw new IllegalArgumentException();
        }
        if (y < 0. || y > 1.) {
            throw new IllegalArgumentException();
        }

    /* check that  'y == 1 - x' : */
        z = x + y - 0.5 - 0.5;

        if (abs(z) > eps * 3.) {
            throw new IllegalArgumentException();
        }

        if (x == 0.) {
            w = log_p ? Double.NEGATIVE_INFINITY : 0;
            w1 = log_p ? 0 : 1;
            return bratio_return(false, w, w1);
        }
        if (y == 0.) {
            if (b == 0.) {
                throw new IllegalArgumentException();
            }
            w = log_p ? 0 : 1;
            w1 = log_p ? Double.NEGATIVE_INFINITY : 0;
            return bratio_return(false, w, w1);
        }

        if (a == 0.) {
            w = log_p ? 0 : 1;
            w1 = log_p ? Double.NEGATIVE_INFINITY : 0;
            return bratio_return(false, w, w1);
        }
        if (b == 0.) {
            w = log_p ? Double.NEGATIVE_INFINITY : 0;
            w1 = log_p ? 0 : 1;
            return bratio_return(false, w, w1);
        }

        eps = max(eps, 1e-15);
        boolean a_lt_b = (a < b);
        if (/* max(a,b) */ (a_lt_b ? b : a) < eps * .001) { /* procedure for a and b < 0.001 * eps */
            // L230:  -- result *independent* of x (!)
            // *w  = a/(a+b)  and  w1 = b/(a+b) :
            if (log_p) {
                if (a_lt_b) {
                    w = log1p(-a / (a + b)); // notably if a << b
                    w1 = log(a / (a + b));
                } else { // b <= a
                    w = log(b / (a + b));
                    w1 = log1p(-b / (a + b));
                }
            } else {
                w = b / (a + b);
                w1 = a / (a + b);
            }

            return new double[]{w, w1};
        }

        if (min(a, b) <= 1.) { /*------------------------ a <= 1  or  b <= 1 ---- */

            do_swap = (x > 0.5);
            if (do_swap) {
                a0 = b;
                x0 = y;
                b0 = a;
                y0 = x; // SET_0_swap.
            } else {
                a0 = a;
                x0 = x;
                b0 = b;
                y0 = y; // SET_0_noswap.
            }
	/* now have  x0 <= 1/2 <= y0  (still  x0+y0 == 1) */

            if (b0 < min(eps, eps * a0)) { /* L80: */
                w = fpser(a0, b0, x0, eps, log_p);
                w1 = log_p ? R_Log1_Exp(w) : 0.5 - w + 0.5;
                return bratio_return(do_swap, w, w1);
            }

            if (a0 < min(eps, eps * b0) && b0 * x0 <= 1.) { /* L90: */
                w1 = apser(a0, b0, x0, eps);
                return bratio_return_from_w1(log_p, do_swap, w1);
            }

            boolean did_bup = false;
            if (max(a0, b0) > 1.) { /* L20:  min(a,b) <= 1 < max(a,b)  */
                if (b0 <= 1.) return bratio_return_bpser_w(a0, b0, x0, eps, log_p, do_swap);

                if (x0 >= 0.29) {
                    return bratio_return_bpser_w1(b0, a0, y0, eps, log_p, do_swap);
       /* was 0.3, PR#13786 */
                }

                if (x0 < 0.1 && pow(x0 * b0, a0) <= 0.7)
                    return bratio_return_bpser_w(a0, b0, x0, eps, log_p, do_swap);

                if (b0 > 15.) {
                    w1 = 0.;
                    n = -1; // irrelevant, just added to make it compile.
                } else {
                    n = 20; /* goto L130; */
                    w1 = bup(b0, a0, y0, x0, n, eps, false);
                    did_bup = true;
                    b0 += n;
                }
            } else { /*  a, b <= 1 */
                if (a0 >= min(0.2, b0)) return bratio_return_bpser_w(a0, b0, x0, eps, log_p, do_swap);

                if (pow(x0, a0) <= 0.9) return bratio_return_bpser_w(a0, b0, x0, eps, log_p, do_swap);

                if (x0 >= 0.3) {
                    return bratio_return_bpser_w1(b0, a0, y0, eps, log_p, do_swap);
                }
                n = 20; /* goto L130; */
                w1 = bup(b0, a0, y0, x0, n, eps, false);
                did_bup = true;
                b0 += n;

            }
            w1 = bgrat(b0, a0, y0, x0, w1, 15 * eps, false);
            if (w1 == 0 || (0 < w1 && w1 < DBL_MIN)) { // w1=0 or very close:
                // "almost surely" from underflow, try more: [2013-03-04]
// FIXME: it is even better to do this in bgrat *directly* at least for the case
//  !did_bup, i.e., where *w1 = (0 or -Inf) on entry
                if (did_bup) { // re-do that part on log scale:
                    w1 = bup(b0 - n, a0, y0, x0, n, eps, true);
                } else {
                    w1 = Double.NEGATIVE_INFINITY; // = 0 on log-scale
                }
                w1 = bgrat(b0, a0, y0, x0, w1, 15 * eps, true);
                if (log_p) {
                    w = R_Log1_Exp(w1);
                } else {
                    w = /* 1 - exp(*w1) */ -expm1(w1);
                    w1 = exp(w1);
                }
                return bratio_return(do_swap, w, w1);
            }
            // else
            return bratio_return_from_w1(log_p, do_swap, w1);
        } else { /* L30: -------------------- both  a, b > 1  {a0 > 1  &  b0 > 1} ---*/

	/* lambda := a y - b x  =  (a + b)y  =  a - (a+b)x    {using x + y == 1},
	 * ------ using the numerically best version : */
            lambda = Double.isFinite(a + b)
                    ? ((a > b) ? (a + b) * y - b : a - (a + b) * x)
                    : a * y - b * x;
            do_swap = (lambda < 0.);
            if (do_swap) {
                lambda = -lambda;
                a0 = b;
                x0 = y;
                b0 = a;
                y0 = x; // SET_0_swap.
            } else {
                a0 = a;
                x0 = x;
                b0 = b;
                y0 = y; // SET_0_noswap.
            }

            if (b0 < 40.) {
                if (b0 * x0 <= 0.7
                        || (log_p && lambda > 650.)) // << added 2010-03; svn r51327
                    return bratio_return_bpser_w(a0, b0, x0, eps, log_p, do_swap);
            } else if (a0 > b0) { /* ----  a0 > b0 >= 40  ---- */
                if (b0 <= 100. || lambda > b0 * 0.03)
                    return bratio_return_bfrac(a0, b0, x0, y0, lambda, eps, log_p, do_swap);
                else {
                    w = basym(a0, b0, lambda, eps * 100., log_p);
                    w1 = log_p ? R_Log1_Exp(w) : 0.5 - w + 0.5;
                    return bratio_return(do_swap, w, w1);
                }

            } else if (a0 <= 100.) {
                return bratio_return_bfrac(a0, b0, x0, y0, lambda, eps, log_p, do_swap);
            } else if (lambda > a0 * 0.03) {
                return bratio_return_bfrac(a0, b0, x0, y0, lambda, eps, log_p, do_swap);
            } else {
                w = basym(a0, b0, lambda, eps * 100., log_p);
                w1 = log_p ? R_Log1_Exp(w) : 0.5 - w + 0.5;
                return bratio_return(do_swap, w, w1);
            }


        } /* else: a, b > 1 */

/*            EVALUATION OF THE APPROPRIATE ALGORITHM */

    /* b0 := fractional_part( b0 )  in (0, 1]  */
        n = (int) b0;
        b0 -= n;
        if (b0 == 0.) {
            --n;
            b0 = 1.;
        }

        w = bup(b0, a0, y0, x0, n, eps, false);

        if (w < DBL_MIN && log_p) { /* do not believe it; try bpser() : */
            //        R_ifDEBUG_printf(" L140: bup(b0=%g,..)=%.15g < DBL_MIN - not used; ", b0, * w);
	/*revert: */
            b0 += n;
	/* which is only valid if b0 <= 1 || b0*x0 <= 0.7 */
            return bratio_return_bpser_w(a0, b0, x0, eps, log_p, do_swap);
        }
//        R_ifDEBUG_printf(" L140: *w := bup(b0=%g,..) = %.15g; ", b0, * w);
        if (x0 <= 0.7) {
	/* log_p :  TODO:  w = bup(.) + bpser(.)  -- not so easy to use log-scale */
            w += bpser(a0, b0, x0, eps, /* log_p = */ false);
//            R_ifDEBUG_printf(" x0 <= 0.7: *w := *w + bpser(*) = %.15g\n", * w);
            return bratio_return_from_w(log_p, do_swap, w);
        }
    /* L150: */
        if (a0 <= 15.) {
            n = 20;
            w += bup(a0, b0, x0, y0, n, eps, false);
//            R_ifDEBUG_printf("\n a0 <= 15: *w := *w + bup(*) = %.15g;", * w);
            a0 += n;
        }
        w = bgrat(a0, b0, x0, y0, w, 15 * eps, false);
        return bratio_return_from_w(log_p, do_swap, w);

    } /* bratio */


    private static double[] bratio_return_bpser_w(final double a0, final double b0, final double x0, final double eps, final boolean log_p, final boolean do_swap) {
        final double w = bpser(a0, b0, x0, eps, log_p);
        final double w1 = log_p ? R_Log1_Exp(w) : 0.5 - w + 0.5;
        return bratio_return(do_swap, w, w1);
    }

    private static double[] bratio_return_bpser_w1(final double b0, final double a0, final double y0, final double eps, final boolean log_p, final boolean do_swap) {
        final double w1 = bpser(b0, a0, y0, eps, log_p);
        final double w = log_p ? R_Log1_Exp(w1) : 0.5 - w1 + 0.5;
        return bratio_return(do_swap, w, w1);
    }

    private static double[] bratio_return_bfrac(final double a0, final double b0, final double x0, final double y0, final double lambda, final double eps, final boolean log_p, final boolean do_swap) {
        final double w = bfrac(a0, b0, x0, y0, lambda, eps * 15., log_p);
        final double w1 = log_p ? R_Log1_Exp(w) : 0.5 - w + 0.5;
        return bratio_return(do_swap, w, w1);
    }

    private static double[] bratio_return_from_w1(final boolean log_p, final boolean do_swap, final double w1) {
        double actualW;
        double actualW1;
        if (log_p) {
            actualW = log1p(-w1);
            actualW1 = log(w1);
        } else {
            actualW = 0.5 - w1 + 0.5;
            actualW1 = w1;
        }
        return do_swap ? new double[]{actualW1, actualW} : new double[]{actualW, actualW1};
    }

    private static double[] bratio_return_from_w(final boolean log_p, final boolean do_swap, final double w) {

        double actualW;
        double actualW1;
        if (log_p) {
            actualW1 = log1p(-w);
            actualW = log(w);
        } else {
            actualW1 = 0.5 - w + 0.5;
            actualW = w;
        }
        return do_swap ? new double[]{actualW1, actualW} : new double[]{actualW, actualW1};
    }

    private static double[] bratio_return(final boolean do_swap, final double w, final double w1) {
        return do_swap ? new double[]{w1, w} : new double[]{w, w1};
    }


    static double fpser(double a, double b, double x, double eps, boolean log_p) {
/* ----------------------------------------------------------------------- *
 *                 EVALUATION OF I (A,B)
 *                                X
 *          FOR B < MIN(EPS, EPS*A) AND X <= 0.5
 * ----------------------------------------------------------------------- */

        double ans, c, s, t, an, tol;

    /* SET  ans := x^a : */
        if (log_p) {
            ans = a * log(x);
        } else if (a > eps * 0.001) {
            t = a * log(x);
            if (t < exparg(1)) { /* exp(t) would underflow */
                return 0.;
            }
            ans = exp(t);
        } else
            ans = 1.;

/*                NOTE THAT 1/B(A,B) = B */

        if (log_p)
            ans += log(b) - log(a);
        else
            ans *= b / a;

        tol = eps / a;
        an = a + 1.;
        t = x;
        s = t / an;
        do {
            an += 1.;
            t = x * t;
            c = t / an;
            s += c;
        } while (abs(c) > tol);

        if (log_p)
            ans += log1p(a * s);
        else
            ans *= a * s + 1.;
        return ans;
    }

    static double bup(double a, double b, double x, double y, int n, double eps,
                      boolean give_log) {
/* ----------------------------------------------------------------------- */
/*     EVALUATION OF I_x(A,B) - I_x(A+N,B) WHERE N IS A POSITIVE INT. */
/*     EPS IS THE TOLERANCE USED. */
/* ----------------------------------------------------------------------- */

        double ret_val;
        int i, k, mu;
        double d, l;

// Obtain the scaling factor exp(-mu) and exp(mu)*(x^a * y^b / beta(a,b))/a

        double apb = a + b,
                ap1 = a + 1.;
        if (n > 1 && a >= 1. && apb >= ap1 * 1.1) {
            mu = (int) abs(exparg(1));
            k = (int) exparg(0);
            if (mu > k)
                mu = k;
            d = exp(-(double) mu);
        } else {
            mu = 0;
            d = 1.;
        }

    /* L10: */
        ret_val = give_log
                ? brcmp1(mu, a, b, x, y, true) - log(a)
                : brcmp1(mu, a, b, x, y, false) / a;
        if (n == 1 ||
                (give_log && ret_val == Double.NEGATIVE_INFINITY) || (!give_log && ret_val == 0.))
            return ret_val;

        int nm1 = n - 1;
        double w = d;

/*          LET K BE THE INDEX OF THE MAXIMUM TERM */

        k = 0;
        if (b > 1.) {
            if (y > 1e-4) {
                double r = (b - 1.) * x / y - a;
                if (r >= 1.)
                    k = (r < nm1) ? (int) r : nm1;
            } else
                k = nm1;

//          ADD THE INCREASING TERMS OF THE SERIES - if k > 0
/* L30: */
            for (i = 0; i < k; ++i) {
                l = (double) i;
                d *= (apb + l) / (ap1 + l) * x;
                w += d;
            }
        }

// L40:     ADD THE REMAINING TERMS OF THE SERIES

        for (i = k; i < nm1; ++i) {
            l = (double) i;
            d *= (apb + l) / (ap1 + l) * x;
            w += d;
            if (d <= eps * w) /* relativ convergence (eps) */
                break;
        }

        // L50: TERMINATE THE PROCEDURE
        if (give_log) {
            ret_val += log(w);
        } else
            ret_val *= w;

        return ret_val;
    }

    static double exparg(int l) {
/* --------------------------------------------------------------------
 *     If l = 0 then  exparg(l) = The largest positive W for which
 *     exp(W) can be computed. With 0.99999 fuzz  ==> exparg(0) =   709.7756  nowadays
 *     if l = 1 (nonzero) then  exparg(l) = the largest negative W for
 *     which the computed value of exp(W) is nonzero.
 *     With 0.99999 fuzz			  ==> exparg(1) =  -709.0825  nowadays
 *     Note... only an approximate value for exparg(L) is needed.
 * -------------------------------------------------------------------- */

        final double lnb = .69314718055995;
        int m = (l == 0) ? DBL_MAX_EXP : DBL_MIN_EXP - 1;

        return m * lnb * .99999;
    } /* exparg */

    static double brcmp1(final int mu, final double a, final double b, final double x, final double y, final boolean give_log) {
/* -----------------------------------------------------------------------
 *          Evaluation of    exp(mu) * x^a * y^b / beta(a,b)
 * ----------------------------------------------------------------------- */

        final double const__ = .398942280401433; /* == 1/sqrt(2*pi); */
    /* R has  M_1_SQRT_2PI */

    /* Local variables */
        double c, t, u, v, z, a0, b0, apb;

        a0 = min(a, b);
        if (a0 < 8.) {
            double lnx, lny;
            if (x <= .375) {
                lnx = log(x);
                lny = alnrel(-x);
            } else if (y > .375) {
                // L11:
                lnx = log(x);
                lny = log(y);
            } else {
                lnx = alnrel(-y);
                lny = log(y);
            }

            // L20:
            z = a * lnx + b * lny;
            if (a0 >= 1.) {
                z -= betaln(a, b);
                return esum(mu, z, give_log);
            }
            // else :
	/* ----------------------------------------------------------------------- */
	/*              PROCEDURE FOR A < 1 OR B < 1 */
	/* ----------------------------------------------------------------------- */
            // L30:
            b0 = max(a, b);
            if (b0 >= 8.) {
	/* L80:                  ALGORITHM FOR b0 >= 8 */
                u = gamln1(a0) + algdiv(a0, b0);
                return give_log
                        ? log(a0) + esum(mu, z - u, true)
                        : a0 * esum(mu, z - u, false);

            } else if (b0 <= 1.) {
                //                   a0 < 1, b0 <= 1
                double ans = esum(mu, z, give_log);
                if (ans == (give_log ? Double.NEGATIVE_INFINITY : 0.))
                    return ans;

                apb = a + b;
                if (apb > 1.) {
                    // L40:
                    u = a + b - 1.;
                    z = (gam1(u) + 1.) / apb;
                } else {
                    z = gam1(apb) + 1.;
                }
                // L50:
                c = give_log
                        ? log1p(gam1(a)) + log1p(gam1(b)) - log(z)
                        : (gam1(a) + 1.) * (gam1(b) + 1.) / z;
                return give_log
                        ? ans + log(a0) + c - log1p(a0 / b0)
                        : ans * (a0 * c) / (a0 / b0 + 1.);
            }
            // else:               algorithm for	a0 < 1 < b0 < 8
            // L60:
            u = gamln1(a0);
            int n = (int) (b0 - 1.);
            if (n >= 1) {
                c = 1.;
                for (int i = 1; i <= n; ++i) {
                    b0 += -1.;
                    c *= b0 / (a0 + b0);
		/* L61: */
                }
                u += log(c); // TODO?: log(c) = log( prod(...) ) =  sum( log(...) )
            }
            // L70:
            z -= u;
            b0 += -1.;
            apb = a0 + b0;
            if (apb > 1.) {
                // L71:
                t = (gam1(apb - 1.) + 1.) / apb;
            } else {
                t = gam1(apb) + 1.;
            }
            // L72:
            return give_log
                    ? log(a0) + esum(mu, z, true) + log1p(gam1(b0)) - log(t) // TODO? log(t) = log1p(..)
                    : a0 * esum(mu, z, false) * (gam1(b0) + 1.) / t;

        } else {

/* ----------------------------------------------------------------------- */
/*              PROCEDURE FOR A >= 8 AND B >= 8 */
/* ----------------------------------------------------------------------- */
            // L100:
            double h, x0, y0, lambda;
            if (a > b) {
                // L101:
                h = b / a;
                x0 = 1. / (h + 1.);// => lx0 := log(x0) = 0 - log1p(h)
                y0 = h / (h + 1.);
                lambda = (a + b) * y - b;
            } else {
                h = a / b;
                x0 = h / (h + 1.);  // => lx0 := log(x0) = - log1p(1/h)
                y0 = 1. / (h + 1.);
                lambda = a - (a + b) * x;
            }
            double lx0 = -log1p(b / a); // in both cases

            // L110:
            double e = -lambda / a;
            if (abs(e) > 0.6) {
                // L111:
                u = e - log(x / x0);
            } else {
                u = rlog1(e);
            }

            // L120:
            e = lambda / b;
            if (abs(e) > 0.6) {
                // L121:
                v = e - log(y / y0);
            } else {
                v = rlog1(e);
            }

            // L130:
            z = esum(mu, -(a * u + b * v), give_log);
            return give_log
                    ? log(const__) + (log(b) + lx0) / 2. + z - bcorr(a, b)
                    : const__ * sqrt(b * x0) * z * exp(-bcorr(a, b));
        }

    }

    static double alnrel(final double a) {
/* -----------------------------------------------------------------------
 *            Evaluation of the function ln(1 + a)
 * ----------------------------------------------------------------------- */

        if (abs(a) > 0.375)
            return log(1. + a);
        // else : |a| <= 0.375
        final double
                p1 = -1.29418923021993,
                p2 = .405303492862024,
                p3 = -.0178874546012214,
                q1 = -1.62752256355323,
                q2 = .747811014037616,
                q3 = -.0845104217945565;
        final double
                t = a / (a + 2.),
                t2 = t * t,
                w = (((p3 * t2 + p2) * t2 + p1) * t2 + 1.) /
                        (((q3 * t2 + q2) * t2 + q1) * t2 + 1.);
        return t * 2. * w;

    }

    static double betaln(double a0, double b0) {
/* -----------------------------------------------------------------------
 *     Evaluation of the logarithm of the beta function  ln(beta(a0,b0))
 * ----------------------------------------------------------------------- */

        final double e = .918938533204673;/* e == 0.5*LN(2*PI) */

        double
                a = min(a0, b0),
                b = max(a0, b0);

        if (a < 8.) {
            if (a < 1.) {
/* ----------------------------------------------------------------------- */
//                    		A < 1
/* ----------------------------------------------------------------------- */
                if (b < 8.)
                    return gamln(a) + (gamln(b) - gamln(a + b));
                else
                    return gamln(a) + algdiv(a, b);
            }
	/* else */
/* ----------------------------------------------------------------------- */
//				1 <= A < 8
/* ----------------------------------------------------------------------- */
            double w;
            if (a < 2.) {
                if (b <= 2.) {
                    return gamln(a) + gamln(b) - gsumln(a, b);
                }
	    /* else */

                if (b < 8.) {
                    w = 0.;
                    int n = (int) (b - 1.);
                    double z = 1.;
                    for (int i = 1; i <= n; ++i) {
                        b += -1.;
                        z *= b / (a + b);
                    }
                    return w + log(z) + (gamln(a) + (gamln(b) - gsumln(a, b)));
                }
                return gamln(a) + algdiv(a, b);
            }
            // else L30:    REDUCTION OF A WHEN B <= 1000

            if (b <= 1e3) {
                int n = (int) (a - 1.);
                w = 1.;
                for (int i = 1; i <= n; ++i) {
                    a += -1.;
                    double h = a / b;
                    w *= h / (h + 1.);
                }
                w = log(w);

                if (b >= 8.)
                    return w + gamln(a) + algdiv(a, b);

                // else
                L40:
                // 	1 < A <= B < 8 :  reduction of B
                n = (int) (b - 1.);
                double z = 1.;
                for (int i = 1; i <= n; ++i) {
                    b += -1.;
                    z *= b / (a + b);
                }
                return w + log(z) + (gamln(a) + (gamln(b) - gsumln(a, b)));
            } else { // L50:	reduction of A when  B > 1000
                int n = (int) (a - 1.);
                w = 1.;
                for (int i = 1; i <= n; ++i) {
                    a += -1.;
                    w *= a / (a / b + 1.);
                }
                return log(w) - n * log(b) + (gamln(a) + algdiv(a, b));
            }

        } else {
/* ----------------------------------------------------------------------- */
            // L60:			A >= 8
/* ----------------------------------------------------------------------- */

            double
                    w = bcorr(a, b),
                    h = a / b,
                    u = -(a - 0.5) * log(h / (h + 1.)),
                    v = b * alnrel(h);
            if (u > v)
                return log(b) * -0.5 + e + w - v - u;
            else
                return log(b) * -0.5 + e + w - u - v;
        }

    }

    static double gamln(double a) {
/* -----------------------------------------------------------------------
 *            Evaluation of  ln(gamma(a))  for positive a
 * ----------------------------------------------------------------------- */
/*     Written by Alfred H. Morris */
/*          Naval Surface Warfare Center */
/*          Dahlgren, Virginia */
/* ----------------------------------------------------------------------- */

        final double d = .418938533204673;/* d == 0.5*(LN(2*PI) - 1) */

        final double c0 = .0833333333333333;
        final double c1 = -.00277777777760991;
        final double c2 = 7.9365066682539e-4;
        final double c3 = -5.9520293135187e-4;
        final double c4 = 8.37308034031215e-4;
        final double c5 = -.00165322962780713;

        if (a <= 0.8)
            return gamln1(a) - log(a); /* ln(G(a+1)) - ln(a) == ln(G(a+1)/a) = ln(G(a)) */
        else if (a <= 2.25)
            return gamln1(a - 0.5 - 0.5);

        else if (a < 10.) {
            int i, n = (int) (a - 1.25);
            double t = a;
            double w = 1.;
            for (i = 1; i <= n; ++i) {
                t += -1.;
                w *= t;
            }
            return gamln1(t - 1.) + log(w);
        } else { /* a >= 10 */
            double t = 1. / (a * a);
            double w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a;
            return d + w + (a - 0.5) * (log(a) - 1.);
        }
    }

    private static double gamln1(double a) {
/* ----------------------------------------------------------------------- */
/*     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 <= A <= 1.25 */
/* ----------------------------------------------------------------------- */

        double w;
        if (a < 0.6) {
            final double p0 = .577215664901533;
            final double p1 = .844203922187225;
            final double p2 = -.168860593646662;
            final double p3 = -.780427615533591;
            final double p4 = -.402055799310489;
            final double p5 = -.0673562214325671;
            final double p6 = -.00271935708322958;
            final double q1 = 2.88743195473681;
            final double q2 = 3.12755088914843;
            final double q3 = 1.56875193295039;
            final double q4 = .361951990101499;
            final double q5 = .0325038868253937;
            final double q6 = 6.67465618796164e-4;
            w = ((((((p6 * a + p5) * a + p4) * a + p3) * a + p2) * a + p1) * a + p0) /
                    ((((((q6 * a + q5) * a + q4) * a + q3) * a + q2) * a + q1) * a + 1.);
            return -(a) * w;
        } else { /* 0.6 <= a <= 1.25 */
            final double r0 = .422784335098467;
            final double r1 = .848044614534529;
            final double r2 = .565221050691933;
            final double r3 = .156513060486551;
            final double r4 = .017050248402265;
            final double r5 = 4.97958207639485e-4;
            final double s1 = 1.24313399877507;
            final double s2 = .548042109832463;
            final double s3 = .10155218743983;
            final double s4 = .00713309612391;
            final double s5 = 1.16165475989616e-4;
            double x = a - 0.5 - 0.5;
            w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) /
                    (((((s5 * x + s4) * x + s3) * x + s2) * x + s1) * x + 1.);
            return x * w;
        }
    }

    private static double algdiv(final double a, final double b) {
/* ----------------------------------------------------------------------- */

/*     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B >= 8 */

/*                         -------- */

/*     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY */
/*     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X). */

/* ----------------------------------------------------------------------- */

    /* Initialized data */

        final double c0 = .0833333333333333;
        final double c1 = -.00277777777760991;
        final double c2 = 7.9365066682539e-4;
        final double c3 = -5.9520293135187e-4;
        final double c4 = 8.37308034031215e-4;
        final double c5 = -.00165322962780713;

        double c, d, h, t, u, v, w, x, s3, s5, x2, s7, s9, s11;

/* ------------------------ */
        if (a > b) {
            h = b / a;
            c = 1. / (h + 1.);
            x = h / (h + 1.);
            d = a + (b - 0.5);
        } else {
            h = a / b;
            c = h / (h + 1.);
            x = 1. / (h + 1.);
            d = b + (a - 0.5);
        }

/* Set s<n> = (1 - x^n)/(1 - x) : */

        x2 = x * x;
        s3 = x + x2 + 1.;
        s5 = x + x2 * s3 + 1.;
        s7 = x + x2 * s5 + 1.;
        s9 = x + x2 * s7 + 1.;
        s11 = x + x2 * s9 + 1.;

/* w := Del(b) - Del(a + b) */

        t = 1. / (b * b);
        w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 *
                s3) * t + c0;
        w *= c / b;

/*                    COMBINE THE RESULTS */

        u = d * alnrel(a / b);
        v = a * (log(b) - 1.);
        if (u > v)
            return w - v - u;
        else
            return w - u - v;
    }

    private static double gsumln(final double a, final double b) {
/* ----------------------------------------------------------------------- */
/*          EVALUATION OF THE FUNCTION LN(GAMMA(A + B)) */
/*          FOR 1 <= A <= 2  AND  1 <= B <= 2 */
/* ----------------------------------------------------------------------- */

        final double x = a + b - 2.;/* in [0, 2] */

        if (x <= 0.25)
            return gamln1(x + 1.);

    /* else */
        if (x <= 1.25)
            return gamln1(x) + alnrel(x);
    /* else x > 1.25 : */
        return gamln1(x - 1.) + log(x * (x + 1.));

    }

    private static double bcorr(double a0, double b0) {
/* ----------------------------------------------------------------------- */

/*     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE */
/*     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A). */
/*     IT IS ASSUMED THAT A0 >= 8 AND B0 >= 8. */

/* ----------------------------------------------------------------------- */
    /* Initialized data */

        final double c0 = .0833333333333333;
        final double c1 = -.00277777777760991;
        final double c2 = 7.9365066682539e-4;
        final double c3 = -5.9520293135187e-4;
        final double c4 = 8.37308034031215e-4;
        final double c5 = -.00165322962780713;

    /* System generated locals */
        double ret_val, r1;

    /* Local variables */
        double a, b, c, h, t, w, x, s3, s5, x2, s7, s9, s11;
/* ------------------------ */
        a = min(a0, b0);
        b = max(a0, b0);

        h = a / b;
        c = h / (h + 1.);
        x = 1. / (h + 1.);
        x2 = x * x;

/*                SET SN = (1 - X^N)/(1 - X) */

        s3 = x + x2 + 1.;
        s5 = x + x2 * s3 + 1.;
        s7 = x + x2 * s5 + 1.;
        s9 = x + x2 * s7 + 1.;
        s11 = x + x2 * s9 + 1.;

/*                SET W = DEL(B) - DEL(A + B) */

/* Computing 2nd power */
        r1 = 1. / b;
        t = r1 * r1;
        w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 *
                s3) * t + c0;
        w *= c / b;

/*                   COMPUTE  DEL(A) + W */

/* Computing 2nd power */
        r1 = 1. / a;
        t = r1 * r1;
        ret_val = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a +
                w;
        return ret_val;
    }

    private static double rlog1(final double x) {
/* -----------------------------------------------------------------------
 *             Evaluation of the function  x - ln(1 + x)
 * ----------------------------------------------------------------------- */

        final double a = .0566749439387324;
        final double b = .0456512608815524;
        final double p0 = .333333333333333;
        final double p1 = -.224696413112536;
        final double p2 = .00620886815375787;
        final double q1 = -1.27408923933623;
        final double q2 = .354508718369557;

        double h, r, t, w, w1;
        if (x < -0.39 || x > 0.57) { /* direct evaluation */
            w = x + 0.5 + 0.5;
            return x - log(w);
        }
    /* else */
        if (x < -0.18) { /* L10: */
            h = x + .3;
            h /= .7;
            w1 = a - h * .3;
        } else if (x > 0.18) { /* L20: */
            h = x * .75 - .25;
            w1 = b + h / 3.;
        } else { /*		Argument Reduction */
            h = x;
            w1 = 0.;
        }

/* L30:              	Series Expansion */

        r = h / (h + 2.);
        t = r * r;
        w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.);
        return t * 2. * (1. / (1. - r) - r * w) + w1;

    } /* rlog1 */

    private static double esum(final int mu, final double x, boolean give_log) {
/* ----------------------------------------------------------------------- */
/*                    EVALUATION OF EXP(MU + X) */
/* ----------------------------------------------------------------------- */

        if (give_log)
            return x + (double) mu;

        // else :
        final double w;
        if (x > 0.) { /* L10: */
            if (mu > 0) return exp((double) mu) * exp(x);
            w = mu + x;
            if (w < 0.) return exp((double) mu) * exp(x);
        } else { /* x <= 0 */
            if (mu < 0) return exp((double) mu) * exp(x);
            w = mu + x;
            if (w > 0.) return exp((double) mu) * exp(x);
        }
        return exp(w);

    }

    private static double gam1(final double a) {
/*     ------------------------------------------------------------------ */
/*     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 <= A <= 1.5 */
/*     ------------------------------------------------------------------ */

        double d, t, w, bot, top;

        t = a;
        d = a - 0.5;
        // t := if(a > 1/2)  a-1  else  a
        if (d > 0.)
            t = d - 0.5;
        if (t < 0.) { /* L30: */
            final double[]
                    r = new double[]{-.422784335098468, -.771330383816272,
                    -.244757765222226, .118378989872749, 9.30357293360349e-4,
                    -.0118290993445146, .00223047661158249, 2.66505979058923e-4,
                    -1.32674909766242e-4};

            final double s1 = .273076135303957;
            final double s2 = .0559398236957378;

            top = (((((((r[8] * t + r[7]) * t + r[6]) * t + r[5]) * t + r[4]
            ) * t + r[3]) * t + r[2]) * t + r[1]) * t + r[0];
            bot = (s2 * t + s1) * t + 1.;
            w = top / bot;
            if (d > 0.)
                return t * w / a;
            else
                return a * (w + 0.5 + 0.5);

        } else if (t == 0) { // L10: a in {0, 1}
            return 0.;

        } else { /* t > 0;  L20: */
            final double
                    p[] = {.577215664901533, -.409078193005776,
                    -.230975380857675, .0597275330452234, .0076696818164949,
                    -.00514889771323592, 5.89597428611429e-4},
                    q[] = {1., .427569613095214, .158451672430138,
                            .0261132021441447, .00423244297896961};

            top = (((((p[6] * t + p[5]) * t + p[4]) * t + p[3]) * t + p[2]
            ) * t + p[1]) * t + p[0];
            bot = (((q[4] * t + q[3]) * t + q[2]) * t + q[1]) * t + 1.;
            w = top / bot;
            if (d > 0.) /* L21: */
                return t / a * (w - 0.5 - 0.5);
            else
                return a * w;
        }
    }

    private static double bpser(double a, double b, double x, double eps, final boolean log_p) {
/* -----------------------------------------------------------------------
 * Power SERies expansion for evaluating I_x(a,b) when
 *	       b <= 1 or b*x <= 0.7.   eps is the tolerance used.
 * NB: if log_p is TRUE, also use it if   (b < 40  & lambda > 650)
 * ----------------------------------------------------------------------- */

        int i, m;
        double ans, c, t, u, z, a0, b0, apb;

        if (x == 0.) {
            return log_p ? Double.NEGATIVE_INFINITY : 0;
        }
/* ----------------------------------------------------------------------- */
/*	      compute the factor  x^a/(a*ExponentialUtils(a,b)) */
/* ----------------------------------------------------------------------- */
        a0 = min(a, b);
        if (a0 >= 1.) { /*		 ------	 1 <= a0 <= b0  ------ */
            z = a * log(x) - betaln(a, b);
            ans = log_p ? z - log(a) : exp(z) / a;
        } else {
            b0 = max(a, b);

            if (b0 < 8.) {

                if (b0 <= 1.) { /*	 ------	 a0 < 1	 and  b0 <= 1  ------ */

                    if (log_p) {
                        ans = a * log(x);
                    } else {
                        ans = pow(x, a);
                        if (ans == 0.) /* once underflow, always underflow .. */
                            return ans;
                    }
                    apb = a + b;
                    if (apb > 1.) {
                        u = a + b - 1.;
                        z = (gam1(u) + 1.) / apb;
                    } else {
                        z = gam1(apb) + 1.;
                    }
                    c = (gam1(a) + 1.) * (gam1(b) + 1.) / z;

                    if (log_p) /* FIXME ? -- improve quite a bit for c ~= 1 */
                        ans += log(c * (b / apb));
                    else
                        ans *= c * (b / apb);

                } else { /* 	------	a0 < 1 < b0 < 8	 ------ */

                    u = gamln1(a0);
                    m = (int) (b0 - 1.);
                    if (m >= 1) {
                        c = 1.;
                        for (i = 1; i <= m; ++i) {
                            b0 += -1.;
                            c *= b0 / (a0 + b0);
                        }
                        u += log(c);
                    }

                    z = a * log(x) - u;
                    b0 += -1.; // => b0 in (0, 7)
                    apb = a0 + b0;
                    if (apb > 1.) {
                        u = a0 + b0 - 1.;
                        t = (gam1(u) + 1.) / apb;
                    } else {
                        t = gam1(apb) + 1.;
                    }

                    if (log_p) /* FIXME? potential for improving log(t) */
                        ans = z + log(a0 / a) + log1p(gam1(b0)) - log(t);
                    else
                        ans = exp(z) * (a0 / a) * (gam1(b0) + 1.) / t;
                }

            } else { /* 		------  a0 < 1 < 8 <= b0  ------ */

                u = gamln1(a0) + algdiv(a0, b0);
                z = a * log(x) - u;

                if (log_p)
                    ans = z + log(a0 / a);
                else
                    ans = a0 / a * exp(z);
            }
        }
        if (ans == (log_p ? Double.NEGATIVE_INFINITY : 0) || (!log_p && a <= eps * 0.1)) {
            return ans;
        }

/* ----------------------------------------------------------------------- */
/*		       COMPUTE THE SERIES */
/* ----------------------------------------------------------------------- */
        double tol = eps / a,
                n = 0.,
                sum = 0., w;
        c = 1.;
        do { // sum is alternating as long as n < b (<==> 1 - b/n < 0)
            n += 1.;
            c *= (0.5 - b / n + 0.5) * x;
            w = c / (a + n);
            sum += w;
        } while (n < 1e7 && abs(w) > tol);
//        if(abs(w) > tol) { // the series did not converge (in time)
//            // warn only when the result seems to matter:
//            if(( log_p && !(a*sum > -1. && abs(log1p(a * sum)) < eps*abs(ans))) ||
//                    (!log_p && abs(a*sum + 1.) != 1.))
//                MATHLIB_WARNING5(
//                        " bpser(a=%g, b=%g, x=%g,...) did not converge (n=1e7, |w|/tol=%g > 1; A=%g)",
//                        a,b,x, abs(w)/tol, ans);//
//        }
//        R_ifDEBUG_printf("  -> n=%.0f iterations, |w|=%g %s %g=tol:=eps/a ==> a*sum=%g\n",
//                n, fabs(w), (fabs(w) > tol) ? ">!!>" : "<=",
//                tol, a*sum);
        if (log_p) {
            if (a * sum > -1.) ans += log1p(a * sum);
            else {
//                if(ans > Double.NEGATIVE_INFINITY)
//                    MATHLIB_WARNING3(
//                            "pbeta(*, log.p=TRUE) -> bpser(a=%g, b=%g, x=%g,...) underflow to -Inf",
//                            a,b,x);
                ans = Double.NEGATIVE_INFINITY;
            }
        } else if (a * sum > -1.)
            ans *= (a * sum + 1.);
        else // underflow to
            ans = 0.;
        return ans;
    }

    private static double bfrac(double a, double b, double x, double y, double lambda,
                                double eps, boolean log_p) {
/* -----------------------------------------------------------------------
       Continued fraction expansion for I_x(a,b) when a, b > 1.
       It is assumed that  lambda = (a + b)*y - b.
   -----------------------------------------------------------------------*/

        double c, e, n, p, r, s, t, w, c0, c1, r0, an, bn, yp1, anp1, bnp1,
                beta, alpha, brc;

        if (!Double.isFinite(lambda)) return Double.NaN;// TODO: can return 0 or 1 (?)
        brc = brcomp(a, b, x, y, log_p);
        if (Double.isNaN(brc)) { // e.g. from   L <- 1e308; pnbinom(L, L, mu = 5)
            return Double.NaN; // TODO: could we know better?
        }
        if (!log_p && brc == 0.) {
            return 0.;
        }

        c = lambda + 1.;
        c0 = b / a;
        c1 = 1. / a + 1.;
        yp1 = y + 1.;

        n = 0.;
        p = 1.;
        s = a + 1.;
        an = 0.;
        bn = 1.;
        anp1 = 1.;
        bnp1 = c / c1;
        r = c1 / c;

/*        CONTINUED FRACTION CALCULATION */

        do {
            n += 1.;
            t = n / a;
            w = n * (b - n) * x;
            e = a / s;
            alpha = p * (p + c0) * e * e * (w * x);
            e = (t + 1.) / (c1 + t + t);
            beta = n + w / s + e * (c + n * yp1);
            p = t + 1.;
            s += 2.;

	/* update an, bn, anp1, and bnp1 */

            t = alpha * an + beta * anp1;
            an = anp1;
            anp1 = t;
            t = alpha * bn + beta * bnp1;
            bn = bnp1;
            bnp1 = t;

            r0 = r;
            r = anp1 / bnp1;
            if (abs(r - r0) <= eps * r)
                break;

	/* rescale an, bn, anp1, and bnp1 */

            an /= bnp1;
            bn /= bnp1;
            anp1 = r;
            bnp1 = 1.;
        } while (n < 10000);// arbitrary; had '1' --> infinite loop for  lambda = Inf
        return (log_p ? brc + log(r) : brc * r);
    } /* bfrac */

    private static double bgrat(double a, double b, double x, double y, double w,
                                double eps, boolean log_w) {
/* -----------------------------------------------------------------------
*     Asymptotic Expansion for I_x(a,b)  when a is larger than b.
*     Compute   w := w + I_x(a,b)
*     It is assumed a >= 15 and b <= 1.
*     eps is the tolerance used.
*     ierr is a variable that reports the status of the results.
*
* if(log_w),  *w  itself must be in log-space;
*     compute   w := w + I_x(a,b)  but return *w = log(w):
*          *w := log(exp(*w) + I_x(a,b)) = logspace_add(*w, log( I_x(a,b) ))
* ----------------------------------------------------------------------- */

        final int n_terms_bgrat = 30;
        double[] c = new double[n_terms_bgrat], d = new double[n_terms_bgrat];
        double bm1 = b - 0.5 - 0.5,
                nu = a + bm1 * 0.5, /* nu = a + (b-1)/2 =: T, in (9.1) of
			     * Didonato & Morris(1992), p.362 */
                lnx = (y > 0.375) ? log(x) : alnrel(-y),
                z = -nu * lnx; // z =: u in (9.1) of D.&M.(1992)

        if (b * z == 0.) { // should not happen, but does, e.g.,
            // for  pbeta(1e-320, 1e-5, 0.5)  i.e., _subnormal_ x,
            // Warning ... bgrat(a=20.5, b=1e-05, x=1, y=9.99989e-321): ..
            throw new IllegalArgumentException();
            //           MATHLIB_WARNING5(
            //                   "bgrat(a=%g, b=%g, x=%g, y=%g): z=%g, b*z == 0 underflow, hence inaccurate pbeta()",
            //                   a,b,x,y, z);
	/* L_Error:    THE EXPANSION CANNOT BE COMPUTED */
            //           *ierr = 1; return;
        }

/*                 COMPUTATION OF THE EXPANSION */
        double
	/* r1 = b * (gam1(b) + 1.) * exp(b * log(z)),// = b/gamma(b+1) z^b = z^b / gamma(b)
	 * set r := exp(-z) * z^b / gamma(b) ;
	 *          gam1(b) = 1/gamma(b+1) - 1 , b in [-1/2, 3/2] */
                // exp(a*lnx) underflows for large (a * lnx); e.g. large a ==> using log_r := log(r):
                // r = r1 * exp(a * lnx) * exp(bm1 * 0.5 * lnx);
                // log(r)=log(b) + log1p(gam1(b)) + b * log(z) + (a * lnx) + (bm1 * 0.5 * lnx),
                log_r = log(b) + log1p(gam1(b)) + b * log(z) + nu * lnx,
                // FIXME work with  log_u = log(u)  also when log_p=FALSE  (??)
                // u is 'factored out' from the expansion {and multiplied back, at the end}:
                log_u = log_r - (algdiv(b, a) + b * log(nu)),// algdiv(b,a) = log(gamma(a)/gamma(a+b))
	/* u = (log_p) ? log_r - u : exp(log_r-u); // =: M  in (9.2) of {reference above} */
	/* u = algdiv(b, a) + b * log(nu);// algdiv(b,a) = log(gamma(a)/gamma(a+b)) */
                // u = (log_p) ? log_u : exp(log_u); // =: M  in (9.2) of {reference above}
                u = exp(log_u);

        if (log_u == Double.NEGATIVE_INFINITY) {
            throw new IllegalArgumentException();
            //           R_ifDEBUG_printf(" bgrat(*): underflow log_u = -Inf  = log_r -u', log_r = %g ",
//                    log_r);
//	/* L_Error:    THE EXPANSION CANNOT BE COMPUTED */ *ierr = 2; return;
        }

        boolean u_0 = (u == 0.); // underflow --> do work with log(u) == log_u !
        double l = // := *w/u .. but with care: such that it also works when u underflows to 0:
                log_w
                        ? ((w == Double.NEGATIVE_INFINITY) ? 0. : exp(w - log_u))
                        : ((w == 0.) ? 0. : exp(log(w) - log_u));

//        R_ifDEBUG_printf(" bgrat(a=%g, b=%g, x=%g, *)\n -> u=%g, l='w/u'=%g, ",
//                a,b,x, u, l);
        double
                q_r = grat_r(b, z, log_r, eps), // = q/r of former grat1(b,z, r, &p, &q)
                v = 0.25 / (nu * nu),
                t2 = lnx * 0.25 * lnx,
                j = q_r,
                sum = j,
                t = 1., cn = 1., n2 = 0.;
        for (int n = 1; n <= n_terms_bgrat; ++n) {
            double bp2n = b + n2;
            j = (bp2n * (bp2n + 1.) * j + (z + bp2n + 1.) * t) * v;
            n2 += 2.;
            t *= t2;
            cn /= n2 * (n2 + 1.);
            int nm1 = n - 1;
            c[nm1] = cn;
            double s = 0.;
            if (n > 1) {
                double coef = b - n;
                for (int i = 1; i <= nm1; ++i) {
                    s += coef * c[i - 1] * d[nm1 - i];
                    coef += b;
                }
            }
            d[nm1] = bm1 * cn + s / n;
            double dj = d[nm1] * j;
            sum += dj;
            if (sum <= 0.) {
                throw new IllegalArgumentException();
//                R_ifDEBUG_printf(" bgrat(*): sum_n(..) <= 0; should not happen (n=%d)\n", n);
//	    /* L_Error:    THE EXPANSION CANNOT BE COMPUTED */ *ierr = 3; return;
            }
            if (abs(dj) <= eps * (sum + l)) {
//               *ierr = 0;
                break;
            } else if (n == n_terms_bgrat) { // never? ; please notify R-core if seen:
                throw new IllegalArgumentException();
//                *ierr = 4;
//                MATHLIB_WARNING5(
//                        "bgrat(a=%g, b=%g, x=%g) *no* convergence: NOTIFY R-core!\n dj=%g, rel.err=%g\n",
//                        a,b,x, dj, fabs(dj) /(sum + l));
            }
        } // for(n .. n_terms..)

/*                    ADD THE RESULTS TO W */

        if (log_w) // *w is in log space already:
            w = logsumexp(w, log_u + log(sum));
        else
            w += (u_0 ? exp(log_u + log(sum)) : u * sum);
        return w;
    }

    private static double logsumexp(final double a, final double b) {
        if (a == Double.NEGATIVE_INFINITY) {
            return b;
        } else if (a == Double.NEGATIVE_INFINITY) {
            return a;
        } else {
            final double max = Math.max(a, b);
            return Math.log(Math.exp(a - max) + Math.exp(b - max)) + max;
        }
    }

    static double apser(double a, double b, double x, double eps) {
/* -----------------------------------------------------------------------
 *     apser() yields the incomplete beta ratio  I_{1-x}(b,a)  for
 *     a <= min(eps,eps*b), b*x <= 1, and x <= 0.5,  i.e., a is very small.
 *     Use only if above inequalities are satisfied.
 * ----------------------------------------------------------------------- */

        final double g = .577215664901533;

        double tol, c, j, s, t, aj;
        double bx = b * x;

        t = x - bx;
        if (b * eps <= 0.02)
            c = log(x) + psi(b) + g + t;
        else // b > 2e13 : psi(b) ~= log(b)
            c = log(bx) + g + t;

        tol = eps * 5. * abs(c);
        j = 1.;
        s = 0.;
        do {
            j += 1.;
            t *= x - bx / j;
            aj = t / j;
            s += aj;
        } while (abs(aj) > tol);

        return -a * (c + s);
    }

    public static double logBeta(final double a, final double b) {
        if(Double.isNaN(a) || Double.isNaN(b))
            return Double.NaN;

        double p = a, q = a;
        if(b < p) p = b;/* := min(a,b) */
        if(b > q) q = b;/* := max(a,b) */

    /* both arguments must be >= 0 */
        if (p < 0)
            return Double.NaN;
        else if (p == 0) {
            return Double.POSITIVE_INFINITY;
        }
        else if (!Double.isFinite(q)) { /* q == +Inf */
            return Double.NEGATIVE_INFINITY;
        }

        final double corr;
        if (p >= 10) {
	/* p and q are big. */
            corr = Gamma.logGammaCorrection(p) + Gamma.logGammaCorrection(q) - Gamma.logGammaCorrection(p + q);
            return log(q) * -0.5 + MathConstants.LN_SQRT_2PI + corr
                    + (p - 0.5) * log(p / (p + q)) + q * log1p(-p / (p + q));
        }
        else if (q >= 10) {
	/* p is small, but q is big. */
            corr = Gamma.logGammaCorrection(q) - Gamma.logGammaCorrection(p + q);
            return Gamma.logGamma(p) + corr + p - p * log(p + q)
                    + (q - 0.5) * log1p(-p / (p + q));
        }
        else {
	/* p and q are small: p <= q < 10. */
	/* R change for very small args */
            if (p < 1e-306) return Gamma.logGamma(p) + (Gamma.logGamma(q) - Gamma.logGamma(p+q));
            else return log(Gamma.gamma(p) * (Gamma.gamma(q) / Gamma.gamma(p + q)));
        }
    }

    private static double brcomp(double a, double b, double x, double y, boolean log_p) {
/* -----------------------------------------------------------------------
 *		 Evaluation of x^a * y^b / ExponentialUtils(a,b)
 * ----------------------------------------------------------------------- */

        final double const__ = .398942280401433; /* == 1/sqrt(2*pi); */
    /* R has  M_1_SQRT_2PI , and M_LN_SQRT_2PI = ln(sqrt(2*pi)) = 0.918938.. */
        int i, n;
        double c, e, u, v, z, a0, b0, apb;

        if (x == 0. || y == 0.) {
            return log_p ? Double.NEGATIVE_INFINITY : 0;
        }
        a0 = min(a, b);
        if (a0 < 8.) {
            double lnx, lny;
            if (x <= .375) {
                lnx = log(x);
                lny = alnrel(-x);
            } else {
                if (y > .375) {
                    lnx = log(x);
                    lny = log(y);
                } else {
                    lnx = alnrel(-y);
                    lny = log(y);
                }
            }

            z = a * lnx + b * lny;
            if (a0 >= 1.) {
                z -= betaln(a, b);
                return log_p ? (z) : exp(z);
            }

/* ----------------------------------------------------------------------- */
/*		PROCEDURE FOR a < 1 OR b < 1 */
/* ----------------------------------------------------------------------- */

            b0 = max(a, b);
            if (b0 >= 8.) { /* L80: */
                u = gamln1(a0) + algdiv(a0, b0);

                return (log_p ? log(a0) + (z - u) : a0 * exp(z - u));
            }
	/* else : */

            if (b0 <= 1.) { /*		algorithm for max(a,b) = b0 <= 1 */

                double e_z = log_p ? (z) : exp(z);

                if (!log_p && e_z == 0.) /* exp() underflow */
                    return 0.;

                apb = a + b;
                if (apb > 1.) {
                    u = a + b - 1.;
                    z = (gam1(u) + 1.) / apb;
                } else {
                    z = gam1(apb) + 1.;
                }

                c = (gam1(a) + 1.) * (gam1(b) + 1.) / z;
	    /* FIXME? log(a0*c)= log(a0)+ log(c) and that is improvable */
                return (log_p
                        ? e_z + log(a0 * c) - log1p(a0 / b0)
                        : e_z * (a0 * c) / (a0 / b0 + 1.));
            }

	/* else : 		  ALGORITHM FOR 1 < b0 < 8 */

            u = gamln1(a0);
            n = (int) (b0 - 1.);
            if (n >= 1) {
                c = 1.;
                for (i = 1; i <= n; ++i) {
                    b0 += -1.;
                    c *= b0 / (a0 + b0);
                }
                u = log(c) + u;
            }
            z -= u;
            b0 += -1.;
            apb = a0 + b0;
            double t;
            if (apb > 1.) {
                u = a0 + b0 - 1.;
                t = (gam1(u) + 1.) / apb;
            } else {
                t = gam1(apb) + 1.;
            }

            return (log_p
                    ? log(a0) + z + log1p(gam1(b0)) - log(t)
                    : a0 * exp(z) * (gam1(b0) + 1.) / t);

        } else {
/* ----------------------------------------------------------------------- */
/*		PROCEDURE FOR A >= 8 AND B >= 8 */
/* ----------------------------------------------------------------------- */
            double h, x0, y0, lambda;
            if (a <= b) {
                h = a / b;
                x0 = h / (h + 1.);
                y0 = 1. / (h + 1.);
                lambda = a - (a + b) * x;
            } else {
                h = b / a;
                x0 = 1. / (h + 1.);
                y0 = h / (h + 1.);
                lambda = (a + b) * y - b;
            }

            e = -lambda / a;
            if (abs(e) > .6)
                u = e - log(x / x0);
            else
                u = rlog1(e);

            e = lambda / b;
            if (abs(e) <= .6)
                v = rlog1(e);
            else
                v = e - log(y / y0);

            z = log_p ? -(a * u + b * v) : exp(-(a * u + b * v));

            return (log_p
                    ? -LN_SQRT_2PI + .5 * log(b * x0) + z - bcorr(a, b)
                    : const__ * sqrt(b * x0) * z * exp(-bcorr(a, b)));
        }
    }

    static double psi(double x) {
/* ---------------------------------------------------------------------
 *                 Evaluation of the Digamma function psi(x)
 *                           -----------
 *     Psi(xx) is assigned the value 0 when the digamma function cannot
 *     be computed.
 *     The main computation involves evaluation of rational Chebyshev
 *     approximations published in Math. Comp. 27, 123-127(1973) by
 *     Cody, Strecok and Thacher. */

/* --------------------------------------------------------------------- */
/*     Psi was written at Argonne National Laboratory for the FUNPACK */
/*     package of special function subroutines. Psi was modified by */
/*     A.H. Morris (NSWC). */
/* --------------------------------------------------------------------- */

        final double piov4 = .785398163397448; /* == pi / 4 */
/*     dx0 = zero of psi() to extended precision : */
        final double dx0 = 1.461632144968362341262659542325721325;

/* --------------------------------------------------------------------- */
/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF */
/*     PSI(X) / (X - X0),  0.5 <= X <= 3. */
        final double p1[] = {.0089538502298197, 4.77762828042627,
                142.441585084029, 1186.45200713425, 3633.51846806499,
                4138.10161269013, 1305.60269827897};
        final double q1[] = {44.8452573429826, 520.752771467162,
                2210.0079924783, 3641.27349079381, 1908.310765963,
                6.91091682714533e-6};
/* --------------------------------------------------------------------- */


/* --------------------------------------------------------------------- */
/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF */
/*     PSI(X) - LN(X) + 1 / (2*X),  X > 3. */

        final double p2[] = {-2.12940445131011, -7.01677227766759,
                -4.48616543918019, -.648157123766197};
        final double q2[] = {32.2703493791143, 89.2920700481861,
                54.6117738103215, 7.77788548522962};
/* --------------------------------------------------------------------- */

        int i, m, n, nq;
        double d2;
        double w, z;
        double den, aug, sgn, xmx0, xmax1, upper, xsmall;

/* --------------------------------------------------------------------- */


/*     MACHINE DEPENDENT CONSTANTS ... */

/* --------------------------------------------------------------------- */
/*	  XMAX1	 = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
		   WITH ENTIRELY INT REPRESENTATION.  ALSO USED
		   AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
		   ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
		   PSI MAY BE REPRESENTED AS LOG(X).
 * Originally:  xmax1 = amin1(ipmpar(3), 1./spmpar(1))  */
        xmax1 = (double) Integer.MAX_VALUE;
        d2 = 1 / DBL_EPSILON; /*= 0.5 / (0.5 * DBL_EPS) = 1/DBL_EPSILON = 2^52 */
        if (xmax1 > d2) xmax1 = d2;

/* --------------------------------------------------------------------- */
/*        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X) */
/*                 MAY BE REPRESENTED BY 1/X. */
        xsmall = 1e-9;
/* --------------------------------------------------------------------- */
        aug = 0.;
        if (x < 0.5) {
/* --------------------------------------------------------------------- */
/*     X < 0.5,  USE REFLECTION FORMULA */
/*     PSI(1-X) = PSI(X) + PI * COTAN(PI*X) */
/* --------------------------------------------------------------------- */
            if (abs(x) <= xsmall) {

                if (x == 0.) {
                    return 0.;
                }
/* --------------------------------------------------------------------- */
/*     0 < |X| <= XSMALL.  USE 1/X AS A SUBSTITUTE */
/*     FOR  PI*COTAN(PI*X) */
/* --------------------------------------------------------------------- */
                aug = -1. / x;
            } else { /* |x| > xsmall */
/* --------------------------------------------------------------------- */
/*     REDUCTION OF ARGUMENT FOR COTAN */
/* --------------------------------------------------------------------- */
	    /* L100: */
                w = -x;
                sgn = piov4;
                if (w <= 0.) {
                    w = -w;
                    sgn = -sgn;
                }
/* --------------------------------------------------------------------- */
/*     MAKE AN ERROR EXIT IF |X| >= XMAX1 */
/* --------------------------------------------------------------------- */
                if (w >= xmax1) {
                    return 0.;
                }
                nq = (int) w;
                w -= (double) nq;
                nq = (int) (w * 4.);
                w = (w - (double) nq * 0.25) * 4.;
/* --------------------------------------------------------------------- */
/*     W IS NOW RELATED TO THE FRACTIONAL PART OF  4. * X. */
/*     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST */
/*     QUADRANT AND DETERMINE SIGN */
/* --------------------------------------------------------------------- */
                n = nq / 2;
                if (n + n != nq) {
                    w = 1. - w;
                }
                z = piov4 * w;
                m = n / 2;
                if (m + m != n) {
                    sgn = -sgn;
                }
/* --------------------------------------------------------------------- */
/*     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X) */
/* --------------------------------------------------------------------- */
                n = (nq + 1) / 2;
                m = n / 2;
                m += m;
                if (m == n) {
/* --------------------------------------------------------------------- */
/*     CHECK FOR SINGULARITY */
/* --------------------------------------------------------------------- */
                    if (z == 0.) {
                        return 0.;
                    }
/* --------------------------------------------------------------------- */
/*     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND */
/*     SIN/COS AS A SUBSTITUTE FOR TAN */
/* --------------------------------------------------------------------- */
                    aug = sgn * (cos(z) / sin(z) * 4.);

                } else { /* L140: */
                    aug = sgn * (sin(z) / cos(z) * 4.);
                }
            }

            x = 1. - x;

        }
    /* L200: */
        if (x <= 3.) {
/* --------------------------------------------------------------------- */
/*     0.5 <= X <= 3. */
/* --------------------------------------------------------------------- */
            den = x;
            upper = p1[0] * x;

            for (i = 1; i <= 5; ++i) {
                den = (den + q1[i - 1]) * x;
                upper = (upper + p1[i]) * x;
            }

            den = (upper + p1[6]) / (den + q1[5]);
            xmx0 = x - dx0;
            return den * xmx0 + aug;
        }

/* --------------------------------------------------------------------- */
/*     IF X >= XMAX1, PSI = LN(X) */
/* --------------------------------------------------------------------- */
        if (x < xmax1) {
/* --------------------------------------------------------------------- */
/*     3. < X < XMAX1 */
/* --------------------------------------------------------------------- */
            w = 1. / (x * x);
            den = w;
            upper = p2[0] * w;

            for (i = 1; i <= 3; ++i) {
                den = (den + q2[i - 1]) * w;
                upper = (upper + p2[i]) * w;
            }

            aug = upper / (den + q2[3]) - 0.5 / x + aug;
        }
        return aug + log(x);

/* --------------------------------------------------------------------- */
/*     ERROR RETURN */
/* --------------------------------------------------------------------- */

    }

    static double grat_r(double a, double x, double log_r, double eps) {
/* -----------------------------------------------------------------------
 *        Scaled complement of incomplete gamma ratio function
 *                   grat_r(a,x,r) :=  Q(a,x) / r
 * where
 *               Q(a,x) = pgamma(x,a, lower.tail=FALSE)
 *     and            r = e^(-x)* x^a / Gamma(a) ==  exp(log_r)
 *
 *     It is assumed that a <= 1.  eps is the tolerance to be used.
 * ----------------------------------------------------------------------- */

        if (a * x == 0.) { /* L130: */
            if (x <= a) {
	    /* L100: */
                return exp(-log_r);
            } else {
	    /* L110:*/
                return 0.;
            }
        } else if (a == 0.5) { // e.g. when called from pt()
	/* L120: */
            if (x < 0.25) {
                double p = erf__(sqrt(x));
                return (0.5 - p + 0.5) * exp(-log_r);

            } else { // 2013-02-27: improvement for "large" x: direct computation of q/r:
                double sx = sqrt(x),
                        q_r = erfc1(1, sx) / sx * SQRT_PI;
                return q_r;
            }

        } else if (x < 1.1) { /* L10:  Taylor series for  P(a,x)/x^a */

            double an = 3.,
                    c = x,
                    sum = x / (a + 3.),
                    tol = eps * 0.1 / (a + 1.), t;
            do {
                an += 1.;
                c *= -(x / an);
                t = c / (a + an);
                sum += t;
            } while (abs(t) > tol);

            double j = a * x * ((sum / 6. - 0.5 / (a + 2.)) * x + 1. / (a + 1.)),
                    z = a * log(x),
                    h = gam1(a),
                    g = h + 1.;

            if ((x >= 0.25 && (a < x / 2.59)) || (z > -0.13394)) {
                // L40:
                double l = rexpm1(z),
                        q = ((l + 0.5 + 0.5) * j - l) * g - h;
                if (q <= 0.) {
		/* L110:*/
                    return 0.;
                } else {
                    return q * exp(-log_r);
                }

            } else {
                double p = exp(z) * g * (0.5 - j + 0.5);
                return /* q/r = */ (0.5 - p + 0.5) * exp(-log_r);
            }

        } else {
	/* L50: ----  (x >= 1.1)  ---- Continued Fraction Expansion */

            double a2n_1 = 1.,
                    a2n = 1.,
                    b2n_1 = x,
                    b2n = x + (1. - a),
                    c = 1., am0, an0;

            do {
                a2n_1 = x * a2n + c * a2n_1;
                b2n_1 = x * b2n + c * b2n_1;
                am0 = a2n_1 / b2n_1;
                c += 1.;
                double c_a = c - a;
                a2n = a2n_1 + c_a * a2n;
                b2n = b2n_1 + c_a * b2n;
                an0 = a2n / b2n;
            } while (abs(an0 - am0) >= eps * an0);

            return /* q/r = (r * an0)/r = */ an0;
        }
    } /* grat_r */

    static double rexpm1(double x) {
/* ----------------------------------------------------------------------- */
/*            EVALUATION OF THE FUNCTION EXP(X) - 1 */
/* ----------------------------------------------------------------------- */

        final double p1 = 9.14041914819518e-10;
        final double p2 = .0238082361044469;
        final double q1 = -.499999999085958;
        final double q2 = .107141568980644;
        final double q3 = -.0119041179760821;
        final double q4 = 5.95130811860248e-4;

        if (abs(x) <= 0.15) {
            return x * (((p2 * x + p1) * x + 1.) /
                    ((((q4 * x + q3) * x + q2) * x + q1) * x + 1.));
        } else { /* |x| > 0.15 : */
            double w = exp(x);
            if (x > 0.)
                return w * (0.5 - 1. / w + 0.5);
            else
                return w - 0.5 - 0.5;
        }

    }

    static double erf__(double x) {
/* -----------------------------------------------------------------------
 *             EVALUATION OF THE REAL ERROR FUNCTION
 * ----------------------------------------------------------------------- */

    /* Initialized data */

        final double c = .564189583547756;
        final double a[] = {7.7105849500132e-5, -.00133733772997339,
                .0323076579225834, .0479137145607681, .128379167095513};
        final double b[] = {.00301048631703895, .0538971687740286,
                .375795757275549};
        final double p[] = {-1.36864857382717e-7, .564195517478974,
                7.21175825088309, 43.1622272220567, 152.98928504694,
                339.320816734344, 451.918953711873, 300.459261020162};
        final double q[] = {1., 12.7827273196294, 77.0001529352295,
                277.585444743988, 638.980264465631, 931.35409485061,
                790.950925327898, 300.459260956983};
        final double r[] = {2.10144126479064, 26.2370141675169,
                21.3688200555087, 4.6580782871847, .282094791773523};
        final double s[] = {94.153775055546, 187.11481179959,
                99.0191814623914, 18.0124575948747};

    /* Local variables */
        double t, x2, ax, bot, top;

        ax = abs(x);
        if (ax <= 0.5) {
            t = x * x;
            top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.;
            bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.;

            return x * (top / bot);
        }

        // else:  |x| > 0.5

        if (ax <= 4.) { //  |x| in (0.5, 4]
            top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
                    + p[5]) * ax + p[6]) * ax + p[7];
            bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
                    + q[5]) * ax + q[6]) * ax + q[7];
            double R = 0.5 - exp(-x * x) * top / bot + 0.5;
            return (x < 0) ? -R : R;
        }

        // else:  |x| > 4

        if (ax >= 5.8) {
            return x > 0 ? 1 : -1;
        }

        // else:  4 < |x| < 5.8
        x2 = x * x;
        t = 1. / x2;
        top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
        bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.;
        t = (c - top / (x2 * bot)) / ax;
        double R = 0.5 - exp(-x2) * t + 0.5;
        return (x < 0) ? -R : R;
    } /* erf */

    static double erfc1(int ind, double x) {
/* ----------------------------------------------------------------------- */
/*         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION */

/*          ERFC1(IND,X) = ERFC(X)            IF IND = 0 */
/*          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE */
/* ----------------------------------------------------------------------- */

    /* Initialized data */

        final double c = .564189583547756;
        final double a[] = {7.7105849500132e-5, -.00133733772997339,
                .0323076579225834, .0479137145607681, .128379167095513};
        final double b[] = {.00301048631703895, .0538971687740286,
                .375795757275549};
        final double p[] = {-1.36864857382717e-7, .564195517478974,
                7.21175825088309, 43.1622272220567, 152.98928504694,
                339.320816734344, 451.918953711873, 300.459261020162};
        final double q[] = {1., 12.7827273196294, 77.0001529352295,
                277.585444743988, 638.980264465631, 931.35409485061,
                790.950925327898, 300.459260956983};
        final double r[] = {2.10144126479064, 26.2370141675169,
                21.3688200555087, 4.6580782871847, .282094791773523};
        final double s[] = {94.153775055546, 187.11481179959,
                99.0191814623914, 18.0124575948747};

        double ret_val;
        double e, t, w, bot, top;

        double ax = abs(x);
        //				|X| <= 0.5 */
        if (ax <= 0.5) {
            t = x * x;
            top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.;
            bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.;
            ret_val = 0.5 - x * (top / bot) + 0.5;
            if (ind != 0) {
                ret_val = exp(t) * ret_val;
            }
            return ret_val;
        }
        // else (L10:):		0.5 < |X| <= 4
        if (ax <= 4.) {
            top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
                    + p[5]) * ax + p[6]) * ax + p[7];
            bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
                    + q[5]) * ax + q[6]) * ax + q[7];
            ret_val = top / bot;

        } else { //			|X| > 4
            // L20:
            if (x <= -5.6) {
                // L50:            	LIMIT VALUE FOR "LARGE" NEGATIVE X
                ret_val = 2.;
                if (ind != 0) {
                    ret_val = exp(x * x) * 2.;
                }
                return ret_val;
            }
            if (ind == 0 && (x > 100. || x * x > -exparg(1))) {
                // LIMIT VALUE FOR LARGE POSITIVE X   WHEN IND = 0
                // L60:
                return 0.;
            }

            // L30:
            t = 1. / (x * x);
            top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
            bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.;
            ret_val = (c - t * top / bot) / ax;
        }

        // L40:                 FINAL ASSEMBLY
        if (ind != 0) {
            if (x < 0.)
                ret_val = exp(x * x) * 2. - ret_val;
        } else {
            // L41:  ind == 0 :
            w = x * x;
            t = w;
            e = w - t;
            ret_val = (0.5 - e + 0.5) * exp(-t) * ret_val;
            if (x < 0.)
                ret_val = 2. - ret_val;
        }
        return ret_val;

    }

    static double basym(double a, double b, double lambda, double eps, boolean log_p) {
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR I_x(A,B) FOR LARGE A AND B. */
/*     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED. */
/*     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT */
/*     A AND B ARE GREATER THAN OR EQUAL TO 15. */
/* ----------------------------------------------------------------------- */


/* ------------------------ */
/*     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP */
/*            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN. */
        final int num_IT = 20;
/*            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1. */

        final double e0 = 1.12837916709551;/* e0 == 2/sqrt(pi) */
        final double e1 = .353553390593274;/* e1 == 2^(-3/2)   */
        final double ln_e0 = 0.120782237635245; /* == ln(e0) */

        final double[] a0 = new double[num_IT + 1], b0 = new double[num_IT + 1],
                c = new double[num_IT + 1], d = new double[num_IT + 1];

        double f = a * rlog1(-lambda / a) + b * rlog1(lambda / b), t;
        if (log_p)
            t = -f;
        else {
            t = exp(-f);
            if (t == 0.) {
                return 0; /* once underflow, always underflow .. */
            }
        }
        double z0 = sqrt(f),
                z = z0 / e1 * 0.5,
                z2 = f + f,
                h, r0, r1, w0;

        if (a < b) {
            h = a / b;
            r0 = 1. / (h + 1.);
            r1 = (b - a) / b;
            w0 = 1. / sqrt(a * (h + 1.));
        } else {
            h = b / a;
            r0 = 1. / (h + 1.);
            r1 = (b - a) / a;
            w0 = 1. / sqrt(b * (h + 1.));
        }

        a0[0] = r1 * .66666666666666663;
        c[0] = a0[0] * -0.5;
        d[0] = -c[0];
        double j0 = 0.5 / e0 * erfc1(1, z0),
                j1 = e1,
                sum = j0 + d[0] * w0 * j1;

        double s = 1.,
                h2 = h * h,
                hn = 1.,
                w = w0,
                znm1 = z,
                zn = z2;
        for (int n = 2; n <= num_IT; n += 2) {
            hn *= h2;
            a0[n - 1] = r0 * 2. * (h * hn + 1.) / (n + 2.);
            int np1 = n + 1;
            s += hn;
            a0[np1 - 1] = r1 * 2. * s / (n + 3.);

            for (int i = n; i <= np1; ++i) {
                double r = (i + 1.) * -0.5;
                b0[0] = r * a0[0];
                for (int m = 2; m <= i; ++m) {
                    double bsum = 0.;
                    for (int j = 1; j <= m - 1; ++j) {
                        int mmj = m - j;
                        bsum += (j * r - mmj) * a0[j - 1] * b0[mmj - 1];
                    }
                    b0[m - 1] = r * a0[m - 1] + bsum / m;
                }
                c[i - 1] = b0[i - 1] / (i + 1.);

                double dsum = 0.;
                for (int j = 1; j <= i - 1; ++j) {
                    dsum += d[i - j - 1] * c[j - 1];
                }
                d[i - 1] = -(dsum + c[i - 1]);
            }

            j0 = e1 * znm1 + (n - 1.) * j0;
            j1 = e1 * zn + n * j1;
            znm1 = z2 * znm1;
            zn = z2 * zn;
            w *= w0;
            double t0 = d[n - 1] * w * j0;
            w *= w0;
            double t1 = d[np1 - 1] * w * j1;
            sum += t0 + t1;
            if (abs(t0) + abs(t1) <= eps * sum) {
                break;
            }
        }

        if (log_p)
            return ln_e0 + t - bcorr(a, b) + log(sum);
        else {
            double u = exp(-bcorr(a, b));
            return e0 * t * u * sum;
        }

    }
}
