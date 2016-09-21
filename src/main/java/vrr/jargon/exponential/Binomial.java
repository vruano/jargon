package vrr.jargon.exponential;

import vrr.jargon.utils.DoubleConstants;

/**
 * Created by valentin on 9/20/16.
 */
public class Binomial {

    public static double deviance(final double x, final double np) {
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
}
