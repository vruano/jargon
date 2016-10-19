package vrr.jargon.exponential;

import com.google.caliper.Benchmark;
import com.google.caliper.Param;
import com.google.caliper.api.Macrobenchmark;
import com.google.caliper.api.VmOptions;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineEvalException;
import org.rosuda.REngine.REngineException;

import java.lang.annotation.Annotation;
import java.util.Random;

/**
 * Created by valentin on 10/12/16.
 */
@VmOptions("-XX:-TieredCompilation")
public class GammaInverseCDFBenchmark {

    @Param({"0.01", "0.1"})
    double shape;

    @Param({"0.01", "0.1"})
    double scale;

    private static final double[] randomProbs = new double[1000];

    static {
        final Random rdn = new Random(131313);
        for (int i = 0; i < randomProbs.length; i++)
            randomProbs[i] = rdn.nextDouble() * 0.999 + 0.0005;
    }

    @Macrobenchmark
    public double jargon(final int reps) {
        double sum = 0;
        for (int r = 0; r < reps; r++) {
            for (int i = 0; i < randomProbs.length; i++)
                sum += Gamma.inverseCDF(randomProbs[i], shape, scale);
        }
        return sum;
    }

    @Macrobenchmark
    public double apache(final int reps) {
        final GammaDistribution dist = new GammaDistributionImpl(shape, scale);
        double sum = 0;
        for (int r = 0; r < reps; r++) {
            for (int i = 0; i < randomProbs.length; i++) {
                try {
                    sum += dist.inverseCumulativeProbability(randomProbs[i]);
                } catch (final MathException ex) {
                    System.err.println("ex : " + randomProbs[i] + " " + shape + " " + scale);
                }
            }
        }
        return sum;
    }

    @Macrobenchmark
    public double rCall(final int reps) throws REXPMismatchException, REngineException {
        double sum = 0;
        new RTest().setUp();
        for (int r = 0; r < reps; r++) {
            final double[] values = new RTest.RCall("qgamma").withArray(randomProbs).with(shape).with("scale", scale).evalAsDoubleArray();
            for (final double v : values)
                sum += v;
        }
        return sum;
    }
}
