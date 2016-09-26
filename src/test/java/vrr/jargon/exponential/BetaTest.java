package vrr.jargon.exponential;

import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import org.testng.Assert;
import org.testng.annotations.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by valentin on 9/22/16.
 */
public class BetaTest extends RTest {

    @BeforeClass
    public void setUp() throws REngineException {
        super.setUp();
    }

    @AfterClass
    public void tearDown() {
        super.tearDown();
    }

    @Test(dataProvider= "inverseCdfData")
    public void testInverseCDF(final double alpha, final double beta, final double p) throws REXPMismatchException, REngineException {
        final double actual = Beta.inverseCDF(p, alpha, beta);
        final double expected = rCall("qbeta")
                .with(p, alpha, beta)
                .evalAsDouble();
        assertEquals(actual, expected);
    }

    @Test(dataProvider= "inverseCdfData")
    public void testLogInverseCDF(final double alpha, final double beta, final double p) throws REXPMismatchException, REngineException {
        final double actual = Beta.logInverseCDF(Math.log(p), alpha, beta);
        final double expected = rCall("qbeta")
                .with(Math.log(p), alpha, beta)
                .with("log", true)
                .evalAsDouble();
        if (expected == 1.112536929253566E-308) { // This is an special case in where R's qbeta uses a 'exp' operation that is
            // an approximation to the actual 'exp' resulting in a different outcome for very small numbers...
            // it does not even agree with the 'exp' operation when called directly in R command-line.
            Assert.assertEquals(actual, 5.562684646268137E-309, "");
        } else
            assertEquals(actual, expected);
    }

    @Test(dataProvider= "betaData")
    public void testBeta(final double alpha, final double beta) throws REXPMismatchException, REngineException {
        final double actual = Beta.beta(alpha, beta);
        final double expected = rCall("beta")
                .with(alpha, beta)
                .evalAsDouble();
        assertEquals(actual, expected);
    }

    @Test(dataProvider= "betaData")
    public void testLogBeta(final double alpha, final double beta) throws REXPMismatchException, REngineException {
        final double actual = Beta.logBeta(alpha, beta);
        final double expected = rCall("lbeta")
                .with(alpha, beta)
                .evalAsDouble();
        assertEquals(actual, expected);
    }

    @Test(dataProvider= "cdfData")
    public void testCDF(final double alpha, final double beta, final double x) throws REXPMismatchException, REngineException {
        final double actual = Beta.CDF(x, alpha, beta);
        final double expected = rCall("pbeta")
                .with(x, alpha, beta)
                .evalAsDouble();
        assertEquals(actual, expected);
    }

    @Test(dataProvider= "cdfData")
    public void testLogCDF(final double alpha, final double beta, final double x) throws REXPMismatchException, REngineException {
        final double actual = Beta.logCDF(x, alpha, beta);
        final double expected = rCall("pbeta")
                .with(x, alpha, beta)
                .with("log.p", true)
                .evalAsDouble();
        assertEquals(actual, expected);
    }

    @DataProvider(name="inverseCdfData")
    public Object[][] inverseCdfData() throws IOException {
        final List<Object[]> result = new ArrayList<>();

        final double[] alphas = new double[] { Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0.01, 0.1, 0.5, 1, 1.5, 3, 10, 30};
        final double[] betas = new double[] { Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0.01, 0.1, 0.5, 1, 1.5, 3, 10, 30, 300};
        final double[] probs = new double[] { Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 0.9999, 0.99999999999999, 1 };
        for (double a: alphas)
            for (double b: betas)
                for (double p: probs) {
                    result.add(new Object[] { a, b, p});
                }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="cdfData")
    public Object[][] cdfData() throws IOException {
        final List<Object[]> result = new ArrayList<>();

        final double[] alphas = new double[] { Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0.01, 0.1, 0.5, 1, 1.5, 3, 10, 30};
        final double[] betas = new double[] { Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0.01, 0.1, 0.5, 1, 1.5, 3, 10, 30, 300};
        final double[] xs = new double[] { Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 0.9999, 0.99999999999999, 1, 2, 10, 20, -1, 200 };
        for (double a: alphas)
            for (double b: betas)
                for (double x: xs) {
                    result.add(new Object[] { a, b, x});
                }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="betaData")
    public Object[][] betaData() throws IOException {
        final List<Object[]> result = new ArrayList<>();

        final double[] alphas = new double[] { Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0.01, 0.1, 0.5, 1, 1.5, 3, 10, 30};
        final double[] betas = new double[] { Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0.01, 0.1, 0.5, 1, 1.5, 3, 10, 30, 300};
        for (double a: alphas)
            for (double b: betas)
                result.add(new Object[] { a, b});
        return result.toArray(new Object[result.size()][]);
    }
}
