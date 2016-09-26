package vrr.jargon.exponential;

import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import org.testng.Assert;
import org.testng.annotations.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by valentin on 9/25/16.
 */
public class BinomialTest extends RTest {

    @BeforeTest
    public void setUp () throws REngineException {
        super.setUp();
    }

    @AfterTest
    public void tearDown() {
        super.tearDown();
    }

    @Test(dataProvider="densityData")
    public void testDensity(final double p, final int n, final int x) throws REXPMismatchException, REngineException {
        final double actual = Binomial.density(x, n, p);
        final double expected = rCall("dbinom")
                .with(x, n, p)
                .evalAsDouble();
        assertEquals(actual, expected);
    }

    @Test(dataProvider="densityData")
    public void testLogDensity(final double p, final int n, final int x) throws REXPMismatchException, REngineException {
        final double actual = Binomial.logDensity(x, n, p);
        final double expected = rCall("dbinom")
                .with(x, n, p)
                .with("log", true)
                .evalAsDouble();
        assertEquals(actual, expected);
    }

    @Test(dataProvider = "densityData")
    public void testCDF(final double p, final int n, final int x) throws REXPMismatchException, REngineException {
        final double actual = Binomial.CDF(x, n, p);
        final double expected = rCall("pbinom")
                .with(x, n, p)
                .evalAsDouble();
        assertEquals(actual, expected);
    }

    @Test(dataProvider = "inverseCdfData")
    public void testInverseCDF(final double p, final int n, final double q) throws REXPMismatchException, REngineException {
        final int actual = Binomial.inverseCDF(q, n, p);
        final double expected = rCall("qbinom")
                .with(q, n, p)
                .evalAsDouble();
        if (Double.isNaN(expected))
            Assert.assertEquals(actual, -1);
        else
            Assert.assertEquals(actual, (int) Math.round(expected));
    }

    @Test(dataProvider = "inverseCdfData")
    public void testLogInverseCDF(final double p, final int n, final double q) throws REXPMismatchException, REngineException {
        final int actual = Binomial.logInverseCDF(Math.log(q), n, p);
        final double expected = rCall("qbinom")
                .with(Math.log(q), n, p)
                .with("log.p", true)
                .evalAsDouble();
        if (Double.isNaN(expected))
            Assert.assertEquals(actual, -1);
        else
            Assert.assertEquals(actual, (int) Math.round(expected));
    }

    @Test(dataProvider = "logInverseCdfExtraData")
    public void testLogInverseCDFExtra(final double p, final int n, final double q) throws REXPMismatchException, REngineException {
        final int actual = Binomial.logInverseCDF(q, n, p);
        final double expected = rCall("qbinom")
                .with(q, n, p)
                .with("log.p", true)
                .evalAsDouble();
        if (Double.isNaN(expected))
            Assert.assertEquals(actual, -1);
        else
            Assert.assertEquals(actual, (int) Math.round(expected));
    }

    @Test(dataProvider = "densityData")
    public void testLogCDF(final double p, final int n, final int x) throws REXPMismatchException, REngineException {
        final double actual = Binomial.logCDF(x, n, p);
        final double expected = rCall("pbinom")
                .with(x, n, p)
                .with("log", true)
                .evalAsDouble();
        assertEquals(actual, expected);
    }

    @DataProvider(name="densityData")
    public Object[][] densityData() {
        final List<Object[]> result = new ArrayList<>();
        final int[] ns = new int[]{0, 1, 2, 3, 10, 100};
        final double[] probs = new double[]{Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0, 1e-100, 1e-20, 1e-10, 1e-5, 0.1, .99, 0.99999999999};
        for (final double p : probs)
            for (final int n : ns)
                for (int x = 0; x < n; x++) {
                    result.add(new Object[]{p, n, x});
                }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="inverseCdfData")
    public Object[][] inverseCdfData() {
        final List<Object[]> result = new ArrayList<>();
        final int[] ns = new int[]{0, 1, 2, 3, 10, 100};
        final double[] probs = new double[]{Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0, 1e-100, 1e-20, 1e-10, 1e-5, 0.1, .99, 0.99999999999, 1.0};
        for (final double p : probs)
            for (final int n : ns)
                for (int x = 0; x < probs.length; x++) {
                    result.add(new Object[]{p, n, probs[x]});
                }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="logInverseCdfExtraData")
    public Object[][] logInverseCdfExtraData() {
        final List<Object[]> result = new ArrayList<>();
        final int[] ns = new int[]{0, 1, 2, 3, 10, 100};
        final double[] probs = new double[]{Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0, 1e-100, 1e-20, 1e-10, 1e-5, 0.1, .99, 0.99999999999};
        final double[] qs = new double[]{Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, -100, -10, -0.1, -0.0000000000001, 0.0};
        for (final double p : probs)
            for (final int n : ns)
                for (final double q : qs)
                    result.add(new Object[]{p, n, q});

        return result.toArray(new Object[result.size()][]);
    }

}
