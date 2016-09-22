package vrr.jargon.exponential;

import org.rosuda.REngine.*;
import org.rosuda.REngine.JRI.JRIEngine;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import vrr.jargon.utils.DoubleConstants;

import java.util.ArrayList;
import java.util.List;

/**
 * Unit test for {@link Gamma}.
 */
public class GammaTest extends RTest {

    @BeforeClass
    public void setUp() throws REngineException {
        super.setUp();
    }

    @AfterClass
    public void tearDown() {
        super.tearDown();
    }

    @Test(dataProvider = "cdfData")
    public void testCDF(final double shape, final double scale, final double x) throws REXPMismatchException, REngineException {

        final double expected = rCall("pgamma")
                .with(x, shape)
                .with("scale", scale)
                .evalAsDouble();

        final double actual = Gamma.CDF(x, shape, scale);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(expected * 0.000000001));
    }

    @Test(dataProvider = "cdfData")
    public void testDensity(final double shape, final double scale, final double x) throws REXPMismatchException, REngineException {

        final double expected = rCall("dgamma")
                .with(x, shape)
                .with("scale", scale)
                .evalAsDouble();

        final double actual = Gamma.density(x, shape, scale);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(expected * 0.000000001));
    }

    @Test(dataProvider = "cdfData")
    public void testLogDensity(final double shape, final double scale, final double x) throws REXPMismatchException, REngineException {

        final double expected = rCall("dgamma")
                .with(x, shape)
                .with("scale", scale)
                .with("log", true)
                .evalAsDouble();

        final double actual = Gamma.logDensity(x, shape, scale);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(expected * 0.000000001));
    }


    @Test(dataProvider = "cdfData")
    public void testLogCDF(final double shape, final double scale, final double x) throws REXPMismatchException, REngineException {

        final double expected = rCall("pgamma")
                .with(x, shape)
                .with("scale", scale)
                .with("log.p", true)
                .evalAsDouble();

        final double actual = Gamma.logCDF(x, shape, scale);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(expected * 0.000000001));
    }

    @Test(dataProvider = "inverseCdfData")
    public void testInverseCDF(final double shape, final double scale, final double p) throws REXPMismatchException, REngineException {

        final double expected = rCall("qgamma")
                .with(p, shape)
                .with("scale", scale)
                .evalAsDouble();

        final double actual = Gamma.inverseCDF(p, shape, scale);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(expected * 0.000000001));
    }

    @Test(dataProvider = "inverseCdfData")
    public void testLogInverseCDF(final double shape, final double scale, final double p) throws REXPMismatchException, REngineException {

        final double expected = rCall("qgamma")
                .with(Math.log(p), shape)
                .with("scale", scale)
                .with("log.p", true)
                .evalAsDouble();

        final double actual = Gamma.logInverseCDF(Math.log(p), shape, scale);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(expected * 0.000000001));
    }

    @Test(dataProvider = "logInverseCdfExtraData")
    public void testLogInverseCDFExtra(final double shape, final double scale, final double logP) throws REXPMismatchException, REngineException {

        final double expected = rCall("qgamma")
                .with(logP, shape)
                .with("scale", scale)
                .with("log.p", true)
                .evalAsDouble();

        final double actual = Gamma.logInverseCDF(logP, shape, scale);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(expected * 0.000000001));
    }

    @Test(dataProvider = "gammaData")
    public void testGamma(final double v) throws REXPMismatchException, REngineException {
        final double expected = rCall("gamma")
                .with(v)
                .evalAsDouble();
        final double actual = Gamma.gamma(v);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(expected * 0.000000001));
    }

    @Test(dataProvider = "gammaData")
    public void testLogGamma(final double v) throws REXPMismatchException, REngineException {
        final double expected = rCall("lgamma")
                .with(v)
                .evalAsDouble();
        final double actual = Gamma.logGamma(v);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(expected * 0.000000001));
    }

    @DataProvider(name = "cdfData")
    public Object[][] cdfData() {
        final List<Object[]> result = new ArrayList<>();
        result.add(new Object[] { 0.0, 0.0, 0.0 });
        result.add(new Object[] { Double.NaN, 0.0, 0.0 });
        result.add(new Object[] { 0.0, Double.NaN, 0.0 });
        result.add(new Object[] { 0.0, 0.0, Double.NaN });
        result.add(new Object[] { 1, 2, 0.5 });
        result.add(new Object[] { .1, .2, 0.5 });

        final double[] shapes = new double[] { 1e-100, 0.000001, 0.0001, 0.001, 0.1, 1, 1.1, 1.3, 2.4, 5.0, 12.1, 1241.0, 121241.0 };
        final double[] scales = new double[] { 1e-100, 0.000001, 0.0001, 0.001, 0.1, 1, 1.1, 1.3, 2.4, 5.0, 12.1, 12313.0, 121133411.0};
        final double[] xs = new double[] { Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0, 0.1, 10, 100, 12212.0};
        for (final double shape : shapes)
            for (final double scale : scales)
                for (final double x : xs) {
                    result.add(new Object[] { shape, scale, x});
                }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name = "inverseCdfData")
    public Object[][] inverseCdfData() {
        final List<Object[]> result = new ArrayList<>();
        result.add(new Object[] { 0.0, 0.0, 0.0 });
        result.add(new Object[] { Double.NaN, 0.0, 0.0 });
        result.add(new Object[] { 0.0, Double.NaN, 0.0 });
        result.add(new Object[] { 0.0, 0.0, Double.NaN });
        result.add(new Object[] { 1, 2, 0.25 });
        result.add(new Object[] { .1, .2, 0.75 });

        final double[] shapes = new double[] { 1e-100, 0.000001, 0.0001, 0.001, 0.1, 1, 1.1, 1.3, 2.4, 5.0, 12.1, 1241.0, 121241.0 };
        final double[] scales = new double[] { 1e-100, 0.000001, 0.0001, 0.001, 0.1, 1, 1.1, 1.3, 2.4, 5.0, 12.1, 12313.0, 121133411.0};
        final double[] probs = new double[] { -10000.0, -1.0, -1.0e-10, -1.3e-23};
        for (final double shape : shapes)
            for (final double scale : scales)
                for (final double p : probs) {
                    result.add(new Object[] { shape, scale, p});
                }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name = "logInverseCdfExtraData")
    public Object[][] logInverseCdfExtraData() {
        final List<Object[]> result = new ArrayList<>();
        result.add(new Object[] { 0.0, 0.0, Double.NEGATIVE_INFINITY });
        result.add(new Object[] { Double.NaN, 0.0, Double.NEGATIVE_INFINITY});
        result.add(new Object[] { 0.0, Double.NaN, Double.NEGATIVE_INFINITY });
        result.add(new Object[] { 1, 2, 0.0 });
        result.add(new Object[] { .1, .2, 0.0 });

        final double[] shapes = new double[] { 1e-100, 0.000001, 0.0001, 0.001, 0.1, 1, 1.1, 1.3, 2.4, 5.0, 12.1, 1241.0, 121241.0 };
        final double[] scales = new double[] { 1e-100, 0.000001, 0.0001, 0.001, 0.1, 1, 1.1, 1.3, 2.4, 5.0, 12.1, 12313.0, 121133411.0};
        final double[] probs = new double[] { 0, 1e-100, 3.2e-3, 0.1, 0.5, 0.654, 0.9, 0.9999999999};
        for (final double shape : shapes)
            for (final double scale : scales)
                for (final double p : probs) {
                    result.add(new Object[] { shape, scale, p});
                }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name = "gammaData")
    public Object[][] gammaData() {
        return new Object[][] {
                new Object[] { Double.NaN },
                new Object[] { Double.POSITIVE_INFINITY },
                new Object[] { Double.NEGATIVE_INFINITY },
                new Object[] { .00001 },
                new Object[] { 1e-100 },
                new Object[] { .1 },
                new Object[] { 3 },
                new Object[] { 20 },
                new Object[] { 2000.1 },
                new Object[] { -.1 },
                new Object[] { -10.0 },
                new Object[] { -20 },
                new Object[] { -2000.1 },
                new Object[] { -10.0000000000001 },
                new Object[] { -0.9999999999999 }
        };
    }
}
