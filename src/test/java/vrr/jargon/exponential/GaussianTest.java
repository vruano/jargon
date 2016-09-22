package vrr.jargon.exponential;

import org.rosuda.REngine.*;
import org.rosuda.REngine.JRI.JRIEngine;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Unit tests for {@link Gaussian}
 */

public class GaussianTest {


    private static JRIEngine R;

    @BeforeClass
    public void setUp() throws REngineException {
        R = new JRIEngine();
    }

    @AfterClass
    public void tearDown() {
        R.close();
    }

    @Test(dataProvider = "densityData")
    public void testDensity(final double mu, final double sigma, final double x) throws REngineException, REXPMismatchException {
        final REXP dnormCall = REXP.asCall("dnorm", new REXP[] { new REXPDouble(x), new REXPDouble(mu), new REXPDouble(sigma), new REXPLogical(false)});
        final double expected = R.eval(dnormCall, null, true).asDouble();
        final double actual = Gaussian.density(x, mu, sigma, false);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, 0.000000001);
    }

    @Test(dataProvider = "densityData")
    public void testLogDensity(final double mu, final double sigma, final double x) throws REngineException, REXPMismatchException {
        final REXP dnormCall = REXP.asCall("dnorm", new REXP[] { new REXPDouble(x), new REXPDouble(mu), new REXPDouble(sigma), new REXPLogical(true)});
        final double expected = R.eval(dnormCall, null, true).asDouble();
        final double actual = Gaussian.density(x, mu, sigma, true);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(0.000000001 * expected));
    }


    @Test(dataProvider = "densityData")
    public void testCDF(final double mu, final double sigma, final double x) throws REXPMismatchException, REngineException {
        final REXP pnormCall = REXP.asCall("pnorm", new REXP[] { new REXPDouble(x), new REXPDouble(mu), new REXPDouble(sigma)});
        final double expected = R.eval(pnormCall, null, true).asDouble();
        final double actual = Gaussian.CDF(x, mu, sigma, false);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(0.000000001 * expected));
    }

    @Test(dataProvider = "densityData")
    public void testLogCDF(final double mu, final double sigma, final double x) throws REXPMismatchException, REngineException {
        final REXP pnormCall = REXP.asCall("pnorm", new REXP[] { new REXPDouble(x), new REXPDouble(mu), new REXPDouble(sigma),
           new REXPLogical(true), new REXPLogical(true) });
        final double expected = R.eval(pnormCall, null, true).asDouble();
        final double actual = Gaussian.CDF(x, mu, sigma, true);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(0.000000001 * expected));
    }

    @Test(dataProvider = "inverseCDFData")
    public void testInverseCDF(final double mu, final double sigma, final double p) throws REXPMismatchException, REngineException {
        final REXP qnormCall = REXP.asCall("qnorm", new REXP[] { new REXPDouble(p), new REXPDouble(mu), new  REXPDouble(sigma)});
        final double expected = R.eval(qnormCall, null, true).asDouble();
        final double actual = Gaussian.inverseCDF(p, mu, sigma, false);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(0.000000001 * expected));
    }

    @Test(dataProvider = "inverseCDFData")
    public void testLogInverseCDF(final double mu, final double sigma, final double p) throws REXPMismatchException, REngineException {
        final REXP qnormCall = REXP.asCall("qnorm", new REXP[] { new REXPDouble(Math.log(p)), new REXPDouble(mu), new  REXPDouble(sigma),
           new REXPLogical(true), new REXPLogical(true)});
        final double expected = R.eval(qnormCall, null, true).asDouble();
        final double actual = Gaussian.inverseCDF(Math.log(p), mu, sigma, true);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, 0.000000001);
    }

    @DataProvider(name = "densityData")
    public Object[][] densityData() {
        final List<Object[]> result = new ArrayList<>();
        result.add(new Object[] { 0.0, 1.0, 0.5});
        result.add(new Object[] { 0.0, 1.0, Double.NaN });
        result.add(new Object[] { 0.0, Double.NaN, 0.5 });
        result.add(new Object[] { Double.NaN, 1.0, 0.5 });
        result.add(new Object[] { 0.0, 1.0, Double.POSITIVE_INFINITY });
        result.add(new Object[] { 0.0, Double.POSITIVE_INFINITY, 0.5 });
        result.add(new Object[] { Double.POSITIVE_INFINITY, 1.0, 0.5 });
        result.add(new Object[] { 0.0, 1.0, Double.NEGATIVE_INFINITY });
        result.add(new Object[] { 0.0, Double.NEGATIVE_INFINITY, 0.5 });
        result.add(new Object[] { Double.NEGATIVE_INFINITY, 1.0, 0.5 });
        result.add(new Object[] { 0.0, 0.0, 0.5 });
        result.add(new Object[] { 0.0, 0.0, 0.0 });
        result.add(new Object[] { 0.0, Double.POSITIVE_INFINITY, 0.5 });

        final double[] mus = new double[] { 0.1 , 10.0, 0.000001, -1.3, -1000212.1, -0.000001 };
        final double[] sigmas = new double[] { 0.0, Double.POSITIVE_INFINITY, -0.1, 100, 0.000001, 2.0 };
        final double[] zscores = new double[] { 0.0 , -1.0, -3.5, -100.0, 1.1, 0.00001, 10. };

        for (final double mu : mus) {
            for (final double sigma : sigmas) {
                for (final double zscore : zscores) {
                    result.add(new Object[] { mu, sigma, sigma * zscore });
                }
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name = "inverseCDFData")
    public Object[][] inverseCDFData() {
        final List<Object[]> result = new ArrayList<>();
        result.add(new Object[] { 0.0, 1.0, 0.25});
        result.add(new Object[] { 0.0, 1.0, Double.NaN });
        result.add(new Object[] { 0.0, Double.NaN, 0.25 });
        result.add(new Object[] { Double.NaN, 1.0, 0.25 });
        result.add(new Object[] { 0.0, 1.0, Double.POSITIVE_INFINITY });
        result.add(new Object[] { 0.0, Double.POSITIVE_INFINITY, 0.25 });
        result.add(new Object[] { Double.POSITIVE_INFINITY, 1.0, 0.25 });
        result.add(new Object[] { 0.0, 1.0, Double.NEGATIVE_INFINITY });
        result.add(new Object[] { 0.0, Double.NEGATIVE_INFINITY, 0.25 });
        result.add(new Object[] { Double.NEGATIVE_INFINITY, 1.0, 0.25 });
        result.add(new Object[] { 0.0, 0.0, 0.25 });
        result.add(new Object[] { 0.0, 0.0, 0.0 });
        result.add(new Object[] { 0.0, Double.POSITIVE_INFINITY, 0.25 });

        final double[] mus = new double[] { 0.1 , 10.0, 0.000001, -1.3, -1000212.1, -0.000001 };
        final double[] sigmas = new double[] { 0.0, Double.POSITIVE_INFINITY, -0.1, 100, 0.000001, 2.0 };
        final double[] p = new double[] { 0.0 , 1.0, .5, 1e-100, 1e-10, 5.2e-2, 0.9999999999, Double.POSITIVE_INFINITY, 1.1 };

        for (final double mu : mus) {
            for (final double sigma : sigmas) {
                for (final double zscore : p) {
                    result.add(new Object[] { mu, sigma, sigma * zscore });
                }
            }
        }
        return result.toArray(new Object[result.size()][]);
    }


}
