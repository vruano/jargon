package vrr.jargon.exponential;

import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Unit test cases for {@link Poisson}
 */
public class PoissonTest extends RTest {

    @BeforeClass
    public void setUp() throws REngineException {
        super.setUp();
    }

    @AfterClass
    public void tearDown() {
        super.tearDown();
    }

    @Test(dataProvider = "densityData")
    public void testDensity(final double lambda, final long x) throws REXPMismatchException, REngineException {
        final double expected = rCall("dpois")
                .with(x, lambda)
                .evalAsDouble();
        final double actual = Poisson.density(x, lambda);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(expected * 0.000000001));
    }

    @Test(dataProvider = "densityData")
    public void testLogDensity(final double lambda, final double x) throws REXPMismatchException, REngineException {
        final double expected = rCall("dpois")
                .with(x, lambda)
                .with("log", true)
                .evalAsDouble();
        final double actual = Poisson.logDensity(x, lambda);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual, expected, Math.abs(expected * 0.000000001));
    }

    @DataProvider(name = "densityData")
    public Object[][] densityData() {
        final List<Object[]> result = new ArrayList<>();
        final double[] lambdas = new double[] {
                Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
                Double.MIN_NORMAL , Double.MIN_VALUE, Double.MAX_VALUE, 0, 0.00001, 1.1, 1.2, 2.0, 3.0, 10.21, -12, 10000.212, 2021210211.1};
        final long[] values = new long[] { 0, -1, 10, 13, Long.MAX_VALUE };
        for (final double lambda : lambdas)
            for (final long x : values)
                result.add(new Object[] { lambda, x});
        return result.toArray(new Object[result.size()][]);
    }
}
