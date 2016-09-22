package vrr.jargon.exponential;

import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

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
        final double expected = rCall("qbeta")
                .with(p, alpha, beta)
                .evalAsDouble();
        final double actual = Beta.inverseCDF(p, alpha, beta);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(expected));
        else
            Assert.assertEquals(actual, expected, Math.abs(expected * 0.0000000001));
    }

    @Test(dataProvider= "inverseCdfData")
    public void testLogInverseCDF(final double alpha, final double beta, final double p) throws REXPMismatchException, REngineException {
        final double expected = rCall("qbeta")
                .with(Math.log(p), alpha, beta)
                .with("log", true)
                .evalAsDouble();
        final double actual = Beta.logInverseCDF(Math.log(p), alpha, beta);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(expected));
        else
            Assert.assertEquals(actual, expected, Math.abs(expected * 0.0000000001));
    }

    @DataProvider(name="inverseCdfData")
    public Object[][] inverseCdfData() throws IOException {
        final List<Object[]> result = new ArrayList<>();
        //final File outputFile = createTempFile("R-output",".tab");
        final File outputFile = File.createTempFile("R-output",".tab");
        //final File scriptFile = createTempFile("R-script", ".R");
        final File scriptFile = File.createTempFile("R-script", ".R");

        final double[] alphas = new double[] { 0.01, 0.1, 0.5, 1, 1.5, 3, 10, 30};
        final double[] betas = new double[] { 0.01, 0.1, 0.5, 1, 1.5, 3, 10, 30, 300};
        final double[] probs = new double[] { 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 0.9999, 0.99999999999999, 1 };
        for (double a: alphas)
            for (double b: betas)
                for (double p: probs) {
                    result.add(new Object[] { a, b, p});
                }
        return result.toArray(new Object[result.size()][]);
    }
}
