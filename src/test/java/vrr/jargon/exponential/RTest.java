package vrr.jargon.exponential;

import org.apache.commons.math.util.DoubleArray;
import org.rosuda.REngine.*;
import org.rosuda.REngine.JRI.JRIEngine;
import org.testng.Assert;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

/**
 * Created by valentin on 9/21/16.
 */
public class RTest {
    protected static JRIEngine R;

    public synchronized void setUp() throws REngineException {
        if (R == null) {
            R = new JRIEngine();
        }
    }

    public synchronized void tearDown() {
        //R.close(); R = null;
    }

    protected RCall rCall(final String method) {
        return new RCall(method);
    }

    protected void assertEquals(final double actual, final double expected, final String message) {
        if (Double.isNaN(actual))
            Assert.assertTrue(Double.isNaN(expected), message);
        else
            Assert.assertEquals(actual, expected, Math.abs(actual * epsilon()), message);
    }

    protected void assertEquals(final double actual, final double expected) {
        if (Double.isNaN(actual))
            Assert.assertTrue(Double.isNaN(expected), null);
        else
            Assert.assertEquals(actual, expected, Math.abs(actual * epsilon()), null);
    }

    protected double epsilon() {
        return 1e-10;
    }

    public static class RCall {
        private final String method;

        private final List<REXP> parameters;

        private final List<String> parameterNames;

        protected RCall(final String method) {
            this.method = Objects.requireNonNull(method);
            parameterNames = new ArrayList<>();
            parameters = new ArrayList<>();
        }

        public RCall withArray(final double ... values) {
            parameters.add(new REXPDouble(values));
            parameterNames.add(null);
            return this;
        }

        public RCall with(final double ... values) {
            for (final double v : values) {
                parameters.add(new REXPDouble(v));
                parameterNames.add(null);
            }
            return this;
        }

        public RCall with(final boolean ... values) {
            for (final boolean v : values) {
                parameters.add(new REXPLogical(v));
                parameterNames.add(null);
            }
            return this;
        }

        public RCall with(final String name, final double v) {
            parameterNames.add(name);
            parameters.add(new REXPDouble(v));
            return this;
        }

        private REXP eval() throws REXPMismatchException, REngineException {
            final REXP[] callArray = new REXP[parameters.size() + 1];
            callArray[0] = new REXPSymbol(method);
            for (int i = 0; i < parameters.size(); i++) {
                callArray[i + 1] = parameters.get(i);
            }
            final RList callItems = new RList(callArray);
            for (int i = 0; i < parameters.size(); i++) {
                final String name = parameterNames.get(i);
                if (name != null) {
                    callItems.setKeyAt(i + 1, name);
                }
            }
            final REXP call = new REXPLanguage(callItems);
            return R.eval(call, null, true);
        }

        public double evalAsDouble() throws REXPMismatchException, REngineException {
            return eval().asDouble();
        }

        public double[] evalAsDoubleArray() throws REXPMismatchException, REngineException {
            return eval().asDoubles();
        }

        public String toString() {
            final REXP[] callArray = new REXP[parameters.size() + 1];
            callArray[0] = new REXPSymbol(method);
            for (int i = 0; i < parameters.size(); i++) {
                callArray[i + 1] = parameters.get(i);
            }
            final RList callItems = new RList(callArray);
            for (int i = 0; i < parameters.size(); i++) {
                final String name = parameterNames.get(i);
                if (name != null) {
                    callItems.setKeyAt(i + 1, name);
                }
            }
            final REXP call = new REXPLanguage(callItems);
            return call.toDebugString();
        }

        public RCall with(final String name, final boolean v) {
            parameterNames.add(name);
            parameters.add(new REXPLogical(v));
            return this;
        }

    }
}
