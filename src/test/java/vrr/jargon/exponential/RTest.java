package vrr.jargon.exponential;

import org.rosuda.REngine.*;
import org.rosuda.REngine.JRI.JRIEngine;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

/**
 * Created by valentin on 9/21/16.
 */
public class RTest {
    protected JRIEngine R;

    @BeforeClass
    public void setUp() throws REngineException {
        R = new JRIEngine();
    }

    @AfterClass
    public void tearDown() {
        R.close();
    }

    protected RCall rCall(final String method) {
        return new RCall(method);
    }

    protected class RCall {
        private final String method;

        private final List<REXP> parameters;

        private final List<String> parameterNames;

        protected RCall(final String method) {
            this.method = Objects.requireNonNull(method);
            parameterNames = new ArrayList<>();
            parameters = new ArrayList<>();
        }

        public RCall with(final double ... values) {
            for (final double v : values) {
                parameters.add(new REXPDouble(v));
                parameterNames.add(null);
            }
            return this;
        }

        public RCall with(final String name, final double v) {
            parameterNames.add(name);
            parameters.add(new REXPDouble(v));
            return this;
        }

        public double evalAsDouble() throws REXPMismatchException, REngineException {
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
            return R.eval(call, null, true).asDouble();
        }


        public RCall with(final String name, final boolean v) {
            parameterNames.add(name);
            parameters.add(new REXPLogical(v));
            return this;
        }
    }
}
