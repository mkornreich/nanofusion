import com.google.gson.Gson;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.nanofusion.simulation.Lattice;
import org.nanofusion.simulation.PVector;
import org.nanofusion.simulation.Simulator;
import org.nanofusion.simulation.SimulatorProperties;

/**
 * Runs the original Java SimulationEngine headlessly and writes the result in
 * the same format as the original SaveProgressDialog. Used only to produce
 * reference output for checking the Python port.
 *
 * Usage: java ReferenceRun out.json key=value ...
 * Keys: name x y z iterations precision timeStep ex ey ez threshold lorentz
 *       logTime transition damping
 */
public class ReferenceRun {

    public static void main(String[] args) throws Exception {
        Map<String, String> a = new HashMap<String, String>();
        for (int i = 1; i < args.length; i++) {
            String[] kv = args[i].split("=", 2);
            a.put(kv[0], kv.length > 1 ? kv[1] : "");
        }
        SimulatorProperties p = new SimulatorProperties();
        if (a.containsKey("name")) p.setName(a.get("name"));
        if (a.containsKey("x")) p.setLatticeXSize(Integer.parseInt(a.get("x")));
        if (a.containsKey("y")) p.setLatticeYSize(Integer.parseInt(a.get("y")));
        if (a.containsKey("z")) p.setLatticeZSize(Integer.parseInt(a.get("z")));
        if (a.containsKey("iterations")) p.setIterations(Integer.parseInt(a.get("iterations")));
        if (a.containsKey("precision")) p.setPrecision(Integer.parseInt(a.get("precision")));
        if (a.containsKey("timeStep")) p.setTimeStep(new BigDecimal(a.get("timeStep")));
        if (a.containsKey("threshold")) p.setThreshold(new BigDecimal(a.get("threshold")));
        if (a.containsKey("ex") || a.containsKey("ey") || a.containsKey("ez")) {
            PVector e = new PVector();
            if (a.containsKey("ex")) e.setX(new BigDecimal(a.get("ex")));
            if (a.containsKey("ey")) e.setY(new BigDecimal(a.get("ey")));
            if (a.containsKey("ez")) e.setZ(new BigDecimal(a.get("ez")));
            p.setElectricField(e);
        }
        p.setIncludeLorentz("true".equals(a.get("lorentz")));
        p.setUseLogTime("true".equals(a.get("logTime")));
        if (a.containsKey("transition")) p.setTransitionStep(Integer.parseInt(a.get("transition")));
        if (a.containsKey("damping")) p.setDampingFactor(new BigDecimal(a.get("damping")));

        // Same loop as SimulationProgressDialog.ProgressThread.run()
        Simulator simulator = new Simulator(p);
        List<Lattice> result = new ArrayList<Lattice>();
        result.add(simulator.getLattice());
        for (int ii = 0; ii < p.getIterations(); ii++) {
            result.add(simulator.step(ii));
        }

        // Same output as SaveProgressDialog.ProgressThread.run()
        BufferedWriter writer = new BufferedWriter(new FileWriter(args[0]));
        Gson gson = new Gson();
        writer.write("{\n");
        writer.write("properties: ");
        gson.toJson(p, writer);
        writer.write(",\n");
        writer.write("results:\n");
        writer.write("[\n");
        for (Lattice lattice : result) {
            gson.toJson(lattice, writer);
            writer.write(",\n");
        }
        writer.write("]\n");
        writer.write("}");
        writer.close();
        boolean fused = simulator.getPairs() != null && simulator.getPairs().size() > 0;
        System.err.println("fusion=" + fused);
        System.exit(0);
    }
}
