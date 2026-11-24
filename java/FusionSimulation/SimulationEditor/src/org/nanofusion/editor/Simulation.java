 /*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.editor;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.nanofusion.simulation.FusionPair;
import org.nanofusion.simulation.Lattice;
import org.nanofusion.simulation.SimulatorProperties;

/**
 * 
 *    Main bean of the NanoFusion Simulation program that contains the parameters
 * used during and the results from running a simulation execution.
 *
 * @author pkornrei
 */
public class Simulation {

    private List<Lattice> results = new ArrayList<Lattice>();
    private SimulatorProperties properties = null;
    private transient boolean changed = false;
    private transient Map<FusionPair, String> pairs;
    private transient List<Double> minDistances = new ArrayList<Double>();

    /**
     * 
     *    Retrieves the results of running the simulation.
     * 
     * @return 
     */
    public List<Lattice> getResults() {
        return results;
    }

    /**
     * 
     *    Sets the results of running the simulation.
     * 
     * @param results 
     */
    public void setResults(List<Lattice> results) {
        this.results = results;
    }

    /**
     *     Retrieves the parameters used for the results.
     * 
     * @return 
     */
    public SimulatorProperties getProperties() {
        return properties;
    }

    /**
     *      Sets the parameters used for the results.
     * 
     * @param properties 
     */
    public void setProperties(SimulatorProperties properties) {
        this.properties = properties;
    }

    public Map<FusionPair, String> getPairs() {
        return pairs;
    }

    public void setPairs(Map<FusionPair, String> pairs) {
        this.pairs = pairs;
    }

    public List<Double> getMinDistances() {
        return minDistances;
    }

    public void setMinDistances(List<Double> minDistances) {
        this.minDistances = minDistances;
    }

    /**
     *   Retrieves whether the represented fusion results has changed.
     * 
     * @return 
     */
    public boolean isChanged() {
        return changed;
    }

    /**
     * 
     *    Retrieves whether the represented fusion results has changed.
     * 
     * @param changed 
     */
    public void setChanged(boolean changed) {
        this.changed = changed;
    }

    @Override
    public int hashCode() {
        if (properties!=null) {
            int hash = properties.getIterations() + properties.getLatticeSize() + properties.getPrecision();
            if (results!=null) {
                hash += results.size();
            }
            return hash;
        }
        return -1;
    }
    
    @Override
    public boolean equals(Object o) {
        if (o instanceof Simulation) {
            if (o!=null) {
                Simulation simulation = (Simulation)o;
                if (this.getProperties().equals(simulation.getProperties())) {
                    return true;
                }
            }
        }
        return false;
    }
    
}
