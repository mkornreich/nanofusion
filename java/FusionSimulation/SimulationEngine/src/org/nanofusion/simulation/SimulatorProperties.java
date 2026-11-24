/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.simulation;

import org.nanofusion.simulation.PVector;
import java.math.BigDecimal;

/**
 * 
 *     Bean that contains the properties for running a simulation.
 *
 * @author pkornrei
 */
public class SimulatorProperties {
    
    private String name = "Simulation 1";
    
    private int latticeXSize = 10;
    
    private int latticeYSize = 10;
    
    private int latticeZSize = 10;
    
    private int iterations = 5;
    
    private BigDecimal timeStep = new BigDecimal("2.0E-17");
    
    private int precision = 6;
    
    private PVector electricField = null;
    
    private BigDecimal threshold = new BigDecimal(2.0E-15);
    
    private boolean includeLorentz = false;
    
    private boolean useLogTime = false;
    
    private int transitionStep = 0;
    
    private BigDecimal dampingFactor = null;

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public int getLatticeXSize() {
        return latticeXSize;
    }

    public void setLatticeXSize(int latticeXSize) {
        this.latticeXSize = latticeXSize;
    }

    public int getLatticeYSize() {
        return latticeYSize;
    }

    public void setLatticeYSize(int latticeYSize) {
        this.latticeYSize = latticeYSize;
    }

    public int getLatticeZSize() {
        return latticeZSize;
    }

    public void setLatticeZSize(int latticeZSize) {
        this.latticeZSize = latticeZSize;
    }

    public int getLatticeSize() {
        return latticeXSize + latticeYSize + latticeZSize;
    }

    public PVector getElectricField() {
        return electricField;
    }

    public void setElectricField(PVector electricField) {
        this.electricField = electricField;
    }

    public int getIterations() {
        return iterations;
    }

    public void setIterations(int iterations) {
        this.iterations = iterations;
    }

    public BigDecimal getTimeStep() {
        return timeStep;
    }

    public void setTimeStep(BigDecimal timeStep) {
        this.timeStep = timeStep;
    }

    public int getPrecision() {
        return precision;
    }

    public void setPrecision(int precision) {
        this.precision = precision;
    }

    public BigDecimal getThreshold() {
        return threshold;
    }

    public void setThreshold(BigDecimal threshold) {
        this.threshold = threshold;
    }

    public boolean isIncludeLorentz() {
        return includeLorentz;
    }

    public void setIncludeLorentz(boolean includeLorentz) {
        this.includeLorentz = includeLorentz;
    }

    public boolean isUseLogTime() {
        return useLogTime;
    }

    public void setUseLogTime(boolean useLogTime) {
        this.useLogTime = useLogTime;
    }

    public int getTransitionStep() {
        return transitionStep;
    }

    public void setTransitionStep(int transitionStep) {
        this.transitionStep = transitionStep;
    }

    public BigDecimal getDampingFactor() {
        return dampingFactor;
    }

    public void setDampingFactor(BigDecimal dampingFactor) {
        this.dampingFactor = dampingFactor;
    }
    
    
    
    @Override
    public boolean equals (Object o) {
        if (o instanceof SimulatorProperties) {
            if (o != null ) {
                SimulatorProperties p = (SimulatorProperties)o;
                if (this.latticeXSize == p.latticeXSize && this.latticeYSize == p.latticeYSize && this.latticeZSize == p.latticeZSize
                        && this.iterations == p.iterations && this.precision == p.precision && this.name.equals(p.name)) {
                    return true;
                }
            }
        }
        return false;
    }
    
    
}
