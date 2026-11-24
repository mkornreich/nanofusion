/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.simulation;

import java.io.Serializable;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author pkornrei
 */
public class Lattice implements Serializable {
    
    private List<Particle> particles = new ArrayList<Particle>();
    private long executionTime = 0;
    
    private BigDecimal time = BigDecimal.ZERO;
    
    public Lattice() {
    }
    
    public Lattice (LVector size, int precision) {
        this(size, precision, null);
    }

    
    public Lattice (LVector size, int precision, PVector electricField) {
        boolean deut = true;
        BigDecimal mass = null;
        BigDecimal invMass = null;
        BigDecimal charge = null;
        BigDecimal z = BigDecimal.ZERO;
        int id = 0;
        String type = "";
        for (int layer=0; layer<size.getZ(); layer++) {
            BigDecimal y = BigDecimal.ZERO;
            deut = !deut;
            for (int column=0; column<size.getY(); column++) {
                BigDecimal x = BigDecimal.ZERO;
                deut = !deut;
                for (int row=0; row<size.getX(); row++) {
                    if (deut) {
                        mass = Constants.DEUTERIUM_MASS;
                        invMass = Constants.INVERSE_DEUTERIUM_MASS;
                        charge = Constants.DEUTERIUM_CHARGE;
                        deut = false;
                        type = "D";
                    }
                    else {
                        mass = Constants.LITHIUM_MASS;
                        invMass = Constants.INVERSE_LITHIUM_MASS;
                        charge = Constants.LITHIUM_CHARGE;
                        deut = true;
                        type = "Li";
                    }
                    PVector position = new PVector(x, y, z);
                    LVector lvector = new LVector(row, column, layer);
                    Particle particle = new Particle(position, charge, electricField);
                    particle.setMass(mass);
                    particle.setInverseMass(invMass);
                    particle.setLatticePosition(lvector);
                    particle.setId(id);
                    particle.setType(type);
                    particle.setPrecision(precision);
                    particles.add(particle);
                    x = x.add(Constants.INITIAL_LATTICE_SPACING);
                    id++;
                }
                y = y.add(Constants.INITIAL_LATTICE_SPACING);
            }
            z = z.add(Constants.INITIAL_LATTICE_SPACING);
        }
    }

    public BigDecimal getTime() {
        return time;
    }

    public void setTime(BigDecimal time) {
        this.time = time;
    }
    
    

    public List<Particle> getParticles() {
        return particles;
    }

    public void setParticles(List<Particle> particles) {
        this.particles = particles;
    }

    public long getExecutionTime() {
        return executionTime;
    }

    public void setExecutionTime(long executionTime) {
        this.executionTime = executionTime;
    }
    
    

    @Override
    public Lattice clone() {
        Lattice clola = new Lattice();
        List<Particle> cloneParticles = new ArrayList<Particle>();
        for (Particle p: particles) {
            cloneParticles.add(p.clone());
        }
        clola.setParticles(cloneParticles);
        return clola;
    }
    
    
    
}
