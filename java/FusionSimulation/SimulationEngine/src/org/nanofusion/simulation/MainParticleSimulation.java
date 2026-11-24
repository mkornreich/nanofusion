/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.simulation;

import com.google.gson.Gson;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.List;

/**
 *
 * @author pkornrei
 */
public class MainParticleSimulation {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        /*double mu = 1.660538921E-27;
        double dm = 2*mu;
        double lm = 7*mu;
        double dc = 1.0;
        double lc = 3.0;
        double ke = 8.9875517873681765E+9;
        double r = 2.04E-10;
        double uc = 1.60217657E-19;
        double num = dc*lc*uc*uc*ke;
        double den = r*r;
        System.out.println(num);
        System.out.println(den);
        double f = num/den;
        System.out.println(f);
        double da = f/dm;
        double la = f/lm;
        System.out.println("Ad = " + da);
        System.out.println("Al = " + la);*/
        /*Simulator simulator = new Simulator();
        List<Lattice> simulation = simulator.simulate(new SimulatorProperties());
        Gson gson = new Gson();
        String result = gson.toJson(simulation);
        System.out.println(result);*/
        
        /*for (Lattice lattice: simulation) {
            List<Particle> particles = lattice.getParticles();
            Particle corner = particles.get(0);
            PVector position = corner.getPosition();
            System.out.println(position.getX().toString());
        }
        for (Lattice lattice: simulation) {
            List<Particle> particles = lattice.getParticles();
            Particle corner = particles.get(1);
            PVector position = corner.getPosition();
            System.out.println(position.getX().toString());
        }*/
        System.exit(0);
    }
}
