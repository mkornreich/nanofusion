/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.simulation;

import java.math.BigDecimal;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 *
 * @author pkornrei
 */
public class Simulator {
    
    private Lattice past = null;
    private Lattice lattice = null;
    private SimulatorProperties props;
    private BigDecimal time = BigDecimal.ZERO;
    private Map<FusionPair,String> pairs = null;
    private BigDecimal minDistance = BigDecimal.ONE;
    private Date executionStart = null;
    private int threadsReturned = 0;

    public Simulator(SimulatorProperties props) {
        executionStart = new Date();
        this.props = props;
        LVector size = new LVector(props.getLatticeXSize(), props.getLatticeYSize(), props.getLatticeZSize());
        past = new Lattice(size, props.getPrecision(), props.getElectricField());
        lattice = past;
        lattice.setExecutionTime(0);
    }
    
    private BigDecimal getTimeStep(int step) {
      BigDecimal timeStep = props.getTimeStep();
      if (props.isUseLogTime()) {
          double expnt = 1.*(step - props.getTransitionStep())/props.getDampingFactor().doubleValue();
          BigDecimal decay = new BigDecimal(1.0/(1.0 + Math.exp(expnt)));
          timeStep = timeStep.multiply(decay);
      }
      return timeStep;
    };
    
    private void reportReturn() {
        threadsReturned++;
    }

    
    public Lattice step(int iteration) {
        BigDecimal timeStep = getTimeStep(iteration);
        time=time.add(timeStep);
        lattice = null;
        lattice = past.clone();
        lattice.setTime(time);
        List<Particle> particles = lattice.getParticles();
        // TODO Calculate force;
        int p = 0;
        int maxThreads = 4;
        ExecutorService executor = Executors.newFixedThreadPool(maxThreads);
        while (p<particles.size()) {
            ParticleThread thread = new ParticleThread();
            thread.setParticleIndex(p);
            executor.execute(thread);
            p++;
        }
        executor.shutdown();
        while (!executor.isTerminated()) {}
        
        /*for (Particle particle: particles) {
            PVector force = new PVector();
            PVector rkd1 = new PVector();
            PVector rkd2 = new PVector();
            PVector rkd3 = new PVector();
            PVector rkd4 = new PVector();
            PVector rkv1 = new PVector();
            PVector rkv2 = new PVector();
            PVector rkv3 = new PVector();
            PVector rkv4 = new PVector();
            for (Particle actor: particles) {
                if (particle.getLatticePosition().isNearNeighbor(actor.getLatticePosition())) {
                    ParticleInteraction xaction = CommonCalculatorUtility.getInteraction(actor, particle, props.getPrecision());
                    //ParticleInteraction xaction = CommonCalculatorUtility.calculateRungeKuttaIncrements(actor, particle, props.getPrecision(), props.getTimeStep());
                    force.add(xaction.getForce());
                    //rkd1.add(xaction.getRkd1());
                    //rkd2.add(xaction.getRkd2());
                    //rkd3.add(xaction.getRkd3());
                    //rkd4.add(xaction.getRkd4());
                    //rkv1.add(xaction.getRkv1());
                    //rkv2.add(xaction.getRkv2());
                    //rkv3.add(xaction.getRkv3());
                    //rkv4.add(xaction.getRkv4());
                    if (xaction.getDistance().subtract(props.getThreshold()).signum()<=0) {
                        if (pairs==null) {
                            pairs = new HashMap<FusionPair, String>();
                        }
                        FusionPair pair = new FusionPair();
                        pair.setParticle1(actor);
                        pair.setParticle2(particle);
                        pairs.put(pair,"Threshold");
                    }
                    LVector partlv = particle.getLatticePosition();
                    LVector actorlv = actor.getLatticePosition();
                    PVector partpos = particle.getPosition();
                    PVector actpos = actor.getPosition();
                    // Particles have passed
                    //if (partlv.getX()<actorlv.getX() && partlv.getY()==actorlv.getY() && partlv.getZ()==actorlv.getZ() && (partpos.getX().subtract(actpos.getX()).signum()>0)) {
                    //    if (pairs==null) {
                    //        pairs = new HashMap<FusionPair, String>();
                    //    }
                    //    FusionPair pair = new FusionPair();
                    //    pair.setParticle1(actor);
                    //    pair.setParticle2(particle);
                    //    pairs.put(pair, "Passed");
                    //}
                }
            }
            particle.setCurrentForce(force, time, props.getPrecision(), props.isIncludeLorentz());
            //particle.setRungeKuttaParameters(rkd1, rkd2, rkd3, rkd4, rkv1, rkv2, rkv3, rkv4);
        }*/
        for (Particle particle: particles) {
            particle.update(timeStep);
            //particle.updateRungeKutta();
        }
        past = null;
        past = lattice;
        long now = (new Date()).getTime();
        lattice.setExecutionTime(now  - executionStart.getTime());
        return lattice;
    }

    public Lattice getLattice() {
        return lattice;
    }

    public void setLattice(Lattice lattice) {
        this.lattice = lattice;
    }

    public Map<FusionPair, String> getPairs() {
        return pairs;
    }

    public BigDecimal getMinDistance() {
        return minDistance;
    }

    private class ParticleThread  extends Thread {
        
        private int particleIndex = -1;

        public int getParticleIndex() {
            return particleIndex;
        }

        public void setParticleIndex(int particleIndex) {
            this.particleIndex = particleIndex;
        }

        public void run() {
            List<Particle> particles = lattice.getParticles();
            Particle particle = particles.get(particleIndex);
            PVector force = new PVector();
            for (Particle actor: particles) {
                if (particle.getLatticePosition().isNearNeighbor(actor.getLatticePosition())) {
                    ParticleInteraction xaction = CommonCalculatorUtility.getInteraction(actor, particle, props.getPrecision());
                    //ParticleInteraction xaction = CommonCalculatorUtility.calculateRungeKuttaIncrements(actor, particle, props.getPrecision(), props.getTimeStep());
                    force.add(xaction.getForce());
                    //rkd1.add(xaction.getRkd1());
                    //rkd2.add(xaction.getRkd2());
                    //rkd3.add(xaction.getRkd3());
                    //rkd4.add(xaction.getRkd4());
                    //rkv1.add(xaction.getRkv1());
                    //rkv2.add(xaction.getRkv2());
                    //rkv3.add(xaction.getRkv3());
                    //rkv4.add(xaction.getRkv4());
                    if (xaction.getDistance().subtract(props.getThreshold()).signum()<=0) {
                        if (pairs==null) {
                            pairs = new HashMap<FusionPair, String>();
                        }
                        FusionPair pair = new FusionPair();
                        pair.setParticle1(actor);
                        pair.setParticle2(particle);
                        pairs.put(pair,"Threshold");
                    }
                    LVector partlv = particle.getLatticePosition();
                    LVector actorlv = actor.getLatticePosition();
                    PVector partpos = particle.getPosition();
                    PVector actpos = actor.getPosition();
                    // Particles have passed
                    //if (partlv.getX()<actorlv.getX() && partlv.getY()==actorlv.getY() && partlv.getZ()==actorlv.getZ() && (partpos.getX().subtract(actpos.getX()).signum()>0)) {
                    //    if (pairs==null) {
                    //        pairs = new HashMap<FusionPair, String>();
                    //    }
                    //    FusionPair pair = new FusionPair();
                    //    pair.setParticle1(actor);
                    //    pair.setParticle2(particle);
                    //    pairs.put(pair, "Passed");
                    //}
                }
            }
            particle.setCurrentForce(force, time, props.getPrecision(), props.isIncludeLorentz());
            reportReturn();
        }
    }
    
}
