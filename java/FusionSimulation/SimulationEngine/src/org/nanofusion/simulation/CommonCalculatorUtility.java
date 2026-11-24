/*
 * BS"D
 */
package org.nanofusion.simulation;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

/**
 * 
 *    Provides utility numerical calculations for the simulation.
 *
 * @author pkornrei
 */
public class CommonCalculatorUtility {
    
    /**
     *    Calculates the square root of a BigDecimal encoded number
     * with a given scale.
     * 
     * @param in
     * @param scale
     * @return 
     */
    public static BigDecimal sqrt(BigDecimal in, int scale){
        BigDecimal sqrt = new BigDecimal(1);
        sqrt.setScale(scale + 3, RoundingMode.FLOOR);
        BigDecimal store = new BigDecimal(in.toString());
        boolean first = true;
        do{
            if (!first){
                store = new BigDecimal(sqrt.toString());
            }
            else first = false;
            store.setScale(scale + 3, RoundingMode.FLOOR);
            sqrt = in.divide(store, scale + 3, RoundingMode.FLOOR).add(store).divide(
                    BigDecimal.valueOf(2), scale + 3, RoundingMode.FLOOR);
        }
        while (!store.equals(sqrt));
        return sqrt.setScale(scale, RoundingMode.FLOOR);
    }
    
    /**
     * 
     *     Simple inexpensive square root calculation of a BigDecimal number.
     * Does not work for values near 0, but for this simulation, since the
     * square root is used for calculating distance, it does not provide such a problem.
     * 
     * 
     * @param value
     * @return 
     */
    public static BigDecimal sqrt(BigDecimal value) {
        BigDecimal x = new BigDecimal(Math.sqrt(value.doubleValue()));
        return x.add(new BigDecimal(value.subtract(x.multiply(x)).doubleValue() / (x.doubleValue() * 2.0)));
        //return x;
    }
    
    /**
     * 
     *     Calculates the square magnitude of the difference between two vectors.
     * 
     * @param p1
     * @param p2
     * @param precision
     * @return 
     */
    public static BigDecimal getSquareDistance(PVector p1, PVector p2, int precision) {
        BigDecimal distance = BigDecimal.ZERO;
        MathContext mc = new MathContext(precision, RoundingMode.HALF_UP);
        BigDecimal x = p1.getX().subtract(p2.getX());
        BigDecimal y = p1.getY().subtract(p2.getY());
        BigDecimal z = p1.getZ().subtract(p2.getZ());
        distance = x.multiply(x, mc).add(y.multiply(y, mc)).add(z.multiply(z, mc));
        return distance;
    }
    
    /**
     *     Calculates the distance between the given vectors.
     * 
     * @param p1
     * @param p2
     * @param precision
     * @return 
     */
    public static BigDecimal getDistance(PVector p1, PVector p2, int precision) {
        return sqrt(getSquareDistance(p1, p2, precision));
    }
    
    /**
     *     Calculates the difference between the given vectors in the direction
     * from the 1st to the 2nd vector.
     * 
     * @param from
     * @param to
     * @return 
     */
    public static PVector getDirection(PVector from, PVector to) {
        return new PVector(to.getX().subtract(from.getX()), to.getY().subtract(from.getY()), to.getZ().subtract(from.getZ()));
    }
    
    public static PVector getForce(Particle p1, Particle p2, int precision) {
        MathContext mc = new MathContext(precision, RoundingMode.HALF_UP);
        BigDecimal potential = Constants.POTENTIAL.multiply(p1.getCharge(), mc).multiply(p2.getCharge(), mc);
        BigDecimal sqdist = getSquareDistance(p1.getPosition(), p2.getPosition(), precision);
        BigDecimal distance = sqrt(sqdist, precision);
        BigDecimal unitMagnitude = potential.divide(distance.multiply(sqdist, mc), mc);
        PVector force = getDirection(p1.getPosition(), p2.getPosition());
        force.multiply(unitMagnitude);
        //System.out.println("Coulomb Force: " + force.getX().toString());
        return force;
    }
    
    public static ParticleInteraction getInteraction(Particle p1, Particle p2, int precision) {
        MathContext mc = new MathContext(precision, RoundingMode.HALF_UP);
        BigDecimal potential = Constants.POTENTIAL.multiply(p1.getCharge(), mc).multiply(p2.getCharge(), mc);
        BigDecimal sqdist = getSquareDistance(p1.getPosition(), p2.getPosition(), precision);
        BigDecimal distance = sqrt(sqdist);
        BigDecimal unitMagnitude = potential.divide(distance.multiply(sqdist, mc), mc);
        PVector force = getDirection(p1.getPosition(), p2.getPosition());
        force.multiply(unitMagnitude);
        ParticleInteraction xaction = new ParticleInteraction();
        xaction.setDistance(distance);
        xaction.setForce(force);
        return xaction;
    }
    
    public static BigDecimal applyElapsedTimeChange(BigDecimal initial, BigDecimal rate, int precision, BigDecimal timeStep) {
        MathContext mc = new MathContext(precision, RoundingMode.HALF_UP);
        return initial.add(rate.multiply(timeStep), mc);
    }
    
    public static PVector getElapsdTimeChange(PVector rate, int precision, BigDecimal timeStep) {
        MathContext mc = new MathContext(precision, RoundingMode.HALF_UP);
        PVector result = rate.clone();
        result.multiply(timeStep);
        return result;
    }
    
    public static ParticleInteraction calculateRungeKuttaIncrements(Particle p1, Particle p2, 
            ParticleActor actor, BigDecimal time, BigDecimal timeStep, int precision) {
        ParticleInteraction xaction = new ParticleInteraction();
        PVector ptmp = p2.getPosition().clone();
        PVector vtmp = p2.getVelocity().clone();
        BigDecimal c1 = p1.getCharge();
        BigDecimal c2 = p2.getCharge();
        PVector pos2 = p1.getPosition().clone();
        /*System.out.println("Particle 1 " + p1.getLatticePosition().getX());
        System.out.println("Particle 2 " + p2.getLatticePosition().getX());
        System.out.println(ptmp.getX().toString());
        System.out.println(pos2.getX().toString());*/
        ParticleInteraction x1 = actor.getInteraction(pos2, ptmp, c1, c2, time, precision);
        xaction.setDistance(x1.getDistance());
        xaction.setRkd1(getElapsdTimeChange(vtmp, precision, timeStep));
        xaction.setRkv1(getElapsdTimeChange(x1.getForce(), precision, timeStep));
        BigDecimal half = new BigDecimal(0.5);
        ptmp = xaction.getRkd1().clone();
        ptmp.multiply(half);
        ptmp.add(p2.getPosition());
        vtmp = xaction.getRkv1().clone();
        vtmp.multiply(half);
        vtmp.add(p2.getVelocity());
        x1 = actor.getInteraction(pos2, ptmp, c1, c2, time, precision);
        xaction.setRkd2(getElapsdTimeChange(vtmp, precision, timeStep));
        xaction.setRkv2(getElapsdTimeChange(x1.getForce(), precision, timeStep));
        ptmp = xaction.getRkd2().clone();
        ptmp.multiply(half);
        ptmp.add(p2.getPosition());
        vtmp = xaction.getRkv2().clone();
        vtmp.multiply(half);
        vtmp.add(p2.getVelocity());
        x1 = actor.getInteraction(pos2, ptmp, c1, c2, time, precision);
        xaction.setRkd3(getElapsdTimeChange(vtmp, precision, timeStep));
        xaction.setRkv3(getElapsdTimeChange(x1.getForce(), precision, timeStep));
        ptmp = xaction.getRkd3().clone();
        ptmp.add(p2.getPosition());
        vtmp = xaction.getRkv3().clone();
        vtmp.add(p2.getVelocity());
        x1 = actor.getInteraction(pos2, ptmp, c1, c2, time, precision);
        xaction.setRkd4(getElapsdTimeChange(vtmp, precision, timeStep));
        xaction.setRkv4(getElapsdTimeChange(x1.getForce(), precision, timeStep));
        return xaction;
    }
    
    public static BigDecimal updateRungeKutta(BigDecimal initialValue, BigDecimal rk1, BigDecimal rk2, BigDecimal rk3, BigDecimal rk4, int precision) {
        //System.out.println("Initial Value: " + initialValue.toString());
        MathContext mc = new MathContext(precision, RoundingMode.HALF_UP);
        BigDecimal r1 = rk1.add(rk2);
        r1 = r1.multiply(new BigDecimal (0.5), mc);
        BigDecimal r2 = rk3.add(rk4);
        BigDecimal result = r1.add(r2);
        result = result.divide(new BigDecimal(3.0), mc);
        result = initialValue.add(result);
        //System.out.println("Result: " + result.toString());
        return result;
    }

    
}
