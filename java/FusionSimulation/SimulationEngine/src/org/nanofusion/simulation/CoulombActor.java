/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.simulation;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import static org.nanofusion.simulation.CommonCalculatorUtility.getDirection;
import static org.nanofusion.simulation.CommonCalculatorUtility.getSquareDistance;
import static org.nanofusion.simulation.CommonCalculatorUtility.sqrt;

/**
 *
 * @author pkornrei
 */
public class CoulombActor implements ParticleActor {

   public ParticleInteraction getInteraction(PVector d1, PVector d2, BigDecimal c1, BigDecimal c2, BigDecimal time, int precision) {
        MathContext mc = new MathContext(precision, RoundingMode.HALF_UP);
        BigDecimal potential = Constants.POTENTIAL.multiply(c1, mc).multiply(c2, mc);
        BigDecimal sqdist = getSquareDistance(d1, d2, precision);
        BigDecimal distance = sqrt(sqdist);
        BigDecimal unitMagnitude = potential.divide(distance.multiply(sqdist, mc), mc);
        PVector force = getDirection(d1, d2);
        force.multiply(unitMagnitude);
        ParticleInteraction xaction = new ParticleInteraction();
        xaction.setDistance(distance);
        xaction.setForce(force);
        return xaction;
    }
}
