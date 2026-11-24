/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.simulation;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

/**
 *
 * @author pkornrei
 */
public class FieldActor implements ParticleActor {
    
    PVector fieldForce;
    
    public FieldActor(BigDecimal charge, PVector electricField) {
        if (electricField!=null) {
            this.fieldForce = new PVector(electricField.getX().multiply(charge).multiply(Constants.CHARGE),
                    electricField.getY().multiply(charge).multiply(Constants.CHARGE),
                    electricField.getZ().multiply(charge).multiply(Constants.CHARGE));
        }
    }

    @Override
    public ParticleInteraction getInteraction(PVector d1, PVector d2, BigDecimal c1, BigDecimal c2, BigDecimal time, int precision) {
        ParticleInteraction result = new ParticleInteraction();
        MathContext mc = new MathContext(precision, RoundingMode.HALF_UP);
        //fieldForce = 
        BigDecimal angle = Constants.OMEGA.multiply(time, mc);
        double wave = Math.sin(angle.doubleValue());
        PVector currentFieldForce = fieldForce.clone();
        currentFieldForce.multiply(new BigDecimal(wave));
        //System.out.println("Field Force (X): " + currentFieldForce.getX().toString());
        result.setForce(currentFieldForce);
        return result;
    }

    
}
