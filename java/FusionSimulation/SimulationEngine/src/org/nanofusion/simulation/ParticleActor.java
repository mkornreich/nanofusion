/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.simulation;

import java.math.BigDecimal;

/**
 *
 * @author pkornrei
 */
public interface ParticleActor {
    public ParticleInteraction getInteraction(PVector d1, PVector d2, BigDecimal c1, BigDecimal c2, BigDecimal time, int precision);
}
