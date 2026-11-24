/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.simulation;

import java.math.BigDecimal;
import java.math.RoundingMode;

/**
 *
 * @author pkornrei
 */
public class Constants {
    
    public static int SCALE = 10;
    public static final BigDecimal CHARGE = new BigDecimal("1.60217657E-19");
    public static final BigDecimal TIME_INTERVAL = new BigDecimal("2.0E-17");
    //public static final BigDecimal TIME_INTERVAL = new BigDecimal("1.0E-17");
    public static final BigDecimal INITIAL_LATTICE_SPACING = new BigDecimal("2.04E-10");
    public static final BigDecimal ATOMIC_MASS_UNIT = new BigDecimal("1.660538921E-27");
    public static final BigDecimal DEUTERIUM_MASS = new BigDecimal("2.01410178").multiply(ATOMIC_MASS_UNIT);
    public static final BigDecimal INVERSE_DEUTERIUM_MASS = BigDecimal.ONE.divide(DEUTERIUM_MASS, 10, RoundingMode.HALF_UP);
    public static final BigDecimal DEUTERIUM_CHARGE = BigDecimal.ONE;
    public static final BigDecimal LITHIUM_MASS = new BigDecimal("7.01600455").multiply(ATOMIC_MASS_UNIT);
    public static final BigDecimal INVERSE_LITHIUM_MASS = BigDecimal.ONE.divide(LITHIUM_MASS, 10, RoundingMode.HALF_UP);
    public static final BigDecimal LITHIUM_CHARGE = new BigDecimal("3.0");
    public static final BigDecimal K_E = new BigDecimal("8.9875517873681765E+9");
    public static final BigDecimal POTENTIAL = K_E.multiply(CHARGE).multiply(CHARGE);
    public static final BigDecimal C = new BigDecimal("3.0E8");
    public static final BigDecimal INVERSE_C = BigDecimal.ONE.divide(C, 10, RoundingMode.HALF_UP);
    public static final BigDecimal LAMBDA = new BigDecimal("3.6E-6");
    public static final BigDecimal OMEGA = new BigDecimal("2.0").multiply(new BigDecimal(Math.PI)).multiply(C).divide(LAMBDA, 18, RoundingMode.HALF_UP);
    
    
    private Constants() {
    }
    
}
