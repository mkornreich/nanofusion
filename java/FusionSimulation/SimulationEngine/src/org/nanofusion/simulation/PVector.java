/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.simulation;

import java.io.Serializable;
import java.math.BigDecimal;

/**
 * 
 *    Bean object representing a three-dimensional vector, such as
 * position or velocity. This vector is not immutable.
 *
 * @author pkornrei
 */
public class PVector implements Serializable{
    
    /** The value of the vector along the x-axis. */
    private BigDecimal x = BigDecimal.ZERO;
    /** The value of the vector along the y-axis. */
    private BigDecimal y = BigDecimal.ZERO;
    /** The value of the vector along the z-axis. */
    private BigDecimal z = BigDecimal.ZERO;
    
    public PVector () {
    }

    public PVector (BigDecimal x, BigDecimal y, BigDecimal z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public BigDecimal getX() {
        return x;
    }

    public void setX(BigDecimal x) {
        this.x = x;
    }

    public BigDecimal getY() {
        return y;
    }

    public void setY(BigDecimal y) {
        this.y = y;
    }

    public BigDecimal getZ() {
        return z;
    }

    public void setZ(BigDecimal z) {
        this.z = z;
    }

    @Override
    public PVector clone() {
        PVector clove = new PVector(x, y, z);
        return clove;
    }
    
    /**
     * 
     *    Adds the given vector to this vector.
     * 
     * @param pv 
     */
    public void add(PVector pv) {
        x = x.add(pv.x);
        y = y.add(pv.y);
        z = z.add(pv.z);
    }
    
    /**
     * 
     *    Performs scalar multiplication of this vector.
     * 
     * @param scalar 
     */
    public void multiply(BigDecimal scalar) {
        x = x.multiply(scalar);
        y = y.multiply(scalar);
        z = z.multiply(scalar);
    }
}
