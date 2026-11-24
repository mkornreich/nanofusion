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
public class ParticleInteraction {
    
    private BigDecimal distance;
    
    private PVector force;
    
    private PVector rkd1;
    private PVector rkd2;
    private PVector rkd3;
    private PVector rkd4;

    private PVector rkv1;
    private PVector rkv2;
    private PVector rkv3;
    private PVector rkv4;

    public BigDecimal getDistance() {
        return distance;
    }

    public void setDistance(BigDecimal distance) {
        this.distance = distance;
    }

    public PVector getForce() {
        return force;
    }

    public void setForce(PVector force) {
        this.force = force;
    }

    public PVector getRkd1() {
        return rkd1;
    }

    public void setRkd1(PVector rkd1) {
        this.rkd1 = rkd1;
    }

    public PVector getRkd2() {
        return rkd2;
    }

    public void setRkd2(PVector rkd2) {
        this.rkd2 = rkd2;
    }

    public PVector getRkd3() {
        return rkd3;
    }

    public void setRkd3(PVector rkd3) {
        this.rkd3 = rkd3;
    }

    public PVector getRkd4() {
        return rkd4;
    }

    public void setRkd4(PVector rkd4) {
        this.rkd4 = rkd4;
    }

    public PVector getRkv1() {
        return rkv1;
    }

    public void setRkv1(PVector rkv1) {
        this.rkv1 = rkv1;
    }

    public PVector getRkv2() {
        return rkv2;
    }

    public void setRkv2(PVector rkv2) {
        this.rkv2 = rkv2;
    }

    public PVector getRkv3() {
        return rkv3;
    }

    public void setRkv3(PVector rkv3) {
        this.rkv3 = rkv3;
    }

    public PVector getRkv4() {
        return rkv4;
    }

    public void setRkv4(PVector rkv4) {
        this.rkv4 = rkv4;
    }


    
    
}
