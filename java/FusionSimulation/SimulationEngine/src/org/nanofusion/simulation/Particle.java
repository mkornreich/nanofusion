/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.simulation;

import com.google.gson.Gson;
import java.io.Serializable;
import java.io.StringWriter;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

/**
 *
 * @author pkornrei
 */
public class Particle implements Serializable {
    
    int id;
    
    private PVector position;
    private PVector velocity = new PVector();
    private String type;
    private transient int precision = 6;
    private transient BigDecimal mass;
    private transient BigDecimal charge;
    private transient PVector currentForce;
    private transient LVector latticePosition;
    private transient BigDecimal inverseMass;
    private transient PVector fieldForce = null;
    private transient PVector rkdx1 = null;
    private transient PVector rkdx2 = null;
    private transient PVector rkdx3 = null;
    private transient PVector rkdx4 = null;
    private transient PVector rkdv1 = null;
    private transient PVector rkdv2 = null;
    private transient PVector rkdv3 = null;
    private transient PVector rkdv4 = null;
    private transient BigDecimal lorentzMultiplier;

    public Particle(PVector position) {
        this(position, null, null);
    }

    public Particle(PVector position, BigDecimal charge) {
        this(position, charge, null);
    }

    public Particle(PVector position, BigDecimal charge, PVector electricField) {
        this.position = position;
        this.charge = charge;
        if (electricField!=null) {
            this.fieldForce = new PVector(electricField.getX().multiply(charge).multiply(Constants.CHARGE),
                    electricField.getY().multiply(charge).multiply(Constants.CHARGE),
                    electricField.getZ().multiply(charge).multiply(Constants.CHARGE));
            this.lorentzMultiplier=this.charge.multiply(Constants.CHARGE).multiply(Constants.INVERSE_C);
        }
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public PVector getPosition() {
        return position;
    }

    public void setPosition(PVector position) {
        this.position = position;
    }

    public PVector getVelocity() {
        return velocity;
    }

    public void setVelocity(PVector velocity) {
        this.velocity = velocity;
    }

    public BigDecimal getMass() {
        return mass;
    }

    public void setMass(BigDecimal mass) {
        this.mass = mass;
    }

    public BigDecimal getInverseMass() {
        return inverseMass;
    }

    public void setInverseMass(BigDecimal inverseMass) {
        this.inverseMass = inverseMass;
    }

    
    public BigDecimal getCharge() {
        return charge;
    }

    public void setCharge(BigDecimal charge) {
        this.charge = charge;
    }

    public PVector getCurrentForce() {
        return currentForce;
    }

    public PVector getFieldForce() {
        return fieldForce;
    }
    
    

    public void setCurrentForce(PVector currentForce, BigDecimal t, int precision, boolean includeLorentz) {
        //System.out.println("Particle: "  + this.getLatticePosition().toString());
        MathContext mc = new MathContext(precision, RoundingMode.HALF_UP);
        this.currentForce = currentForce;
        //System.out.println("Coulomb Force (X): " + currentForce.getX().toString());
        if (this.fieldForce!=null) {
            BigDecimal angle = Constants.OMEGA.multiply(t, mc);
            double wave = Math.sin(angle.doubleValue());
            PVector currentFieldForce = fieldForce.clone();
            currentFieldForce.multiply(new BigDecimal(wave));
            //System.out.println("Field Force (X): " + currentFieldForce.getX().toString());
            this.currentForce.add(currentFieldForce);
            //System.out.println("Current Force (X): " + currentForce.getX().toString());
       }
       if (includeLorentz) {
           BigDecimal lx = this.lorentzMultiplier.multiply(this.velocity.getZ());
           BigDecimal lz = BigDecimal.ZERO.subtract(this.lorentzMultiplier.multiply(this.velocity.getX()));
           PVector lorentzForce = new PVector(lx, BigDecimal.ZERO, lz);
           this.currentForce.add(lorentzForce);
       }
    }

    public LVector getLatticePosition() {
        return latticePosition;
    }

    public void setLatticePosition(LVector latticePosition) {
        this.latticePosition = latticePosition;
    }

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }

    public int getPrecision() {
        return precision;
    }

    public void setPrecision(int precision) {
        this.precision = precision;
    }
    
    
    public Particle clone() {
        Particle clopa = new Particle(position.clone());
        clopa.setVelocity(velocity.clone());
        clopa.setMass(mass);
        clopa.setInverseMass(inverseMass);
        clopa.setCharge(charge);
        clopa.setId(id);
        clopa.setType(type);
        clopa.fieldForce = fieldForce;
        clopa.setLatticePosition(latticePosition);
        clopa.setPrecision(precision);
        clopa.lorentzMultiplier = lorentzMultiplier;
        return clopa;
    }
    
    public void update(BigDecimal timeStep) {
        //if (this.getLatticePosition().getY()==3 && this.getLatticePosition().getZ()==3) {
        //    System.out.println("Particle: "  + this.getLatticePosition().toString());
        //}
        position.setX(CommonCalculatorUtility.applyElapsedTimeChange(position.getX(), velocity.getX(), precision, timeStep));
        position.setY(CommonCalculatorUtility.applyElapsedTimeChange(position.getY(), velocity.getY(), precision, timeStep));
        position.setZ(CommonCalculatorUtility.applyElapsedTimeChange(position.getZ(), velocity.getZ(), precision, timeStep));
        //if (this.getLatticePosition().getY()==3 && this.getLatticePosition().getZ()==3) {
        //    System.out.println("Position (X): " + position.getX().toString());
        //}
        
        currentForce.multiply(inverseMass);
        /*if (this.getLatticePosition().getY()==3 && this.getLatticePosition().getZ()==3) {
            System.out.println("Acceleration (X): " + currentForce.getX().toString());
        }*/
        velocity.setX(CommonCalculatorUtility.applyElapsedTimeChange(velocity.getX(), currentForce.getX(), precision, timeStep));
        velocity.setY(CommonCalculatorUtility.applyElapsedTimeChange(velocity.getY(), currentForce.getY(), precision, timeStep));
        velocity.setZ(CommonCalculatorUtility.applyElapsedTimeChange(velocity.getZ(), currentForce.getZ(), precision, timeStep));
        /*if (this.getLatticePosition().getY()==3 && this.getLatticePosition().getZ()==3) {
            System.out.println("Velocity (X): " + velocity.getX().toString());
        }*/
    }
    
    public void updateRungeKutta() {
        //System.out.println("Particle No. " + latticePosition.getX());
        position.setX(CommonCalculatorUtility.updateRungeKutta(position.getX(), rkdx1.getX(), rkdx2.getX(), rkdx3.getX(), rkdx4.getX(), precision));
        position.setY(CommonCalculatorUtility.updateRungeKutta(position.getY(), rkdx1.getY(), rkdx2.getY(), rkdx3.getY(), rkdx4.getY(), precision));
        position.setZ(CommonCalculatorUtility.updateRungeKutta(position.getZ(), rkdx1.getZ(), rkdx2.getZ(), rkdx3.getZ(), rkdx4.getZ(), precision));
        //if (id == 0) {
        //    System.out.println("Position (X): " + position.getX().toString());
        //}
        
        //System.out.println("Acceleration (X): " + currentForce.getX().toString());
        velocity.setX(CommonCalculatorUtility.updateRungeKutta(velocity.getX(), rkdv1.getX(), rkdv2.getX(), rkdv3.getX(), rkdv4.getX(), precision));
        velocity.setY(CommonCalculatorUtility.updateRungeKutta(velocity.getY(), rkdv1.getY(), rkdv2.getY(), rkdv3.getY(), rkdv4.getY(), precision));
        velocity.setZ(CommonCalculatorUtility.updateRungeKutta(velocity.getZ(), rkdv1.getZ(), rkdv2.getZ(), rkdv3.getZ(), rkdv4.getZ(), precision));
        //System.out.println("Velocity (X): " + velocity.getX().toString());
    }
    
    public void setRungeKuttaParameters(PVector rkd1, PVector rkd2, PVector rkd3, PVector rkd4,  PVector rkv1, PVector rkv2, PVector rkv3, PVector rkv4) {
        this.rkdx1 = rkd1;
        this.rkdx2 = rkd2;
        this.rkdx3 = rkd3;
        this.rkdx4 = rkd4;
        this.rkdv1 = rkv1;
        this.rkdv2 = rkv2;
        this.rkdv3 = rkv3;
        this.rkdv4 = rkv4;
    }
    
    @Override
    public boolean equals (Object o) {
        if (o!=null && o instanceof Particle) {
            Particle p = (Particle)o;
            return this.id == p.id;
        }
        return false;
    }
    
    @Override
    public String toString () {
        Gson gson = new Gson();
        StringWriter writer = new StringWriter();
        gson.toJson(this, writer);
        return writer.toString();
    }
    
}
