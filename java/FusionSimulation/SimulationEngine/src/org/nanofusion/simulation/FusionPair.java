/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.simulation;

import com.google.gson.Gson;
import java.io.Serializable;
import java.io.StringWriter;

/**
 *
 * @author pkornrei
 */
public class FusionPair implements Serializable {
    
    Particle particle1 = null;
    Particle particle2 = null;

    public Particle getParticle1() {
        return particle1;
    }

    public void setParticle1(Particle particle1) {
        this.particle1 = particle1;
    }

    public Particle getParticle2() {
        return particle2;
    }

    public void setParticle2(Particle particle2) {
        this.particle2 = particle2;
    }
    
    @Override
    public boolean equals(Object o) {
        if (o!=null && o instanceof FusionPair) {
            FusionPair fp = (FusionPair)o;
            return (
                    (this.particle1.equals(fp.particle1) && this.particle2.equals(fp.particle2)) || 
                    (this.particle1.equals(fp.particle2) && this.particle2.equals(fp.particle1)));
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
