/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.editor;

import com.google.gson.Gson;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import org.nanofusion.simulation.Lattice;
import org.netbeans.spi.actions.AbstractSavable;
import org.openide.filesystems.FileChooserBuilder;

/**
 *
 * @author pkornrei
 */
public class SimulationSavable extends AbstractSavable implements SaveFinishListener {
    
    private Simulation simulation;
    private SimulationSavedListener listener;
    private SaveProgressDialog progressDialog = null;
    
    SimulationSavable(Simulation simulation,SimulationSavedListener listener) {
        this.simulation = simulation;
        this.listener = listener;
        register();
    }
    
    @Override
    protected String findDisplayName() {
        return simulation.getProperties().getName();
    }

    @Override
    protected void handleSave() throws IOException {
        if (progressDialog==null) {
            progressDialog = new SaveProgressDialog();
        }
        progressDialog.saveSimulation(simulation, this);
        progressDialog.setVisible(true);
    }

    @Override
    public boolean equals(Object o) {
        if (o instanceof SimulationSavable) {
            SimulationSavable s = (SimulationSavable)o;
            if (this.simulation.equals(o)) {
                return true;
            }
        }
        return false;
    }

    @Override
    public int hashCode() {
        return simulation.hashCode();
    }

    @Override
    public void notifySaveFinished() {
        listener.onSimulationSaved(this);
        unregister();
    }
    
}
