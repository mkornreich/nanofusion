/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.visualization.threed;

import java.util.List;
import javax.swing.JPanel;

/**
 *
 * @author pkornrei
 */
public interface Plotter {
    
    public JPanel plot(List<PlotSet> plots);
    
    public JPanel plot(List<PlotSet> plots, PlotBounds bounds);
    
    public void setBounds (PlotBounds bounds);
    
    
}
