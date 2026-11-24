/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.visualization.threed;

/**
 *
 * @author pkornrei
 */
public class PlotBounds {
    
    private double[] minima;
    private double[] maxima;
    
    public PlotBounds() {
        
    }
    
    public PlotBounds(double[] minima, double[] maxima) {
        this.maxima = maxima;
        this.minima = minima;
    }

    public double[] getMinima() {
        return minima;
    }

    public void setMinima(double[] minima) {
        this.minima = minima;
    }

    public double[] getMaxima() {
        return maxima;
    }

    public void setMaxima(double[] maxima) {
        this.maxima = maxima;
    }
    
    
}
