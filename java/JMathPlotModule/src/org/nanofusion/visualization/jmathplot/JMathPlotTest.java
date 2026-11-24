/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.visualization.jmathplot;

import java.awt.Color;
import javax.swing.JPanel;
import org.math.plot.canvas.Plot3DCanvas;

/**
 *
 * @author pkornrei
 */
public class JMathPlotTest {
    
    public JPanel getTestPlot(double[] x, double[] y, double[] z) {
        Plot3DCanvas canvas = new Plot3DCanvas();
        canvas.addScatterPlot("Test", Color.red, x, y, z);
        return canvas;
    }
}
