/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.visualization.threed;

import java.awt.Color;
import java.util.List;
import javax.swing.JPanel;
import org.math.plot.canvas.Plot3DCanvas;

/**
 *
 * @author pkornrei
 */
public class ThreeDPlotter implements Plotter {
    
    private PlotBounds bounds;
    
    private Plot3DCanvas generatePlot(List<PlotSet> plots) {
        Plot3DCanvas canvas = new Plot3DCanvas();
        if (plots!=null && plots.size()>0) {
            for (PlotSet plot: plots) {
                canvas.addScatterPlot("Test", plot.getColor(), plot.pointsToData());
            }
            
        }
        return canvas;
    }

    @Override
    public JPanel plot(List<PlotSet> plots) {
        Plot3DCanvas canvas = generatePlot(plots);
        if (this.bounds!=null) {
            canvas.setFixedBounds(bounds.getMinima(), bounds.getMaxima());
        }
        return canvas;
    }

    @Override
    public JPanel plot(List<PlotSet> plots, PlotBounds bounds) {
        Plot3DCanvas canvas = generatePlot(plots);
        if (bounds!=null) {
            canvas.setFixedBounds(bounds.getMinima(), bounds.getMaxima());
        }
        return canvas;
    }

    @Override
    public void setBounds(PlotBounds bounds) {
        this.bounds = bounds;
    }

    
}
