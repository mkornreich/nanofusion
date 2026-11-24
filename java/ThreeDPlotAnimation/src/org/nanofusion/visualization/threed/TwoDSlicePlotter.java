/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.visualization.threed;

import java.awt.Color;
import java.awt.Shape;
import java.util.List;
import javax.swing.JPanel;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.util.ShapeUtilities;

/**
 *
 * @author pkornrei
 */
public class TwoDSlicePlotter implements Plotter {
    
    private PlotBounds bounds;

    @Override
    public JPanel plot(List<PlotSet> plots) {
        XYSeriesCollection seriesCollection = new XYSeriesCollection();
        if (bounds!=null) {
            double[] minima = bounds.getMinima();
            double[] maxima = bounds.getMinima();
            NumberAxis dAxis = new NumberAxis();
            dAxis.setLowerBound(minima[0]);
            dAxis.setUpperBound(maxima[0]);
            NumberAxis rAxis = new NumberAxis();
            rAxis.setLowerBound(minima[0]);
            rAxis.setUpperBound(maxima[0]);
        }
        for (PlotSet points: plots) {
            seriesCollection.addSeries(points.pointsToSliceData());
        }
        JFreeChart chart = ChartFactory.createScatterPlot("Plot", "X", "Y", seriesCollection);
        XYPlot plot = chart.getXYPlot();
        int p = 0;
        for (PlotSet points: plots) {
            Shape shape = ShapeUtilities.createDiamond(2);
            XYItemRenderer renderer = plot.getRenderer();
            renderer.setSeriesShape(p, shape);
            renderer.setSeriesPaint(p, points.getColor());
            p++;
        }
        ChartPanel panel = new ChartPanel(chart);
        return panel;
    }

    @Override
    public JPanel plot(List<PlotSet> plots, PlotBounds bounds) {
        setBounds(bounds);
        return plot(plots);
    }

    @Override
    public void setBounds(PlotBounds bounds) {
        this.bounds = bounds;
    }
    
}
