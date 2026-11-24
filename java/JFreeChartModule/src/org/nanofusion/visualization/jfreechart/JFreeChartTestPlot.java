/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.visualization.jfreechart;

import java.awt.Color;
import java.awt.Shape;
import javax.swing.JPanel;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.util.ShapeUtilities;

/**
 *
 * @author pkornrei
 */
public class JFreeChartTestPlot {
    
    public JPanel getTestPlot(XYDataset data) {
        JFreeChart chart = ChartFactory.createScatterPlot("Test Plot", "X", "", data);
        NumberAxis dAxis = new NumberAxis();
        dAxis.setLowerBound(-1.5);
        dAxis.setUpperBound(1.5);
        NumberAxis rAxis = new NumberAxis();
        rAxis.setLowerBound(-1.5);
        XYPlot plot = chart.getXYPlot();
        rAxis.setUpperBound(1.5);
        plot.setDomainAxis(dAxis);
        plot.setRangeAxis(rAxis);
        Shape shape = ShapeUtilities.createDiamond(3);
        XYItemRenderer renderer = plot.getRenderer();
        renderer.setSeriesShape(0, shape);
        renderer.setSeriesPaint(0, Color.BLUE);
        ChartPanel panel = new ChartPanel(chart);
        return panel;
    }
}
