/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.visualization.threed;

import java.awt.Color;
import java.util.List;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author pkornrei
 */
public class PlotSet {
    
    private List<PlotPoint> points;
    
    private Color color;
    
    private String name;

    public List<PlotPoint> getPoints() {
        return points;
    }

    public void setPoints(List<PlotPoint> points) {
        this.points = points;
    }

    public Color getColor() {
        return color;
    }

    public void setColor(Color color) {
        this.color = color;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }


    
    public double[][] pointsToData() {
        if (points!=null && points.size()>0) {
            double[][] result = new double[points.size()][3];
            for (int p=0; p<points.size(); p++) {
                PlotPoint point=points.get(p);
                result[p][0]=point.getX();
                result[p][1]=point.getY();
                result[p][2]=point.getZ();
            }
            return result;
        }
        return null;
    }
    
    public XYSeries pointsToSliceData() {
         //XYSeriesCollection seriesCollection = new XYSeriesCollection();
         XYSeries series = new XYSeries(name);
         if (points!=null && points.size()>0) {
             for (PlotPoint point: points) {
                 series.add(point.getX(), point.getY());
             }
             return series;
         }
         return null;
    }
    
}
