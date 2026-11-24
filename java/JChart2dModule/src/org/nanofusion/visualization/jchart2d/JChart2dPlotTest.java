/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.visualization.jchart2d;

import info.monitorenter.gui.chart.Chart2D;
import info.monitorenter.gui.chart.ITrace2D;
import info.monitorenter.gui.chart.TracePoint2D;
import info.monitorenter.gui.chart.pointpainters.PointPainterDisc;
import info.monitorenter.gui.chart.traces.Trace2DLtd;
import info.monitorenter.gui.chart.traces.Trace2DSimple;
import info.monitorenter.gui.chart.traces.painters.TracePainterDisc;
import info.monitorenter.gui.chart.views.ChartPanel;
import java.awt.Color;
import java.util.List;
import javax.swing.JLayeredPane;

/**
 *
 * @author pkornrei
 */
public class JChart2dPlotTest {
    
    public JLayeredPane getTestPlot(List<TracePoint2D> points) {
        Chart2D chart = new Chart2D();
        ITrace2D trace = new Trace2DSimple();
        chart.addTrace(trace);
        trace.setTracePainter(new TracePainterDisc(4));
        trace.setColor(Color.DARK_GRAY);
        for (TracePoint2D point: points) {
            trace.addPoint(point);
        }
        JLayeredPane pane = new ChartPanel(chart);
        return pane;
    }
    
}
