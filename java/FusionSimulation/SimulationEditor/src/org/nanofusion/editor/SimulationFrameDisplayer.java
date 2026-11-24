/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nanofusion.editor;

import java.awt.Color;
import java.awt.Dimension;
import java.util.ArrayList;
import java.util.List;
import javax.swing.JLabel;
import javax.swing.JPanel;
import org.nanofusion.simulation.Constants;
import org.nanofusion.simulation.Lattice;
import org.nanofusion.simulation.Particle;
import org.nanofusion.visualization.threed.PlotBounds;
import org.nanofusion.visualization.threed.PlotPoint;
import org.nanofusion.visualization.threed.PlotSet;
import org.nanofusion.visualization.threed.ThreeDPlotter;
import org.nanofusion.visualization.threed.TwoDSlicePlotter;

/**
 *
 * @author pkornrei
 */
public class SimulationFrameDisplayer {
    
    public static void displayFrame(List<Lattice> simulation, JPanel panel, int frameNo, JLabel lblFrameNo) {
        Lattice lattice = simulation.get(frameNo);
        List<Particle> particles = lattice.getParticles();
        PlotSet deut = new PlotSet();
        PlotSet lith = new PlotSet();
        List<PlotPoint> deutPoints = new ArrayList<PlotPoint>();
        List<PlotPoint> lithPoints = new ArrayList<PlotPoint>();
        double lowerLimit = -1.* Constants.INITIAL_LATTICE_SPACING.doubleValue();
        double middleLimit = 12.* Constants.INITIAL_LATTICE_SPACING.doubleValue();
        double upperLimit = 25.* Constants.INITIAL_LATTICE_SPACING.doubleValue();
        double[] lowerLimits = new double[] {lowerLimit, lowerLimit, lowerLimit};
        double[] upperLimits = new double[] {upperLimit, middleLimit, middleLimit};
        PlotBounds bounds = new PlotBounds(lowerLimits, upperLimits);
        for (Particle particle: particles) {
            PlotPoint point = new PlotPoint();
            point.setX(particle.getPosition().getX().doubleValue());
            point.setY(particle.getPosition().getY().doubleValue());
            point.setZ(particle.getPosition().getZ().doubleValue());
            if (particle.getCharge().doubleValue()>1.5) {
                lithPoints.add(point);
            }
            else {
                deutPoints.add(point);
            }
        }
        deut.setPoints(deutPoints);
        lith.setPoints(lithPoints);
        deut.setColor(Color.BLUE);
        lith.setColor(Color.RED);
        List<PlotSet> sets = new ArrayList<PlotSet>();
        sets.add(deut);
        sets.add(lith);
        ThreeDPlotter plotter = new ThreeDPlotter();
        javax.swing.JPanel plot = plotter.plot(sets, bounds);
        /*for (int ff=0; ff<simulation.size(); ff++) {
            Particle particle = simulation.get(ff).getParticles().get(particleNo);
            PlotPoint point = new PlotPoint();
            PVector position = particle.getPosition();
            point.setX(position.getX().doubleValue());
            point.setY(position.getY().doubleValue());
            point.setZ(position.getY().doubleValue());
            deutPoints.add(point);
        }
        deut.setPoints(deutPoints);
        deut.setColor(Color.BLUE);
        List<PlotSet> sets = new ArrayList<PlotSet>();
        sets.add(deut);
        ThreeDPlotter plotter = new ThreeDPlotter();
        javax.swing.JPanel plot = plotter.plot(sets);*/
        Dimension parentSize = panel.getSize();
        //Dimension pnlSize = new Dimension();
        plot.setSize(parentSize);
        panel.removeAll();
        panel.repaint();
        panel.add(plot);
        panel.repaint();
        int displayFrameNo = frameNo + 1;
        lblFrameNo.setText(displayFrameNo + "");
    }
    
    public static void displayFrame2D(List<Lattice> simulation, JPanel panel, int frameNo, JLabel lblFrameNo, int z) {
        Lattice lattice = simulation.get(frameNo);
        List<Particle> particles = lattice.getParticles();
        PlotSet deut = new PlotSet();
        PlotSet lith = new PlotSet();
        List<PlotPoint> deutPoints = new ArrayList<PlotPoint>();
        List<PlotPoint> lithPoints = new ArrayList<PlotPoint>();
        double lowerLimit = -5.* Constants.INITIAL_LATTICE_SPACING.doubleValue();
        double upperLimit = 20.* Constants.INITIAL_LATTICE_SPACING.doubleValue();
        double[] lowerLimits = new double[] {lowerLimit, lowerLimit, lowerLimit};
        double[] upperLimits = new double[] {upperLimit, upperLimit, upperLimit};
        PlotBounds bounds = new PlotBounds(lowerLimits, upperLimits);
        for (Particle particle: particles) {
            if (particle.getLatticePosition().getZ() == z) {
                PlotPoint point = new PlotPoint();
                point.setX(particle.getPosition().getX().doubleValue());
                point.setY(particle.getPosition().getY().doubleValue());
                if (particle.getCharge().doubleValue()>1.5) {
                    lithPoints.add(point);
                }
                else {
                    deutPoints.add(point);
                }
            }
        }
        deut.setPoints(deutPoints);
        lith.setPoints(lithPoints);
        deut.setColor(Color.BLUE);
        deut.setName("Deuterium");
        lith.setColor(Color.RED);
        lith.setName("Lithium");
        List<PlotSet> sets = new ArrayList<PlotSet>();
        sets.add(deut);
        sets.add(lith);
        TwoDSlicePlotter plotter = new TwoDSlicePlotter();
        javax.swing.JPanel plot = plotter.plot(sets, bounds);
        /*for (int ff=0; ff<simulation.size(); ff++) {
            Particle particle = simulation.get(ff).getParticles().get(particleNo);
            PlotPoint point = new PlotPoint();
            PVector position = particle.getPosition();
            point.setX(position.getX().doubleValue());
            point.setY(position.getY().doubleValue());
            point.setZ(position.getY().doubleValue());
            deutPoints.add(point);
        }
        deut.setPoints(deutPoints);
        deut.setColor(Color.BLUE);
        List<PlotSet> sets = new ArrayList<PlotSet>();
        sets.add(deut);
        ThreeDPlotter plotter = new ThreeDPlotter();
        javax.swing.JPanel plot = plotter.plot(sets);*/
        Dimension parentSize = panel.getSize();
        //Dimension pnlSize = new Dimension();
        plot.setSize(parentSize);
        panel.removeAll();
        panel.repaint();
        panel.add(plot);
        panel.repaint();
        int displayFrameNo = frameNo + 1;
        lblFrameNo.setText(displayFrameNo + "");
    }
    
}
