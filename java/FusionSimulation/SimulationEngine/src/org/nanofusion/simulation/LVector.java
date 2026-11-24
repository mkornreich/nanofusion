//BS"D
package org.nanofusion.simulation;

/**
 *
 * @author pkornrei
 */
public class LVector {
    
    private int x;
    
    private int y;
    
    private int z;
    
    public LVector(int x, int y, int z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public int getX() {
        return x;
    }

    public void setX(int x) {
        this.x = x;
    }

    public int getY() {
        return y;
    }

    public void setY(int y) {
        this.y = y;
    }

    public int getZ() {
        return z;
    }

    public void setZ(int z) {
        this.z = z;
    }
    
    @Override
    public String toString() {
        String s="X: " + this.x + ", Y: " + this.y + ", Z: " + this.z;
        return s;
    }
    
    /**
     * 
     * @param lvec
     * @return 
     */
    public boolean isNearNeighbor(LVector lvec) {
        int ldist = Math.abs(x - lvec.x) + Math.abs(y - lvec.y) + Math.abs(z - lvec.z);
        return ldist > 0 && ldist < 3;
    }
    
    @Override
    public boolean equals(Object o) {
        if (o!=null && o instanceof LVector) {
            LVector lv = (LVector)o;
            return ((this.x == lv.x) && (this.y == lv.y) && (this.z == lv.z));
        }
        return false;
    }
}
