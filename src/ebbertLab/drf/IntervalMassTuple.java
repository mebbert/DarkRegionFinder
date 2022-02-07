/**
 * 
 */
package ebbertLab.drf;

import htsjdk.samtools.util.Interval;

/**
 * @author markebbert
 *
 */
public class IntervalMassTuple {
	  public final Interval interval; // The interval object
	  public final double mass;// The mass as a percentage

	  public IntervalMassTuple(Interval interval, double mass) { 
	    this.interval = interval; 
	    this.mass = mass; 
	  }
	  
	  public String toString(){
		  return this.interval.toString() + "->" + this.mass;
	  }
}
