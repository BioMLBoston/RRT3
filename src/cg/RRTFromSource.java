package cg;

import io.DCDIO;
import io.Writer;

import java.io.DataOutputStream;
import java.io.File;
import java.util.ArrayList;
import java.util.Calendar;

/**
 * @author Dong Luo 2012-08-21
 *
 */

public class RRTFromSource implements Runnable {
	private RRT rrt;
	private String base;
	
	public RRTFromSource(RRT rrt, String base) {
		this.rrt = rrt;
		this.base = base+"/";
	}

	@Override
	public void run() {
		String time = "\nRRT Start time: "+Calendar.getInstance().getTime()+"\n";
		String energy = "";
		ArrayList<CGMolecule> path = rrt.searchPathFromSource();
		time += "RRT End time: "+Calendar.getInstance().getTime()+"\n";
		System.out.println(time);
		if (path != null) {
			System.out.printf("\nA path found is saved to disk.\n");
			DataOutputStream out = DCDIO.getOutHandle(base+"path.dcd");
			DCDIO.writeHeader(out, path.get(0), "path.dcd");
			int i = 0;
			for (CGMolecule m:path) {
				DCDIO.writeData(out, m);
				i++;
				energy += String.format("%d %.2f %.2f %.2f %.2f %.2f\n", i, m.getEnergy(), m.getEnergy(GordonEnergy.angleType), m.getEnergy(GordonEnergy.dihedralType), m.getEnergy(GordonEnergy.vanDeWallsType), m.getEnergy(GordonEnergy.hydrogenBondType));
			}
			DCDIO.closeOutHandle(out);
			DCDIO.writeFrameNumber(base+"path.dcd", i);
			Writer.writeStringToFile(time+rrt.getLog(), base+"log.txt");
			Writer.writeStringToFile(energy, base+"energy.txt");
		} else System.out.println("No path is found.");
	                  }

	public static void main(String[] args) {
            System.out.println("Working Directory = " +
              System.getProperty("user.dir"));
		if (args.length < 5) {
			System.out.println("Usage: java -jar testRRT work-dir start-pdb-file goal-pdb-file output-dir score-threshold");
			System.out.println("\t\tStart-pdb-file and goal-pdb-file are assumed be at work-dir");
			System.out.println("\t\tOutput-dir is relative to work-dir");
			System.exit(0);
		                      }
		CGMolecule m1 = (CGMolecule) CGMoleculeComplex.newCGMoleculeComplex(args[0]+"/"+args[1]).getMolecule("A");
		CGMolecule m2 = (CGMolecule) CGMoleculeComplex.newCGMoleculeComplex(args[0]+"/"+args[2]).getMolecule("A");
		String base = args[0]+"/"+args[3];
		new File(base).mkdirs();
		RRT rrt = new RRT(m1, m2, base, Integer.parseInt(args[4]));
               
		Thread ti = new Thread(new RRTFromSource(rrt, base));
               
		ti.start();
		try {
			ti.join();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
                                                                
	}
}
