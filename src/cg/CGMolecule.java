package cg;

/**
 * @author Dong Luo 2012-05-23
 *
 */

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.StringCharacterIterator;
import java.util.ArrayList;
import java.util.ListIterator;
import java.util.Random;

import cg.GordonEnergy.DihedralTypes;
import cg.GordonEnergy.HydrogenBondCapability;

import pdb.PDBAtom;
import pdb.PDBMolecule;
import utilities.Vector;

public class CGMolecule extends PDBMolecule {
	private ArrayList<DihedralTypes> dihedralTypeList=null;
	private ArrayList<Element> elementList=null;	// rigid elements
	private ArrayList<Integer> freeDihedralAngleList=null;	// free dihedral angle indexes, depend on elements
	private ArrayList<Integer> freeAngleList=null;	// free 3-atom angle indexes
	private ArrayList<Double> angleVector=null, distVector=null;	// debug use
	private ArrayList<Double> flexIndex=null;
	

	private CGMolecule(ArrayList<PDBAtom> list, String chainID) {
		super(list, chainID);
	}

	/**
	 * A deep copy is implemented
	 * @param molecule is the one to be cloned
	 * @return a new coarse grained molecule contains the C Alpha atoms only from given molecule
	 */
	public static CGMolecule newCGMolecule(PDBMolecule molecule) {
		ArrayList<PDBAtom> list = new ArrayList<PDBAtom>();
		ListIterator<PDBAtom> iter = molecule.extractBackbone().listIterator();
		while (iter.hasNext()) {
			PDBAtom a = iter.next();
			if (a != null) {
                            
                                String xyzStr = a.getCoordinatesString();
				double x = Double.parseDouble(xyzStr.substring(0, 7).trim());
				double y = Double.parseDouble(xyzStr.substring(8, 15).trim());
				double z = Double.parseDouble(xyzStr.substring(16, 23).trim());
				list.add(new PDBAtom(a.getIndex(), x, y, z, a.getAtomType(), a.getAminoAcidType(), a.getChainID(), a.getResidueIndex()));
			}
		}
		CGMolecule inst = new CGMolecule(list, molecule.getChainID());
		if (molecule instanceof CGMolecule) {
			CGMolecule old = (CGMolecule)molecule;
			if (old.dihedralTypeList != null)
				inst.dihedralTypeList = new ArrayList<DihedralTypes>(old.dihedralTypeList);
			if (old.freeAngleList != null)
				inst.freeAngleList = new ArrayList<Integer>(old.freeAngleList);
			if (old.freeDihedralAngleList != null)
				inst.freeDihedralAngleList = new ArrayList<Integer>(old.freeDihedralAngleList);
			if (old.getElementList() != null) {
				inst.elementList = new ArrayList<Element>();
				for (Element e:old.getElementList())
					inst.elementList.add(new Element(inst, e));
			}
			if (old.flexIndex != null)
				inst.flexIndex = new ArrayList<Double>(old.flexIndex);
		} else {
			inst.initDihedralTypes();
			inst.initElements();
		}
		return inst;
	}
	
	/* 
	 * any continuous "E" or "H" are taken as an element
	 * assume dihedralTypeList is the same size of pdbAtomlist
	 */
	private void initElements() {		
		elementList=new ArrayList<Element>();
		DihedralTypes t0 = dihedralTypeList.get(0), t;
		int offset = 0;
		for (int i=1;i<dihedralTypeList.size();i++) {
			t = dihedralTypeList.get(i);
			if ((t0 != DihedralTypes.H) && (t0 != DihedralTypes.E)) {
				t0 = t; offset = i; continue;
			}
			if ((t != t0) && (i == (offset+1))) {
				t0 = t; offset = i; continue;
			}
			if ((t == t0) && (i != dihedralTypeList.size()-1)) continue;
			if (t != t0) {
				if ((i-offset)>=Element.minimumElementSize) elementList.add(new Element(this, t0, offset, i-1));
				t0 = t; offset = i;
			} else if ((i-offset+1)>=Element.minimumElementSize)
				elementList.add(new Element(this, t0, offset, i));
		}
		// free angle/dihedral angle lists and score vector need to be reset
		freeAngleList = null; freeDihedralAngleList = null;
	}
	
	public ArrayList<Element> getElementList() {
		if (elementList == null) initElements();
		return elementList;
	}
	
	public void setElements(ArrayList<Element> list) {
		elementList = list;
		freeAngleList = null; freeDihedralAngleList = null;
	}

	/**
	 * Copy the element assignment from another molecule
	 * @param mol is the molecule to be copied
	 */
	public void setElements(CGMolecule mol) {
		if (elementList != null)
			elementList.clear();
		else
			elementList = new ArrayList<Element>();
		for (Element e:mol.getElementList())
			elementList.add(new Element(this, e));
		freeAngleList = null; freeDihedralAngleList = null;
	}

	/*
	 * This function offer the opportunity to manually assign rigid element
	 * Remove any element that is overlapping with new element
	 * add new element and keep in sorted order
	 */
	public void setElement(DihedralTypes t, int start, int end) {
		if (start < 0 || end > dihedralTypeList.size()) return;
		if ((end-start+1) < Element.minimumElementSize) return;
		
		getElementList(); // initialize it if not yet
		// remove overlapping elements
		// until find the point to insert the new element
		int i = 0;
		while (i < elementList.size()) {
			Element e = elementList.get(i);
			if (end <= e.getStart()) break;
			if (start < e.getEnd()) elementList.remove(i);
			else i++;
		}
		elementList.add(i, new Element(this, t, start, end));
		
		freeAngleList = null; freeDihedralAngleList = null;
	}
	
	public void removeElement(int offset) {
		if (offset < 0 || offset >= dihedralTypeList.size()) return;
		getElementList(); // initialize it if not yet
		for (int i=0;i<elementList.size();i++)
			if (elementList.get(i).getStart() == offset) {
				elementList.remove(i);
				break;
			}
		freeAngleList = null; freeDihedralAngleList = null;
	}
	
	public void setElement(String filePath) {
		try{
			BufferedReader reader = new BufferedReader(new FileReader(filePath));
			//elementList = new ArrayList<Element>();
			String pdbLine; int gid = 0, offset = 0, cid = -1;
			while( (pdbLine = reader.readLine()) != null){
				if (pdbLine.length() < 64) continue;
				if (pdbLine.substring(13,15).equals("CA")) {
					cid++;
					int id = Integer.parseInt(pdbLine.substring(60,63).trim());
					if (gid != id) {
						if ((cid-offset)>=Element.minimumElementSize)
							setElement(DihedralTypes.H, offset, cid-1);
							//elementList.add(new Element(this, DihedralTypes.H, offset, cid-1));
						gid = id; offset = cid;
					}
				}
			}
		}
		catch(FileNotFoundException e){
			System.out.println("Cannot find file " + filePath);
		}
		catch(IOException e){
			System.out.println("Error occured while reading from file " + filePath);
		}
	}
	
	public void setFlexIndex(String filePath) {
		try{
			BufferedReader reader = new BufferedReader(new FileReader(filePath));
			flexIndex = new ArrayList<Double>();
			String pdbLine;
			while( (pdbLine = reader.readLine()) != null){
				if (pdbLine.length() < 66) continue;
				
                                if (pdbLine.substring(13,15).equals("CA")||pdbLine.substring(13,14).equals("N")||pdbLine.substring(13,14).equals("O")||pdbLine.substring(13,14).equals("C") )
					flexIndex.add(Double.parseDouble(pdbLine.substring(60,66).trim()));
			}
		}
		catch(FileNotFoundException e){
			System.out.println("Cannot find file " + filePath);
		}
		catch(IOException e){
			System.out.println("Error occured while reading from file " + filePath);
		}
		if (flexIndex!=null && flexIndex.size()!=getPDBAtomList().size()) {
			flexIndex = null;
			System.out.println("Setting flexIndex failed, size not agree.");
		}
	}
	
	public ArrayList<Double> getFlexIndex() {
		return flexIndex;
	}
	
	public int getRandomLoopIndex(Random g) {
		ArrayList<Integer> list = getFreeDihedralAngleList();
		if (list == null || list.size() == 0)
			return 0;
		else
			return list.get(g.nextInt(list.size()));
	}
	
	public int getRandomAngleIndex(Random g) {
		if (freeAngleList == null || freeAngleList.size() == 0)
			return 0;
		else
			return freeAngleList.get(g.nextInt(freeAngleList.size()));
	}

	public ArrayList<Integer> getFreeAngleList() {
		if (freeAngleList == null && getElementList() != null) {
			freeAngleList = new ArrayList<Integer>();
			int i = getPDBAtomList().size();
			for (Element e:elementList) {
				for (;i<=e.getStart();i++)
					freeAngleList.add(i);
				i = e.getEnd();
			}
		}
		return freeAngleList;
	}

	public void setFreeAngleList(ArrayList<Integer> list) {
		freeAngleList = list;
	}

	public ArrayList<Integer> getFreeDihedralAngleList() {
		if (freeDihedralAngleList == null && getElementList() != null) {
			freeDihedralAngleList = new ArrayList<Integer>();
			int i = getPDBAtomList().size();
			for (Element e:elementList) {
				for (;i<=e.getStart();i++)
					freeDihedralAngleList.add(i);
				i = e.getEnd();
			}
		}
		return freeDihedralAngleList;
	}
	
	public void setFreeDihedralAngleList(ArrayList<Integer> list) {
		freeDihedralAngleList = list;
	}

	/**
	 * assume m1 and m2 are same molecule but in different configurations
	 * This function will set the same rigid elements for m1 and m2
	 * and their freeAngleList, freeDihedralAngleList as those differed larger than threshold value
	 * @param m1
	 * @param m2
	 * @return the summation of delta of angles
	 */
	public static ArrayList<Double> setElementsAndFreeAngles(CGMolecule m1, CGMolecule m2, double threshold) {
		if (m1 == null || m2 == null) return null;
		Element.setSharedElements(m1, m2);
		ArrayList<PDBAtom> l1 = m1.getPDBAtomList();
		ArrayList<PDBAtom> l2 = m2.getPDBAtomList();
		if (l1 == null || l2 == null) return null;
		if (l1.size()<4 || l2.size()<4) return null;
		double angleDeltaThreshold = threshold*java.lang.Math.PI/180;
		ArrayList<Double>deltaAngle = new ArrayList<Double>();
		double sum = 0;
		ArrayList<Integer> al1 = new ArrayList<Integer>(), dal1 = new ArrayList<Integer>();
		ArrayList<Integer> al2 = new ArrayList<Integer>(), dal2 = new ArrayList<Integer>();
		PDBAtom a1, a2, a3, a4, b1, b2, b3, b4;
		for (int i:m1.getFreeAngleList()) {
			a1=l1.get(i-1);a2=l1.get(i);a3=l1.get(i+1);
			b1=l2.get(i-1);b2=l2.get(i);b3=l2.get(i+1);
			double angle1 = a1.subtract(a2).angle(a3.subtract(a2));
			double angle2 = b1.subtract(b2).angle(b3.subtract(b2));
			double d = java.lang.Math.abs(angle1-angle2);
			if (d>java.lang.Math.PI) d = java.lang.Math.PI*2 - d;
			if (d>angleDeltaThreshold) {
				al1.add(i);
				al2.add(i);
				sum += d;
				deltaAngle.add(sum);
			}
		}
		for (int i:m1.getFreeDihedralAngleList()) {
			a1=l1.get(i-1);a2=l1.get(i);a3=l1.get(i+1);a4=l1.get(i+2);
			b1=l2.get(i-1);b2=l2.get(i);b3=l2.get(i+1);b4=l2.get(i+2);
			double angle1 = calculateDihedralAngle(a1, a2, a3, a4);
			double angle2 = calculateDihedralAngle(b1, b2, b3, b4);
			double d = java.lang.Math.abs(angle1-angle2);
			if (d>java.lang.Math.PI) d = java.lang.Math.PI*2 - d;
			if (d>angleDeltaThreshold) {
				dal1.add(i);
				dal2.add(i);
				sum += d;
				deltaAngle.add(sum);
			}
		}
		m1.setFreeAngleList(al1);m1.setFreeDihedralAngleList(dal1);
		m2.setFreeAngleList(al2);m2.setFreeDihedralAngleList(dal2);
		return deltaAngle;
	}
	
	public void outputDihedralAngles() {
		ArrayList<PDBAtom> list = getPDBAtomList();
		int size=list.size();
		if (size<4) return;
		System.out.println("Dihedral angles:");
		PDBAtom a1, a2, a3, a4;
		a2=list.get(0);a3=list.get(1);a4=list.get(2);
		for (int i=3;i<size;i++) {
			a1=a2;a2=a3;a3=a4;
			a4=list.get(i);
			if (a1!=null && a2!=null && a3!=null && a4!=null)
				System.out.printf("%.2f ", 180*calculateDihedralAngle(a1, a2, a3, a4)/java.lang.Math.PI);
		}
	}
	
	/**
	 * dihedral type is defined by initial structure
	 */
	private void initDihedralTypes() {
		ArrayList<PDBAtom> list = getPDBAtomList();
		int size=list.size();
		if (size<4) return;	// no dihedral angle exists
		// define the dihedral type of first atom to be T
		dihedralTypeList=new ArrayList<DihedralTypes>();
		dihedralTypeList.add(DihedralTypes.T);
		PDBAtom a1, a2, a3, a4;
		a2=list.get(0);a3=list.get(1);a4=list.get(2);
		for (int i=3;i<size;i++) {
			a1=a2;a2=a3;a3=a4;
			a4=list.get(i);
			if (a1!=null && a2!=null && a3!=null && a4!=null)
				dihedralTypeList.add(GordonEnergy.getDihedralType(calculateDihedralAngle(a1, a2, a3, a4)));
		}
		// define the dihedral type of last 2 atoms to be T
		dihedralTypeList.add(DihedralTypes.T);
		dihedralTypeList.add(DihedralTypes.T);
	}

	/**
	 * 
	 * @param typeStr
	 * @param offset
	 */
	public void setDihedralTypes(String typeStr, int offset) {
		if (typeStr==null || offset<0 || (typeStr.length()+offset)>dihedralTypeList.size()) {
			System.out.println("Invalid input for setting dihedral types.");
			return;
		}
		StringCharacterIterator it = new StringCharacterIterator(typeStr);
		// Iterate over the characters in the forward direction
		for (char ch=it.first(); ch != StringCharacterIterator.DONE; ch=it.next()) {
			dihedralTypeList.set(offset, GordonEnergy.getDihedralType(ch));
			offset++;
		}
		elementList = null;
	}
	
	public void setDihedralTypes(ArrayList<DihedralTypes> list) {
		if (list.size() != dihedralTypeList.size()) {
			System.out.println("Invalid input size for setting dihedral types.");
			return;
		}
		dihedralTypeList = list;
		elementList = null;
	}
	
	public ArrayList<DihedralTypes> getDihedralTypes() {
		return dihedralTypeList;
	}
	
	/**
	 * 
	 * @return total energy
	 */
	public double getEnergy() {
		return getEnergy(GordonEnergy.angleType)+
				getEnergy(GordonEnergy.dihedralType)+
				getEnergy(GordonEnergy.vanDeWallsType)+
				getEnergy(GordonEnergy.hydrogenBondType);
	}

	public boolean isAngleEnergyValid() {
		ArrayList<PDBAtom> list = getPDBAtomList();
		int size=list.size();
		if (size<3 || size!=dihedralTypeList.size()) return true;
		DihedralTypes type;
		PDBAtom a1, a2, a3;
		a2=list.get(0);a3=list.get(1);
		for (int i=2;i<size;i++) {
			a1=a2;a2=a3;
			a3=list.get(i);
			type=dihedralTypeList.get(i-2);
			if (a1!=null && a2!=null && a3!=null) {
				double e = GordonEnergy.angleEnergy(a1.subtract(a2).angle(a3.subtract(a2)), type);
				if (e > GordonEnergy.angleEnergyThreshold) return false;
			}
		}
		return true;
	}
	
	private double getAngleEnergy() {
		double energy = 0;
		ArrayList<PDBAtom> list = getPDBAtomList();
		int size=list.size();
		if (size<3 || size!=dihedralTypeList.size()) return energy;
		DihedralTypes type;
		PDBAtom a1, a2, a3;
		a2=list.get(0);a3=list.get(1);
		for (int i=2;i<size;i++) {
			a1=a2;a2=a3;
			a3=list.get(i);
			type=dihedralTypeList.get(i-2);
			if (a1!=null && a2!=null && a3!=null) {
				double e = GordonEnergy.angleEnergy(a1.subtract(a2).angle(a3.subtract(a2)), type);
				energy += e;
			}
		}
		return energy;
	}
	
	public static double calculateDihedralAngle(PDBAtom a1, PDBAtom a2, PDBAtom a3, PDBAtom a4) {
		Vector v1v2 = a2.subtract(a1);
		Vector v2v3 = a3.subtract(a2);
		Vector v3v4 = a4.subtract(a3);

		Vector surface1Normal = v2v3.crossProduct(v3v4);
		Vector surface2Normal = v1v2.crossProduct(v2v3);

		//return surface2Normal.angle(surface1Normal);

		// The sign of the dihedral angle is determined by the angle between 
		// the normal to the plane v1v2 and the vector v3v4 (which is v3-v4). 
		// If the angle is < 90 degrees (or the dot product between the two 
		// vectors is negative), reverse the sign of the dihedral angle.
		double angle = surface2Normal.angle(surface1Normal);
		double sign = surface2Normal.dotProduct(v3v4);
		if (sign < 0)
		return (-1*angle);
		else 
		return angle;
		
	}
/*
	private static double calculateDihedralAngle2(PDBAtom a, PDBAtom b, PDBAtom c, PDBAtom d) {
		Vector ba=b.subtract(a);ba=ba.multiply(1.0/ba.magnitude());
		Vector cb=c.subtract(b);cb=cb.multiply(1.0/cb.magnitude());
		Vector dc=d.subtract(c);dc=dc.multiply(1.0/dc.magnitude());

		Vector n1=ba.crossProduct(cb);n1=n1.multiply(1.0/n1.magnitude());
		Vector n2=cb.crossProduct(dc);n2=n2.multiply(1.0/n2.magnitude());
		double ret=n2.angle(n1);
		if (n1.dotProduct(dc)<0) ret=-ret;
		return ret;
	}
	*/
	
	private double getDihedralEnergy() {
		double energy = 0;
		ArrayList<PDBAtom> list = getPDBAtomList();
		int size=list.size();
		if (size<4 || size!=dihedralTypeList.size()) return energy;
		DihedralTypes type;
		PDBAtom a1, a2, a3, a4;
		a2=list.get(0);a3=list.get(1);a4=list.get(2);
		for (int i=3;i<size;i++) {
			a1=a2;a2=a3;a3=a4;
			a4=list.get(i);
			type=dihedralTypeList.get(i-2);
			if (a1!=null && a2!=null && a3!=null && a4!=null) {
				energy += GordonEnergy.dihedralEnergy(calculateDihedralAngle(a1, a2, a3, a4), type);
			}
		}
		return energy;
	}
	
	public boolean isVanDeWallsEnergyValid() {
		ArrayList<PDBAtom> list = getPDBAtomList();
		int size=list.size();
		if (size<4) return true;
		PDBAtom a1, a2;
		for (int i=0;i<size-3;i++) {
			a1 = list.get(i);
			for (int j=i+3;j<size;j++) {
				a2 = list.get(j);
				double e = GordonEnergy.vanDeWallsEnergy(a1.distance(a2),
											GordonEnergy.getAtomType(a1.getAminoAcidType()),
											GordonEnergy.getAtomType(a2.getAminoAcidType()));
				if (e > GordonEnergy.vanDeWallsEnergyThreshold) return false;
			}
		}
		return true;
	}
	
	private double getVanDeWallsEnergy() {
		double energy = 0;
		ArrayList<PDBAtom> list = getPDBAtomList();
		int size=list.size();
		if (size<4) return energy;
		PDBAtom a1, a2;
		for (int i=0;i<size-3;i++) {
			a1 = list.get(i);
			for (int j=i+3;j<size;j++) {
				a2 = list.get(j);
				double e = GordonEnergy.vanDeWallsEnergy(a1.distance(a2),
											GordonEnergy.getAtomType(a1.getAminoAcidType()),
											GordonEnergy.getAtomType(a2.getAminoAcidType()));
				energy += e;
			}
		}
		return energy;
	}
	
	private double getHydrogenBondEnergy() {
		double energy = 0;
		ArrayList<PDBAtom> list = getPDBAtomList();
		int size=list.size();
		if (size<6 || size!=dihedralTypeList.size()) return energy;
		HydrogenBondCapability c1, c2;
		PDBAtom pre_a1, a1, post_a1, a2;
		a1=list.get(0);post_a1=list.get(1);
		for (int i=2;i<size-4;i++) {
			pre_a1=a1;a1=post_a1;
			post_a1 = list.get(i);
			c1 = GordonEnergy.getHydrogenBondCapability(dihedralTypeList.get(i-1));
			if (c1 == HydrogenBondCapability.A) {
				c2 = GordonEnergy.getHydrogenBondCapability(dihedralTypeList.get(i+2));
				if (c2 != HydrogenBondCapability.A) continue;
				a2 = list.get(i+2);
				PDBAtom pre_a2=list.get(i+1), post_a2=list.get(i+3);
				energy+=GordonEnergy.HBEnergy(pre_a1, a1, post_a1, pre_a2, a2, post_a2, c1);
			} else if (c1 == HydrogenBondCapability.B) {
				for (int j=i;j<size-1;j++) {
					c2 = GordonEnergy.getHydrogenBondCapability(dihedralTypeList.get(j));
					if (c2 != HydrogenBondCapability.B) continue;
					a2 = list.get(j);
					if (a1.distance(a2)>GordonEnergy.HBSheetSearchDistance) continue;
					PDBAtom pre_a2=list.get(j-1), post_a2=list.get(j+1);
					energy+=GordonEnergy.HBEnergy(pre_a1, a1, post_a1, pre_a2, a2, post_a2, c1);
				}
			}
		}
		return energy;
	}
	
	public double getEnergy(String type) {
		double energy=0;
		if (type.equalsIgnoreCase(GordonEnergy.angleType)) {
			energy = getAngleEnergy();
		} else if (type.equalsIgnoreCase(GordonEnergy.dihedralType)) {
			energy = getDihedralEnergy();
		} else if (type.equalsIgnoreCase(GordonEnergy.vanDeWallsType)) {
			energy = getVanDeWallsEnergy();
		} else if (type.equalsIgnoreCase(GordonEnergy.hydrogenBondType)) {
			energy = getHydrogenBondEnergy();
		}
		return energy;
	}

	@Override
	public String toString() {
		String nameStr = "", atomTypeStr = "", dihedralTypeStr = "", elementStr = "";
		ListIterator<PDBAtom> iter = getPDBAtomList().listIterator();
		while (iter.hasNext()) {
			PDBAtom a = iter.next();
			if (a != null) {
				nameStr += utilities.Definitions.getAminoAcidOneLetterName(a.getAminoAcidType());
				atomTypeStr += GordonEnergy.getAtomType(a.getAminoAcidType());
			}
		}
		ListIterator<DihedralTypes> iter2 = dihedralTypeList.listIterator();
		while (iter2.hasNext()) dihedralTypeStr += iter2.next();
		ListIterator<Element> iter3 = getElementList().listIterator();
		while (iter3.hasNext()) elementStr += iter3.next()+"; ";
		return "Chain "+getChainID()+":\n"+nameStr+"\n"+atomTypeStr+"\n"+dihedralTypeStr+"\n"+elementStr+"\n";
	}

	@Override
	public void translate(double deltaX, double deltaY, double deltaZ)
	{
		super.translate(deltaX, deltaY, deltaZ);
		for (Element e:getElementList()) e.inValidateVectors();
	}

	/**
	 * This function is used to correct the distance between two adjacent atoms
	 * caused by numeric error accumulated during sampling
	 * @param index is the index of first atom
	 * @param delta is the distance difference to be corrected
	 */
	public void translate(int index, double delta) {
		ArrayList<PDBAtom> l = getPDBAtomList();
		PDBAtom a1 = l.get(index), a2 = l.get(index+1);
		double ratio = delta/a1.distance(a2);
		double deltaX = (a2.getX()-a1.getX())*ratio;
		double deltaY = (a2.getY()-a1.getY())*ratio;
		double deltaZ = (a2.getZ()-a1.getZ())*ratio;
		for (int i=index+1;i<l.size();i++) l.get(i).translate(deltaX, deltaY, deltaZ);
		for (Element e:getElementList()) e.inValidateVectors();
	}
	
	/*
	 * only reset the vectors of elements after index
	 */
	public void resetElementVectors(int index) {
		for (Element e:getElementList())
			if ((e.getStart() >= index) || (e.getEnd() > index)) e.inValidateVectors();
	}
	
	/**
	 * 3-atoms angle rotation by a value
	 * @param index is where rotate is done
	 * @param angle is in degree
	 */
	public void rotateAroundPoint(int index, double angle)
	{
		ArrayList<PDBAtom> l = getPDBAtomList();
		if (index <= 0 || index >= l.size()-1) return;
		PDBAtom a1=l.get(index-1), a2=l.get(index), a3=l.get(index+1);
		Vector surfaceNormal = a2.subtract(a1).crossProduct(a3.subtract(a2));
		ListIterator<PDBAtom> iter = l.listIterator(index+1);
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotate(a2, a2.add(surfaceNormal), angle);
		}
		resetElementVectors(index);
	}
	
	/**
	 * 3-atoms angle rotation to a value
	 * @param index is where rotate is done
	 * @param angle is in degree
	 */
	public void rotateAroundPointTo(int index, double angle)
	{
		ArrayList<PDBAtom> l = getPDBAtomList();
		if (index <= 0 || index >= l.size()-1) return;
		PDBAtom a1=l.get(index-1), a2=l.get(index), a3=l.get(index+1);
		rotateAroundPoint(index, angle - a1.subtract(a2).angle(a3.subtract(a2))*180/java.lang.Math.PI);
	}
	
	/**
	 * @param bondBeginAtomIndex - the index of the atom at the beginning edge of the bond
 	 * @param bondEndAtomIndex - the index of the atom at the ending edge of the bond
	 * @param angle - the angle by which atoms will be rotated around the bond
	 */
	@Override
	public void rotateAroundBond(int bondBeginAtomIndex, int bondEndAtomIndex, double angle)
	{
		super.rotateAroundBond(bondBeginAtomIndex, bondEndAtomIndex, angle);
		resetElementVectors(bondEndAtomIndex);
	}
	
	/**
	 * @param bondBegin - the atom at the beginning edge of the bond
	 * @param bondEnd - the atom at the ending edge of the bond
	 * @param angle -
	 */
	@Override
	public  void rotateAroundBond(PDBAtom bondBegin, PDBAtom bondEnd, double angle)
	{
		super.rotateAroundBond(bondBegin, bondEnd, angle);
		resetElementVectors(getPDBAtomList().indexOf(bondEnd));
	}
	
	@Override
	public void rotateRigidBodyAroundArbitraryAxis(Vector begin, Vector end, double angle)
	{
		super.rotateRigidBodyAroundArbitraryAxis(begin, end, angle);
		resetElementVectors(0);
	}
	
	@Override
	public void rotateRigidBodyAroundArbitraryAxis(double alpha, double beta, double theta)
	{
		super.rotateRigidBodyAroundArbitraryAxis(alpha, beta, theta);
		resetElementVectors(0);
	}
	
	@Override
	public void rotateAroundXAxis(double angle)
	{
		super.rotateAroundXAxis(angle);
		resetElementVectors(0);
	}
	
	@Override
	public void rotateAroundYAxis(double angle)
	{
		super.rotateAroundYAxis(angle);
		resetElementVectors(0);
	}
	
	@Override
	public void rotateAroundZAxis(double angle)
	{
		super.rotateAroundZAxis(angle);
		resetElementVectors(0);
	}
	
	@Override
	public void rotateRigidBodyAroundXAxis(double angle)
	{
		super.rotateRigidBodyAroundXAxis(angle);
		// super method called other methods
		// nothing need do here for elements
	}
	
	@Override
	public void rotateRigidBodyAroundYAxis(double angle)
	{
		super.rotateRigidBodyAroundYAxis(angle);
		// super method called other methods
		// nothing need do here for elements
	}
	
	@Override
	public void rotateRigidBodyAroundZAxis(double angle)
	{
		super.rotateRigidBodyAroundZAxis(angle);
		// super method called other methods
		// nothing need do here for elements
	}

	public Vector getCentroid(int start, int end)
	{
		double totalX = 0;
		double totalY = 0;
		double totalZ = 0;

		ArrayList<PDBAtom> list = getPDBAtomList();
		if (start<0 || end>list.size() || end<start) return null;
		for (int i=start;i<=end;i++){
			PDBAtom a = list.get(i);
			totalX = totalX + a.getX();
			totalY = totalY + a.getY();
			totalZ = totalZ + a.getZ();
		}
		
		int chainSize = end-start+1;
		totalX /= chainSize;
		totalY /= chainSize;
		totalZ /= chainSize;
		
		return new Vector(totalX, totalY, totalZ);
	}

	/**
	 * Derive a new molecule with flexible 3-atoms and dihedral angles randomly sampled
	 * @return
	 */
	public CGMolecule newRandomMolecule() {
		CGMolecule molecule = null;
		Random g1=new Random(), g2=new Random();
		double minAngle = 50, range = 150 - 50, diRange = 360;
		do {
			molecule = CGMolecule.newCGMolecule(this);
			for (int i:getFreeAngleList())
				molecule.rotateAroundPointTo(i, g1.nextDouble()*range+minAngle);
			for (int i:getFreeDihedralAngleList())
				molecule.rotateAroundBond(i, i+1, g2.nextDouble()*diRange);
		} while (!molecule.isVanDeWallsEnergyValid());
		return molecule;
	}

	/**
	 * Derive a new molecule with flexible 3-atoms and dihedral angles randomly sampled
	 * @param deltaAngle contains the accumulated difference of angles
	 * @return
	 */
	public CGMolecule newRandomMolecule(ArrayList<Double> deltaAngle) {
		ArrayList<Integer> l1 = getFreeAngleList(), l2 = getFreeDihedralAngleList();
		int s1 = l1.size(), s2 = l2.size();
		if (deltaAngle == null || deltaAngle.size() != (s1+s2)) return null;
		CGMolecule molecule = null;
		Random g0=new Random(), g1=new Random(), g2=new Random();
		double minAngle = 50, range = 150 - 50, diRange = 360;
		ArrayList<Integer> toBeSampled = new ArrayList<Integer>();
		do {
			molecule = CGMolecule.newCGMolecule(this);
			//	pick angles to be sampled depending on their difference values stored in deltaAngle
			toBeSampled.clear();
			for (int i=0;i<s1+s2;i++) {
				double d = g0.nextDouble()*deltaAngle.get(s1+s2-1);
				int index = 0;
				for (;index<(s1+s2);index++)
					if (d < deltaAngle.get(index)) break;
				if (!toBeSampled.contains(index)) toBeSampled.add(index);
			}
			for (int i:toBeSampled) {
				/*
				if (i<s1)
					molecule.rotateAroundPointTo(l1.get(i), g1.nextDouble()*range+minAngle);
				else
					molecule.rotateAroundBond(l2.get(i-s1), l2.get(i-s1)+1, g2.nextDouble()*diRange);
					*/
				if (i<s1)
					molecule.rotateAroundPoint(l1.get(i), g1.nextDouble()*10-5);
				else
					molecule.rotateAroundBond(l2.get(i-s1), l2.get(i-s1)+1, g2.nextDouble()*10-5);
			}
		} while (!molecule.isVanDeWallsEnergyValid());
		/*
		String info = "";
		for (int i:toBeSampled) info += i+" ";
		System.out.println(info); */
		return molecule;
	}

	/**
	 * Derive a new molecule with flexible 3-atoms and dihedral angles randomly sampled
	 * @return
	 */
	public CGMolecule newDerivedMolecule() {
		CGMolecule molecule = null;
		ArrayList<Integer> l = getFreeAngleList();
		int s1 = l.size();
		Random g1=new Random(), g2=new Random(), g3=new Random();
		int iter = 0;
		do {
			if (iter++>1000) System.out.println("No new derived molecule can be generate after 1000 trys");
			molecule = CGMolecule.newCGMolecule(this);
			int i = g1.nextInt(s1+getFreeDihedralAngleList().size());
			if (i < s1)
				molecule.rotateAroundPoint(l.get(i), g2.nextDouble()*GordonEnergy.sampleAngleRange-GordonEnergy.sampleAngleRange/2);
			else
				molecule.rotateAroundBond(getFreeDihedralAngleList().get(i-s1), getFreeDihedralAngleList().get(i-s1)+1, g3.nextDouble()*GordonEnergy.sampleAngleRange-GordonEnergy.sampleAngleRange/2);
		} while (!molecule.isAngleEnergyValid() || !molecule.isVanDeWallsEnergyValid());
		if (iter > 1000) return null;
		return molecule;
	}
/*
	public ArrayList<Double> getAngleVector() {
		return angleVector;
	}
	
	public ArrayList<Double> getDistVector() {
		return distVector;
	}
*/	
	/**
	 * get score vector relative to given molecule
	 * @param reference
	 * @return
	 */
	public ArrayList<Double> getScoreVector(CGMolecule reference) {
		ArrayList<Double> scoreVector = null;
		if (reference == null || reference.getElementList() == null || getElementList() == null) return null;
		ArrayList<Element> l1 = reference.elementList, l2 = elementList;
		if (l1.size() != l2.size()) return null;
		int len = l1.size(), size = len*(len-1)/2;
		double angles1[]=new double[size], angles2[]=new double[size];
		double dists1[]=new double[size], dists2[]=new double[size];
		int index = 0;
		for (int i=0;i<len;i++) {
			Element ei1=l1.get(i), ei2=l2.get(i);
			for (int j=i+1;j<len;j++) {
				Element ej1=l1.get(j), ej2=l2.get(j);
				angles1[index] = Element.Angle(ei1, ej1);
				angles2[index] = Element.Angle(ei2, ej2);
				dists1[index] = Element.distance(ei1, ej1);
				dists2[index] = Element.distance(ei2, ej2);
				index++;
			}
		}
		scoreVector = new ArrayList<Double>();
		//angleVector = new ArrayList<Double>(); distVector = new ArrayList<Double>();
		for (int i=0;i<len;i++) {
			double s = 0;
			//double sa = 0, sd = 0;
			for (int j=0;j<len;j++) {
				if (j==i) continue;
				if (j<i)
					index = (len-1)*j-j*(j-1)/2+i-j-1;
				else
					index = (len-1)*i-i*(i-1)/2+j-i-1;
				double d = java.lang.Math.abs(angles1[index]-angles2[index])*GordonEnergy.weightAngle;
				s += d;
				//sa += d;
				d = java.lang.Math.abs(dists1[index]-dists2[index])*GordonEnergy.weightDistance;
				s += d;
				//sd += d;
			}
			scoreVector.add(s);
			//angleVector.add(sa); distVector.add(sd);
		}
		return scoreVector;
	}

	public static double getScore(ArrayList<Double> scoreVector) {
		if (scoreVector != null && scoreVector.size() > 0) {
			double sum = 0;
			for (double s:scoreVector) sum += s*s;
			return java.lang.Math.sqrt(sum);
		}
		return -1;
	}
	
	public double getScore(CGMolecule reference) {
		return getScore(getScoreVector(reference));
	}

	public static double squareDistance(ArrayList<Double> sv1, ArrayList<Double> sv2) {
		if (sv1 == null || sv2 == null || sv1.size() != sv2.size()) return -1;
		double sum = 0;
		for (int i=0;i<sv1.size();i++)
			sum += (sv1.get(i)-sv2.get(i))*(sv1.get(i)-sv2.get(i));
		return sum;
	}
	
	public static double distance(ArrayList<Double> sv1, ArrayList<Double> sv2) {
		if (sv1 == null || sv2 == null || sv1.size() != sv2.size()) return -1;
		return java.lang.Math.sqrt(squareDistance(sv1, sv2));
	}
	
	public static double squareDistance(CGMolecule m1, CGMolecule m2, CGMolecule reference) {
		if (m1 == null || m2 == null) return -1;
		return squareDistance(m1.getScoreVector(reference), m2.getScoreVector(reference));
	}

	public static double distance(CGMolecule m1, CGMolecule m2, CGMolecule reference) {
		if (m1 == null || m2 == null) return -1;
		return distance(m1.getScoreVector(reference), m2.getScoreVector(reference));
	}
	
	public double getAngle(int index) {
		ArrayList<PDBAtom> l = getPDBAtomList();
		PDBAtom a1, a2, a3;
		a1=l.get(index-1);a2=l.get(index);a3=l.get(index+1);
		return a1.subtract(a2).angle(a3.subtract(a2));
	}
	
	public double getDihedralAngle(int index) {
		ArrayList<PDBAtom> l = getPDBAtomList();
		return CGMolecule.calculateDihedralAngle(l.get(index-1), l.get(index), l.get(index+1), l.get(index+2));
	}
}
