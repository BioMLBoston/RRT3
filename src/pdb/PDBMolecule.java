package pdb;

import java.util.ArrayList;
import java.util.ListIterator;
import utilities.Definitions;
import utilities.Vector;

/**
 * 
 * @author Bahar Akbal-Delibas (abakbal@cs.umb.edu)
 * This class represents a single protein chain
 *
 */
public class PDBMolecule 
{
	///////////////////////////////////////////
	//// Constructors of PDBMolecule class ////
	///////////////////////////////////////////
	public PDBMolecule(String theChainID)
	{
		this.chainID = theChainID;
		this.pdbAtomList = new ArrayList<PDBAtom>(); // create an empty pdbAtomList
	}
	
	public PDBMolecule(ArrayList<PDBAtom> theAtomList, String theChainID)
	{
		this.pdbAtomList = theAtomList;
		this.chainID = theChainID;
	}
	
	//////////////////////////////////////////////////
	//// Setters and getters of PDBMolecule class ////
	//////////////////////////////////////////////////
	
	/**
	 * @param theAtomList
	 */
	public void setPDBAtomList(ArrayList<PDBAtom> theAtomList) 
	{
		this.pdbAtomList = theAtomList;
	}
	
	/**
	 * @return the pdbAtomList
	 */
	public ArrayList<PDBAtom> getPDBAtomList() {
		return this.pdbAtomList;
	}
	
	/**
	 * 
	 * @return the chain ID of the PDB molecule
	 */
	public String getChainID()
	{
		return this.chainID;
	}
	
	//////////////////////////////////////////
	//// Operations of PDBMolecule class  ////
	//////////////////////////////////////////
	
	/**
	 * @return the backbone atoms of this molecule
	 */
	public ArrayList<PDBAtom> extractBackbone()
	{
		// create an iterator on this chain
		ListIterator<PDBAtom> iter = this.getPDBAtomList().listIterator(0); 
		int chainSize = this.getPDBAtomList().size();
		// create a list of backbone atoms
		ArrayList<PDBAtom> backbone = new ArrayList<PDBAtom>();		
		while(iter.hasNext() && iter.nextIndex() <= chainSize)
		{
			PDBAtom a = iter.next();
			if(a.isCA() || a.isN() || a.isC() || a.isO())
				backbone.add(a);
		}
		return backbone;
	}
	
	/**
	 * This method returns the list of C-alpha atoms in the molecule. 
	 * The list of C-alpha atoms can be used as a coarse-grained representation of the molecule.
	 * @return the C-alpha atoms of this molecule
	 */
	public ArrayList<PDBAtom> extractCAlphas()
	{
		// create an iterator on this chain
		ListIterator<PDBAtom> iter = this.getPDBAtomList().listIterator(0); 
		int chainSize = this.getPDBAtomList().size();
		// create a list of C-alpha atoms
		ArrayList<PDBAtom> cAlphas = new ArrayList<PDBAtom>();		
		while(iter.hasNext() && iter.nextIndex() <= chainSize)
		{
			PDBAtom a = iter.next();
			if(a.isCA())
				cAlphas.add(a);
		}
		return cAlphas;
	}
	
	/**
	 * This method computes the centroid of a 3D point set (i.e. a molecule)
	 * @return a Vector which holds the x, y, z coordinates of the centroid
	 */
	public Vector getCentroid()
	{
		double totalX = 0;
		double totalY = 0;
		double totalZ = 0;
		Vector centroid = new Vector();
		
		ListIterator<PDBAtom> iter = this.getPDBAtomList().listIterator(0); 
		int chainSize = this.getPDBAtomList().size();	
		while(iter.hasNext() && iter.nextIndex() <= chainSize)
		{
			PDBAtom a = iter.next();
			totalX = totalX + a.getX();
			totalY = totalY + a.getY();
			totalZ = totalZ + a.getZ();
		}
		
		centroid.setX(totalX/chainSize);
		centroid.setY(totalY/chainSize);
		centroid.setZ(totalZ/chainSize);
		
		return centroid;
	}
	
	/**
	 * This method translates this molecule (i.e. chain) by calling the 
	 * translate() method on each of its atoms with the given delta values
	 * on x-coordinate, y-coordinate and z-coordinate. 
	 */
	public void translate(double deltaX, double deltaY, double deltaZ)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.translate(deltaX, deltaY, deltaZ);
		}
	}
	
	
	/**
	 * @param bondBeginAtomIndex - the index of the atom at the beginning edge of the bond
 	 * @param bondEndAtomIndex - the index of the atom at the ending edge of the bond
	 * @param angle - the angle by which atoms will be rotated around the bond
	 */
	public void rotateAroundBond(int bondBeginAtomIndex, int bondEndAtomIndex, double angle)
	{
		// Specify the bond by the indices of two atoms at the beginning and ending edges
		// This is needed to perform the rotation on atoms to be rotated.
		PDBAtom bondBegin = this.pdbAtomList.get(bondBeginAtomIndex);
		PDBAtom bondEnd = this.pdbAtomList.get(bondEndAtomIndex);
		
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator(bondEndAtomIndex+1);
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotate(bondBegin, bondEnd, angle);
		}
	}
	
	
	/**
	 * @param bondBegin - the atom at the beginning edge of the bond
	 * @param bondEnd - the atom at the ending edge of the bond
	 * @param angle -
	 */
	public  void rotateAroundBond(PDBAtom bondBegin, PDBAtom bondEnd, double angle)
	{
		// Specify the bond by two atoms at the beginning and ending edges
		// This is needed to perform the rotation on atoms to be rotated.
		int bondEndAtomIndex = this.pdbAtomList.indexOf(bondEnd);
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator(bondEndAtomIndex+1);
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotate(bondBegin, bondEnd, angle);
		}
	}
	
	public void rotateRigidBodyAroundArbitraryAxis(Vector begin, Vector end, double angle)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotate(begin, end, angle);
		}
	}
	
	public void rotateRigidBodyAroundArbitraryAxis(double alpha, double beta, double angle)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotate(alpha, beta, angle);
		}
	}
	
	public void rotateAroundXAxis(double angle)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotateAroundXAxis(angle);
		}
	}
	
	public void rotateAroundYAxis(double angle)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotateAroundYAxis(angle);
		}
	}
	
	public void rotateAroundZAxis(double angle)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotateAroundZAxis(angle);
		}
	}
	
	public void rotateRigidBodyAroundXAxis(double angle)
	{
		// First translate the centroid of the molecule to the origin
		Vector centroid = this.getCentroid();
		this.translate(-1*centroid.getX(), -1*centroid.getY(), -1*centroid.getZ());
		
		// Then rotate the molecule around x axis
		this.rotateAroundXAxis(angle);
		
		// Finally translate the centroid back to original position
		this.translate(centroid.getX(), centroid.getY(), centroid.getZ());
	}
	
	public void rotateRigidBodyAroundYAxis(double angle)
	{
		// First translate the centroid of the molecule to the origin
		Vector centroid = this.getCentroid();
		this.translate(-1*centroid.getX(), -1*centroid.getY(), -1*centroid.getZ());
		
		// Then rotate the molecule around y axis
		this.rotateAroundYAxis(angle);
		
		// Finally translate the centroid back to original position
		this.translate(centroid.getX(), centroid.getY(), centroid.getZ());
	}
	
	public void rotateRigidBodyAroundZAxis(double angle)
	{
		// First translate the centroid of the molecule to the origin
		Vector centroid = this.getCentroid();
		this.translate(-1*centroid.getX(), -1*centroid.getY(), -1*centroid.getZ());
		
		// Then rotate the molecule around z axis
		this.rotateAroundZAxis(angle);
		
		// Finally translate the centroid back to original position
		this.translate(centroid.getX(), centroid.getY(), centroid.getZ());
	}
	
	/**
	 * Checks if there are any steric clashes within the given protein conformation. 
	 * @return true if clash exists between atoms of the molecule
	 * @return false if no clash exists between atoms of the molecule
	 */
	public boolean stericClashExists()
	{
		PDBAtom a1, a2;
		boolean stericClashOccurred = false;
		double energy = 0;
		
		for(int i = 0 ; i < pdbAtomList.size() ; i++)
		{	
			a1 = pdbAtomList.get(i);
			
			// check if pdbAtom at index i clashes with any other atom in the list.
			for(int j = i+4 ; j < pdbAtomList.size() ; j++)
			{
				a2 = pdbAtomList.get(j);
				
				double dist = a1.distance(a2);
				double a1Radius = a1.getRmin2();
				double a2Radius = a2.getRmin2();
				energy = energy + a1.pairwiseLJPotential(a2);
				if( (dist < a1Radius + a2Radius) && !isClashIgnorable(a1, a2) )
				{
					System.out.println("Steric clash occurred between atom " + i + " and atom " + j);
					System.out.println("Distance between these two atoms: " + dist);
					stericClashOccurred = true;
				}
			}
		}
					
		// if the control reaches here, this means that
		// no clashes detected.
		return stericClashOccurred;
	}
	
	/**
	 * This method checks if the steric clash between a pair of atoms is small enough to be ignored.
	 * @param a1 first atom of the clashing pair
	 * @param a2 second atom of the clashing pair
	 * @return true if the distance between clashing atoms is greater than 0.7 Angstroms
	 * @return false otherwise
	 */
	private static boolean isClashIgnorable(PDBAtom a1, PDBAtom a2){
		// The threshold value of 0.7 for detecting clashes is set by experimenting 
		// with the backbone of the 1COA structure of CI2, which is provided with 
		// Lydia Kavraki's Assignment-2.
		// This distance value can be improved by more trials and observations.
		double dist = a1.distance(a2);
		if(dist > 0.7)
			return true;
		return false;
	}
	
	/**
	 * This method gets the atom with the specified atom type from the specified amino acid.
	 * @param atomTypeStr the type of the desired atom, e.g "C", "CA", "N", "O", etc.
	 * @param theResidueIndex the index of the amino acid within the molecule, which the desired atom belongs to. 
	 * @return PDBAtom with the specified name and residue index
	 * @return null if the molecule does not contain such atom
	 */
	public PDBAtom getPDBAtomInMolecule(String atomTypeStr, int theResidueIndex)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			if(atom.getAtomType() == Definitions.getAtomType(atomTypeStr)
					&& atom.getResidueIndex() == theResidueIndex)
				return atom;
		}
		// if the control reaches here, this means that pdbAtomList
		// does not contain the given atom so return null
		return null;
	}
	
	/**
	 * This method gets the list of atoms with the specified atom type from the specified amino acid type.
	 * @param atomTypeStr the desired type of the list of atoms, e.g "C", "CA", "N", "O", etc.
	 * @param aminoTypeStr the desired type of the amino acids that contain the atoms 
	 * @return ArrayList<PDBAtom> an array list of atoms with the specified name and residue type
	 * @return null if the molecule does not contain such atom
	 */
	public ArrayList<PDBAtom> getPDBAtomsInMolecule(String atomTypeStr, String aminoTypeStr)
	{
		ArrayList<PDBAtom> outputList = new ArrayList<PDBAtom>();
		
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			if(atom.getAtomType() == Definitions.getAtomType(atomTypeStr)
					&& atom.getAminoAcidType() == Definitions.getAminoAcidType(aminoTypeStr))
				outputList.add(atom);
		}
		// if the control reaches here, this means that pdbAtomList
		// does not contains the given atom so return null
		return outputList;
	}
	
	/**
	 * This method gets the list of atoms with the specified AMBER atom type.
	 * @param amberTypeStr the desired AMBER type of the list of atoms.
	 * @return ArrayList<AmberPDBAtom> an array list of atoms with the specified 
	 * name and residue type (null if the molecule does not contain such atom)
	 */
	public ArrayList<AmberPDBAtom> getAmberPDBAtomsInMolecule(String amberTypeStr)
	{
		ArrayList<AmberPDBAtom> outputList = new ArrayList<AmberPDBAtom>();
		
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			if(atom instanceof AmberPDBAtom)
			{
				AmberPDBAtom amberAtom = (AmberPDBAtom) atom;
				if(amberAtom.getAmberAtomType() == Definitions.getAmberAtomType(amberTypeStr)){
					outputList.add(amberAtom);
				}
			}
		}
		// if the control reaches here, this means that amberPDBAtomList
		// does not contains the given atom so return null
		return outputList;
	}
	
	/**
	 * This method adds a new pdbAtom to the pdbAtomList of this PDBMolecule.
	 * @param pdbAtom the new atom to be added
	 */
	public void addPDBAtom(PDBAtom pdbAtom)
	{
		this.pdbAtomList.add(pdbAtom);
	}
	
	/**
	 * This method gives the index of a specific PDBAtom that belongs to this PDBMolecule.
	 * @param pdbAtom
	 * @return index of the given PDBAtom; -1 if the given pdbAtom does not exist in this molecule
	 */
	public int getAtomIndex(PDBAtom pdbAtom)
	{
		for(int i = 0; i < this.pdbAtomList.size(); i++)	
		{
			PDBAtom atom = this.pdbAtomList.get(i);
			if(atom.equals(pdbAtom))
				return i;
		}
		// otherwise pdbAtom does not exist in this molecule
		// so return -1
		return -1;
	}
	
	/**
	 * This method checks if the given pair atoms are bonded or not by looking at their indices.
	 * @param a1 first PDBAtom of the given pair	
	 * @param a2 second PDBAtom of the given pair
	 * @return true if the given pair of atoms are bonded; false otherwise.
	 */
	public boolean isBonded(PDBAtom a1, PDBAtom a2)
	{
		int indexOfA1 = this.getAtomIndex(a1);
		int indexOfA2 = this.getAtomIndex(a2);
		
		if(indexOfA1 == -1 || indexOfA2 == -1)
		{
			return false;
		}
		else
		{
			if(Math.abs(indexOfA1 - indexOfA2) < 4)
			{
				return true;	
			}
			else
			{
				return false;
			}
		}
	}
	
	/**
	 * This method checks if the given pair of atoms are bonded with a peptide bond.
	 * @param a1 first PDBAtom of the given pair
	 * @param a2 first PDBAtom of the given pair
	 * @return true if the given pair of atoms share a peptide bond; false otherwise.
	 */
	public boolean isPeptideBond(PDBAtom a1, PDBAtom a2)
	{
		if(a1.getChainID().equalsIgnoreCase(this.chainID) &&
				a2.getChainID().equalsIgnoreCase(this.chainID))
		{
			if (a1.getResidueIndex() == a2.getResidueIndex()+1)
			{
				if( a1.isN() && a2.isC() )
					return true;
				return false;
			}
			else if (a2.getResidueIndex() == a1.getResidueIndex()+1)
			{
				if( a1.isC() && a2.isN() )
					return true;
				return false;
			}
			else
				return false;
		}
		else
			return false;
	}
	
	public void setTerminalTypes()
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			if(this.containsResidueWithSmallerIndex(atom.getResidueIndex()) &&
					this.containsResidueWithLargerIndex(atom.getResidueIndex()))
			{
				// then the chain contains residue(s) both with an index
				// larger than the atom's residue index and smaller than
				// the atom's residue index
				atom.setTerminalType(Definitions.TerminalTypes.NONTERMINAL);
			}
			else 
			{
				if(this.containsResidueWithSmallerIndex(atom.getResidueIndex()))
				{
					// then the chain contains residue(s) with an index
					// smaller than the atom's residue index but it does not
					// contain any residues with an index larger than the 
					// atom's residue index
					atom.setTerminalType(Definitions.TerminalTypes.CTERMINAL);
				}
				else
				{
					// then the chain contains residue(s) with an index
					// larger than the atom's residue index but it does not
					// contain any residues with an index smaller than the 
					// atom's residue index
					atom.setTerminalType(Definitions.TerminalTypes.NTERMINAL);
				}
			}
		}
	}
	
	private boolean containsResidueWithSmallerIndex(int theResidueIndex)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			if(atom.getResidueIndex() < theResidueIndex)
				return true;
		}
		return false;
	}
	
	private boolean containsResidueWithLargerIndex(int theResidueIndex)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			if(atom.getResidueIndex() > theResidueIndex)
				return true;
		}
		return false;
	}

	
	/////////////////////////////////////////////
	//// The attributes of PDBMolecule class ////
	/////////////////////////////////////////////
	
	/**
	 * The list of PDB atoms in the molecule
	 */
	private ArrayList<PDBAtom> pdbAtomList; 
	
	/**
	 * The chain ID of the PDBMolecule
	 */
	private String chainID;

}
