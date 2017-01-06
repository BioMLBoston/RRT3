package cg;

/**
 * @author Dong Luo 2012-05-24
 *
 */

import io.Parser;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.ListIterator;

import cg.GordonEnergy.DihedralTypes;

import pdb.PDBMolecule;
import pdb.PDBMoleculeComplex;

public class CGMoleculeComplex extends PDBMoleculeComplex {

	private CGMoleculeComplex(ArrayList<PDBMolecule> theList) {
		super(theList);
	}
	
	public static CGMoleculeComplex newCGMoleculeComplex(PDBMoleculeComplex moleculeComplex) {
		ArrayList<PDBMolecule> moleculeList=new ArrayList<PDBMolecule>();
		ListIterator<PDBMolecule> iter = moleculeComplex.getMolecules().listIterator();
		while(iter.hasNext())
		{
			PDBMolecule pdbMolecule = iter.next();
			if(pdbMolecule != null)
				moleculeList.add(CGMolecule.newCGMolecule(pdbMolecule));
		}
		return new CGMoleculeComplex(moleculeList);
	}
	
	public static CGMoleculeComplex newCGMoleculeComplex(String filePath) {
		CGMoleculeComplex inst = newCGMoleculeComplex(Parser.parsePDBFile(filePath));
		inst.secondaryStructureFromPDBFile(filePath);
		return inst;
	}
	
	/**
	 * dihedral types corrected by PDB SS section but ignore first and last residues
	 * rigid element re-assign is automatically done
	 * TURN SS is ignored
	 * PDB SS Format reference:
	 * http://www.wwpdb.org/documentation/format23/sect5.html
	 * @param filePath
	 */
	public void secondaryStructureFromPDBFile(String filePath){
		try{
			BufferedReader reader = new BufferedReader(new FileReader(filePath));
			String pdbLine;
			CGMolecule old = null, molecule = null;
			ArrayList<DihedralTypes> dtl = null;
			while( (pdbLine = reader.readLine()) != null){
				if (pdbLine.length() < 4) continue;
				if (pdbLine.substring(0,6).equals("HELIX ")) {
					molecule = (CGMolecule) getMolecule(pdbLine.substring(19, 20));
					if (molecule == null) continue;
					if (molecule != old) {
						if (old != null && dtl != null) old.setDihedralTypes(dtl);
						old = molecule;
						dtl = molecule.getDihedralTypes();
					}
					int offset = Integer.parseInt(pdbLine.substring(21,25).trim())
							-molecule.getPDBAtomList().get(0).getResidueIndex();
					int length = Integer.parseInt(pdbLine.substring(71,76).trim());
					if (length >= Element.minimumElementSize) {
						for (int i=offset+1;i<offset+length-1;i++)
						//for (int i=offset;i<offset+length-1;i++)
							dtl.set(i, DihedralTypes.H);
					}
				} else if (pdbLine.substring(0,6).equals("SHEET ")) {
					molecule = (CGMolecule) getMolecule(pdbLine.substring(21, 22));
					if (molecule == null) continue;
					if (molecule != old) {
						if (old != null && dtl != null) old.setDihedralTypes(dtl);
						old = molecule;
						dtl = molecule.getDihedralTypes();
					}
					int offset = Integer.parseInt(pdbLine.substring(22,26).trim())-1;
					int length = Integer.parseInt(pdbLine.substring(33,37).trim())-offset;
					offset -= molecule.getPDBAtomList().get(0).getResidueIndex()-1;
					if (length >= Element.minimumElementSize) {
						for (int i=offset+1;i<offset+length-1;i++)
						//for (int i=offset;i<offset+length-1;i++)
							dtl.set(i, DihedralTypes.E);
					}
				}
			}
			// last molecule redefine dihedral types
			if (old != null && dtl != null) old.setDihedralTypes(dtl);
			
			reader.close();
		}
		catch(FileNotFoundException e){
			System.out.println("Cannot find file " + filePath);
		}
		catch(IOException e){
			System.out.println("Error occured while reading from file " + filePath);
		}
	}
	
	/**
	 * pair wise energy is only calculated within each chain for now
	 * ToDo: define the way to calculate inter-chain energy
	 */
	public double getEnergy() {
		double energy=0;
		ListIterator<PDBMolecule> iter = getMolecules().listIterator();
		while(iter.hasNext())
		{
			CGMolecule molecule = (CGMolecule) iter.next();
			if(molecule != null) energy+=molecule.getEnergy();
		}
		return energy;
	}

	/**
	 * pair wise energy is only calculated within each chain for now
	 * ToDo: define the way to calculate inter-chain energy
	 */
	public double getEnergy(String type) {
		double energy=0;
		ListIterator<PDBMolecule> iter = getMolecules().listIterator();
		while(iter.hasNext())
		{
			CGMolecule molecule = (CGMolecule) iter.next();
			if(molecule != null) energy+=molecule.getEnergy(type);
		}
		return energy;
	}

	public double getMoleculeEnergy(String chainID) {
		CGMolecule molecule=(CGMolecule) getMolecule(chainID);
		if (molecule == null) return 0;
		return molecule.getEnergy();
	}

	public double getMoleculeEnergy(String type, String chainID) {
		CGMolecule molecule=(CGMolecule) getMolecule(chainID);
		if (molecule == null) return 0;
		return molecule.getEnergy(type);
	}

	@Override
	public String toString() {
		String info = new String();
		ListIterator<PDBMolecule> iter = getMolecules().listIterator();
		while (iter.hasNext()) info += (CGMolecule) iter.next();
		return info;
	}
}
