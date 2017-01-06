package pdb;
import utilities.Definitions.AminoAcidTypes;
import utilities.Definitions.AtomTypes;
import utilities.Definitions.TerminalTypes;
import utilities.Definitions;

/**
 * @author Bahar Akbal-Delibas (abakbal@cs.umb.edu)
 *
 */
public class PDBAtom extends Atom {

	///////////////////////////////////////
	//// Constructors of PDBAtom class ////
	///////////////////////////////////////
	
	/**
	 * Default constructor of PDBAtom class
	 */
	public PDBAtom(){
		super();
		this.setAtomType(AtomTypes.NONE);
		this.setAminoAcidType(AminoAcidTypes.UNK);
		this.setResidueIndex(0);
		// set eps and rmin2
		assignLJParams();
	}
	
	/**
	 * Constructor of PDBAtom class
	 * @param theIndex the pdb file atom index of this atom
	 * @param theX the X-coordinate of the atom in space
	 * @param theY the Y-coordinate of the atom in space
	 * @param theZ the Z-coordinate of the atom in space
	 * @param theAtomTypeStr String value of the type of this atom (e.g. "CA", "C", "N", "O", etc.)
	 * @param theAminoTypeStr String value of the type of the residue that this atom is part of (e.g. "ALA", "LEU", "MET", etc.)
	 * @param theResidueIndex the pdb file residue index of the residue that this atom is part of
	 */
	public PDBAtom(int theIndex, double theX, double theY, double theZ, 
			String theAtomTypeStr, String theAminoTypeStr, String theChainID, int theResidueIndex) {
		super(theIndex, theX, theY, theZ);
		this.setAtomType(Definitions.getAtomType(theAtomTypeStr));
		this.setAminoAcidType(Definitions.getAminoAcidType(theAminoTypeStr));
		this.setChainID(theChainID);
		this.setResidueIndex(theResidueIndex);
		// set eps and rmin2
		assignLJParams();
	}
	
	/**
	 * Constructor of PDBAtom class
	 * @param theIndex the pdb file atom index of this atom
	 * @param theX the X-coordinate of the atom in space
	 * @param theY the Y-coordinate of the atom in space
	 * @param theZ the Z-coordinate of the atom in space
	 * @param theAtomType the type of this atom (e.g. CA, C, N, O, etc.)
	 * @param theAminoType the type of the residue that this atom is part of (e.g. ALA, LEU, MET, etc.)
	 * @param theChainID the String value of the chain ID that this atom belongs to (e.g. "A", "B", "C", etc.)
	 * @param theResidueIndex the pdb file residue index of the residue that this atom is part of
	 */
	public PDBAtom(int theIndex, double theX, double theY, double theZ, 
			AtomTypes theAtomType, AminoAcidTypes theAminoType, String theChainID, int theResidueIndex) {
		super(theIndex, theX, theY, theZ);
		this.setAtomType(theAtomType);
		this.setAminoAcidType(theAminoType);
		this.setChainID(theChainID);
		this.setResidueIndex(theResidueIndex);
		// set eps and rmin2
		assignLJParams();
	}
	
	public PDBAtom(double theX, double theY, double theZ){
		super();
		this.setX(theX);
		this.setY(theY);
		this.setZ(theZ);
		this.setAtomType(AtomTypes.NONE);
		this.setAminoAcidType(AminoAcidTypes.UNK);
		this.setResidueIndex(0);
		// set eps and rmin2
		assignLJParams();
	}
	
	///////////////////////////////////////////
	//// Setters and getters PDBAtom class ////
	///////////////////////////////////////////
	
	/**
	 * @param rmin2
	 */
	public void setRmin2(double rmin2) {
		this.rmin2 = rmin2;
	}
	
	/**
	 * @return the rmin2
	 */
	public double getRmin2() {
		return rmin2;
	}

	/**
	 * @param atomType the atomType to set
	 */
	public void setAtomType(AtomTypes atomType) {
		this.atomType = atomType;
	}

	/**
	 * @return the atomType
	 */
	public AtomTypes getAtomType() {
		return atomType;
	}

	/**
	 * @param aminoAcidType the aminoAcidType to set
	 */
	public void setAminoAcidType(AminoAcidTypes aminoAcidType) {
		this.aminoAcidType = aminoAcidType;
	}

	/**
	 * @return the aminoAcidType
	 */
	public AminoAcidTypes getAminoAcidType() {
		return aminoAcidType;
	}
	
	/**
	 * @param theChainID the chain ID to set
	 */
	public void setChainID(String theChainID)
	{
		this.chainID = theChainID;
	}
	
	/**
	 * @return the chainID
	 */
	public String getChainID()
	{
		return this.chainID;
	}
	
	/**
	 * @param eps
	 */
	public void setEps(double eps) {
		this.eps = eps;
	}
	
	/**
	 * @return the eps
	 */
	public double getEps() {
		return eps;
	}

	/**
	 * @param i the residueIndex to set
	 */
	public void setResidueIndex(int i) {
		this.residueIndex = i;
	}

	/**
	 * @return the residueIndex
	 */
	public int getResidueIndex() {
		return residueIndex;
	}
	
	public void setTerminalType(TerminalTypes theTerminalType)
	{
		this.terminalType = theTerminalType;
	}
	
	public TerminalTypes getTerminalType()
	{
		return this.terminalType;
	}
	
	/////////////////////////////////////////
	//// Other methods of PDB Atom class ////
	/////////////////////////////////////////
	
	/**
	 * @return true if the atom is Oxygen; false otherwise
	 */
	public boolean isO(){
		if(this.atomType == AtomTypes.O || this.atomType == AtomTypes.OXT || 
				this.atomType == AtomTypes.OT1 || this.atomType == AtomTypes.OT2)
			return true;
		else
			return false;
	}
	
	/**
	 * @return true if the atom is Nitrogen; false otherwise
	 */
	public boolean isN(){
		if(this.atomType == AtomTypes.N)
			return true;
		else
			return false;
	}
	
	/**
	 * @return true if the atom is backbone Carbon; false otherwise
	 */
	public boolean isC(){
		if(this.atomType == AtomTypes.C)
			return true;
		else
			return false;
	}
	
	/**
	 * @return true if the atom is backbone C-alpha; false otherwise
	 */
	public boolean isCA(){
		if(this.atomType == AtomTypes.CA)
			return true;
		else
			return false;
	}
	
	/**
	 * @return true if the atom is a Hydrogen atom; false otherwise
	 */
	public boolean isH(){
		if(Definitions.getAtomTypeString(this.atomType).trim().startsWith("H"))
			return true;
		else
			return false;
	}
	
	/**
	 * @return true if the atom is backbone C-beta; false otherwise
	 */
	public boolean isCB(){
		if(this.atomType == AtomTypes.CB)
			return true;
		else
			return false;
	}
	
	/**
	 * @return true if the atom is backbone Carbon or C-beta atom; false otherwise
	 */
	public boolean extBackBone(){
		return (isO() || isN() || isC() || isCA() || isCB());
	}
	
	/**
	 * @return true if 
	 */
	public boolean isSideChainConsideredForHBond(){
		if(this.aminoAcidType == AminoAcidTypes.SER || this.aminoAcidType == AminoAcidTypes.THR || 
				this.aminoAcidType == AminoAcidTypes.ASN || this.aminoAcidType == AminoAcidTypes.ASP ||
				this.aminoAcidType == AminoAcidTypes.GLN || this.aminoAcidType == AminoAcidTypes.GLU)
			return true;
		else
			return false;
	}

	/**
	 * The parameters are taken from the parm99 Amber force field
	 * The epsilons are sqrt epsilons to save tons of calls to sqrt
	 */
	public void assignLJParams()
	{
	  if(isN()) {
	    rmin2 = 1.824;
	    // eps = 0.412310563; // 0.17;
	    eps = 0.17;
	  }
	  if(isO()) {
	    rmin2 = 1.6612;
	    //eps = 0.458257569; //0.21;
	    eps = 0.21;
	  }
	  if(isC()) {
	    rmin2 = 1.908;
	    //eps = 0.293257566; // 0.086;
	    eps =  0.086;
	  }
	  if(isCA() || isCB()) {
	    rmin2 = 1.908;
	    //eps = 0.330756708; //0.1094;
	    eps = 0.1094;
	  }
	}
	
	public String toString2()
	{
		return this.pdbLine();
	}
	
	public String toString()
	{
		return this.pdbLine();
	}
	
	public String pdbLine()
	{
		String pdbLine = "ATOM  ";
		int numberOfSpaces;
		
		// Append Atom serial number substring to pdbLine
		// According to PDB standard it is a 4 character string
		// that is located between index 7 and index 10 of PDB entry
		int serialNumber = this.getIndex();
		numberOfSpaces = 0;
		if(serialNumber < 10)
			numberOfSpaces = 4;
		else if(serialNumber < 100)
			numberOfSpaces = 3;
		else if(serialNumber < 1000)
			numberOfSpaces = 2;
		else if(serialNumber < 10000)
			numberOfSpaces = 1;
		else
			numberOfSpaces = 0;
		
		String serialNumberStr = "";
		// Prepend required number of spaces to serialNumber substring of PDB entry
		for(int i = 0; i < numberOfSpaces; i++)
		{
			serialNumberStr = " " + serialNumberStr;
		}
		serialNumberStr = serialNumberStr + serialNumber;
		
		pdbLine = pdbLine + serialNumberStr;
		
		// According to PDB standard, index 11 of PDB entry is a space
		// so append a space to the PDB entry
		pdbLine = pdbLine + " ";
		
		// Append Atom name substring to pdbLine
		// According to PDB standard it is a 4 character string
		// that is located between index 12 and index 15 of PDB entry
		String atomName = Definitions.getAtomTypeString(this.getAtomType());
		pdbLine = pdbLine + atomName;
		
		// According to PDB standard, index 16 of PDB entry is a space
		// so append a space to the PDB entry
		pdbLine = pdbLine + " ";
		
		// Append Residue name substring to pdbLine
		// According to PDB standard it is a 4 character string
		// that is located between index 17 and index 19 of PDB entry
		String residueName = Definitions.getAminoAcidTypeString(this.getAminoAcidType());
		pdbLine = pdbLine + residueName;
		
		// According to PDB standard, index 20 of PDB entry is a space
		// so append a space to the PDB entry
		pdbLine = pdbLine + " ";
		
		// According to PDB standard, index 21 of PDB entry contains
		// chain ID, so append the chain ID to the PDB entry
		pdbLine = pdbLine + this.getChainID();
		
		// Append Residue index number substring to pdbLine
		// According to PDB standard it is a 4 character string
		// that is located between index 22 and index 25 of PDB entry
		int resIndex = this.getResidueIndex();
		numberOfSpaces = 0;
		if(resIndex < 10)
			numberOfSpaces = 3;
		else if(resIndex < 100)
			numberOfSpaces = 2;
		else if(resIndex < 1000)
			numberOfSpaces = 1;
		else
			numberOfSpaces = 0;
		
		String residueIndexStr = "";
		// Prepend required number of spaces to serialNumber substring of PDB entry
		for(int i = 0; i < numberOfSpaces; i++)
		{
			residueIndexStr = " " + residueIndexStr;
		}
		residueIndexStr = residueIndexStr + resIndex;
		
		pdbLine = pdbLine + residueIndexStr;
		
		// According to PDB standard, substring of PDB entry between 
		// index 26 and index 29 contains a string that can be ignored for us for now
		// so append 4 spaces to the PDB entry
		pdbLine = pdbLine + "    ";
		
		// Append Atom Coordinates to the PDB entry
		pdbLine = pdbLine + this.getCoordinatesString();
		
		return pdbLine;
	}
	
	public boolean equals(Object o){
		if(o instanceof PDBAtom)
		{
			PDBAtom a1 = (PDBAtom)o;
			if(this.atomType == a1.getAtomType() &&
					this.aminoAcidType == a1.getAminoAcidType() &&
					this.residueIndex == a1.getResidueIndex() &&
					this.chainID.equalsIgnoreCase(a1.getChainID()))
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			return false;
		}
		
	}

	public double pairwiseLJPotential(PDBAtom a2){
		// This function computes the Lennard-Jones potential between two atoms.
		double rMin = (this.getRmin2() + a2.getRmin2()) / 2;
		double eps = this.getEps() * a2.getEps();
		eps = Math.sqrt(eps);
		double dist = this.distance(a2);
		
		double ljTerm1 = Math.pow(rMin / dist, 12);
		double ljTerm2 = Math.pow(rMin / dist, 6);
		double vLJr = eps * (ljTerm1 - 2 * ljTerm2);
		//double ljTerm1 = Math.pow(rMin / dist, 9);
		//double ljTerm2 = Math.pow(rMin / dist, 6);
		//double vLJr = eps * (ljTerm1 - 2*ljTerm2);
		
		return vLJr;	
	}

	
	/////////////////////////////////////////
	//// The attributes of PDBAtom class ////
	/////////////////////////////////////////
	
	
	/**
	 * The type of the atom
	 */
	private AtomTypes atomType;
	
	/**
	 * The type of the amino acid that contains the atom
	 */
	private AminoAcidTypes aminoAcidType;
	
	/**
	 * The chain ID of the atom
	 */
	private String chainID;
	
	/**
	 * vdW parameter - well depth
	 */
	private double eps;
	
	/**
	 * vdW parameter - radius
	 */
	private double rmin2;
	
	/**
	 * Residue index
	 */
	private int residueIndex;
	
	private TerminalTypes terminalType;
}
