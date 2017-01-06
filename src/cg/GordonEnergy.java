package cg;

import pdb.PDBAtom;
import utilities.Vector;
import utilities.Definitions.AminoAcidTypes;

/**
 * @author Dong Luo 2012-08-27
 *
 */

public class GordonEnergy {
	/**
	 * Reference:
	 * A coarse-grained α-carbon protein model with anisotropic hydrogen-bonding
	 * Yap EH, Fawzi NL, Head-Gordon T.
	 * Proteins. 2008 Feb 15;70(3):626-38.
	 */
	
	/**
	 * The atom types
	 */
	public static enum AtomTypes{
		// B is for strong attraction
		// V is for weak attraction
		// N is for weak repulsion
		// L is for strong repulsion
		B, V, N, L 
	}
	/**
	 * The dihedral angle types
	 */
	public static enum DihedralTypes{
		// H is for helical
		// E is for extended
		// T is for turn other than P, U, Q
		// P is for +90° turn
		// U is for 0° turn
		// Q is for -90° turn
		H, E, T, P, U, Q 
	}
	/**
	 * The hydrogen bond types
	 */
	public static enum HydrogenBondCapability{
		// A is for helical dihedral
		// B is for extended dihedral
		// C is for other dihedral
		A, B, C 
	}
	
	public static double defaultEnergyScale=1.0;
	public static double helicalOptimalAngle=95;
	public static double optimalAngle=105;
	public static double angleConstant=20, angleEnergyThreshold=15;
	public static double minAngle=50, maxAngle=150, minDihedralAngle=-180, maxDihedralAngle=180;
	public static double HA=0, HB=1.2, HC=1.2, HD=1.2, HAngle=0.17;
	public static double EA=0.45, EB=0, EC=0.6, ED=0, EAngle=-0.35;
	public static double TA=0.2, TB=0.2, TC=0.2, TD=0.2, TAngle=0;
	public static double PA=0.36, PB=0, PC=0.48, PD=0, PAngle=1.57;
	public static double QA=0.36, QB=0, QC=0.48, QD=0, QAngle=-1.57;
	public static double UA=0.36, UB=0, UC=0.48, UD=0, UAngle=3.14;
	public static double BBS1=1.4, BBS2=1.0, BVS1=0.7, BVS2=1.0, VVS1=0.35, VVS2=1.0;
	public static double LLS1=0.333, LLS2=-1.0, LVS1=0.333, LVS2=-1.0, LBS1=0.333, LBS2=-1.0;
	public static double NXS1=1.0, NXS2=0.0;
	public static double vanDeWallsSum=1.16, vanDeWallsEnergyThreshold=3;
	public static double HBlowScale=0.7, HBHighScale=0.98;
	public static double lengthUnit=3.8, HBAlphaDistance=1.35*lengthUnit, HBSheetDistance=1.25*lengthUnit;
	public static double HBSheetSearchDistance=3.0*lengthUnit;
	public static double HBWidth=0.5;
	public static double sampleRatio=0.05, sampleAngleRange=10;
	public static double weightAngle=1, weightDistance=5;
	
	public static AtomTypes getAtomType(AminoAcidTypes aminoType)
	{
		//B-strong attraction
		if(aminoType == AminoAcidTypes.TRP) return AtomTypes.B;
		else if(aminoType == AminoAcidTypes.CYS) return AtomTypes.B;
		else if(aminoType == AminoAcidTypes.LEU) return AtomTypes.B;
		else if(aminoType == AminoAcidTypes.ILE) return AtomTypes.B;
		else if(aminoType == AminoAcidTypes.PHE) return AtomTypes.B;
		else if(aminoType == AminoAcidTypes.MET) return AtomTypes.B;
		else if(aminoType == AminoAcidTypes.VAL) return AtomTypes.B;
		// V-weak attraction
		else if(aminoType == AminoAcidTypes.ALA) return AtomTypes.V;
		else if(aminoType == AminoAcidTypes.TYR) return AtomTypes.V;
		// L-strong repulsion
		else if(aminoType == AminoAcidTypes.GLU) return AtomTypes.L;
		else if(aminoType == AminoAcidTypes.ASP) return AtomTypes.L;
		else if(aminoType == AminoAcidTypes.ASN) return AtomTypes.L;
		else if(aminoType == AminoAcidTypes.HIS) return AtomTypes.L;
		else if(aminoType == AminoAcidTypes.GLN) return AtomTypes.L;
		else if(aminoType == AminoAcidTypes.LYS) return AtomTypes.L;
		else if(aminoType == AminoAcidTypes.ARG) return AtomTypes.L;
		// N-weak repulsion
		else if(aminoType == AminoAcidTypes.PRO) return AtomTypes.N;
		else if(aminoType == AminoAcidTypes.GLY) return AtomTypes.N;
		else if(aminoType == AminoAcidTypes.SER) return AtomTypes.N;
		else if(aminoType == AminoAcidTypes.THR) return AtomTypes.N;
		else return AtomTypes.N;
	}

	public static DihedralTypes getDihedralType(char c) {
		if (c=='H') return DihedralTypes.H;
		if (c=='E') return DihedralTypes.E;
		if (c=='T') return DihedralTypes.T;
		if (c=='P') return DihedralTypes.P;
		if (c=='Q') return DihedralTypes.Q;
		if (c=='U') return DihedralTypes.U;
		System.out.println("Dihedral type "+c+" not recognized, will use the default type T.");
		return DihedralTypes.T;
	}
	
	/**
	 * 
	 * @param angle is the dihedral angle in radians
	 * @return corresponding dihedral type, only global minimum considered
	 */
	public static DihedralTypes getDihedralType(double angle) {
		double H=50, E=-160, T=180, P=90, Q=-90, U=0;
		angle=angle*180/java.lang.Math.PI;
		DihedralTypes ret=DihedralTypes.H;
		double min=java.lang.Math.abs(angle-H);
		double delta=java.lang.Math.abs(angle-E);
		if (min>delta) {
			min=delta;
			ret=DihedralTypes.E;
		}
		delta=java.lang.Math.abs(angle-T);
		if (min>delta) {
			min=delta;
			ret=DihedralTypes.T;
		}
		delta=java.lang.Math.abs(angle-P);
		if (min>delta) {
			min=delta;
			ret=DihedralTypes.P;
		}
		delta=java.lang.Math.abs(angle-Q);
		if (min>delta) {
			min=delta;
			ret=DihedralTypes.Q;
		}
		delta=java.lang.Math.abs(angle-U);
		if (min>delta) {
			min=delta;
			ret=DihedralTypes.U;
		}
		return ret;
	}

	/**
	 * 
	 * @param v is the value
	 * @param r1 is the first reference
	 * @param r2 is the second reference
	 * @param r3 is the third reference
	 * @return the minimum delta between value and references
	 */
	private static double minDelta(double v, double r1, double r2, double r3) {
		double min=java.lang.Math.abs(v-r1);
		double delta=java.lang.Math.abs(v-r2);
		if (min>delta) min=delta;
		delta=java.lang.Math.abs(v-r3);
		if (min>delta) min=delta;
		return min;
	}
	
	/**
	 * 
	 * @param angle is the diheral angle in radians
	 * @return the corresponding dihedral type, both global and local minimums considered
	 */
	public static DihedralTypes getDihedralType2(double angle) {
		double H1=-65, H2=50, H3=165;
		double E1=-160, E2=-45, E3=85;
		double T1=-60, T2=60, T3=180;
		double P1=-155, P2=-25, P3=90;
		double Q1=-90, Q2=25, Q3=155;
		double U1=-115, U2=0, U3=115;
		angle=angle*180/java.lang.Math.PI;
		DihedralTypes ret=DihedralTypes.H;
		double min=minDelta(angle, H1, H2, H3);
		double delta=minDelta(angle, E1, E2, E3);
		if (min>delta) {
			min=delta;
			ret=DihedralTypes.E;
		}
		delta=minDelta(angle, T1, T2, T3);
		if (min>delta) {
			min=delta;
			ret=DihedralTypes.T;
		}
		delta=minDelta(angle, P1, P2, P3);
		if (min>delta) {
			min=delta;
			ret=DihedralTypes.P;
		}
		delta=minDelta(angle, Q1, Q2, Q3);
		if (min>delta) {
			min=delta;
			ret=DihedralTypes.Q;
		}
		delta=minDelta(angle, U1, U2, U3);
		if (min>delta) {
			min=delta;
			ret=DihedralTypes.U;
		}
		return ret;
	}

	public static HydrogenBondCapability getHydrogenBondCapability(DihedralTypes type) {
		if (type == DihedralTypes.H) return HydrogenBondCapability.A;
		if (type == DihedralTypes.E) return HydrogenBondCapability.B;
		return HydrogenBondCapability.C;
	}
	
	/**
	 * 
	 * @param angle is in radians
	 * @param type is dihedral propensity of previous atom
	 * @return angle energy located at current atom
	 */
	public static double angleEnergy(double angle, DihedralTypes type) {
		double detaAngle;
		if (type == DihedralTypes.H)
			detaAngle=angle-helicalOptimalAngle*java.lang.Math.PI/180.0;
		else
			detaAngle=angle-optimalAngle*java.lang.Math.PI/180.0;
		return 0.5*angleConstant*detaAngle*detaAngle;
	}
	
	/**
	 * 
	 * @param angle is in radians
	 * @param type is the dihedral angle type
	 * @return dihedral angle energy
	 */
	public static double dihedralEnergy(double angle, DihedralTypes type) {
		double A, B, C, D, tAngle;
		if (type == DihedralTypes.H) {
			A=HA; B=HB; C=HC; D=HD; tAngle=angle+HAngle;
		} else if (type == DihedralTypes.E) {
			A=EA; B=EB; C=EC; D=ED; tAngle=angle+EAngle;
		} else if (type == DihedralTypes.P) {
			A=PA; B=PB; C=PC; D=PD; tAngle=angle+PAngle;
		} else if (type == DihedralTypes.Q) {
			A=QA; B=QB; C=QC; D=QD; tAngle=angle+QAngle;
		} else if (type == DihedralTypes.U) {
			A=UA; B=UB; C=UC; D=UD; tAngle=angle+UAngle;
		} else {
			A=TA; B=TB; C=TC; D=TD; tAngle=angle+TAngle;
		}
		return (A*(1+java.lang.Math.cos(tAngle))+
				B*(1-java.lang.Math.cos(tAngle))+
				C*(1+java.lang.Math.cos(3.0*tAngle))+
				D*(1+java.lang.Math.cos(tAngle+java.lang.Math.PI/4.0))
				);
	}
	
	/**
	 * 
	 * @param r is distance between two atoms
	 * @param t1 is CG atom type of first atom
	 * @param t2 is CG atom type of second atom
	 * @return van de Walls energy between the two atoms
	 */
	public static double vanDeWallsEnergy(double r, AtomTypes t1, AtomTypes t2) {
		double S1, S2;
		if (t1 == AtomTypes.B && t2 == AtomTypes.B) {
			S1=BBS1; S2=BBS2;
		} else if ((t1 == AtomTypes.B && t2 == AtomTypes.V) ||
				(t1 == AtomTypes.V && t2 == AtomTypes.B)) {
			S1=BVS1; S2=BVS2;
		} else if (t1 == AtomTypes.V && t2 == AtomTypes.V) {
			S1=VVS1; S2=VVS2;
		} else if (t1 == AtomTypes.L && t2 == AtomTypes.L) {
			S1=LLS1; S2=LLS2;
		} else if ((t1 == AtomTypes.L && t2 == AtomTypes.V) ||
				(t1 == AtomTypes.V && t2 == AtomTypes.L)) {
			S1=LVS1; S2=LVS2;
		} else if ((t1 == AtomTypes.L && t2 == AtomTypes.B) ||
				(t1 == AtomTypes.B && t2 == AtomTypes.L)) {
			S1=LBS1; S2=LBS2;
		} else {
			S1=NXS1; S2=NXS2;
		}
		return 4.0*S1*(java.lang.Math.pow(vanDeWallsSum/r, 12)-S2*java.lang.Math.pow(vanDeWallsSum/r, 6));
	}
	
	/**
	 * 
	 * @param a0, a1, a2 are the i-1, i, i+1 atoms
	 * @param a3, a4, a5 are the j-1, j, j+1 atoms
	 * @param type describes where the hydrogen bond locates, can be A(helical) or (B)sheet
	 * @return hydrogen bond energy for atom pair i, j
	 * rij is the vector between atom i and j
	 * tHBi is the unit vector normal to plane defined by atoms i-1, i, i+1
	 * tHBj is the unit vector normal to plane defined by atoms j-1, j, j+1
	 */
	public static double HBEnergy(PDBAtom a0, PDBAtom a1, PDBAtom a2,
					PDBAtom a3, PDBAtom a4, PDBAtom a5, HydrogenBondCapability type) {
		Vector rij = a4.subtract(a1);
		Vector tHBi = a1.subtract(a0).crossProduct(a2.subtract(a1));
		tHBi = tHBi.multiply(1.0/tHBi.magnitude());
		Vector tHBj = a4.subtract(a3).crossProduct(a5.subtract(a4));
		tHBj = tHBj.multiply(1.0/tHBj.magnitude());

		double r=rij.magnitude();
		Vector rijUnit=rij.multiply(1/r);
		if (type == HydrogenBondCapability.A) r-=HBAlphaDistance;
		else if (type == HydrogenBondCapability.B) r-=HBSheetDistance;
		else return 0;
		double HBdist=HBWidth*HBWidth;
		double F=java.lang.Math.exp(-r*r/HBdist);
		double G=java.lang.Math.exp((java.lang.Math.abs(tHBi.dotProduct(rijUnit))-1)/HBdist);
		double H=java.lang.Math.exp((java.lang.Math.abs(tHBj.dotProduct(rijUnit))-1)/HBdist);
		return F*G*H;
	}
	
	public static String angleType="Angle";
	public static String dihedralType="Dihedral";
	public static String vanDeWallsType="vanDeWalls";
	public static String hydrogenBondType="HydrogenBond";

}
