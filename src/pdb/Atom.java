package pdb;

import utilities.Vector;
import javax.vecmath.GMatrix;
import java.lang.Math;
import java.text.DecimalFormat;

/**
 * @author Bahar Akbal-Delibas (abakbal@cs.umb.edu)
 * 
 */

public abstract class Atom extends Vector{
	
	////////////////////////////////////
	//// Constructors of Atom class ////
	////////////////////////////////////
	
	/**
	 * Default constructor of Atom class
	 * Sets energy to 0.0 and index to 0
	 */
	public Atom(){
		super();
		this.energy = 0.0;
		this.index = 0;
	}
	
	/**
	 * Constructor of Atom class.
	 * Sets 
	 * 	 the energy of the atom to the given energy value, 
	 * 	 the index of the atom to the given index value, 
	 * 	 3D coordinates of the atom to the given double values.
	 * @param theIndex is a integer value that is set as the index of the created atom
	 * @param theX position of the atom on the x-coordinate.
	 * @param theY position of the atom on the y-coordinate.
	 * @param theZ position of the atom on the z-coordinate.
	 */
	public Atom(int theIndex, double theX, double theY, double theZ){
		super(theX, theY, theZ);
		this.energy = 0.0;
		this.index = theIndex;
	}
	
	///////////////////////////////////////////
	//// Setters and getters of Atom class ////
	///////////////////////////////////////////
	
	/**
	 * The setter of energy attribute
	 * @param theEnergy is a double value that is set as the energy of this atom
	 */
	public void setEnergy(double theEnergy){
		this.energy = theEnergy;
	}
	
	/**
	 * The setter of index attribute
	 * @param theIndex is an integer value that is set as the index of this atom
	 */
	public void setIndex(int theIndex){
		this.index = theIndex;
	}
	
	/**
	 * The setter of coordinates attribute
	 * @param theCoordinates is an instance of utilities.Vector class that is 
	 * set as the 3D coordinates of this atom.
	 */
	/*public void setCoordinates(Vector theCoordinateVector){
		this.coordinateVector = theCoordinateVector;
	}*/
	
	/**
	 * The getter of energy attribute
	 * @return the energy of this atom.
	 */
	public double getEnergy(){
		return this.energy;
	}
	
	/**
	 * The getter of index attribute
	 * @return the index of this atom
	 */
	public int getIndex(){
		return this.index;
	}
	
	/**
	 * The getter of coordinates attribute
	 * @return the 3D coordinates of this atom
	 */
	/*public Vector getCoordinateVector(){
		return this.coordinateVector;
	}*/
	
	/////////////////////////////////////
	//// Other methods of Atom class ////
	/////////////////////////////////////
	
	/**
	 * The distance between this atom and the atom p
	 * @param p The atom whose distance from this atom is calculated by this method.
	 * @return The distance between this atom and the atom @param p
	 */
	/*public double distanceFrom(Atom p){
		// get the distance between coordinate vectors of these two atoms 
		double dist = this.coordinateVector.distance(p.getCoordinateVector());
		// Calculate the distance between this atom and the atom p
		return dist;
	}*/
	
	/**
	 * The method isNeighborWith checks if the given atom p
	 * is closer than a given threshold, maxDistance, to this atom.
	 * @param p The atom that is checked if it is a neighbor of this atom.
	 * @param maxDistance 
	 * @return true if p is a neighbor of this atom; false otherwise.
	 */
	public boolean isNeighborWith(Atom p, double maxDistance)
	{
		double d = this.distance(p);
		return (d <= maxDistance * 1.01);
	}
	
	/**
	 * This method translates the atom on the x, y and z coordinates 
	 * as far as the given delta values.
	 * @param deltaX the distance value that this atom will be translated on x-coordinate
	 * @param deltaY the distance value that this atom will be translated on y-coordinate
	 * @param deltaZ the distance value that this atom will be translated on z-coordinate
	 */
	public void translate(double deltaX, double deltaY, double deltaZ)
	{
		this.setX(this.getX() + deltaX);
		this.setY(this.getY() + deltaY);
		this.setZ(this.getZ() + deltaZ);
		/*
		GMatrix translationMatrix = new GMatrix(4, 4);
		translationMatrix.setColumn(0, new double[]{1, 0, 0, 0});
		translationMatrix.setColumn(1, new double[]{0, 1, 0, 0});
		translationMatrix.setColumn(2, new double[]{0, 0, 1, 0});
		translationMatrix.setColumn(3, new double[]{deltaX, deltaY, deltaZ, 1});
		
		GMatrix atomPositionVector = new GMatrix(4, 1);
		atomPositionVector.setColumn(0, new double[]{this.x, this.y, this.z, 1});
		
		GMatrix resultingAtomPositionVector = new GMatrix(4, 1);
		resultingAtomPositionVector.mul(translationMatrix, atomPositionVector);
		
		// Set the new position of this atom after the rotation
		double newX = resultingAtomPositionVector.getElement(0, 0);
		double newY = resultingAtomPositionVector.getElement(1, 0);
		double newZ = resultingAtomPositionVector.getElement(2, 0);
		
		this.setX(newX);
		this.setY(newY);
		this.setZ(newZ);
		*/
	}
	
	/**
	 * This method performs rotation of an atom around an arbitrary bond,
	 * which is specified by two other atoms at the beginning and ending edges. 
	 * This method uses axis-angle representation to perform rotation.
	 * @param begin the starting edge of the bond around which this atom is rotated
	 * @param end the ending edge of the bond around which this atom is rotated
	 * @param theta rotation angle
	 * */
	public void rotate(Vector begin, Vector end, double theta)
	{
		// In order to find the resulting position vector of this atom after 
		// rotation, represented by transpose( [x', y', z', 1] ), we should multiply 
		// the rotation matrix of the bond R(bond, theta) by the position vector of 
		// this atom, transpose( [x, y, z, 1]. 
		//
		// transpose( [x', y', z', 1] ) = R(bond, theta) * transpose( [x, y, z, 1] )
		//
		// R(bond, theta) = T(begin) * R(axis, theta) * T(-begin)
		// where 
		// 		T(x) is a translation by the vector x
		//		R(axis, theta) is a rotation around an axis that goes through
		//						the origin of the specified coordinate system.
		//
		// This rotation around an arbitrary point is realized by translating the point 
		// to the origin, rotating the target atom around the axis through the origin, 
		// and then translating it back (the composition of these 3 transformations yields 
		// a unique 4x4 homogeneous matrix that achieves the same effect).
		//
		// The axis of rotation can easily be computed from the positions of atoms begin 
		// and end, and must have unit norm.
		
		// IMPORTANT NOTE: java.Math class expects radian values as parameters for 
		// trigonometric functions. 
		// Thus, we need to convert the theta value in degrees to radian.
		theta = Math.toRadians(theta);
		
		// First, create the forward translation matrix, T(begin)
		GMatrix forwardTranslationMatrix = new GMatrix(4, 4);
		forwardTranslationMatrix.setColumn(0, new double[]{1, 0, 0, 0});
		forwardTranslationMatrix.setColumn(1, new double[]{0, 1, 0, 0});
		forwardTranslationMatrix.setColumn(2, new double[]{0, 0, 1, 0});
		forwardTranslationMatrix.setColumn(3, new double[]{begin.getX(), begin.getY(), begin.getZ(), 1});
		
		// Then create the backward translation matrix, T(-begin) 
		GMatrix backwardTranslationMatrix = new GMatrix(4, 4);
		backwardTranslationMatrix.setColumn(0, new double[]{1, 0, 0, 0});
		backwardTranslationMatrix.setColumn(1, new double[]{0, 1, 0, 0});
		backwardTranslationMatrix.setColumn(2, new double[]{0, 0, 1, 0});
		backwardTranslationMatrix.setColumn(3, new double[]{-1*begin.getX(), -1*begin.getY(), -1*begin.getZ(), 1});
		
		// Create the axis of rotation by subtracting the coordinates 
		// of atoms begin and end, and then dividing by its norm.
		double axisX = end.getX() - begin.getX();
		double axisY = end.getY() - begin.getY();
		double axisZ = end.getZ() - begin.getZ();
		double axisNorm = Math.sqrt(axisX * axisX + axisY * axisY + axisZ * axisZ);		
		Vector axis = new Vector(axisX/axisNorm, axisY/axisNorm, axisZ/axisNorm);
		
		// Create the axis rotation matrix, R(axis, theta)
		// Rotation Matrix is created as described by John J. Craig in
		// "Intoduction to Robotics Mechanics & Control" page: 47, Equation 2.77
		GMatrix axisRotationMatrix = new GMatrix(4, 4);
		axisRotationMatrix.setColumn(0, 
				new double[]{axis.getX() * axis.getX() * (1-Math.cos(theta)) + Math.cos(theta),
							 axis.getX() * axis.getY() * (1-Math.cos(theta)) - axis.getZ() * Math.sin(theta),
							 axis.getX() * axis.getZ() * (1-Math.cos(theta)) + axis.getY() * Math.sin(theta),
							 0});
		
		axisRotationMatrix.setColumn(1, 
				new double[]{axis.getX() * axis.getY() * (1-Math.cos(theta)) + axis.getZ() * Math.sin(theta),
							 axis.getY() * axis.getY() * (1-Math.cos(theta)) + Math.cos(theta),
							 axis.getY() * axis.getZ() * (1-Math.cos(theta)) - axis.getX() * Math.sin(theta),
							 0});
		
		axisRotationMatrix.setColumn(2, 
				new double[]{axis.getX() * axis.getZ() * (1-Math.cos(theta)) - axis.getY() * Math.sin(theta),
							 axis.getY() * axis.getZ() * (1-Math.cos(theta)) + axis.getX() * Math.sin(theta),
							 axis.getZ() * axis.getZ() * (1-Math.cos(theta)) + Math.cos(theta),
							 0});
		
		axisRotationMatrix.setColumn(3, new double[]{0, 0, 0, 1});
		
		// Now create the bond rotation matrix, R(bond, theta)
		GMatrix bondRotationMatrix = (GMatrix) forwardTranslationMatrix.clone();
		bondRotationMatrix.mul(axisRotationMatrix);
		bondRotationMatrix.mul(backwardTranslationMatrix);
		
		// Finally, multiply the position vector of this atom, 
		// transpose( [x, y, z, 1] ), by the bond rotation matrix to 
		// calculate the resulting position vector after the rotation.
		GMatrix atomPositionVector = new GMatrix(4, 1);
		atomPositionVector.setColumn(0, new double[]{this.x, this.y, this.z, 1});
		
		// GMatrix resultingAtomPositionVector = (GMatrix) bondRotationMatrix.clone();
		GMatrix resultingAtomPositionVector = new GMatrix(4, 1);
		resultingAtomPositionVector.mul(bondRotationMatrix, atomPositionVector);
		
		// Set the new position of this atom after the rotation
		double newX = resultingAtomPositionVector.getElement(0, 0);
		double newY = resultingAtomPositionVector.getElement(1, 0);
		double newZ = resultingAtomPositionVector.getElement(2, 0);
		
		this.setX(newX);
		this.setY(newY);
		this.setZ(newZ);
	}
	
	public void rotate(double alpha, double beta, double theta)
	{
		// In order to find the resulting position vector of this atom after 
		// rotation, represented by transpose( [x', y', z', 1] ), we should multiply 
		// the rotation matrix of the bond R(bond, theta) by the position vector of 
		// this atom, transpose( [x, y, z, 1]. 
		//
		// transpose( [x', y', z', 1] ) = R(axis, theta) * transpose( [x, y, z, 1] )
		//
		// where 
		//		R(axis, theta) is a rotation around an axis that goes through
		//						the origin of the specified coordinate system and makes
		//						alpha degrees angle with x-axis, 
		//						90-alpha degrees angle with y-axis
		//						beta degrees angle with z-axis.
		//
		// The axis of rotation is computed using alpha and beta angles 
		// and must have unit norm.
		
		// IMPORTANT NOTE: java.Math class expects radian values as parameters for 
		// trigonometric functions. 
		// Thus, we need to convert all angle values in degrees to radian.
		alpha = Math.toRadians(alpha);
		beta = Math.toRadians(beta);
		theta = Math.toRadians(theta);
		double pi = Math.toRadians(180);
		
		// Calculate the rotation of axis using alpha and beta angles
		double axisX = Math.cos(alpha);
		double axisY = Math.cos(pi/2 - alpha);
		double axisZ = Math.cos(beta);
		Vector axis = new Vector(axisX, axisY, axisZ);
		
		// Create the axis rotation matrix, R(axis, theta)
		// Rotation Matrix is created as described by John J. Craig in
		// "Intoduction to Robotics Mechanics & Control" page: 47, Equation 2.77
		GMatrix axisRotationMatrix = new GMatrix(4, 4);
		axisRotationMatrix.setColumn(0, 
				new double[]{axis.getX() * axis.getX() * (1-Math.cos(theta)) + Math.cos(theta),
							 axis.getX() * axis.getY() * (1-Math.cos(theta)) - axis.getZ() * Math.sin(theta),
							 axis.getX() * axis.getZ() * (1-Math.cos(theta)) + axis.getY() * Math.sin(theta),
							 0});
		
		axisRotationMatrix.setColumn(1, 
				new double[]{axis.getX() * axis.getY() * (1-Math.cos(theta)) + axis.getZ() * Math.sin(theta),
							 axis.getY() * axis.getY() * (1-Math.cos(theta)) + Math.cos(theta),
							 axis.getY() * axis.getZ() * (1-Math.cos(theta)) - axis.getX() * Math.sin(theta),
							 0});
		
		axisRotationMatrix.setColumn(2, 
				new double[]{axis.getX() * axis.getZ() * (1-Math.cos(theta)) - axis.getY() * Math.sin(theta),
							 axis.getY() * axis.getZ() * (1-Math.cos(theta)) + axis.getX() * Math.sin(theta),
							 axis.getZ() * axis.getZ() * (1-Math.cos(theta)) + Math.cos(theta),
							 0});
		
		axisRotationMatrix.setColumn(3, new double[]{0, 0, 0, 1});
		
		// Finally, multiply the position vector of this atom, 
		// transpose( [x, y, z, 1] ), by the axis rotation matrix to 
		// calculate the resulting position vector after the rotation.
		GMatrix atomPositionVector = new GMatrix(4, 1);
		atomPositionVector.setColumn(0, new double[]{this.x, this.y, this.z, 1});
		
		// GMatrix resultingAtomPositionVector = (GMatrix) bondRotationMatrix.clone();
		GMatrix resultingAtomPositionVector = new GMatrix(4, 1);
		resultingAtomPositionVector.mul(axisRotationMatrix, atomPositionVector);
		
		// Set the new position of this atom after the rotation
		double newX = resultingAtomPositionVector.getElement(0, 0);
		double newY = resultingAtomPositionVector.getElement(1, 0);
		double newZ = resultingAtomPositionVector.getElement(2, 0);
		
		this.setX(newX);
		this.setY(newY);
		this.setZ(newZ);
	}
	
	public void rotateAroundXAxis(double theta)
	{
		theta = Math.toRadians(theta);
		
		// Create the Euler rotation matrix for x axis.
		// Rotation Matrix is created as described by John J. Craig in
		// "Introduction to Robotics Mechanics & Control" page: 47, Equation 2.74
		GMatrix axisRotationMatrix = new GMatrix(4, 4);
		axisRotationMatrix.setColumn(0, new double[]{1, 0, 0, 0});
		axisRotationMatrix.setColumn(1, new double[]{0, Math.cos(theta), Math.sin(theta), 0});
		axisRotationMatrix.setColumn(2, new double[]{0, -1 * Math.sin(theta), Math.cos(theta), 0});
		axisRotationMatrix.setColumn(3, new double[]{0, 0, 0, 1});
		
		// Finally, multiply the position vector of this atom, 
		// transpose( [x, y, z, 1] ), by the bond rotation matrix to 
		// calculate the resulting position vector after the rotation.
		GMatrix atomPositionVector = new GMatrix(4, 1);
		atomPositionVector.setColumn(0, new double[]{this.x, this.y, this.z, 1});
		
		GMatrix resultingAtomPositionVector = new GMatrix(4, 1);
		resultingAtomPositionVector.mul(axisRotationMatrix, atomPositionVector);
		
		// Set the new position of this atom after the rotation
		double newX = resultingAtomPositionVector.getElement(0, 0);
		double newY = resultingAtomPositionVector.getElement(1, 0);
		double newZ = resultingAtomPositionVector.getElement(2, 0);
		
		this.setX(newX);
		this.setY(newY);
		this.setZ(newZ);
		
	}
	
	public void rotateAroundYAxis(double theta)
	{
		theta = Math.toRadians(theta);
		
		// Create the Euler rotation matrix for y axis.
		// Rotation Matrix is created as described by John J. Craig in
		// "Introduction to Robotics Mechanics & Control" page: 47, Equation 2.75
		GMatrix axisRotationMatrix = new GMatrix(4, 4);
		axisRotationMatrix.setColumn(0, new double[]{Math.cos(theta), 0, -1 * Math.sin(theta), 0});
		axisRotationMatrix.setColumn(1, new double[]{0, 1, 0, 0});
		axisRotationMatrix.setColumn(2, new double[]{Math.sin(theta), 0, Math.cos(theta), 0});
		axisRotationMatrix.setColumn(3, new double[]{0, 0, 0, 1});
		
		// Finally, multiply the position vector of this atom, 
		// transpose( [x, y, z, 1] ), by the bond rotation matrix to 
		// calculate the resulting position vector after the rotation.
		GMatrix atomPositionVector = new GMatrix(4, 1);
		atomPositionVector.setColumn(0, new double[]{this.x, this.y, this.z, 1});
		
		GMatrix resultingAtomPositionVector = new GMatrix(4, 1);
		resultingAtomPositionVector.mul(axisRotationMatrix, atomPositionVector);
		
		// Set the new position of this atom after the rotation
		double newX = resultingAtomPositionVector.getElement(0, 0);
		double newY = resultingAtomPositionVector.getElement(1, 0);
		double newZ = resultingAtomPositionVector.getElement(2, 0);
		
		this.setX(newX);
		this.setY(newY);
		this.setZ(newZ);
		
	}
	
	public void rotateAroundZAxis(double theta)
	{
		theta = Math.toRadians(theta);
		
		// Create the Euler rotation matrix for z axis.
		// Rotation Matrix is created as described by John J. Craig in
		// "Introduction to Robotics Mechanics & Control" page: 47, Equation 2.76
		GMatrix axisRotationMatrix = new GMatrix(4, 4);
		axisRotationMatrix.setColumn(0, new double[]{Math.cos(theta), Math.sin(theta), 0, 0});
		axisRotationMatrix.setColumn(1, new double[]{-1 * Math.sin(theta), Math.cos(theta), 0, 0});
		axisRotationMatrix.setColumn(2, new double[]{0, 0, 1, 0});
		axisRotationMatrix.setColumn(3, new double[]{0, 0, 0, 1});
		
		// Finally, multiply the position vector of this atom, 
		// transpose( [x, y, z, 1] ), by the bond rotation matrix to 
		// calculate the resulting position vector after the rotation.
		GMatrix atomPositionVector = new GMatrix(4, 1);
		atomPositionVector.setColumn(0, new double[]{this.x, this.y, this.z, 1});
		
		GMatrix resultingAtomPositionVector = new GMatrix(4, 1);
		resultingAtomPositionVector.mul(axisRotationMatrix, atomPositionVector);
		
		// Set the new position of this atom after the rotation
		double newX = resultingAtomPositionVector.getElement(0, 0);
		double newY = resultingAtomPositionVector.getElement(1, 0);
		double newZ = resultingAtomPositionVector.getElement(2, 0);
		
		this.setX(newX);
		this.setY(newY);
		this.setZ(newZ);	
	}
	
	
	/*
	 * This method rotates atom by using homogeneous coordinate 
	 * transformation matrix
	 * @param theta rotation around x-axis
	 * @param phi rotation around y-axis
	 * @param psi rotation around z-axis
	 * */
	/*
	 public void rotate(double theta, double phi, double psi)
	{
		// First, create this atom's position vector.
		// Since we are applying a homogeneous transformation, and
		// homogeneous coordinate transformation matrices operate on 
		// 4D homogeneous coordinate vector representation,
		// we add a fourth dimension as 1 to the coordinate vector.
		GMatrix atomPositionVector = new GMatrix(4, 1);
		atomPositionVector.setColumn(0, new double[]{this.x, this.y, this.z, 1});
		
		// Create the 4X4 homogeneous transformation matrix for rotation
		GMatrix transformationMatrix = new GMatrix(4, 4);
		transformationMatrix.setColumn(0, 
				new double[]{Math.cos(psi) * Math.cos(phi) + Math.sin(psi) * Math.sin(theta) * Math.sin(phi),
							-1 * Math.sin(psi) * Math.cos(theta),
							Math.sin(psi) * Math.sin(theta) * Math.cos(phi) - Math.cos(psi) * Math.sin(phi),
							0});
		
		transformationMatrix.setColumn(1, 
				new double[]{Math.sin(psi) * Math.cos(phi) - Math.cos(psi) * Math.sin(theta) *Math.sin(phi),
							 Math.cos(psi) * Math.cos(theta),
							 -1 * Math.cos(psi) * Math.sin(theta) * Math.cos(phi) - Math.sin(psi) * Math.sin(phi),
							 0});
		
		transformationMatrix.setColumn(2, 
				new double[]{Math.cos(theta) * Math.sin(phi),
							 Math.sin(theta),
							 Math.cos(theta) * Math.cos(phi),
							 0});
		
		transformationMatrix.setColumn(3, 
				new double[]{0, 0, 0, 1});
		
		// Now, multiply the transformation matrix with the coordinate vector
		// to obtain the new coordinate vector.
		// GMatrix.mul(GMatrix m1) sets the value of a matrix m to the result of 
		// multiplying itself with matrix m1 (m = m * m1).
		transformationMatrix.mul(atomPositionVector);
		
		// Finally, set new x, y, z values of this atom after transformation.
		double[] newCoordinates = new double[4];
		transformationMatrix.getColumn(0, newCoordinates);
		this.x = newCoordinates[0];
		this.y = newCoordinates[1];
		this.z = newCoordinates[2];
	}
	
	*/
	public String getCoordinatesString()
	{
		return normalizeCoordinateString(this.x) + normalizeCoordinateString(this.y) + normalizeCoordinateString(this.z); 
	}
	
	private static String normalizeCoordinateString(double coord)
	{
		// BigDecimal bd = new BigDecimal(coord);
	    // bd = bd.setScale(3);
	    // double roundedCoord = bd.doubleValue();
		DecimalFormat roundFormatter = new DecimalFormat("###0.000");
		double roundedCoord = Double.parseDouble(roundFormatter.format(coord));
	    
	    int numberOfSpacesToPrepend;
	    if(Math.abs(roundedCoord) < 10)
	    	numberOfSpacesToPrepend = 3;
	    else if(Math.abs(roundedCoord) < 100)
	    	numberOfSpacesToPrepend = 2;
	    else if(Math.abs(roundedCoord) < 1000)
	    	numberOfSpacesToPrepend = 1;
	    else
	    	numberOfSpacesToPrepend = 0;
	    
	    // If the rounded coordinate value is a negative number
	    // then need to prepend "-" to the number
	    // so decrement number of spaces to prepend by 1
	    if(roundedCoord < 0)
	    	numberOfSpacesToPrepend = numberOfSpacesToPrepend - 1;
	    
	    // When rounded Coordinate is zero but the actual coordinate is 
	    // a very small negative number but not zero (e.g, -0.0000001)
	    // roundFormatter.format() returns "-0.000", which becomes a problem in the PDB file.
	    // In order to avoid this problem, make sure that normalizedCoordinateString is always
	    // set to "0.000" when rounded Coordinate is zero.
	    String normalizedCoordinateString = "";
	    if(roundedCoord == 0)
	    {
	    	normalizedCoordinateString = "0.000";
	    }
	    else
	    {
	    	normalizedCoordinateString = roundFormatter.format(coord);
	    }
	    for(int i = 0; i < numberOfSpacesToPrepend; i++)
	    	normalizedCoordinateString = " " + normalizedCoordinateString;
	    
	    // normalizedCoordinateString = normalizedCoordinateString + roundedCoord;
	    
	    return normalizedCoordinateString;
	}
	
	//////////////////////////////////////
	//// The attributes of Atom class ////
	//////////////////////////////////////
	
	/**
	 * The calculated energy of the atom
	 */
	private double energy;
	
	/**
	 * The serial number of the atom
	 */
	private int index; 	
}

