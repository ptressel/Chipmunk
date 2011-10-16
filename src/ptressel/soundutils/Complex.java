package ptressel.soundutils;

/**
 * Class representing complex numbers.  This is merely a container with a
 * few operations.
 *
 * Unlike the Java wrapper classes for primitive values (e.g. Integer),
 * Complex provides operations that alter the current instance, for
 * efficiency.
 *
 * For a full-featured Complex class, see JSci.maths.Complex in the JSci
 * library at: http://jsci.sourceforge.net/ .
 */

public class Complex implements Cloneable {

  // Constants

  /** Unit imaginary number, i.e., a Complex representing i. */
  public static final Complex I = new Complex( 0, 1 );

  /** Zero */
  public static final Complex ZERO = new Complex( 0, 0 );

  /** One */
  public static final Complex ONE = new Complex( 1, 0 );

  /** This system's newline */
  public static final String newline = System.getProperty( "line.separator" );

  // Class methods
  // These mainly provide arithmetic that does not alter the source values.
  // There are two forms of each -- one creates a new Complex for the
  // result, and the other (having one extra Complex argument) writes the
  // result into the supplied Complex.  The latter is what would be used
  // by code that is trying to avoid object creation, but also needs non-
  // destructive operations.

  /**
   * Make an array filled with separate Complex(0,0) instances.  This is
   * intended for applications that reuse their Complex arrays.
   *
   * @param size length of the array
   */
  public static Complex[] makeArray( int size ) {
    Complex[] array = new Complex[size];
    for ( int i = 0; i < size; i++ ) {
      array[i] = new Complex(0,0);
    }
    return array;
  }

  /**
   * Compute the unit vector at the given phase from the positive real
   * axis, measured counterclockwise.  That is, for an angle w, find
   * exp( i w ) = cos w + i sin w.  Return this in a new Complex.
   *
   * <p>Note this function is significantly inaccurate.  Where possible
   * (e.g. for simple angles), use known values.
   *
   * @param w phase in radians
   * @return a new Complex representing exp( i w )
   */
  public static Complex unitVector( double w ) {
    return new Complex( Math.cos(w), Math.sin(w) );
  }
  
  /**
   * <p>Compute the unit vector at the given phase from the positive real
   * axis, measured counterclockwise.  That is, for an angle w, find
   * exp( i w ) = cos w + i sin w.
   *
   * <p>Store the resulting real and imag parts in the supplied Complex and
   * return that.
   *
   * <p>Note this function is significantly inaccurate.  Where possible
   * (e.g. for simple angles), use known values.
   *
   * @param w phase in radians
   * @param out place for results
   * @return the supplied Complex, now containing exp( i w )
   */
  public static Complex unitVector( double w, Complex out ) {
    return out.set( Math.cos(w), Math.sin(w) );
  }
  
  /**
   * Compute the conjugate of the given Complex.  Return it in a new Complex.
   *
   * @param x Complex to conjugate
   * @return the conjugate in a new Complex
   */
  public static Complex conjugate( Complex x ) {
    return new Complex( x.real, -x.imag );
  }

  /**
   * Compute the conjugate of the given Complex.  Store it in the supplied
   * Complex and return that.
   *
   * @param x Complex to conjugate
   * @param out place for results
   * @return the supplied Complex, now containing the conjugate
   */
  public static Complex conjugate( Complex x, Complex out ) {
    return out.set( x.real, -x.imag );
  }

  /**
   * Add two Complex and return a new Complex with the result.
   *
   * @param a addend
   * @param b addend
   * @return a + b in a new Complex
   */
  public static Complex add( Complex a, Complex b ) {
    return new Complex( a.real + b.real, a.imag + b.imag );
  }

  /**
   * Add two Complex, store the result in the supplied Complex, and
   * return that.
   *
   * @param a addend
   * @param b addend
   * @param out place for results
   * @return the supplied Complex, now containing a + b
   */
  public static Complex add( Complex a, Complex b, Complex out ) {
    return out.set( a.real + b.real, a.imag + b.imag );
  }

  /**
   * Subtract two Complex and return a new Complex with the result.
   *
   * @param a minuend
   * @param b subtrahend
   * @return a - b in a new Complex
   */
  public static Complex subtract( Complex a, Complex b ) {
    return new Complex( a.real - b.real, a.imag - b.imag );
  }

  /**
   * Subtract two Complex, store the result in the supplied Complex, and
   * return that.
   *
   * @param a minuend
   * @param b subtrahend
   * @param out place for results
   * @return the supplied Complex, now containing a - b
   */
  public static Complex subtract( Complex a, Complex b, Complex out ) {
    return out.set( a.real - b.real, a.imag - b.imag );
  }

  /**
   * Multiply two Complex and return a new Complex with the result.
   * That is, for arguments a and b, the real part of the result is
   * a.real*b.real - a.imag*b.imag and the imaginary part is
   * a.real*b.imag + a.imag*b.real.
   *
   * @param a multiplicand
   * @param b multiplicand
   * @return product of a and b
   */
  public static Complex multiply( Complex a, Complex b ) {
    return new Complex( a.real*b.real - a.imag*b.imag,
                        a.real*b.imag + a.imag*b.real );
  }

  /**
   * Multiply two Complex, store the result in the supplied Complex, and 
   * return that.
   * That is, for arguments a and b, the real part of the result is
   * a.real*b.real - a.imag*b.imag and the imaginary part is
   * a.real*b.imag + a.imag*b.real.
   *
   * @param a multiplicand
   * @param b multiplicand
   * @param out place for results
   * @return the supplied Complex, now containing the product of a and b
   */
  public static Complex multiply( Complex a, Complex b, Complex out ) {
    return out.set( a.real*b.real - a.imag*b.imag,
                    a.real*b.imag + a.imag*b.real );
  }

  /**
   * Multiply one Complex by the conjugate of another, and return a new
   * Complex with the result.  That is, for arguments a and b, the real
   * part of the result is a.real*b.real + a.imag*b.imag and the imaginary
   * part is -( a.real*b.imag + a.imag*b.real ).
   *
   * @param a multiplicand
   * @param b multiplicand to be conjugated
   * @return product of a and conjugate of b
   */
  public static Complex multiplyConjugate( Complex a, Complex b ) {
    return new Complex(  a.real*b.real + a.imag*b.imag,
                        -a.real*b.imag + a.imag*b.real );
  }

  /**
   * Divide one Complex by another and return a new Complex with the result.
   * That is, for arguments a and b, the quotient is (a * conj(b)) / mag(b).
   *
   * @param a dividend
   * @param b divisor
   * @return quotient of a and b
   */
  public static Complex divide( Complex a, Complex b ) {
    double magB = b.magnitude();
    return new Complex( (  a.real*b.real + a.imag*b.imag ) / magB,
                        ( -a.real*b.imag + a.imag*b.real ) / magB );
  }

  /**
   * Multiply a Complex by a scalar and return a new Complex with the result.
   *
   * @param x multiplicand
   * @param c scalar multiplicand
   * @return a Complex with real and imag parts multiplied by c.
   */
  public static Complex scalarMultiply( Complex x, double c ) {
    return new Complex( c*x.real, c*x.imag );
  }

  /**
   * Helper method for comparing floating point values.  Check whether
   * the given double values are equal within the given number of ulps,
   * using the average of the two numbers' ulps.
   * (See Math.ulps(double).)
   *
   * @param d1 comparand
   * @param d2 comparand
   * @param nUlps allowed error
   * @return true if within nUlps ulps
   */
  public static boolean closeEnough( double d1, double d2, int nUlps ) {
    // Find nUlps times the average of the ulps.
    double okError = nUlps * ( Math.ulp(d1) + Math.ulp(d2) ) / 2;
    // Are they within specified number of ulps of each other?
    return ( Math.abs( d1 - d2 ) <= okError );
  }

  // Instance variables

  /** Real part */
  double real = 0.;

  /** Imaginary part */
  double imag = 0.;

  // Constructors

  /**
   * Make a Complex with the given real and imagninary parts.
   *
   * @param real real part
   * @param imag imagninary part
   */
  public Complex( double real, double imag ) {
    this.real = real;
    this.imag = imag;
  }

  /**
   * Make a Complex with the given magnitude and phase.  If this constructor
   * had only two arguments, its signature would conflict with the
   * constructor that takes the real and imaginary parts as arguments.
   * Rather than make a static method to do what is really the job of a
   * constructor, we add a dummy argument.  It is completely ignored --
   * put in whatever value you want.
   *
   * @param r magnitude
   * @param w phase
   * @param dummy something to distinguish this from Complex(x,y)
   */
  public Complex( double r, double w, int dummy ) {
    this.real = r * Math.cos(w);
    this.imag = r * Math.sin(w);
  }

  // Instance methods

  // Non-destructive methods -- these don't alter the object.

  /** Return the real part. */
  public double getReal() { return real; }

  /** Return the imaginary part. */
  public double getImag() { return imag; }

  /**
   * Compute the magnitude of this Complex = the square root of this value
   * times its conjugate = sqrt( real ^ 2 + imag ^ 2 ).
   *
   * @return the magnitude
   */
  public double magnitude() {
    return Math.hypot( real, imag );
    // return Math.sqrt( real*real + imag*imag );
  }

  /**
   * Return the square of the magnitude.  For use by functions that need
   * the square, so they don't have to roll their own, or call magnitude
   * and endure a useless sqrt call.
   *
   * @return square of the magnitude
   */
  public double magSqr() {
    return real*real + imag*imag;
  }

  /**
   * Return the phase in radians of this Complex, in the range pi to -pi.
   * That is, return the arctan of real / imag.  This calls Math.atan2(x,y)
   * with x as the real part and y the imaginary part, and so correctly
   * handles zero or infinite real or imaginary parts.
   *
   * @return the phase
   */
  public double phase() {
    return Math.atan2( real, imag );
  }

  // Support for Cloneable -- we don't need to override Object.clone() as
  // we have only primitive instance variables.  We provide equals to obey
  // the recommendation that a clone be "equals" to its source.

  /**
   * Check whether this Complex is exactly equal to the given Complex.
   * Overrides Object.equals.
   *
   * @param o comparand
   */
  public boolean equals( Object o ) {
    // Check for a compatible class.
    if ( !Complex.class.isAssignableFrom( o.getClass() ) ) return false;

    // o is a subclass of Complex -- safe to cast it.
    Complex c = (Complex) o;
    return (real == c.real) && (imag == c.imag);
  }

  /**
   * Check whether this Complex is equal to the given Complex within the
   * given number of ulps.  (See Math.ulps(double).)  That is:
   * Compare real and imag parts separately.  If each part of this Complex
   * and the supplied Complex are within an nUlps ulps (using whichever
   * one's ulp is greater), we take them to be the same.
   *
   * @param x comparand
   *
   * @return true if the real and imag parts are each differ by no more than
   * two ulps
   */
  public boolean closeEnough( Complex x, int nUlps ) {
    // Are they within specified number of ulps of each other?
    return closeEnough( real, x.real, nUlps ) &&
           closeEnough( imag, x.imag, nUlps );
  }

  /** String representation of a Complex. */
  public String toString() {
    return "Re: " + real + " Im: " + imag;
  }

  /** String representation for a Complex array. */
  public static String arrayString( Complex[] array ) {
    StringBuffer s = new StringBuffer("");
    for ( int i = 0; i < array.length; i++  ) {
      s.append( new String( i + ": " + array[i] ) + newline );
    }
    return new String(s);
  }

  /**
   * Find the number of ulps by which two Complex differ, separately for the
   * real and complex parts.  Use the average of the ulps for each of the
   * numbers being compared as the scale.
   *
   * @param a minuend
   * @param b subtrahend
   * @return a two element array containing, first, the difference between a
   * and b, measured in ulps, and second, the averaged ulps used as the scale
   */
  public static Complex[] diffUlps( Complex a, Complex b ) {
    double realUlps = ( Math.ulp(a.real) + Math.ulp(b.real) ) / 2;
    double imagUlps = ( Math.ulp(a.imag) + Math.ulp(b.imag) ) / 2;
    Complex[] result = { new Complex( ( a.real - b.real ) / realUlps,
                                      ( a.imag - b.imag ) / imagUlps ),
                         new Complex( realUlps, imagUlps ) };
    return result;
  }

  // Destructive methods:  These alter the current instance, and so are more
  // efficient than the static methods (which create new objects), so long as
  // the current instance doesn't need to be preserved.  Note this is not the
  // behavior of the Java number wrapper classes!!
  //
  // Where appropriate, the methods return the (modified) object itself, to
  // allow chaining operations together.

  // Unary operations and setters

  /**
   * Set both real and complex parts -- more efficient than separate calls
   * to setReal and setImag.
   *
   * @param real real part
   * @param imag imaginary part
   * @return the modified object
   */
  public Complex set( double real, double imag ) {
    this.real = real;
    this.imag = imag;
    return this;
  }

  /**
   * Copy the real and complex parts from an existing Complex.
   *
   * @param x the value to copy
   * @return the modified object
   */
  public Complex set( Complex x ) {
    this.real = x.real;
    this.imag = x.imag;
    return this;
  }

  /**
   * Set the real part.
   *
   * @param real real part
   * @return the modified object
   */
  public Complex setReal( double real ) {
    this.real = real;
    return this;
  }

  /**
   * Set the imaginary part.
   *
   * @param imag imaginary part
   * @return the modified object
   */
  public Complex setImag( double imag ) {
    this.imag = imag;
    return this;
  }

  /**
   * Conjugate this Complex.
   *
   * @return the modified object
   */
  public Complex conjugate() {
    imag = -imag;
    return this;
  }

  /**
   * Rotate counterclockwise by some multiple of 90 degrees.
   *
   * @param nTimes how many 
   * @return the modified object
   */
  public Complex rotate90( int nTimes ) {
    double temp;
    // Reduce nTimes to range 0-3.
    if ( nTimes > 3 ) nTimes %= 4;
    switch (nTimes) {
      case 0:
        return this;
      case 1:
        // 90: re=a, im=b -> re=-b, im=a
        temp = real;
        real = -imag;
        imag = temp;
        return this;
      case 2:
        // 180: re=a, im=b -> re=-a, im=-b
        real = -real;
        imag = -imag;
        return this;
      case 3:
        // 270: re=a, im=b -> re=b, im=-a
        temp = real;
        real = imag;
        imag = -temp;
        return this;
    }
    return this;  // Java doesn't know we can't get here.
  }

  // Binary operations

  /**
   * Add the supplied Complex to this Complex, with result stored
   * in this Complex.
   *
   * @param x addend
   * @return the modified object
   */
  public Complex add( Complex x ) {
    real += x.real;
    imag += x.imag;
    return this;
  }

  /**
   * Subtract the supplied Complex from this Complex, with result stored
   * in this Complex.
   *
   * @param x subtrahend
   * @return the modified object
   */
  public Complex subtract( Complex x ) {
    real -= x.real;
    imag -= x.imag;
    return this;
  }

  /**
   * Multiply this Complex by the supplied Complex, with result stored
   * in this Complex.
   *
   * @param x multiplicand
   * @return the modified object
   */
  public Complex multiply( Complex x ) {
    double newReal = real*x.real - imag*x.imag;
    imag = real*x.imag + imag*x.real;
    real = newReal;
    return this;
  }

  /**
   * Multiply this Complex by the conjugate of the supplied Complex, with
   * result stored in this Complex.
   *
   * @param x multiplicand to be conjugated
   * @return the modified object
   */
  public Complex multiplyConjugate( Complex x ) {
    double newReal = real*x.real + imag*x.imag;
    imag = -real*x.imag + imag*x.real;
    real = newReal;
    return this;
  }

  /**
   * Divide this Complex by the supplied Complex, with result stored in
   * this Complex.
   *
   * @param x divisor
   * @return the modified object
   */
  public Complex divide( Complex x ) {
    double magX = x.magnitude();
    double newReal = ( real*x.real + imag*x.imag ) / magX;
    imag = ( imag*x.real - real*x.imag ) / magX;
    real = newReal;
    return this;
  }

  /**
   * Multiply this Complex by a scalar, with result stored in this Complex
   *
   * @param c scalar multiplicand
   * @return the modified object
   */
  public Complex scalarMultiply( double c ) {
    real *= c;
    imag *= c;
    return this;
  }

  // Unit test

  /**
   * A very brief test of class Complex.
   */
  public static void main( String[] args ) {
    String msg;

    // Check magnitude.
    // A right triangle with integral hypotenuse.
    Complex rightTriangle345 = new Complex( 3, 4 );
    double hyp345 = rightTriangle345.magnitude();
    if ( hyp345 == 5 ) {
      System.out.println( "Magnitude ok." );
    } else {
      System.out.println( "Magnitude not ok: Got " + hyp345 +
                          "; should be 5." );
    }

    // Check angle.
    Complex piOver4Vector = new Complex( 1, 1 );  // A 45 degree angle.
    double ourPiOver4 = piOver4Vector.phase();
    double mathPiOver4 = Math.PI / 4;
    if ( ourPiOver4 == mathPiOver4 ) {
      System.out.println( "Phase exactly equal." );
    }
    if ( closeEnough( ourPiOver4, mathPiOver4, 2 ) ) {
      System.out.println( "Phase close enough." );
    } else {
      System.out.println( "Phase not ok: Got " + ourPiOver4 +
                          "; should be " + mathPiOver4 + "." );
    }

    // Get and set.
    Complex x = new Complex( 4, -1 );  // whatever
    x.setReal(3);
    x.setImag(4);
    if ( x.getReal() == 3 && x.getImag() == 4 ) {
      System.out.println( "Set and get ok." );
    } else {
      System.out.println( "Set and get not ok" );
    }

    // Check equals and (destructive) conjugate.
    Complex cp4 = new Complex( 1, -1 );
    cp4.conjugate();  // flip the sign
    if ( piOver4Vector.equals( cp4 ) ) {
      System.out.println( "Equals and conjugate ok." );
    } else {
      System.out.println( "Equals and conjugate not ok." );
    }
    // cp4 is munged -- don't use it any more.

    // Now that we've tested equals, test clone.
    // Silly compiler *knows* we support clone, but it gripes about not
    // handling the exception anyway.
    Complex cl = null;
    try {
      cl = (Complex) rightTriangle345.clone();
    }
    catch( CloneNotSupportedException e ) {
      System.out.println( "Clone attempt threw CloneNotSupportedException." );
    }
    if ( cl.equals( rightTriangle345 ) ) {
      System.out.println( "Clone ok." );
    } else {
      System.out.println( "Clone not ok." );
    }
    // Now clone an array of Complex.
    Complex[] ar1 = { new Complex( 1, 2 ), new Complex( -3, 4 ) };
    Complex[] ar2 = null;
    // But here, the silly compiler fusses if I *do* include the try/catch.
    // It should make up its mind!!!
    //try {
      ar2 = (Complex[]) ar1.clone();
    //}
    //catch( CloneNotSupportedException e ) {
    //  System.out.println(
    //    "Array clone attempt threw CloneNotSupportedException." );
    //}
    if ( ar1[0].equals( ar2[0] ) && ar1[1].equals( ar2[1] ) ) {
      System.out.println( "Array clone ok." );
    } else {
      System.out.println( "Array clone not ok." );
    }

    // Non-destructive multiply.
    Complex rt = new Complex( 3, 4 );
    Complex rt2 = new Complex( 3, 4 );
    Complex crt = new Complex( 3, -4 );
    // Product of (3,4) w/ self and w/ conjugate.
    Complex prd345 = new Complex( 3*3-4*4, 2*3*4 );
    Complex prd345c = new Complex( 3*3+4*4, 0 );

    Complex p1 = multiply( rt, crt );
    if ( p1.equals( prd345c ) ) {
      System.out.println( "1st non-destructive multiply ok." );
    } else {
      System.out.println( "1st non-destructive multiply not ok. Got: " +
                          p1 + "; should be " + prd345c + "." );
    }
    Complex p2 = multiply( rt, rt );
    if ( p2.equals( prd345 ) ) {
      System.out.println( "2nd non-destructive multiply ok." );
    } else {
      System.out.println( "2nd non-destructive multiply not ok. Got: " +
                          p2 + "; should be " + prd345 + "." );
    }
    Complex p3 = multiplyConjugate( rt, crt );
    if ( p3.equals( prd345 ) ) {
      System.out.println( "3rd non-destructive multiply ok." );
    } else {
      System.out.println( "3rd non-destructive multiply not ok. Got: " +
                          p3 + "; should be " + prd345 + "." );
    }

    // Destructive multiply.
    rt.multiply( crt );
    rt2.multiplyConjugate( crt );
    if ( rt.equals( prd345c ) && rt2.equals( prd345 ) ) {
      System.out.println( "Destructive multiply ok." );
    } else {
      System.out.println( "Destructive multiply not ok." );
    }
    // rt and rt2 are munged.

    // Add, subtract
    Complex a = new Complex( 1, 2 );
    Complex a1 = new Complex( 1, 2 );
    Complex b = new Complex( 2, 1 );
    Complex sum = new Complex( 3, 3 );
    Complex dif = new Complex( -1, 1 );
    if ( add( a, b ).equals(sum) && subtract( a, b ).equals(dif) ) {
      System.out.println( "Non-destructive add & subtract ok." );
    } else {
      System.out.println( "Non-destructive add & subtract not ok." );
    }
    a.add(b);
    a1.subtract(b);
    if( a.equals(sum) && a1.equals(dif) ) {
      System.out.println( "Destructive add & subtract ok." );
    } else {
      System.out.println( "Destructive add & subtract not ok." );
    }
    // a, a1 munged.

    // Divide, scalarMultiply, non-destructive conjugate.
    Complex u = new Complex( 2, 3 );
    Complex v = new Complex( -1, 1 );
    Complex q
      = scalarMultiply( multiply( u, conjugate(v) ), 1 / v.magnitude() );
    Complex qnd = divide( u, v ); // non-destructive
    u.divide(v);  // destructive
    if ( q.equals(qnd) && q.equals(u) ) {
      System.out.println( "Non-destructive and destructive divide ok." );
    } else {
      System.out.println( "Non-destructive and destructive divide not ok." );
    }

    // Rotate
    // Some vectors for comparison.
    Complex r0 = new Complex(2,3);
    Complex r1 = new Complex(-3,2);
    Complex r2 = new Complex(-2,-3);
    Complex r3 = new Complex(3,-2);
    // The ones we'll change.
    Complex z0 = new Complex(2,3);
    Complex z1 = new Complex(2,3);
    Complex z2 = new Complex(2,3);
    Complex z3 = new Complex(2,3);
    Complex z7 = new Complex(2,3);
    // Do the rotates and compare against the expected values.
    if ( z0.rotate90(0).equals( r0 ) && z1.rotate90(1).equals( r1 ) &&
         z2.rotate90(2).equals( r2 ) && z3.rotate90(3).equals( r3 ) &&
         z7.rotate90(7).equals( r3 ) ) {
      System.out.println( "Rotate ok." );
    } else {
      System.out.println( "Rotate not ok." );
    }

    // Compute some unit vectors using unitVector and compare with known
    // values using diffUlps.  Compute the arguments to unitVector exactly
    // as is done by the FFT constructor.
    Complex u0 = Complex.unitVector( 2*Math.PI*0/4 );
    Complex u90 = Complex.unitVector( 2*Math.PI*1/4 );
    Complex u180 = Complex.unitVector( 2*Math.PI*2/4 );
    Complex u270 = Complex.unitVector( 2*Math.PI*3/4 );
    // Now the real values.
    Complex t0 = new Complex( 1, 0 );
    Complex t90 = new Complex( 0, 1 );
    Complex t180 = new Complex( -1, 0 );
    Complex t270 = new Complex( 0, -1 );
    // Get differences.
    Complex d0 = subtract( t0, u0 );
    Complex d90 = subtract( t90, u90 );
    Complex d180 = subtract( t180, u180 );
    Complex d270 = subtract( t270, u270 );
    Complex[] du0 = diffUlps( t0, u0 );
    Complex[] du90 = diffUlps( t90, u90 );
    Complex[] du180 = diffUlps( t180, u180 );
    Complex[] du270 = diffUlps( t270, u270 );

    System.out.println(
      "Errors due to calculating unit vectors using trig functions:" );
    System.out.println( "Angle	Diff	Diff in ulps	Ulps" );
    System.out.println( "0 deg	" + d0 + "	" + du0[0] + "	" + du0[1]);
    System.out.println( "90 deg	" + d90 + "	" + du90[0] + "	" + du90[1]);
    System.out.println(
      "180 deg	" + d180 + "	" + du180[0] + "	" + du180[1]);
    System.out.println(
      "270 deg	" + d270 + "	" + du270[0] + "	" + du270[1]);
  }
}
