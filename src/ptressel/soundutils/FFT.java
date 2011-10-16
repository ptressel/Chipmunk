package ptressel.soundutils;

import java.io.*;
import ptressel.soundutils.SoundUtils.IndexAndValue;

/**
 * <p>An object of class FFT is provides an fft calculator tailored to a
 * particular number of samples, chosen at construction time.
 *
 * <p>Suggested usage:  Make one FFT object for each input processing stream,
 * with the desired sample length.  For each sample to be transformed, call
 * fft(sample) on the FFT object.  After this, do queries on the FFT object
 * using its various convenience methods, or examine its results array.
 *
 * <p>FFT objects should be reused, else the efficiency of avoiding the setup
 * is wasted.
 *
 * <p>The FFT results persist until another fft call is done, when the
 * results array is overwritten.  So if results must be preserved past the
 * next use of the FFT object, they should be copied out or cloned.  (The
 * array element type, Complex, supports cloning.  A method is provided that
 * returns a copy of the results, as well as a method that merely returns
 * a reference to the internal results array.)  However, if it will often be
 * required to have more than one set of FFT results available at once, it
 * will be more efficient to create enough FFT objects.
 *
 * <p>References on the fft algorithm used here:
 *
 * <p>Schilling & Harris, "Fundamentals of Signal Processing", 2005,
 * section 3.5.  See especially the pseudo-code on p. 200.  An attempt has
 * been made to use notation similar to that in this reference.
 *
 * <p>Lynn & Fuerst, "Introductory Digital Signal Processing with Computer
 * Applications", 1989 (there's also a 2nd edition in 1998), section 7.3.
 */

// TODO Currently this and other SoundUtils classes only support 8-bit and
// 16-bit sound depths.  24-bit is needed for general use.  Support 24-bit
// depth throughout.  Find out if 32-bit or deeper sound is ever used, if so,
// support it, if not, fail if anyone requests it.

public class FFT
{
  // Class constants

  /**
   * Default half-width at half height of the smoothing filter.
   */
  public static final int defaultSigma = 2;

  /** The local system's newline. */
  public static final String newline = System.getProperty( "line.separator" );

  // Instance variables

  /**
   * Number of samples in the signal and result.
   */
  private int N;  // Yes, I know it's uppercase.

  /**
   * Actual number of elements in the results array, excluding the high
   * end that are just conjgates of the low end.
   */
  private int rN;

  /**
   * Log base 2 of the number of samples.  This is the number of stages in
   * the fft calculation.
   */
  private int logN;

  /** Sample rate in Hz. */
  private double sampleRate;

  /**
   * Storage for results.  This is allocated once when the FFT is created,
   * and reused.  The results array is length N, with elements of type
   * Complex.  However, meaningful results lie only in element 0 to N/2,
   * below the Nyquist frequency.  Query methods will not return values
   * above that.
   */
  private Complex[] results;

  /**
   * Internal storage for intermediate results.  Two arrays are used during
   * computation, alternately serving as the output of each stage.  Their
   * references are swapped between scratch and results each time, so
   * whichever array is pointed to by results is the output.
   */
  private Complex[] scratch;

  /** Storage for power density spectrum. */
  private double[] power;

  /** Maximum value in the power spectrum. */
  private double maxPower;

  /** Average power. */
  private double avgPower;

  /** Storage for a smoothed power density spectrum. */
  private double[] smoothedPower;

  /** Maximum value in the smoothed power spectrum. */
  private double maxSmoothed;

  // There's no average smoothed power because it's equal to the average
  // power:  Our smoothing filter has total weight one.

  /** Storage for the smoothing filter. */
  private double[] filter;

  /** Midpoint of the filter. */
  private int fMid;

  /**
   * The half-width at half height of the smoothing filter, in terms of
   * indices in the results and power arrays.
   */
  private int sigma;

  /**
   * The unique roots of unity needed in the fft calculation for sample and
   * results length N.  These are exp( -i 2 pi n/N ) where 0 <= n < N.
   */
  private Complex[] roots;

  /**
   * Table of indices for copying the signal in the appropriately scrambled
   * order as input to the first stage of the decimation-in-time fft
   * calculation.  Note these are merely the linear indices with the bit
   * order reversed.
   */
  private int[] bitrev;

  // Constructor

  /**
   * Make an FFT for length N signals.  This can be reused for other signals
   * of the same length.
   *
   * Currently N must be a positive integral power of 2.
   *
   * The sample rate isn't needed for the FFT operation itself, but rather
   * allows answering queries in terms of actual times and frequencies, not
   * abstract units.
   *
   * A Gaussian with the given sigma (measured in terms of bins in the
   * spectrum, and truncated at two sigma out) is used to produce a smoothed
   * version of the power spectrum.
   *
   * @param N the number of samples in the signal and result
   * @param sampleRate sample frequency in Hz
   * @param sigma half-width at half height of the smoothing filter
   * @throws IllegalArgumentException if N is not a power of 2
   */
  public FFT( int N, double sampleRate, int sigma )
      throws IllegalArgumentException {
    init( N, sampleRate, sigma );
  }

  /**
   * Make an FFT for length N signals.  This can be reused for other signals
   * of the same length.
   *
   * Currently N must be a positive integral power of 2.
   *
   * The sample rate isn't needed for the FFT operation itself, but rather
   * allows answering queries in terms of actual times and frequencies, not
   * abstract units.
   *
   * A default filter (a Gaussian with sigma of two bins in the spectrum,
   * truncated at two sigma out) is used to produce a smoothed version of
   * the power spectrum.
   *
   * @param N the number of samples in the signal and result
   * @param sampleRate sample frequency in Hz
   * @throws IllegalArgumentException if N is not a power of 2
   */
  public FFT( int N, double sampleRate )
      throws IllegalArgumentException {
    init( N, sampleRate, defaultSigma );
  }

  /** Worker method for the constructors. */
  private void init( int N, double sampleRate, int sigma )
      throws IllegalArgumentException {

    /*debug*/
    /*
    System.out.println();
    System.out.println( "FFT(" + N + ")" );
    */

    // Check that N is a power of 2 -- see if more than one bit is on.
    // Note the spiffy new Integer methods introduced in Java 1.5.  Not
    // using 1.5?  Hah -- write yer own.
    if ( Integer.highestOneBit(N) != Integer.lowestOneBit(N) || N <= 0 ) {
      throw new IllegalArgumentException( "Size must be power of 2" );
    }
    // Survived -- save the number as N, and get log N, which is the number
    // of zeros to the right of the single 1 bit.
    this.N = N;
    logN = Integer.numberOfTrailingZeros(N);
    // Meaningful results only go run from index zero to N/2, i.e. length of
    // useful part of results is N/2 + 1.
    rN = N/2 + 1;

    /*debug*/
    /*
    System.out.println( "N = " + N + ", logN = " + logN + ", rN = " + rN );
    */

    // Save rate.
    this.sampleRate = sampleRate;

    // Create the scratch and results array.  The alert code reader will
    // note that we allocate a full N entries in these arrays.  This is
    // because it is much simpler, and avoids lots of conditionals and
    // method calls on the critical path, if we just let fft fill the entire
    // N entries, rather than trying to catch those above N/2 and substitute
    // the conjugate of the mirror image below N/2.
    scratch = Complex.makeArray(N);
    results = Complex.makeArray(N);

    // Arrays of extracted values are only N/2+1 long, however.  We fill
    // power as soon as we calculate an fft, but until smooth is called,
    // smoothedPower is just a duplicate of power.
    power = new double[rN];
    smoothedPower = new double[rN];

    // Compute the required roots of unity, exp( -i2pi n/N), n = 0 to N-1.
    // Note for simple angles (i.e. the compass points), we will replace the
    // computed values with known values, to avoid at least in obvious cases
    // possible the wildly inaccurate trig functions (off by ~1E-15! sounds
    // small, but it's ~1E15 ulps).
    //
    // At the same time...
    // Fill in the array of indices for rearranging the signal for the first
    // stage of the fft calculation.  The nth entry in this table is n with
    // its bits backwards.  Java's Integer class provides a reverse method
    // (since 1.5).  Hope it is doing the reverse in some nifty fashion,
    // making use of the local machine's instruction set if possible...
    roots = new Complex[N];
    bitrev = new int[N];
    for ( int n = 0; n < N; n++ ) {

      // Ask for vectors with N equal angles around the unit circle.
      // (Caution to code modifiers:  Do not rearrange the order of
      // operations such that any of these do integer arithmetic!)
      roots[n] = Complex.unitVector( -2*Math.PI*n/N );

      // Reverse the bits in n and put that in the index table.  We can use
      // the new Integer.reverse method, which reverses the whole 32 bit int.
      // But our indices are down in the low end:  bits 0 through log N - 1.
      // So we need to shift the reversed int back down.  The zero bit lands
      // in bit # 31, and we want it to be in bit # log N - 1, so we shift it
      // down by 31 - (logN - 1) = 32 - logN.
      bitrev[n] = Integer.reverse(n) >>> (32-logN);
      /*debug*/
      /*
      System.out.println( "n = " + n + ", roots[n] = " + roots[n] +
        ", rev = " + Integer.reverse(n) +
        ", shift = " + (32-logN) + ", bitrev = " + bitrev[n] );
      */
    }

    // Fix up the compass points in roots, 0, -pi/2, -pi, -3pi/2.
    roots[0].set( 1, 0 );        // 0
    roots[N/2].set( -1, 0 );     // -pi
    if ( N >= 4 ) {
      roots[N/4].set( 0, -1 );   // -pi/2, same as 3pi/2
      roots[3*N/4].set( 0, 1 );  // -3pi/2, same as pi/2
    }

    // Construct Gaussian filter for smoothing.  We do not bother with the
    // normalization constant up front since we're going to chop off the
    // tails, and would just need to re-normalize anyway -- the weights we
    // actually use must add to one.  Width of the whole filter is
    // 4 * sigma + 1, i.e. from the center point at the peak of the Gaussian,
    // the filter extends for 2*sigma more points on either side.  Range of
    // filter indices is 0 to 4*sigma, and the peak is at 2*sigma.

    // Make sure the filter will fit in the power spectrum array!  Not a
    // problem for real data, but our debugging cases have N=2 and N=4...
    // Adjust sigma downward to fit.
    int fLen = 4*sigma+1;   // filter length with current sigma
    if ( fLen > rN ) {
      // Too big -- make filter length the largest odd # < = length of power
      // array, as filter should be symmetric and have a single peak.  Min
      // length is 1, which is what we'll get with a minimum length FFT of 2.
      fLen = (rN % 2 == 0) ? rN-1 : rN;
      // Fix up sigma -- round up if fractional.
      sigma = (int) Math.ceil( ((double)(fLen-1))/4 );
    }

    filter = new double[ fLen ];  // filter storage
    fMid = (fLen-1)/2;            // filter midpoint, i.e. mean
    filter[ fMid ] = 1;     // unnormalized peak = exp(0)
    double sum = 1;         // add up unnormalized values here, start w/ peak
    // The rest are symmetric about the peak.
    for ( int i = 0; i < fMid; i++ ) {
      double delta = i - fMid;  // point minus mean
      filter[i] = filter[ fLen-1-i ]
        = Math.exp(-delta*delta/(2*sigma*sigma));
      sum += 2*filter[i];
    }
    // Now divide through by the sum to normalize to 1.
    for ( int i = 0; i < fLen; i++ ) {
      filter[i] /= sum;
    }
  }

  // Instance methods

  // Information available after construction, that doesn't depend on results

  /**
   * Return the number of signal samples for which this FFT was
   * created.
   *
   * @return the length of the input signal
   */
  public int size() { return N; }

  /**
   * Return the number of meaningful values in the results array, i.e. the
   * number of independent spectrum values, i.e. those at 0 to the Nyquist
   * frequency.  This is equal to half the number of samples in the input
   * signal plus 1.
   *
   * @return the number of points in the results, from 0 to the Nyquist
   * frequency
   */
  public int resultsSize() { return rN; }

  /**
   * Return the Nyquist frequency, i.e. the frequency corresponding to the
   * highest results value, at index N/2.
   *
   * @return the Nyquist frequency
   */
  public double nyquist() { return sampleRate / 2; }

  /**
   * Return the frequency corresponding to a given index n, i.e.
   * n s / N where s is the sample rate.  No attempt is made to constrain
   * the index to the range actually available in the results.
   *
   * @param n index
   * @return frequency corresponding to n
   */
  public double frequencyAtIndex( int n ) {
    return n * sampleRate / N;
  }

  /**
   * <p>Return the index nearest to the given frequency.  Note the index is
   * not constrained to be in the span of the results array -- it is simply
   * the nearest integer after conversion from the frequency according to:
   * n = round( f*N/s ) where f is the given frequency and s is the sample
   * rate.  So receiving a value from indexForFrequency does not guarantee
   * that it can be used to look up an fft result.  The caller may wish to
   * use this conversion for other purposes.
   *
   * <p>The above formula can be obtained as follows:  The fft produces the
   * same number of points as in the supplied signal array.  The values in
   * the top half of the results, above the Nyquist frequency, which is half
   * the sample rate, are a mirror image of those below.  The frequencies in
   * the whole results array range from 0 to twice the Nyquist frequency, i.e
   * from 0 to the sample frequency.  This range is split into the same
   * number of points as there are in the signal array.  So the spacing
   * between frequency points is:
   * <pre>
   *   frequency spacing = (s/2) / (N/2) = s/N
   * </pre>
   * where s is the sample rate and n is the buffer size.  So the frequency
   * f(n) at index n is n*s/N, or n = f(n)*N/s.  For frequencies not
   * exactly at one of the f(n) points, we round to the nearest integer.
   *
   * @param f frequency
   * @return index nearest to the given frequency
   * @throws ArrayIndexOutOfBoundsException if resulting index cannot be
   * converted to an int without under- or overflow
   */
  public int indexForFrequency( double f )
      throws IndexOutOfBoundsException {

    // Range check the raw value before we convert it.
    double x = f * N / sampleRate;
    if ( x > Integer.MAX_VALUE || x < Integer.MIN_VALUE )
      throw new IndexOutOfBoundsException(
        "Resulting index not within int range" );
    // Safe to convert.
    return (int) Math.round(x);
  }

  // Raw results

  /**
   * <p>Return the value in the results at index n.  No attempt is made to
   * prevent reading out values before the first fft call.  Values are only
   * returned for the meaningfuo part of the results, i.e. at indices
   * 0 through N/2.
   *
   * <p>This is the officially sanctioned means of examining the results
   * (as opposed to calling getResultsArray to get a reference to the
   * array itself) but is still not safe -- it returns references to the
   * actual Complex objects residing in the results array.  So the caller
   * should be careful not to change their contents.
   *
   * @param n index into results
   * @return the Complex value at index n in the results
   * @throws ArrayIndexOutOfBoundsException if n is > N/2
   */
  public Complex getResult( int n ) throws ArrayIndexOutOfBoundsException {
    return results[n];
  }

  /**
   * <p>Get a reference to the internal complex spectrum array.
   *
   * <p>No check is made for
   * whether this array contains actual results of a computation, as
   * that would only be relevant at the outset:  Once fft has been called
   * successfully, this FFT contains results.  The array is provided in case
   * the caller needs to extract some info not provided by FFT methods.  It
   * is suggested, instead, that the caller add the required methods to FFT...
   *
   * <p>Note that this array is the personal property of the FFT object, and,
   * no matter what the caller may be using it for, the next fft call will
   * overwrite its contents.  Further note that if the caller messes with the
   * contents of the array, subsequent information extraction calls to FFT
   * will not return results corresponding to the spectrum derived in the
   * prior fft call.  This is the price we pay for efficiency.  The caller
   * Has Been Warned.
   *
   * @return a reference to the complex spectrum array
   */
  public Complex[] getResultsArray() { return results; }

  /**
   * Make a copy of the results array.  This is appropriate if either the
   * results have to be preserved, and this FFT object will be reused -or-
   * if one wants to do some manipulation of the results, without damaging
   * the internal FFT values.
   */
 public Complex[] copyResultsArray() { return results.clone(); }

  /**
   * <p>Got your own spectrum, and want to use the FFT query methods, but
   * also want to keep the current FFT object contents?  Trade the arrays
   * temporarily -- put yours in with swapResultsArray and save the returned
   * reference.  When done, put the real one back in using swapResultsArray.
   *
   * <p>Besides inserting the caller's spectrum, swapResultsArray recomputes
   * the power and smoothed power from the new spectrum.
   *
   * <p>Because the power and smoothed power would need to be recomputed if
   * a saved spectrum is reinserted, it would be more efficient to make an
   * extra FFT just for the purpose of examining the caller's own spectra.
   *
   * @param newArray alternate array (containing a spectrum)
   * @throws IllegalArgumentException if the array is the wrong size
   */
  public Complex[] swapResultsArray( Complex[] newArray )
      throws IllegalArgumentException {

    // Make sure it's appropriate.  No, they can't have a longer array...yeesh.
    if ( newArray.length != N )
      throw new IllegalArgumentException(
        "Array length doesn't match this FFT length." );

    // Trade the arrays.
    Complex[] temp = results;
    results = newArray;
    // Populate the power and smoothedPower arrays.
    computePower();
    smooth();

    return temp;
  }

  /**
   * <p>Get a reference to the internal power array.
   *
   * <p>No check is made for 
   * whether this array contains actual results of a computation, as 
   * that would only be relevant at the outset:  Once fft has been called 
   * successfully, this FFT contains results.  The array is provided in case 
   * the caller needs to extract some info not provided by FFT methods.  It
   * is suggested, instead, that the caller add the required methods to FFT...
   *
   * <p>Note that this array is the personal property of the FFT object, and,
   * no matter what the caller may be using it for, the next fft call will
   * overwrite its contents.  Further note that if the caller messes with the
   * contents of the array, subsequent information extraction calls to FFT
   * will not return results corresponding to the spectrum derived in the
   * prior fft call.  This is the price we pay for efficiency.  The caller
   * Has Been Warned.
   *
   * @return a reference to the power array
   */
  public double[] getPowerArray() { return power; }

  /**
   * <p>Get a reference to the internal smoothed power array.
   *
   * <p>No check is made for
   * whether this array contains actual results of a computation, as
   * that would only be relevant at the outset:  Once fft has been called
   * successfully, this FFT contains results.  The array is provided in case
   * the caller needs to extract some info not provided by FFT methods.  It
   * is suggested, instead, that the caller add the required methods to FFT...
   * 
   * <p>Note that this array is the personal property of the FFT object, and,
   * no matter what the caller may be using it for, the next fft call will
   * overwrite its contents.  Further note that if the caller messes with the
   * contents of the array, subsequent information extraction calls to FFT
   * will not return results corresponding to the spectrum derived in the
   * prior fft call.  This is the price we pay for efficiency.  The caller
   * Has Been Warned.
   *
   * @return a reference to the smoothed power array
   */
  public double[] getSmoothedPowerArray() { return smoothedPower; }

  /**
   * <p>Get a reference to the internal filter array.
   *
   * <p>Note that this array is the personal property of the FFT object.
   * If the caller messes with the contents of the array, subsequent fft
   * calls will not smooth the power spectrum according to the stated scheme.
   * The caller Has Been Warned.
   *
   * @return a reference to the filter array
   */
  public double[] getFilterArray() { return filter; }

  // Private helpers, mainly used during the fft itself to populate the
  // secondary results -- the power and smoothed power arrays.

  /**
   * Power at the given results index.  Power at a single index is the
   * squared magnitude of the result value there, divided by N.  This
   * method is used only internally by the fft method, to populate the
   * power array.  External callers use powerAtIndex.
   *
   * @param n index
   * @return power at that index
   * @throws ArrayIndexOutOfBoundsException if index is outside span
   * of results
   */
  private double powerAtIndexFromResults( int n )
      throws ArrayIndexOutOfBoundsException {

    // This is the only place we actually compute the power.  All other
    // power queries are serviced from the power array.
    return results[n].magSqr() / N;
  }

  /**
   * Compute the power for the spectrum up to the Nyquist frequency.
   * Also find the average and maximum power.
   */
  private void computePower() {
    maxPower = 0;
    avgPower = 0;
    for ( int n = 0; n < power.length; n++ ) {
      power[n] = powerAtIndexFromResults(n);
      if ( power[n] > maxPower ) maxPower = power[n];
      avgPower += power[n];
    }
    avgPower /= power.length;
  }

  /**
   * Smooth the power using a Gaussian filter with sigma chosen via the
   * constructor.  Also find the maximum of the smoothed power.
   */
  private void smooth() {
    maxSmoothed = 0;

    // With the filter centered at each point in the power array, add the
    // contributions of the neighborhood and put them in the smoothed array.
    // Here, do only the points that will use the whole filter.  The lowest
    // smoothed point that gets to use the whole filter is the one at 2*sigma,
    // and it gets contributions from points 0 through 4*sigma in power.
    // Those are matched up with points 0 through 4*sigma in the filter.
    // The smoothed element at n gets contributions from points n-2*sigma
    // through n+2*sigma in power, again matched up with filter points 0
    // through 4*sigma.  The last smoothed point that uses the whole filter
    // is that at 2*sigma points prior to the end of the array.  But do
    // everything in terms of array lengths.
    for ( int n = fMid; n < power.length-fMid; n++ ) {
      // Sum up the contributions weighted by the filter.  Yes, smoothedPower
      // is already zero.  So sue me.
      smoothedPower[n] = 0;
      for ( int i = -fMid; i <= fMid; i++ ) {
        smoothedPower[n] += power[n+i] * filter[fMid+i];
      }
      if ( smoothedPower[n] > maxSmoothed ) maxSmoothed = smoothedPower[n];
    }

    // Now those pesky edge values.  Add up the weighted filter values that
    // actually get used, then fix up the normalization.  Do both ends at
    // once, as the normalization factor is symmetric.
    for ( int n = 0; n < fMid; n++ ) {
      // The two smoothed locations we're working on:
      smoothedPower[n] = 0;
      smoothedPower[power.length-1-n] = 0;
      double sum = 0;  // add up filter values here
      // Cut off the filter on the "left" side for the low end of the power
      // array, or on the right for the high end.  When n=0, the filter peak
      // is positioned over 0 in the power array, and points 0 through 2*sigma
      // in the power array are matched with filter points 2*sigma (the peak)
      // through 4*sigma.
      for ( int i = -n; i <= fMid; i++ ) {
        smoothedPower[n] += power[n+i] * filter[fMid+i];
        smoothedPower[power.length-n-1]
          += power[power.length-1-n-i] * filter[fMid+i];
        sum += filter[fMid+i];
      }
      smoothedPower[n] /= sum;
      smoothedPower[power.length-1-n] /= sum;
    }
  }

  // Possibly interesting queries

  /**
   * Average of the power.  This is the same for either the smoothed or
   * unsmoothed power.
   *
   * @return average power
   */
  public double averagePower() { return avgPower; }

  /**
   * Maximum value in the (unsmoothed) power.
   *
   * @return maximum power
   */
  public double maximumPower() { return maxPower; }

  /**
   * Maximum value in the smoothed power.
   *
   * @return maximum smoothed power
   */
  public double maximumSmoothedPower() { return maxSmoothed; }

  /**
   * Power at the given index.  Power at a single index is the
   * squared magnitude of the result value there, divided by N.
   *
   * @param n index
   * @return power at that index
   * @throws ArrayIndexOutOfBoundsException if index is outside span
   * of meaningful data, 0 to Nyquist frequency
   */
  public double powerAtIndex ( int n )
      throws ArrayIndexOutOfBoundsException {

    return power[n];
  }

  /**
   * Power at the given frequency.  If frequency does not correspond
   * precisely to an index in the result, it is rounded to the nearest
   * index.  Power at a single index is the squared magnitude of the
   * result value there, divided by N.
   *
   * @param f frequency
   * @return power at index nearest to that frequency
   * @throws ArrayIndexOutOfBoundsException if nearest index is outside span
   * of meaningful data, 0 to Nyquist frequency
   */
  public double powerAtFrequency( double f )
      throws ArrayIndexOutOfBoundsException {

    return powerAtIndex( indexForFrequency( f ) );
  }

  /**
   * <p>Get the power in a given index range.  The range is inclusive of the
   * low endpoint, but exclusive of the high endpoint.  This simplifies
   * calling powerInIndexRange repeatedly with successive ranges.  The power
   * in the range is just the sum of the power at each index in the range.
   *
   * <p>See indexForFrequency for an explanation of conversion from the
   * index (bin) number in the complex results or power arrays and the
   * frequency.
   *
   * <p>The scale of the units in which the power is expressed depends on
   * what the input signal scale was.  E.g. if the signal started out as
   * 8 bits, then its values ranged from either -2^7 to 2^7-1 if signed,
   * or 0 to 2^8-1 if unsigned.
   *
   * @param n1 low (inclusive) end of index range
   * @param n2 high (exclusive) end of index range
   * @return power in that range
   * @throws ArrayIndexOutOfBoundsException if range exceeds span of
   * meaningful data, 0 to Nyquist frequency
   */
  public double powerInIndexRange( int n1, int n2 )
      throws ArrayIndexOutOfBoundsException {

    // Note we do not reorder the indices if n1 isn't the low index --
    // caller might be doing some procedure in which indices out of order
    // occurs in an end case.

    // Loop over the indices, adding up the power at each.
    double powerSum = 0;
    for ( int n = n1; n < n2; n++ ) {
      powerSum += powerAtIndex(n);
    }

    return powerSum;
  }

  /**
   * <p>Get the power in a given frequency range.  The range is inclusive of
   * the low endpoint, but exclusive of the high endpoint.  This simplifies
   * calling powerInFreqRange repeatedly with successive ranges.  If
   * frequencies do not correspond precisely to indices in the results, they
   * are rounded to the nearest indices.  The power in the range is just
   * the sum of the power at each index in the range.
   *
   * <p>The scale of the units in which the power is expressed depends on
   * what the input signal scale was.  E.g. if the signal started out as
   * 8 bits, then its values ranged from either -2^7 to 2^7-1 if signed,
   * or 0 to 2^8-1 if unsigned.
   *
   * @param f1 low (inclusive) end of frequency range
   * @param f2 high (exclusive) end of frequency range
   * @return power in that range
   * @throws ArrayIndexOutOfBoundsException if range exceeds span of
   * meaningful data, 0 to Nyquist frequency
   */
  public double powerInFreqRange( double f1, double f2 )
      throws ArrayIndexOutOfBoundsException {

    // Get the nearest indices.
    int n1 = indexForFrequency( f1 );
    int n2 = indexForFrequency( f2 );

    // Let powerInIndexRange do the rest of the work.
    return powerInIndexRange( n1, n2 );
  }

  // TODO
  // Q: Have powerPeaks and smoothedPowerPeaks return IndexAndValue arrays?
  // If so, change any callers and tests.
  /**
   * <p>Find maxima in the (unsmoothed) power and return an array with their
   * indices.  Maxima are points that are larger than their neighbors, so
   * if the power has a lot of noise, there may be many.  If there is a
   * flat region in the data, and points on either side are lower, then
   * all points in the flat region are reported as maxima.
   *
   * <p>A noise margin can be specified -- points are not regarded as
   * different unless they differ by more than this amount.  A caution:
   * This is not a substitute for smoothing, and the margin should be
   * small.  No attempt is made to guard against missing a peak whose
   * slopes are shallow enough that adjacent points differ by less than
   * the margin.
   *
   * <p>A threshold can be specified -- peaks at less than this minimum are
   * not recorded.  Note that any margin applies to the test for points
   * being above the threshold, as otherwise the definition of equal within
   * the margin would be violated for points just on either size of the
   * threshold.
   *
   * @param margin Values are considered equal if they differ by no more than
   * this.
   * @param threshold Peaks are not reported unless at least this high.
   * @param plateauCenterOnly If true, only the middle index of a plateau
   * peak is returned (or middle-1 if the number of indices in the plateau
   * is even).  The value reported is the largest in the plateau.
   *
   * @return indices and values of maxima in the power
   */
  public IndexAndValue[] powerPeaks(
          double margin, double threshold, boolean plateauCenterOnly ) {

    return SoundUtils.peaks( power, margin, threshold, plateauCenterOnly );
  }

  /**
   * <p>Find maxima in the smoothed power and return an array with their
   * indices.  Maxima are points that are larger than their neighbors, so
   * if the power has a lot of noise, there may be many.  If there is a
   * flat region in the data, and points on either side are lower, then
   * all points in the flat region are reported as maxima.
   *
   * <p>A noise margin can be specified -- points are not regarded as
   * different unless they differ by more than this amount.  A caution:
   * This is not a substitute for smoothing, and the margin should be
   * small.  No attempt is made to guard against missing a peak whose
   * slopes are shallow enough that adjacent points differ by less than
   * the margin.
   *
   * <p>A threshold can be specified -- peaks at less than this minimum are
   * not recorded.  Note that any margin applies to the test for points
   * being above the threshold, as otherwise the definition of equal within
   * the margin would be violated for points just on either size of the
   * threshold.
   *
   * @param margin Values are considered equal if they differ by no more than
   * this.
   * @param threshold Peaks are not reported unless at least this high.
   * @param plateauCenterOnly If true, only the middle index of a plateau
   * peak is returned (or middle-1 if the number of indices in the plateau
   * is even).  The value reported is the largest in the plateau.
   *
   * @return indices and values of maxima in the smoothed power
   */
  public IndexAndValue[] smoothedPowerPeaks(
          double margin, double threshold, boolean plateauCenterOnly ) {

    return SoundUtils.peaks(
            smoothedPower, margin, threshold, plateauCenterOnly );
  }

  // Things that do actual work

  /**
   * <p>Compute the FFT of a sampled time-domain signal.  It is assumed that
   * samples are taken at equal intervals.
   *
   * <p>The function does not care about the sample rate, so the array
   * positions in the resulting spectrum are abstract units.  The frequency
   * corresponding to offset n in the results array is f = n s / N where
   * s is the sample rate.  That implies one will get frequencies up to the
   * sample rate.  This is a delusion.  If the input is real (as it is here),
   * then the points above N/2 are just the complex conjugates of those
   * below -- the same sample values that fit a sinusoid with frequency
   * n s / N where n <= N/2 also fit a sinusoid with frequency (N-n) s / N.
   * So any results for frequencies above (N/2) s / N = s/2, the Nyquist
   * frequency, should be ignored.  Query methods of FFT will not return
   * values above N/2.
   *
   * <p>Typical use of fft:<br>
   * <pre>
   * // Make an FFT with the desired sample length and sampling rate.
   * FFT fftProcessor = new FFT( N, s );
   * // Acquire a signal.
   * double[] signal = ...
   * // Do the fft.
   * fft( signal );
   * // Perform queries.
   * double avg = fftProcessor.averagePower();
   * int[] peaks = fftProcessor.smoothedPowerPeaks( avg/10, avg );
   * </pre>
   *
   * <p>If you're reading the code, variable names roughly follow the
   * notation of Schilling & Harris, "Fundamentals of Signal Processing",
   * section 3.5, p. 200 in the 2005 edition.
   *
   * @param signal an array of (real) signal amplitudes
   * @throws IllegalArgumentException if the signal array is not the right
   * length
   */
  public void fft( double signal[] ) {

    /*debug*/
    /*
    System.out.println();
    System.out.println( "fft" );
    */

    // Check that the signal is length N...
    if ( signal.length != N )
      throw new IllegalArgumentException(
        "Signal length must match that for which this FFT object was made." );

    // Copy the signal into the real part of the input array, in bit-reversed
    // order.  Actually, we copy into what will be the output array, because
    // our first action in each stage is to swap the arrays.
    for ( int j = 0; j < N; j++ ) {
      results[j].set( signal[bitrev[j]], 0 );
      /*debug*/
      /*
      System.out.println( "in[" + j + "] = " + results[j] );
      */
    }

    // Some scratch variables
    Complex[] temp;   // for swapping the array references
    int s, g, b, n;      // assorted values -- see below
    Complex y = new Complex(0,0);  // holds the second term in each butterfly
    int i, ks, m, mg; // bunch of loop indices -- don't know if declaring
                      // 'em out here lets Java avoid a some work in the
                      // loop...

    // Loop over log N stages.
    for ( i = 1; i <= logN; i++ ) {

      // Swap the arrays to get the results of the previous stage in position
      // as the input for this stage.
      temp = results;
      results = scratch;
      scratch = temp;

      // Stage i has spacing between "groups" in the arrays = s = 2^i.
      s = 1<<i;

      // So the number of groups is N/s.
      g = N/s;

      // Elements of the input array are combined pairwise -- the number of
      // pairs ("butterflies") in each group is s/2.
      b = s/2;
      /*debug*/
      /*
      System.out.println();
      System.out.println( "i = " + i + ", s = " + s + ", g = " + g +
                              ", b = " + b );
      */

      // Loop over groups.  Diverging from S&H, the index increments by the
      // group size each time, to avoid a multiply.  The index at each
      // iteration will be the base index of the group.  The final value of
      // ks will be (g-1)*s = N-s.
      for ( ks = 0; ks < N; ks += s ) {
        /*debug*/
        /*
        System.out.println( "ks = " + ks );
        */

        // Loop over butterflies in this group.  We'll also need the product
        // m*g, so include that as a loop variable.
        for ( m = 0, mg = 0; m < b; m++, mg += g ) {

          // Get the index of the first element in the butterfly, i.e. the
          // one with no multiplier.
          n = ks + m;
          /*debug*/
          /*
          System.out.println( "m = " + m + ", n = " + n +
                              ", n+b = " + (n+b) + ", mg = " + mg );
          */

          // => Note that in the following steps, the recipient of the
          // computations is a pre-existing Complex object -- the contents
          // of the recipient are changed in place.  Note also that
          // non-destructive arithmetic is used, as the contents of the
          // source Complex objects need to be preserved.

          // Compute the second element in the butterfly.  Its multiplier is
          // the unit vector exp( -i2pi mg / N ), i.e. the value at index
          // m*g in our roots table.
          Complex.multiply( scratch[n+b], roots[mg], y );
          /*debug*/
          /*
          System.out.println( "in[n] = " + scratch[n] + ", in[n+b] = " +
                              scratch[n+b] );
          System.out.println( "roots[mg] = " + roots[mg] + ", y = " + y );
          */

          // Combine with the first butterfly element.
          Complex.add( scratch[n], y, results[n] );
          Complex.subtract( scratch[n], y, results[n+b] );
          /*debug*/
          /*
          System.out.println( "out[n] = " + results[n] + ", out[n+b] = " +
                              results[n+b] );
          */
        }
      }
    }

    // When we get here, the results should be in the spectrum array.
    // Before we leave, compute the power density, as that will surely be
    // more commonly used than the raw complex results.  Note we only
    // compute the power out to the Nyquist frequency.
    computePower();

    // And perhaps more useful than the power at each point will be the
    // smoothed power, as peak-finding will be simpler and less subject
    // to noise.
    smooth();
  }

  /** DFT, used for testing. */
  public static Complex[] dft( double[] signal ) {
    int N = signal.length;
    Complex[] results = Complex.makeArray(N);

    // Get (and clean up) roots of unity.  See FFT constructor.
    Complex[] roots = new Complex[N];
    for ( int n = 0; n < N; n++ ) {
      roots[n] = Complex.unitVector( -2*Math.PI*n/N );
    }
    // Fix up the compass points in roots, 0, -pi/2, -pi, -3pi/2.
    roots[0].set( 1, 0 );        // 0
    roots[N/2].set( -1, 0 );     // -pi
    if ( N >= 4 ) {
      roots[N/4].set( 0, -1 );   // -pi/2, same as 3pi/2
      roots[3*N/4].set( 0, 1 );  // -3pi/2, same as pi/2
    }

    // Compute the DFT.  To get the nth output, sum over input (index m),
    // times exp( -i 2pi nm / N ).  Loop over output index:
    Complex w = new Complex(0,0);
    for ( int n = 0; n < N; n++ ) {
      results[n].set(0,0);  // Clear the output
      // Loop over input index:
      for ( int m = 0; m < N; m++ ) {
        // Look up the root of unity in our cleaned up table.  The angle
        // wraps around to 0 at 2pi.  Beware using destructive operations
        // of Complex.  The following statement originally read:
        //   Complex w = roots[ (n*m) % N ];
        // But that copies the reference to the Complex object in the
        // roots array into w.  So in the w.scalarMultiply(...) below,
        // the entry in roots gets changed.  Grrr...
        w.set( roots[ (n*m) % N ] );
        // Complex w = Complex.unitVector( -2*Math.PI*n*m/N );
        results[n].add( w.scalarMultiply( signal[m] ) );
      }
    }
    return results;
  }

  /**
   * Test fft.  The results here have been verified against Matlab's fft.
   * The values obtained by this fft show roundoff error compared to the
   * Matlab results.
   */
  public static void main( String[] args ) {
    // Ask user where to write output files (avoid chooser as this likely
    // will be invoked from Netbeans).
    BufferedReader in
      = new BufferedReader( new InputStreamReader( System.in ) );
    String path = "";
    System.out.println( "Directory for output files " +
                        "(if empty, no files will be written):" );
    try {
      path = in.readLine().trim();
    }
    catch( IOException e ) {
      System.out.println( "readLine threw IOException: " + e.getMessage() );
      System.exit(0);
    }

    // Try a little bitty input -- N = 2.
    int N1 = 2;
    double[] s1 = { 1, 2 };
    System.out.println( "Input array: 1, 2" );

    // With this input, both hand calculation and Matlab's fft get { 3, -1 }.
    Complex[] m1 = { new Complex( 3, 0 ), new Complex( -1, 0 ) };
    System.out.println( "True (by hand & Matlab) fft results: " );
    System.out.println( Complex.arrayString( m1 ) );

    // Compute the DFT.
    Complex[] r1d = dft( s1 );
    // Show the DFT results.
    System.out.println( "DFT results: " );
    System.out.println( Complex.arrayString( r1d ) );

    // Make an FFT and call it.
    // We put in some random choice for the sample frequency -- say, 1 Hz.
    FFT fftLen2 = new FFT(N1,1);
    fftLen2.fft( s1 );
    Complex[] r1f = fftLen2.getResultsArray();
    System.out.println( "FFT results: " );
    System.out.println( Complex.arrayString( r1f ) );
    // Same as true results?
    int diff = 0;
    int far = 0;  // Here, "far" means > 2 ulps.  Should probably be 2 log 2,
                  // i.e. 3ish, to be ~ # operations.
    for ( int n = 0; n < N1; n++ ) {
      if ( !r1f[n].equals( m1[n] ) ) {
        diff++;
        if ( !r1f[n].closeEnough( m1[n], 2 ) ) far++;
      }
    }
    if ( diff != 0 ) {
      System.out.println( "# entries that differ = " + diff );
      if ( far != 0 ) {
        System.out.println( "# entries that differ by > 2 ulps = " + far );
      }
    }
    System.out.println();

    // Next up -- N = 4.

    int N2 = 4;
    double[] s2 = { 1, 2, 3, 4 };
    System.out.println( "Input array: 1, 2, 3, 4" );

    // With this input, both hand calculation and Matlab's fft get
    // { 10, -2+2i, -2, -2-2i }.
    Complex[] m2 = { new Complex( 10, 0 ), new Complex( -2, 2 ),
                     new Complex( -2, 0 ), new Complex( -2, -2 ) };
    System.out.println( "True (by hand & Matlab) fft results: " );
    System.out.println( Complex.arrayString( m2 ) );

    // Compute the DFT.
    Complex[] r2d = dft( s2 );
    // Just for fun, show the results.
    System.out.println( "DFT results: " );
    System.out.println( Complex.arrayString( r2d ) );

    // Make an FFT and call it.
    FFT fftLen4 = new FFT(N2,1);
    fftLen4.fft( s2 );
    Complex[] r2f = fftLen4.getResultsArray();
    System.out.println( "FFT results: " );
    System.out.println( Complex.arrayString( r2f ) );
    // Same?
    diff = 0;
    far = 0;  // This time, we'll take far to be 4 ulps off.
    for ( int n = 0; n < N2; n++ ) {
      if ( !r2f[n].equals( m2[n] ) ) {
        diff++;
        if ( !r2f[n].closeEnough( m2[n], 4 ) ) far++;
      }
    }
    if ( diff != 0 ) {
      System.out.println( "# entries that differ = " + diff );
      if ( far != 0 ) {
        System.out.println( "# entries that differ by > 4 ulps = " + far );
      }
    }
    System.out.println();

    // Now use a CompositeWave to make a known set of pure frequencies.
    // The specific frequencies don't matter -- only their ratio to the
    // sampling rate.

    System.out.println( "Four pure tones" );
    System.out.println();

    double[] freq = { 400, 420, 1000, 7999 };  // Note near Nyquist freq
    int[] amp = { 1, 3, 5, 2 };
    double[] phase = { 0, 1, 2, 3 };  // say anything
    double sampleRate = 16000;
    int nBytes = 2;
    boolean bigEndian = false;
    boolean signed = true;
    
    Wave waves = new CompositeWave( freq, amp, phase, sampleRate, nBytes,
                                    signed, bigEndian );

    // Get a power of 2 worth of samples.  No need to end on period
    // boundaries, and we're actually more interested in how few samples
    // we can get away with using, so don't give it a huge buffer.  We
    // get out half the size we put in, and the range is 0 to the Nyquist
    // frequency, so the spacing between frequency bins is:
    //   bin width = (s/2) / (n/2) = s/n
    // where s is the sample rate and n is the buffer size.  If we want to
    // distinguish our 400Hz and 420Hz signals, we would probably like them
    // to be at least two bins apart.  So:
    //   bin width < 20/2  -> s/n < 10 -> n > 1600
    // Nearest power of 2 above that is 2048.  Ask our Wave for that many
    // samples and unpack the encoding.

    double[] samples = new double[2048];
    waves.getNextRawSamples( samples );
    // Make sure we got something...
    System.out.println( "A few samples:" );
    for ( int i = 0; i < 20; i++ ) {
      System.out.println( i + ": " + samples[i] );
    }
    System.out.println();

    // Make an appropriate FFT and do the transform.
    FFT wavesFFT = new FFT( samples.length, sampleRate );
    wavesFFT.fft( samples );
    Complex[] spectrum = wavesFFT.getResultsArray();
    double[] power = wavesFFT.getPowerArray();
    double[] smoothed = wavesFFT.getSmoothedPowerArray();

    // Show the filter.
    System.out.println( "Filter:" );
    double[] filter = wavesFFT.getFilterArray();
    for( int i = 0; i < filter.length; i++ ) {
      System.out.println( i + ": " + filter[i] );
    }
    System.out.println();

    // Show a few pieces of the result -- complex spectrum, power, and
    // smoothed power around one of the frequencies.
    int i420 = wavesFFT.indexForFrequency( 420 );
    int fMid = (filter.length-1)/2;
    int low = i420 - fMid - 2;  // show a few beyond the smoothing range
    int hi = i420 + fMid + 2;
    low = (low < 0) ? 0 : low;
    hi = (hi > power.length-1) ? power.length-1 : hi;
    System.out.println( "420 Hz is at index " + i420 );
    for( int i = low; i <= hi ; i++ ) {
      System.out.println( i + ": complex: " + spectrum[i] +
                          ", power: " + power[i] +
                          ", smoothed: " + smoothed[i] );
    }
    System.out.println();

    // Find significantly non-zero entries.  Look for 1/1000 of the maximum
    // power.  First, find the max power.
    double max = 0;
    for ( int i = 0; i < power.length; i++ ) {
      if ( power[i] > max ) max = power[i];
    }
    double thousandth = max/1000;
    System.out.println( "Spectrum entries with power > .001 max:" );
    for ( int i = 0; i < power.length; i++ ) {
      if ( power[i] >= thousandth ) {
        System.out.println( i + ": " + spectrum[i] + " Power: " + power[i] );
      }
    }
    System.out.println();

    // Try the peak finders.  Won't really know what to set for the threshold
    // and margin til we have a way to measure the background noise (e.g.
    // histogramming or order statistics).  Choose threshold and margin based
    // on the maximum and/or average.
    // (Would rather use median than average here...  This assumes that the
    // peaks aren't yanking the average up too badly.)
    double avgPower = wavesFFT.averagePower();
    // double maxPower = wavesFFT.maximumPower();
    // double maxSmoothed = wavesFFT.maximumSmoothedPower();
    IndexAndValue[] powerPeaks =
            wavesFFT.powerPeaks( avgPower/10, avgPower, false );
    IndexAndValue[] smoothedPeaks =
            wavesFFT.smoothedPowerPeaks( avgPower/10, avgPower, false );

    System.out.println( "Power peaks, using margin = " + (avgPower/10)
                        + " and threshold = " + avgPower );
    for ( int i = 0; i < powerPeaks.length; i++ ) {
      System.out.println( powerPeaks[i] );
    }
    System.out.println();
    System.out.println( "Smoothed power peaks, using margin = " + (avgPower/10)
                        + " and threshold = " + avgPower );
    for ( int i = 0; i < smoothedPeaks.length; i++ ) {
      System.out.println( smoothedPeaks[i] );
    }
    System.out.println();

    // If the user wants, write out the signal, spectrum, power, and smoothed
    // power in a form appropriate for importing into Matlab.  Show the index
    // of each element and the frequency of the low end of its frequency bin
    // (after the Matlab "..." line continuation).
    if ( path.length() != 0 ) {
      if ( !path.endsWith( "\\" ) ) {
        path = path + "\\";
      }
      try {
        File spectrumFile = new File( path + "waves_spectrum.m" );
        File signalFile = new File( path + "waves_signal.m" );
        File powerFile = new File( path + "waves_power.m" );
        File smoothedFile = new File( path + "waves_smoothed.m" );
        FileWriter spectrumOut = new FileWriter( spectrumFile );
        FileWriter signalOut = new FileWriter( signalFile );
        FileWriter powerOut = new FileWriter( powerFile );
        FileWriter smoothedOut = new FileWriter( smoothedFile );
        spectrumOut.write( "wavesSpectrum = [ " );
        signalOut.write( "wavesSignal = [ " );
        powerOut.write( "wavesPower = [ " );
        smoothedOut.write( "wavesSmoothed = [ " );
        for ( int i = 0; i < power.length; i++ ) {
          String indexAndFrequency =
            Integer.toString(i) + ") " +
            Double.toString( wavesFFT.frequencyAtIndex(i) ) + "Hz";
          spectrumOut.write(
            Double.toString( spectrum[i].getReal() ) + "+" +
            Double.toString( spectrum[i].getImag() ) + "i" + " ... " +
            indexAndFrequency + newline );
          signalOut.write( Double.toString( samples[i] ) + " ... " +
            i + newline );
          powerOut.write( Double.toString( power[i] ) + " ... " +
            indexAndFrequency + newline );
          smoothedOut.write( Double.toString( smoothed[i] ) + " ... " +
            indexAndFrequency + newline );
        }
        // That was only half the signal.  Write out the rest of the spectrum
        // too.
        for ( int i = power.length; i < samples.length; i++ ) {
          String indexAndFrequency =
            Integer.toString(i) + ") " +
            Double.toString( wavesFFT.frequencyAtIndex(i) ) + "Hz";
          spectrumOut.write(
            Double.toString( spectrum[i].getReal() ) + "+" +
            Double.toString( spectrum[i].getImag() ) + "i" + " ... " +
            indexAndFrequency + newline );
          signalOut.write( Double.toString( samples[i] ) + " ... " +
            i + newline );
        }
        spectrumOut.write( "];" + newline );
        signalOut.write( "];" + newline );
        powerOut.write( "];" + newline );
        smoothedOut.write( "];" + newline );
        spectrumOut.close();
        signalOut.close();
        powerOut.close();
        smoothedOut.close();
      }
      catch( IOException e ) {
        System.out.println(
          "Attempt to write to waves*.m failed with IOException:" );
        System.out.println( e.getMessage() );
      }
    }

    // With the default sigma, the smoothed power did not contain the lowest
    // peak.  Try again with sigma = 1.  May need to go to a non-integer
    // sigma to get better control of the filter.  On the other hand, this
    // could just show we need to put the frequencies further apart for this
    // sampling rate.
    System.out.println(
      "Redo the smoothed peak finding with sigma half as large." );
    FFT wavesFFT1 = new FFT( samples.length, sampleRate, 1 );
    wavesFFT1.fft( samples );
    double[] smoothed1 = wavesFFT1.getSmoothedPowerArray();

    // Show the filter.
    System.out.println( "Filter:" );
    double[] filter1 = wavesFFT1.getFilterArray();
    for( int i = 0; i < filter1.length; i++ ) {
      System.out.println( i + ": " + filter1[i] );
    }
    System.out.println();

    double avgPower1 = wavesFFT1.averagePower();
    IndexAndValue[] smoothedPeaks1
            = wavesFFT1.smoothedPowerPeaks( avgPower1/10, avgPower1, false );

    System.out.println( "Smoothed power peaks, using margin = " + (avgPower1/10)
                        + " and threshold = " + avgPower1 );
    for ( int i = 0; i < smoothedPeaks1.length; i++ ) {
      System.out.println( smoothedPeaks1[i] );
    }
    System.out.println();
  }
}
