package ptressel.soundutils;

import java.util.Arrays;

/**
 * Assorted methods that grew up in the other Sound* and related classes,
 * and have been harvested for shared use.  Thus their calling conventions
 * are somewhat idiosyncratic.
 * TODO Make usage and parameter names uniform.
 */

// TODO Currently, except for maxAmplitudePerBitDepth, this and other
// SoundUtils classes only support 8-bit and 16-bit sound depths.  24-bit is
// needed for general use.  Support 24-bit depth throughout.  Find out if
// 32-bit or deeper sound is ever used, if so, support it, if not, fail if
// anyone requests it.

public class SoundUtils {
	
	public static final int SAMPLE_BYTES = 2;

  // Helper classes

  /**
   * Provide a convenience holder for an array index and the value at that
   * index, for methods such as topNPeaks that need to return both the
   * locations and values of items in an array.
   * 
   * These are only holders for data, thus no accessors.
   */
  public static class IndexAndValue implements Comparable<IndexAndValue>
  {
    /** Index of this item in some array or collection in which it was found. */
    public int index;

    /** Value of this item. */
    public double value;

    public IndexAndValue() {
      this.index = 0;
      this.value = 0.0;
    }

    public IndexAndValue( int index, double value ) {
      this.index = index;
      this.value = value;
    }
    
    /**
     * In order to use JUnit Assert#assertArrayEquals on an array of
     * IndexAndValue instances, need to provide an equals() that is true
     * whenever two instances have the same frequency and power.
     * 
     * @param obj object to compare with this IndexAndValue
     * @return true if obj is an IndexAndValue instance with the same
     * frequency and power as this instance
     */
    @Override
    public boolean equals(Object obj) {
      if (obj == null) {
        return false;
      }
      if (getClass() != obj.getClass()) {
        return false;
      }
      final IndexAndValue other = (IndexAndValue) obj;
      if (this.index != other.index) {
        return false;
      }
      if (this.value != other.value) {
        return false;
      }
      return true;
    }

    /**
     * Provide a hashCode to go along with {@link IndexAndValue#equals()},
     * so that equal objects have the same hashCode.
     *
     * @return a hash for this IndexAndValue
     */
    @Override
    public int hashCode() {
      final Double valueObj = new Double(value);
      return index ^ valueObj.hashCode();
    }

    /**
     * Implement Comparable so an array of IndexAndValue elements can be sorted,
     * first, in descending order by value, then within equal values in
     * increasing order of index (as they would be in a stable sort, but
     * include the ordering explicitly).  Order is descending by value since
     * this will be used to pick out the highest peaks, and we'll deliver them
     * with highest peaks first.
     *
     * @param o comparand
     * @return -1 if {@literal this.value > o.value}, or if values are equal
     * but {@literal this.index < o.index}
     * @return 0 if value and index are both equal
     * @return 1 if {@literal this.value < o.value} or values are equal but
     * {@literal this.index > o.index}
     */
    public int compareTo(IndexAndValue o) {
      if ( this.value > o.value ) {
        return -1;
      } else if ( this.value == o.value ) {
        if ( this.index < o.index ) {
          return -1;
        } else if ( this.index == o.index ) {
          return 0;
        } else {
          return 1;
        }
      } else {
        return 1;
      }
    }

    /** For printing, show index: value. */
    @Override
    public String toString() {
      return "" + index + ": " + value;
    }
  }

  // Class data

  // Constants used for audio format endianness manipulations

  /** Mask for a byte at the low-order end of an integer */
  public static final int BYTE_MASK = (1 << Byte.SIZE) - 1;

  /** Mask for two bytes at the low-order end of an integer */
  public static final int TWO_BYTE_MASK = (1 << (2*Byte.SIZE)) - 1;

  /**
   * Provide conversion for byte depth to maximum absolute amplitude of a sound
   * sample of that depth, which is the absolute value of the largest negative
   * number of the given number of bytes.  For n bytes that's 2 to the n-1.
   * Element k of the following array holds the conversion for k+1 bytes.
   * Used by {@link maxAmplitudePerBitDepth(int)}.
   */
  public static final float[] maxAmplitudePerByteDepth =
          { (float) Math.abs( Byte.MIN_VALUE ),
            (float) Math.abs( Short.MIN_VALUE ),
            (float) Math.pow( 2, 23 )
          };

  // Class methods

  /**
   * Provide conversion for bit depth to maximum absolute amplitude of a sound
   * sample of that depth.
   */
  public static float maxAmplitudePerBitDepth( int bits ) {
    if ( bits == 8 || bits == 16 || bits == 24 ) {
      return maxAmplitudePerByteDepth[ ( bits / Byte.SIZE ) - 1 ];
    }
    throw new IllegalArgumentException(
            "Only 8, 16, and 24 bit sound depths are supported." );
  }
  
  /**
   * <p>Helper that extracts a sample from a PCM buffer, obeying the specified
   * sample size, endianness, and signed or not.
   *
   * <p>Only supports sample widths of 1 or 2 bytes.  (What?!  Your sound
   * card can do more than 16 bits?  I don't HEAR you, na na na na...)
   *
   * <p>Requests for 2 byte samples are not required to be aligned on a 2 byte
   * boundary within the buffer.
   *
   * @throws ArrayIndexOutOfBoundsException if requested sample bytes are
   * outside the buffer.
   * @throws IllegalArgumentException if sample width is other than 1 or 2.
   *
   * @param buffer containing samples
   * @param byteOffset of earliest byte in the desired sample
   * @return the reassembled sample value
   */
  public static int extractSample( byte[] buffer, int byteOffset )
      throws IllegalArgumentException { 
  
    int result = 0;
    short temp = 0;
  
    // 1 and 2 bytes are handled separately for efficiency.
    // TODO Check that sample byteOffset is a multiple of the byte size.
    // If performance is an issue, make this an assert; convert other
    // sanity checks to asserts.
    if ( SAMPLE_BYTES == 1 ) {
      // Get the byte -- let it sign-extend.
      result = buffer[byteOffset];
      // Based on signed, clear the high bits or not.
      if ( !signed ) {
        result &= BYTE_MASK;
      }
    } else if ( SAMPLE_BYTES == 2 ) {
      // Get the bytes out in the appropriate order.
      byte low   = 0; 
      byte high  = 0; 
      if ( bigEndian ) {
        // Here, the first byte is the high byte.
        high = buffer[byteOffset];
        low  = buffer[byteOffset+1];
      } else {
        // Here, the second byte is the high byte.
        high = buffer[byteOffset+1];
        low  = buffer[byteOffset];
      }
      // Copy in the low byte let it sign-extend, then clear all but that byte.
      result = low;
      result &= BYTE_MASK;
      // Include the high byte.  Yup, this sign-extends if the high bit of
      // high is on...
      result |= high << Byte.SIZE;
      // Low two bytes are in position.  Nuke the sign if we don't want it.
      if ( !signed ) {
        result &= TWO_BYTE_MASK;
      }
    } else {
      throw new IllegalArgumentException(
        "SoundUtils.extractSample: SAMPLE_BYTES must be 1 or 2, got "
        + SAMPLE_BYTES );
    }

    return result;
  }

  /**
   * <p>Convert the given sample value (assumed signed and scaled as intended)
   * according to the given audio format values, and store it in successive
   * bytes starting at byteOffset in the supplied buffer.  If the value is out of
   * range for the audio format, it will be clipped to the limiting value.
   *
   * <p>Currently only supports signed, mono, 1 or 2 bytes.
   *
   * @param sampleValue the value to format and store
   * @param buffer the array to store it in
   * @param byteOffset the lowest index byte to put it in
   * @param SAMPLE_BYTES
   * @param bigEndian true if sampleValue is big endian.
   * @param signed true if sampleValue is signed
   *
   * @throws ArrayIndexOutOfBoundsException if byteOffset + SAMPLE_BYTES extends
   * past the end of the buffer
   * @throws IllegalArgumentException if the audio format is not supported
   */
  public static void insertSample( double sampleValue, byte[] buffer,
                                   int byteOffset, int SAMPLE_BYTES,
                                   boolean signed, boolean bigEndian )
      throws ArrayIndexOutOfBoundsException, IllegalArgumentException {

    /*debug*/
    /*
    System.out.println( "byteOffset = " + byteOffset );
    */

    // Assume the value is in range -- get it as an integer.
    long integerValue = Math.round( sampleValue );
    /*debug*/
    /*
    System.out.println( "insertSample: sampleValue = " + sampleValue +
                        ", integerValue = " + integerValue );
    */

    switch ( SAMPLE_BYTES ) {

      case 1:
        // One byte is the easy case -- clip to fit the allowed range and get
        // the low byte.
        if ( sampleValue > Byte.MAX_VALUE )
          integerValue = Byte.MAX_VALUE;
        else if ( sampleValue < Byte.MIN_VALUE )
          integerValue = Byte.MIN_VALUE;
        /*debug*/
        /*
        System.out.println( "insertSample: after clipping = " + integerValue +
          ", in hex = " + Long.toHexString( integerValue ) );
        */
        buffer[byteOffset] = (byte) integerValue;
        /*debug*/
        /*
        System.out.println( "insertSample: stored value = " + buffer[byteOffset] +
          ", in hex = " + Integer.toHexString( BYTE_MASK & buffer[byteOffset] ) );
        */
        break;

      case 2:
        // For two bytes, need to split the bytes and order them according
        // to endianness.
        if ( sampleValue > Short.MAX_VALUE )
          integerValue = Short.MAX_VALUE;
        else if ( sampleValue < Short.MIN_VALUE )
          integerValue = Short.MIN_VALUE;
        /*debug*/
        /*
        System.out.println( "insertSample: after clipping = " + integerValue +
          ", in hex = " + Long.toHexString( integerValue ) );
        */
        // Separate the high and low bytes.
        byte low = (byte) integerValue;
        byte high = (byte) ( integerValue >> Byte.SIZE );
        /*debug*/
        /*
        System.out.println( "insertSample: low = " + low + " (hex " +
          Integer.toHexString(BYTE_MASK&low) + "), high = " + high +
          " (hex " + Integer.toHexString(BYTE_MASK&high) + ")" );
        */
        if ( bigEndian ) {
          // Here, the high-order byte goes first.
          buffer[byteOffset] = high;
          buffer[byteOffset+1] = low;
          /*debug*/
          /*
          System.out.println( "insertSample: big endian, 0 byte = " +
            buffer[byteOffset] + " (hex " +
            Integer.toHexString(BYTE_MASK&buffer[byteOffset]) +
            "), 1 byte = " + buffer[byteOffset+1] + " (hex " +
            Integer.toHexString(BYTE_MASK&buffer[byteOffset+1]) + ")" );
          */
        } else {
          // Low byte goes first.
          buffer[byteOffset] = low;
          buffer[byteOffset+1] = high;
          /*debug*/
          /*
          System.out.println( "insertSample: little endian, 0 byte = " +
            buffer[byteOffset] + " (hex " +
            Integer.toHexString(BYTE_MASK&buffer[byteOffset]) +
            "), 1 byte = " + buffer[byteOffset+1] +
            Integer.toHexString(BYTE_MASK&buffer[byteOffset+1]) + ")" );
          */
        }
        break;

      default:
        throw new IllegalArgumentException(
          "Number of bytes per sample must be 1 or 2." );
    }
  }

  /**
   * <p>Find maxima in the given array and return an array with their
   * indices.  Maxima are points that are larger than their neighbors, so
   * if the power has a lot of noise, there may be many.  If there is a
   * flat region in the data, and points on either side are lower, then
   * all points in the flat region are reported as maxima.  Points at the
   * ends of the array are regarded as higher than their non-existent
   * neighbors outside the array.
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
   * @param array in which to find maxima
   * @param margin Values are considered equal if they differ by no more than
   * this.
   * @param threshold Peaks are not reported unless at least this high.
   * @param plateauCenterOnly If true, only the middle index of a plateau
   * peak is returned (or middle-1 if the number of indices in the plateau
   * is even).  The value reported is the largest in the plateau.
   * 
   * @return indices and values of maxima in the array
   */
  // TODO Provide a generic that allows T[] array, for use with Complex
  // array elements.  (Note Complex will need a compareTo that compares
  // amplitudes.  That can't be reconciled with its equals().)
  public static IndexAndValue[] peaks(
          double[] array, double margin, double threshold,
          boolean plateauCenterOnly ) {
   
    // We can't find more peaks than the entire array (and that would
    // only happen if it were flat).  We'll prune unused entries later.
    IndexAndValue[] tempPeaks = new IndexAndValue[array.length];
    int numPeaks = 0;   // How far along we are in tempPeaks.
    boolean up = true;  // Remember if we've seen a rising edge (true) or a
                        // falling edge (false) most recently before a run of
                        // equal values.  We pretend we saw a rising edge at
                        // the outset.
    int plateau = -1;   // Remember the first point in a plateau (run of values
                        // in our list of peaks that are within margin of each
                        // other, but are above the points that surround them).
                        // We won't know if the flat region is a peak or just a
                        // non-strictly monotonic slope til we get to its far
                        // edge, so we may need to back these points out of the
                        // list of peaks.
    double plateauMaxValue = Double.NEGATIVE_INFINITY;
                        // Highest value seen in the plateau.

    margin = Math.abs(margin);  // Beware silly callers.
    
    for ( int i = 0; i < array.length-1; i++ ) {
                        
      // Split into greater, equal, less than next neighbor
      double diff = array[i] - array[i+1];
     
      if ( diff > margin ) {

        // Current point is above next.  We record it as a peak unless
        // we're already on a slope headed down, or it's too small.
        if ( up && ( array[i] > threshold - margin ) ) {
          // We passed a rising edge previously, and are about to head down --
          // this point is a peak.
          tempPeaks[ numPeaks++ ] = new IndexAndValue( i, array[i] );
          // If we were in a plateau, this point is still part of it.
          if ( ( plateau > 0 ) && ( array[i] > plateauMaxValue ) ) {
            plateauMaxValue = array[i];
          }
        }

        // If we were on a plateau, we've reached its end -- do we want all of
        // it or just the center?
        if ( plateau >= 0 && plateauCenterOnly ) {
          // Just the center -- overwrite the first plateau entry with the
          // center index in the plateau and the max value.  Note we have
          // already created an IndexAndValue there.  Note numPeaks points
          // to the peaks entry after the end of the plateau.
          tempPeaks[plateau].index =
            ( tempPeaks[plateau].index + tempPeaks[numPeaks-1].index ) / 2;
          tempPeaks[plateau].value = plateauMaxValue;
          // Reset the next peak index so we'll fill in the next peak after
          // the plateau center.
          // TODO perhaps.  This will cause some already-made IndexAndValue
          // instances to be abandoned, but may not be worth checking for
          // and reusing them.
          numPeaks = plateau + 1;
        }

        // End any plateau check.
        plateau = -1;
        plateauMaxValue = Double.NEGATIVE_INFINITY;

        // When we advance to the next point, we will have seen a falling edge.
        up = false;

      } else if ( diff < -margin ) {

        // Current point is below next, so isn't a peak.  We are (now) on a
        // rising edge.  If we had provisional peaks on a plateau, we need
        // to cancel them.
        if ( plateau >= 0 ) {
          // Yes, we were in a plateau.  The first point was at index
          // plateau in tempPeaks.  Reset tempPeak's next recording location
          // to that so we overwrite the plateau points.  Clear the plateau
          // check.
          numPeaks = plateau;
          plateau = -1;
          plateauMaxValue = Double.NEGATIVE_INFINITY;
        }

        // Remember the rising edge.
        up = true;

      } else {

        // Current point is (approximately) equal to next.
        // If we saw a rising edge most recently, then we record this as a
        // tentative peak so long as it's big enough.
        if ( up && array[i] > threshold - margin ) {
          tempPeaks[ numPeaks ] = new IndexAndValue( i, array[i] );
          // If we aren't already in a plateau, record the start of one.
          if ( plateau == -1 ) {
            plateau = numPeaks;
          }
          if ( array[i] > plateauMaxValue ) {
            plateauMaxValue = array[i];
          }
          numPeaks++;
        }
      }
    }

    // We stopped at the point before the end -- check the end.  We consider
    // the point after it to be lower, so we are in the "current point is
    // above next" case -- this is a peak if we saw a rising edge most
    // recently.
    if ( up ) {
      tempPeaks[ numPeaks++ ] =
              new IndexAndValue( array.length-1, array[ array.length-1 ] );
      // Handle end of plateau (see diff > margin case above).
      if ( plateau > 0 && array[ array.length-1 ] > plateauMaxValue ) {
        plateauMaxValue = array[ array.length-1 ];
      }
      if ( plateau >= 0 && plateauCenterOnly ) {
        tempPeaks[plateau].index =
          ( tempPeaks[plateau].index + tempPeaks[numPeaks-1].index ) / 2;
        tempPeaks[plateau].value = plateauMaxValue;
        numPeaks = plateau + 1;
      }
    }

    // If we didn't fill tempPeaks (hah!), copy it into an array of the right
    // size.
    IndexAndValue[] peaks = new IndexAndValue[numPeaks];
    System.arraycopy( tempPeaks, 0, peaks, 0, numPeaks );

    return peaks;
  }

  /**
   * Return the indices and values of (at most) the highest N peaks with the
   * supplied array, sorted in decreasing order by value.  If there are fewer
   * fewer than N peaks, the returned array will only contain as many entries
   * as there are peaks.  The sort used is stable -- peaks of the same
   * amplitude will be returned in order of their index in the original array.
   * 
   * For noisy data, it is recommended to smooth the data or average over
   * mulitple samples.
   * 
   * See method {@link #peak} for how to use the margin, threshold, and
   * plateauCenterOnly parameters.
   *
   * @param N maximum number of peaks to return
   * @param array in which to find peaks
   * @param margin values are considered equal if they differ by no more than
   * this
   * @param threshold peaks are not reported unless at least this high
   * @param plateauCenterOnly If true, only the middle index of a plateau
   * peak is returned (or middle-1 if the number of indices in the plateau
   * is even).  The value reported is the largest in the plateau.
   *
   * @return indices and values of maxima in the array
   */
  public static IndexAndValue[] topNPeaks( int N, double[] array,
                                           double margin, double threshold,
                                           boolean plateauCenterOnly ) {
    // TODO NEXT
    // Call peaks(), sort result, truncate to N.
    IndexAndValue allPeaks[] =
            peaks( array, margin, threshold, plateauCenterOnly );
    Arrays.sort(allPeaks);
    return Arrays.copyOf( allPeaks, Math.min( N, allPeaks.length ) );
  }

  /**
   * <p>Make one cycle worth of pure tone values, starting at maximum amplitude
   * (like a cosine), with the given amplitude, frequency, and sample rate.
   * The resulting double array of values is unformatted and unclipped.
   * (This is appropriate for combining into a composite waveform, where any
   * individual component should not be individually clipped even if it would
   * exceed the audio format's range, as it may combine destructively with
   * another tone, so that the sum is in range.  Clipping first would give the
   * wrong sum, i.e. would differ from what would happen to physical sound
   * waves.)
   *
   * @param frequency frequency of the signal to be generated
   * @param amplitude amplitude of the signal to be generated
   * @param sampleRate audio format: sample rate
   *
   * @return samples values representing one cycle of the signal
   */
  // TODO Here and elsewhere:  Should cyles start with zero amplitude (sine)?
  // Might be better for avoiding a pop at the beginning of playing the sound,
  // but not so much use later so long as contiguous samples are extracted.
  public static double[] generateCycle( double frequency, int amplitude,
                                        double sampleRate ) {

    // Make a buffer for one cycle of the wave.  The number of samples in
    // one cycle is the sample frequency / wave frequency.
    int samplesPerCycle = (int) (sampleRate / frequency);
    double[] cycle = new double[ samplesPerCycle ];

    // Generate samples out of a cosine for one cycle.
    for ( int i = 0; i < samplesPerCycle; i++ ) {

      // Compute the sample value for sample i.  This is:
      //   2 pi i / samplesPerCycle
      // (Don't precompute 2 * Math.PI / samplesPerCycle as that may
      // compound any roundoff error.  We only need to do this setup once...)
      cycle[i] = amplitude * Math.cos( 2 * i * Math.PI / samplesPerCycle );
    }

    return cycle;
  }

  /**
   * <p>Make one cycle worth of samples, starting at phase zero (i.e. a cosine),
   * with the given amplitude, frequency, and sample rate.  Encode them in a
   * short array according to the given audio format.  Clip to fit the sample
   * size if necessary.
   *
   * <p>The returned array of samples can be copied repeatedly to fill a
   * buffer for output to a playback device configured with the given audio
   * format.
   *
   * <p>Supported audio formats are mono, 2 bytes per sample, signed, PCM.
   *
   * @param frequency frequency of the signal to be generated
   * @param amplitude amplitude of the signal to be generated
   *
   * @param sampleRate audio format: sample rate
   * @param channels audio format: number of channels (only 1 supported)
   *
   * @return encoded samples representing one cycle of the signal
   * @throws IllegalArgumentException if specified audio format is not
   * supported
   */
  // TODO Here and elsewhere:  Pass in an AudioFormatParams instead of
  // individual signed, endian, bytes, #channels values.
  public static short[] generateCycleFormatted(
          double frequency, int amplitude,
          double sampleRate, int channels ) {

    // Get the raw (unformatted and unclipped) values.
    double[] rawCycle = generateCycle( frequency, amplitude, sampleRate );

    // Make a buffer for the formatted samples.
    short[] cycle = new short[ rawCycle.length ];

    // Get the amplitude limits for the various sample sizes.  For 2s
    // complement, max is ( 2^#bits / 2 ) -1 , and min is - ( 2^#bits / 2 ).
    int bits = Byte.SIZE * SAMPLE_BYTES;
    long max = 1 << (bits-1);  // Hmm, will Java do long arithmetic if need be?
    long min = -max;
    --max;
   
    // Clip and format the values.
    for ( int i = 0; i < rawCycle.length; i++ ) {

      double value = rawCycle[i];

      // Is that out of range?  If so, constrain it to the legal range for
      // the given sample size.
      if ( value > max ) value = max;
      else if ( value < min ) value = min;

      // Deal with the audio format and store this in the cycle buffer.
      // We've already done an amplitude check.  Note insertSample wants the
      // byte offset to store the result at, not the sample number.
      insertSample( value, cycle, i*SAMPLE_BYTES, SAMPLE_BYTES, signed,
                    bigEndian );
      /*debug*/
      /*
      String s0 = toHexString( cycle[i] );
      System.out.print( "PureWave: " + i + ": before format = " + value +
        ", after format = " + cycle[i] + " (hex " + s0 + ")" );
      if ( SAMPLE_BYTES > 1 ) {
        String s1 = toHexString( cycle[i+1] );
        System.out.print( ", " + cycle[i+1] + " (hex " + s1 + ")" );
      }
      System.out.println();
      */
    }

    return cycle;
  }

  /**
   * Helper for debug messages that want to print bytes in hex.  Unfortunately
   * the Byte class doesn't provide toHexString, so we make one here.
   * (Integer#toHexString sign-extends, so a negative byte gets left-filled
   * with 6 unnecessary "f"s.)
   *
   * @param b byte to format as hex
   * @return the two-hex-digit String representation of b
   */
  public static String toHexString( byte b ) {
    // Don't let the silly thing sign-extend...
    // TODO Should this left-fill positive numbers with zeros, so all bytes
    // are shown as two chars?  Integer#toHexString doesn't.
    return Integer.toHexString( b & BYTE_MASK );
  }

  /**
   * Helper for debug messages that want to print shorts in hex.  Unfortunately
   * the Short class doesn't provide toHexString, so we make one here.
   *
   * @param s short to format as hex
   * @return the four-hex-digit String representation of s
   */
  public static String toHexString( short s ) {
    // Don't let the silly thing sign-extend...
    // TODO Should this left-fill positive numbers with zeros, so all shorts
    // are shown as four chars?  Integer#toHexString doesn't.
    return Integer.toHexString( s & TWO_BYTE_MASK );
  }

  /**
   * <p>Helper for command line options with arguments.  Finds the specified
   * option string among the given args array and returns the following
   * argument, if any.
   *
   * <p>If there are multiple instances, this returns the first.
   *
   * <p>Official POSIX options with arguments have the form -<name> <arg>
   * where <name> is a string not containing - and <arg> is a string, quoted
   * if it contains whitespace.  We don't place any restriction on the
   * format of the option, but assume the argument is a single string and
   * that it does not start with -.  Quoting a string starting with - does
   * not protect it from being regarded as a new option, not an argument.
   *
   * @param args the array of String arguments, usually from main
   * @param option the string indicating the option.
   * @return the argument string if there was one, null if no matching option
   * was found, or the empty string if the following item is an option or
   * there is no following item.
   */
  public static String getArgForOption( String[] args, String option ) {

    // Look for it.
    for ( int i = 0; i < args.length; i++ ) {
      if ( option.equals( args[i] ) ) {
        // Found the option.  Is there an arg?
        if ( args.length > i+1 && args[i+1].charAt(0) != '-' ) {
          // There's another token on the line and it doesn't start with -,
          // ergo it's our arg.
          return args[i+1];
        } else {
          // Nothing else on the line or another option follows -- return empty.
          return "";
        }
      }
    }

    return null;  // Didn't find the option
  }

  /**
   * <p>Helper for command line options without arguments.  Checks for the
   * specified option string among the given args array and returns true
   * if present.
   *
   * @param args tokenized arguments
   * @param option the option string
   * @return true if option is present
   */
  public static boolean checkOption( String[] args, String option ) {

    // Look for it.
    for ( int i = 0; i < args.length; i++ ) {
      if ( option.equals( args[i] ) ) return true; // Found the option.
    }
    return false;  // Didn't find the option
  }

  /** Helper for main:  run one set of tests on extractSample. */
  private static boolean oneTest( byte[] input, int[] results, int SAMPLE_BYTES,
                                 boolean signed, boolean bigEndian ) {

    boolean allOk = true;  // Assume the best.

    // Loop over the results.  Index into input is result index * SAMPLE_BYTES.
    for ( int r = 0, i = 0; r < results.length; r++, i += SAMPLE_BYTES ) {
      // Call extractSample.
      int convSays = extractSample( input, i, SAMPLE_BYTES );
      // Compare with expected result.
      if ( convSays != results[r] ) {
        System.out.println( "Input: " + input[i]
          + (SAMPLE_BYTES == 2 ? (", " + input[i+1]) : "")
          + " Expected: " + results[r] + " Got: " + convSays );
        allOk = false;
      }
    }

    return allOk;
  }

  /** Tests */
  public static void main( String[] args ) {

    // Test extractSample.
    System.out.println();
    System.out.println( "Test extractSample: " );
    System.out.println();

    // Input for all tests.
    byte[] input = { 1, -1, -2, 2, 0, -3, -4, 0 };

    // Test all combinations of # bytes, endianness, signed.

    // 1, big, signed
    int[] result1 = { 1, -1, -2, 2, 0, -3, -4, 0 };
    System.out.println( "1 byte, signed, big" );
    if ( oneTest( input, result1, 1, true, true ) ) {
      System.out.println( "Ok." );
    }

    // 1, signed, little
    int[] result2 = { 1, 0xFF, 0xFE, 2, 0, 0xFD, 0xFC, 0 };
    System.out.println( "1 byte, signed, little" );
    if ( oneTest( input, result2, 1, true, false ) ) {
      System.out.println( "Ok." );
    }

    // Should get same results for little as big if 1 byte.
    // 1, unsigned, big
    int[] result3 = result1;
    System.out.println( "1 byte, unsigned, big" );
    if ( oneTest( input, result3, 1, false, true ) ) {
      System.out.println( "Ok." );
    }

    // 1, unsigned, little
    int[] result4 = result2;
    System.out.println( "1 byte, unsigned, little" );
    if ( oneTest( input, result4, 1, false, false ) ) {
      System.out.println( "Ok." );
    }

    // 2, signed, big
    int[] result5 = { 0x1FF, 0xFFFFFE02, 0xFD, 0xFFFFFC00 };
    System.out.println( "2 bytes, signed, big" );
    if ( oneTest( input, result5, 2, true, true ) ) {
      System.out.println( "Ok." );
    }

    // 2, signed, little
    int[] result6 = { 0x1FF, 0xFE02, 0xFD, 0xFC00 };
    System.out.println( "2 bytes, signed, little" );
    if ( oneTest( input, result6, 2, true, false ) ) {
      System.out.println( "Ok." );
    }

    // 2, unsigned, bi
    int[] result7 = { 0xFFFFFF01, 0x2FE, 0xFFFFFD00, 0xFC };
    System.out.println( "2 bytes, unsigned, big" );
    if ( oneTest( input, result7, 2, false, true ) ) {
      System.out.println( "Ok." );
    }

    // 2, unsigned, little
    int[] result8 = { 0xFF01, 0x2FE, 0xFD00, 0xFC };
    System.out.println( "2 bytes, unsigned, little" );
    if ( oneTest( input, result8, 2, false, false ) ) {
      System.out.println( "Ok." );
    }

    // insertSample was tested during testing of PureWave (its former home)
    // but should have tests added here.

    // Show output of toHexString(byte).
    System.out.println();
    System.out.println( "Test toHexString: " );
    System.out.println();

    // What's the issue?  Only want to see the bits in the actual byte, not
    // what it gets promoted to when used as an arg to Integer.toHexString.
    // Negative values get sign-extended and end up with a bunch of extraneous
    // Fs.  E.g.:
    //   int x1 = 0xEE;
    // is fine as an arg to Integer.toHexString -- the literal is already an
    // int value, and is positive.
    //   byte b2 = 0xEE;
    //   int x2 = b2;
    // Whoops -- b2 is negative, so the assignment to x2 sign-extends.

    // negative inputs
    byte b = (byte) 0xEE;
    System.out.println( "With byte b = 0xEE:" );
    System.out.println( "Integer.toHexString( b ) = " +
                        Integer.toHexString( b ) );
    System.out.println( "SoundUtils.toHexString( b ) = " +
                        SoundUtils.toHexString( b ) );
    short s = (short) 0xEEEE;
    System.out.println( "With short s = 0xEEEE:" );
    System.out.println( "Integer.toHexString( s ) = " +
                        Integer.toHexString( s ) );
    System.out.println( "SoundUtils.toHexString( s ) = " +
                        SoundUtils.toHexString( s ) );

    // positive inputs
    b = 1;
    System.out.println( "With byte b = 1:" );
    System.out.println( "Integer.toHexString( b ) = " +
                        Integer.toHexString( b ) );
    System.out.println( "SoundUtils.toHexString( b ) = " +
                        SoundUtils.toHexString( b ) );
    System.out.println( "With short s = 1:" );
    System.out.println( "Integer.toHexString( s ) = " +
                        Integer.toHexString( s ) );
    System.out.println( "SoundUtils.toHexString( s ) = " +
                        SoundUtils.toHexString( s ) );

    // Show samples produced by generateCycle.
    System.out.println();
    System.out.println( "Test generateCycle: " );
    System.out.println();

    // Choose params so there are only a few samples through the wave.
    // Amplitude corresponds to the 2 byte clipped test for formatted
    // samples, below.
    double[] rawCycle
      = generateCycle( 400, Short.MAX_VALUE+Short.MAX_VALUE/3, 8*400 );

    // Show what we got.  Since we asked for 1 byte samples, we don't
    // need any byte swapping or repacking.
    System.out.println( "Unformatted samples:" );
    for ( int i = 0; i < rawCycle.length; i++ ) {
      System.out.println( i + ": " + rawCycle[i] );
    }
    System.out.println();

    // Show samples produced by generateCycleFormatted.
    System.out.println();
    System.out.println( "Test generateCycleFormatted: " );
    System.out.println();

    // First 1 byte, no clipping.  Choose params so there are only a few
    // samples through the wave.
    byte[] cycle
      = generateCycleFormatted( 400, Byte.MAX_VALUE/2, 8*400, 1,
                                true, true, 1 );

    // Show what we got.  Since we asked for 1 byte samples, we don't
    // need any byte swapping or repacking.
    System.out.println( "1-byte samples:" );
    for ( int i = 0; i < cycle.length; i++ ) {
      System.out.println( i + ": " + cycle[i] );
    }
    System.out.println();

    // Repeat for two-byte samples.  Again, no clipping.
    cycle = generateCycleFormatted( 400, Short.MAX_VALUE/2, 8*400, 2,
                                    true, true, 1 );

    System.out.println( "2-byte samples:" );
    for ( int i = 0; i < cycle.length; i += 2 ) {
      String s0 = toHexString( cycle[i] );
      String s1 = toHexString( cycle[i+1] );
      System.out.println( i + ": " +
        cycle[i] + " (hex " + s0 + "), " +
        cycle[i+1] + " (hex " + s1 + ")" );
    }
    System.out.println();


    // Repeat for two-byte samples, this time, with amplitude over the limit.
    cycle = generateCycleFormatted(400, Short.MAX_VALUE+Short.MAX_VALUE/3,
                                   8*400, 2, true, true, 1);

    System.out.println( "2-byte samples, clipped:" );
    for ( int i = 0; i < cycle.length; i += 2 ) {
      String s0 = toHexString( cycle[i] );
      String s1 = toHexString( cycle[i+1] );
      System.out.println( i + ": " +
        cycle[i] + " (hex " + s0 + "), " +
        cycle[i+1] + " (hex " + s1 + ")" );
    }
    System.out.println();

    // Try out the peak finder.  Make some little arrays with all the cases.
    System.out.println( "Checking peak finder." );

    // With a margin of 2, this one has peaks only at the ends.  It exercises
    // plateau rejection, margin, and peak at ends cases.  Use margin = 2,
    // threshold = 0;
    double[] array1 = { 100, 90, 80, 82, 80, 78, 80, 70, 80, 90, 90, 90, 100 };
    IndexAndValue[] peaks1 = {
      new IndexAndValue( 0, 100 ),
      new IndexAndValue( array1.length-1, 100 )  };
    IndexAndValue[] reported1 = peaks( array1, 2, 0, false );
    // Compare peaks.
    boolean peaksOk1 = false;
    if ( peaks1.length == reported1.length ) {
      for ( int i = 0; i < peaks1.length; i++ ) {
        if ( peaks1[i] != reported1[i] ) break;  // Wah.
      }
      peaksOk1 = true;  // Whew.
    }
    if ( !peaksOk1 ) {
      System.out.println( "1st set of peaks should have been:" );
      for ( int i = 0; i < peaks1.length; i++ ) {
        System.out.print( peaks1[i] + " " );
      }
      System.out.println();
      System.out.println( "Got instead:" );
      for ( int i = 0; i < reported1.length; i++ ) {
        System.out.print( reported1[i] + " " );
      }
    } else {
      System.out.println( "1st set of peaks ok." );
    }
    System.out.println();

    // With a margin of 2 and threshold of 25, this has no peaks at the ends,
    // two plateau peaks, and a single point peak.  It tests ignoring differences
    // within margin and ignoring values below threshold.
    double[] array2 = {
      80, 90, 92, 90, 80, 10, 12, 20, 10, 90, 88, 90, 80, 83, 80 };
    IndexAndValue[] peaks2 = {
      new IndexAndValue( 1, 90 ),
      new IndexAndValue( 2, 92 ),
      new IndexAndValue( 3, 90 ),
      new IndexAndValue( 9, 90 ),
      new IndexAndValue( 10, 88 ),
      new IndexAndValue( 11, 90 ),
      new IndexAndValue( 13, 83 ),
    };
    IndexAndValue[] reported2 = peaks( array2, 2, 25, false );
    // Compare peaks.
    boolean peaksOk2 = false;
    if ( peaks2.length == reported2.length ) {
      for ( int i = 0; i < peaks2.length; i++ ) {
        if ( peaks2[i] != reported2[i] ) break;  // Wah.
      }
      peaksOk2 = true;  // Whew.
    }
    if ( !peaksOk2 ) {
      System.out.println( "2nd set of peaks should have been:" );
      for ( int i = 0; i < peaks2.length; i++ ) {
        System.out.print( peaks2[i] + " " );
      }
      System.out.println();
      System.out.println( "Got instead:" );
      for ( int i = 0; i < reported2.length; i++ ) {
        System.out.print( reported2[i] + " " );
      }
    } else {
      System.out.println( "2nd set of peaks ok." );
    }
    System.out.println();

    // Using the same array of values, find peaks but ask for only centers of
    // plateaus.
    IndexAndValue[] peaks2a = {
      new IndexAndValue( 2, 92 ),
      new IndexAndValue( 10, 90 ),
      new IndexAndValue( 13, 83 ),
    };
    IndexAndValue[] reported2a = peaks( array2, 2, 25, true );
    // Compare peaks.
    boolean peaksOk2a = false;
    if ( peaks2a.length == reported2a.length ) {
      for ( int i = 0; i < peaks2a.length; i++ ) {
        if ( peaks2a[i] != reported2a[i] ) break;  // Wah.
      }
      peaksOk2a = true;  // Whew.
    }
    if ( !peaksOk2a ) {
      System.out.println( "Peaks with plateau centers only should have been:" );
      for ( int i = 0; i < peaks2a.length; i++ ) {
        System.out.print( peaks2a[i] + " " );
      }
      System.out.println();
      System.out.println( "Got instead:" );
      for ( int i = 0; i < reported2a.length; i++ ) {
        System.out.print( reported2a[i] + " " );
      }
    } else {
      System.out.println( "Peaks with plateau centers only ok." );
    }
    
    // Show hash code for an IndexAndValue, and values it's constructed from.
    IndexAndValue iv = new IndexAndValue( 10, 100.0 );
    long longValueOfValue = Double.doubleToLongBits( iv.value );
    int hashCodeOfValue = new Double( iv.value ).hashCode();
    System.out.println(
            "IndexAndValue(" + iv.index + "," + iv.value + ").hashCode() = " +
            iv.hashCode() );
    System.out.println(
            "Double.doubleToLongBits(" + iv.value + ") = " + longValueOfValue );
    System.out.println(
            "Double(" + iv.value + ").hashCode() = " + hashCodeOfValue );
    System.out.println();
  }
}
