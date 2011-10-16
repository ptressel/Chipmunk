package ptressel.soundutils;

/**
 * <p>Encapsulates a set of frequency, amplitude, and phase parameters for a
 * CompositeWave.  These can be provided to a CompositeWave constructor, or
 * can represent a chord-based symbol.  A set of CompositeWaveParams is
 * provided to Receive to serve as the state symbols for StateRotation.
 *
 * <p>CompositeWaveParams does not include audio format parameters.  The
 * constructors that accept CompositeWaveParams, if they need audio format
 * info, should accept AudioFormatParams.
 *
 * <p>CompositeWaveParams includes arrays for frequency, amplitude, and
 * initial phase for each pure tone included in the chord.  See CompositeWave
 * for a description of these values.
 */
public class CompositeWaveParams {

  /** The system's newline, for toString. */
  public static String newline = System.getProperty( "line.separator" );

  /** Number of pure tones in this chord */
  private int nWaves;

  /** Frequency */
  private double[] frequency;

  /** Initial phase in radians */
  private double[] initialPhase;

  /** Amplitude */
  private int[] amplitude;

  /**
   * <p>Make a CompositeWaveParams to hold the supplied frequency, phase, and
   * amplitude information.  See CompositeWave for a description of these
   * values.
   *
   * <p>The initial phase argument may be left null, in which case the initial
   * phase will be taken to be zero.  Initial phase is not useful unless one
   * is doing phase shift keying.
   *
   * <p>The contents of the caller's arrays are copied: CompositeWaveParams
   * does not store references to the caller's arrays.  The caller may change
   * their arrays later without worry.  However, the accessors return
   * references to the actual CompositeWaveParams internal data, not copies.
   * The caller is trusted not to damage these.
   *
   * @param frequency frequencies for each wave in the symbol
   * @param amplitude the maximum amplitude for each wave in the symbol
   * @param initialPhase for each wave, the phase in radians that will appear
   * at the outset of the symbol
   * @throws IllegalArgumentException if the numbers of frequencies,
   * amplitudes, and initial phases (if supplied) are not the same
   */
  public CompositeWaveParams( double[] frequency, int[] amplitude,
           double[] initialPhase ) throws IllegalArgumentException {

    // Sanity check -- all arrays must be same length.
    nWaves = frequency.length;
    if ( amplitude.length != nWaves ||
         ( initialPhase != null && initialPhase.length != nWaves ) ) {
      throw new IllegalArgumentException(
        "Same number of values must be supplied for each wave parameter." );
    }

    // Store the user's waveform values.  We make our own copies of the
    // arrays, so the user can do what they want with their arrays.
    this.frequency = frequency.clone();
    this.amplitude = amplitude.clone();
    // If the user didn't supply initial phases, make some.
    if ( initialPhase != null ) {
      this.initialPhase = initialPhase.clone();
    } else {
      this.initialPhase = new double[nWaves];
      for ( int i = 0; i < nWaves; i++ ) {
        this.initialPhase[i] = 0;
      }
    }
  }

  /** Get the frequency array. */
  public double[] getFrequency() { return frequency; }

  /** Get the amplitude array. */
  public int[] getAmplitude() { return amplitude; }

  /** Get the initial phase array. */
  public double[] getInitialPhase() { return initialPhase; }

  /** Provide a string describing the waves in the CompositeWaveParams. */
  public String toString() {
    StringBuffer s = new StringBuffer( "Composite wave:" + newline );
    // Add each wave's info on a line.
    for ( int i = 0; i < nWaves; i++ ) {
      s.append( " Wave " + i + ": " + frequency[i] + " Hz, " + amplitude[i] +
                " amplitude units, " + initialPhase[i] + " radians" + newline );
    }

    return new String(s);
  }
}
