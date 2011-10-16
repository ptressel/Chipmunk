package ptressel.soundutils;

/**
 * <p>Provide WAVE (PCM) data with a single tone.  Data can be extracted in
 * a stream of single samples, or by filling a buffer.  Once initialized
 * with wave (specified by frequency, amplitude, initial phase) and audio
 * format information (sample rate, number of bytes per sample, signed or
 * not, endianness), the PureWave can be repeatedly called to get the next
 * sample or sequence of samples.
 *
 * <p>PureWave currently only supports mono, signed, one or two bytes
 * per sample.
 *
 * <p>In order to get low latency, a Wave can be used in a separate thread,
 * with a pool of buffers to fill and place on a queue.  Then when a new
 * buffer of samples is needed, the consumer thread can get it from the
 * queue of pre-filled buffers.  After the consumer is done with it, the
 * buffer can be placed back in the empty pool.
 *
 * <p>Performance caution:
 * PureWave must, for each sample, compute the unformatted value of the
 * current amplitude, then format it.  This involves calling trig functions.
 * If this is too slow, then a specialized version of PureWave can be written,
 * that accepts only frequencies that are integral submultiples of the 
 * sample rate (so that there are an integral number of samples per cycle)
 * Then it can produce one set of formatted samples covering the entire cycle,
 * and repeatedly loop through that array, copying out samples as needed.
 * System.arraycopy can be used when multiple samples are requested.
 */
public class PureWave implements Wave {

  // Instance data

  // Audio format information -- we don't use a Java Sound AudioFormat
  // because there's no actual connection to Java Sound, or even, really
  // to sound...  But we do use the same types for this data as the
  // parameters in the AudioFormat constructors.

  /** Sample rate (note this is a float in an actual AudioFormat) */
  private final double sampleRate;

  /** Number of bytes per sample */
  private final int sampleBytes;

  /** Is sample value signed or unsigned? */
  private final boolean signed;

  /** Is sample value big Endian or small? */
  private final boolean bigEndian;

  // Wave information supplied by user

  /** Frequency */
  private final double frequency;

  /** Initial phase in radians */
  private final double initialPhase;

  /** Amplitude */
  private final int amplitude;

  // Internal data

  /** Phase of next sample */
  private double nextSamplePhase = 0;

  /** Phase increment per sample */
  private double phaseIncrement;

  // Constructor

  /**
   * <p>Set up a generator for a single-frequency wave.
   *
   * <p>Waves are specified here by their frequency, amplitude, and phase.  
   *
   * <p>Amplitude is the magnitude of the wave's peaks -- currently only
   * signed PCM is supported, so if the supplied amplitude is A, the sample
   * values will range from -A to +A.  The scale of the amplitude should
   * be appropriate to the number of bytes per sample.  That is, for one
   * byte, the amplitude should be within Byte.MIN_VALUE to Byte.MAX_VALUE.
   * For two bytes, it should be within Short.MIN_VALUE to Short.MAX_VALUE.
   * Note if this wave is to be combined with another, the sum of their
   * amplitudes should fit in the appropriate range, else the combined
   * samples will be clipped.
   *
   * <p>Phase is measured in radians starting from 0 at maximum amplitude
   * (like a cosine).  The initial phase is likely irrelevant if not doing
   * phase shift keying.
   *
   * <p>The very first sample returned will have the initial phase.
   * Subsequent samples will continue the cycle until the phase is set
   * or reset to its initial value.  Note there is no reason to reset the
   * phase when doing frequency or amplitude modulation, or simply denoting
   * symbols by frequencies.  In fact, doing so can introduce unwanted
   * high frequency components due to any abrupt amplitude change.
   *
   * <p>Only single-channel, signed PCM is currently supported, and only
   * one or two byte sample sizes.  If other values are specified, the
   * constructor will throw an exception.
   *
   * @param frequency the frequency
   * @param initialPhase the phase in radians of the very first sample that
   * will be returned
   * @param amplitude the maximum amplitude
   *
   * @param sampleRate the sample rate
   * @param sampleBytes the number of bytes per sample (currently must be 1
   * or 2)
   * @param signed true if samples are signed (currently only signed is
   * supported)
   * @param bigEndian true if samples should be provided in big endian order
   * (relevant only if sampleBytes is > 1)
   *
   * @throws IllegalArgumentException if an unsupported value of the sample
   * bytes or signed is given, or if the amplitude is out of range for the
   * number of sample bytes
   */
  public PureWave( double frequency, int amplitude, double initialPhase, 
                   double sampleRate, int sampleBytes, boolean signed,
                   boolean bigEndian ) throws IllegalArgumentException {

    // Check for supported options and amplitude range.
    if ( !signed )
      throw new IllegalArgumentException(
        "Only signed values are currently supported." );
    if ( sampleBytes == 1 ) {
      if ( amplitude < Byte.MIN_VALUE || amplitude > Byte.MAX_VALUE )
        throw new IllegalArgumentException(
          "Amplitude out of range for a one byte sample size." );
    } else if ( sampleBytes == 2 ) {
      if ( amplitude < Short.MIN_VALUE || amplitude > Short.MAX_VALUE )
        throw new IllegalArgumentException(
          "Amplitude out of range for a two byte sample size." );
    } else {
      throw new IllegalArgumentException(
        "Only sample sizes of 1 or 2 bytes are currently supported." );
    }

    // Store the specified wave parameters.  (Frequency is not used past
    // construction.
    this.frequency = frequency;
    this.amplitude = amplitude;
    this.initialPhase = initialPhase;

    // We do use the audio format values (or at least, the ones we currently
    // support) past construction.
    this.sampleRate = sampleRate;
    this.sampleBytes = sampleBytes;
    this.signed = signed;
    this.bigEndian = bigEndian;

    // What is the phase change per sample?  The period of one cycle of
    // frequency f is 1/f.  For sampling frequency s, there is one sample
    // every 1/s, so the fraction of a cycle by which the wave advances on
    // each sample is (1/s) / (1/f) = f/s.  One cycle is 2pi radians, so the
    // number of radians by which the wave advances per sample is 2pi f/s.
    // Note that f/s is under no obligation to come out an exact value -- it
    // may have roundoff error.  If so, then the roundoff error will
    // accumulate gradually.  An alternative would be to count samples and
    // (re)compute the phase on each sample k as 2pi k f/s, leaving the
    // division by s til last.  Unfortunately, that would eventually overflow.
    // For now, tolerate the roundoff.
    phaseIncrement = 2 * Math.PI * frequency / sampleRate;

    // Note also that we are going to have to call trig functions on each
    // sample.  *If* we knew that the sampleRate was a multiple of the
    // frequency, we could precompute full cycles of samples, and just
    // cycle through that array.  If the user wanted a buffer full of
    // samples, we could copy the cycle array into the user's buffer
    // repeatedly.  This is how the original version was implemented,
    // until it became clear that one could not require users to choose
    // only frequencies that were a submultiple of the sample rate,
    // especially not for musical scales...
  }

  // Public instance methods

  /**
   * Get the next single, unformatted sample and return it as a double.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * @return the next unformatted sample
   */
  public double getNextRawSample() {
    // Compute the current amplitude of the wave.
    double value = Math.cos( nextSamplePhase ) * amplitude;
    // Update the phase.  Wrap at 2pi.
    nextSamplePhase = ( nextSamplePhase + phaseIncrement ) % (2*Math.PI);
    return value;
  }

  /**
   * <p>Get the next (single) formatted sample and store it at the specified
   * location in the supplied buffer.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * <p>If there aren't enough bytes in the caller's buffer to accomodate the
   * formatted sample, then insertNextSample will throw
   * ArrayIndexOutOfBoundsException.
   * The wave's phase will advance to the next sample regardless.
   *
   * @param buffer to put the sample in
   * @param byteOffset starting index for the sample
   *
   * @throws ArrayIndexOutOfBoundsException if the sample extends past the
   * end of the buffer
   * @throws IllegalArgumentException if the audio format is not supported
   */
  public void insertNextSample( byte[] buffer, int byteOffset )
      throws ArrayIndexOutOfBoundsException, IllegalArgumentException {

    // Get an unformatted sample, update phase.
    double value = getNextRawSample();
    // Format it and put it in the user's buffer.
    SoundUtils.insertSample( value, buffer, byteOffset, sampleBytes, signed,
                             bigEndian );
  }

  /**
   * <p>Get the next specified number of unformatted samples from the wave
   * and store them in the supplied double array at the given offset.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * @param rawBuffer storage for the unformatted samples
   * @param offset starting index (in bytes) for new samples
   * @param nSamples number of samples to provide
   *
   * @throws ArrayIndexOutOfBoundsException if the samples would extend past
   * the end of the buffer
   */
  public void getNextRawSamples( double[] rawBuffer, int offset, int nSamples )
      throws ArrayIndexOutOfBoundsException {

    if ( ( rawBuffer.length - offset ) < nSamples )
      throw new IndexOutOfBoundsException(
        "Requested number of samples will not fit." );

    for ( int i = offset; i < offset + nSamples; i++ ) {
      rawBuffer[i] = getNextRawSample();
    }
  }

  /**
   * <p>Get the next specified number of samples (note: samples, not bytes)
   * from the wave and store them in the supplied byte array at the given
   * offset.
   *
   * <p>The starting offset is in bytes, and there is no requirement that the
   * offset be aligned on a sampleBytes boundary, nor that the buffer length
   * be an integral number of samples, as the caller may have other information
   * stored in the buffer (e.g. a header for the audio format).  For the same
   * reason, this does not wrap if the number of samples will not fit in the
   * tail of the buffer.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * <p>If there is not enough room in the user's buffer for the requested
   * number of samples, then as many samples as will fit will be written in,
   * and insertNextSamples will throw ArrayIndexOutOfBoundsException.
   *
   * @param buffer to put the samples in
   * @param byteOffset starting index for the samples
   * @param nSamples number of samples to provide
   *
   * @throws ArrayIndexOutOfBoundsException if the samples would extend past
   * the end of the buffer
   * @throws IllegalArgumentException if the audio format is not supported
   */
  public void insertNextSamples( byte[] buffer, int byteOffset, int nSamples )
      throws ArrayIndexOutOfBoundsException, IllegalArgumentException {

    for ( int i = byteOffset; i < byteOffset + nSamples*sampleBytes;
          i += sampleBytes ) {
      // Get a sample, format it, put it in the buffer, increment the phase.
      insertNextSample( buffer, i );
    }
  }

  /**
   * <p>Get the next specified number of bytes worth of samples, format the
   * samples, and store them in the supplied byte array at the given offset.
   * If the number of bytes is not a multiple of the sample size, the number
   * of bytes provided will be rounded down to the nearest whole sample.
   * (This is for the benefit of callers that don't know about the audio format,
   * mainly MockTargetDataLine.)
   *
   * <p>The starting offset is in bytes, and there is no requirement that the
   * offset be aligned on a sampleBytes boundary, nor that the buffer length
   * be an integral number of samples, as the caller may have other information
   * stored in the buffer (e.g. a header for the audio format).  For the same
   * reason, this does not wrap if the number of samples will not fit in the
   * tail of the buffer.  We do, however, require that the number of bytes
   * requested be a multiple of the sample size.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * @param buffer storage for the samples
   * @param byteOffset starting index in buffer at which to store the
   * samples
   * @param nBytes number of bytes to insert
   */
  public void insertNextBytes( byte[] buffer, int byteOffset, int nBytes ) {
    if ( ( nBytes % sampleBytes ) != 0 ) {
      throw new IllegalArgumentException(
              "Requested number of bytes is not a multiple of frame size." );
    }
    insertNextSamples( buffer, byteOffset, nBytes / sampleBytes );
  }

  /**
   * <p>Fill the supplied buffer with unformatted samples from the wave.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * @param rawBuffer storage for samples
   */
  public void getNextRawSamples( double[] rawBuffer ) {
    // This won't throw exceptions due to length, so don't make the caller
    // deal with that.
    try {
      getNextRawSamples( rawBuffer, 0, rawBuffer.length );
    }
    catch ( ArrayIndexOutOfBoundsException e ) {}
  }

  /**
   * <p>Fill the supplied buffer with samples from the wave.
   *
   * <p>This version of insertNextSamples is not intended for buffers that include
   * non-sample data such as headers.  Instead, use
   * {@link insertNextSamples(byte[], int, int)}.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * <p>If the buffer is not an integral multiple of the sample size, this
   * will fill in as many samples as it can, leaving the remaining bytes
   * unmodified.
   *
   * @param buffer the buffer to fill
   *
   * @throws IllegalArgumentException if the audio format is not supported
   */
  public void insertNextSamples( byte[] buffer )
      throws IllegalArgumentException {

    // This won't throw exceptions due to length, so don't make the caller
    // deal with that.
    try {
      insertNextSamples( buffer, 0, buffer.length / sampleBytes );
    }
    catch ( ArrayIndexOutOfBoundsException e ) {}
  }

  /**
   * Set phase for next sample to the specified value, in radians.  This is
   * not useful unless doing phase shift keying.
   *
   * @param nextPhase phase in radians of next sample to return
   */
  public void setPhase( double nextPhase ) {
    // Don't bother constraining to 0..2pi.
    nextSamplePhase = nextPhase;
  }

  /**
   * Reset phase for next sample to originally-specified initial phase(s).
   */
  public void resetPhase() {
    nextSamplePhase = initialPhase;
  }

  /** Test */
  public static void main( String[] args ) {

    // Choose params so there are only a few samples through the wave.
    // Copy out some samples and print them.
    Wave wave = new PureWave( 400, Byte.MAX_VALUE/2, 0, 8*400, 1, true, true );

    // Expect 8 samples per cycle, so get 2.5 cycles.
    byte[] buffer = new byte[20];

    // Request some samples.
    wave.insertNextSamples( buffer );

    // Show what's in there.  Since we asked for 1 byte samples, we don't
    // need any byte swapping or repacking.
    System.out.println( "1-byte samples:" );
    for ( int i = 0; i < buffer.length; i++ ) {
      System.out.println( i + ": " + buffer[i] );
    }
    System.out.println();

    // Repeat for two-byte samples.
    wave = new PureWave( 400, Short.MAX_VALUE/2, 0, 8*400, 2, true, true );
    // Need buffer that's twice as large.
    buffer = new byte[40];
    wave.insertNextSamples( buffer );
    System.out.println( "2-byte samples:" );
    for ( int i = 0; i < buffer.length; i += 2 ) {
      int v = SoundUtils.extractSample( buffer, i, 2, true, true );
      String s0 = SoundUtils.toHexString( buffer[i] );
      String s1 = SoundUtils.toHexString( buffer[i+1] );
      System.out.println( i + ": " + v +
        ", hi: " + buffer[i] + " (hex " + s0 + "), " + 
        ", lo: " + buffer[i+1] + " (hex " + s1 + ")" );
    }
    System.out.println();

    // Get some more samples one at a time.  These should follow on from the
    // 2-byte samples above.
    System.out.println( "More 2-byte samples, requested one at a time:" );
    // Although we don't have to, place the samples in successive locations
    // in the buffer, to be sure the offset is working.
    for ( int i = 0; i < buffer.length; i += 2 ) {
      wave.insertNextSample( buffer, i );
      int v = SoundUtils.extractSample( buffer, i, 2, true, true );
      String s0 = SoundUtils.toHexString( buffer[i] );
      String s1 = SoundUtils.toHexString( buffer[i+1] );
      System.out.println( i + ": " + v +
        ", hi: " + buffer[i] + " (hex " + s0 + "), " + 
        ", lo: " + buffer[i+1] + " (hex " + s1 + ")" );
    }
    System.out.println();

    // Get a buffer full of raw samples.
    System.out.println( "Still more -- a buffer full of raw samples:" );
    double[] dbuffer = new double[20];
    wave.getNextRawSamples( dbuffer );
    for ( int i = 0; i < dbuffer.length; i++ ) {
      System.out.println( i + ": " + dbuffer[i] );
    }
    System.out.println();

    // Get some raw samples, one at a time.
    System.out.println( "Yet more samples, raw, requested one at a time:" );
    double value = 0;
    for ( int i = 0; i < 20; i++ ) {
      value = wave.getNextRawSample();
      System.out.println( i + ": " + value );
    }
  }
}
