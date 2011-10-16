package ptressel.soundutils;

/**
 * <p>Provide WAVE (PCM) data with a mixture of tones.  Data can be extracted
 * in a stream of single samples, or by filling a buffer.  Once initialized
 * with wave (specified by frequency, amplitude, initial phase) and audio
 * format information (sample rate, number of bytes per sample, signed or
 * not, endianness), the CompositeWave can be repeatedly called to get the
 * next sample or sequence of samples.
 *
 * <p>The user supplies arrays of the same length for each wave parameter,
 * where each wave's parameters are at the same index in each array.
 *
 * <p>CompositeWave currently only supports mono, signed, one or two bytes
 * per sample.
 *
 * <p>In order to get low latency, a Wave can be used in a separate thread,
 * with a pool of buffers to fill and place on a queue.  Then when a new
 * buffer of samples is needed, the consumer thread can get it from the
 * queue of pre-filled buffers.  After the consumer is done with it, the
 * buffer can be placed back in the empty pool.
 *
 * <p>Performance caution:
 * CompositeWave must, for each sample, compute and sum the unformatted
 * values of each pure tone, and then format the sum.  If this is too slow,
 * then a specialized version of CompositeWave can be written, that accepts
 * only a set of waves whose frequencies are integral submultiples of the
 * sample rate (so that there are an integral number of samples per cycle
 * for each wave), and whose periods are not mutually prime (and whose
 * least common multiple is not excessive).  Then it can produce one set of
 * formatted samples covering the entire least common multiple of all
 * periods (more precisely, of the number of samples over a period), and
 * repeatedly loop through that array, copying out samples as needed.
 * System.arraycopy can be used when multiple samples are requested.
 *
 * Implementation note:
 * Although it would be possible to have CompositeWave construct a number
 * of PureWaves and request and combine their samples, that would only be
 * useful for raw samples, as formatted samples can't be directly added,
 * and it would hurt performance.
 */
public class CompositeWave implements Wave {

  // Instance data

  // Audio format information -- we don't use a Java Sound AudioFormat
  // because there's no actual connection to Java Sound, or even, really
  // to sound...  But we do use the same types for this data as the
  // parameters in the AudioFormat constructors.

  /** Sample rate (note this is a float in an actual AudioFormat) */
  private double sampleRate;

  /** Number of bytes per sample */
  private final int sampleBytes = 2;

  // A composite wave is specified as a collection of pure tones.  Each
  // pure tone has its own set of values equivalent to those in PureTone.
  // One tone's values are at the same array index throughout the following
  // arrays.

  // Wave information supplied by user

  /** Frequency */
  private double[] frequency;

  /** Initial phase in radians */
  private double[] initialPhase;

  /** Amplitude */
  private int[] amplitude;

  // Internal data

  /** Number of pure tone waves in the composite */
  int nWaves;

  /** Phases in each wave for next sample to return */
  private double[] nextSamplePhase;

  /** Phase increment per sample */
  private double[] phaseIncrement;

  // Constructor

  // TODO Perhaps default phase such that initial amplitude is zero,
  // not maximum.  An abrupt change in amplitude causes a click in
  // the audio.
  
  /**
   * <p>Set up a generator for a composite wave.
   *
   * <p>Each component of the composite wave is a pure tone, specified by
   * its frequency, amplitude, and phase at any point in time.
   *
   * <p>Amplitude is the magnitude of the wave's peaks -- currently only
   * signed PCM is supported, so if the supplied amplitudes sum to A, the
   * sample values will range from -A to +A, so long as that is within
   * the range possible for the number of bytes per sample, otherwise,
   * samples that are out of range will be clipped.
   *
   * To avoid clipping, the sum of the specified amplitudes should be within
   * Byte.MIN_VALUE to Byte.MAX_VALUE if there is one byte per sample, or
   * for two bytes, it should be within Short.MIN_VALUE to Short.MAX_VALUE.
   *
   * <p>Phase is measured in radians starting from 0 at maximum amplitude
   * (like a cosine).  The available initial phases are those at sample times
   * k/s where s is the sample rate and k is an integer.  That is, the phases
   * of the samples that can be produced are 2 pi k f / s where f is the
   * frequency of the wave.  The supplied initial phase will be rounded to
   * the sample time with the nearest phase.
   *
   * <p>The initial phase is likely irrelevant unless doing phase shift keying.
   * That argument may be left null to indicate initial phase zero for all.
   *
   * <p>The very first sample returned will have the specified initial phases
   * for each tone.  Subsequent samples will continue the cycle until the
   * phase is reset.  Note there is no reason to reset the phase when doing
   * frequency or amplitude modulation, or simply denoting symbols by
   * frequencies.
   *
   * <p>Only single-channel, signed PCM is currently supported, and only
   * one or two byte sample sizes.  If other values are specified, the
   * constructor will throw an IllegalArgumentException.
   *
   * @param frequency frequencies for each wave
   * @param initialPhase for each wave, the phase in radians of the very
   * first sample that will be returned
   * @param amplitude the maximum amplitude for each wave
   *
   * @param sampleRate the sample rate
   *
   * @throws IllegalArgumentException if an unsupported audio format is
   * specified, or if the numbers of frequencies, amplitudes, and initial
   * phases (if supplied) are not the same
   */
  public CompositeWave( double[] frequency, int[] amplitude,
           double[] initialPhase, double sampleRate )
      throws IllegalArgumentException {

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
    nextSamplePhase = this.initialPhase.clone();

    // We do use the audio format values (or at least, the ones we currently
    // support) past construction.
    this.sampleRate = sampleRate;

    // Make a place for the phase increments per sample.
    phaseIncrement = new double[nWaves];

    // Compute the phase change per sample for each wave.  (See long
    // explanation in PureWave.)
    for ( int i = 0; i < nWaves; i++ ) {
      phaseIncrement[i] = 2 * Math.PI * frequency[i] / sampleRate;
    }
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
    // Add up the values of all the waves at the current sample.
    double value = 0;
    for ( int w = 0; w < nWaves; w++ ) {
      // Compute the current amplitude of this wave.
      value += Math.cos( nextSamplePhase[w] ) * amplitude[w];
      // Update the phase.  Wrap at 2pi.
      nextSamplePhase[w]
        = ( nextSamplePhase[w] + phaseIncrement[w] ) % (2*Math.PI);
    }
    return value;
  }

  /**
   * <p>Get the next (single) formatted sample and store it at the given offset
   * in the supplied buffer.
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

    // Get the next unformatted sample and advance the sample counters.
    double value = getNextRawSample();
    // Clip, format and store the composite value.
    SoundUtils.insertSample( value, buffer, byteOffset, sampleBytes );
  }

  /**
   * <p>Get the next specified number of unformatted samples from the wave
   * and store them in the supplied double array at the given offset.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * @param buffer to put the samples in
   * @param offset starting index for the samples
   * @param nSamples number of samples to provide
   *
   * @throws ArrayIndexOutOfBoundsException if the samples would extend past
   * the end of the buffer
   */
  public void getNextRawSamples( double[] buffer, int offset, int nSamples )
      throws ArrayIndexOutOfBoundsException, IllegalArgumentException {

    if ( ( buffer.length - offset ) < nSamples )
      throw new IndexOutOfBoundsException(
        "Requested number of samples will not fit." );

    for ( int i = offset; i < offset + nSamples; i++ ) {
      buffer[i] = getNextRawSample();
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
   * @param nSamples number of samples to insert
   *
   * @throws ArrayIndexOutOfBoundsException if the samples would extend past
   * the end of the buffer
   * @throws IllegalArgumentException if the audio format is not supported
   */
  public void insertNextSamples( byte[] buffer, int byteOffset, int nSamples )
      throws ArrayIndexOutOfBoundsException, IllegalArgumentException {

    for ( int i = byteOffset; i < byteOffset + sampleBytes * nSamples;
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
   * <p>Fill the given buffer with unformatted samples from the wave.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * @param buffer the buffer to fill
   */
  public void getNextRawSamples( double[] buffer ) {
    // This won't throw exceptions due to length, so don't make the caller
    // deal with that.
    try {
      getNextRawSamples( buffer, 0, buffer.length );
    }
    catch ( IndexOutOfBoundsException e ) {}
  }

  /**
   * <p>Fill the given buffer with samples from the wave.
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
  public void insertNextSamples( byte[] buffer ) throws IllegalArgumentException {

    // This won't throw exceptions due to length, so don't make the caller
    // deal with that.
    try {
      insertNextSamples( buffer, 0, buffer.length / sampleBytes );
    }
    catch ( ArrayIndexOutOfBoundsException e ) {}
  }

  /**
   * Set phases of each wave for next sample to the specified values, in
   * radians.  This is not useful unless doing phase shift keying.
   *
   * @param nextPhase phases in radians for each wave in next sample to return
   */
  public void setPhase( double[] nextPhase ) {
    // Don't bother constraining to 0..2pi.
    for ( int i = 0; i < nWaves; i++ )
      nextSamplePhase[i] = nextPhase[i];
  }

  /**
   * Reset phase for next sample to originally-specified initial phase(s).
   */
  public void resetPhase() {
    // It's surely a tiny array so don't use System.arraycopy.
    for ( int i = 0; i < nWaves; i++ )
      nextSamplePhase[i] = initialPhase[i];
  }

  /** Test */
  public static void main( String[] args ) {

    // Choose params so there are only a few samples through the wave.

    // First try something simple -- two waves an octave apart.
    double freq2Tones[] = { 400, 800 };
    int max = Byte.MAX_VALUE / 2;
    int[] amp2Tones = { max, max };
    int oneByte = 1;
    boolean signed = true;
    boolean bigEndian = true;
    Wave wave = new CompositeWave( freq2Tones, amp2Tones, null, 8*800,
                                   oneByte, signed, bigEndian );
    System.out.println( "1 byte, signed, big endian" );
    System.out.println( "Freq: " + freq2Tones[0] + ", " + freq2Tones[1]
      + ", Amp: " + amp2Tones[0] + ", " + amp2Tones[1] );
    System.out.println();

    // Get enough samples to see one cycle of the 400Hz tone, with 8 sample
    // over a cycle for the 800Hz tone, i.e. 16 samples.  No, get 3x that.
    byte[] buffer = new byte[3*16];

    // Request some samples.
    wave.insertNextSamples( buffer );

    // Show what's in there.  Since we asked for 1 byte samples, we don't
    // need any byte swapping or repacking.
    System.out.println( "1-byte samples, two tones an octave apart:" );
    for ( int i = 0; i < buffer.length; i++ ) {
      System.out.println( i + ": " + buffer[i] );
    }
    System.out.println();

    // Same but little endian.
    boolean notBigEndian = false;
    wave = new CompositeWave( freq2Tones, amp2Tones, null, 8*800,
                              oneByte, signed, notBigEndian );
    System.out.println( "1 byte, signed, little endian" );
    System.out.println( "Freq: " + freq2Tones[0] + ", " + freq2Tones[1]
      + ", Amp: " + amp2Tones[0] + ", " + amp2Tones[1] );
    System.out.println();

    // Request some samples.
    wave.insertNextSamples( buffer );

    // Show what's in there.  Since we asked for 1 byte samples, we don't
    // need any byte swapping or repacking.
    System.out.println( "1-byte samples, two tones an octave apart:" );
    for ( int i = 0; i < buffer.length; i++ ) {
      System.out.println( i + ": " + buffer[i] );
    }
    System.out.println();

    // This time, non-zero phase.
    double[] phase2 = { 1, 3 };
    wave = new CompositeWave( freq2Tones, amp2Tones, phase2, 8*800,
                              oneByte, signed, bigEndian );
    System.out.println( "1 byte, signed, big endian, non-zero phase" );
    System.out.println( "Freq: " + freq2Tones[0] + ", " + freq2Tones[1]
      + ", Amp: " + amp2Tones[0] + ", " + amp2Tones[1] );
    System.out.println();

    // Request some samples.
    wave.insertNextSamples( buffer );

    // Show what's in there.  Since we asked for 1 byte samples, we don't
    // need any byte swapping or repacking.
    System.out.println( "1-byte samples, two tones an octave apart:" );
    for ( int i = 0; i < buffer.length; i++ ) {
      System.out.println( i + ": " + buffer[i] );
    }
    System.out.println();

    // This time, one zero phase, one non-zero phase.
    double[] phaseZnZ = { 0, 1 };
    wave = new CompositeWave( freq2Tones, amp2Tones, phaseZnZ, 8*800,
                              oneByte, signed, bigEndian );
    System.out.println( "1 byte, signed, big endian, one zero phase" );
    System.out.println( "Freq: " + freq2Tones[0] + ", " + freq2Tones[1]
      + ", Amp: " + amp2Tones[0] + ", " + amp2Tones[1] );
    System.out.println();

    // Request some samples.
    wave.insertNextSamples( buffer );

    // Show what's in there.  Since we asked for 1 byte samples, we don't
    // need any byte swapping or repacking.
    System.out.println( "1-byte samples, two tones an octave apart:" );
    for ( int i = 0; i < buffer.length; i++ ) {
      System.out.println( i + ": " + buffer[i] );
    }
    System.out.println();

    double[] freq3Tones = { 400, 600, 800 };
    // Stay just a titch down from the max value for 3 equal amplitude waves.
    int max1Byte = Byte.MAX_VALUE / 3 - 3;
    int max2Byte = Short.MAX_VALUE / 3 - 3;
    int[] amp1Byte = { max1Byte, max1Byte, max1Byte };
    int[] amp2Byte = { max2Byte, max2Byte, max2Byte };

    // Need a sampling rate that will give us an integral number of samples
    // over the period of each wave, so we can see them come back into phase.
    // If f is the freq of a wave, period T = 1/f, and fs the sampling freq,
    // period Ts = 1/fs, then there are T/Ts = fs/f samples per period.  So
    // we need an fs into which all 3 f's divide evenly, e.g. the least common
    // multiple.  Here, that's 3*800 = 2400.  But that only gives us 3 samples
    // over the 800Hz signal -- want an even # so we see both peaks, so get
    // 6*800.  The number of samples over the 800Hz wave's period is 6, that
    // for 600Hz is 8, for 400Hz, it's 12.  These will come back into phase
    // after the LCM of the periods, i.e. 24 samples.  Get several times that,
    // in case something isn't matching up at the end.
    //
    // This case was verified as follows using Matlab:
    // n=[1:24];             % Get 24 sample points
    // s800 = (n-1)*2*pi/6;  % 6 samples per cycle
    // s600 = (n-1)*2*pi/8;  % 8 per cycle
    // s400 = (n-1)*2*pi/12; % 12 per cycle
    // v800 = 39*cos(s800);  % values for the s800 samples
    // v600 = 39*cos(s600);
    // v400 = 39*cos(s400);
    // v = v800 + v600 + v400;

    int samp3 = 6*800;
    wave = new CompositeWave( freq3Tones, amp1Byte, null, samp3,
                              oneByte, signed, bigEndian );
    System.out.println( "1 byte, signed, big endian" );
    System.out.println( "Freq: " + freq3Tones[0] + ", " + freq3Tones[1] + ", "
      + freq3Tones[2] + ", Amp: " + amp1Byte[0] + ", " + amp1Byte[1] + ", "
      + amp1Byte[2] + ", Sampling: " + samp3 );
    System.out.println();

    buffer = new byte[3*24];

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
    int twoBytes = 2;
    wave = new CompositeWave( freq3Tones, amp2Byte, null, samp3,
                              twoBytes, signed, bigEndian);
    System.out.println( "2 bytes, signed, big endian" );
    System.out.println( "Freq: " + freq3Tones[0] + ", " + freq3Tones[1] + ", "
      + freq3Tones[2] + ", Amp: " + amp2Byte[0] + ", " + amp2Byte[1] + ", "
      + amp2Byte[2] + ", Sampling: " + samp3 );
    System.out.println();
    // Need buffer that's twice as large.
    buffer = new byte[2*3*24];
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

    // Try out the case we're using for testing FFT.
    System.out.println();
    System.out.println( "Four pure tones" );

    double[] freq = { 400, 420, 1000, 7999 };  // Note near Nyquist freq
    int[] amp = { 1, 3, 5, 2 };
    double[] phase = { 0, 1, 2, 3 };  // say anything
    int sampleRate = 16000;

    System.out.println( "2 bytes, signed, little endian, some non-zero phases");
    System.out.println( "Freq: " + freq[0] + ", " + freq[1] + ", "
      + freq[2] + ", " + freq[3] + ", Amp: " + amp[0] + ", " + amp[1] + ", "
      + amp[2] + ", " + amp[3] + ", Sampling: " + sampleRate );

    wave = new CompositeWave( freq, amp, phase, sampleRate, twoBytes,
                              signed, notBigEndian );

    // Get a power of 2 worth of samples.
    double[] samples = new double[2048];
    wave.getNextRawSamples( samples );
    // Make sure we got something...
    System.out.println( "A few samples:" );
    for ( int i = 0; i < 20; i++ ) {
      System.out.println( i + ": " + samples[i] );
    }
    System.out.println();

    // Same, but with no phase specified.
    System.out.println( "2 bytes, signed, little endian, no phases" );
    System.out.println( "Freq: " + freq[0] + ", " + freq[1] + ", "
      + freq[2] + ", " + freq[3] + ", Amp: " + amp[0] + ", " + amp[1] + ", "
      + amp[2] + ", " + amp[3] + ", Sampling: " + sampleRate );

    wave = new CompositeWave( freq, amp, null, sampleRate, twoBytes,
                              signed, notBigEndian );

    // Get a power of 2 worth of samples.
    samples = new double[2048];
    wave.getNextRawSamples( samples );
    // Make sure we got something...
    System.out.println( "A few samples:" );
    for ( int i = 0; i < 20; i++ ) {
      System.out.println( i + ": " + samples[i] );
    }
    System.out.println();

  }
}
