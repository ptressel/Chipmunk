package ptressel.soundutils;

/**
 * <p>Provide WAVE (PCM) data with a specified mix of tones.  Data can be
 * extracted in a stream of single samples, or by filling a buffer.  Once
 * initialized with a wave or set of waves to combine per channel (where each
 * wave is specified by frequency, amplitude, initial phase) and audio format
 * information (sample rate, number of bytes per sample, signed or not,
 * endianness), the Wave can be repeatedly called to get the next unformatted
 * sample or sequence of samples, or to format the next sample or samples and
 * store them in a buffer.
 * 
 * <p>The get* methods produce unformatted amplitudes, and are intended for use
 * when the caller merely wants raw wave amplitudes.  The insert* methods
 * produce formatted samples, and are intended for filling a buffer for output.
 * Calls to any get or insert method will advance the sample position (phase)
 * of the wave, so that successive calls produce contiguous samples.  It is
 * unlikely that the caller will want to mix calls to the get (unformatted)
 * and insert (formatted) methods.
 *
 * <p>Classes implementing this interface may provide more methods that are
 * specific to the type of wave they provide.  In particular, they may
 * provide methods to set the phase to a specified value for each single
 * frequency component of the waveform.
 *
 * <p>In order to get low latency, a Wave can be used in a separate thread,
 * with a pool of buffers to fill and place on a queue.  Then when a new
 * buffer of samples is needed, the consumer thread can get it from the
 * queue of pre-filled buffers.  After the consumer is done with it, the
 * buffer can be placed back in a pool of empty buffers.
 */
public interface Wave {

  /**
   * Get the next single, unformatted sample and return it as a double.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * @return the next unformatted sample
   */
  public double getNextRawSample();

  /**
   * <p>Get the next (single) sample and store it at the given offset
   * in the supplied buffer.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * @param buffer to put the sample in
   * @param byteOffset starting index for the sample
   *
   * @throws IndexOutOfBoundsException if the sample extends past the
   * end of the buffer
   */
  public void insertNextSample( byte[] buffer, int byteOffset );

  /**
   * <p>Get the next specified number of unformatted samples from the wave
   * and store them in the supplied double array at the given offset.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * @param rawBuffer storage for the unformatted samples
   * @param offset starting index for the samples
   * @param nSamples number of samples to provide
   *
   * @throws IndexOutOfBoundsException if the samples would extend past the
   * end of the buffer
   */
  public void getNextRawSamples( double[] rawBuffer, int offset, int nSamples )
      throws IndexOutOfBoundsException;

  /**
   * <p>Fill the supplied buffer with unformatted samples from the wave.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * @param rawBuffer buffer to fill
   */
  public void getNextRawSamples( double[] rawBuffer );

  /**
   * <p>Get the next specified number of samples (note: samples, not bytes)
   * from the wave, format them, and store them in the supplied byte array at
   * the given offset.
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
   * @param buffer storage for the samples
   * @param byteOffset starting index in buffer at which to store the
   * samples
   * @param nSamples number of samples to insert
   *
   * @throws IndexOutOfBoundsException if the samples would extend past the
   * end of the buffer
   */
  public void insertNextSamples( byte[] buffer, int byteOffset, int nSamples )
      throws IndexOutOfBoundsException;

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
   * tail of the buffer.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   *
   * @param buffer storage for the samples
   * @param byteOffset starting index in buffer at which to store the
   * samples
   * @param nBytes number of bytes to insert
   */
  public void insertNextBytes( byte[] buffer, int byteOffset, int nBytes )
      throws IndexOutOfBoundsException;

  /**
   * <p>Fill the supplied buffer with formatted samples from the wave.
   * 
   * <p>This version of insertNextSamples is not intended for buffers that
   * include non-sample data such as headers, as it will insert samples
   * starting at index 0.  Instead, use
   * {@link insertNextSamples(byte[], int, int)}.
   *
   * <p>If the buffer is not an integral multiple of the sample size, this
   * will fill in as many samples as it can, leaving the remaining bytes
   * unmodified.
   *
   * <p>The next call to any of the get or insert methods will return samples
   * following the samples returned by this call.
   * 
   * @param buffer storage for samples
   * 
   * @throws IllegalArgumentException if the buffer size is not an integral
   * multiple of the audio format's sample size
   */
  public void insertNextSamples( byte[] buffer )
      throws IllegalArgumentException;

  /**
   * Reset phase for next sample to originally-specified initial phase(s).
   */
  public void resetPhase();
}
