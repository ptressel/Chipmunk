package ptressel.soundutils;

import java.util.concurrent.*;

/**
 * Fill buffers with samples from a Wave and queue them up.
 *
 * @see PlayWave for more info.
 */

public class WaveGen extends Thread {

  // Instance data

  /** Flag that tells us when to quit. */
  private boolean quitFlag = false;

  /** Sample provider */
  private final Wave wave;

  /** Transfer queue */
  private final LinkedBlockingQueue<byte[]> transferQueue;

  /** Buffer pool */
  private final LinkedBlockingQueue<byte[]> bufferPool;

  // Debugging info -- note we only need audio format info for debug
  // messages.

  /** Want debugging statistics printed? */
  private boolean debug;

  /** Want msgs every time the output thread might block? */
  private boolean debugBlocking;

  /** Print debug stats after this many loop passes */
  private int statLoopInterval;

  /** Forgetting rate for moving averages */
  private double statForgettingRate;

  /** Audio format: Bytes per sample */
  private int sampleBytes;

  /** Audio format: True if signed */
  private boolean signed;

  /** Audio format: True if big endian */
  private boolean bigEndian;

  /** Counter to tell when to print debugging statistics. */
  private int statCtr = 0;

  /** Moving average value of bytes received. */
  private double statAvgValue = 0;

  /** Moving average of transfer queue length. */
  private double statAvgQueueLen = 0;

  /** Moving average of buffer pool length. */
  private double statAvgPoolLen = 0;

  // Constructor

  /**
   * Make a WaveGen thread.
   *
   * @param wave the source for samples
   * @param transferQueue the queue into which full buffers will be placed
   * @param bufferPool a source of empty buffers
   */
  public WaveGen( Wave wave,
                  LinkedBlockingQueue<byte[]> transferQueue,
                  LinkedBlockingQueue<byte[]> bufferPool ) {

    this.wave = wave;
    this.transferQueue = transferQueue;
    this.bufferPool = bufferPool;

    debug = false;
    debugBlocking = false;
  }

  /**
   * Make a WaveGen thread with the given debug options, including audio
   * format info.
   *
   * @param debug want debug messages?
   * @param debugBlocking want grossly excessive messages around operations
   * that might block?
   * @param statLoopInterval print debug stats after this many loop passes
   * @param statForgettingRate fraction of old "average" to keep in moving
   * averages
   * @param sampleBytes audio format: bytes per sample
   * @param signed audio format: true if signed
   * @param bigEndian audio format: true if big endian
   *
   * @param wave the source for samples
   * @param transferQueue the queue into which full buffers will be placed
   * @param bufferPool a source of empty buffers
   */
  public WaveGen( Wave wave,
                  LinkedBlockingQueue<byte[]> transferQueue,
                  LinkedBlockingQueue<byte[]> bufferPool,
                  boolean debug, boolean debugBlocking, int statLoopInterval,
                  double statForgettingRate, int sampleBytes,
                  boolean signed, boolean bigEndian ) {

    // Standard setup.
    this( wave, transferQueue, bufferPool );

    // Now set the debug options.
    this.debug = debug;
    this.debugBlocking = debugBlocking;
    this.statLoopInterval = statLoopInterval;
    this.statForgettingRate = statForgettingRate;
    this.sampleBytes = sampleBytes;
    this.signed = signed;
    this.bigEndian = bigEndian;

    // Start the buffer pool stats off with the right value.
    if ( debug ) statAvgPoolLen = bufferPool.size();
  }

  // Thread methods

  /**
   * The run loop:  Get a buffer from the buffer pool, if there are any, else
   * wait.  Fill it from the Wave.  Queue it up on the transfer queue.  Repeat
   * until our quit flag goes high.
   */
  @Override
  public void run() {

    // The buffer we're moving from capture to transfer queue.
    byte[] buffer = null;

    runLoop:
    while(true) {
      // First check if we've been told to quit.
      if ( quitFlag ) break runLoop;

      // Get a buffer from the buffer pool, if possible, else wait.
      // Since we're not reading from a capture device that might drop
      // sound if we're not waiting on it, we can block on the buffer
      // pool.
      takeLoop:
      while(true) {
        try {
          if ( debugBlocking ) {
            System.out.println( "WaveGen: about to take from buffer pool." );
          }
          buffer = bufferPool.take();
          if ( debugBlocking ) {
            System.out.println( "WaveGen: take finished." );
          }
          break takeLoop;
        }
        catch( InterruptedException e ) {
          // Check for quit.
          if ( quitFlag ) break runLoop;
        }
      }

      // Ask our Wave to fill the buffer.
      try {
        wave.insertNextSamples( buffer );
      }
      catch( IllegalArgumentException e ) {
        // Bug in whatever set this up.
        System.out.println(
          "WaveGen: buffer size must be multiple of sample size" );
        break runLoop;
      }

      // Put it in the transfer queue
      putLoop:
      while(true) {
        try {
          if ( debugBlocking ) {
            System.out.println( "WaveGen: about to put to transfer queue." );
          }
          transferQueue.put( buffer );
          if ( debugBlocking ) {
            System.out.println( "WaveGen: put finished." );
          }
          break putLoop;
        }
        catch( InterruptedException e ) {
          // Check for quit.  
          if ( quitFlag ) break runLoop;
        }
      }

      // Debugging statistics:  Get average of samples in this buffer, and
      // update the moving average.  If we're at the stats printing interval,
      // do it and reset the counter.
      if ( debug ) {
        // Get average of this buffer.
        int avg = 0;
        for ( int i = 0; i < buffer.length; i += sampleBytes ) {
          avg += SoundUtils.extractSample( buffer, i, sampleBytes,
                   signed, bigEndian );
        }
        avg /= ( buffer.length / sampleBytes );
        // Update moving average of sound values.
        statAvgValue
          = statForgettingRate * statAvgValue +
            ( 1 - statForgettingRate ) * avg;
        // Update moving average of queue lengths.
        statAvgQueueLen
          = statForgettingRate * statAvgQueueLen +
            ( 1 - statForgettingRate ) * transferQueue.size();
        statAvgPoolLen
          = statForgettingRate * statAvgPoolLen +
            ( 1 - statForgettingRate ) * bufferPool.size();
        // Want to print them now?
        if ( statCtr >= statLoopInterval ) {
          // System.out.println does a good job of not dumping one thread's
          // line right in the middle of another.
          System.out.println( "WaveGen: # passes = " + statCtr );
          System.out.println( "WaveGen: current sound value = " + avg );
          System.out.println( "WaveGen: avg sound value = " + statAvgValue );
          System.out.println( "WaveGen: avg queue len = " + statAvgQueueLen);
          System.out.println( "WaveGen: avg pool len = " + statAvgPoolLen);
          statCtr = 0;
        } else {
          statCtr++;
        }
      }

      // Forget we had this buffer.
      buffer = null;

      // Give the poor output thread a break.
      yield();
    }

    // We've been told to quit -- fall out of the loop.
  }

  // Other methods

  /**
   * Turn on the quit flag.  This is the officially sanctioned method for
   * stopping this thread.  Do not call stop() -- you have been warned.
   */
  public void quit() {
    quitFlag = true;
  }
}
