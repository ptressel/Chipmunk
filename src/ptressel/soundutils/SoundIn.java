package ptressel.soundutils;

import android.media.AudioFormat;
import android.media.AudioRecord;
import android.media.MediaRecorder.AudioSource;
import android.util.Log;

import java.util.concurrent.*;

/**
 * Copy sound from the capture line and queue it up.
 * 
 * @see SoundCopy for more info.
 */

public class SoundIn extends Thread {
	private static final String TAG = "Chipmunk";

	// The current system's newline string ("\n" for Unix, or "\r\n" for
	// Windows).
	private final static String NEWLINE = System.getProperty("line.separator");

	// Instance data

	/** Flag that tells us when to quit. */
	private boolean quitFlag = false;

	/** Capture line */
	// private final TargetDataLine capture;

	/** Transfer queue */
	private final LinkedBlockingQueue<short[]> transferQueue;

	/** Buffer pool */
	private final LinkedBlockingQueue<short[]> bufferPool;

	/** Buffer size */
	private final int bufferSize;

	private AudioRecord recorder;

	// Constructor

	/**
	 * Make a SoundIn thread.
	 * 
	 * @param capture
	 *            the input line from which to acquire sound
	 * @param transferQueue
	 *            the queue into which full buffers will be placed
	 * @param bufferPool
	 *            a source of empty buffers
	 * @param bufferSize
	 *            size that buffers should be (in case more are needed)
	 */
	public SoundIn(int sampleRate, int bufferSize,
			LinkedBlockingQueue<short[]> transferQueue,
			LinkedBlockingQueue<short[]> bufferPool) {

		this.transferQueue = transferQueue;
		this.bufferPool = bufferPool;
		this.bufferSize = bufferSize;

		try { // ... initialise

			recorder = new AudioRecord(AudioSource.MIC, sampleRate,
					AudioFormat.CHANNEL_IN_MONO,
					AudioFormat.ENCODING_PCM_16BIT, bufferSize * 2);

		} catch (Throwable x) {
			Log.v(TAG, "Error reading voice audio", x);
		}
	}

	// Thread methods

	/**
	 * The run loop: Get a buffer from the buffer pool, if there are any, else
	 * make a new buffer. Fill it from the capture line. Queue it up on the
	 * transfer queue. Repeat until our quit flag goes high.
	 */
	@Override
	public void run() {

		// The buffer we're moving from capture to transfer queue.
		short[] buffer = null;

		// Turn the line on. If we can't, we quit.
		try {
			on();
		} catch (Exception e) {
			quitFlag = true;
		}

		runLoop: while (true) {
			// First check if we've been told to quit.
			if (quitFlag)
				break runLoop;

			// Get a buffer from the buffer pool, if possible, else make a new
			// one.
			// The poll method doesn't block -- it returns null if nothing is
			// available.
			buffer = bufferPool.poll();
			if (buffer == null) {
				// Tried recycling -- now be greedy.
				buffer = new short[bufferSize];
			}

			// Read something from the capture line. Nothing fancy is done about
			// partial buffers. Might be good to zero the buffer before reading,
			// or wrap the short array in something that stores the actual
			// number
			// of samples read.
			recorder.read(buffer, 0, bufferSize);

			// Put it in the transfer queue
			putLoop: while (true) {
				try {
					transferQueue.put(buffer);
					break putLoop;
				} catch (InterruptedException e) {
					// Check for quit.
					if (quitFlag)
						break runLoop;
				}
			}

			// Forget we had this buffer.
			buffer = null;

			// Give the poor output thread a break.
			yield();
		}

		// We've been told to quit. Shut down the line and fall out of the loop.
		off();
	}

	// Other methods

	/**
	 * Turn on the quit flag. This is the officially sanctioned method for
	 * stopping this thread. Do not call stop() -- you have been warned.
	 */
	public void quit() {
		quitFlag = true;
	}

	/**
	 * Start the capture line.
	 */
	public void on() {
		recorder.startRecording();
	}

	/**
	 * Close the capture line.
	 */
	public void off() {
		// We aren't going to copy over any more info, so we flush anything the
		// line has in its internal buffers.
		recorder.stop();
		recorder.release();
	}
}
