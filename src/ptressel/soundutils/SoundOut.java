package ptressel.soundutils;

import java.util.concurrent.LinkedBlockingQueue;

import android.media.AudioFormat;
import android.media.AudioManager;
import android.media.AudioTrack;

/**
 * Get data from the transfer queue and write it to the output.
 * 
 * @see SoundCopy for more info.
 */

public class SoundOut extends Thread {

	// Instance data

	private Effects effects;

	/** Flag that tells us when to quit. */
	private boolean quitFlag = false;

	/** Transfer queue */
	private final LinkedBlockingQueue<short[]> transferQueue;

	/** Buffer pool */
	private final LinkedBlockingQueue<short[]> bufferPool;

	private AudioTrack track;

	// Constructor

	/**
	 * Make a SoundOut thread.
	 * 
	 * @param output
	 *            the output line to which to write data
	 * @param transferQueue
	 *            the queue from which full buffers will be obtained
	 * @param bufferPool
	 *            queue into which empty buffers will be placed
	 */
	public SoundOut(int sampleRate, LinkedBlockingQueue<short[]> transferQueue,
			LinkedBlockingQueue<short[]> bufferPool, Effects effects) {

		this.transferQueue = transferQueue;
		this.bufferPool = bufferPool;

		this.effects = effects;

		// http://www.badlogicgames.com/wordpress/?p=228

		int minSize = AudioTrack.getMinBufferSize(sampleRate,
				AudioFormat.CHANNEL_CONFIGURATION_MONO,
				AudioFormat.ENCODING_PCM_16BIT);
		track = new AudioTrack(AudioManager.STREAM_MUSIC, sampleRate,
				AudioFormat.CHANNEL_CONFIGURATION_MONO,
				AudioFormat.ENCODING_PCM_16BIT, minSize, AudioTrack.MODE_STREAM);
	}

	// Thread methods

	/**
	 * The run loop: Get a buffer from the transfer queue -- block til one is
	 * available. Write it to the output. Put the (spent) buffer in the buffer
	 * pool. Repeat until our quit flag goes high.
	 */
	@Override
	public void run() {

		// The buffer we're moving from transfer queue to output.
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

			// Get a buffer from the transfer queue. The take method blocks
			// until
			// an item is available.
			takeLoop: while (true) {
				try {
					buffer = transferQueue.take();
					break takeLoop;
				} catch (InterruptedException e) {
					// Check for quit.
					if (quitFlag)
						break runLoop;
				}
			}

			effects.processBuffer(buffer);

			// Write to the output.
			track.write(buffer, 0, buffer.length);

			// Put buffer back in the buffer pool. Note the SourceDataLine write
			// had better not return till it's copied out the buffer, 'cause
			// we're
			// about to reuse it...

			putLoop: while (true) {
				try {
					bufferPool.put(buffer);
					break putLoop;
				} catch (InterruptedException e) {
					// Check for quit.
					if (quitFlag)
						break runLoop;
				}
			}

			// Forget we had this buffer.
			buffer = null;

			// The input thread does a yield for us, so in fairness, we should
			// do
			// one for it. But we're the ones that can't keep up (under Linux
			// with Sun's Java). So consider being unfair. On the other hand,
			// it's the input side that we want to not miss anything. This side
			// we could always turn off when we don't want to hear it.
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
	 * Start the output line.
	 */
	public void on() {
		track.play();
	}

	/**
	 * Close the output line.
	 */
	public void off() {
		// Tell the output to drain its internal buffers, then close it.
		track.flush();
		track.release();
	}
}
