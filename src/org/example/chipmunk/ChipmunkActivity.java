package org.example.chipmunk;

import java.util.concurrent.LinkedBlockingQueue;

import ptressel.soundutils.Effects;
import ptressel.soundutils.SoundIn;
import ptressel.soundutils.SoundOut;
import android.app.Activity;
import android.media.AudioFormat;
import android.media.AudioRecord;
import android.os.Bundle;

public class ChipmunkActivity extends Activity {

	private static final int SAMPLE_RATE = 8000;
	private int bufferSize;
	/** Capture thread */
	private SoundIn soundIn;

	/** Playback thread */
	private SoundOut soundOut;

	/**
	 * Buffer pool: Input side gets buffers to fill from here unless there are
	 * no more, in which case it allocates one. Output side puts spent buffers
	 * back in. If all goes well, the length of this queue should not grow
	 * indefinitely. The intent of the buffer pool is to reuse buffers rather
	 * than allocating new ones for each input, then throwing them away.
	 * Avoiding object creation and destruction is usually the most significant
	 * optimization in this sort of application. The buffer pool is prefilled
	 * with a fair number of buffers.
	 */
	private LinkedBlockingQueue<short[]> transferBufferPool = new LinkedBlockingQueue<short[]>();

	/**
	 * Transfer queue: Input side puts filled buffers on one end. Output side
	 * dequeues them from the other.
	 */
	private LinkedBlockingQueue<short[]> transferQueue = new LinkedBlockingQueue<short[]>();
	private Effects effects;

	/** Called when the activity is first created. */
	@Override
	public void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.main);

		// Make our in and out threads. We pass both of them the transfer queue
		// and buffer pool, along with their separate lines, debug info, etc.

		bufferSize = AudioRecord.getMinBufferSize(8000,
				AudioFormat.CHANNEL_IN_MONO, AudioFormat.ENCODING_PCM_16BIT);
		effects = new Effects(bufferSize);
		soundIn = new SoundIn(SAMPLE_RATE, bufferSize, transferQueue,
				transferBufferPool);
		soundOut = new SoundOut(SAMPLE_RATE, transferQueue, transferBufferPool,
				effects);
	}
}