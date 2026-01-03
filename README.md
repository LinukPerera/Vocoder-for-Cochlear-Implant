### **Project Title:** Shannon Vocoder for Cochlear Implant Simulation

### **Project Description:**

This project implements an **18-channel Shannon vocoder** in MATLAB to simulate cochlear implant sound processing. The vocoder encodes speech signals into frequency bands, extracts amplitude envelopes, and modulates noise carriers to mimic the way cochlear implants stimulate the auditory nerve. The project not only generates vocoded speech but also provides a **comprehensive analysis** of every stage in the processing pipeline.

### **Key Objectives:**

1. To implement an **18-channel bandpass filter bank** based on Shannon’s frequency mapping.
2. To extract **temporal envelopes** from each band using half-wave rectification and low-pass filtering.
3. To modulate white noise carriers with the extracted envelopes to simulate cochlear implant stimulation.
4. To perform detailed analysis, including waveform, spectrogram, frequency response, and filter pole-zero visualizations.
5. To ensure robustness by handling Nyquist limits, frequency clamping, and filter stability.

### **MATLAB Implementation Highlights:**

* **Input Signal Preprocessing:**

  * Audio resampled to 22.05 kHz and converted to mono for consistency.

* **Filter Bank Design:**

  * 18 channels with logarithmically spaced cut-off frequencies using Shannon’s formula.
  * 3rd-order elliptical bandpass filters for analysis.
  * Clamping of filter frequencies to prevent exceeding Nyquist limit.

* **Envelope Extraction:**

  * Half-wave rectification followed by a 1st-order low-pass filter at 160 Hz.
  * Captures speech rhythm and intelligibility cues relevant to cochlear implant users.

* **Noise Carrier Modulation:**

  * Each envelope modulates a white noise carrier, simulating auditory nerve stimulation.
  * Output signals from all channels are summed to form the final vocoded signal.

* **Analysis & Visualization:**

  * Time-domain and frequency-domain plots for input and output.
  * Step-by-step channel analysis (example: Channel 8) showing BPF, rectified signal, envelope, modulation, and final output.
  * Filter bank characteristics: cut-off frequencies, bandwidths, and magnitude responses.
  * Comparison of filter orders (Shannon 3rd-order, Rosen 6th-order, optimal ellipord).
  * Low-pass envelope filter response verification.
  * Pole-zero plot for a representative analysis filter.

* **Output:**

  * Vocoded speech played back using `soundsc`.
  * Visual and numerical analysis suitable for report submission.

### **Applications:**

* Simulating cochlear implant hearing for research and education.
* Testing and visualizing speech processing algorithms.
* Understanding temporal and spectral processing in auditory prostheses.
