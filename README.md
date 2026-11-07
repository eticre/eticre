## Hi there 👋

# GPX Elevation Cleaner (FFT-based)

A Python command-line tool to remove elevation spikes and GPS noise from `.gpx` tracks using the Fast Fourier Transform (FFT). 

## 🧠 Why not use standard filters?

Traditional smoothing filters (like moving averages or median filters) often fail on GPS elevation data. When these filters "smooth" a spike, they create an artificial step in the altitude profile. Because the slope (gradient) is the mathematical derivative of altitude, the derivative of a step is an infinite spike. As a result, standard filters often replace a positive spike with a massive, fake negative spike in the gradient.

This tool solves the problem entirely by transforming the altitude data into the frequency domain. It filters out the high-frequency noise (GPS errors and spikes) while preserving the low-frequency terrain profile. The result is an inherently smooth, infinitely differentiable curve that yields realistic gradient calculations.

## ✨ Features

- **Spatial Resampling**: Interpolates the track at a fixed distance step (default 1 meter) to prevent mathematical explosions caused by points that are too close to each other.
- **Reflection Padding**: Prevents the "Gibbs phenomenon" (artificial edge spikes at the start and end of the track) by mirroring the signal before the FFT and truncating it afterward.
- **Recursive Low-Pass Filter**: Simulates a steep Butterworth filter by applying the frequency gain multiple times, resulting in an ultra-smooth curve.
- **Sliding Window Gradient**: Computes the slope percentage (`%`) over a customizable spatial window (e.g., every 50 meters) instead of point-by-point.
- **Visualization**: Generates a high-quality overlay plot of the original vs. cleaned gradient.
- **Export**: Optionally saves a new `.gpx` file with the corrected elevation data.

## 📦 Installation

You need Python 3.7+ and the following libraries:

```bash
pip install numpy matplotlib gpxpy

parameters:
Parameter	Default	            Description
input	(required)	            Path to the input GPX file.
output	(required)	            Path to save the output image (.png).
--step	1.0	                    Resampling distance in meters.
--cutoff-m	250.0	            FFT cutoff wavelength in meters. Variations shorter than this will be treated as noise and removed.
--passes	3	                Number of recursive filter passes. Higher values = steeper cutoff and smoother curve.
--slope-m	50.0	            Spatial window in meters used to calculate the slope (gradient). A larger window smooths out the  final gradient line.
--ylim-percentile	98.0	    Y-axis limit based on the percentile of the cleaned slope. Prevents massive original spikes from ruining the plot scale. Set to 0 to disable.
--gpx-clean	None	            (Optional) Path to save a new GPX file with the cleaned elevations.


basic usage:
only view
python smooth-gpx.py input.gpx output.png

generate image and new gpx file:
python smooth-gpx.py input.gpx output.png --gpx-clean output_cleaned.gpx

🔬How it works (Under the hood)
Interpolation: The track is resampled at regular intervals (e.g., every 1m) to ensure uniform distance between points.

Detrending: A linear trend (a straight line from start to finish) is subtracted from the elevation profile to prevent FFT edge artifacts.

Reflection Padding: The signal is mirrored at both ends (20% of its length) to make the signal continuous. This is the secret to eliminating edge spikes.

Frequency Domain Filtering: The FFT is computed. A Butterworth-like gain curve is applied, effectively zeroing out high frequencies (noise) while keeping low frequencies (actual hills). This is raised to the power of --passes for recursive smoothing.

Inverse FFT & Trend Addition: The signal is reconstructed, the trend is added back, and the padding is removed.

Slope Calculation: The gradient is calculated by comparing the elevation at i - slope-m/2 and i + slope-m/2, creating a smooth and realistic profile.