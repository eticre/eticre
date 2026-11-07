#!/usr/bin/env python3
"""
Pulizia definitiva GPX via FFT con padding e passaggi ricorsivi.
Esporta il grafico di pendenza ed eventualmente il file GPX pulito.
"""

import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import gpxpy
except ImportError:
    print("Manca gpxpy. Installa con: pip install gpxpy", file=sys.stderr)
    sys.exit(1)

def haversine(lat1, lon1, lat2, lon2):
    R = 6371000.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp = p2 - p1
    dl = np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))

def load_gpx(path):
    with open(path, "r", encoding="utf-8") as f:
        gpx = gpxpy.parse(f)
    lats, lons, eles = [], [], []
    for track in gpx.tracks:
        for seg in track.segments:
            for pt in seg.points:
                if pt.elevation is not None:
                    lats.append(pt.latitude)
                    lons.append(pt.longitude)
                    eles.append(pt.elevation)
    return np.array(lats), np.array(lons), np.array(eles)

def cumulative_distance_m(lats, lons):
    d = np.zeros(len(lats))
    if len(lats) > 1:
        d[1:] = haversine(lats[:-1], lons[:-1], lats[1:], lons[1:])
    return np.cumsum(d)

def resample_track(lats, lons, eles, dist_orig, step=1.0):
    max_dist = dist_orig[-1]
    dist_uniform = np.arange(0, max_dist, step)
    if dist_uniform[-1] < max_dist:
        dist_uniform = np.append(dist_uniform, max_dist)
        
    lat_u = np.interp(dist_uniform, dist_orig, lats)
    lon_u = np.interp(dist_uniform, dist_orig, lons)
    ele_u = np.interp(dist_uniform, dist_orig, eles)
    return lat_u, lon_u, ele_u, dist_uniform

def fft_lowpass_padded(ele, step_m, cutoff_period_m=250.0, passes=3):
    n = len(ele)
    if n < 2 or cutoff_period_m <= 0 or passes <= 0:
        return ele.copy()

    # Padding Riflesso per evitare artefatti ai bordi
    pad_len = max(10, int(n * 0.2))
    padded_ele = np.pad(ele, pad_len, mode='reflect')
    n_padded = len(padded_ele)

    # Rimozione trend lineare
    t = np.arange(n_padded)
    p = np.polyfit(t, padded_ele, 1)
    trend = np.polyval(p, t)
    detrended = padded_ele - trend

    # FFT e filtro Butterworth
    fft_vals = np.fft.fft(detrended)
    freqs = np.fft.fftfreq(n_padded, d=step_m)

    f_cutoff = 1.0 / cutoff_period_m
    gain = 1.0 / (1.0 + (np.abs(freqs) / f_cutoff)**4)
    gain = gain ** passes

    # Antitrasformata
    fft_filtered = fft_vals * gain
    filtered_detrended = np.real(np.fft.ifft(fft_filtered))
    filtered_padded = filtered_detrended + trend

    # Rimozione del padding
    return filtered_padded[pad_len:-pad_len]

def compute_slope_percent(ele, dist, window_m=50.0, step=1.0):
    half_pts = int((window_m / 2.0) / step)
    if half_pts < 1:
        half_pts = 1
        
    n = len(ele)
    slope = np.full(n, np.nan)
    
    for i in range(n):
        lo = max(0, i - half_pts)
        hi = min(n - 1, i + half_pts)
        d_dist = dist[hi] - dist[lo]
        if d_dist > 0:
            slope[i] = ((ele[hi] - ele[lo]) / d_dist) * 100.0
    return slope

def write_clean_gpx(src_path, dst_path, dist_orig, ele_uniform_cleaned, dist_uniform):
    """Interpola l'altitudine pulita sulla traccia originale e salva il file."""
    ele_cleaned_orig = np.interp(dist_orig, dist_uniform, ele_uniform_cleaned)
    with open(src_path, "r", encoding="utf-8") as f:
        gpx = gpxpy.parse(f)
    
    idx = 0
    for track in gpx.tracks:
        for seg in track.segments:
            for pt in seg.points:
                if pt.elevation is not None:
                    pt.elevation = float(ele_cleaned_orig[idx])
                    idx += 1
                    
    with open(dst_path, "w", encoding="utf-8") as f:
        f.write(gpx.to_xml())

def main():
    ap = argparse.ArgumentParser(description="Pulizia GPX definitiva con FFT")
    ap.add_argument("input",  help="File GPX di ingresso")
    ap.add_argument("output", help="File immagine di uscita (PNG)")
    ap.add_argument("--step",     type=float, default=1.0,  help="Passo interpolazione (m)")
    ap.add_argument("--cutoff-m", type=float, default=250.0, help="Periodo di taglio FFT in metri (default 250)")
    ap.add_argument("--passes",   type=int,   default=3,    help="Passaggi ricorsivi (default 3)")
    ap.add_argument("--slope-m",  type=float, default=50.0, help="Finestra calcolo pendenza (m, default 50)")
    ap.add_argument("--ylim-percentile", type=float, default=98.0, help="Limite asse Y (0 = nessun limite)")
    ap.add_argument("--gpx-clean", metavar="FILE", help="(opzionale) salva il GPX con le altitudini pulite")
    args = ap.parse_args()

    # 1. Caricamento
    lats, lons, eles = load_gpx(args.input)
    if len(eles) < 3:
        print("Punti insufficienti nel file GPX.", file=sys.stderr)
        sys.exit(2)
        
    dist_orig = cumulative_distance_m(lats, lons)
    lats_u, lons_u, eles_u, dist_u = resample_track(lats, lons, eles, dist_orig, step=args.step)

    # 2. Filtro FFT
    cleaned_u = fft_lowpass_padded(
        eles_u, step_m=args.step, 
        cutoff_period_m=args.cutoff_m, 
        passes=args.passes
    )

    # 3. Calcolo pendenze
    slope_orig_u    = compute_slope_percent(eles_u,    dist_u, window_m=args.slope_m, step=args.step)
    slope_cleaned_u = compute_slope_percent(cleaned_u, dist_u, window_m=args.slope_m, step=args.step)
    dist_km = dist_u / 1000.0

    residual = eles_u - cleaned_u
    outliers_u = np.abs(residual) > 5.0

    # 4. Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(dist_km, slope_orig_u, color="#1f77b4", lw=1.0, alpha=0.4, label="Originale")
    ax.plot(dist_km, slope_cleaned_u, color="#d62728", lw=2.0, 
            label=f"Pulita (FFT {args.passes} passaggi, cutoff {args.cutoff_m}m)")

    if outliers_u.any():
        ax.scatter(dist_km[outliers_u], np.zeros(np.sum(outliers_u)),
                   color="orange", marker="|", s=100, zorder=5, 
                   label=f"Spike rilevati ({int(outliers_u.sum())})")

    ax.axhline(0, color="black", lw=0.6, alpha=0.5)
    ax.set_xlabel("Distanza (km)")
    ax.set_ylabel(f"Pendenza (%) calcolata su {args.slope_m}m")
    ax.set_title("Pendenza percentuale: Originale vs Pulita (FFT con Padding)")
    ax.grid(True, alpha=0.3)

    if args.ylim_percentile and 0 < args.ylim_percentile < 100:
        finite_vals = slope_cleaned_u[np.isfinite(slope_cleaned_u)]
        if finite_vals.size:
            p_hi = np.percentile(np.abs(finite_vals), args.ylim_percentile)
            ax.set_ylim(-max(p_hi * 1.2, 10.0), max(p_hi * 1.2, 10.0))

    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    plt.close(fig)

    print(f"Grafico salvato in: {args.output}")
    
    # 5. Salvataggio GPX (opzionale)
    if args.gpx_clean:
        write_clean_gpx(args.input, args.gpx_clean, dist_orig, cleaned_u, dist_u)
        print(f"GPX pulito salvato in: {args.gpx_clean}")

if __name__ == "__main__":
    main()