import os
from pathlib import Path
import numpy as np
import pandas as pd
import healpy as hp
from astropy.io import fits


def detect_healpix_ordering_and_nside(fits_path: str):
    """
    Returns (nside, nest_bool, prob_hdu_index_or_name, prob_colname_or_None, is_image_map)
    Tries common GW skymap layouts:
      - HEALPix map stored as IMAGE in an extension named 'PROB' or primary HDU
      - HEALPix map stored as BINTABLE with a 'PROB' column
    """
    with fits.open(fits_path, memmap=True) as hdul:
        # Try to find an extension named 'PROB'
        prob_hdu = None
        for i, h in enumerate(hdul):
            if (h.name or "").upper() == "PROB":
                prob_hdu = i
                break

        # If not found, fall back to first extension that looks like healpix
        if prob_hdu is None:
            # Common in ligo.skymap: probability in HDU 1 (BINTABLE)
            # We'll just use 1 if exists
            prob_hdu = 1 if len(hdul) > 1 else 0

        hdr = hdul[prob_hdu].header

        # ORDERING can appear in various places
        ordering = (hdr.get("ORDERING") or hdr.get("INDXSCHM") or "").strip().upper()
        nest = True if ordering == "NESTED" else False if ordering == "RING" else None

        # NSIDE might be present
        nside = hdr.get("NSIDE")

        # Determine if image or table
        is_image = hdul[prob_hdu].data is not None and hdul[prob_hdu].header.get("XTENSION") != "BINTABLE"
        if hdul[prob_hdu].header.get("XTENSION") == "BINTABLE":
            is_image = False

        # If BINTABLE, find a probability-like column
        prob_col = None
        if not is_image and hdul[prob_hdu].data is not None:
            cols = [c.upper() for c in hdul[prob_hdu].columns.names]
            for candidate in ["PROB", "PROBABILITY"]:
                if candidate in cols:
                    prob_col = candidate
                    break

        return nside, nest, prob_hdu, prob_col, is_image


def read_prob_map(fits_path: str):
    """
    Returns (prob, nest_bool) where prob is a 1D array of length npix.
    """
    nside_hdr, nest_hdr, prob_hdu, prob_col, is_image = detect_healpix_ordering_and_nside(fits_path)

    # healpy can usually read both image maps and table maps via hp.read_map
    # but specifying nest helps avoid wrong RA/Dec mapping.
    # If nest_hdr is None, hp.read_map will still read; we then infer later if needed.
    prob = hp.read_map(
        fits_path,
        hdu=prob_hdu,
        field=0,
        nest=nest_hdr if nest_hdr is not None else False,
        verbose=False
    )
    prob = np.asarray(prob, dtype=np.float64)

    # Determine ordering used for pix->ang conversion:
    # If header had it, trust it; otherwise assume RING (common default) unless you know otherwise.
    nest = nest_hdr if nest_hdr is not None else False

    # Sanity: infer nside from length if header missing
    npix = prob.size
    nside = hp.npix2nside(npix)
    if nside_hdr is not None and int(nside_hdr) != int(nside):
        print(f"[warn] {fits_path}: header NSIDE={nside_hdr} but inferred nside={nside}. Using inferred.")

    return prob, nest


def fits_to_csv_one(fits_path: str, csv_path: str):
    prob, nest = read_prob_map(fits_path)
    npix = prob.size
    nside = hp.npix2nside(npix)

    pix_id = np.arange(npix, dtype=np.int64)
    theta, phi = hp.pix2ang(nside, pix_id, nest=nest)

    ra_deg = np.degrees(phi)
    dec_deg = 90.0 - np.degrees(theta)

    df = pd.DataFrame({
        "pix_id": pix_id,
        "ra_deg": ra_deg,
        "dec_deg": dec_deg,
        "probability": prob
    })
    df.to_csv(csv_path, index=False)
    print(f"Wrote {csv_path} (npix={npix}, nside={nside}, ordering={'NESTED' if nest else 'RING'})")

def fits_to_csv_one(fits_path: str, csv_path: str):
    prob, nest = read_prob_map(fits_path)
    npix = prob.size
    nside = hp.npix2nside(npix)

    pix_id = np.arange(npix, dtype=np.int64)
    theta, phi = hp.pix2ang(nside, pix_id, nest=nest)

    ra_deg = np.degrees(phi)
    dec_deg = 90.0 - np.degrees(theta)

    df = pd.DataFrame({
        "pix_id": pix_id,
        "ra_deg": ra_deg,
        "dec_deg": dec_deg,
        "probability": prob
    })

    ordering_str = "NESTED" if nest else "RING"

    # Write a metadata line first, then the CSV header+data
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(f"# nside {nside} {ordering_str}\n")
        df.to_csv(f, index=False)

    print(f"Wrote {csv_path} (npix={npix}, nside={nside}, ordering={ordering_str})")


def convert_folder(input_dir: str, output_dir: str, pattern="*.fits"):
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # include .fits and .fits.gz
    fits_files = sorted(list(in_path.glob("*.fits")) + list(in_path.glob("*.fits.gz")))
    if not fits_files:
        print(f"No FITS files found in {input_dir}")
        return

    for f in fits_files:
        out_csv = out_path / (f.stem.replace(".fits", "") + "_pixels.csv")
        fits_to_csv_one(str(f), str(out_csv))


if __name__ == "__main__":
    convert_folder(
        input_dir="fits_files",
        output_dir="csv_files"
    )
