import os
from pathlib import Path
import numpy as np
import healpy as hp
from astropy.io import fits


def detect_healpix_ordering_and_nside(fits_path: str):
    """
    Returns (nside, nest_bool, prob_hdu_index, prob_colname_or_None, is_image_map)
    """
    with fits.open(fits_path, memmap=True) as hdul:
        prob_hdu = None
        for i, h in enumerate(hdul):
            if (h.name or "").upper() == "PROB":
                prob_hdu = i
                break

        if prob_hdu is None:
            prob_hdu = 1 if len(hdul) > 1 else 0

        hdr = hdul[prob_hdu].header
        ordering = (hdr.get("ORDERING") or hdr.get("INDXSCHM") or "").strip().upper()
        nest = True if ordering == "NESTED" else False if ordering == "RING" else None

        nside = hdr.get("NSIDE")

        is_image = hdul[prob_hdu].data is not None and hdr.get("XTENSION") != "BINTABLE"
        if hdr.get("XTENSION") == "BINTABLE":
            is_image = False

        prob_col = None
        if not is_image and hdul[prob_hdu].data is not None and hasattr(hdul[prob_hdu], "columns"):
            cols = [c.upper() for c in hdul[prob_hdu].columns.names]
            for candidate in ["PROB", "PROBABILITY"]:
                if candidate in cols:
                    prob_col = candidate
                    break

        return nside, nest, prob_hdu, prob_col, is_image


def read_prob_map(fits_path: str):
    """
    Returns (prob, nest_bool)
    """
    nside_hdr, nest_hdr, prob_hdu, prob_col, is_image = detect_healpix_ordering_and_nside(fits_path)

    # For most GW FITS handled by healpy, field=0 works for probability.
    prob = hp.read_map(
        fits_path,
        hdu=prob_hdu,
        field=0,
        nest=nest_hdr if nest_hdr is not None else False,
        verbose=False
    )
    prob = np.asarray(prob, dtype=np.float64)

    nest = nest_hdr if nest_hdr is not None else False

    npix = prob.size
    nside = hp.npix2nside(npix)
    if nside_hdr is not None and int(nside_hdr) != int(nside):
        print(f"[warn] {fits_path}: header NSIDE={nside_hdr} but inferred nside={nside}. Using inferred.")

    return prob, nest


def fits_to_txt_one(fits_path: str, txt_path: str, float_decimals: int = 3):
    prob, nest = read_prob_map(fits_path)
    npix = prob.size
    nside = hp.npix2nside(npix)

    pix_id = np.arange(npix, dtype=np.int64)
    theta, phi = hp.pix2ang(nside, pix_id, nest=nest)

    lon_deg = np.degrees(phi)                  # RA / longitude
    lat_deg = 90.0 - np.degrees(theta)         # Dec / latitude

    # You requested header to say "NEST" not "NESTED"
    ordering_str = "NEST" if nest else "RING"

    # Write space-separated text
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write(f"# nside {nside} {ordering_str}\n")
        f.write("# pix_id lat lon prob\n")

        # Formatting:
        # pix_id as int
        # lat/lon with fixed decimals 
        # prob in scientific notation
        lat_fmt = f"{{:.{float_decimals}f}}"
        lon_fmt = f"{{:.{float_decimals}f}}"

        for i in range(npix):
            f.write(
                f"{int(pix_id[i])} "
                f"{lat_fmt.format(lat_deg[i])} "
                f"{lon_fmt.format(lon_deg[i])} "
                f"{prob[i]:.16e}\n"
            )

    print(f"Wrote {txt_path} (npix={npix}, nside={nside}, ordering={ordering_str})")


def convert_folder_to_txt(input_dir: str, output_dir: str):
    in_path = Path(input_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    fits_files = sorted(list(in_path.glob("*.fits")) + list(in_path.glob("*.fits.gz")))
    if not fits_files:
        print(f"No FITS files found in {input_dir}")
        return

    for f in fits_files:
        # make a nice name even if it was .fits.gz
        stem = f.name.replace(".fits.gz", "").replace(".fits", "")
        out_txt = out_path / f"{stem}_pixels.txt"
        fits_to_txt_one(str(f), str(out_txt))


if __name__ == "__main__":
    convert_folder_to_txt(
        input_dir="fits_files",
        output_dir="txt_files"
    )
