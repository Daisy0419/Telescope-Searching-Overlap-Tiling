import math
import healpy as hp
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.coordinates import ICRS
from typing import Literal
from astropy_healpix import HEALPix
from regions import (
    CircleSkyRegion,
    PolygonSkyRegion,
    RectangleSkyRegion,
    Region,
    Regions,
)

def _generate_rectangular_tiling_core(
    width: u.Quantity,
    height: u.Quantity,
    *,
    overlap: float,
    stagger: bool,
    pole_eps: float,
    ra_reference: Literal["center", "equatorward_edge"],
    cos_floor: float,
):

    """
    Generate centers for a grid of RA/Dec-aligned rectangular pointings.

    The center-to-center spacings are
        declination spacing = height * (1 - overlap)
        local RA spacing    = width  * (1 - overlap)

    The reference declination can be selected in two ways:
    ``"center"``
        Evaluate the cosine correction at the center of the row. 

    ``"equatorward_edge"``
        Evaluate the correction at the edge of each tile that is closer
        to the equator. 

    Parameters
    ----------
    width : astropy.units.Quantity
        Angular width of one rectangular field of view.

    height : astropy.units.Quantity
        Angular height of one rectangular field of view.

    overlap : float
        Fraction by which the nominal center spacing is reduced in both
        directions. It must satisfy ``0 <= overlap < 1``.

    stagger : bool
        If True, shift every other declination row by half of its target
        RA spacing, producing a brick-like arrangement.

    pole_eps : float
        Small angular tolerance in degrees used when positioning the
        extreme declination rows.

    ra_reference : {"center", "equatorward_edge"}
        Latitude at which to calculate the cosine correction for the
        target RA spacing.

    cos_floor : float
        Positive lower bound placed on ``cos(reference_dec)`` to avoid
        division by zero near the poles.

    Returns
    -------
    tile_indices : list[int]
        Consecutive zero-based tile identifiers.

    tile_coords : list[astropy.coordinates.SkyCoord]
        ICRS coordinates of the tile centers. Coordinates are ordered
        row by row, from south to north, and then by RA within each row.

    """

    if not width.unit.is_equivalent(u.deg):
        raise ValueError("width must be an angular quantity.")
    if not height.unit.is_equivalent(u.deg):
        raise ValueError("height must be an angular quantity.")
    if not 0.0 <= overlap < 1.0:
        raise ValueError("overlap must be in [0, 1).")

    width_deg = float(width.to_value(u.deg))
    height_deg = float(height.to_value(u.deg))

    if width_deg <= 0.0 or height_deg <= 0.0:
        raise ValueError("width and height must be positive.")

    width_step = width_deg * (1.0 - overlap)
    height_step = height_deg * (1.0 - overlap)

    dec_min = -90.0 + height_deg / 2.0 - pole_eps
    dec_max = +90.0 - height_deg / 2.0 + pole_eps

    if dec_min > dec_max + 1e-12:
        return [], []

    num_dec = (
        int(np.floor((dec_max + 1e-12 - dec_min) / height_step))
        + 1
    )
    dec_centers = dec_min + np.arange(num_dec) * height_step

    ra_rows = []
    dec_rows = []

    for row_index, dec in enumerate(dec_centers):
        if ra_reference == "equatorward_edge":
            if dec >= 0.0:
                reference_dec = dec - height_deg / 2.0
            else:
                reference_dec = dec + height_deg / 2.0

            reference_dec = np.clip(
                reference_dec,
                -90.0 + pole_eps,
                +90.0 - pole_eps,
            )
        else:
            reference_dec = dec

        cos_dec = max(
            cos_floor,
            math.cos(math.radians(reference_dec)),
        )
        target_ra_step = width_step / cos_dec

        ra_offset = (
            0.5 * target_ra_step
            if stagger and row_index % 2 == 1
            else 0.0
        )

        if target_ra_step >= 360.0 - 1e-9:
            ra_centers = np.array([ra_offset % 360.0])
        else:
            num_ra = int(np.ceil(360.0 / target_ra_step))
            base = np.linspace(
                0.0,
                360.0,
                num_ra,
                endpoint=False,
            )
            ra_centers = (base + ra_offset) % 360.0

        ra_rows.append(ra_centers)
        dec_rows.append(
            np.full(ra_centers.size, dec, dtype=float)
        )

    all_ra = np.concatenate(ra_rows)
    all_dec = np.concatenate(dec_rows)

    coordinates = SkyCoord(
        ra=all_ra * u.deg,
        dec=all_dec * u.deg,
        frame=ICRS(),
    )

    tile_indices = list(range(all_ra.size))
    tile_coords = list(coordinates)

    return tile_indices, tile_coords


def generate_rectangular_tiling(
    width: u.Quantity,
    height: u.Quantity,
    stagger: bool = True,
    pole_eps: float = 1e-4,
):
    """Center-latitude RA spacing with no intentional overlap."""
    return _generate_rectangular_tiling_core(
        width,
        height,
        overlap=0.0,
        stagger=stagger,
        pole_eps=pole_eps,
        ra_reference="center",
        cos_floor=1e-9,
    )


def generate_rectangular_tiling_overlap(
    width: u.Quantity,
    height: u.Quantity,
    overlap: float = 0.12,
    stagger: bool = True,
    pole_eps: float = 1e-5,
):
    """Center-latitude RA spacing with configurable overlap."""
    return _generate_rectangular_tiling_core(
        width,
        height,
        overlap=overlap,
        stagger=stagger,
        pole_eps=pole_eps,
        ra_reference="center",
        cos_floor=1e-6,
    )


def generate_rectangular_tiling_equatorward_edge_ra_spacing(
    width: u.Quantity,
    height: u.Quantity,
    stagger: bool = True,
    pole_eps: float = 1e-4,
):
    """RA spacing evaluated at each row's equatorward edge."""
    return _generate_rectangular_tiling_core(
        width,
        height,
        overlap=0.0,
        stagger=stagger,
        pole_eps=pole_eps,
        ra_reference="equatorward_edge",
        cos_floor=1e-9,
    )



def generate_healpix_pixel_centers(
    nside: int = 64,
    order: str = "nested",
):
    """
    Return all HEALPix pixel centers at the given NSIDE.

    Parameters
    ----------
    nside : int
        HEALPix NSIDE.
    order : {"nested","ring"}
        HEALPix ordering.

    Returns
    -------
    tile_ids : list[int]
        HEALPix pixel indices (0..12*nside^2-1) in the chosen ordering.
    tile_coords_list : list[SkyCoord]
        Pixel-center coordinates (ICRS).
    """
    if order not in ("nested", "ring"):
        raise ValueError("`order` must be 'nested' or 'ring'.")
    if not (isinstance(nside, int) and nside > 0):
        raise ValueError("`nside` must be a positive integer.")

    hpx = HEALPix(nside=nside, order=order, frame=ICRS())
    npix = hpx.npix

    ipix = np.arange(npix, dtype=np.int64)
    coords = hpx.healpix_to_skycoord(ipix)  # SkyCoord array

    tile_ids = ipix.tolist()
    tile_coords_list = list(coords)  # list[SkyCoord], one per pixel

    return tile_ids, tile_coords_list