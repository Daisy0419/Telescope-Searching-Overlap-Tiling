import ligo.skymap
import ligo.skymap.moc
import ligo.skymap.io.fits
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import os
import re


def convert_mac_healpix(input_fits_file, output_fits_file):
    map, header = ligo.skymap.io.fits.read_sky_map(input_fits_file, nest=True)

    target_nside = 64

    nside_sqr = len(map)/12
    nside = np.sqrt(nside_sqr).astype(int)
    assert nside * nside == nside_sqr, f"Imported map is not a valid nside!"
    print(f"nside: {nside}")

    scale = nside / target_nside
    scale = (scale * scale).astype(int)

    level = np.log2(scale).astype(int)
    assert 2**level == scale, f"nside {nside} note a power of 2 multiple of target_nside {target_nside}"

    print(f"Scaling down {scale}X")

    # map_reduced = np.array(map).reshape(-1,scale).sum(axis=-1)
    # map_ring = hp.pixelfunc.reorder(map_reduced, inp="NESTED", out="RING", n2r=True)
    # hp.write_map(output_fits_file, np.array(map_reduced), overwrite=True)
    m_reduced = np.asarray(map).reshape(-1, scale).sum(axis=-1)
    hp.write_map(output_fits_file, m_reduced, nest=True, overwrite=True)


def visualize_converted_map(output_fits_file):
    hmap = hp.read_map(output_fits_file)
    hp.mollview(hmap, cmap="viridis")
    fig = plt.gcf() 
    fig.set_size_inches(20, 12) 

    plt.show()

def visualize_converted_map_nested(output_fits_file):
    hmap = hp.read_map(output_fits_file, nest=True, verbose=False)
    hp.mollview(hmap, nest=True, cmap="viridis")
    plt.gcf().set_size_inches(20, 12)
    plt.show()



def get_all_files(path):
    files = []
    for entry in os.listdir(path):
        full_path = os.path.join(path, entry)
        if os.path.isfile(full_path):
            files.append(full_path)
    return files


def getEvents(folder_path):
    # Pattern to match and extract the GW event name
    pattern = re.compile(r"(GW\d{6}_\d{6})")

    # List to store extracted event IDs
    event_ids = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".fits"):
            match = pattern.search(filename)
            if match:
                event_ids.append(match.group(1))

    print(event_ids)
    return event_ids


path = "ligo_data"
# map_files = get_all_files(path)
# for file in map_files:
#     convert_mac_healpix(file, output_fits_file)
#     print(file)
#     visualize_converted_map(output_fits_file)


files = getEvents(path)

for file in files:
    input_fits_file = f"ligo_data/IGWN-GWTC3p0-v2-{file}_PEDataRelease_cosmo_reweight_C01_Mixed.fits"
    output_fits_file = f"fits_files/{file}.fits"
    convert_mac_healpix(input_fits_file, output_fits_file)
    visualize_converted_map_nested(output_fits_file)

