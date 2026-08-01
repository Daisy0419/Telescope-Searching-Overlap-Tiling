import matplotlib.pyplot as plt
import healpy as hp
import matplotlib.ticker as mticker


def visualize_converted_map(fits_file, save_path=None):
    hmap = hp.read_map(fits_file, verbose=False)
    # hp.mollview(hmap, cmap="viridis")
    hp.mollview(hmap,
                cmap="viridis",
                title='',
                notext=True,
                # norm='hist',
                cbar=True)

    fig = plt.gcf()
    fig.set_size_inches(3, 4)  # Paper-friendly size
    # fig.set_size_inches(2, 3)

    if len(fig.axes) > 1:
        cbar_ax = fig.axes[1]  # healpy creates it as the second axis
        cbar_ax.set_xlabel("Probability", fontsize=13, labelpad=20)
        cbar_ax.xaxis.set_label_position('bottom')
        cbar_ax.xaxis.set_ticks_position('bottom')
        cbar_ax.tick_params(labelsize=12)
        cbar_ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.4f'))

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()



# file = "GW200311_115853.fits"
# save_path = "map_GW200311_115853.png"

# # file = "GW200216_220804.fits"
# # save_path = "map_GW200216_220804.png"
# visualize_converted_map(file, save_path)

skymaps = ['GW200129_065458', 'GW200225_060421', 'GW200220_124850', 'GW200322_091133', 'GW200302_015811', 
           'GW191216_213338', 'GW191129_134029', 'GW200316_215756', 'GW200128_022011', 'GW200306_093714', 
           'GW191204_110529', 'GW191105_143521', 'GW200210_092254', 'GW200112_155838', 'GW191103_012549', 
           'GW200105_162426', 'GW191127_050227', 'GW200219_094415', 'GW200224_222234', 'GW200209_085452', 
           'GW191126_115259', 'GW191230_180458', 'GW200115_042309', 'GW191109_010717', 'GW191222_033537', 
           'GW200202_154313', 'GW191113_071753', 'GW200216_220804', 'GW200208_222617', 'GW191215_223052', 
           'GW200308_173609', 'GW200311_115853', 'GW200208_130117', 'GW230529_181500', 'GW200220_061928', 
           'GW191219_163120', 'GW191204_171526']

for map in skymaps:
    file = f"{map}.fits"
    save_path = f"map_{map}.png"
    visualize_converted_map(file, save_path)

