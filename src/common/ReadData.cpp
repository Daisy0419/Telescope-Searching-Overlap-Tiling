#include "ReadData.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <stdexcept>

// #include <highfive/H5File.hpp>
// #include <highfive/H5DataSet.hpp>
// #include <highfive/H5DataSpace.hpp>


//Read Tiling Data
void readTiles(const std::string& filename, 
             std::vector<int>& ranks, 
             std::vector<double>& ras, 
             std::vector<double>& decs, 
             std::vector<std::vector<int>>& member_pixels) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    // Skip the header line
    if (!std::getline(file, line)) {
        std::cerr << "Error reading header" << std::endl;
        return;
    }

    while (std::getline(file, line)) {
        if (line.empty()) continue; // Skip empty lines

        std::istringstream ss(line);
        std::string token;

        // Read ID (rank)
        if (!std::getline(ss, token, ',')) continue;
        int rank;
        try {
            rank = std::stoi(token);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing rank: " << token << std::endl;
            continue;
        }
        ranks.push_back(rank);

        // Read RA
        if (!std::getline(ss, token, ',')) continue;
        double ra;
        try {
            ra = std::stod(token);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing RA: " << token << std::endl;
            continue;
        }
        ras.push_back(ra);

        // Read DEC
        if (!std::getline(ss, token, ',')) continue;
        double dec;
        try {
            dec = std::stod(token);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing DEC: " << token << std::endl;
            continue;
        }
        decs.push_back(dec);

        // Read the remaining line as HEALPixels (space-separated)
        std::string heal_str;
        std::getline(ss, heal_str);
        std::istringstream heal_ss(heal_str);
        std::vector<int> pixels;
        int pixel;
        while (heal_ss >> pixel) {
            pixels.push_back(pixel);
        }
        member_pixels.push_back(pixels);
    }

    file.close();

}



//Read Probability Map - csv
bool readHEALPixelsCSV(const std::string& filename,
                       std::vector<double>& pixel_probs,
                       int nside) {
    const int pixel_count = 12 * nside * nside;
    pixel_probs.assign(pixel_count, 0.0);

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << "\n";
        return false;
    }

    auto trim = [](const std::string& s) -> std::string {
        size_t b = 0, e = s.size();
        while (b < e && std::isspace((unsigned char)s[b])) ++b;
        while (e > b && std::isspace((unsigned char)s[e - 1])) --e;
        return s.substr(b, e - b);
    };

    std::vector<unsigned char> seen(pixel_count, 0);
    std::string line;
    size_t lineno = 0, loaded = 0, dup = 0, oob = 0, bad = 0, nan_inf = 0;

    bool header_skipped = false;

    while (std::getline(file, line)) {
        ++lineno;
        line = trim(line);
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        // Skip exactly one header line
        if (!header_skipped) {
            header_skipped = true;
            continue;
        }

        // Parse CSV: pix_id,ra_deg,dec_deg,probability
        std::istringstream ss(line);
        std::string token;

        int pix_id;
        double ra, dec, prob;

        try {
            if (!std::getline(ss, token, ',')) throw std::runtime_error("missing pix_id");
            pix_id = std::stoi(trim(token));

            if (!std::getline(ss, token, ',')) throw std::runtime_error("missing ra");
            ra = std::stod(trim(token)); (void)ra;

            if (!std::getline(ss, token, ',')) throw std::runtime_error("missing dec");
            dec = std::stod(trim(token)); (void)dec;

            if (!std::getline(ss, token, ',')) throw std::runtime_error("missing prob");
            prob = std::stod(trim(token));
        } catch (...) {
            std::cerr << "Parse error at line " << lineno << ": " << line << "\n";
            ++bad;
            continue;
        }

        if (pix_id < 0 || pix_id >= pixel_count) {
            std::cerr << "Out-of-bounds pix_id " << pix_id
                      << " at line " << lineno << " (pixel_count=" << pixel_count << ")\n";
            ++oob;
            continue;
        }

        if (!std::isfinite(prob)) {
            std::cerr << "Non-finite prob for pix_id " << pix_id
                      << " at line " << lineno << " -> skipped\n";
            ++nan_inf;
            continue;
        }

        if (seen[pix_id]) ++dup;   // keep last value
        seen[pix_id] = 1;
        pixel_probs[pix_id] = prob;
        ++loaded;
    }

    const size_t seen_count = std::accumulate(seen.begin(), seen.end(), size_t{0});
    const size_t missing = (size_t)pixel_count > seen_count ? ((size_t)pixel_count - seen_count) : 0;
    const double sum_prob = std::accumulate(pixel_probs.begin(), pixel_probs.end(), 0.0);

    std::cerr << "Loaded " << loaded << " entries "
              << "(unique=" << seen_count
              << ", duplicates=" << dup
              << ", missing=" << missing
              << ", out-of-bound=" << oob
              << ", non-finite=" << nan_inf
              << ", parse_errors=" << bad << ")\n";
    std::cerr << "Sum(prob) = " << sum_prob << "\n";

    return true;
}



//Read Probability Map - txt(space seperated)
bool readHEALPixelsTXT(const std::string& filename,
                    std::vector<double>& pixel_probs,
                    int nside) {
    int pixel_count = 12 * nside * nside;
    pixel_probs.assign(pixel_count, 0.0);

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << "\n";
        return false;
    }

    std::vector<unsigned char> seen(pixel_count, 0);
    std::string line;
    size_t lineno = 0, loaded = 0, dup = 0, oob = 0, bad = 0, nan_inf = 0;

    while (std::getline(file, line)) {
        ++lineno;
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        std::istringstream ss(line);
        int pix_id;
        double lat, lon, prob;
        if (!(ss >> pix_id >> lat >> lon >> prob)) {
            std::cerr << "Parse error at line " << lineno << ": " << line << "\n";
            ++bad;
            continue;
        }

        if (pix_id < 0 || pix_id >= pixel_count) {
            std::cerr << "Out-of-bounds pix_id " << pix_id
                      << " at line " << lineno << " (pixel_count=" << pixel_count << ")\n";
            ++oob;
            continue;
        }

        if (!std::isfinite(prob)) {
            std::cerr << "Non-finite prob for pix_id " << pix_id
                      << " at line " << lineno << " -> skipped\n";
            ++nan_inf;
            continue;
        }

        if (seen[pix_id]) ++dup; // keep last value
        seen[pix_id] = 1;
        pixel_probs[pix_id] = prob;
        ++loaded;
    }

    const size_t seen_count = std::accumulate(seen.begin(), seen.end(), size_t{0});
    const size_t missing = pixel_count > seen_count ? (pixel_count - seen_count) : 0;
    const double finite_sum = std::accumulate(pixel_probs.begin(), pixel_probs.end(), 0.0);

    std::cerr << "Loaded " << loaded << " entries "
              << "(unique=" << seen_count
              << ", duplicates=" << dup
              << ", missing=" << missing
              << ", out bound=" << oob
              << ", non-finite=" << nan_inf
              << ", parse_errors=" << bad << ")\n";
    std::cerr << "Sum(prob) over finite entries = " << finite_sum << "\n";

    return true;
}


// bool readHEALPixelsHDF5_Compound(const std::string& filename,
//                                   std::vector<double>& pixel_probs,
//                                   int nside,
//                                   const std::string& dataset_name = "likelihood_map") {
//     int pixel_count = 12 * nside * nside;
//     pixel_probs.assign(pixel_count, 0.0);
//     try {
//         HighFive::File file(filename, HighFive::File::ReadOnly);

//         struct MapEntry {
//             int64_t pixel;
//             double lat;
//             double lon;
//             double probability;
//         };

//         auto dataset = file.getDataSet(dataset_name);
//         size_t n = dataset.getSpace().getDimensions()[0];
//         std::vector<MapEntry> entries(n);

//         hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(MapEntry));
//         H5Tinsert(memtype, "pixel",       offsetof(MapEntry, pixel),       H5T_NATIVE_INT64);
//         H5Tinsert(memtype, "lat",         offsetof(MapEntry, lat),         H5T_NATIVE_DOUBLE);
//         H5Tinsert(memtype, "lon",         offsetof(MapEntry, lon),         H5T_NATIVE_DOUBLE);
//         H5Tinsert(memtype, "probability", offsetof(MapEntry, probability), H5T_NATIVE_DOUBLE);

//         herr_t status = H5Dread(dataset.getId(), memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, entries.data());
//         H5Tclose(memtype);

//         if (status < 0) {
//             std::cerr << "H5Dread failed for dataset " << dataset_name << " in " << filename << "\n";
//             return false;
//         }

//         std::vector<unsigned char> seen(pixel_count, 0);
//         size_t loaded = 0, dup = 0, oob = 0, nan_inf = 0;

//         for (size_t i = 0; i < entries.size(); ++i) {
//             int pix_id = static_cast<int>(entries[i].pixel);
//             double prob = entries[i].probability;

//             if (pix_id < 0 || pix_id >= pixel_count) {
//                 std::cerr << "Out-of-bounds pix_id " << pix_id
//                           << " at index " << i << " (pixel_count=" << pixel_count << ")\n";
//                 ++oob;
//                 continue;
//             }
//             if (!std::isfinite(prob)) {
//                 std::cerr << "Non-finite prob for pix_id " << pix_id
//                           << " at index " << i << " -> skipped\n";
//                 ++nan_inf;
//                 continue;
//             }
//             if (seen[pix_id]) ++dup;
//             seen[pix_id] = 1;
//             pixel_probs[pix_id] = prob;
//             ++loaded;
//         }

//         const size_t seen_count = std::accumulate(seen.begin(), seen.end(), size_t{0});
//         const size_t missing = pixel_count > seen_count ? (pixel_count - seen_count) : 0;
//         const double finite_sum = std::accumulate(pixel_probs.begin(), pixel_probs.end(), 0.0);
//         std::cerr << "Loaded " << loaded << " entries "
//                   << "(unique=" << seen_count
//                   << ", duplicates=" << dup
//                   << ", missing=" << missing
//                   << ", out_of_bounds=" << oob
//                   << ", non-finite=" << nan_inf << ")\n";
//         std::cerr << "Sum(prob) over finite entries = " << finite_sum << "\n";
//         return true;

//     } catch (const HighFive::Exception& e) {
//         std::cerr << "HDF5 error reading file " << filename << ": " << e.what() << "\n";
//         return false;
//     }
// }



// //Read Probability Map - HDF5 (using HighFive)
// bool readHEALPixelsHDF5(const std::string& filename,
//                         std::vector<double>& pixel_probs,
//                         int nside,
//                         const std::string& pix_id_dataset,
//                         const std::string& prob_dataset) {
//     int pixel_count = 12 * nside * nside;
//     pixel_probs.assign(pixel_count, 0.0);

//     try {
//         HighFive::File file(filename, HighFive::File::ReadOnly);

//         // Read datasets
//         std::vector<int> pix_ids;
//         std::vector<double> probs;

//         file.getDataSet(pix_id_dataset).read(pix_ids);
//         file.getDataSet(prob_dataset).read(probs);

//         if (pix_ids.size() != probs.size()) {
//             std::cerr << "Error: pix_id and prob datasets have different sizes ("
//                       << pix_ids.size() << " vs " << probs.size() << ")\n";
//             return false;
//         }

//         std::vector<unsigned char> seen(pixel_count, 0);
//         size_t loaded = 0, dup = 0, oob = 0, nan_inf = 0;

//         for (size_t i = 0; i < pix_ids.size(); ++i) {
//             int pix_id = pix_ids[i];
//             double prob = probs[i];

//             if (pix_id < 0 || pix_id >= pixel_count) {
//                 std::cerr << "Out-of-bounds pix_id " << pix_id
//                           << " at index " << i << " (pixel_count=" << pixel_count << ")\n";
//                 ++oob;
//                 continue;
//             }
//             if (!std::isfinite(prob)) {
//                 std::cerr << "Non-finite prob for pix_id " << pix_id
//                           << " at index " << i << " -> skipped\n";
//                 ++nan_inf;
//                 continue;
//             }
//             if (seen[pix_id]) ++dup; // keep last value
//             seen[pix_id] = 1;
//             pixel_probs[pix_id] = prob;
//             ++loaded;
//         }

//         const size_t seen_count = std::accumulate(seen.begin(), seen.end(), size_t{0});
//         const size_t missing = pixel_count > seen_count ? (pixel_count - seen_count) : 0;
//         const double finite_sum = std::accumulate(pixel_probs.begin(), pixel_probs.end(), 0.0);

//         std::cerr << "Loaded " << loaded << " entries "
//                   << "(unique=" << seen_count
//                   << ", duplicates=" << dup
//                   << ", missing=" << missing
//                   << ", out_of_bounds=" << oob
//                   << ", non-finite=" << nan_inf << ")\n";
//         std::cerr << "Sum(prob) over finite entries = " << finite_sum << "\n";

//         return true;

//     } catch (const HighFive::Exception& e) {
//         std::cerr << "HDF5 error reading file " << filename << ": " << e.what() << "\n";
//         return false;
//     }
// }


// Helper to extract and lowercase the file extension
std::string getExtension(const std::string& filename) {
    size_t dot_pos = filename.rfind('.');
    if (dot_pos == std::string::npos || dot_pos == filename.length() - 1) {
        return "";
    }
    std::string ext = filename.substr(dot_pos + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return ext;
}

// map reading wrapper
bool readHEALPixels(const std::string& filename,
                    std::vector<double>& pixel_probs,
                    int nside) {
    std::string ext = getExtension(filename);

    if (ext == "csv") {
        return readHEALPixelsCSV(filename, pixel_probs, nside);
    } else if (ext == "txt" || ext == "dat") {
        return readHEALPixelsTXT(filename, pixel_probs, nside);
    // } else if (ext == "h5" || ext == "hdf5" || ext == "he5") {
    //     return readHEALPixelsHDF5(filename, pixel_probs, nside);
    //     // return readHEALPixelsHDF5_Compound(filename, pixel_probs,nside);
    } else {
        std::cerr << "Unknown file extension '" << ext
                  << "' for file: " << filename << "\n"
                  << "Supported extensions: .csv, .txt, .dat, .h5, .hdf5, .he5\n";
        return false;
    }
}



void computeTileProbs(const std::vector<std::vector<int>>& member_pixels,
                    const std::vector<double>& pixel_probs,
                    std::vector<double>& tile_probs) {
    tile_probs.resize(member_pixels.size(), 0.0);

    double sum = std::accumulate(pixel_probs.begin(), pixel_probs.end(), 0.0);
    std::cout << "init sum: " << sum << "\n";
    std::cout << "num pixels: " << pixel_probs.size() << "\n";


    for (size_t i = 0; i < member_pixels.size(); ++i) {
        tile_probs[i] = 0.0;
        for (int p : member_pixels[i]) {
                if(p > pixel_probs.size()) 
                    throw std::runtime_error("Error: pixel counts do not match tiles!");
                tile_probs[i] += pixel_probs[p];

        }
    }
}


int getTopTileIndices(std::vector<double>& tile_probs,
                      std::vector<int>& ranks,
                      std::vector<double>& ras,
                      std::vector<double>& decs,
                      std::vector<std::vector<int>>& member_pixels,
                      double coverage) {
    std::cout << "tile num: " << member_pixels.size() << "\n";
    std::cout << "total tile prob: " << std::accumulate(tile_probs.begin(), tile_probs.end(), 0.0) << "\n";
    if (tile_probs.empty()) return 0;

    // // Compute total sum to treat coverage as a fraction
    // double total_sum = 0.0;
    // for (double p : tile_probs) total_sum += p;
    // double threshold = coverage * total_sum;

    std::vector<size_t> indices(tile_probs.size());
    std::iota(indices.begin(), indices.end(), 0);

    // Sort indices in descending order of tile_probs
    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        return tile_probs[a] > tile_probs[b];
    });

    // Create sorted versions
    std::vector<double> sorted_tile_probs(tile_probs.size());
    std::vector<int> sorted_ranks(ranks.size());
    std::vector<double> sorted_ras(ras.size());
    std::vector<double> sorted_decs(decs.size());
    std::vector<std::vector<int>> sorted_member_pixels(member_pixels.size());

    for (size_t i = 0; i < indices.size(); ++i) {
        sorted_tile_probs[i] = tile_probs[indices[i]];
        sorted_ranks[i] = ranks[indices[i]];
        sorted_ras[i] = ras[indices[i]];
        sorted_decs[i] = decs[indices[i]];
        sorted_member_pixels[i] = std::move(member_pixels[indices[i]]);
    }

    // Assign back to originals
    tile_probs = std::move(sorted_tile_probs);
    ranks = std::move(sorted_ranks);
    ras = std::move(sorted_ras);
    decs = std::move(sorted_decs);
    member_pixels = std::move(sorted_member_pixels);

    // Accumulate until >= threshold using the now sorted tile_probs
    double eps = 10e-5;
    double cum_sum = 0.0;
    size_t n = 0;
    for (; n < tile_probs.size(); ++n) {
        // if (tile_probs[n] < eps) break;
        cum_sum += tile_probs[n];
        if (cum_sum >= coverage) break;
    }

    // Include the last one if we broke out
    if (n < tile_probs.size()) ++n;

    // Resize all to keep only the top N
    tile_probs.resize(n);
    ranks.resize(n);
    ras.resize(n);
    decs.resize(n);
    member_pixels.resize(n);
    std::cout << "total tile prob after trim: " << std::accumulate(tile_probs.begin(), tile_probs.end(), 0.0) << "\n";

    return static_cast<int>(n);
}

int getNonZeroTileIndices(std::vector<double>& tile_probs,
                      std::vector<int>& ranks,
                      std::vector<double>& ras,
                      std::vector<double>& decs,
                      std::vector<std::vector<int>>& member_pixels) {
    std::cout << "tile num: " << member_pixels.size() << "\n";
    std::cout << "total tile prob: " << std::accumulate(tile_probs.begin(), tile_probs.end(), 0.0) << "\n";
    if (tile_probs.empty()) return 0;

    // Find indices of tiles with prob > 0
    std::vector<size_t> keep;
    for (size_t i = 0; i < tile_probs.size(); ++i) {
        if (tile_probs[i] > 0.0) keep.push_back(i);
    }

    // Sort kept indices in descending order of tile_probs
    std::sort(keep.begin(), keep.end(), [&](size_t a, size_t b) {
        return tile_probs[a] > tile_probs[b];
    });

    size_t n = keep.size();

    // Create compacted versions
    std::vector<double> new_tile_probs(n);
    std::vector<int> new_ranks(n);
    std::vector<double> new_ras(n);
    std::vector<double> new_decs(n);
    std::vector<std::vector<int>> new_member_pixels(n);

    for (size_t i = 0; i < n; ++i) {
        new_tile_probs[i] = tile_probs[keep[i]];
        new_ranks[i] = ranks[keep[i]];
        new_ras[i] = ras[keep[i]];
        new_decs[i] = decs[keep[i]];
        new_member_pixels[i] = std::move(member_pixels[keep[i]]);
    }

    tile_probs = std::move(new_tile_probs);
    ranks = std::move(new_ranks);
    ras = std::move(new_ras);
    decs = std::move(new_decs);
    member_pixels = std::move(new_member_pixels);

    std::cout << "tiles with prob > 0: " << n << "\n";
    std::cout << "total tile prob after trim: " << std::accumulate(tile_probs.begin(), tile_probs.end(), 0.0) << "\n";

    return static_cast<int>(n);
}

// int getTopTileOverlap(std::vector<double>& tile_probs,
//                       std::vector<int>& ranks,
//                       std::vector<double>& ras,
//                       std::vector<double>& decs,
//                       std::vector<std::vector<int>>& member_pixels,
//                       const std::vector<double>& pixel_probs,
//                       double coverage) {
//     if (tile_probs.empty()) return 0;
//     std::cout << "tile num: " << member_pixels.size() << "\n";
//     std::cout << "total tile prob: " << std::accumulate(tile_probs.begin(), tile_probs.end(), 0.0) << "\n";

//     std::vector<size_t> indices(tile_probs.size());
//     std::iota(indices.begin(), indices.end(), 0);

//     // Sort indices in descending order of tile_probs
//     std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
//         return tile_probs[a] > tile_probs[b];
//     });

//     // Create sorted versions
//     std::vector<double> sorted_tile_probs(tile_probs.size());
//     std::vector<int> sorted_ranks(ranks.size());
//     std::vector<double> sorted_ras(ras.size());
//     std::vector<double> sorted_decs(decs.size());
//     std::vector<std::vector<int>> sorted_member_pixels(member_pixels.size());

//     for (size_t i = 0; i < indices.size(); ++i) {
//         sorted_tile_probs[i] = tile_probs[indices[i]];
//         sorted_ranks[i] = ranks[indices[i]];
//         sorted_ras[i] = ras[indices[i]];
//         sorted_decs[i] = decs[indices[i]];
//         sorted_member_pixels[i] = std::move(member_pixels[indices[i]]);
//     }

//     // Assign back to originals
//     tile_probs = std::move(sorted_tile_probs);
//     ranks = std::move(sorted_ranks);
//     ras = std::move(sorted_ras);
//     decs = std::move(sorted_decs);
//     member_pixels = std::move(sorted_member_pixels);

//     // Accumulate until >= coverage using the sorted tile_probs
//     std::vector<bool> visited(pixel_probs.size(), false);
//     double cum_sum = 0.0;
//     size_t n = 0;
//     for (; n < tile_probs.size(); ++n) {
//         for(int pix : member_pixels[n]) {
//             if(!visited[pix]) {
//                 cum_sum += pixel_probs[pix];
//                 visited[pix] = true;
//             }
//         }
//         if (cum_sum >= coverage) break;
//     }

//     // Include the last one if we broke out
//     if (n < tile_probs.size()) ++n;

//     // Resize all to keep only the top N
//     tile_probs.resize(n);
//     ranks.resize(n);
//     ras.resize(n);
//     decs.resize(n);
//     member_pixels.resize(n);
//     std::cout << "total tile prob after trim: " << std::accumulate(tile_probs.begin(), tile_probs.end(), 0.0) << "\n";

//     return static_cast<int>(n);
// }


// int getTopTileOverlapAwareGreedy(
//     std::vector<double>& tile_probs,
//     std::vector<int>& ranks,
//     std::vector<double>& ras,
//     std::vector<double>& decs,
//     std::vector<std::vector<int>>& member_pixels,
//     const std::vector<double>& pixel_probs,
//     double coverage ) {
//     const size_t T = tile_probs.size();
//     if (T == 0) return 0;

//     std::cout << "tile num: " << member_pixels.size() << "\n";
//     std::cout << "total tile prob (standalone sum): "
//               << std::accumulate(tile_probs.begin(), tile_probs.end(), 0.0) << "\n";

//     // Interpret coverage as fraction of total pixel mass if coverage <= 1.
//     double total_mass = std::accumulate(pixel_probs.begin(), pixel_probs.end(), 0.0);
//     double threshold  = (coverage <= 1.0 ? coverage * total_mass : coverage);

//     std::vector<char> covered(pixel_probs.size(), 0); // 0/1 covered HEALPixels
//     std::vector<char> picked(T, 0);

//     auto marginal_gain = [&](size_t tile_idx) -> double {
//         double g = 0.0;
//         for (int pix : member_pixels[tile_idx]) {
//             // Safety: ensure indices are valid (remove assert if you prefer speed)
//             assert(pix >= 0 && static_cast<size_t>(pix) < pixel_probs.size());
//             if (!covered[pix]) g += pixel_probs[pix];
//         }
//         return g;
//     };

//     struct Node {
//         double key;     // last known marginal gain (lazy)
//         size_t idx;     // tile index
//     };
//     struct Cmp {
//         bool operator()(const Node& a, const Node& b) const { return a.key < b.key; } // max-heap
//     };
//     std::priority_queue<Node, std::vector<Node>, Cmp> pq;

//     // Initialize heap with true marginal gains under empty coverage.
//     // (This is O(total pixels across all tiles) once.)
//     for (size_t i = 0; i < T; ++i) {
//         pq.push(Node{marginal_gain(i), i});
//     }

//     std::vector<size_t> selected;
//     selected.reserve(T);

//     double cum = 0.0;
//     const double EPS = 1e-15;

//     while (cum + EPS < threshold) {
//         // Remove already-picked tiles from top
//         while (!pq.empty() && picked[pq.top().idx]) pq.pop();
//         if (pq.empty()) break;

//         Node top = pq.top();
//         pq.pop();
//         if (picked[top.idx]) continue;

//         // Recompute true marginal gain under *current* covered set
//         double true_gain = marginal_gain(top.idx);

//         // Peek next best (lazy key) among unpicked
//         while (!pq.empty() && picked[pq.top().idx]) pq.pop();
//         double next_key = pq.empty() ? -std::numeric_limits<double>::infinity() : pq.top().key;

//         // Lazy-greedy check: if it's still best, accept; otherwise update and retry later
//         if (true_gain + EPS >= next_key) {
//             // Accept this tile
//             picked[top.idx] = 1;
//             selected.push_back(top.idx);

//             // Update coverage + cum by newly covered pixels
//             for (int pix : member_pixels[top.idx]) {
//                 if (!covered[pix]) {
//                     covered[pix] = 1;
//                     cum += pixel_probs[pix];
//                 }
//             }
//         } else {
//             // Not actually best anymore; push back with updated key
//             pq.push(Node{true_gain, top.idx});
//         }
//     }

//     // If nothing selected, early exit
//     if (selected.empty()) return 0;

//     // Reorder + trim vectors to keep only selected tiles in the selection order.
//     std::vector<double> new_tile_probs; new_tile_probs.reserve(selected.size());
//     std::vector<int>    new_ranks;      new_ranks.reserve(selected.size());
//     std::vector<double> new_ras;        new_ras.reserve(selected.size());
//     std::vector<double> new_decs;       new_decs.reserve(selected.size());
//     std::vector<std::vector<int>> new_member_pixels; new_member_pixels.reserve(selected.size());

//     for (size_t idx : selected) {
//         new_tile_probs.push_back(tile_probs[idx]);
//         new_ranks.push_back(ranks[idx]);
//         new_ras.push_back(ras[idx]);
//         new_decs.push_back(decs[idx]);
//         new_member_pixels.push_back(std::move(member_pixels[idx]));
//     }

//     tile_probs     = std::move(new_tile_probs);
//     ranks          = std::move(new_ranks);
//     ras            = std::move(new_ras);
//     decs           = std::move(new_decs);
//     member_pixels  = std::move(new_member_pixels);

//     std::cout << "covered mass (union over pixels): " << cum << " / " << total_mass << "\n";
//     std::cout << "selected tiles: " << selected.size() << "\n";

//     return static_cast<int>(selected.size());
// }


int getTopTileOverlap(
    std::vector<double>& tile_probs,
    std::vector<int>& ranks,
    std::vector<double>& ras,
    std::vector<double>& decs,
    std::vector<std::vector<int>>& member_pixels,
    const std::vector<double>& pixel_probs,
    double coverage) {

    const std::size_t tile_count = tile_probs.size();

    if (ranks.size() != tile_count ||
        ras.size() != tile_count ||
        decs.size() != tile_count ||
        member_pixels.size() != tile_count) {
        throw std::invalid_argument(
            "getTopTileOverlap: inconsistent tile-array sizes.");
    }

    if (tile_count == 0) {
        return 0;
    }

    if (pixel_probs.empty()) {
        throw std::invalid_argument(
            "getTopTileOverlap: pixel probability map is empty.");
    }

    // coverage is an absolute probability, e.g. 0.99.
    if (!std::isfinite(coverage) || coverage <= 0.0) {
        throw std::invalid_argument(
            "getTopTileOverlap: coverage must be positive and finite.");
    }

    for (double probability : pixel_probs) {
        if (!std::isfinite(probability) || probability < 0.0) {
            throw std::invalid_argument(
                "getTopTileOverlap: pixel probabilities must be "
                "finite and nonnegative.");
        }
    }

    // Normalize each tile's pixel list and recompute its standalone
    // probability. Pixels may appear in multiple different tiles.
    for (std::size_t tile = 0; tile < tile_count; ++tile) {
        auto& pixels = member_pixels[tile];

        std::sort(pixels.begin(), pixels.end());
        pixels.erase(
            std::unique(pixels.begin(), pixels.end()),
            pixels.end());

        double probability = 0.0;
        for (int pixel : pixels) {
            if (pixel < 0 ||
                static_cast<std::size_t>(pixel) >= pixel_probs.size()) {
                throw std::out_of_range(
                    "getTopTileOverlap: invalid HEALPixel index.");
            }

            probability += pixel_probs[pixel];
        }

        tile_probs[tile] = probability;
    }

    // Rank tiles by standalone probability. This ranking is used only
    // to construct the target HEALPixel set.
    std::vector<std::size_t> order(tile_count);
    std::iota(order.begin(), order.end(), 0);

    std::sort(
        order.begin(), order.end(),
        [&](std::size_t lhs, std::size_t rhs) {
            if (tile_probs[lhs] != tile_probs[rhs]) {
                return tile_probs[lhs] > tile_probs[rhs];
            }

            // Deterministic tie-breaking.
            return lhs < rhs;
        });

    // Stage 1: construct the target HEALPixel set. Each pixel contributes
    // to target_probability at most once, regardless of tile overlap.
    std::vector<char> target_pixel(pixel_probs.size(), 0);

    double target_probability = 0.0;
    std::size_t seed_tile_count = 0;
    constexpr double epsilon = 1e-12;

    for (std::size_t tile : order) {
        ++seed_tile_count;

        for (int pixel : member_pixels[tile]) {
            if (!target_pixel[pixel]) {
                target_pixel[pixel] = 1;
                target_probability += pixel_probs[pixel];
            }
        }

        if (target_probability + epsilon >= coverage) {
            break;
        }
    }

    if (target_probability + epsilon < coverage) {
        std::cerr
            << "getTopTileOverlap: requested coverage "
            << coverage
            << " is unreachable; coverable probability is "
            << target_probability << ".\n";
    }

    // Stage 2: retain every tile that covers at least one target pixel.
    // This includes alternative overlapping tiles beyond the seed cutoff.
    std::vector<std::size_t> retained_indices;
    retained_indices.reserve(tile_count);

    std::vector<char> retained_pixel(pixel_probs.size(), 0);
    double retained_union_probability = 0.0;

    for (std::size_t tile : order) {
        bool touches_target = false;

        for (int pixel : member_pixels[tile]) {
            if (target_pixel[pixel]) {
                touches_target = true;
                break;
            }
        }

        if (!touches_target) {
            continue;
        }

        retained_indices.push_back(tile);

        // Diagnostic union probability of all retained candidates.
        for (int pixel : member_pixels[tile]) {
            if (!retained_pixel[pixel]) {
                retained_pixel[pixel] = 1;
                retained_union_probability += pixel_probs[pixel];
            }
        }
    }

    if (retained_indices.empty()) {
        throw std::runtime_error(
            "getTopTileOverlap: no tile covers a target HEALPixel.");
    }

    // Reorder and trim all tile-indexed arrays together.
    std::vector<double> retained_tile_probs;
    std::vector<int> retained_ranks;
    std::vector<double> retained_ras;
    std::vector<double> retained_decs;
    std::vector<std::vector<int>> retained_member_pixels;

    retained_tile_probs.reserve(retained_indices.size());
    retained_ranks.reserve(retained_indices.size());
    retained_ras.reserve(retained_indices.size());
    retained_decs.reserve(retained_indices.size());
    retained_member_pixels.reserve(retained_indices.size());

    for (std::size_t tile : retained_indices) {
        retained_tile_probs.push_back(tile_probs[tile]);
        retained_ranks.push_back(ranks[tile]);
        retained_ras.push_back(ras[tile]);
        retained_decs.push_back(decs[tile]);
        retained_member_pixels.push_back(
            std::move(member_pixels[tile]));
    }

    tile_probs = std::move(retained_tile_probs);
    ranks = std::move(retained_ranks);
    ras = std::move(retained_ras);
    decs = std::move(retained_decs);
    member_pixels = std::move(retained_member_pixels);

    std::cout
        << "seed tiles scanned: " << seed_tile_count << "\n"
        << "target unique pixel probability: "
        << target_probability << "\n"
        << "tiles retained: " << retained_indices.size() << "\n"
        << "retained-tile union probability: "
        << retained_union_probability << "\n";

    return static_cast<int>(retained_indices.size());
}

// read rank, ra, dec and prize from Shaon tiling file
void read_data_from_file(const std::string& filename, std::vector<double>& probability, 
                        std::vector<int>& ranks, std::vector<double>& ra, std::vector<double>& dec) {

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line, token;
    getline(file, line);  // Skip header

    while (getline(file, line)) {
        std::stringstream ss(line);
        std::vector<std::string> row;

        while (getline(ss, token, ',')) {
            row.push_back(token);
        }

        ranks.push_back(std::stoi(row[0])); 
        ra.push_back(std::stod(row[2]));  
        dec.push_back(std::stod(row[3]));  
        probability.push_back(std::stod(row[4]));  
    }
    file.close();

}


