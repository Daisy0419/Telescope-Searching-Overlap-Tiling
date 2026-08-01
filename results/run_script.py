import os
import subprocess


def run(executable, dataset, tiling, budgets,
        w_max, w_acc, dwell_time,
        is_deepslow, result_file):

    for budget in budgets:
        print(
            f"\nRunning {executable} with "
            f"dataset={dataset}, "
            f"tiling={tiling}, "
            f"budget={budget}, "
            f"w_max={w_max}, "
            f"w_acc={w_acc}, "
            f"dwell_time={dwell_time}, "
            f"is_deepslow={is_deepslow}"
        )

        # Must match the C++ argument order.
        cmd = [
            executable,
            dataset,
            tiling,
            str(budget),
            str(w_max),
            str(w_acc),
            str(dwell_time),
            str(int(is_deepslow)),
            result_file,
        ]

        try:
            result = subprocess.run(
                cmd,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )

            if result.stdout:
                print(result.stdout)

            if result.stderr:
                print(result.stderr)

        except subprocess.CalledProcessError as error:
            print(f"Error running command: {cmd}")

            if error.stdout:
                print(error.stdout)

            if error.stderr:
                print(error.stderr)


def get_files(folder_path):
    return sorted(
        os.path.splitext(filename)[0]
        for filename in os.listdir(folder_path)
        if filename.endswith(".txt")
    )


def get_tile_files(folder_path):
    return sorted(
        filename
        for filename in os.listdir(folder_path)
        if filename.endswith(".csv")
    )


if __name__ == "__main__":
    executable = "../build/ts"
    map_path = "../data/test_maps"
    tile_path = "../data/tilings"
    result_path = "../results"

    os.makedirs(result_path, exist_ok=True)

    skymaps = get_files(map_path)
    tilefiles = get_tile_files(tile_path)

    budgets = [100]
    w_max = 10
    w_acc = 10
    dwell_time = 1
    is_deepslow = 0

    for tiling_name in tilefiles:
        tiling = os.path.join(tile_path, tiling_name)
        tile_stem = os.path.splitext(tiling_name)[0]

        # One result file per tiling.
        result_file = os.path.join(
            result_path,
            f"maxp_{tile_stem}.csv",
        )

        for skymap in skymaps:
            dataset = os.path.join(
                map_path,
                f"{skymap}.txt",
            )

            print(f"dataset: {dataset}")

            run(
                executable,
                dataset,
                tiling,
                budgets,
                w_max,
                w_acc,
                dwell_time,
                is_deepslow,
                result_file,
            )