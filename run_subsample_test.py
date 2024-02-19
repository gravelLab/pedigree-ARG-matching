import os
import shutil
from run_alignment import *
from utility import *

os.chdir("pedigrees")
simulation = get_directory_path("Specify the simulation directory in the pedigrees folder:")
result_directory = input("Specify the name of the directory where the results will be saved:")
os.chdir(simulation)
subsample_directory_average_results = dict()
dictionary_result_file = input("Specify the name of the resulting dictionary:")
result_file = open(dictionary_result_file, 'w')
for subsample_test_directory in sorted(os.listdir(), key=len):
    os.chdir(subsample_test_directory)
    total_results = 0
    total_samples = 0
    for separate_simulation in os.listdir():
        total_samples += 1
        if os.path.exists(result_directory) and os.path.isdir(result_directory):
            shutil.rmtree(result_directory)
        print(f"Running the alignment on {subsample_test_directory}/{separate_simulation}")
        run_alignment_and_save_results(directory=separate_simulation, result_directory_name=result_directory)
        os.chdir(separate_simulation)
        clade_names = [file for file in os.listdir() if os.path.isfile(file) and not file.endswith('.pedigree')]
        assert len(clade_names) == 1
        clade_name = clade_names[0]
        os.chdir(result_directory)
        os.chdir(clade_name)
        clade_roots = [file for file in os.listdir() if os.path.isdir(file) and file.isdigit()]
        assert len(clade_roots) == 1
        clade_root = clade_roots[0]
        os.chdir(clade_root)
        simulation_results = len(os.listdir()) - 1
        assert simulation_results > 1
        total_results += simulation_results
        print(f"{simulation_results} alignments")
        os.chdir("../../../..")
        print("Done")
    subsample_directory_average_results[subsample_test_directory] = total_results / total_samples
    print(f"{subsample_test_directory}: {total_results / total_samples}")
    result_file.write(f"{subsample_test_directory}: {total_results / total_samples}\n")
    result_file.flush()
    os.chdir("..")

for key, value in subsample_directory_average_results.items():
    print(f"{key}: {value}")
result_file.close()
# save_dictionary_to_file(dictionary=subsample_directory_average_results, dictionary_filename=dictionary_result_file)
