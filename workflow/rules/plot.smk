sys.path.append("workflow/scripts")
from evaluate_utils import cal_accuracy, cal_accuracy_multiple
from plotting import plot_prc, plot_confusion_matrices


data_folder = main_config["test_data_folder"]
scenario_list = list(scenario_configs.keys())
threshold_list = sprime_config["threshold_list"]

plot_params = dict()
plot_params["extra_text"] = main_config["extra_text"]
plot_params["models_for_plotting"] = main_config["models_for_plotting"]


import os

base_path = main_config["results_folder"]
plot_inputs = [
    f"{base_path}/{model}/{model}_accuracy_aggregated.txt"
    for model in main_config["models_for_plotting"]
]


plot_inputs_confusion_matrices = [
    f"{base_path}/{model}/{model}_metrics_aggregated.txt"
    for model in main_config["models_for_plotting"]
]



rule plot_prc:
    input:
        f"{base_path}/{{model}}/{{model}}_accuracy_aggregated.txt"
    output:
        f"{base_path}/plots/{{model}}_precision_recall_curve.pdf"
    run:
        plot_prc(input[0], output[0], plot_params)

rule plot_confusion_matrices:
    input:
        f"{base_path}/{{model}}/{{model}}_metrics_aggregated.txt"
    output:
        f"{base_path}/plots/{{model}}_confusion_matrix.pdf"
    run:
        plot_confusion_matrices(input[0], output[0])


rule plot_all:
    input:
        expand(
            f"{base_path}/plots/{{model}}_precision_recall_curve.pdf",
            model=main_config["models_for_plotting"]
        ),
        expand(
        f"{base_path}/plots/{{model}}_confusion_matrix.pdf",
        model=main_config["models_for_plotting"]
    )