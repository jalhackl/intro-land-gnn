import yaml
import os

# Load all configuration files
def load_configs():
    """Load all configuration files and return as dict"""
    configs = {}

    with open("config/main.yaml", "r") as f:
        configs["main"] = yaml.safe_load(f)

    with open("config/hmmix/config_hmmix.yaml", "r") as f:
        configs["hmmix"] = yaml.safe_load(f)

    with open("config/sprime/config_sprime.yaml", "r") as f:
        configs["sprime"] = yaml.safe_load(f)
    
    with open("config/sstar/config_sstar.yaml", "r") as f:
        configs["sstar"] = yaml.safe_load(f)


    with open("config/slurm/config.yaml", "r") as f:
        configs["slurm"] = yaml.safe_load(f)

    if "scenario_subfolder" in configs["main"]:
        for scenario_config in configs["main"]["scenarios"]:
            with open(os.path.join("config", "scenarios", configs["main"]["scenario_subfolder"], scenario_config + ".yaml")) as f:
                configs[scenario_config] = yaml.safe_load(f)
    else:
        for scenario_config in configs["main"]["scenarios"]:
            with open(os.path.join("config", "scenarios", scenario_config + ".yaml")) as f:
                configs[scenario_config] = yaml.safe_load(f)


    return configs


def load_scenario_configs(main_config):
    scenario_configs = {}
    for scenario_config in main_config["scenarios"]:
        with open(os.path.join("config", "scenarios", scenario_config + ".yaml")) as f:
            scenario_configs[scenario_config] = yaml.safe_load(f)

    return scenario_configs



# Load configs
configs = load_configs()
main_config = configs["main"]
hmmix_config = configs["hmmix"]
sprime_config = configs["sprime"]
slurm_config = configs["slurm"]

sstar_config = configs["sstar"]

#scenario_configs = [configss[config] for config in configs if config.startswith("config_")]
scenario_configs = load_scenario_configs(main_config)
