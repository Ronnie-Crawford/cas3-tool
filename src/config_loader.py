import json

def read_config(config_path : str) -> dict:

    """
    Reads a JSON configuration file and returns a dictionary of configuration
    parameters.

    Arguments:
        - config_path (str): The path to the JSON configuration file.

    Returns:
        - config (dict): A dictionary containing the configuration parameters.
    """

    with open(config_path, "r") as config_file:

        config = json.load(config_file)

    return config

# Load the conifg file once
config_path = "./config.json"
config = read_config(config_path)
