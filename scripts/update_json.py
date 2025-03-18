#!/usr/bin/env python3

import json
import argparse
import os

parser = argparse.ArgumentParser(prog="UpdateJson", 
                                 description="Updates the drs_SPARCED.json during replication of Figure 2 (D-G).")

parser.add_argument('-p', '--path', metavar='path', default='../sim_configs/drs_SPARCED.json', 
                    help='Path to the JSON configuration file.')
parser.add_argument('--drug', metavar='drug', required=True, help='Drug to perturb cell population with.')
parser.add_argument('--dose', metavar='dose', required=True, help='Drug dose to be set.')

def change_drug_dose_settings(path, drug, dose):
    """Modifies the JSON file at path to the user-provided drug and dose arguments."""

    # Ensure file exists
    if not os.path.exists(path):
        print(f"Error: JSON file '{path}' not found.")
        exit(1)

    # Load JSON safely
    try:
        with open(path, encoding='utf-8', mode='r') as file:
            config = json.load(file)
    except json.JSONDecodeError:
        print(f"Error: The JSON file '{path}' is not properly formatted.")
        exit(1)

    # Ensure 'drs' key exists
    if "drs" not in config:
        print("Error: 'drs' key not found in the JSON configuration.")
        exit(1)

    # Update values
    config["drs"]["drug"] = drug
    config["drs"]["dose_um"] = dose

    # Save updated JSON
    with open(path, encoding='utf-8', mode='w') as file:
        json.dump(config, file, indent=2)

    print(f"JSON file updated: Drug = {drug}, Dose = {dose}")


if __name__ == '__main__':
    args = parser.parse_args()
    change_drug_dose_settings(args.path, args.drug, args.dose)
