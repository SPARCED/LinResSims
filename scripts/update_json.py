#!/bin/bash python3

import json 
import argparse

parser = argparse.ArgumentParser(prog="UpdateJson", description="Intended only for Updating the drs_SPARCED.json\
                                during replication of figure 2, panels D, E, F, and G.")

parser.add_argument('-p', '--path', metavar='path', 
                    default='../sim_configs/drs_SPARCED.json', help='Path to the json configurartion file.')
parser.add_argument('--drug', metavar='drug', help='Drug to perturb cell population with.')
parser.add_argument('--dose', metavar='dose', help='Drug dose to be set.')

def change_drug_dose_settings(path, drug, dose):
    """Modifies the json file at Path to the user-provided drug and dose arguements provided."""\
    
    with open(path, encoding='utf-8', mode='r') as file:
        config = json.load(file)

    config['drs']['drug'] = drug
    config['drs']['dose'] = dose

    with open(path, encoding='utf-8', mode='w') as file:
        json.dump(config, file, indent=2)

    print(f"JSON file was updated to drug {drug} and dose {dose}")


if __name__ == '__main__':

    args = parser.parse_args()

    change_drug_dose_settings(args.path, args.drug, args.dose)