#!/usr/bin/env python
import os
import earthaccess
import re
import json
import argparse


def convert_to_aws_json_array(env_var_dict):
    json_array = [{"name": key, "value": value} for key, value in env_var_dict.items()]
    return json_array


def get_aws_credentials():

    if auth_object is None:
        print("Failed to authenticate with Earthdata.")
        return None

    try:
        cmr_token = os.environ["CMR_TOKEN"]
        auth_object = earthaccess.Auth()
        auth_object.token = {"access_token": cmr_token}
    except KeyError:
        auth_object = earthaccess.login()
        cmr_token = auth_object.token
        env_var_credentials = {"CMR_TOKEN": cmr_token["access_token"]}

    credentials = earthaccess.get_s3_credentials("ASF")
    for key, value in credentials.items():
        # env_var_key = re.sub("([a-z0-9])([A-Z])", r"\1_\2", key).upper()
        # env_var_key = "AWS_" + re.sub("([a-z0-9])([A-Z])", r"\1_\2", key).upper()
        env_var_key = key
        env_var_credentials[env_var_key] = value

    return env_var_credentials


def update_job_definition(input_filename, output_filename):
    with open(input_filename, "r") as f:
        job_definition = json.load(f)

    env_var_credentials = get_aws_credentials()

    if env_var_credentials:
        new_env_vars = convert_to_aws_json_array(env_var_credentials)
        if "containerOverrides" not in job_definition:
            job_definition["containerOverrides"] = {}
        if "environment" not in job_definition["containerOverrides"]:
            job_definition["containerOverrides"]["environment"] = []
        job_definition["containerOverrides"]["environment"].extend(new_env_vars)

        with open(output_filename, "w") as f:
            json.dump(job_definition, f, indent=2)

        print(f"Updated job definition saved to {output_filename}")
    else:
        print("Failed to update job definition due to missing AWS credentials.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add S3 environment variables to AWS Batch job definition."
    )
    parser.add_argument("input_file", help="Input AWS Batch job definition JSON file.")
    parser.add_argument(
        "output_file",
        help=(
            "Output AWS Batch job definition JSON file with updated environment"
            " variables."
        ),
    )

    args = parser.parse_args()
    update_job_definition(args.input_file, args.output_file)
